# ============================================
# Step 1: PyClone Input 생성 (완성 버전)
# ============================================

library(VariantAnnotation)
library(GenomicRanges)
library(dplyr)

# ============================================
# 경로 설정
# ============================================

vcf_dir <- "/BiO/Hyein/M1M2VCF/somatic/coding/"
cns_dir <- "/BiO/Hyein/cnv/results/"
output_dir <- "pyclone_inputs/"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================
# Patient 리스트 생성
# ============================================

vcf_files <- list.files(vcf_dir, 
                        pattern = "\\.PASS\\.filtered3\\.vep\\.snpeff\\.vcf$",
                        full.names = FALSE)

# Patient ID 추출 함수
extract_patient_id <- function(filename) {
  id <- gsub("\\.PASS\\.filtered3.*", "", filename)
  return(id)
}

all_patient_ids <- sapply(vcf_files, extract_patient_id, USE.NAMES = FALSE)

print(paste("Total patients:", length(all_patient_ids)))

# ============================================
# 한 환자 처리 함수
# ============================================

process_patient <- function(patient_id) {
  
  cat("\n========================================\n")
  cat("Processing:", patient_id, "\n")
  cat("========================================\n")
  
  tryCatch({
    
    # 파일 찾기
    vcf_file <- paste0(vcf_dir, patient_id, ".vcf")
    cns_file <- paste0(cns_dir, patient_id, ".markdup.sorted.call.cns")
    
    # 파일 존재 확인
    if(!file.exists(vcf_file)) {
      stop(paste("VCF not found:", vcf_file))
    }
    if(!file.exists(cns_file)) {
      stop(paste("CNS not found:", cns_file))
    }
    
    # VCF 읽기
    vcf <- readVcf(vcf_file, genome="hg38")
    
    sample_names <- colnames(geno(vcf)$AD)
    tumor_idx <- grep("TD", sample_names)
    if(length(tumor_idx) == 0) tumor_idx <- 2
    tumor_idx <- tumor_idx[1]
    
    # Mutations 추출
    chrom <- as.character(seqnames(vcf))
    pos <- start(vcf)
    ref <- as.character(ref(vcf))
    alt <- as.character(unlist(alt(vcf)))
    
    ad <- geno(vcf)$AD[, tumor_idx]
    ref_counts <- sapply(ad, function(x) x[1])
    alt_counts <- sapply(ad, function(x) x[2])
    vaf <- alt_counts / (ref_counts + alt_counts)
    
    mutations <- data.frame(
      mutation_id = paste0(chrom, "_", pos, "_", ref, "_", alt),
      chrom = chrom,
      pos = pos,
      ref_counts = as.numeric(ref_counts),
      var_counts = as.numeric(alt_counts),
      vaf = vaf,
      stringsAsFactors = FALSE
    )
    
    cat("  VCF mutations:", nrow(mutations), "\n")
    
    # CNS 읽기
    cns <- read.table(cns_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    cat("  CNS segments:", nrow(cns), "\n")
    
    # Chromosome 맞추기
    if(grepl("^chr", cns$chromosome[1]) && !grepl("^chr", mutations$chrom[1])) {
      cns$chromosome <- gsub("^chr", "", cns$chromosome)
    } else if(!grepl("^chr", cns$chromosome[1]) && grepl("^chr", mutations$chrom[1])) {
      cns$chromosome <- paste0("chr", cns$chromosome)
    }
    
    # CN 정보 추가
    mutations_gr <- GRanges(
      seqnames = mutations$chrom,
      ranges = IRanges(start = mutations$pos, width = 1)
    )
    
    cns_gr <- GRanges(
      seqnames = cns$chromosome,
      ranges = IRanges(start = cns$start, end = cns$end),
      cn = cns$cn
    )
    
    overlaps <- findOverlaps(mutations_gr, cns_gr)
    
    mutations$total_cn <- 2
    mutations$total_cn[queryHits(overlaps)] <- cns_gr$cn[subjectHits(overlaps)]
    
    mutations <- mutations %>%
      mutate(
        minor_cn = case_when(
          total_cn == 0 ~ 0,
          total_cn == 1 ~ 0,
          total_cn >= 2 ~ pmin(1, floor(total_cn / 2))
        ),
        major_cn = total_cn - minor_cn
      )
    
    # Purity 추정
    diploid <- mutations %>%
      filter(minor_cn == 1, major_cn == 1, 
             var_counts >= 3,
             ref_counts + var_counts >= 10)
    
    purity_diploid <- ifelse(nrow(diploid) > 5,
                              min(quantile(diploid$vaf, 0.75, na.rm=TRUE) * 2, 1.0),
                              NA)
    
    loh <- mutations %>%
      filter(minor_cn == 0, major_cn == 1,
             var_counts >= 3,
             ref_counts + var_counts >= 10)
    
    purity_loh <- ifelse(nrow(loh) > 3,
                         min(quantile(loh$vaf, 0.75, na.rm=TRUE), 1.0),
                         NA)
    
    purities <- c(purity_diploid, purity_loh)
    purities <- purities[!is.na(purities)]
    purity <- ifelse(length(purities) > 0, median(purities), 0.7)
    purity <- min(max(purity, 0.3), 1.0)
    
    cat("  Purity:", round(purity, 3), "\n")
    
    # PyClone Input
    pyclone_input <- mutations %>%
      filter(var_counts >= 3, ref_counts + var_counts >= 10) %>%
      mutate(sample_id = patient_id,
             normal_cn = 2,
             tumour_content = purity) %>%
      dplyr::select(mutation_id, sample_id, ref_counts, var_counts,
                    normal_cn, minor_cn, major_cn, tumour_content)
    
    cat("  Final mutations:", nrow(pyclone_input), "\n")
    
    # 저장
    # Patient ID에서 파일명으로 안전한 문자만 사용
    safe_id <- gsub("[^A-Za-z0-9_-]", "_", patient_id)
    output_file <- paste0(output_dir, safe_id, "_input.tsv")
    
    write.table(pyclone_input, output_file, sep = "\t", 
                quote = FALSE, row.names = FALSE)
    
    cat("  ✓ Saved:", output_file, "\n")
    
    return(list(
      patient_id = patient_id,
      n_mutations = nrow(pyclone_input),
      purity = purity,
      status = "SUCCESS"
    ))
    
  }, error = function(e) {
    cat("  ✗ ERROR:", e$message, "\n")
    return(list(
      patient_id = patient_id,
      n_mutations = NA,
      purity = NA,
      status = paste("ERROR:", e$message)
    ))
  })
}

# ============================================
# 테스트: 첫 환자
# ============================================

cat("\n=== TESTING FIRST PATIENT ===\n")
test_result <- process_patient(all_patient_ids[1])
print(test_result)

# ============================================
# 전체 실행 (주석 해제하여 사용)
# ============================================

results <- list()
for(i in 1:length(all_patient_ids)) {
   result <- process_patient(all_patient_ids[i])
   results[[i]] <- result
 }
 
 # Summary
 results_df <- do.call(rbind, lapply(results, as.data.frame))
 write.csv(results_df, "pyclone_input_summary.csv", row.names = FALSE)
 
 cat("\n=== SUMMARY ===\n")
 cat("Total patients:", nrow(results_df), "\n")
 cat("Success:", sum(results_df$status == "SUCCESS"), "\n")
 cat("Errors:", sum(results_df$status != "SUCCESS"), "\n")

#Pyclone 분석 
pyclone-vi  fit -i YON30-TD-180427_input.tsv -o ./results/ -c 10 -d beta-binomial -r 5

pyclone-vi write-results-file -i ./results -o ./results.tsv
