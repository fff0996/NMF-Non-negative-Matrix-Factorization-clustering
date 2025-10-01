library(ComplexHeatmap); library(circlize); library(matrixStats); library(RColorBrewer) ; library(IntNMF)
library(maftools)
library(data.table)
library(dplyr)
library(tidyr)

merged_maf <- readRDS("merged_maf.rds")
# merged_maf: maftools::merge_mafs()로 만든 MAF 객체라고 가정
d <- as.data.table(merged_maf@data)

# 비침묵 변이만 사용 (Silent/UTR/Intron 등 제외)
drop_classes <- c("Silent","Intron","IGR","RNA","3'UTR","5'UTR")
nonsyn <- d[!(d$Variant_Classification %in% drop_classes), , drop = FALSE]

# 샘플×유전자 0/1 매트릭스
mut_mat <- nonsyn %>%
  transmute(SampleID = Tumor_Sample_Barcode, Gene = Hugo_Symbol, val = 1L) %>%
  distinct() %>%
  pivot_wider(names_from = Gene, values_from = val, values_fill = 0) %>%
  as.data.frame()

rownames(mut_mat) <- mut_mat$SampleID
mut_mat$SampleID <- NULL


# 1. 각 유전자(열)별로 mutation frequency 계산
gene_freq <- colSums(mut_mat) / nrow(mut_mat)

# 2. 0.05 < freq < 0.60 조건 만족하는 유전자 선택
selected_genes <- gene_freq > 0.05 & gene_freq < 0.60

# 3. 필터링된 매트릭스
mut_mat_filtered <- mut_mat[, selected_genes]

# 확인
dim(mut_mat_filtered)
write.table(mut_mat_filtered,"NMF_mutinput.txt",sep="\t",quote=F,row.names=F)


gsva <- read.table("Hallmark GSVA위치")
cnv <- read.table("gistic 결과 위치")
mut <- read.table("NMF_mutinput.txt")
# 샘플 (행) * features (열) 형태로 맞춰주기 


# common 공통 샘플로 통일
common <- Reduce(intersect,list(rownames(gsva), rownames(cnv),rownames(mut)))
#나중에 히트맵 그릴때 이 데이터로 그림 
gsva_raw <- gsva[common,]
cnv_raw <- cnv[common,]
mut_raw <- mut[common,]

# IntNMF 입력: Gain/Loss 분리 (정확성)
cnv_gain <- pmax(as.matrix(cnv_raw), 0)
cnv_loss <- abs(pmin(as.matrix(cnv_raw), 0))
dat2_combined <- cbind(cnv_gain/max(cnv_gain), cnv_loss/max(cnv_loss))

# GSVA는 shift (feature 수 관리)
dat1 <- as.matrix(gsva_raw)
if (!all(dat1 >= 0)) dat1 <- dat1 + abs(min(dat1)) + .Machine$double.eps
dat1 <- dat1/max(dat1)

# Mutation
dat3 <- as.matrix(mut_raw)
p <- colMeans(dat3 == 1)
dat3 <- dat3[, p >= 0.05 & p <= 0.60, drop=FALSE]

# IntNMF
dat <- list(dat1, dat2_combined, dat3)
fit <- nmf.mnnals(dat=dat, k=3, maxiter=200, st.count=20, n.ini=15, 
                  ini.nndsvd=TRUE, seed=TRUE)
SilhouettePlot(fit, cluster.col = NULL)
ConsensusMatPlot(fit,rowLab=TRUE,colLab=TRUE)
## 히트맵: 원본 사용 (간단하게)
common <- Reduce(intersect, list(rownames(gsva_raw), rownames(cnv_raw), 
                                 rownames(mut_raw), rownames(fit$W)))
W <- fit$W[common, , drop=FALSE]
cl <- setNames(paste0("C", fit$clusters[common]), common)
sord <- names(sort(apply(W, 1, max)))
clord <- cl[sord]
lev <- unique(clord)

# 원본 데이터 (Gain/Loss 합치지 않음)
G <- gsva_raw[sord, , drop=FALSE]
C <- cnv_raw[sord, , drop=FALSE]  # signed 그대로
M <- mut_raw[sord, , drop=FALSE]

# Z-score/clip
zcol <- function(X, cap=2){
  X <- as.matrix(X)
  Z <- scale(X)
  Z[!is.finite(Z)] <- 0
  Z[Z > cap] <- cap
  Z[Z < -cap] <- -cap
  Z
}

Gz <- zcol(G, cap=2)
C <- as.matrix(C)
Ccl <- C
Ccl[Ccl > 2] <- 2
Ccl[Ccl < -2] <- -2
Mt <- as.matrix(M)

# 히트맵 (3개만, gain/loss 따로 안 그림)
col_gsva <- colorRamp2(c(-2,0,2), c("#2b6cb0","white","#c53030"))
col_cnv <- colorRamp2(c(-2,0,2), c("#225ea8","white","#e34a33"))
col_mut <- c("0"="white","1"="black")

ha <- HeatmapAnnotation(Cluster=factor(clord, levels=lev))

hts <- list(
  gsva = Heatmap(t(Gz), name="GSVA", col=col_gsva, top_annotation=ha,
                 show_row_names=TRUE, show_column_names=FALSE,
                 cluster_rows=TRUE, cluster_columns=FALSE,
                 column_split=factor(clord, levels=lev)),
  
  cnv = Heatmap(t(Ccl), name="CNV", col=col_cnv,
                show_row_names=TRUE, show_column_names=FALSE,
                cluster_rows=TRUE, cluster_columns=FALSE,
                column_split=factor(clord, levels=lev)),
  
  mutation = Heatmap(t(Mt), name="Mutation", col=col_mut,
                     show_row_names=TRUE, show_column_names=FALSE,
                     cluster_rows=TRUE, cluster_columns=FALSE,
                     column_split=factor(clord, levels=lev))
)

draw(Reduce(`%v%`, hts), merge_legend=TRUE)


#Save
ht <- Reduce(`%v%`, hts)  # 이미 만든 히트맵 결합

## 크기 자동 산출(대충 가독 사이즈)
ns <- length(clord)                                  # 샘플 수(열)
nr <- sum(c(ncol(Mt), ncol(Ccl), ncol(Gz)), na.rm=TRUE)  # 총 피처 수(행)
w_in <- max(6,  min(0.12 * ns, 30))  # 가로(inch)
h_in <- max(8,  min(0.15 * nr, 60))  # 세로(inch)

##anno 데이터가 있다는 전제
## 0) 체크
stopifnot(nrow(Gz)==nrow(Ccl), nrow(Gz)==nrow(Mt))
samp <- rownames(Gz)
annx <- anno[match(samp, anno$sample), ]
stopifnot(identical(annx$sample, samp))
lev <- c("C1","C2","C3")

## 1) 분할·상태
split_fac <- factor(annx$clord, levels=lev)
ord <- order(split_fac)  # 군집 블록 고정
ec_state <- ifelse(annx$ecDNA_called & annx$ecDNA_any, "pos",
             ifelse(annx$ecDNA_called & !annx$ecDNA_any, "neg", "NA"))
ec_state <- factor(ec_state, levels=c("neg","pos","NA"))

annx$categ <- as.character(annx$categ)
annx$categ[is.na(annx$categ)] <- "NA"  # NA를 문자열로

# Factor로 변환 (모든 level 명시)
categ_levels <- c("Signature_3", "Signature_APOBEC", "Signature_clock", 
                  "Signature_17", "Signature_8", "NA")
annx$categ <- factor(annx$categ, levels = categ_levels)

# 색상 (순서 동일하게)
categ_colors <- c(
  "Signature_3" = "#E41A1C",
  "Signature_APOBEC" = "#377EB8",
  "Signature_clock" = "#999999",
  "Signature_17" = "#999999",
  "Signature_8" = "#999999",
  "NA" = "#999999"
)
## 2) 상단 어노테이션
library(ComplexHeatmap); library(grid)
ha_top <- HeatmapAnnotation(
  df = data.frame(
    Cluster = split_fac,
    APOBEC  = annx$categ,  # exclude=NULL 제거
    ecDNA   = ec_state,
    driver  = factor(annx$ecDNA_driver, exclude=NULL),
    row.names = samp
  ),
  TMB    = anno_barplot(annx$TMB, border=FALSE),
  burden = anno_barplot(annx$ecDNA_burden, border=FALSE),
  col = list(
    Cluster = setNames(c("#1f77b4","#2ca02c","#ff7f0e"), lev),
    APOBEC  = categ_colors,
    ecDNA   = c(neg="#BDBDBD", pos="#000000", "NA"="#FFFFFF")
  ),
  annotation_name_side="left",
  show_annotation_name=TRUE
)


## 3) 히트맵(행=샘플이므로 t() 유지)
hts <- list(
  gsva = Heatmap(t(Gz),  name="GSVA",      col=col_gsva, top_annotation=ha_top,
                 show_row_names=TRUE, show_column_names=FALSE,
                 cluster_rows=TRUE,  cluster_columns=FALSE,
                 column_split=split_fac, column_order=ord),
  cnv  = Heatmap(t(Ccl), name="CopyNumber", col=col_cnv,
                 show_row_names=TRUE, show_column_names=FALSE,
                 cluster_rows=TRUE,  cluster_columns=FALSE,
                 column_split=split_fac, column_order=ord),
  mutation = Heatmap(t(Mt), name="Mutation", col=col_mut,
                 show_row_names=TRUE, show_column_names=FALSE,
                 cluster_rows=TRUE,  cluster_columns=FALSE,
                 column_split=split_fac, column_order=ord)
)

draw(Reduce(`%v%`, hts), merge_legend=TRUE,
     heatmap_legend_side="right", annotation_legend_side="right")

## PDF
pdf("intNMF_heatmap.pdf", width=w_in, height=h_in, onefile=FALSE)
draw(ht, merge_legend=TRUE, heatmap_legend_side="right", annotation_legend_side="right")
dev.off()

## PNG(고해상도)
png("intNMF_heatmap.png", units="in", res=300, width=w_in, height=h_in, type="cairo")
draw(ht, merge_legend=TRUE, heatmap_legend_side="right", annotation_legend_side="right")
dev.off()
