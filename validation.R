# ============================================
# External Validation: Centroid-based Classification
# ============================================

library(dplyr)
library(tidyr)

# ============================================
# 1. Centroid 계산 (Internal)
# ============================================

# View 1: GSVA
C1_gsva <- colMeans(internal_dat1[fit$clusters == 1, ])
C2_gsva <- colMeans(internal_dat1[fit$clusters == 2, ])
C3_gsva <- colMeans(internal_dat1[fit$clusters == 3, ])

# View 2: CNV
C1_cnv <- colMeans(internal_dat2[fit$clusters == 1, ])
C2_cnv <- colMeans(internal_dat2[fit$clusters == 2, ])
C3_cnv <- colMeans(internal_dat2[fit$clusters == 3, ])

# View 3: Mutation
C1_mut <- colMeans(internal_dat3[fit$clusters == 1, ])
C2_mut <- colMeans(internal_dat3[fit$clusters == 2, ])
C3_mut <- colMeans(internal_dat3[fit$clusters == 3, ])

# Centroid list
centroids <- list(
  C1 = list(gsva = C1_gsva, cnv = C1_cnv, mut = C1_mut),
  C2 = list(gsva = C2_gsva, cnv = C2_cnv, mut = C2_mut),
  C3 = list(gsva = C3_gsva, cnv = C3_cnv, mut = C3_mut)
)

cat("✅ Centroids calculated\n")
cat("  - GSVA:", length(C1_gsva), "features\n")
cat("  - CNV:", length(C1_cnv), "features\n")
cat("  - Mutation:", length(C1_mut), "features\n")

# ============================================
# 2. Immune pathways
# ============================================

immune_pathways <- c(
  "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
)

# ============================================
# 3. Thresholds 계산
# ============================================

cluster_features <- data.frame(
  sample = rownames(internal_dat1),
  cluster = fit$clusters,
  
  cnv_gain_magnitude = rowSums(internal_dat2[, 1:50]),
  cnv_loss_magnitude = rowSums(internal_dat2[, 51:100]),
  
  mut_burden = rowSums(internal_dat3),
  
  immune_score = rowMeans(internal_dat1[, immune_pathways])
)

thresholds <- cluster_features %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    
    cnv_gain_median = median(cnv_gain_magnitude),
    cnv_gain_Q1 = quantile(cnv_gain_magnitude, 0.25),
    cnv_gain_Q3 = quantile(cnv_gain_magnitude, 0.75),
    
    cnv_loss_median = median(cnv_loss_magnitude),
    cnv_loss_Q1 = quantile(cnv_loss_magnitude, 0.25),
    cnv_loss_Q3 = quantile(cnv_loss_magnitude, 0.75),
    
    mut_median = median(mut_burden),
    mut_Q1 = quantile(mut_burden, 0.25),
    mut_Q3 = quantile(mut_burden, 0.75),
    
    immune_median = median(immune_score),
    immune_Q1 = quantile(immune_score, 0.25),
    immune_Q3 = quantile(immune_score, 0.75)
  )

cat("✅ Thresholds calculated\n")
print(thresholds)

# ============================================
# 4. 분류 함수
# ============================================

classify_by_correlation <- function(new_gsva, new_cnv, new_mut, centroids) {
  
  new_gsva <- as.numeric(new_gsva)
  new_cnv <- as.numeric(new_cnv)
  new_mut <- as.numeric(new_mut)
  
  # Mutation SD = 0이면 mutation correlation 제외
  if(sd(new_mut) == 0) {
    cor_C1 <- cor(new_gsva, centroids$C1$gsva) + cor(new_cnv, centroids$C1$cnv)
    cor_C2 <- cor(new_gsva, centroids$C2$gsva) + cor(new_cnv, centroids$C2$cnv)
    cor_C3 <- cor(new_gsva, centroids$C3$gsva) + cor(new_cnv, centroids$C3$cnv)
  } else {
    cor_C1 <- cor(new_gsva, centroids$C1$gsva) + 
              cor(new_cnv, centroids$C1$cnv) + 
              cor(new_mut, centroids$C1$mut)
    
    cor_C2 <- cor(new_gsva, centroids$C2$gsva) + 
              cor(new_cnv, centroids$C2$cnv) + 
              cor(new_mut, centroids$C2$mut)
    
    cor_C3 <- cor(new_gsva, centroids$C3$gsva) + 
              cor(new_cnv, centroids$C3$cnv) + 
              cor(new_mut, centroids$C3$mut)
  }
  
  correlations <- c(C1 = cor_C1, C2 = cor_C2, C3 = cor_C3)
  assigned_cluster <- names(which.max(correlations))
  
  return(list(
    cluster = assigned_cluster,
    correlations = correlations
  ))
}

validate_by_rules <- function(new_gsva, new_cnv, new_mut, 
                              assigned_cluster, thresholds,
                              immune_pathways, gain_cols, loss_cols) {
  
  new_cnv_gain_mag <- sum(new_cnv[gain_cols], na.rm=TRUE)
  new_cnv_loss_mag <- sum(new_cnv[loss_cols], na.rm=TRUE)
  new_mut_burden <- sum(new_mut, na.rm=TRUE)
  
  # Immune score (as.numeric 필수!)
  new_immune <- mean(as.numeric(new_gsva[immune_pathways]), na.rm=TRUE)
  
  # "C1" → 1
  cluster_idx <- as.numeric(gsub("C", "", assigned_cluster))
  
  cnv_gain_match <- (new_cnv_gain_mag >= thresholds$cnv_gain_Q1[cluster_idx] & 
                     new_cnv_gain_mag <= thresholds$cnv_gain_Q3[cluster_idx])
  
  cnv_loss_match <- (new_cnv_loss_mag >= thresholds$cnv_loss_Q1[cluster_idx] & 
                     new_cnv_loss_mag <= thresholds$cnv_loss_Q3[cluster_idx])
  
  mut_match <- (new_mut_burden >= thresholds$mut_Q1[cluster_idx] & 
                new_mut_burden <= thresholds$mut_Q3[cluster_idx])
  
  immune_match <- (new_immune >= thresholds$immune_Q1[cluster_idx] & 
                   new_immune <= thresholds$immune_Q3[cluster_idx])
  
  confidence <- mean(c(cnv_gain_match, cnv_loss_match, mut_match, immune_match), na.rm=TRUE)
  
  return(list(
    cnv_gain_magnitude = new_cnv_gain_mag,
    cnv_loss_magnitude = new_cnv_loss_mag,
    mut_burden = new_mut_burden,
    immune_score = new_immune,
    confidence = confidence
  ))
}

classify_external_sample <- function(new_gsva, new_cnv, new_mut, 
                                     centroids, thresholds, 
                                     immune_pathways, gain_cols, loss_cols) {
  
  cor_result <- classify_by_correlation(new_gsva, new_cnv, new_mut, centroids)
  rule_result <- validate_by_rules(new_gsva, new_cnv, new_mut, 
                                   cor_result$cluster, thresholds,
                                   immune_pathways, gain_cols, loss_cols)
  
  return(list(
    assigned_cluster = cor_result$cluster,
    correlations = cor_result$correlations,
    features = rule_result,
    confidence = rule_result$confidence
  ))
}

# ============================================
# 5. External 분류
# ============================================

cat("\n🚀 Classifying external samples...\n")

external_results <- list()

for(i in 1:nrow(external_dat1)) {
  
  result <- classify_external_sample(
    new_gsva = external_dat1[i, ],
    new_cnv = external_dat2[i, ],
    new_mut = external_dat3[i, ],
    centroids = centroids,
    thresholds = thresholds,
    immune_pathways = immune_pathways,
    gain_cols = 1:50,
    loss_cols = 51:100
  )
  
  external_results[[i]] <- result
  
  if(i %% 10 == 0) cat("  Processed", i, "/ 88 samples\n")
}

# ============================================
# 6. 결과 정리
# ============================================

external_clusters <- data.frame(
  sample = rownames(external_dat1),
  assigned_cluster = sapply(external_results, function(x) x$assigned_cluster),
  confidence = sapply(external_results, function(x) x$confidence),
  cor_C1 = sapply(external_results, function(x) x$correlations["C1"]),
  cor_C2 = sapply(external_results, function(x) x$correlations["C2"]),
  cor_C3 = sapply(external_results, function(x) x$correlations["C3"])
)

cat("\n=== EXTERNAL CLASSIFICATION RESULTS ===\n")
table(external_clusters$assigned_cluster)
cat("\nConfidence:\n")
summary(external_clusters$confidence)

# ============================================
# 7. Internal vs External 비교
# ============================================

# Internal features
internal_features <- data.frame(
  sample = rownames(internal_dat1),
  cluster = paste0("C", fit$clusters),
  
  cnv_gain = rowSums(internal_dat2[, 1:50]),
  cnv_loss = rowSums(internal_dat2[, 51:100]),
  mut_burden = rowSums(internal_dat3),
  immune_score = rowMeans(internal_dat1[, immune_pathways])
)

internal_summary <- internal_features %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    cnv_gain_median = median(cnv_gain),
    cnv_loss_median = median(cnv_loss),
    mut_median = median(mut_burden),
    immune_median = median(immune_score)
  ) %>%
  arrange(cluster)

# External features
external_features <- data.frame(
  sample = rownames(external_dat1),
  cluster = external_clusters$assigned_cluster,
  
  cnv_gain = rowSums(external_dat2[, 1:50]),
  cnv_loss = rowSums(external_dat2[, 51:100]),
  mut_burden = rowSums(external_dat3),
  immune_score = rowMeans(external_dat1[, immune_pathways])
)

external_summary <- external_features %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    cnv_gain_median = median(cnv_gain),
    cnv_loss_median = median(cnv_loss),
    mut_median = median(mut_burden),
    immune_median = median(immune_score)
  ) %>%
  arrange(cluster)

# 전체 비교
comparison <- bind_rows(
  internal_summary %>% mutate(cohort = "Internal"),
  external_summary %>% mutate(cohort = "External")
)

cat("\n=== FULL COMPARISON ===\n")
print(comparison[, c("cohort", "cluster", "n", "cnv_gain_median", "cnv_loss_median", "mut_median", "immune_median")])

# C2만 비교
cat("\n=== C2 COMPARISON (Immune Cold) ===\n")
cat("Internal C2:\n")
print(internal_summary %>% filter(cluster == "C2"))
cat("\nExternal C2:\n")
print(external_summary %>% filter(cluster == "C2"))

# 저장
saveRDS(external_clusters, "external_clusters_centroid.rds")
saveRDS(list(
  internal_summary = internal_summary,
  external_summary = external_summary,
  comparison = comparison
), "validation_summary.rds")

cat("\n✅ Validation complete!\n")
cat("📊 Results saved:\n")
cat("  - external_clusters_centroid.rds\n")
cat("  - validation_summary.rds\n")
