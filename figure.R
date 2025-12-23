#=============================================================================
# GSVA 그룹별 barplot 그림 
## ===========================================================================
## 저장 없이 화면에만 출력 (print) 버전
## =============================

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

plot_gsva_set_print <- function(GSVA_mat, split_fac, pathways, title=NULL) {

  ## 0) shape 체크
  GSVA_mat <- as.matrix(GSVA_mat)
  stopifnot(!is.null(rownames(GSVA_mat)), !is.null(colnames(GSVA_mat)))
  stopifnot(!is.null(names(split_fac)))
  stopifnot(is.factor(split_fac))

  common <- intersect(rownames(GSVA_mat), names(split_fac))
  stopifnot(length(common) >= 10)

  pathways <- intersect(pathways, colnames(GSVA_mat))
  if (length(pathways) == 0) stop("지정한 pathway가 GSVA_mat colnames에 없음")

  ## 1) long
  df <- GSVA_mat[common, pathways, drop=FALSE] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Sample") %>%
    mutate(Cluster = split_fac[Sample]) %>%
    pivot_longer(-c(Sample, Cluster),
                 names_to="Pathway", values_to="GSVA") %>%
    mutate(
      Cluster = factor(Cluster, levels = levels(split_fac)),
      Pathway = factor(Pathway, levels = pathways),
      Pathway_lab = str_wrap(Pathway, width = 20)
    )

  ## 2) summary
  df_sum <- df %>%
    group_by(Pathway_lab, Cluster) %>%
    summarise(
      mean = mean(GSVA, na.rm=TRUE),
      sem  = sd(GSVA, na.rm=TRUE) / sqrt(sum(!is.na(GSVA))),
      .groups="drop"
    )

  ## 3) plot (저장 X, print만)
  p <- ggplot(df_sum, aes(x=Cluster, y=mean)) +
    geom_col(width=0.7) +
    geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem),
                  width=0.2, linewidth=0.4) +
    geom_jitter(
      data=df,
      aes(x=Cluster, y=GSVA),
      width=0.12, height=0,
      size=0.7, alpha=0.35
    ) +
    facet_wrap(~Pathway_lab, nrow=1, scales="free_y") +
    labs(title=title, x=NULL, y="GSVA score") +
    theme_bw(base_size=10) +
    theme(
      plot.title = element_text(face="bold", size=11),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      strip.text = element_text(size=9),
      axis.text.x = element_text(size=9)
    )

  print(p)
  invisible(list(plot=p, df=df, df_sum=df_sum))
}

## -----------------------------
## 실행 예시 (저장 안 함)
## -----------------------------
# GSVA_mat는 네가 위에서 만들어둔 GSVA_mat 사용 (gsva_raw or Gz 기반)
# split_fac는 C1/C2/C3 factor (names=sample)

pw_immune <- c("INTERFERON GAMMA RESPONSE","INTERFERON ALPHA RESPONSE",
               "TNFA SIGNALING VIA NFKB","IL6 JAK STAT3 SIGNALING","INFLAMMATORY RESPONSE")

pw_metab  <- c("FATTY ACID METABOLISM","BILE ACID METABOLISM",
               "XENOBIOTIC METABOLISM","PEROXISOME","ADIPOGENESIS")
pw_ploi <- c("MYC TARGETS V1","MYC TARGETS V2","E2F TARGETS","G2M CHECKPOINT","MITOTIC SPINDLE")

# 화면 출력
res_immune <- plot_gsva_set_print(GSVA_mat, split_fac, pw_immune, "Immune / inflammatory (C2)")
res_metab  <- plot_gsva_set_print(GSVA_mat, split_fac, pw_metab,  "Metabolic / differentiation (C3)")





#======================================================
#CNV Amp/Del CNV mean 그룹별 
#======================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

## ---- 0) shape check
cnv_cyto <- as.data.frame(cnv_cyto)
stopifnot(!is.null(rownames(cnv_cyto)))
stopifnot(is.factor(split_fac), !is.null(names(split_fac)))

common_samp <- intersect(rownames(cnv_cyto), names(split_fac))
stopifnot(length(common_samp) >= 10)

cnv_cyto  <- cnv_cyto[common_samp, , drop=FALSE]
split_fac <- split_fac[common_samp]

# numeric 강제
cnv_cyto[] <- lapply(cnv_cyto, function(x) as.numeric(as.character(x)))

## ---- 1) stats (Kruskal + BH)
pvals <- sapply(colnames(cnv_cyto), function(band) {
  x <- cnv_cyto[[band]]
  ok <- is.finite(x) & !is.na(split_fac)
  if (length(unique(split_fac[ok])) < 2) return(1)
  suppressWarnings(kruskal.test(x[ok] ~ split_fac[ok])$p.value)
})
qvals <- p.adjust(pvals, "BH")

cnv_stats <- data.frame(
  CNV_peak = colnames(cnv_cyto),
  P_value  = pvals[colnames(cnv_cyto)],
  FDR      = qvals[colnames(cnv_cyto)],
  stringsAsFactors = FALSE
)

sel_bands <- cnv_stats$CNV_peak[is.finite(cnv_stats$FDR) & cnv_stats$FDR < 0.05]
stopifnot(length(sel_bands) >= 1)
cat("FDR<0.05 bands:", length(sel_bands), "\n")

## ---- 2) mean ± SE table for selected bands
calc_mean_se <- function(x, g) {
  df <- data.frame(x=x, g=g)
  df %>%
    group_by(g) %>%
    summarise(
      n    = sum(!is.na(x)),
      mean = mean(x, na.rm=TRUE),
      se   = sd(x, na.rm=TRUE) / sqrt(n),
      .groups = "drop"
    )
}

plot_df <- lapply(sel_bands, function(band){
  ms <- calc_mean_se(cnv_cyto[[band]], split_fac)
  ms$CNV_peak <- band
  ms
}) %>% bind_rows()

plot_df <- plot_df %>%
  dplyr::rename(Cluster = g) %>%
  mutate(
    Cluster = factor(as.character(Cluster), levels = levels(split_fac))
  ) %>%
  dplyr::select(CNV_peak, Cluster, n, mean, se)

## ---- 3) DEL/AMP 분류
mean_overall <- plot_df %>%
  group_by(CNV_peak) %>%
  summarise(mean_all = mean(mean, na.rm=TRUE), .groups="drop") %>%
  mutate(type = ifelse(mean_all < 0, "DEL", "AMP"))

plot_df <- plot_df %>%
  left_join(mean_overall[, c("CNV_peak", "type")], by="CNV_peak")

## ---- 4) DEL/AMP 각각 정렬 (effect size)
order_del <- plot_df %>%
  filter(type == "DEL") %>%
  group_by(CNV_peak) %>%
  summarise(effect = max(mean) - min(mean), .groups="drop") %>%
  arrange(desc(effect))

order_amp <- plot_df %>%
  filter(type == "AMP") %>%
  group_by(CNV_peak) %>%
  summarise(effect = max(mean) - min(mean), .groups="drop") %>%
  arrange(desc(effect))

## ---- 5) DEL plot (가로)
plot_del <- plot_df %>% filter(type == "DEL")
plot_del$CNV_peak <- factor(plot_del$CNV_peak, levels=order_del$CNV_peak)


                     
p_del <- ggplot(plot_del, aes(x=CNV_peak, y=mean, fill=Cluster)) +
  geom_bar(stat="identity", position=position_dodge(0.8), width=0.7) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                position=position_dodge(0.8), width=0.25, linewidth=0.4) +
  scale_fill_manual(values=c("C1"="#8EB3D6", "C2"="#F2C94C", "C3"="#7DBE6C")) +
  labs(x=NULL, y="Mean CNV value", title="Deletion") +
  theme_minimal(base_size=11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1, size=9),
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(face="bold", hjust=0.5)
  )

## ---- AMP plot
p_amp <- ggplot(plot_amp, aes(x=CNV_peak, y=mean, fill=Cluster)) +
  geom_bar(stat="identity", position=position_dodge(0.8), width=0.7) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                position=position_dodge(0.8), width=0.25, linewidth=0.4) +
  scale_fill_manual(values=c("C1"="#8EB3D6", "C2"="#F2C94C", "C3"="#7DBE6C")) +
  labs(x=NULL, y="Mean CNV value", title="Amplification") +
  theme_minimal(base_size=11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1, size=9),
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(face="bold", hjust=0.5)
  )

library(patchwork)
p_combined <- p_amp / p_del + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "top")

print(p_combined)
