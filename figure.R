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
#CNV Amp/Del Frequency 그룹별 
#======================================================
library(dplyr)
library(tidyr)

amp_thr <- 1
del_thr <- -1
top_n_amp <- 18
top_n_del <- 18

cnv_bin_amp <- cnv_cyto >= amp_thr
cnv_bin_del <- cnv_cyto <= del_thr

amp_freq_all <- colMeans(cnv_bin_amp, na.rm=TRUE)
del_freq_all <- colMeans(cnv_bin_del, na.rm=TRUE)

top_amp  <- names(sort(amp_freq_all, decreasing=TRUE))[1:min(top_n_amp, length(amp_freq_all))]
top_del  <- names(sort(del_freq_all, decreasing=TRUE))[1:min(top_n_del, length(del_freq_all))]
top_band <- unique(c(top_amp, top_del))

# cytoband order: (amp+del) 전체 빈도(최대값) 기준 정렬
ord_score <- pmax(amp_freq_all[top_band], del_freq_all[top_band], na.rm=TRUE)
band_levels <- names(sort(ord_score, decreasing=TRUE))

df_amp <- as.data.frame(cnv_bin_amp[, top_band, drop=FALSE]) %>%
  mutate(sample = rownames(.), cluster = as.character(split_fac[sample])) %>%
  pivot_longer(cols = all_of(top_band), names_to="cytoband", values_to="event") %>%
  group_by(cluster, cytoband) %>%
  summarise(freq = mean(event, na.rm=TRUE), .groups="drop") %>%
  mutate(type = "Amp")

df_del <- as.data.frame(cnv_bin_del[, top_band, drop=FALSE]) %>%
  mutate(sample = rownames(.), cluster = as.character(split_fac[sample])) %>%
  pivot_longer(cols = all_of(top_band), names_to="cytoband", values_to="event") %>%
  group_by(cluster, cytoband) %>%
  summarise(freq = mean(event, na.rm=TRUE), .groups="drop") %>%
  mutate(type = "Del")

plot_df <- bind_rows(df_amp, df_del) %>%
  mutate(
    cluster = factor(cluster, levels=c("C1","C2","C3")),
    type    = factor(type, levels=c("Amp","Del")),
    cytoband = factor(cytoband, levels=band_levels)
  )
library(ggplot2)

# 논문용(너무 쨍하지 않게) - 필요하면 바꿔도 됨
cluster_cols <- c(
  "C1" = "#D55E00",  # burnt orange
  "C2" = "#009E73",  # green
  "C3" = "#0072B2"   # blue
)

p <- ggplot(plot_df, aes(x=cytoband, y=freq, fill=cluster)) +
  geom_col(
    width = 0.78,
    position = position_dodge2(width = 0.82, padding = 0.12),
    color = NA
  ) +
  facet_grid(type ~ ., scales="free_y", switch="y") +
  scale_fill_manual(values = cluster_cols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
  labs(x=NULL, y="Frequency") +
  theme_classic(base_size = 11) +
  theme(
    # strip(패널 제목) 최소화
    strip.background = element_blank(),
    strip.placement  = "outside",
    strip.text.y.left = element_text(face="bold", size=11),
    strip.text.x = element_blank(),

    # x축 글자(겹침 방지)
    axis.text.x  = element_text(angle=45, hjust=1, vjust=1, size=9),
    axis.ticks.x = element_blank(),

    # 선/여백 정리
    axis.line.x = element_blank(),
    panel.spacing = unit(8, "mm"),

    # legend
    legend.title = element_blank(),
    legend.position = "right",
    legend.key.height = unit(4, "mm")
  )

p
