# 같은 subclone에서 함께 나타나는 gene pairs
subclonal_data <- df_with_gene[status == "Subclonal" & !is.na(gene)]

# 같은 sample, 같은 cluster
subclonal_cooccur <- merge(subclonal_data[, .(sample, cluster_id, gene, clord)],
                           subclonal_data[, .(sample, cluster_id, gene, clord)],
                           by = c("sample", "cluster_id", "clord"),
                           allow.cartesian = TRUE)

subclonal_cooccur <- subclonal_cooccur[gene.x < gene.y]  # 중복 제거

# 빈도
cooccur_freq <- subclonal_cooccur[, .N, by = .(gene.x, gene.y, clord)][
  order(-N)][1:20]

print(cooccur_freq)

# Network 시각화
library(igraph)
library(ggraph)

edges <- cooccur_freq[clord == "C2" & N >= 3]  # C2에서 3번 이상
g <- graph_from_data_frame(edges[, .(gene.x, gene.y, N)])

ggraph(g, layout = 'fr') +
  geom_edge_link(aes(width = N), alpha = 0.5) +
  geom_node_point(size = 5, color = "steelblue") +
  geom_node_text(aes(label = name), repel = TRUE) +
  labs(title = "C2: Co-occurring Subclonal Mutations") +
  theme_graph()

library(data.table)

# PyClone 결과 파일들이 있는 디렉토리
result_dir <- "pyclone_results/"

# 모든 결과 파일 찾기
result_files <- list.files(result_dir, 
                          pattern = "_results.tsv$",  # 또는 적절한 패턴
                          full.names = TRUE)

cat("Found", length(result_files), "result files\n")

# 모든 파일 읽어서 합치기
all_results <- rbindlist(lapply(result_files, fread), fill = TRUE)

######

library(data.table)
library(ggplot2)

# 1. 기본 확인
table(all_results$clord)

# 2. Clonal/Subclonal 분류
all_results[, status := ifelse(cellular_prevalence > 0.9, "Clonal", "Subclonal")]

# 3. 샘플별 요약 (그룹 정보 포함)
sample_summary <- all_results[, .(
  clord = unique(clord),
  n_mutations = .N,
  n_clonal = sum(status == "Clonal"),
  n_subclonal = sum(status == "Subclonal"),
  clonal_fraction = mean(status == "Clonal"),
  n_clusters = length(unique(cluster_id)),
  mean_ccf = mean(cellular_prevalence),
  median_ccf = median(cellular_prevalence)
), by = sample]

head(sample_summary)

# 4. 그룹별 비교
group_comparison <- sample_summary[, .(
  n_samples = .N,
  mean_clonal_frac = mean(clonal_fraction),
  sd_clonal_frac = sd(clonal_fraction),
  mean_n_clusters = mean(n_clusters),
  sd_n_clusters = sd(n_clusters),
  mean_ccf = mean(mean_ccf)
), by = clord]

print(group_comparison)
library(rstatix)
kruskal.test(clonal_fraction ~ clord, data = sample_summary)

# Pairwise 비교 (3그룹이니까)
pairwise.wilcox.test(sample_summary$clonal_fraction, 
                     sample_summary$clord, 
                     p.adjust.method = "BH")


library(ggplot2)

# 1. Boxplot + p-values
ggplot(sample_summary, aes(x = clord, y = clonal_fraction, fill = clord)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
  labs(title = "Clonal Mutation Fraction by Group",
       subtitle = "Kruskal-Wallis p = 0.259 (NS)",
       x = "Group", y = "Clonal Fraction") +
  theme_bw() +
  theme(legend.position = "none")

# 2. 그룹별 요약 통계
sample_summary[, .(
  N = .N,
  Mean = mean(clonal_fraction),
  SD = sd(clonal_fraction),
  Median = median(clonal_fraction),
  Q1 = quantile(clonal_fraction, 0.25),
  Q3 = quantile(clonal_fraction, 0.75)
), by = clord]

# 3. Violin plot
ggplot(sample_summary, aes(x = clord, y = clonal_fraction, fill = clord)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, 
               fill = "red", color = "red") +
  labs(title = "Clonal Fraction Distribution",
       x = "Group", y = "Clonal Fraction") +
  theme_bw()


# 1. Cluster 개수 비교
kruskal.test(n_clusters ~ clord, data = sample_summary)

ggplot(sample_summary, aes(x = clord, y = n_clusters, fill = clord)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
  labs(title = "Number of Clusters by Group",
       x = "Group", y = "Number of Clusters") +
  theme_bw()

# 2. Subclonal mutations 개수
kruskal.test(n_subclonal ~ clord, data = sample_summary)

# 3. CCF 분포 (모든 mutation)
kruskal.test(cellular_prevalence ~ clord, data = all_results)

ggplot(all_results, aes(x = clord, y = cellular_prevalence, fill = clord)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  labs(title = "CCF Distribution by Group",
       x = "Group", y = "Cellular Prevalence") +
  theme_bw()


# C2와 C3의 cluster 구조 차이?
cluster_complexity <- all_results[, .(
  n_clusters = length(unique(cluster_id)),
  shannon = diversity(cellular_prevalence, index = "shannon"),
  clord = unique(clord)
), by = sample]

# 그룹별 비교
kruskal.test(n_clusters ~ clord, data = cluster_complexity)
kruskal.test(shannon ~ clord, data = cluster_complexity)

# 시각화
library(gridExtra)
p1 <- ggplot(cluster_complexity, aes(x = clord, y = n_clusters, fill = clord)) +
  geom_boxplot() + theme_bw()

p2 <- ggplot(cluster_complexity, aes(x = clord, y = shannon, fill = clord)) +
  geom_boxplot() + theme_bw()

grid.arrange(p1, p2, ncol = 2)
