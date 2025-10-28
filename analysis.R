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


