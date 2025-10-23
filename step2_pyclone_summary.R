
library(dplyr)
patient_files <- list.files("./pyclone_results", pattern="\\.tsv$", full.names=TRUE)

pyclone_summary <- data.frame()

for(patient in patient_files){
  tryCatch({
    
    # 파일 읽기
    results <- read.table(patient, sep="\t", header=TRUE, check.names=FALSE)
    
    # Cluster summary
    cluster_summary <- results %>%
      group_by(cluster_id) %>%
      summarise(
        median_prev = median(cellular_prevalence),
        n_mutations = n()
      ) %>%
      arrange(desc(median_prev))
    
    # Main clone
    main_clone_prev <- cluster_summary$median_prev[1]
    main_clone_size <- cluster_summary$n_mutations[1]
    
    # Subclones
    n_clusters <- nrow(cluster_summary)
    n_subclones <- n_clusters - 1
    
    if(n_clusters > 1) {
      subclone_mutations <- sum(cluster_summary$n_mutations[2:n_clusters])
      subclone_fraction <- subclone_mutations / nrow(results)
    } else {
      subclone_mutations <- 0
      subclone_fraction <- 0
    }
    
    # Shannon diversity
    props <- cluster_summary$n_mutations / sum(cluster_summary$n_mutations)
    shannon <- -sum(props * log(props + 1e-10))
    
    # Patient ID 추출
    patient_id <- gsub("\\.tsv", "", basename(patient))
    
    # Summary
    patient_summary <- data.frame(
      patient_id = patient_id,
      n_total_mutations = nrow(results),
      n_clusters = n_clusters,
      n_subclones = n_subclones,
      main_clone_prev = main_clone_prev,
      subclone_fraction = subclone_fraction,
      shannon_diversity = shannon
    )
    
    pyclone_summary <- rbind(pyclone_summary, patient_summary)
    
    cat("✓", patient_id, "- clusters:", n_clusters, "\n")
    
  }, error = function(e) {
    cat("✗ Error for", basename(patient), ":", e$message, "\n")
  })
}

#cc는 clustering결과 파일 
names(pyclone_summary)[1] <- "sample"

pyclone_summary <- left_join(pyclone_summary,cc,by="sample")

################################################################
# Cluster간 비교
################################################################


library(dplyr)
library(ggplot2)
library(ggpubr)


summary_by_cluster <- pyclone_summary %>%
  group_by(clord) %>%
  summarise(
    n = n(),
    mean_clusters = mean(n_clusters),
    mean_subclones = mean(n_subclones),
    mean_subclone_frac = mean(subclone_fraction),
    mean_shannon = mean(shannon_diversity)
  )

wilcox_subclone <- wilcox.test(
  subclone_fraction ~ (clord == "C2"),
  data = pyclone_summary
)
cat("\nSubclone fraction (C2 vs others):\n")
cat("p-value:", wilcox_subclone$p.value, "\n")

# Shannon diversity
wilcox_shannon <- wilcox.test(
  shannon_diversity ~ (clord == "C2"),
  data = pyclone_summary
)

p1 <- ggplot(pyclone_summary, aes(x=clord, y=subclone_fraction, fill=clord)) +
  geom_boxplot(alpha=0.7, outlier.shape=NA) +
  geom_jitter(width=0.2, alpha=0.4, size=2) +
  stat_compare_means(comparisons = list(c("C1","C2"), c("C2","C3"))) +
  theme_bw() +
  labs(title="Subclonal Fraction by Cluster",
       x="Cluster", y="Subclonal Mutation Fraction") +
  theme(legend.position="none")

# Shannon diversity
p2 <- ggplot(pyclone_summary, aes(x=clord, y=shannon_diversity, fill=clord)) +
  geom_boxplot(alpha=0.7, outlier.shape=NA) +
  geom_jitter(width=0.2, alpha=0.4, size=2) +
  stat_compare_means(comparisons = list(c("C1","C2"), c("C2","C3"))) +
  theme_bw() +
  labs(title="Shannon Diversity by Cluster",
       x="Cluster", y="Shannon Diversity Index") +
  theme(legend.position="none")

# Number of subclones
p3 <- ggplot(pyclone_summary, aes(x=clord, y=n_subclones, fill=clord)) +
  geom_boxplot(alpha=0.7, outlier.shape=NA) +
  geom_jitter(width=0.2, alpha=0.4, size=2) +
  stat_compare_means(comparisons = list(c("C1","C2"), c("C2","C3"))) +
  theme_bw() +
  labs(title="Number of Subclones by Cluster",
       x="Cluster", y="Number of Subclones") +
  theme(legend.position="none")
library(cowplot)
plot_grid(p1, p2, p3, ncol=2)
