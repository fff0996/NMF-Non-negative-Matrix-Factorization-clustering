library(ComplexHeatmap); library(circlize); library(matrixStats); library(RColorBrewer) ; library(IntNMF)

gsva <- read.table("Hallmark GSVA위치")
cnv <- read.table("gistic 결과 위치")
library(maftools)
library(data.table)
library(dplyr)
library(tidyr)
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

# 희귀 노이즈 제거: 너무 희귀한 유전자는 드롭 (예: ≥3 샘플에서 변이)
min_samples <- 3
keep_genes <- colSums(mut_mat) >= min_samples
mut_mat <- mut_mat[, keep_genes, drop = FALSE]




mut <- read.table("somatic mutation calling 결과 위치")
# 샘플 (행) * features (열) 형태로 맞춰주기 
# common 공통 샘플로 통일
common <- Reduce(intersect,list(rownames(gsva), rownames(cnv),rownames(mut)))
#나중에 히트맵 그릴때 이 데이터로 그림 
gsva_raw <- gsva[common,]
cnv_raw <- cnv[common,]
mut_raw <- mut[common,]

dat1 <- gsva_raw
dat2 <- cnv_raw
dat3 <- mut_raw

## Make all data positive by shifting to positive direction.
## Also rescale the datasets so that they are comparable.
if (!all(dat1>=0)) dat1 <- pmax(dat1 + abs(min(dat1)), .Machine$double.eps)
dat1 <- dat1/max(dat1)
if (!all(dat2>=0)) dat2 <- pmax(dat2 + abs(min(dat2)), .Machine$double.eps)
dat2 <- dat2/max(dat2)
#dat3 (mutation) 는 binary
#if (!all(dat3>=0)) dat3 <- pmax(dat3 + abs(min(dat3)), .Machine$double.eps)
#dat3 <- dat3/max(dat3)
# 빈도 필터
M <- dat3
p <- colMeans(M==1)
keep <- which(p >= 0.05 & p <= 0.60)
M1 <- M[, keep, drop=FALSE]
dat3 <- M1

# The function nmf.mnnals requires the samples to be on rows and variables on columns.

dat <- list(dat1,dat2,dat3)
fit <- nmf.mnnals(dat=dat,k=4,maxiter=200,st.count=20,n.ini=15,ini.nndsvd=TRUE,seed=TRUE)
ClusterEntropy(ComputedClusters=fit$clusters, TrueClasses=true.cluster.assignment$cluster.id)
SilhouettePlot(fit, cluster.col = NULL)



#결과 바탕으로 히트맵 그리기.
## 1) 공통 샘플 정렬
common <- Reduce(intersect, list(rownames(gsva_raw), rownames(cnv_raw), rownames(mut_raw), rownames(fit$W)))
stopifnot(length(common) >= 2)

## 2) 정렬
W     <- fit$W[common, , drop=FALSE]
cl    <- setNames(paste0("C", fit$clusters[common]), common)
sord  <- names(sort(apply(W, 1, max)))    # W 최대값 기준 정렬
clord <- cl[sord]; lev <- unique(clord)

## 3) 히트맵 입력 = 원본(raw)
G <- gsva_raw[sord, , drop=FALSE]          # GSVA raw
C <- cnv_raw[sord,  , drop=FALSE]          # CNV  raw (log2)
M <- mut_raw[sord,  , drop=FALSE]          # MUT  raw (0/1)
## 4) 표준화/클립
zcol <- function(X, cap=2){
  X <- as.matrix(X)
  Z <- scale(X)                 # column-wise
  Z[!is.finite(Z)] <- 0
  Z[Z >  cap] <-  cap
  Z[Z < -cap] <- -cap
  Z
}

Gz <- zcol(G, cap=2) 
Ccl  <- pmax(pmin(C, 2), -2)               # CNV는 원본을 [-2,2]로 클립
Mt   <- M                                  # Mutation은 0/1 그대로

## 5) 색
col_gsva <- circlize::colorRamp2(c(-2,0,2), c("#2b6cb0","white","#c53030"))
col_cnv  <- circlize::colorRamp2(c(-2,0,2), c("#225ea8","white","#e34a33"))
col_mut  <- c("0"="white","1"="black")
ha <- ComplexHeatmap::HeatmapAnnotation(Cluster=factor(clord, levels=lev))

## 6) 히트맵
hts <- list(
   gsva     = ComplexHeatmap::Heatmap(t(Gz),  name="GSVA",      col=col_gsva, top_annotation=ha,
                                     show_row_names=TRUE, show_column_names=FALSE,
                                     cluster_rows=TRUE,  cluster_columns=FALSE,
                                     column_split=factor(clord, levels=lev)),
  cnv      = ComplexHeatmap::Heatmap(t(Ccl), name="CopyNumber", col=col_cnv,
                                     show_row_names=TRUE, show_column_names=FALSE,
                                     cluster_rows=TRUE,  cluster_columns=FALSE,
                                     column_split=factor(clord, levels=lev)),
    mutation = ComplexHeatmap::Heatmap(t(Mt),  name="Mutation",  col=col_mut,
                                     show_row_names=TRUE, show_column_names=FALSE,
                                     cluster_rows=TRUE,  cluster_columns=FALSE,
                                     column_split=factor(clord, levels=lev))

)
ComplexHeatmap::draw(Reduce(`%v%`, hts), merge_legend=TRUE)


#Save
ht <- Reduce(`%v%`, hts)  # 이미 만든 히트맵 결합

## 크기 자동 산출(대충 가독 사이즈)
ns <- length(clord)                                  # 샘플 수(열)
nr <- sum(c(ncol(Mt), ncol(Ccl), ncol(Gz)), na.rm=TRUE)  # 총 피처 수(행)
w_in <- max(6,  min(0.12 * ns, 30))  # 가로(inch)
h_in <- max(8,  min(0.15 * nr, 60))  # 세로(inch)

## PDF
pdf("intNMF_heatmap.pdf", width=w_in, height=h_in, onefile=FALSE)
draw(ht, merge_legend=TRUE, heatmap_legend_side="right", annotation_legend_side="right")
dev.off()

## PNG(고해상도)
png("intNMF_heatmap.png", units="in", res=300, width=w_in, height=h_in, type="cairo")
draw(ht, merge_legend=TRUE, heatmap_legend_side="right", annotation_legend_side="right")
dev.off()
