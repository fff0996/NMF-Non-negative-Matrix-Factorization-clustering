#Make input format

#1. rna-seq (Gene set)
d <- read.table("Hallmark_GSVA_results.csv",sep=",",header=T,check.names=F)

rownames(d) <- d[,1]      # 첫 컬럼을 rownames로
d <- d[,-1]               # 첫 컬럼 제거 (실제 값만 남김)

# 이제 d는 행=gene set, 열=샘플
#head(d)[,1:4]

expr <- t(as.matrix(d))   # samples x gene_sets
#비음수 만들기 x축이동; min-shift방식
col_mins <- apply(expr, 2, min, na.rm=TRUE)
#각 열에서 해당 gene set의 최솟값을 빼줌
#예: 만약 HALLMARK_ADIPOGENESIS 열의 최소가 -0.5라면,
#그 열 전체에 -(-0.5) = +0.5 를 더한 효과
#→ 그 열의 최소값이 정확히 0이 됨
#혹시 값이 딱 0이 되면 수치계산에서 곤란할 수 있어서
#0 대신 아주 작은 양수(0.000001) 로 밀어줌
expr_pos <- sweep(expr, 2, col_mins, FUN="-") + 1e-6
