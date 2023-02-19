#!/usr/bin/env Rscript
library("DESeq2")

dir.create(dirname("${outFile}"), recursive=TRUE)

countData <- read.csv("${file}", header=TRUE, row.names=1, sep="\t")
countData <- round(countData, 0)
repNames  <- colnames(countData)
suffixPos <- unlist(gregexpr("_[0-9]+\$", repNames))
groups    <- substr(repNames, 1, suffixPos-1)
groups    <- factor(groups)
groups    <- relevel(groups, ref="${groups[0]}")
metaData  <- data.frame(row.names=repNames, groups)

dds <- DESeqDataSetFromMatrix(countData=countData, colData=metaData,
                              design=~groups)
dds <- dds[ rowSums( counts(dds) ) >= 10, ]
dds <- DESeq(dds)

res  <- results(dds, alpha=${alpha})
coef <- "groups_${groups[1]}_vs_${groups[0]}"
res  <- lfcShrink(dds, res=res, coef=coef, type="apeglm")

skipped       <- setdiff(rownames(countData), rownames(res))
res[skipped,] <- NA

res_df <- cbind(rownames(res), data.frame(res, row.names=NULL))
colnames(res_df)[1] <- "locus_tag"
write.table(res_df, file="${outFile}", sep="\t", row.names=FALSE, quote=FALSE)
