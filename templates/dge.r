#!/usr/bin/env Rscript
# Created by Michal Bukowski (michal.bukowski@tuta.io) under GPL-3.0 license

# Runs DESeq2 DGE analysis for a given samlon quantmerge output and specified params

library("DESeq2")

# Create output directory for the output TSV file
dir.create(dirname("${outFile}"), recursive=TRUE)

# Read salmon quantmerge output (NumReads for replicas), infer group names from
# column names (group_n, where n is a replica number), set the reference group
# as the first group on the list groups passed to the script.
countData <- read.csv("${file}", header=TRUE, row.names=1, sep="\t")
countData <- round(countData, 0)
repNames  <- colnames(countData)
suffixPos <- unlist(gregexpr("_[0-9]+\$", repNames))
groups    <- substr(repNames, 1, suffixPos-1)
groups    <- factor(groups)
groups    <- relevel(groups, ref="${groups[0]}")
# Based on column names and group names inffered from them, create meta data
# assigning seach column name (row.names) to a group (a column)
metaData  <- data.frame(row.names=repNames, groups)

# Use countData and metaData to run DGE analysis. Filter out rows with the sum
# of counts < 10.
dds <- DESeqDataSetFromMatrix(countData=countData, colData=metaData,
                              design=~groups)
dds <- dds[ rowSums( counts(dds) ) >= 10, ]
dds <- DESeq(dds)

# Process the results using the given alpha threshold, perform lfc shrinkage.
res  <- results(dds, alpha=${alpha})
coef <- "groups_${groups[1]}_vs_${groups[0]}"
res  <- lfcShrink(dds, res=res, coef=coef, type="apeglm")

# Add to final results rows skipped before, fill them up with NA values.
skipped       <- setdiff(rownames(countData), rownames(res))
res[skipped,] <- NA

# Rearrange the data structure in order to have row names (locus_tags) as a first
# column. Save the data to a TSV file.
res_df <- cbind(rownames(res), data.frame(res, row.names=NULL))
colnames(res_df)[1] <- "locus_tag"
write.table(res_df, file="${outFile}", sep="\t", row.names=FALSE, quote=FALSE)
