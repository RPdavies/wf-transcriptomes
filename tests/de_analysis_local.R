coldata <- read.csv("../../nextflow/sample_sheet_testPaired.csv",
                    stringsAsFactors = T)
coldata$condition <- relevel(coldata$condition, ref = "control")

gene_cts <- read.delim("../../nextflow/de_analysis/all_gene_counts.tsv")

# Not including an intercept leads to broken results! 
#design <- model.matrix(~ 0 + ID + condition, data = coldata)
#design <- model.matrix(~ 0 + condition, data = coldata)

#design <- model.matrix(~ ID + condition, data = coldata)
design <- model.matrix(~ condition, data = coldata)

library("edgeR")

y <- DGEList(gene_cts)
y <- calcNormFactors(y)
y <- estimateDisp(y,design)
fit <- glmQLFit(y,design)
head(fit$coefficients)
qlf <- glmQLFTest(fit)
edger_res <- topTags(qlf, n=nrow(y), sort.by="PValue")[[1]]

head(edger_res)

length(which(edger_res$PValue < 0.05 & edger_res$logFC > 0))
length(which(edger_res$PValue < 0.05 & edger_res$logFC < 0))

plotMD(qlf)
abline(h=c(-1,1), col="blue")
plotQLDisp(fit)
plotBCV(y)
