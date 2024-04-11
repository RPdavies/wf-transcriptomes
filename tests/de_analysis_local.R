rm(list = ls())
# coldata <- read.csv("../../nextflow/sample_sheet.csv",
#                     row.names="alias", sep=",", stringsAsFactors=TRUE)
coldata <- read.csv("../../nextflow/sample_sheet_testPaired.csv",
                    row.names="alias", sep=",", stringsAsFactors=TRUE)
coldata$sample_id <- rownames(coldata)
#coldata$condition <- factor(coldata$condition, levels=rev(levels(coldata$condition)))

# check if control condition exists, sets as reference 
if(!"control" %in% coldata$condition)
  stop("sample_sheet.csv does not contain 'control' 
       condition - unable to set reference")
coldata$condition <- relevel(coldata$condition, ref = "control")

gene_cts <- read.delim("../../nextflow/de_analysis/all_gene_counts.tsv")

# Not including an intercept leads to broken results! Likely because 
# explicit contrasts matrix is then expected

#design <- model.matrix(~ 0 + ID + condition, data = coldata)
#design <- model.matrix(~ 0 + condition, data = coldata)

# check if ID is specified in sample_sheet, if so include in model
# to get 'paired t-test'
if("ID" %in% colnames(coldata))
  design <- model.matrix(~ ID + condition, data = coldata) else
    design <- model.matrix(~ condition, data = coldata)

library("edgeR")

y <- DGEList(gene_cts)
y <- calcNormFactors(y)
y <- estimateDisp(y,design)
fit <- glmQLFit(y,design)
head(fit$coefficients)
qlf <- glmQLFTest(fit)
edger_res <- topTags(qlf, n=nrow(y), sort.by="PValue")[[1]]

#write.table(as.data.frame(edger_res), file="results_dge_correct.tsv", sep="\t")

# plotMD(qlf)
# abline(h=c(-1,1), col="blue")
# plotQLDisp(fit)
# plotBCV(y)

#### make plots to test LFC direction 
i <- 1
png(paste0("test", Sys.time(), ".png"), width = 800)
par(mfrow = c(1,2))
boxplot(counts ~ condition, 
        data = data.frame(counts = y$counts[rownames(edger_res)[i],],
                          condition = coldata$condition),
        main = rownames(edger_res)[i])

plotMD(qlf)
for(j in 1:10)
points(edger_res$logCPM[i], edger_res$logFC[i],
       col = 3)
text(edger_res$logCPM[i], edger_res$logFC[i],
     labels = rownames(edger_res)[i],col="black",cex=1,pos=3)
dev.off()
