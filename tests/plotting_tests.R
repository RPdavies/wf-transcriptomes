plotMD(qlf)
abline(h=c(-1,1), col="blue")
# circle and label the most signifficant gene for direction verification
for(j in 1:10)
  points(edger_res$logCPM[1], edger_res$logFC[1],
         col = 3)
text(edger_res$logCPM[1], edger_res$logFC[1],
     labels = rownames(edger_res)[1], pos=3)
boxplot(counts ~ condition, 
        data = data.frame(counts = y$counts[rownames(edger_res)[1],],
                          condition = coldata$condition),
        main = rownames(edger_res)[1])