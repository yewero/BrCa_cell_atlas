
suppressMessages({
  library("GSEABase")
})

load("xCell.data.rda")

signatures <- xCell.data$signatures
gs_LS <- geneIds(signatures)

gene_LS <- lapply(names(gs_LS), function(x) {
  cellType <- gsub("\\.txt", "", gsub("%", "_", x))
  y <- data.frame(cellType = cellType, gene = gs_LS[[x]], stringsAsFactors = F)
  return(y)
})
gene_DF <- do.call("rbind", gene_LS)
cat("[info] Cell type number:", length(unique(gene_DF$cellType)), "\n")
cat("[info] Total gene number:", length(gene_DF$gene), "\n")
cat("[info] Unique gene number:", length(unique(gene_DF$gene)), "\n")

write.table(x = gene_DF, file = "signatures_xCell.txt", row.names = F, col.names = T, quote = F, sep = "\t")
