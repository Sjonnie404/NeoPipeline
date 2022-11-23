# Rscript that will be run remotely for Deseq2 normalization
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

cat('Loading packages...\n')
shhh <- suppressPackageStartupMessages # It's a library, so shhh!

shhh(library('DESeq2'))
shhh(library('dplyr'))
shhh(library('tibble'))
shhh(library('tidyverse'))

cat('Reading data files...\n')
path = args[1]
in_path = paste(path,'/readsMatrix.csv', sep="")

matrix = read.csv2(in_path, sep=',')
stripped_matrix = matrix[,-1] # sets gene ID to row name
row.names(stripped_matrix) = matrix[,1]

# Deseq expects a coldata variable to generate a dds object.
dummy_coldata_values <- 1:ncol(stripped_matrix)
dummy_coldata_names <- colnames(stripped_matrix)
dummy_coldata <- data.frame(dummy_coldata_names, dummy_coldata_values)

# dds object normaly does all the normalization, however, since we don't have control goroup data, this could not be
# performed. Because of this we mainly use the dds object to exploit is for it VST function.
dds <- DESeqDataSetFromMatrix(countData = stripped_matrix,
                              colData = dummy_coldata,
                              design = ~1)

# Won't have effect on outcome but speeds up memory
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

cat('Performing normalization calculations...\n')
# if number of patients is below 15?, use normal transformations, the mathematical shortcut can't be used because
# initial prediction will be flawed.
if (ncol(stripped_matrix) <= 15){
  cat('Low number of patients detected!\n')
  cat('Mathematical shortcuts could be flawed, using normal transformation instead!\n')
  vst_res = varianceStabilizingTransformation(dds)
} else {
  vst_res = vst(dds)
}
# Generates new data object with transformed counts.
vst_res = assay(vst_res)
vst_res = data.frame(vst_res)

# calculate row median
vst_res$row_median <- apply(vst_res[,-1], 1, median)
# Saves the top 5% of highest scoring genes.
vst_res_sig = data.frame(filter_at(vst_res, vars(row_median),
                                ~. > quantile(., probs = 0.95)))

out_path = paste(path, '/significant_genes.csv', sep="")
# Write only the gene names without the counts.
vst_res_sig$genes <- rownames(vst_res_sig)
write.table(vst_res_sig$genes, file=out_path, sep=',', row.names=FALSE, col.names=FALSE)
cat('Sucessfully wrote significant genes\n')

