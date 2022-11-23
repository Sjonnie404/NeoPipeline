# Rscript that will be run remotely for Deseq2 normalization
cat('Loading packages...\n')
shhh <- supressPackageStartupMessages 
shhh(library('DESeq2'))
shhh(library('dplyr'))
shhh(library('tibble'))
shhh(library('tidyverse'))

cat('Reading data files...\n')
#in_path = 'breast_project_testing/CMI/readsMatrix.csv' # Debug
path = args[1]
in_path = paste(path,'/readsMatrix.csv', sep="")

matrix = read.csv2(in_path, sep=',')
stripped_matrix = matrix[,-1] # sets gene ID to rowname
row.names(stripped_matrix) = matrix[,1]
dummy_coldata_values <- 1:ncol(stripped_matrix)
dummy_coldata_names <- colnames(stripped_matrix)
dummy_coldata <- data.frame(dummy_coldata_names, dummy_coldata_values)

#dummy_coldata <- data.frame(stripped_matrix[,1])
#dummy_coldata$numbering <- 1:nrow(dummy_coldata)

dds <- DESeqDataSetFromMatrix(countData = stripped_matrix,
                              colData = dummy_coldata,
                              design = ~1)

# Won't have effect on outcome but speeds up memory
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

cat('Performing normalization calculations...\n')
# if number of patients is below 15?, use normal transformations,
# the mathematical shortcut can't be used because initial prediction
# will be flawed.
if (ncol(stripped_matrix) <= 15){
  cat('Low number of patients detected!\n')
  cat('Mathematical shortcuts could be flawed, using normal transformation instead!\n')
  vst_res = varianceStabilizingTransformation(dds)
} else {
  vst_res = vst(dds)
}
vst_res = assay(vst_res)
vst_res = data.frame(vst_res)

vst_res$row_median <- apply(vst_res[,-1], 1, median)
vst_res_sig = data.frame(filter_at(vst_res, vars(row_median),
                                ~. > quantile(., probs = 0.95)))

out_path = paste(path, '/significant_genes.csv', sep="")

vst_res_sig$genes <- rownames(vst_res_sig)
write.table(vst_res_sig$genes, file=out_path, sep=',', row.names=FALSE, col.names=FALSE)
cat('Sucessfully wrote significant genes\n')

