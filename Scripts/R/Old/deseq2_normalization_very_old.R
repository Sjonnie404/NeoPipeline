# Rscript that will be run remotely for Deseq2 normalization
library('DESeq2')
# library('dplyr')

cat('###### Succesfully loaded DESeq2 & Dplyr\n')
cat('Running code....\n')
# cat(getwd())

path = '/home/shane/Documents/pycharm_project_519/Output/Counts/gdc_download_20220422_083506.202427/readsMatrix.csv'
matrix = data.frame(read.csv(path, sep='\t'))

# head(matrix, 2)
cat('\n---\n')
# matrix %>% remove_rownames %>% column_to_rownames(var='Gene')
new_matrix = matrix[,-1]
cat(dim(new_matrix))
cat('\n#######################################\n')

# new_matrix
cat('\n')

# YOU ARE HERE, SOMETHING is going really wrong.

matrix[, 1]
cat('\n\n######\n\n')
rownames(new_matrix) = matrix[,1]

# new_matrix <- data.frame(matrix[,-1], row.names=matrix[,1])
# head(new_matrix, 10)
# samp.with.rownames <- data.frame(samp[,-1], row.names=samp[,1])
# matrix = read.csv2('/Output/Counts/gdc_download_20220422_083506.202427/readsMatrix.csv') # TODO fix dynamic path

# testing
# airway <- get(load('/home/shane/Documents/pycharm_project_519/Output/Counts/gdc_download_20220422_083506.202427/airway.RData'))
# str(airway)
# head(matrix,10)
# se <- airway
# se
# ddsHT <- DESeqDataSet(matrix, design = ~ cell + dex)
#
# dds <- DESeq(ddsHT)
# res <- results(dds)
#
# res



