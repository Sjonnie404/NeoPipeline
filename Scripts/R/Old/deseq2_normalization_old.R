# Rscript that will be run remotely for Deseq2 normalization
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")


library('DESeq2')
library("BiocParallel")
# register(MulticoreParam(4))
library('dplyr')
library('tidyverse')
library('tibble')
# library('pheatmap')
library('RColorBrewer')


# print('###### Succesfully loaded DESeq2 & Dplyr\n')
# print('Running code....\n')
# cat(getwd())



###### NOTE: Need to cite ashr & IHW

path = '/home/shane/Documents/pycharm_project_519/Output/Counts/gdc_download_20220422_083506.202427/'
matrix = read.csv2(paste(path,'readsMatrix.csv', sep=''), sep=',')

sample = read.csv2(paste(path,'gdc_sample_sheet.2022-05-09.tsv', sep=''), sep='\t')

new_matrix = matrix[,-1]
row.names(new_matrix) = matrix[,1]

coldata <- data.frame(matrix(ncol=1,nrow=473, dimnames=list(NULL, "condition")))
row.names(coldata) <- colnames(new_matrix)
coldata$type <- sample$Sample.Type
coldata$type <- factor(coldata$type)
# coldata$canditate <- 1:nrow(coldata)
# coldata$canditate <- factor(coldata$canditate)


dds <- DESeqDataSetFromMatrix(countData = new_matrix,
                              colData = coldata,
                              design = ~ type)

keep <- rowSums(counts(dds)) >= 10 # Not needed but speeds up memory
dds <- dds[keep,]

dds <- DESeq(dds)
plotDispEsts(dds)
res <- results(dds)
# res

# --------------------- example -----------------------------
# keep <- rowSums(counts(dds)) >= 10
# dds <- dds[keep, ]
#
# # varianceStabilizingTransformation
# varst <- vst(dds)
# ## retrieve variance transformed data from varst
# vst_mat <- assay(varst)
# #colnames(vst_mat) <- count_features$Sample
# vst_mat <- t(vst_mat)
# --------------------- example -----------------------------


resAsh <- lfcShrink(dds, coef=2, type="normal")
#resAsh
#plotMA(res, ylim=c(-2,2))
#plotMA(resAsh, ylim=c(-2,2))

res <- resAsh

res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()


rld <- vst(dds, blind=TRUE)
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

#pheatmap(rld_cor)


# Set thresholds
padj_cutoff <- 0.01

# Subset the significant results
sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
  dplyr::arrange(padj)

# Check significant genes output
# sig_res
normalized_counts <- counts(dds,
                            normalized = TRUE)


## Order results by padj values
top20_sig_genes <- sig_res %>%
  dplyr::arrange(padj) %>%
  dplyr::pull(gene) %>%
  head(n=20)

write.csv(top20_sig_genes, file=paste(path,'Deseq_significant_genes.csv'), row.names=F)

cat('Finished determining most significant genes!')
quit()

top20_sig_norm <- data.frame(normalized_counts) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% top20_sig_genes)

gathered_top20_sig <- top20_sig_norm %>%
  gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm))], key = "samplename", value = "normalized_counts")

gathered_top20_sig <- inner_join(ei[, c("sample_id", "group_id" )], gathered_top20_sig, by = c("sample_id" = "samplename"))

## plot using ggplot2
ggplot(gathered_top20_sig) +
  geom_point(aes(x = gene,
                 y = normalized_counts),
             position=position_jitter(w=0.1,h=0)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))


# Extract normalized counts for only the significant genes
sig_norm <- data.frame(normalized_counts) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% sig_res$gene)

# Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")

# Run pheatmap using the metadata data frame for the annotation
pheatmap(sig_norm[ , 2:length(colnames(sig_norm))],
         color = heat_colors,
         cluster_rows = T,
         show_rownames = F,
         border_color = NA,
         fontsize = 10,
         scale = "row",
         fontsize_row = 10,
         height = 20)

## Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 1.5 in either direction
# res_table_thres <- res_tbl %>%
  # mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 0.58)

res_table_thres <- res_tbl %>%
  mutate(threshold = padj < 0.01 & log2FoldChange >= 0.58)


## Volcano plot
ggplot(res_table_thres) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  ggtitle("Volcano plot") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  scale_y_continuous(limits = c(0,10)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))



#resOrdered <- res[order(res$pvalue),]
#resOrdered

#summary(res)

#res05 <- results(dds, alpha=0.05)
#summary(res05)

# (unevaluated code chunk)
#library("IHW")
#resIHW <- results(dds, filterFun=ihw)
#summary(resIHW)
#sum(resIHW$padj < 0.1, na.rm=TRUE)
#metadata(resIHW)$ihwResult



