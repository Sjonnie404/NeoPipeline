# Rscript that will be run remotely for Deseq2 normalization
library('DESeq2')

library("BiocParallel")
# register(MulticoreParam(4))
library('dplyr')
library('tibble')
library('tidyverse')
library('pheatmap')
library('RColorBrewer')


# print('###### Succesfully loaded DESeq2 & Dplyr\n')
# print('Running code....\n')
# cat(getwd())



###### NOTE: Need to cite ashr & IHW

path = 'readsMatrix.csv'
matrix = read.csv2(path, sep=',')

sample = read.csv2('gdc_sample_sheet.2022-05-09_edited_again.tsv', sep='\t')

new_matrix = matrix[,-1] # sets gene ID to rowname
row.names(new_matrix) = matrix[,1]


coldata <- data.frame(sample$Sample.Type)
row.names(coldata) <- colnames(new_matrix)
coldata$patient <- colnames(new_matrix)
coldata$patient <- factor(coldata$patient)

coldata$type <- sample$Sample.Type
coldata$type <- factor(coldata$type)
coldata = coldata[,-1]

# coldata$canditate <- 1:nrow(coldata)
# coldata$canditate <- factor(coldata$canditate)


dds <- DESeqDataSetFromMatrix(countData = new_matrix,
                              colData = coldata, 
                              design = ~1) # TODO Check out the formula This cant be the case.

keep <- rowSums(counts(dds)) >= 10 # Not needed but speeds up memory
dds <- dds[keep,]


tmp = vst(dds)
tmp = assay(tmp)
tmp = data.frame(tmp)

tmp$row_median <- apply(tmp[,-1], 1, median)
tmp$row_sum <- apply(tmp[,-1], 1, sum)


# tmp %>%
#   ggplot(aes(x='Genes', y=row_median, colour='')) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_point(position = position_jitterdodge(), alpha = 0.15) +
#   theme_bw()+scale_fill_grey(start = 0.8, end = 1) +
#   labs(title=paste("All genes, #datapoints:",nrow(tmp)), x="", y="score",fill="")

# tmp2 = data.frame(filter_at(tmp, vars(row_median), ~. > quantile(., probs = 0.75)))
# tmp2 %>%
#   ggplot(aes(x='Genes', y=row_median, colour='')) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_point(position = position_jitterdodge(), alpha = 0.2) +
#   theme_bw()+scale_fill_grey(start = 0.8, end = 1) +
#   labs(title=paste("3rd quartile, #datapoints:",nrow(tmp2)), x="", y="score",fill="")


# tmp3 = data.frame(filter_at(tmp, vars(row_median), ~. > quantile(., probs = 0.95)))
# tmp3 %>%
#   ggplot(aes(x='Genes', y=row_median, colour='')) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_point(position = position_jitterdodge(), alpha = 0.4) +
#   theme_bw()+scale_fill_grey(start = 0.8, end = 1) +
#   labs(title=paste("top 5 %, #datapoints:",nrow(tmp3)), x="", y="score",fill="")


tmp4 = data.frame(filter_at(tmp, vars(row_median), ~. > quantile(., probs = 0.99)))
# tmp4 %>%
#   ggplot(aes(x='Genes', y=row_median, colour='')) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_point(position = position_jitterdodge(), alpha = 0.6) +
#   theme_bw()+scale_fill_grey(start = 0.8, end = 1) +
#   labs(title=paste("top 1 %, #datapoints:",nrow(tmp4)), x="", y="score",fill="")


tmp4$genes <- rownames(tmp4) 


write.csv(tmp4$genes, file='Deseq_significant_genes_test.csv', row.names=F)
cat('Sucessfully wrote significant genes')
stop()


#----------------------
dds <- DESeq(dds)

plotDispEsts(dds)

#idx <- identify(res$baseMean, res$log2FoldChange)
# rowsnames(res[idx])

res <- results(dds)
res

resultsNames(dds)

# NOTE: This is only shrinkage, so no normalization, this is only for better visualization
resApe <- lfcShrink(dds, coef='type_Primary.Tumor_vs_Metastatic', type="apeglm")
#resAsh
plotMA(res, ylim=c(-2,2))
plotMA(resApe, ylim=c(-2,2))

# res <- resAsh

res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()

res_tble_pos_fold <- dplyr::filter(res_tbl, log2FoldChange > 0) %>%
                drop_na()

# sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) #%>%



# res_ash_tbl <- resAsh %>%
#       data.frame() %>%
#       rownames_to_column(var="gene") %>%
#       as_tibble()


rld <- vst(dds, blind=TRUE)
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

pheatmap(rld_cor)


# Set thresholds
padj_cutoff <- 0.05

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

write.csv(top20_sig_genes, file='Deseq_significant_genes.csv', row.names=F)

exit()

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

ggplot(res_tble_pos_fold) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj))) +
  ggtitle("Positive fold change") +
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



