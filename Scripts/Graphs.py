########################################################################################
#  Script to analyse & plot data
#  Made by Shane Pullens, Utrecht University - Theoretical bioinformatics.
#  Version 1.0
########################################################################################
# This function looks at all possible skin genes.
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

# Note to self: This really should have been done with a notebook or different functions......

project_dir = Path.cwd()

#Import libraries
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt
# %matplotlib inline

# venn2_unweighted(subsets = (215, 9441, 6), set_labels = ('Candidate peptides', 'Known peptide DB'))
# plt.show()

import seaborn as sns
#
# Make the Venn diagrams & pie charts.
#define data

#####
# Breast
# "BTTcryp|HLA-1:\t 39435"
# "BTTcryp|HLA-2:\t 62919"
# Skin
# "STTcryp|HLA-1:\t 31204"
# "STTcryp|HLA-2:\t 51205"
####
data = [31204, 51205]
labels = ['HLA-A*01:01\n31204', 'HLA-A*02:01\n51205']

#define Seaborn color palette to use
colors = sns.color_palette('pastel')[0:2]

# create pie chart
plt.pie(data, labels = labels, colors = colors, autopct='%.0f%%')
plt.show()
exit()

## This translates back from the found transcripts to the gene type of the gene.
# mode = 'Canonical'
mode = 'Cryptic'
project_dir = Path.cwd()
new_file_name = 'final_skin_rna_genes'
true_file_name = 'skin_TCGA-SKCM_20221124_144444.636999'
#
TCGA_RNA = true_file_name
TCGA = new_file_name
# # new_file_name = 'final_skin_rna_genes'
# # true_file_name = 'skin_TCGA-SKCM_20221029_151238.792423'
#
if mode == 'Cryptic':
    mhcpan_output = 'MHCpan_output_cryptic.csv'
    fasta_file_name = true_file_name+'_cryptic_sequences_translated.fasta'
else:
    mhcpan_output = 'MHCpan_output_canonical.csv'
    fasta_file_name = true_file_name+'_canonical_sequences_translated.fasta'
#
fasta_file = open(Path(project_dir / 'Output' / 'Counts' / new_file_name / fasta_file_name), 'r')
fasta = "".join(fasta_file.readlines())


peptide_df = pd.read_csv(Path(project_dir / 'Output' / 'Counts' / new_file_name / mhcpan_output))
peptide_df = pd.read_csv(Path(project_dir / 'Output' / 'Counts' / 'RNA_ID.txt'))

# peptide_df = peptide_df[peptide_df['MHC'] == 'HLA-A*02:01']
trans_gene_dict = {}
# translist = set(peptide_df['Identity'])
translist = peptide_df['ID']
# print(translist)
# exit()
# print(translist.shape)



for i, read in enumerate(fasta.replace('>','$$>').split('$$')):
    if read == '':
        continue
    header, seq = read.split('\n', 1)

    if header == '':  # Skips whitelines
        continue
    header = header.replace('>','').split('|')
    trans_gene_dict.update({header[0].split('.')[0]: header[2]})

# print(trans_gene_dict)

# translist = translist.str.replace('.','-')
translist = translist.str.split('.').str[0]
print(translist)
print('-'*80)
print(trans_gene_dict)
# exit()

# genelist = [trans_gene_dict[gene] for gene in translist]

# print(genelist)
# print(len(genelist))
# print(len(set(genelist)))
#
# exit()
# HCMI_RNA = 'BOX_RNA_skin_HCMI-CMDC_20220725_092631.124914'
# HCMI_ALL = 'BOX_ALL_skin_HCMI-CMDC_20220725_092440.108408'
#
# ER_RNA = 'BOX_RNA_skin_EXCEPTIONAL_RESPONDERS-ER_20220725_095103.246184'
# ER_ALL = 'BOX_ALL_skin_EXCEPTIONAL_RESPONDERS-ER_20220725_094942.365046'


## Calculates the amount of genes there are per gene_type
df = pd.read_csv(Path(project_dir / 'Output' / 'Counts' / TCGA / 'Raw_counts' / 'd93b1580-c20d-4611-ac99-c0cb33f37a45.rna_seq.augmented_star_gene_counts.tsv'), sep='\t', skiprows=6,header=None)
genelist = translist

df[[0, 99]] = df[0].astype(str).str.split('.', 1, expand=True)
# with pd.option_context('display.max_rows', 10, 'display.max_columns', None):  # more options can be specified also
#     print(df)
# exit()


# print(df.shape)
df = df[df[0].isin(genelist)]
# with pd.option_context('display.max_rows', 10, 'display.max_columns', None):  # more options can be specified also
#     print(df)
# print(df.shape)
# exit()
df_counts = df[2].value_counts()
df_counts = pd.DataFrame({'Gene_type': df_counts.index, 'Amount': df_counts.values})

print(df_counts)
exit()
# df = df[['Unnamed: 2']]
# df_counts = df_counts[df_counts['Amount'] > 200]
df_counts = df_counts.replace({'processed_pseudogene':'proc_pseudogene',
                               'unprocessed_pseudogene':'unproc_pseudogene',
                               'transcribed_unprocessed_pseudogene':'trans_unproc_pseudogene',
                               'transcribed_processed_pseudogene':'trans_proc_pseudogene'})

# with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
# #     print(df)
# #     print('-'*80)
#     print(df_counts)
# # exit()
# BLUE, ORANGE = sns.color_palette()[:2]
# sns.set_theme()
#

## This makes the bar plot of amount of genes per genetype
# new_df = pd.DataFrame({'Gene_type': ['lncRNA', 'trans_unproc_pseudogene', 'unproc_pseudogene', 'lncRNA', 'trans_unproc_pseudogene', 'unproc_pseudogene'],
#                        'Amount': [194, 46, 9, 195, 44, 9],
#                        'HLA': ['HLA-A*01:01', 'HLA-A*01:01', 'HLA-A*01:01', 'HLA-A*02:01', 'HLA-A*02:01', 'HLA-A*02:01']})
# ax = sns.barplot(data = new_df, x='Amount', y='Gene_type', hue='HLA')
# sns.move_legend(ax, 'lower right')
# # plt.pie(df_counts['Amount'], labels=df_counts['Gene_type'], autopct='%.0f%%')
#
# for i in ax.containers:
#     ax.bar_label(i,)
#
# ax.set(title='Number of genes per gene-type', ylabel=None, xlabel=None)
# plt.tight_layout()
# plt.show()

# HCMI Data
# HCMI_matrix_RNA = pd.read_csv(Path(project_dir / 'Output' / 'Counts' / 'Boxplot_data' / HCMI_RNA / 'readsMatrix.csv'))
# HCMI_matrix_RNA['median'] = HCMI_matrix_RNA.iloc[:, 1:].median(axis=1)  # Calculate median for each row
# This selects bottom 95%, however, we want the top 5%
# HCMI_limiter_RNA = HCMI_matrix_RNA['median'].quantile(0.95)   # selects 95%
# HCMI_mask_RNA = HCMI_matrix_RNA[HCMI_matrix_RNA['median'] <= HCMI_limiter_RNA].index
# HCMI_matrix_RNA_top5 = HCMI_matrix_RNA.drop(HCMI_matrix_RNA.index[HCMI_mask_RNA])
# HCMI_median_RNA = HCMI_matrix_RNA_top5[['Gene_id', 'median']]
# HCMI_median_RNA['dataset'] = 'HCMI'
# HCMI_median_RNA['type'] = 'RNA'
#
# HCMI_matrix_ALL = pd.read_csv(Path(project_dir / 'Output' / 'Counts' / 'Boxplot_data' / HCMI_ALL / 'readsMatrix.csv'))
# HCMI_matrix_ALL['median'] = HCMI_matrix_ALL.iloc[:, 1:].median(axis=1)  # Calculate median for each row
# # This selects bottom 95%, however, we want the top 5%
# HCMI_limiter_ALL = HCMI_matrix_ALL['median'].quantile(0.95)   # selects 95%
# HCMI_mask_ALL = HCMI_matrix_ALL[HCMI_matrix_ALL['median'] <= HCMI_limiter_ALL].index
# HCMI_matrix_ALL_top5 = HCMI_matrix_ALL.drop(HCMI_matrix_ALL.index[HCMI_mask_ALL])
# HCMI_median_ALL = HCMI_matrix_ALL_top5[['Gene_id', 'median']]
# HCMI_median_ALL['dataset'] = 'HCMI'
# HCMI_median_ALL['type'] = 'ALL'
#
#
# # ER Data
# # ER_matrix_RNA = pd.read_csv(Path(project_dir / 'Output' / 'Counts' / 'Boxplot_data' / ER_RNA / 'readsMatrix.csv'))
# # ER_matrix_RNA['median'] = ER_matrix_RNA.iloc[:, 1:].median(axis=1)  # Calculate median for each row
# # # This selects bottom 95%, however, we want the top 5%
# # ER_limiter_RNA = ER_matrix_RNA['median'].quantile(0.95)   # selects 95%
# # ER_mask_RNA = ER_matrix_RNA[ER_matrix_RNA['median'] <= ER_limiter_RNA].index
# # ER_matrix_RNA_top5 = ER_matrix_RNA.drop(ER_matrix_RNA.index[ER_mask_RNA])
# # ER_median_RNA = ER_matrix_RNA_top5[['Gene_id', 'median']]
# # ER_median_RNA['dataset'] = 'ER'
# # ER_median_RNA['type'] = 'RNA'
# #
# # ER_matrix_ALL = pd.read_csv(Path(project_dir / 'Output' / 'Counts'/ 'Boxplot_data' / ER_ALL / 'readsMatrix.csv'))
# # ER_matrix_ALL['median'] = ER_matrix_ALL.iloc[:, 1:].median(axis=1)  # Calculate median for each row
# # # This selects bottom 95%, however, we want the top 5%
# # ER_limiter_ALL = ER_matrix_ALL['median'].quantile(0.95)   # selects 95%
# # ER_mask_ALL = ER_matrix_ALL[ER_matrix_ALL['median'] <= ER_limiter_ALL].index
# # ER_matrix_ALL_top5 = ER_matrix_ALL.drop(ER_matrix_ALL.index[ER_mask_ALL])
# # ER_median_ALL = ER_matrix_ALL_top5[['Gene_id', 'median']]
# # ER_median_ALL['dataset'] = 'ER'
# # ER_median_ALL['type'] = 'ALL'
#
# # TCGA Data
# TCGA_matrix_RNA = pd.read_csv(Path(project_dir / 'Output' / 'Counts' / 'Boxplot_data' / TCGA_RNA / 'readsMatrix.csv'))
# TCGA_matrix_RNA['median'] = TCGA_matrix_RNA.iloc[:, 1:].median(axis=1)  # Calculate median for each row
# # This selects bottom 95%, however, we want the top 5%
# TCGA_limiter_RNA = TCGA_matrix_RNA['median'].quantile(0.95)   # selects 95%
# TCGA_mask_RNA = TCGA_matrix_RNA[TCGA_matrix_RNA['median'] <= TCGA_limiter_RNA].index
# TCGA_matrix_RNA_top5 = TCGA_matrix_RNA.drop(TCGA_matrix_RNA.index[TCGA_mask_RNA])
# TCGA_median_RNA = TCGA_matrix_RNA_top5[['Gene_id', 'median']]
# TCGA_median_RNA['dataset'] = 'TCGA'
# TCGA_median_RNA['type'] = 'RNA'
#
# TCGA_matrix_ALL = pd.read_csv(Path(project_dir / 'Output' / 'Counts' / 'Boxplot_data' / TCGA_ALL / 'readsMatrix.csv'))
# TCGA_matrix_ALL['median'] = TCGA_matrix_ALL.iloc[:, 1:].median(axis=1)  # Calculate median for each row
# # This selects bottom 95%, however, we want the top 5%
# TCGA_limiter_ALL = TCGA_matrix_ALL['median'].quantile(0.95)   # selects 95%
# TCGA_mask_ALL = TCGA_matrix_ALL[TCGA_matrix_ALL['median'] <= TCGA_limiter_ALL].index
# TCGA_matrix_ALL_top5 = TCGA_matrix_ALL.drop(TCGA_matrix_ALL.index[TCGA_mask_ALL])
# TCGA_median_ALL = TCGA_matrix_ALL_top5[['Gene_id', 'median']]
# TCGA_median_ALL['dataset'] = 'TCGA'
# TCGA_median_ALL['type'] = 'ALL'
#
#
# combined_median = pd.concat([HCMI_median_ALL,HCMI_median_RNA,ER_median_ALL,ER_median_RNA,TCGA_median_ALL,TCGA_median_RNA])
# print(combined_median)
# print(combined_median.shape)
#
#
# # combined_median['HCMI_median_all'] = HCMI_median_ALL['median']
# # combined_median['HCMI_median_RNA'] = HCMI_median_RNA['median']
#
# # combined_median['ER_median_all'] = ER_median_ALL['median']
# # combined_median['ER_median_RNA'] = ER_median_RNA['median']
# #
# # combined_median['TCGA_median_all'] = TCGA_median_ALL['median']
# # combined_median['TCGA_median_RNA'] = TCGA_median_RNA['median']
#
#
# # print(combined_median.shape)
# # print(combined_median)
# # exit()
#
# # TODO: Need to add all median series together in a df for the boxplots
# # TODO: Do believe we don't need the indexes.
#
#
# # with pd.option_context('display.max_rows',1000, 'display.max_columns',1000):
# #     print(matrix_RNA)
#
# # Testing
# combined_median = combined_median[combined_median['median'] <= 1000]
# combined_median = combined_median[combined_median['type'] == 'RNA']
# BLUE, ORANGE = sns.color_palette()[:2]
#
# ax = sns.boxplot(x='dataset', y='median', data=combined_median, color=ORANGE)
# ax.set(title='Medians of top 5% expressed genes in skin cancer.\nCutoff = 1.000, RNA only')
# plt.tight_layout()
# plt.show()
