########################################################################################
#  Script to analyse & plot data
#  Made by Shane Pullens, Utrecht University - Theoretical bioinformatics.
#  Version 1.0
########################################################################################

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

project_dir = Path.cwd()
HCMI_RNA = 'BOX_RNA_skin_HCMI-CMDC_20220725_092631.124914'
HCMI_ALL = 'BOX_ALL_skin_HCMI-CMDC_20220725_092440.108408'

ER_RNA = 'BOX_RNA_skin_EXCEPTIONAL_RESPONDERS-ER_20220725_095103.246184'
ER_ALL = 'BOX_ALL_skin_EXCEPTIONAL_RESPONDERS-ER_20220725_094942.365046'

TCGA_RNA = 'BOX_RNA_skin_TCGA-SKCM_20220725_092926.911548'
TCGA_ALL = 'BOX_ALL_skin_TCGA-SKCM_20220725_094019.550812'

# HCMI Data
HCMI_matrix_RNA = pd.read_csv(Path(project_dir / 'Output' / 'Counts' / 'Boxplot_data' / HCMI_RNA / 'readsMatrix.csv'))
HCMI_matrix_RNA['median'] = HCMI_matrix_RNA.iloc[:, 1:].median(axis=1)  # Calculate median for each row
# This selects bottom 95%, however, we want the top 5%
HCMI_limiter_RNA = HCMI_matrix_RNA['median'].quantile(0.95)   # selects 95%
HCMI_mask_RNA = HCMI_matrix_RNA[HCMI_matrix_RNA['median'] <= HCMI_limiter_RNA].index
HCMI_matrix_RNA_top5 = HCMI_matrix_RNA.drop(HCMI_matrix_RNA.index[HCMI_mask_RNA])
HCMI_median_RNA = HCMI_matrix_RNA_top5[['Gene_id', 'median']]
HCMI_median_RNA['dataset'] = 'HCMI'
HCMI_median_RNA['type'] = 'RNA'

HCMI_matrix_ALL = pd.read_csv(Path(project_dir / 'Output' / 'Counts' / 'Boxplot_data' / HCMI_ALL / 'readsMatrix.csv'))
HCMI_matrix_ALL['median'] = HCMI_matrix_ALL.iloc[:, 1:].median(axis=1)  # Calculate median for each row
# This selects bottom 95%, however, we want the top 5%
HCMI_limiter_ALL = HCMI_matrix_ALL['median'].quantile(0.95)   # selects 95%
HCMI_mask_ALL = HCMI_matrix_ALL[HCMI_matrix_ALL['median'] <= HCMI_limiter_ALL].index
HCMI_matrix_ALL_top5 = HCMI_matrix_ALL.drop(HCMI_matrix_ALL.index[HCMI_mask_ALL])
HCMI_median_ALL = HCMI_matrix_ALL_top5[['Gene_id', 'median']]
HCMI_median_ALL['dataset'] = 'HCMI'
HCMI_median_ALL['type'] = 'ALL'


# ER Data
ER_matrix_RNA = pd.read_csv(Path(project_dir / 'Output' / 'Counts' / 'Boxplot_data' / ER_RNA / 'readsMatrix.csv'))
ER_matrix_RNA['median'] = ER_matrix_RNA.iloc[:, 1:].median(axis=1)  # Calculate median for each row
# This selects bottom 95%, however, we want the top 5%
ER_limiter_RNA = ER_matrix_RNA['median'].quantile(0.95)   # selects 95%
ER_mask_RNA = ER_matrix_RNA[ER_matrix_RNA['median'] <= ER_limiter_RNA].index
ER_matrix_RNA_top5 = ER_matrix_RNA.drop(ER_matrix_RNA.index[ER_mask_RNA])
ER_median_RNA = ER_matrix_RNA_top5[['Gene_id', 'median']]
ER_median_RNA['dataset'] = 'ER'
ER_median_RNA['type'] = 'RNA'

ER_matrix_ALL = pd.read_csv(Path(project_dir / 'Output' / 'Counts'/ 'Boxplot_data' / ER_ALL / 'readsMatrix.csv'))
ER_matrix_ALL['median'] = ER_matrix_ALL.iloc[:, 1:].median(axis=1)  # Calculate median for each row
# This selects bottom 95%, however, we want the top 5%
ER_limiter_ALL = ER_matrix_ALL['median'].quantile(0.95)   # selects 95%
ER_mask_ALL = ER_matrix_ALL[ER_matrix_ALL['median'] <= ER_limiter_ALL].index
ER_matrix_ALL_top5 = ER_matrix_ALL.drop(ER_matrix_ALL.index[ER_mask_ALL])
ER_median_ALL = ER_matrix_ALL_top5[['Gene_id', 'median']]
ER_median_ALL['dataset'] = 'ER'
ER_median_ALL['type'] = 'ALL'

# TCGA Data
TCGA_matrix_RNA = pd.read_csv(Path(project_dir / 'Output' / 'Counts' / 'Boxplot_data' / TCGA_RNA / 'readsMatrix.csv'))
TCGA_matrix_RNA['median'] = TCGA_matrix_RNA.iloc[:, 1:].median(axis=1)  # Calculate median for each row
# This selects bottom 95%, however, we want the top 5%
TCGA_limiter_RNA = TCGA_matrix_RNA['median'].quantile(0.95)   # selects 95%
TCGA_mask_RNA = TCGA_matrix_RNA[TCGA_matrix_RNA['median'] <= TCGA_limiter_RNA].index
TCGA_matrix_RNA_top5 = TCGA_matrix_RNA.drop(TCGA_matrix_RNA.index[TCGA_mask_RNA])
TCGA_median_RNA = TCGA_matrix_RNA_top5[['Gene_id', 'median']]
TCGA_median_RNA['dataset'] = 'TCGA'
TCGA_median_RNA['type'] = 'RNA'

TCGA_matrix_ALL = pd.read_csv(Path(project_dir / 'Output' / 'Counts' / 'Boxplot_data' / TCGA_ALL / 'readsMatrix.csv'))
TCGA_matrix_ALL['median'] = TCGA_matrix_ALL.iloc[:, 1:].median(axis=1)  # Calculate median for each row
# This selects bottom 95%, however, we want the top 5%
TCGA_limiter_ALL = TCGA_matrix_ALL['median'].quantile(0.95)   # selects 95%
TCGA_mask_ALL = TCGA_matrix_ALL[TCGA_matrix_ALL['median'] <= TCGA_limiter_ALL].index
TCGA_matrix_ALL_top5 = TCGA_matrix_ALL.drop(TCGA_matrix_ALL.index[TCGA_mask_ALL])
TCGA_median_ALL = TCGA_matrix_ALL_top5[['Gene_id', 'median']]
TCGA_median_ALL['dataset'] = 'TCGA'
TCGA_median_ALL['type'] = 'ALL'


combined_median = pd.concat([HCMI_median_ALL,HCMI_median_RNA,ER_median_ALL,ER_median_RNA,TCGA_median_ALL,TCGA_median_RNA])
print(combined_median)
print(combined_median.shape)


# combined_median['HCMI_median_all'] = HCMI_median_ALL['median']
# combined_median['HCMI_median_RNA'] = HCMI_median_RNA['median']

# combined_median['ER_median_all'] = ER_median_ALL['median']
# combined_median['ER_median_RNA'] = ER_median_RNA['median']
#
# combined_median['TCGA_median_all'] = TCGA_median_ALL['median']
# combined_median['TCGA_median_RNA'] = TCGA_median_RNA['median']


# print(combined_median.shape)
# print(combined_median)
# exit()

# TODO: Need to add all median series together in a df for the boxplots
# TODO: Do believe we don't need the indexes.


# with pd.option_context('display.max_rows',1000, 'display.max_columns',1000):
#     print(matrix_RNA)

# Testing
combined_median = combined_median[combined_median['median'] <= 1000]
combined_median = combined_median[combined_median['type'] == 'RNA']
BLUE, ORANGE = sns.color_palette()[:2]

ax = sns.boxplot(x='dataset', y='median', data=combined_median, color=ORANGE)
ax.set(title='Medians of top 5% expressed genes in skin cancer.\nCutoff = 1.000, RNA only')
plt.tight_layout()
plt.show()
