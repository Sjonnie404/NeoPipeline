########################################################################################
#  Script to fetch gene sequences based on coordinates
#  Input: Count file that contains Ensemble IDs
#  Output: Fasta file with headers, containing Ensemble information and NCBI sequence
#  Made by Shane Pullens, Utrecht University - Theoretical bioinformatics.
########################################################################################

# imports
import pprint
from pathlib import Path
import pandas as pd
import seaborn as sns
import subprocess

#  Predefined variables
pp = pprint.PrettyPrinter(indent=4)  # Used to clearly print json files
server = "https://rest.ensembl.org"
verbose = False
project_dir = Path.cwd()
sns.set()  # Sets seaborn as default style


def main():
    return None


def runDeseq(matrix_path):
    """
    Runs the DESeq function in a R subprocess to determine best genes from the counts data.
    R script an be found at ~/Scripts/R/deseq2_normalization.R
    :param matrix: Matrix that contains counts for each gene of each patient.
    :return: None, saves a list of significant genes instead named 'significant_genes.csv'
    """
    print('Please note: Deseq will break when used with 3 or less patients!')
    print('~'*80)
    print('Starting R script...')

    command = '/home/shane/miniconda3/envs/neoPipe/bin/Rscript'
    args = [str(matrix_path)]
    path2script = '/home/shane/Documents/pycharm_project_519/Scripts/R/deseq2_normalization.R'
    retcode = subprocess.run([command, path2script] + args, universal_newlines=True)
    print(retcode)

    print('Finishing R script...')
    print('~'*80)
    return None


def getReadsMatrix(paths, RNAonly=True):
    """
    Here, all count data from each patient (each patient has a separate file) is merged together in a single matrix.
    :param paths: list of all patient filepaths
    :param RNAonly: Boolean to select all genes, or only RNA-related genes.
    :return: returns the combined reads matrix.
    """
    counter = 0
    grouped_df = pd.DataFrame()
    file_path = ''

    for file_path in paths:
        counter += 1
        df = pd.read_csv(file_path, sep='\t', skiprows=6, header=None)
        df.columns = ['gene_id', 'gene_name', 'gene_type', 'unstranded', 'stranded_first', 'stranded_second',
                           'tpm_unstranded', 'fpkm_unstranded', 'fpkm_uq_unstranded']

        if RNAonly:
            search_criteria = ['RNA', 'unprocessed_pseudogene']
            df = df[df['gene_type'].str.contains('|'.join(search_criteria))]

        df = df[df['gene_id'].str.contains('_PAR_Y') == False]  # remove all _PAR_Y genes due to gene name conflicts
        df = df.iloc[:, [0, 3]]
        df.columns = ['Gene_id', 'Patient_{i}'.format(i=counter)]

        if counter <= 1:
            grouped_df = df
        else:
            grouped_df = pd.concat([grouped_df, df['Patient_{i}'.format(i=counter)]], axis=1)

    print('Trying to save matrix...')
    try:
        grouped_df.to_csv(Path(file_path.parent.parent / 'readsMatrix.csv'), index=False)
    except:
        print(f'Something went wrong trying to save {file_path.parent.parent} / readsMatrix.csv.')
    print('successfully wrote reads matrix')
    return grouped_df


if __name__ == "__main__":
    main()