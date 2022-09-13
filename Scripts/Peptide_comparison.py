########################################################################################
#  Script for all MHCpan output functions: Analytics & comparisions.
#  Made by Shane Pullens, Utrecht University - Theoretical bioinformatics.
#  Version 1.0
########################################################################################



# imports
from pathlib import Path
import collections
import pandas as pd


def main():
    project_dir = Path.cwd()
    new_file_name = 'Old/skin_TCGA-SKCM_20220822_132048.170180'
    mhcpan_output = 'MHCpan_output.csv'
    fasta_file_name = 'skin_TCGA-SKCM_20220822_132048.170180_cryptic_sequences_translated.fasta'
    peptide_df = pd.read_csv(Path(project_dir / 'Output' / 'Counts' / new_file_name / mhcpan_output))
    fasta_file = open(Path(project_dir / 'Output' / 'Counts' / new_file_name / fasta_file_name), 'r')
    fasta = "".join(fasta_file.readlines())

    # with pd.option_context('display.max_rows', 100, 'display.max_columns', 20):
    #     print(peptide_df)
    # exit()
    peptide_amount_prediciton(peptide_df, fasta, ratio_percentage=1)
    print('got here')
    exit()
    return None


def peptide_amount_prediciton(peptide_df, fasta, ratio_percentage = 2, kmer=9):
    """
    This function compares the amount of expected binders (amount of unique 9-mers, sliding windows * ratio) compared to
    the total amount of found peptides in the read.
    :param peptides: MHCpan output with predicted binders
    :param fasta: the fasta file that contains the reads
    :param ratio: the user defined ratio for expected binders in a pool of random 9-mer peptides.
    :return:
    """
    counts_df = peptide_df['Identity'].value_counts().to_frame()
    counts_df.index.name = 'Transcripts'
    counts_df.reset_index(inplace=True)
    # counts_df['Transcripts'] = counts_df.index.values
    # print(counts_df)
    # exit()

    #predicted_count = dict(zip(counts.index.values, list(counts)))

    # with pd.option_context('display.max_rows', 100, 'display.max_columns', 20):
    #     print(peptide_df)
    fasta_df = pd.DataFrame(columns=['Transcripts','Fasta_observed'])
    for i, read in enumerate(fasta.replace('>','$$>').split('$$')):
        header, seq = read.split('\n', 1)

        if header == '':  # Skips whitelines
            continue
        header = header.split('.')[0].replace('>','')
        totalkmers = len(seq.replace('\n','')) - kmer + 1
        estKmers = totalkmers*(ratio_percentage/100)
        fasta_df.loc[i] = [header, estKmers]
        # fasta_count[header] = estKmers
    # print(fasta_df)

    print('MHCpan predicted:')
    print(counts_df)
    print('Observed:')
    print(fasta_df)

    # full_df = counts_df.join(fasta_df)
    full_df = pd.merge(counts_df, fasta_df, how='left', on='Transcripts')


    #TODO Also add the
    print(full_df)

    exit()
    return None


def peptide_comparison(file1, file2):
    """
    This function will compare the found cryptic peptides to the canonical peptides, any overlap will result in removal
    of the cryptic peptide from the list, since it can also be a canonical peptide.
    :param file1:
    :param file2:
    :return:
    """

    return None

if __name__ == "__main__":
    main()