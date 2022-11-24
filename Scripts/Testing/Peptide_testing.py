########################################################################################
#  Script for all MHCpan output functions: Analytics & comparisions.
#  Made by Shane Pullens, Utrecht University - Theoretical bioinformatics.
#  Version 1.0
########################################################################################



# imports
from pathlib import Path
import collections
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from scipy import interpolate
import seaborn as sns
from scipy import stats


def main():
    mode = 'Cryptic'
    # mode = 'Canonical'
    project_dir = Path.cwd()
    new_file_name = 'debug_test_remove_me_toStop'
    true_file_name = 'skin_TCGA-SKCM_20221124_091344.401089'

    # new_file_name = 'final_skin_rna_genes'
    # true_file_name = 'skin_TCGA-SKCM_20221029_151238.792423'

    if mode == 'Cryptic':
        mhcpan_output = 'MHCpan_output_cryptic.csv'
        fasta_file_name = true_file_name+'_cryptic_sequences_translated.fasta'
    else:
        mhcpan_output = 'MHCpan_output_canonical.csv'
        fasta_file_name = true_file_name+'_canonical_sequences_translated.fasta'

    fasta_file = open(Path(project_dir / 'Output' / 'Counts' / new_file_name / fasta_file_name), 'r')
    fasta = "".join(fasta_file.readlines())
    peptide_df = pd.read_csv(Path(project_dir / 'Output' / 'Counts' / new_file_name / mhcpan_output))
    # peptide_df = peptide_df[peptide_df['%Rank_EL'] <= 0.1]


    full_df_HLA1, full_df_HLA2 = peptide_amount_prediciton(peptide_df, fasta, ratio_percentage= 1.0)

    # full_df.rename({'Fasta_observed': 'Fasta_expected'}, axis=1, inplace=True)
    full_df_HLA1.rename({'Fasta_observed': 'Fasta_expected'}, axis=1, inplace=True)

    # full_df['Ratio'] = full_df['MHCpan_observed'].div(full_df['Fasta_expected'])  # divides
    full_df_HLA1['Ratio'] = full_df_HLA1['MHCpan_observed'].div(full_df_HLA1['Fasta_expected'])  # divides

    # new_df = full_df[full_df['Fasta_expected'] <= 250].fillna(0).replace(np.inf, 0)
        # new_df = full_df
        # work_df = full_df.drop(['%Rank_EL'], axis=1)
    # work_df = new_df['Ratio'].to_frame('Ratio')

    # df_mva = work_df.rolling(100).mean()  # moving average with a window size of 30

    # print(work_df)
    # print('-'*80)
    # print(df_mva)
    # exit()

    title = new_file_name.split('_', 1)[1]
    # title = 'Bacterial Validation'
    sns.set_theme()
    # ax = sns.lineplot(data=new_df['Ratio'], linewidth=1, alpha=0.4).set(title=f'Ratio Observed / Expected\n{title} - {mode}')
    # ax1 = sns.lineplot(data=df_mva, linewidth=1, color='red')
    # plt.tight_layout()
    # plt.show()
    #
    # exit()
    ### rank_df = new_df[new_df['%Rank_BA'] <= 0.1]
    # full_df_HLA1 = full_df_HLA1[full_df_HLA1['%Rank_EL'] <= 0.1]
    # print(new_df)
    # print(rank_df)

    # print(new_df['MHCpan_observed'].mean())
    # exit()

    #### slope, intercept, r_value, p_value, std_err = stats.linregress(rank_df['Fasta_expected'], rank_df['MHCpan_observed'])
    # With new_df
    # slope, intercept, r_value, p_value, std_err = stats.linregress(new_df['Fasta_expected'], new_df['MHCpan_observed'])
    # ax2 = sns.scatterplot(data=new_df, x='Fasta_expected', y='MHCpan_observed',
    #                 linewidth=0, alpha=1, hue='%Rank_EL').set(title=f'Observed v.s. expected peptides\n{title} - {mode}')
    # ax3 = sns.regplot(data=new_df, x='Fasta_expected', y='MHCpan_observed', scatter=False)

    # Full df (HLA1)
    slope, intercept, r_value, p_value, std_err = stats.linregress(full_df_HLA1['Fasta_expected'], full_df_HLA1['MHCpan_observed'])

    ax2 = sns.scatterplot(data=full_df_HLA1, x='Fasta_expected', y='MHCpan_observed',
                      linewidth=0, alpha=1, hue='%Rank_EL').set(title=f'Observed v.s. expected peptides\n{title} - {mode}')

    ax3 = sns.regplot(data=full_df_HLA1, x='Fasta_expected', y='MHCpan_observed', scatter=False)

    ax3.text(0.5, 0.9, f"y={slope:.1f}x + {intercept:.1f}", horizontalalignment='right',
             verticalalignment='bottom', transform=ax3.transAxes)
    plt.tight_layout()
    plt.show()

    return None


def peptide_amount_prediciton(peptide_df, fasta, ratio_percentage = 2.0, kmer=9):
    """
    This function compares the amount of expected binders (amount of unique 9-mers, sliding windows * ratio) compared to
    the total amount of found peptides in the read.
    :param peptides: MHCpan output with predicted binders
    :param fasta: the fasta file that contains the reads
    :param ratio: the user defined ratio for expected binders in a pool of random 9-mer peptides.
    :return:
    """
    mhcpan_df_HLA1 = peptide_df[peptide_df['MHC'] == 'HLA-A*01:01']
    mhcpan_df_HLA2 = peptide_df[peptide_df['MHC'] == 'HLA-A*02:01']

    mhcpan_df_HLA1 = mhcpan_df_HLA1.sort_values('%Rank_EL', ascending=False).drop_duplicates(['Peptide'], keep='first')
    mhcpan_df_HLA2 = mhcpan_df_HLA2.sort_values('%Rank_EL', ascending=False).drop_duplicates(['Peptide'], keep='first')

    mhcpan_df_HLA1 = mhcpan_df_HLA1['Identity'].value_counts().to_frame()
    mhcpan_df_HLA1.index.name = 'Transcripts'
    mhcpan_df_HLA1.reset_index(inplace=True)
    mhcpan_df_HLA1.rename({'Identity': 'MHCpan_observed'}, axis=1, inplace=True)

    mhcpan_df_HLA2 = mhcpan_df_HLA2['Identity'].value_counts().to_frame()
    mhcpan_df_HLA2.index.name = 'Transcripts'
    mhcpan_df_HLA2.reset_index(inplace=True)
    mhcpan_df_HLA2.rename({'Identity': 'MHCpan_observed'}, axis=1, inplace=True)

    fasta_df = pd.DataFrame(columns=['Transcripts','Fasta_observed'])
    for i, read in enumerate(fasta.replace('>','$$>').split('$$')):
        if read == '':
            continue
        header, seq = read.split('\n', 1)

        if header == '':  # Skips whitelines
            continue
        header = header.split('.')[0].replace('>','')
        totalkmers = len(seq.replace('\n','')) - kmer + 1

        estKmers = round(totalkmers*(ratio_percentage/100))
        fasta_df.loc[i] = [header, estKmers]

    fasta_df = fasta_df.groupby('Transcripts').sum('Fasta_observerd').reset_index()
    full_df_HLA1 = mhcpan_df_HLA1.merge(fasta_df, how='outer', on='Transcripts').fillna(0)
    full_df_HLA2 = mhcpan_df_HLA2.merge(fasta_df, how='outer', on='Transcripts').fillna(0)
    # with pd.option_context('display.max_rows', 20, 'display.max_columns', None):  # more options can be specified also
    #     print(full_df)

    EL_df = peptide_df[['Identity', '%Rank_EL', '%Rank_BA']].rename({'Identity': 'Transcripts'}, axis=1)
    full_df_HLA1 = full_df_HLA1.merge(EL_df, how='outer', on='Transcripts')
    full_df_HLA2 = full_df_HLA2.merge(EL_df, how='outer', on='Transcripts')

    return full_df_HLA1, full_df_HLA2


def duplicate_remover(df):
    """
    Method to remove all duplicates, and from each peptide group, only keep the peptide with the lowest %EL rank.
    However, all duplicates seem to hold the same EL_rank (which makes sense due to MHCpan). Meaning this is now a
    oneline, kept this function here for future purpose perhaps.
    :param df:
    :return:
    """
    duplicated_df = df[df.duplicated(['Peptide'], keep=False)]

    df_no_dups = df.merge(duplicated_df, how='left', indicator=True)
    df_no_dups = df_no_dups[df_no_dups['_merge'] == 'left_only']

    duplicated_df = duplicated_df.sort_values('%Rank_EL', ascending=False).drop_duplicates(['Peptide'], keep='first')
    df = pd.concat([df_no_dups, duplicated_df]).sort_values('%Rank_EL')

    return df

if __name__ == "__main__":
    main()