########################################################################################
#  Script to rank the peptides & compare them to known peptide files.
#  Made by Shane Pullens, Utrecht University - Theoretical bioinformatics.
#  Version 1.0
########################################################################################


import time
import pandas as pd
from pathlib import Path

project_dir = Path.cwd()


def main():
    mode = 'Cryptic'
    new_file_name = 'wednessday_full_run_skin_with_check'
    true_file_name = 'skin_TCGA-SKCM_20221012_093328.490034'

    if mode == 'Cryptic':
        mhcpan_output = 'MHCpan_output_cryptic.csv'
    else:
        mode = 'Canonical'
        mhcpan_output = 'MHCpan_output_canonical.csv'
    peptide_df = pd.read_csv(Path(project_dir / 'Output' / 'Counts' / new_file_name / mhcpan_output))

    hla_a01_df = peptide_df[peptide_df['MHC'] == 'HLA-A*01:01']
    hla_a02_df = peptide_df[peptide_df['MHC'] == 'HLA-A*02:01']

    # with pd.option_context('display.max_rows', 10, 'display.max_columns', None):  # more options can be specified also
    #     print(hla_a01_top_peptides)
    #     print('-'*80)
    #     print(hla_a02_top_peptides)

    print('Starting HLA1')
    hla_a01_top_peptides = peptide_filtering(hla_a01_df, 0.01, 0, inclusive=True)
    hla_a01_output_peptides = peptide_comparison(hla_a01_top_peptides, hla='HLA-A01')
    print('Starting HLA2')
    hla_a02_top_peptides = peptide_filtering(hla_a02_df, 0.01, 0, inclusive=True)
    hla_a02_output_peptides = peptide_comparison(hla_a02_top_peptides, hla='HLA-A02')


    return None


def peptide_filtering(df, cutoff_percentage, absolute_cutoff=0, inclusive=True):
    """
    This function will compare the found cryptic peptides to the canonical peptides, any overlap will result in removal
    of the cryptic peptide from the list, since it can also be a canonical peptide.
    :param df: Dataframe containing obtained peptides:
    :param: cutoff_percentage: cutoff for top x percentage of the peptides.
    :param: absolute_cutoff: when not null, absolute number of peptides is used instead of percentage.
    :param mode: definition of cryptic or canonical mode, known peptides lists will be used accordingly.
    :return: 3 lists of peptides; 1: list of wrongly classified peptides, 2: list of known peptides, 3: list of newly, not classified peptides.
    """
    top_df = pd.DataFrame()
    percentage_mode = True
    df = df.sort_values(by=['%Rank_EL'])

    if absolute_cutoff != 0:
        percentage_mode = False

    if percentage_mode:
        print('Using percentage cutoff')
        cutoff_quantile = cutoff_percentage / 100
        df_cutoff = df['%Rank_EL'].quantile(cutoff_quantile)
        df_mask = df['%Rank_EL'] <= df_cutoff
        top_df = df[df_mask]

    elif not percentage_mode:
        print('Using absolute cutoff')
        top_df = df[:absolute_cutoff]

    if inclusive:
        print('Inclusive mode enabled, adding all peptides with same EL rank as peptides in selected scope...')
        same_rank = top_df['%Rank_EL'].max()
        print(df)
        top_df = df[df['%Rank_EL'] <= same_rank]
        print('Amount of peptides:')
        print(top_df.shape)
        # with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
        #     print(top_df)
    return top_df


def peptide_comparison(df, hla, mode='cryp'):
    """

    :param df:
    :return:
    """
    text = ''
    if mode == 'can':
        text = 'Canonical'
    elif mode == 'cryp':
        text = 'Cryptic'

    print(f'Comparing {hla} peptides with known, {text} peptide database...')
    global project_dir
    comparison_df = pd.DataFrame()
    confirmed_canonical = pd.DataFrame()
    if hla == 'HLA-A01':
        if mode == 'cryp':
            comparison_df = pd.read_csv(Path(project_dir / 'Data' / 'Known_peptides' / 'HLA-A01_cryptic.csv'))
        elif mode == 'can':
            comparison_df = pd.read_csv(Path(project_dir / 'Data' / 'Known_peptides' / 'HLA-A01_canonical.csv'))
    elif hla == 'HLA-A02':
        if mode == 'cryp':
            comparison_df = pd.read_csv(Path(project_dir / 'Data' / 'Known_peptides' / 'HLA-A02_cryptic.csv'))
        elif mode == 'can':
            comparison_df = pd.read_csv(Path(project_dir / 'Data' / 'Known_peptides' / 'HLA-A02_canonical.csv'))

    try:
        merged_df = pd.merge(df, comparison_df, on=['Peptide'], how='inner')
        print('Canonical peptides have been detected! Removing from list...')
        confirmed_canonical = df[df['Peptide'].isin(merged_df['Peptide'])]
        print(confirmed_canonical.shape)
        best_ranking_peptides = df[~df['Peptide'].isin(merged_df['Peptide'])]  # Removed all known canonical peptides
        df = best_ranking_peptides
    except:
        print('No known peptides have been found, saving peptides...')

    return df, confirmed_canonical

if __name__ == "__main__":
    main()



# with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
#     print(df)