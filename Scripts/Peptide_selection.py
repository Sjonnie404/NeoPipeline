########################################################################################
#  Script to rank the peptides & compare them to known peptide files.
#  Made by Shane Pullens, Utrecht University - Theoretical bioinformatics.
#  Version 1.0
########################################################################################

# Imports
import pandas as pd
from pathlib import Path

project_dir = Path.cwd()


def main():
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
    df.sort_values(by=['%Rank_EL'], inplace=True)

    # When an absolute cutoff has been selected, it will overrule the percentage cutoff.
    if absolute_cutoff != 0:
        percentage_mode = False

    # Save top x% of the best ranking peptides
    if percentage_mode:
        print('Using percentage cutoff')
        cutoff_quantile = cutoff_percentage / 100
        df_cutoff = df['%Rank_EL'].quantile(cutoff_quantile)
        df_mask = df['%Rank_EL'] <= df_cutoff
        top_df = df[df_mask]

    # Save top x number of best ranking peptides.
    elif not percentage_mode:
        print('Using absolute cutoff')
        top_df = df[:absolute_cutoff]

    # Inclusive mode checks if there are peptides that fall out of scope due to the threshold, but have the same rank.
    # When the same rank is detected, they are still added to the candidate list, overruling the threshold.
    if inclusive:
        print('Inclusive mode enabled, adding all peptides with same EL rank as peptides in selected scope...')
        same_rank = top_df['%Rank_EL'].max()
        print(df)
        top_df = df[df['%Rank_EL'] <= same_rank]
        print('Amount of peptides:')
        print(top_df.shape)
    return top_df


def peptide_comparison(df, hla):
    """
    This method compares our candidate peptides with known peptides from our database.
    If they are found, they are moved from the candidate list to the false positive list.
    :param df: list of candidate peptides
    :param hla: HLA molecule to compare to
    :return df: filtered candidate peptide list with solely unknown peptides
    :return confirmed_canonical: peptide list of peptides that have been found in other literature.
    """
    print(f'Comparing {hla} peptides with known peptide database...')
    global project_dir
    comparison_df = pd.DataFrame()
    confirmed_canonical = pd.DataFrame()
    if hla == 'HLA-A01':
        comparison_df = pd.read_csv(Path(project_dir / 'Data' / 'Known_peptides' / 'HLA-A01_canonical.csv'))
    elif hla == 'HLA-A02':
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