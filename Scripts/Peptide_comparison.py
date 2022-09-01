########################################################################################
#  Script for all MHCpan output functions: Analytics & comparisions.
#  Made by Shane Pullens, Utrecht University - Theoretical bioinformatics.
#  Version 1.0
########################################################################################



# imports
from pathlib import Path

import pandas as pd


def main():
    project_dir = Path.cwd()
    new_file_name = 'skin_TCGA-SKCM_20220822_132048.170180'
    mhcpan_output = 'MHCpan_output.csv'
    filepath = Path(project_dir / 'Output' / 'Counts' / new_file_name / mhcpan_output)

    peptide_df = pd.read_csv(filepath)
    print('got here')
    exit()
    return None



def peptide_analystics(data):

    return None


def read_mhcPan_output(path):
    """

    :param path:
    :return:
    """
    df = pd.read_csv(path)

    return df


if __name__ == "__main__":
    main()