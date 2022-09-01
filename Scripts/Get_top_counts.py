########################################################################################
#  Script to fetch gene sequences based on coordinates
#  Input: Count file that contains Ensemble IDs
#  Output: Fasta file with headers, containing Ensemble information and NCBI sequence
#  Made by Shane Pullens, Utrecht University - Theoretical bioinformatics.
# Note, only need to connect to server when using NetMHCpan, and for deployment
########################################################################################


# imports
import time
# from Bio.Seq import Seq
# import ensembl_rest
import pprint
# import requests, sys
# import re
# import string
import os.path
# from os import path
from pathlib import Path
from datetime import datetime
import pandas as pd
# from tqdm.auto import tqdm  # Note This might be removed when we launch to webapp
# import numpy as np
# import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import seaborn as sns
import subprocess
# from sklearn.preprocessing import MaxAbsScaler, MinMaxScaler, StandardScaler, RobustScaler # Not needed since we use R for this.
# from rnalysis import filtering

#  Predefined variables
pp = pprint.PrettyPrinter(indent=4)  # Used to clearly print json files
server = "https://rest.ensembl.org"
verbose = False
project_dir = Path.cwd()
sns.set()  # Sets seaborn as default style


def main():
    global project_dir
    timestamp = datetime.now().strftime("%d-%b-%Y-h%H-m%M")

    # dir_name = 'gdc_download_20220422_083506.202427'  # TODO This should be user defined
    dir_name = 'breast_TCGA-BRCA_20220627_115207.223937'

    matrix_path = Path(project_dir / 'Output' / 'Counts' / dir_name)
    runDeseq(matrix_path)
    exit()

    input_paths = Path(project_dir / 'Output' / 'Counts' / dir_name / 'Raw_counts').glob('*counts.tsv')
    # Note: check if this can be done more nicely, is needed because without length, we can't make a loading bar, and
    #  we need to copy it, because else we destroy the generator
    input_path_len = len(list(Path(project_dir / 'Output' / 'Counts' / dir_name / 'Raw_counts').glob('*counts.tsv')))
    target_path = Path(project_dir / 'Output' / 'Fasta' / dir_name / timestamp)

    getReadsMatrix(input_paths, True, True)
    exit()
    # df = readFullDF(Path(project_dir / 'Output' / 'Counts' / dir_name / 'Full_dataframe.csv'))
    matrix = pd.read_csv(Path(project_dir / 'Output' / 'Counts' / dir_name / 'readsMatrix.csv'))

    exit()
    single_df_size = 60660


    # normalized_matrix = normalize(matrix) # Note: might be deprecated due to DEseq2
    # df_list = splice_df(df, single_df_size)
    # make_df_Plots(df, df_list)

    make_matrix_plots(matrix, Path(project_dir / 'Output' / 'Counts' / dir_name / 'readsMatrix.csv'))

    return None


def runDeseq(matrix_path):
    """
    Runs the DESeq function in a R subprocess to determine best genes from the counts data.
    R script an be found at ~/Scripts/R/deseq2_normalization.R
    :param matrix: Matrix that contains counts for each gene of each patient.
    :return: None, saves a list of significant genes instead named 'significant_genes.csv'
    """
    print('~'*80)
    print('Starting R script...')

    command = '/home/shane/miniconda3/envs/neoPipe/bin/Rscript'
    # args = ['--vanilla', str(matrix_path)] # Note: The vanilla argument seems to breaks R's argument parser.
    args = [str(matrix_path)]
    # path2script = '/home/shane/Documents/pycharm_project_519/Scripts/R/deseq2_normalization.R'
    path2script = '/home/shane/Documents/pycharm_project_519/Scripts/R/deseq2_normalization.R'

    retcode = subprocess.run([command, path2script] + args, universal_newlines=True)
    print(retcode)

    print('Finishing R script...')
    print('~'*80)
    return None


# def normalize(df, mode):
#     """
#     # Note: should be deprecated
#     :param df:
#     :param mode:
#     :return:
#     """
#     if mode == 'abs_max':
#         abs_scaler = MaxAbsScaler()
#         df_normalized = pd.DataFrame(abs_scaler.fit_transform(df), columns=df.columns)
#
#     if mode == 'min_max':
#         minmax_scaler = MinMaxScaler()
#         df_normalized = pd.DataFrame(minmax_scaler.fit_transform(df), columns=df.columns)
#
#     if mode == 'standarized':
#         std_scalar = StandardScaler()
#         df_normalized = pd.DataFrame(std_scalar.fit_transform(df), columns=df.columns)
#
#     return df_normalized


def make_matrix_plots(matrix, path):
    new_matrix = matrix.loc[:, matrix.columns != 'Gene']



    matrix_norm = normalize(new_matrix, mode='min_max')

    matrix_norm = pd.concat([matrix[['Gene']], matrix_norm], axis=1).set_index('Gene')


    # ax = sns.heatmap(matrix.set_index('Gene'))
    # ax.set_title('Heatmap of reads with no normalization (RNA only')
    # plt.tight_layout()
    # plt.show()
    #
    # ax = sns.heatmap(matrix_norm)
    # ax.set_title('Heatmap of reads with min-max normalization (RNA only)')
    # plt.tight_layout()
    # plt.show()
    # exit()

    exit()
    return None


def make_df_Plots(df, df_list):
    global project_dir
    # df = df_list[0] #TODO: remove Temp

    # >>>>> Gene type distribution plots
    # count_data = df.value_counts(['gene_type']).reset_index()
    # count_data.columns = ['gene_type', 'counts']
    #
    # count_data_top = count_data.loc[count_data['counts'] > (1000*473)]
    # count_data_top['gene_type'] = count_data_top['gene_type'].str.replace('_', '\n')
    # types_t = count_data_top['gene_type'].tolist()
    # amounts_t = count_data_top['counts'].tolist()
    #
    # count_data_bot = count_data[(count_data['counts'] < (1000*473)) & (count_data['counts'] > (50*473))]
    #
    # types_b = count_data_bot['gene_type'].tolist()
    # amounts_b = count_data_bot['counts'].tolist()
    # #
    # fig1 = plt.figure()
    # ax1 = fig1.add_subplot(111)
    # ax1.bar(types_t, amounts_t)
    # ax1.set_title('Gene types >1000 *473')
    # plt.xticks(rotation=45, ha='right')
    # plt.tight_layout()
    # plt.savefig(Path(project_dir / 'Figures' / 'Plots' / 'Genetype_dis_top.png'))
    #
    # fig2 = plt.figure()
    # ax2 = fig2.add_subplot(111)
    # ax2.bar(types_b, amounts_b)
    # ax2.set_title('Gene types 50 > 1000 *437')
    # plt.xticks(rotation=90, ha='right')
    # plt.tight_layout()
    # plt.savefig(Path(project_dir / 'Figures' / 'Plots' / 'Genetype_dis_bot.png'))
    # # print(count_data)
    #
    # count_data.to_csv(project_dir / 'Figures' / 'Plots' / 'Genetype_dis.csv')
    # plt.show()

    #<<<<< Gene_type distribution plots <<<<<

    # >>>>> Transcript plots >>>>> #TODO: make transcript plots using fetch transcript script.
    # Note: Still need  to make transcript plots.


    # <<<<< Transcript plots <<<<<
    # >>>>> Scatter plots >>>>>


    # # TODO: rerun to remove axis label
    # print('Running...')
    # fig3, axs = plt.subplots(2,3)
    # fig3.suptitle('Scatterplot of Unstranded counts VS. FPKM_unstranded counts')
    # fig3.supxlabel('Unstranded Counts')
    # fig3.supylabel('FPKM_unstranded Counts')
    # sns.regplot(ax=axs[0, 0], x=df_list[0]['unstranded'], y=df_list[0]['fpkm_unstranded'], line_kws={"color": "red"})
    # axs[0, 0].set_title('Patient #1')
    #
    # sns.regplot(ax=axs[0, 1], x=df_list[99]['unstranded'], y=df_list[99]['fpkm_unstranded'], line_kws={"color": "red"})
    # axs[0, 1].set_title('Patient #100')
    #
    # sns.regplot(ax=axs[0, 2], x=df_list[199]['unstranded'], y=df_list[199]['fpkm_unstranded'], line_kws={"color": "red"})
    # axs[1, 0].set_title('Patient #200')
    #
    # sns.regplot(ax=axs[1, 0], x=df_list[299]['unstranded'], y=df_list[299]['fpkm_unstranded'], line_kws={"color": "red"})
    # axs[1, 1].set_title('Patient #300')
    #
    # sns.regplot(ax=axs[1, 1], x=df_list[399]['unstranded'], y=df_list[399]['fpkm_unstranded'], line_kws={"color": "red"})
    # axs[1, 1].set_title('Patient #400')
    #
    # sns.regplot(ax=axs[1, 2], x=df['unstranded'], y=df['fpkm_unstranded'], line_kws={"color": "red"})
    # axs[1, 2].set_title('All Patients')
    #
    # for ax in axs.flat:
    #     ax.set_xlabel(None)
    #     ax.set_ylabel(None)
    # plt.tight_layout()
    # plt.show()
    # exit()


###
    # start_time = time.time()
    # print('0/7 Running...')
    # ax4 = sns.regplot(x=df['unstranded'], y=df['fpkm_unstranded'],
    #                   line_kws={"color": "red"})
    # ax4.set_title('Unstranded counts VS. FPKM_unstranded counts \n(in all patients)')
    # plt.tight_layout()
    # plt.savefig(Path(project_dir / 'Figures' / 'Plots' / 'unstranded_vs_fpkm_unstranded.png'))
    # plt.close()
    # print('1/7 saved unstranded_vs_fpkm_unstranded')
    # # plt.show()
    #
    # ax5 = sns.regplot(x=df['unstranded'], y=df['tpm_unstranded'],
    #                   line_kws={"color": "red"})
    # ax5.set_title('Unstranded counts VS. TPM_unstranded counts \n(in all patients)')
    # plt.tight_layout()
    # plt.savefig(Path(project_dir / 'Figures' / 'Plots' / 'unstranded_vs_tpm_unstranded.png'))
    # plt.close()
    # print('2/7 saved: unstranded_vs_tpm_unstranded')
    #
    # ax6 = sns.regplot(x=df['fpkm_unstranded'], y=df['tpm_unstranded'],
    #                   line_kws={"color": "red"})
    # ax6.set_title('FPKM_unstranded counts VS. TPM_unstranded counts \n(in all patients)')
    # plt.tight_layout()
    # plt.savefig(Path(project_dir / 'Figures' / 'Plots' / 'fpkm_unstranded_vs_tpm_unstranded.png'))
    # plt.close()
    # print('3/7 saved: fpkm_unstranded_vs_tpm_unstranded')
    #
    # ax7 = sns.regplot(x=df['fpkm_unstranded'], y=df['fpkm_uq_unstranded'],
    #                   line_kws={"color": "red"})
    # ax7.set_title('FPKM_unstranded counts VS. FPKM_uq_unstranded counts \n(in all patients)')
    # plt.tight_layout()
    # plt.savefig(Path(project_dir / 'Figures' / 'Plots' / 'FPKM_vs_FPKM_uq.png'))
    # plt.close()
    # print('4/7 saved: FPKM_vs_FPKM_uq')
    #
    # ax8 = sns.regplot(x=df['unstranded'], y=df['stranded_first'],
    #                   line_kws={"color": "red"})
    # ax8.set_title('unstranded counts VS. stranded_first counts \n(in all patients)')
    # plt.tight_layout()
    # plt.savefig(Path(project_dir / 'Figures' / 'Plots' / 'unstranded_vs_stranded_first.png'))
    # plt.close()
    # print('5/7 saved: unstranded_vs_stranded_first')
    # #
    # ax9 = sns.regplot(x=df['stranded_first'], y=df['stranded_second'],
    #               line_kws={"color": "red"})
    # ax9.set_title('stranded_first counts VS. stranded_second counts \n(in all patients)')
    # plt.tight_layout()
    # plt.savefig(Path(project_dir / 'Figures' / 'Plots' / 'stranded_second_vs_stranded_first.png'))
    # plt.close()
    # print('6/7 Saved: stranded_second_vs_stranded_first')
    #
    # ax10 = sns.regplot(x=df['unstranded'], y=df['stranded_second'],
    #                   line_kws={"color": "red"})
    # ax10.set_title('unstranded counts VS. stranded_second counts \n(in all patients)')
    # plt.tight_layout()
    # plt.savefig(Path(project_dir / 'Figures' / 'Plots' / 'unstranded_vs_stranded_second.png'))
    # plt.close()
    # print('7/7 saved: unstranded_vs_stranded_second')
    # print('Finished!')
    # plt.show()

###

    # <<<<< Scatter plots <<<<<
    # >>>>> Violin plots >>>>>
    print('Running...')
    #
    # # print(df_list[0]['unstranded'])
    # # print('min:\t', min(df_list[0]['unstranded']))
    # # print('max:\t', max(df_list[0]['unstranded']))
    starttime = time.time()



    #
    # df_list = df_list[0] # TEMP
    # df = df_list
    df = df[df['gene_type'].str.contains('RNA')]
    # df = df[df['unstranded'] > 1000*437]
    # count_data_bot = count_data[(count_data['counts'] < (1000*473)) & (count_data['counts'] > (50*473))]

    df_abs_max_scaled = normalize(df[['unstranded', 'stranded_first', 'stranded_second',
                                      'tpm_unstranded', 'fpkm_unstranded', 'fpkm_uq_unstranded']], mode='abs_max')
    df_min_max_scaled = normalize(df[['unstranded', 'stranded_first', 'stranded_second',
                                     'tpm_unstranded', 'fpkm_unstranded', 'fpkm_uq_unstranded']], mode='min_max')
    df_standarized = normalize(df[['unstranded', 'stranded_first', 'stranded_second',
                                  'tpm_unstranded', 'fpkm_unstranded', 'fpkm_uq_unstranded']], mode='standarized')

    fig9, axs1 = plt.subplots(2,2)
    fig9.suptitle('Violin plot of fpkm_unstranded counts with different \ndistribution methods. (RNA only)')
    # fig9.supxlabel('Distributions')

    sns.violinplot(ax=axs1[0, 0], x=df['fpkm_unstranded'])
    axs1[0, 0].set_title('No transformations')

    sns.violinplot(ax=axs1[0, 1], x=df_abs_max_scaled['fpkm_unstranded'])
    axs1[0, 1].set_title('Maximum absolute scaling')

    sns.violinplot(ax=axs1[1, 0], x=df_min_max_scaled['fpkm_unstranded'])
    axs1[1, 0].set_title('Min - Max scaling')

    sns.violinplot(ax=axs1[1, 1], x=df_standarized['fpkm_unstranded'])
    axs1[1, 1].set_title('Z-score scaling')
    plt.tight_layout()
    plt.show()

    # sns.boxplot(x=df_standarized['unstranded'])
    plt.tight_layout()
    plt.show()


    stoptime = time.time()

    print('runtime:\t', str(stoptime-starttime))

    # ax13 = sns.displot(df_list[0]["stranded_first"])
    # ax13.fig.suptitle('Distribution off "stranded_first" \n(in all patients)')
    # plt.tight_layout()
    # plt.show()
    exit()
    # <<<<< Violin plots <<<<<


    # stop_time = time.time()
    # print('Runtime:\t', stop_time-start_time)

    return None


def splice_df(df, df_size):
    amounts = df.shape[0] / df_size
    df.set_index(df.columns[0])
    # print(amounts)

    df_list = []
    for i in range(int(amounts)):
        j = i+1
        subdf = df[i*df_size:j*df_size]
        subdf = subdf.set_index(subdf['Unnamed: 0'])
        subdf = subdf[1:]
        # if i > 2:
        #     exit()
        df_list.append(subdf)
    return df_list


def getReadsMatrix(paths, RNAonly=True, saveMatrix=False):
    """

    :param paths:
    :param RNAonly:
    :param saveMatrix:
    :return:
    """
    print('RNA only:\t', RNAonly)
    counter = 0
    grouped_df = pd.DataFrame()

 #   file_path = ''
    for file_path in paths:
        counter += 1
        df = pd.read_csv(file_path, sep='\t', skiprows=6, header=None)
        df.columns = ['gene_id', 'gene_name', 'gene_type', 'unstranded', 'stranded_first', 'stranded_second',
                           'tpm_unstranded', 'fpkm_unstranded', 'fpkm_uq_unstranded']

        if RNAonly: # TODO: Got the task to add 'unprocessed_pseudogenes, now we also have 'transcribed_x & translated_x' added, these are however in low numbers.
                    # NOTE: Should those also be added to the pool?
            search_criteria = ['RNA', 'unprocessed_pseudogene']
            df = df[df['gene_type'].str.contains('|'.join(search_criteria))]

        df = df[df['gene_id'].str.contains('_PAR_Y') == False]  # remove all _PAR_Y genes due to gene name conflicts
        df = df.iloc[:, [0, 3]]
        df.columns = ['Gene_id', 'Patient_{i}'.format(i=counter)]
        # with pd.option_context('display.max_rows',1000, 'display.max_columns',1000):
        #     print(df)

        if counter <= 1:
            grouped_df = df
        else:
            grouped_df = pd.concat([grouped_df, df['Patient_{i}'.format(i = counter)]], axis=1)

    if saveMatrix:
        print('Trying to save matrix')
        grouped_df.to_csv(Path(file_path.parent.parent / 'readsMatrix.csv'), index=False)
        print('successfully wrote reads matrix')

    # Note: this should maybe be deprecated
    # if makeFullDF:
    #     full_df = pd.DataFrame()
    #     for file_path in paths:
    #         df = pd.read_csv(file_path, sep='\t', skiprows=6, header=None)
    #         full_df = pd.concat([full_df, df])
    #         full_df.columns = ['gene_id', 'gene_name', 'gene_type', 'unstranded', 'stranded_first', 'stranded_second',
    #                            'tpm_unstranded', 'fpkm_unstranded', 'fpkm_uq_unstranded']
    #         full_df.to_csv(Path(file_path.parent / 'Full_dataframe.csv'), index=False)
    #     print('successfully wrote full df.')
    return grouped_df


def readFullDF(path):
    df = pd.read_csv(path)
    # with pd.option_context('display.max_rows',1000, 'display.max_columns',1000):
    #     print(df)
    return df


def read_data(path):
    file = open(path, 'r')
    file = file.readlines()
    gene_score_dict = {}

    for line in file[6:]:  # Skip the first 6 header lines
        gene_id, gene_name, gene_type, score, _, _, _, _, _, = line.replace('\n', '').split('\t')
        gene_id = gene_id.split('.')[0]
        gene_score_dict[gene_id] = [score, gene_name, gene_type]

    return gene_score_dict


def write_file(text, filename, path, timestamp, add_timestamp=True, overwrite_check=False):
    """
    # TODO Add documentation
    :param text:
    :param filename:
    :param path:
    :param timestamp:
    :param add_timestamp:
    :param overwrite_check:
    :return:
    """
    if add_timestamp:  # Check if timestamp should be implemented
        filename = timestamp + filename  #+'.fasta' This was only here for testing right?

    # Check if the user already added an extension, skip adding another.
    filename = Path(filename)
    if filename.suffix not in ('.fasta', '.fa'):
        filename = str(filename) + '.fasta'

    absolute_path = Path(path / filename)

    # Check if the user enabled overwrite check. This checks if the user is going to overwrite a file and stops it.
    if overwrite_check:
        if os.path.isfile(absolute_path):
            print('>!! File: ' + filename + ' already exists, aborted writing.')
            return None

    try:
        with open(absolute_path, 'w') as f:
            f.write(text)
        if verbose:
            print('>>> succesfully written '+filename)
    except:
        print('! An error occurred while trying to write '+filename)
    return None


def comparelists(a,b,c):
    z = [x for x in b
         if x not in a and x not in c] #list comprehension
    print(z)


# Note: This might need to be changed when adding the scripts to a pipeline
if __name__ == "__main__":
    main()

######## Code snippets #########
# if os.path.isdir(target_path):
#     if verbose:
#         print('Path already exist.')
# else:
#     os.makedirs(target_path, exist_ok=True)

# outer_pbar = tqdm(total=input_path_len, desc='File: ', leave=True)





    # full_df_dict = {}
    # gene_df_dict = {}
    # full_gene_df = pd.DataFrame()

    # for file_path in input_paths:
        # print(file_path.stem)

        # gene_dict = read_data(file_path)  # Testing


        # full_gene_df = pd.concat([full_gene_df, df[1]])
        # gene_df_dict[file_path.stem] = df[1]
        # full_df_dict[file_path.stem] = df  # TODO: consider removing '.rna_seq...... from key'


        # with pd.option_context('display.max_rows',1000, 'display.max_columns',1000):
        #     print(df)

        # filestem = Path(file_path).stem

        # for gene_id, dict_values in gene_dict.items():  #in tqdm(gene_dict.items(), desc='Gene ID: ', leave=False, position=0+1):
        #     print(gene_id)
        #     print(dict_values)


    # full_df.reset_index(drop=True, inplace=True)

# with pd.option_context('display.max_rows',1000, 'display.max_columns',1000):
#     print(full_df)