# Note: this is a test environment, here we test the possibility for using subprocess on mhcpan.
import pandas as pd
import subprocess
from pathlib import Path
import os
from tqdm import tqdm
import shutil
from queue import Queue
from threading import Thread
from time import time


def main():
    project_dir = Path.cwd()
    new_file_name = 'skin_TCGA-SKCM_20220822_132048.170180'
    target_fasta = 'testing_sequence_mhcPan.fasta'
    filepath = Path(project_dir / 'Output' / 'Counts' / new_file_name)

    alleles = ['HLA-A01:01', 'HLA-A02:01']
    alleles = ['HLA-A01:01']
    run_mhcPan(alleles, filepath, target_fasta, canonical=False, safe_mode='all')
    exit()
    return None


def mhcPan_thread_ripper(path, file_name, n_threads):
    """
    This function rips apart the fasta file that will be used for MHCpan. Because MHC pan only allows fasta FILES
    as input and not fasta STRINGs, we need to split the fasta file up into blocks.
    This function splits up the original fasta file into the same amount as selected threads, and safes them to a
    temporary folder, whilst returning the fasta chunk file names.
    in a temp
    :param path: Path to data
    :param file_name: name of imput fasta file
    :return: list of fasta chunk file names.
    """
    file = open(Path(path / file_name))
    fasta_list = file.read().replace('>', '$$>').split('$$')
    list_len = len(fasta_list)
    list_mod = list_len % n_threads

    chunk_size = int((list_len - list_mod) / n_threads)
    fasta_chunks = [fasta_list[x:x + chunk_size] for x in range(0, len(fasta_list), chunk_size)]

    fasta_chunks[-2] = fasta_chunks[-2] + fasta_chunks[-1]  # adds the final bit of reads to the last full block
    del fasta_chunks[-1]

    os.mkdir(Path(path / 'Tmp'))
    filenames = []
    for i, fasta_chunk in enumerate(fasta_chunks):
        fasta_chunk = ''.join(str(n) for n in fasta_chunk)
        fasta_chunk_filename = file_name.replace('.fasta', '_{}.fasta'.format(i+1))
        with open(Path(path / 'Tmp' / fasta_chunk_filename), 'w') as f:
            f.write(fasta_chunk)
        filenames.append(fasta_chunk_filename)
    return filenames


def mhcPan_thread_assembly(path, rm_tmp_dir=True):
    """

    :param path:
    :param rm_tmp_dir:
    :return:
    """
    whole_df = pd.DataFrame()
    for peptide_filename in tqdm(path.glob('*_cryptic_peptides.csv')):
        peptide_df = pd.read_csv(peptide_filename)
        whole_df = pd.concat([whole_df, peptide_df])

    whole_df.to_csv(Path(path.parent / 'MHCpan_output.csv'))
    if rm_tmp_dir:
        try:
            shutil.rmtree(path)
            print('Successfully removed temporary directory.')
        except OSError as e:
            print("Error: %s - %s." % (e.filename, e.strerror))

    return None


def run_mhcPan_threading(alleles, filepath, filename, canonical, safe_mode='SB'):
    """

    :param alleles:
    :param filepath:
    :param fasta:
    :param canonical:
    :param safe_mode:
    :return:
    """
    if safe_mode == 'SB' or safe_mode == 'sb':  # strong binders
        safe_mode = 'SB'
        # print('Selected mode:\t', safe_mode)
    elif safe_mode == 'WB' or safe_mode == 'wb':  # weak binders
        safe_mode = 'WB'
        # print('Selected mode:\t', safe_mode)
    elif safe_mode == 'AB' or safe_mode == 'ab':  # all binders
        safe_mode = 'AB'
        # print('Selected mode:\t', safe_mode)
    elif safe_mode == 'all' or safe_mode == 'ALL' or safe_mode == 'a':  # Whole output
        safe_mode = 'all'
        # print('Selected mode:\t', safe_mode)
    else:
        # print("Didn't select a proper safe mode, automatically sleected SB")
        safe_mode = 'SB'

    inner_netMHCpan_output_list = []
    netMHCpan_output_df = ''
    for x, allele in enumerate(alleles):
        # print('Working on binding prediction on allele:\t', allele)
        command = 'netMHCpan'
        args = [f"-f {filepath / filename} -s 1 -l 9 -a {allele}"]  # >> {filepath / outfile}"]
        path2script = '/home/shane/Tools/netMHCpan-4.1/'
        retcode = subprocess.check_output([command, path2script] + args, universal_newlines=True)

        # Takes the raw output, splits every line, filters out the lines that start with '#', then removes the empty
        # items in the list (these originate from the white lines), ## then joins it all together to a single string.
        netMHCpan_output_lines = list(filter(None, [s for s in retcode.split('\n') if '#' not in s]))

        # Removes leading and footer lines from standard netMHCpan output.
        del netMHCpan_output_lines[:2]
        del netMHCpan_output_lines[1]

        netMHCpan_list = []
        for i, j in enumerate(netMHCpan_output_lines):
            if x == 0:
                if i == 0:
                    netMHCpan_list.append('\t'.join(j.split()).split('\t'))
            else:
                line = '\t'.join(j.split()).split('\t')
                if not j.startswith('Pos'):  # removes unwanted headers.
                    if line[0].isdigit():  # removes extra footers
                        if len(line) == 18:
                            del line[-2]  # Removes the '<=' from the list
                        netMHCpan_list.append(line)

        inner_netMHCpan_output_list.append(netMHCpan_list)
    netMHCpan_output_df = pd.DataFrame(inner_netMHCpan_output_list[1:], columns=inner_netMHCpan_output_list[0])
        # with pd.option_context('display.max_rows', 100, 'display.max_columns', 20):
        #     print(netMHCpan_output_df)

    timestamp = ''
    if canonical:
        out_file_name = filename.replace('.fasta', '_canonical_peptides.csv')
    else:
        out_file_name = filename.replace('.fasta', '_cryptic_peptides.csv')
    if safe_mode == 'all':
        # print('Safe all mode enabled, saving full netMHCpan output...')
        netMHCpan_output_df.to_csv(Path(filepath / out_file_name), index=False)
        # write_file(full_netMHCpan_output, outfile, filepath, timestamp, False, True, False)
        print('Successfully saved full netMHCpan output')
        return netMHCpan_output_df
    elif safe_mode == 'AB':
        netMHCpan_output_AB = netMHCpan_output_df[netMHCpan_output_df['BindLevel'].notnull()]
        netMHCpan_output_AB.to_csv(Path(filepath / out_file_name), index=False)
        if netMHCpan_output_AB.isnull:
            # print('No binders found')
            return netMHCpan_output_AB

    elif safe_mode == 'WB':
        netMHCpan_output_WB = netMHCpan_output_df[netMHCpan_output_df['BindLevel'] == 'WB']
        netMHCpan_output_WB.to_csv(Path(filepath / out_file_name), index=False)
        if netMHCpan_output_WB.isnull:
            print('No Weak binders found')
        return netMHCpan_output_WB

    elif safe_mode == 'SB':
        netMHCpan_output_SB = netMHCpan_output_df[netMHCpan_output_df['BindLevel'] == 'SB']

        netMHCpan_output_SB.to_csv(Path(filepath / out_file_name), index=False)
        if netMHCpan_output_SB.isnull:
            print('No Strong binders found')
        return netMHCpan_output_SB

    # print('Finishing MHCpan...')
    # print('~'*80)

    return None


def run_mhcPan(alleles, filepath, filename, canonical, safe_mode='SB'):
    """

    :param alleles:
    :param filepath:
    :param fasta:
    :param canonical:
    :param safe_mode:
    :return:
    """

    if safe_mode == 'SB' or safe_mode == 'sb':  # strong binders
        safe_mode = 'SB'
        print('Selected mode:\t', safe_mode)
    elif safe_mode == 'WB' or safe_mode == 'wb':  # weak binders
        safe_mode = 'WB'
        print('Selected mode:\t', safe_mode)
    elif safe_mode == 'AB' or safe_mode == 'ab':  # all binders
        safe_mode = 'AB'
        print('Selected mode:\t', safe_mode)
    elif safe_mode == 'all' or safe_mode == 'ALL' or safe_mode == 'a':  # Whole output
        safe_mode = 'all'
        print('Selected mode:\t', safe_mode)
    else:
        print("Didn't select a proper safe mode, automatically sleected SB")
        safe_mode = 'SB'

    print('~' * 80)
    print('Starting netMHCpan...')
    outfile = ''
    if canonical:
        print('Detecting binders in canonical proteins')
        outfile = 'canonical_peptides.txt'
    elif not canonical:
        print('Detecting binders in cryptic proteins')
        outfile = 'cryptic_peptides.txt'

    full_netMHCpan_output = ''
    netMHCpan_output_df = ''
    for allele in alleles:
        print('Working on binding prediction on allele:\t', allele)
        command = 'netMHCpan'
        print("Running command:\t", f"-f{filepath / filename} -s 1 -l 9 -a {allele} -BA 1 >> {filepath / outfile}")
        args = [f"-f {filepath / filename} -s 1 -l 9 -a {allele} -BA 1 >> {filepath / outfile}"]
        path2script = '/home/shane/Tools/netMHCpan-4.1/'
        retcode = subprocess.check_output([command, path2script] + args, universal_newlines=True)

        # Takes the raw output, splits every line, filters out the lines that start with '#', then removes the empty
        # items in the list (these originate from the white lines), ## then joins it all together to a single string.
        netMHCpan_output_lines = list(filter(None, [s for s in retcode.split('\n') if '#' not in s]))

        # Removes leading lines from standard netMHCpan output.
        del netMHCpan_output_lines[:2]
        del netMHCpan_output_lines[1]

        netMHCpan_list = []
        for i, j in enumerate(netMHCpan_output_lines):
            if i == 0:
                netMHCpan_list.append('\t'.join(j.split()).split('\t'))
            else:
                line = '\t'.join(j.split()).split('\t')
                if not j.startswith('Pos'):  # removes unwanted headers.
                    if line[0].isdigit():  # removes extra footers
                        if len(line) == 18:
                            del line[-2]  # Removes the '<=' from the list
                        netMHCpan_list.append(line)

        netMHCpan_output_df = pd.DataFrame(netMHCpan_list[1:], columns=netMHCpan_list[0])
        # with pd.option_context('display.max_rows', 900, 'display.max_columns', 4):
        #     print(netMHCpan_output_df)

    timestamp = ''
    if safe_mode == 'all':
        print('Safe all mode enabled, saving full netMHCpan output...')
        netMHCpan_output_df.to_csv(Path(filepath / outfile), index=False)
        # write_file(full_netMHCpan_output, outfile, filepath, timestamp, False, True, False)
        print('Successfully saved full netMHCpan output')
    elif safe_mode == 'AB':
        netMHCpan_output_AB = netMHCpan_output_df[netMHCpan_output_df['BindLevel'].notnull()]
        netMHCpan_output_AB.to_csv(Path(filepath / outfile), index=False)
        if netMHCpan_output_AB.isnull:
            print('No binders found')

    elif safe_mode == 'WB':
        netMHCpan_output_WB = netMHCpan_output_df[netMHCpan_output_df['BindLevel'] == 'WB']
        netMHCpan_output_WB.to_csv(Path(filepath / outfile), index=False)
        if netMHCpan_output_WB.isnull:
            print('No Weak binders found')

    elif safe_mode == 'SB':
        netMHCpan_output_SB = netMHCpan_output_df[netMHCpan_output_df['BindLevel'] == 'SB']
        netMHCpan_output_SB.to_csv(Path(filepath / outfile), index=False)
        if netMHCpan_output_SB.isnull:
            print('No Strong binders found')

    print('Finishing MHCpan...')
    print('~' * 80)

    return None


if __name__ == "__main__":
    main()
