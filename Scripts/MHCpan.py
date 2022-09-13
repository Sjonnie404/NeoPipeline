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
    new_file_name = 'breast_TCGA-BRCA_20220902_080722.028812'
    target_fasta = 'testing_mhcpan.fasta'
    filepath = Path(project_dir / 'Output' / 'Counts' / new_file_name)

    alleles = ['HLA-A01:01', 'HLA-A02:01']
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
    fasta_list = list(filter(None, file.read().replace('>', '$$>').split('$$')))
    fasta_list = [x for x in fasta_list if x.startswith('>')]
    list_len = len(fasta_list)

    fasta_chunks = []  # When we have a nice amount of reads, we can simply split on the headers
    if n_threads > list_len:
        print('Workload is smaller than #threads, lowering #threads to ', str(list_len))
        n_threads = list_len
        fasta_chunks = fasta_list

    else:  # If we don;t have a nice list, we need to spread the reads over the threads.
        list_mod = list_len % n_threads
        chunk_size = int((list_len - list_mod) / n_threads)
        fasta_chunks = [fasta_list[x:x + chunk_size] for x in range(0, list_len, chunk_size)]

        fasta_chunks[-2] = fasta_chunks[-2] + fasta_chunks[-1]  # adds the final bit of reads to the last full block
        del fasta_chunks[-1]
    try:
        os.mkdir(Path(path / 'Tmp'))
    except:
        print('Tmp dir already exist')
    filenames = []
    for i, fasta_chunk in enumerate(fasta_chunks):
        fasta_chunk = ''.join(str(n) for n in fasta_chunk)
        fasta_chunk_filename = file_name.replace('.fasta', '_{}.fasta'.format(i+1))
        with open(Path(path / 'Tmp' / fasta_chunk_filename), 'w') as f:
            f.write(fasta_chunk)
        filenames.append(fasta_chunk_filename)
    return filenames, n_threads


def mhcPan_thread_assembly(path, canonical):
    """

    :param path:
    :param canonical:
    :return:
    """
    file_search = ''
    if canonical:
        file_search = '*_canonical_peptides.csv'
    if not canonical:
        file_search = '*_cryptic_peptides.csv'

    whole_df = pd.DataFrame()
    for peptide_filename in path.glob(file_search):
        # print(peptide_filename)
        peptide_df = pd.read_csv(peptide_filename)
        whole_df = pd.concat([whole_df, peptide_df])

    if canonical:
        whole_df.to_csv(Path(path.parent / 'MHCpan_output_canonical.csv'))
    if not canonical:
        whole_df.to_csv(Path(path.parent / 'MHCpan_output_cryptic.csv'))
    return None


def remove_tmp_dir(path):
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

    # print('Selected safe mode:\t', safe_mode)
    netMHCpan_list = []
    netMHCpan_output_df = pd.DataFrame()
    for x, allele in enumerate(alleles):
        # print('Working on binding prediction on allele:\t', allele)
        command = 'netMHCpan'
        args = [f"-f {filepath / filename} -s 1 -l 9 -a {allele} -BA 1"]  # >> {filepath / outfile}"]
        path2script = '/home/shane/Tools/netMHCpan-4.1/'
        retcode = subprocess.check_output([command, path2script] + args, universal_newlines=True)

        # Takes the raw output, splits every line, filters out the lines that start with '#', then removes the empty
        # items in the list (these originate from the white lines), ## then joins it all together to a single string.
        netMHCpan_output_lines = list(filter(None, [s for s in retcode.split('\n') if '#' not in s]))
        # Removes leading and footer lines from standard netMHCpan output.

        del netMHCpan_output_lines[:2]
        del netMHCpan_output_lines[1]

        for i, j in enumerate(netMHCpan_output_lines):
            if i == 0:
                if x == 0:
                    netMHCpan_list.append('\t'.join(j.split()).split('\t'))

            elif i != 0:
                line = '\t'.join(j.split()).split('\t')
                if not j.startswith('Pos'):  # removes unwanted headers.
                    if line[0].isdigit():  # removes extra footers
                        if len(line) == 18:
                            del line[-2]  # Removes the '<=' from the list
                        elif len(line) == 16:  # Note: this is al hardcoded and will break when MHCpan changes
                            line.append('NB')
                        try:
                            netMHCpan_list.append(line)
                        except:
                            print("Something went wrong when trying to append the following line:")
                            print("!"*80)
                            print(line)
                            print("!"*80)  # Note: debugging

    try:
        netMHCpan_output_df = pd.DataFrame(netMHCpan_list[1:], columns=netMHCpan_list[0])
    except:
        print('Something went wrong with the dataframe!')
        for item in netMHCpan_list[1:]:
            if len(item) != 17:
                print(item)
        print('---')
        print(netMHCpan_list[0])
        print('Something went wrong with the dataframe...')
        exit()

    timestamp = ''
    if canonical:
        out_file_name = filename.replace('.fasta', '_canonical_peptides.csv')
    else:
        out_file_name = filename.replace('.fasta', '_cryptic_peptides.csv')
    if safe_mode == 'all':
        # print('Safe all mode enabled, saving full netMHCpan output...')
        netMHCpan_output_df.to_csv(Path(filepath / out_file_name), index=False)
        print('Saved output')
        # write_file(full_netMHCpan_output, outfile, filepath, timestamp, False, True, False)
        # print('Successfully saved full netMHCpan output')
        return netMHCpan_output_df
    elif safe_mode == 'AB':
        netMHCpan_output_AB = netMHCpan_output_df[netMHCpan_output_df['BindLevel'] != 'NB']
        if netMHCpan_output_AB.empty:
            print(netMHCpan_output_AB)
            print('No binders found')
            return netMHCpan_output_AB
        netMHCpan_output_AB.to_csv(Path(filepath / out_file_name), index=False)
        print('Saved output')

    elif safe_mode == 'WB':
        netMHCpan_output_WB = netMHCpan_output_df[netMHCpan_output_df['BindLevel'] == 'WB']
        if netMHCpan_output_WB.empty:
            print('No Weak binders found')
            return netMHCpan_output_WB
        netMHCpan_output_WB.to_csv(Path(filepath / out_file_name), index=False)
        print('Saved output')

    elif safe_mode == 'SB':
        netMHCpan_output_SB = netMHCpan_output_df[netMHCpan_output_df['BindLevel'] == 'SB']
        if netMHCpan_output_SB.empty:
            print('No Strong binders found')
            return netMHCpan_output_SB
        netMHCpan_output_SB.to_csv(Path(filepath / out_file_name), index=False)
        print('Saved output')

    print('Finishing MHCpan...')
    print('~'*80)

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
        outfile = 'canonical_peptides_test.txt'
    elif not canonical:
        print('Detecting binders in cryptic proteins')
        outfile = 'cryptic_peptides_test.txt'

    netMHCpan_list = []
    for allele in alleles:
        print('Working on binding prediction on allele:\t', allele)
        command = 'netMHCpan'
        # print("Running command:\t", f"-f{filepath / filename} -s 1 -l 9 -a {allele} -BA 1 >> {filepath / outfile}")
        args = [f"-f {filepath / filename} -s 1 -l 9 -a {allele} -BA 1 >> {filepath / outfile}"]
        path2script = '/home/shane/Tools/netMHCpan-4.1/'
        retcode = subprocess.check_output([command, path2script] + args, universal_newlines=True)
        # Takes the raw output, splits every line, filters out the lines that start with '#', then removes the empty
        # items in the list (these originate from the white lines), ## then joins it all together to a single string.
        netMHCpan_output_lines = list(filter(None, [s for s in retcode.split('\n') if '#' not in s]))

        # Removes leading lines from standard netMHCpan output.
        del netMHCpan_output_lines[:2]
        del netMHCpan_output_lines[1]

        for i, j in enumerate(netMHCpan_output_lines):
            if i == 0:
                netMHCpan_list.append('\t'.join(j.split()).split('\t'))
            else:
                line = '\t'.join(j.split()).split('\t')
                if not j.startswith('Pos'):  # removes unwanted headers.
                    if line[0].isdigit():  # removes extra footers
                        if len(line) == 15:
                            del line[-2]  # Removes the '<=' from the list
                        elif len(line) == 13:
                            line.append('NB')
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
        netMHCpan_output_AB = netMHCpan_output_df[netMHCpan_output_df['BindLevel'] != 'NB']
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






