########################################################################################
#  Script to fetch gene sequences based on coordinates
#  Input: Count file that contains Ensemble IDs
#  Output: Fasta file with headers, containing Ensemble information and NCBI sequence
#  Made by Shane Pullens, Utrecht University - Theoretical bioinformatics.
# Note, only need to connect to server when using NetMHCpan, and for deployment
########################################################################################
# TODO: Add sys arguments
# Nice to have: GIU style?
# TODO: Add logs, use logging module
# Nice to have: Summary output, which could be used in the logs.
# TODO: This needs to be rewritten to take gene_list input instead of tsv (I think this is done?)

#  Imports
# from Bio import Entrez, SeqIO
import itertools
import time

from Bio.Seq import Seq
# import ensembl_rest
import pprint
import requests, sys
import re
import string
import os.path
# from os import path
from pathlib import Path
from datetime import datetime
from tqdm.auto import tqdm  # Note This might be removed when we launch to webapp

#  Predefined variables
pp = pprint.PrettyPrinter(indent=4)  # Used to clearly print json files
server = "https://rest.ensembl.org"
verbose = False


def main(verbose=False):
    # TODO: CLEANUP
    timestamp = datetime.now().strftime("%d-%b-%Y-h%H-m%M")
    project_dir = Path.cwd()
    #dir_name = 'gdc_download_20220404_113548.750166'  # TODO This should be user defined

    dir_name = 'testing_breast'
    file_name = 'overlapping_genes.txt'

    # input_path = Path(project_dir / 'Output' / 'Counts' / dir_name).glob('*.tsv')
    # input_path_len = len(list(Path(project_dir / 'Output' / 'Counts' / dir_name).glob('*tsv'))) # Note: check if this can be done more nicely
    # input_path = Path(project_dir / 'Output' / 'Counts' / dir_name).glob('*script.tsv')
    # input_path_len = len(list(Path(project_dir / 'Output' / 'Counts' / dir_name).glob('*script.tsv')))

    # TODO: The output should be make more clear than what it is now.
    input_path = Path(project_dir / 'Data' / dir_name / file_name)
 #   input_path_len = len(list(Path(project_dir / 'Data' / file_name)))
    # target_path = Path(project_dir / 'Output' / 'Fasta' / dir_name / timestamp)
    target_path = Path(project_dir / 'Output' / dir_name / 'Fasta')

    if os.path.isdir(target_path):
        if verbose:
            print('Path already exist.')
    else:
        os.makedirs(target_path, exist_ok=True)

    file_path = input_path

    ## Old code ###
    # Nice to have: make this run in parallel

    # outer_pbar = tqdm(total=input_path_len, desc='File: ', leave=True)
    # for file_path in input_path:
    #     #print(file_path)
    #     gene_dict = read_data(file_path)  # Testing
    #     # gene_scores = read_data(Path(project_dir / 'Counts' / ))
    #     backbone_fasta_output = ''
    #     trans_fasta_output = ''
    #     # filename = Path(file_path).name
    #     filestem = Path(file_path).stem
    #
    #     for gene_id, dict_values in tqdm(gene_dict.items(), desc='Gene ID: ', leave=False, position=0+1):
    #         if verbose:
    #             print('>>> Started process for:\t' + str(gene_id))
    #
    #         backbone_trans_id, extend_5, extend_3, trans_id_list = get_Ensenble_data(gene_id, verbose)
    #         backbone_fasta_seq = get_backbone_sequence(backbone_trans_id, extend_5, extend_3, gene_id, dict_values, verbose)
    #         trans_fasta_seqs = get_sequence(trans_id_list, gene_id, dict_values, verbose)
    #
    #         backbone_fasta_output = backbone_fasta_output +backbone_fasta_seq
    #         trans_fasta_output = trans_fasta_output +trans_fasta_seqs
    #         if verbose:
    #             print('>>> Finished process for:\t' + str(gene_id))
    #
    #
    #     write_file(backbone_fasta_output, filestem+'_backbone', target_path, timestamp, True, True)
    #     write_file(trans_fasta_output, filestem+'_transcripts', target_path, timestamp, True, True)
    #     outer_pbar.update(1)


    #print(file_path)
    #gene_dict = read_data(file_path)  # Testing
    file = open(file_path, 'r')
    file = file.readlines()
    gene_list = []
    for line in file:
        # print(line.replace('"','').replace("\n","").split('.')[0])
        gene_id = line.replace('"','').replace("\n","").split('.')[0]
        gene_list.append(gene_id)

    # gene_scores = read_data(Path(project_dir / 'Counts' / ))
    backbone_fasta_output = ''
    trans_fasta_output = ''
    # filename = Path(file_path).name
    filestem = Path(file_path).stem
    gene_dict = {}  # tmp
    for gene_id in gene_list:
        if verbose:
            print('>>> Started process for:\t' + str(gene_id))

        backbone_trans_id, extend_5, extend_3, trans_id_list = get_Ensenble_data(gene_id, verbose=False)
        backbone_fasta_seq = get_backbone_sequence(backbone_trans_id, extend_5, extend_3, gene_id, verbose=False)
        trans_fasta_seqs = get_sequence(trans_id_list, gene_id, verbose=False)

        backbone_fasta_output = backbone_fasta_output +backbone_fasta_seq
        trans_fasta_output = trans_fasta_output +trans_fasta_seqs
        if verbose:
            print('>>> Finished process for:\t' + str(gene_id))

    write_file(backbone_fasta_output, filestem+'_backbone_TEST', target_path, timestamp, True, True)
    write_file(trans_fasta_output, filestem+'_transcripts_TEST', target_path, timestamp, True, True)
    return None


def read_data(path):
    file = open(path, 'r')
    file = file.readlines()
    gene_score_dict = {}

    for line in file[6:]:  # Skip the first 6 header lines
        gene_id, gene_name, gene_type, score, _, _, _, _, _, = line.replace('\n', '').split('\t')
        gene_id = gene_id.split('.')[0]
        gene_score_dict[gene_id] = [score, gene_name, gene_type]

    return gene_score_dict


def get_Ensenble_data(gene_id='ENSG00000157764', verbose=False):
    global server  # Add this as global variable since it won't change.
    ext = "/lookup/id/"+str(gene_id)+"?expand=1"
    gene_r = requests.get(server + ext, headers={"Content-Type": "application/json"})

    if not gene_r.status_code == 200:
        print('Error occurred whilst fetching url:\t', ext)
        #print('Adding to errornous genelist...')
        # raise Exception('Bad response')
        # TODO: find way to save 'bad' genes.
        return '', 1, 1, ['']

    decoded = gene_r.json()
    trans_dict = decoded.get('Transcript')
    trans_dict_own = {}

    if verbose:
        print('  > Found '+str(len(trans_dict))+' transcripts for '+str(gene_id)+'!')

    for trans in trans_dict:
        start = trans.get('start')
        end = trans.get('end')
        id = trans.get('id')
        length = end-start
        trans_dict_own.update({id: [length, start, end]})

    trans_id_list = list(trans_dict_own.keys())
    trans_list_des = sorted(trans_dict_own.items(), key=lambda x: x[1], reverse=True)
    backbone_trans_id = trans_list_des[0][0]  # Extracts actual Id from dict list.
    expected_start = trans_dict_own.get(backbone_trans_id)[1]
    expected_stop = trans_dict_own.get(backbone_trans_id)[2]

    if verbose:
        print('  > Selected backbode ID: '+str(backbone_trans_id))
        print('  > Initial start location: '+str(expected_start))
        print('  > Initial stop location: ' + str(expected_stop))
        print('  > Checking for better start&stop locations...')

    real_start = expected_start
    real_stop = expected_stop
    for key, value in trans_dict_own.items():
        if value[1] < real_start:   # If we have a earlier start coordinate than our longest read, keep that.
            real_start = value[1]
        if value[2] > real_stop:
            real_stop = value[2]

    if verbose:
        print('    > Initial length: '+ str(expected_stop-expected_start))
        if expected_start > real_start:
            print('    > Found better start coordinate, changing from: '+str(expected_start)+' to: '+str(real_start))
        else:
            print('    > No better start coordinate has been found, keeping: '+str(expected_start))
        if expected_stop < real_stop:
            print('    > Found better stop coordinate, changing from: ' +str(expected_stop)+ ' to: ' +str(real_stop))
        else:
            print('    > No better stop coordinate has been found, keeping: ' + str(expected_stop))
        print('    > New length: '+ str(real_stop-real_start))

    prime_extender_5 = expected_start - real_start
    prime_extender_3 = real_stop - expected_stop

    return backbone_trans_id, prime_extender_5, prime_extender_3, trans_id_list


def get_backbone_sequence(trans_id, prime_extender_5, prime_extender_3, gene_id, verbose=False):
    global server  # Add this as global variable since it won't change.

    ext = "/sequence/id/"+str(trans_id)+"?type=cdna;mask_feature=True;expand_5prime="+str(prime_extender_5)+";expand_3prime="+str(prime_extender_3)
    # mask=soft;
    # also: mask_feature=True
    # also "/sequence/id/"+finds[0]+"?expand=1;type=cdna;mask=soft"
    backbone_r = requests.get(server + ext, headers={"Content-Type": "text/x-fasta"})

    # TODO: check if this needs a better try except
    if verbose:
        print('  > Fetching backbone sequence from database...')

    if not backbone_r.status_code == 200:
        print('Error occurred whilst fetching url:\t', ext)
        #print('Adding to errornous genelist...')
        # raise Exception('Bad response')
        # TODO: find way to save 'bad' genes.
        return ''

    # Cuts the 5- & 3-prime UTR from the sequence, makes it lowercase and puts it back, mimicking the UTRs.
    fasta = backbone_r.text
    header, seq = fasta.split('\n', 1)
    seq = seq.upper()
    header = header+" backbone "+gene_id
    prime5 = seq[:prime_extender_5].lower()
    prime3 = seq[-prime_extender_3-1:].lower()
    cds = seq[prime_extender_5:-prime_extender_3-1]
    seq = prime5+cds+prime3
    header = header.replace(' ','|') # TODO: Fix this shouldn't be needed.
    fasta = header+'\n'+seq

    if verbose:
        print('    > Successfully fetched sequence!')
    return fasta


def get_sequence(trans_id_list, gene_id, verbose=False):
    is_coding = True
    global server  # Add this as global variable since it won't change.
    fasta_output = ''
    non_coding_num = 0

    for trans_id in trans_id_list:
        if verbose:
            print('  > Fetching transcript: '+trans_id+ ' from database...')

        ext = "/sequence/id/"+str(trans_id)+"?;type=cds"
        cds_r = requests.get(server + ext, headers={"Content-Type": "text/x-fasta"})
        cdna_r = cds_r  # Added to fix referenced before assigning error.

        if not cds_r.ok:
            if verbose:
                print("     > Seems like: "+trans_id+" doesn't have a known CDS, running checks...")

            # print('Error occurred whilst fetching url:\t', ext)
            # print('Retrying with cdna parameter')
            ext_cdna = "/sequence/id/" + str(trans_id) + "?;type=cdna"
            cdna_r = requests.get(server + ext_cdna, headers={"Content-Type": "text/x-fasta"})  # This is done to make sure the file didn't error out, but there really is no CSD #TODO explain better
            is_coding = False
            non_coding_num += 1
            if verbose:
                print('     > '+trans_id + " has no CDS confirmed, skipping...")
            continue

        if not cdna_r.ok:
            print('An error occurred, see the error code below.')
            cdna_r.raise_for_status()
            sys.exit()


        # Cuts the 5- & 3-prime UTR from the sequence, makes it lowercase and puts it back, mimicking the UTRs.
        fasta = cds_r.text
        header, seq = fasta.split('\n', 1)
        header = header+' '+gene_id
        header = header.replace(' ', '|')
        fasta = header+'\n'+seq

        if is_coding:  # Only add sequence to fasta when there's coding DNA.
            fasta_output = fasta_output + fasta

        # TODO Discuss if we want the headers for non CDS files.
        # NOTE: your in the canonical translation function
        # if is_non_coding:
        #     #       print(header+' NoCDS')
        #     fasta_output = fasta_output + '\n' + header + '|hasNoCDS'
        # elif not is_non_coding:
        #     fasta_output = fasta_output + '\n' + fasta


    print("  > Finished fetching transcript sequences.")
    print("  > Analytics:")
    print("    > Total number of transcripts: "+str(len(trans_id_list)) + '\t|\t' +
          "Number of transcripts without CDS: "+str(non_coding_num) + '\t|\t' +
          str(round(100*non_coding_num/len(trans_id_list), 2))+'%')

    return fasta_output


def cleanup_fasta(seq):
    """
    This function takes sequences that have gaps and weird stops in the sequence and removes them.
    These get joined together again in fasta format.
    :param (str) seq: Sequence that needs to be cleaned
    :return (str): Cleaned sequence in fasta format.
    """
    # Nice to have: See if this can be optimized using biopython
    seq = seq.replace(' ', '').replace('\n', '')  # TODO: Make it so it also works with single line strings.
    temp_list = [seq[index: index + 60] for index in range(0, len(seq), 60)]
    cleaned_seq = '\n'.join(temp_list)

    return cleaned_seq


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


# Note: This might need to be changed when adding the scripts to a pipeline
if __name__ == "__main__":
    main()
