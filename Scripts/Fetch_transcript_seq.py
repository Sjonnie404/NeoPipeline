########################################################################################
#  Script to fetch gene sequences based on coordinates
#  Input: Count file that contains Ensemble IDs
#  Output: Fasta file with headers, containing Ensemble information and NCBI sequence
#  Made by Shane Pullens, Utrecht University - Theoretical bioinformatics.
########################################################################################
# TODO: Add sys arguments
# Nice to have: GIU style?
# TODO: Add logs, use logging module
# Nice to have: Summary output, which could be used in the logs.

#  Imports
# from Bio import Entrez, SeqIO
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

#  Predefined variables
pp = pprint.PrettyPrinter(indent=4)  # Used to clearly print json files
server = "https://rest.ensembl.org"
verbose = True

def main():
    project_dir = Path.cwd()
    # path = "C:\\Users\\shane\\PycharmProjects\\NeoPipeline"
    # gene_scores = read_data(path+'\\Data\\Testing_ids.txt')
    gene_scores = read_data(Path(project_dir / 'Data' / 'Testing_ids.txt'))
    backbone_fasta_output = ''
    trans_fasta_output = ''

    for gene_id in gene_scores.keys():
        print('>>> Started process for:\t' + str(gene_id))
        gene_score = gene_scores.get(gene_id)
        backbone_trans_id, extend_5, extend_3, trans_id_list = get_Ensenble_data(gene_id, verbose)
        backbone_fasta_seq = get_backbone_sequence(backbone_trans_id, extend_5, extend_3, gene_id, gene_score, verbose)
        trans_fasta_seqs = get_sequence(trans_id_list, gene_id, gene_score, verbose)

        backbone_fasta_output = backbone_fasta_output +backbone_fasta_seq
        trans_fasta_output = trans_fasta_output +trans_fasta_seqs
        print('>>> Finished process for:\t' + str(gene_id))

    print('...')
    output_file_name = 'Testing_api'
    # write_file(backbone_fasta_output, output_file_name+'_backbone', path+'\\Output\\Fasta', True, True)
    # write_file(trans_fasta_output, output_file_name+'_transcripts', path+'\\Output\\Fasta', True, True)

    write_file(backbone_fasta_output, output_file_name+'_backbone', Path(project_dir / 'Output' / 'Fasta'), True, True)
    write_file(trans_fasta_output, output_file_name+'_transcripts', Path(project_dir / 'Output' / 'Fasta'), True, True)
    return None


def read_data(path):
    file = open(path, 'r')
    file = file.readlines()
    gene_score_dict = {}

    for line in file:
        gene_id, score = line.replace('\n', '').split('\t')
        gene_id = gene_id.split('.')[0]
        gene_score_dict[gene_id] = score

    return gene_score_dict


def get_Ensenble_data(gene_id='ENSG00000157764', verbose=False):
    global server  # Add this as global variable since it won't change.
    ext = "/lookup/id/"+str(gene_id)+"?expand=1"
    gene_r = requests.get(server + ext, headers={"Content-Type": "application/json"})

    if not gene_r.status_code == 200:
        print('Error occurred whilst fetching url:\t', ext)
        raise Exception('Bad response')

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


def get_backbone_sequence(trans_id, prime_extender_5, prime_extender_3, gene_id, htseq, verbose=False):
    global server  # Add this as global variable since it won't change.

    ext = "/sequence/id/"+str(trans_id)+"?mask_feature=True;expand_5prime="+str(prime_extender_5)+";expand_3prime="+str(prime_extender_3)
    # mask=soft;
    # also: mask_feature=True
    # also "/sequence/id/"+finds[0]+"?expand=1;type=cdna;mask=soft"
    backbone_r = requests.get(server + ext, headers={"Content-Type": "text/x-fasta"})

    # TODO: check if this needs a better try except
    if verbose:
        print('  > Fetching backbone sequence from database...')

    if not backbone_r.status_code == 200:
        print('Error occurred whilst fetching url:\t', ext)
        raise Exception('Bad response')

    #  Cuts the 5- & 3-prime from the sequence, makes it lowercase and puts it back, mimicking the introns.
    fasta = backbone_r.text
    header, seq = fasta.split('\n', 1)
    header = header+" backbone "+str(gene_id)+' ht_seq:'+str(htseq)
    prime5 = seq[:prime_extender_5].lower()
    prime3 = seq[-prime_extender_3-1:].lower()
    cds = seq[prime_extender_5:-prime_extender_3-1]
    seq = prime5+cds+prime3
    header = header.replace(' ','|') # TODO: Fix this shouldn't be needed.
    fasta = header+'\n'+seq

    if verbose:
        print('    > Successfully fetched sequence!')

    return fasta


def get_sequence(trans_id_list, gene_id, htseq, verbose=False):
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


        #  Cuts the 5- & 3-prime from the sequence, makes it lowercase and puts it back, mimicking the introns.
        fasta = cds_r.text
        header, seq = fasta.split('\n', 1)
        header = header+' '+str(gene_id)+' ht_seq:'+str(htseq)
        header = header.replace(' ', '|')
        fasta = header+'\n'+seq

        if is_coding:  # Only add sequence to fasta when there's coding DNA.
            fasta_output = fasta_output + fasta

        # TODO Discuss if we want the headers for non CDS files.
        # if is_non_coding:
        #     #       print(header+' NoCDS')
        #     fasta_output = fasta_output + '\n' + header + '|hasNoCDS'
        # elif not is_non_coding:
        #     fasta_output = fasta_output + '\n' + fasta

    if verbose:
        print("  > Finished fetching transcript sequences.")
        print("  > Analytics:")
        print("    > Total number of transcripts: "+str(len(trans_id_list)) + '\t|\t' +
              "Number of transcripts without CDS: "+str(non_coding_num) + '\t|\t' +
              str(round(100*non_coding_num/len(trans_id_list), 2))+'%')

    return fasta_output


def sequence_analysis(fasta):
    # DEPRECATED
    rm_lower = str.maketrans('', '', string.ascii_lowercase)
    rm_upper = str.maketrans('', '', string.ascii_uppercase)

    # TODO: Check what is the best way to characterize the UTRS, (now going for a added '|' before and after UTR)
    header, sequence = fasta.split('\n', 1)
    cds = sequence.translate(rm_lower).replace(' ', '').replace('\n', '')

    # Nice to have: find better solution for splitting and saving the split string.
    UTR_5, temp_seq = cds.replace('ATG', '1ATG', 1).split('1', 1)

    print(cds)
    print(len(cds))
    print('-'*80)
    # print(len(cds) % 3)
    testing = Seq(cds)

    # n = 3
    # split_strings = [cds[index: index + n] for index in range(0, len(cds), n)]
    # print(split_strings)
    # print(split_strings.index('ATG'))
    # print('-'*80)

    # print('Simple translation:')
    print(testing.translate())
    print(len(testing.translate()))
    print(testing.translate().find('M'))
    exit()
    print('cds translation:')
    print(len(testing.translate(cds=True)))
    print(testing.translate(cds=True))
    exit()

    # print(len(UTR_5))
    # print(len(UTR_5) % 3)
    print(UTR_5)
    print('-----')
    print(cleanup_fasta(UTR_5))

    temp_seq = temp_seq.translate(rm_lower).replace(' ', '').replace('\n', '')

    refurbished_seq = cleanup_fasta(temp_seq)
    print(len(refurbished_seq))
    print(len(refurbished_seq) % 3)

    # seq = N_parser(seq)
    # seq = seq.replace('\n', '',)

    # print('!!!!!!!!')
    # print(seq)
    # print(len(seq) % 3)
    # gene = Seq(seq)

    # print(len(gene))
    # print(gene.translate())
    # print(test)
    # try:
    #     print(gene.translate(cds=True))
    #     print('Found gene in:\t' + str(trans_id))
    #     print('adding to file...')
    # except:
    #     print('No gene found in:\t' + str(trans_id))

    # cds = seq[prime_extender_5:-prime_extender_3-1]
    # seq = prime5+cds+prime3
    # fasta = header+'\n'+seq

    return None


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



def write_file(text, filename, path=Path(), add_timestamp=True, overwrite_check=False):
    """
    # TODO Add documentation
    :param text:
    :param filename:
    :param path:
    :param add_timestamp:
    :param overwrite_check:
    :return:
    """
    if add_timestamp:  # Check if timestamp should be implemented
        time = datetime.now().strftime("%d-%b-%Y-h%H-m%M_")
        filename = time + filename  #+'.fasta' This was only here for testing right?

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
        print('>>> succesfully written '+filename)
    except:
        print('! An error occurred while trying to write '+filename)
    return None

# Note: This might need to be changed when adding the scripts to a pipeline
if __name__ == "__main__":
    main()