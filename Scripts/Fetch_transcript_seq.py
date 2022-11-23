########################################################################################
#  Script to fetch gene sequences based on coordinates
#  Input: Count file that contains Ensemble IDs
#  Output: Fasta file with headers, containing Ensemble information and NCBI sequence
#  Made by Shane Pullens, Utrecht University - Theoretical bioinformatics.
########################################################################################

#  Imports
import time
import pprint
import requests
import os.path
from pathlib import Path

#  Predefined variables
pp = pprint.PrettyPrinter(indent=4)  # Used to clearly print json files
server = "https://rest.ensembl.org"
verbose = False


def main():
    return None


def threading_transcript_sequences(gene_list, verbose=False):
    backbone_fasta_output = ''
    trans_fasta_output = ''

    for i, gene_id in enumerate(gene_list):
        backbone_trans_id, extend_5, extend_3, trans_id_list = get_Ensemble_data(gene_id, verbose=False)
        backbone_fasta_seq = get_backbone_sequence(backbone_trans_id, extend_5, extend_3, gene_id, verbose=False)
        trans_fasta_seqs = get_sequence(trans_id_list, gene_id, verbose=False)

        backbone_fasta_output = backbone_fasta_output + backbone_fasta_seq
        trans_fasta_output = trans_fasta_output + trans_fasta_seqs

    return backbone_fasta_output, trans_fasta_output


def read_data(path):
    file = open(path, 'r')
    file = file.readlines()
    gene_score_dict = {}

    for line in file[6:]:  # Skip the first 6 header lines
        gene_id, gene_name, gene_type, score, _, _, _, _, _, = line.replace('\n', '').split('\t')
        gene_id = gene_id.split('.')[0]
        gene_score_dict[gene_id] = [score, gene_name, gene_type]

    return gene_score_dict


def get_Ensemble_data(gene_id='ENSG00000157764', verbose=False):
    """
    Here we fetch the transcript sequences for all gene names using the Ensemble database
    :param gene_id: string containing the gene id
    :param verbose: boolean for verbose mode
    :return: backbone_trans_id: gives the ID for the transcript that has the highest overlap with the gene.
    :return: prime_extender_5: gives the amount of based that the start should be moved to cover the whole gene.
    :return: prime_extender_3: gives the amount of based that the stop should be moved to cover the whole gene.
    :return: trans_id_list: returns a list of transcript IDs
    """
    global server  # Add this as global variable since it won't change.
    # Web-scraping magic
    ext = "/lookup/id/"+str(gene_id)+"?expand=1"
    gene_r = requests.get(server + ext, headers={"Content-Type": "application/json"})

    # If the API doesn't return a 200 status code, we now something broke on Ensemble's site.
    if not gene_r.status_code == 200:
        if verbose:
            print('[data] Error occurred whilst fetching url:\t', ext)
        #print('Adding to errornous genelist...')
        # raise Exception('Bad response')
        return '', 1, 1, ['']

    # Extract data from fetched Json file.
    decoded = gene_r.json()
    trans_dict = decoded.get('Transcript')
    trans_dict_own = {}

    if verbose:
        print('  > Found '+str(len(trans_dict))+' transcripts for '+str(gene_id)+'!')

    # Extracts data dict for each transcript key.
    for trans in trans_dict:
        start = trans.get('start')
        end = trans.get('end')
        id = trans.get('id')
        length = end-start
        trans_dict_own.update({id: [length, start, end]})

    trans_id_list = list(trans_dict_own.keys())
    # Sorts the id list in descending order
    trans_list_des = sorted(trans_dict_own.items(), key=lambda x: x[1], reverse=True)
    backbone_trans_id = trans_list_des[0][0]  # Extracts actual Id from dict list.
    # Extracts the start and stop location of selected transcript.
    expected_start = trans_dict_own.get(backbone_trans_id)[1]
    expected_stop = trans_dict_own.get(backbone_trans_id)[2]

    if verbose:
        print('  > Selected backbode ID: '+str(backbone_trans_id))
        print('  > Initial start location: '+str(expected_start))
        print('  > Initial stop location: ' + str(expected_stop))
        print('  > Checking for better start&stop locations...')

    # Here we check if we find a transcript (from the same gene) with an earlier or later stop codon.
    # This is done to see if we can get better coverage of the gene.
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

    # Calculations for better start- and stop coordinates.
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
        if verbose:
            print('[backbone] Error occurred whilst fetching url:\t', ext)
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
    header = header.replace(' ','|')
    fasta = header+'\n'+seq

    if verbose:
        print('    > Successfully fetched sequence!')
    # print('Crypt - Success') #REMOVE
    return fasta


def get_sequence(trans_id_list, gene_id, verbose=False, max_tries=2, delay=5):
    """
    Extracts all sequences from found transcripts that contain a known coding region.
    :param trans_id_list: list of transcript IDs
    :param gene_id: Gene ID that's currently selected.
    :param verbose: Boolean for verbose mode.
    :param max_tries: Maximum amount of tries for server call after error.
    :param delay: Delay in seconds after retry to reduce server traffic.
    :return: fasta list of all found sequences.
    """
    is_coding = True
    global server  # Add this as global variable since it won't change.
    fasta_output = ''
    non_coding_num = 0

    for trans_id in trans_id_list:
        if verbose:
            print('  > Fetching transcript: '+trans_id+' from database...')

        ext = "/sequence/id/"+str(trans_id)+"?;type=cds"
        cds_r = requests.get(server + ext, headers={"Content-Type": "text/x-fasta"})
        cdna_r = cds_r  # Added to fix referenced before assigning error.

        if not cds_r.ok:    # if we get an error on gathering the CDS, we need to check if the error originates from
            # there not being a CDS to begin with, or another error all along.
            # To tackle this, when this error occurs we check if can download the dna without the CDS.
            # If this is the case, we know the initial error originates from there not being a CDS and nothing else.

            if verbose:
                print("     > Seems like: "+trans_id+" doesn't have a known CDS, running checks...")
            ext_cdna = "/sequence/id/" + str(trans_id) + "?;type=cdna"

            tries = 0
            while True:
                tries += 1
                if tries <= max_tries:
                    try:
                        cdna_r = requests.get(server + ext_cdna, headers={"Content-Type": "text/x-fasta"})
                    except requests.exceptions.ConnectionError:
                        print(f'Received connection timeout error, waiting {delay} seconds and trying again...')
                        print(f'This was try {tries}, maximum number of tries is set to {max_tries}.')
                        time.sleep(delay)
                    break

            is_coding = False
            non_coding_num += 1
            if verbose:
                print('     > '+trans_id + " has no CDS confirmed, skipping...")
            continue

        if not cdna_r.ok:
            if verbose:
                print('[Transcript] An error occurred, see the error code below.')
            cdna_r.raise_for_status()
            return ''

        # Cuts the 5- & 3-prime UTR from the sequence, makes it lowercase and puts it back, mimicking the UTRs.
        fasta = cds_r.text
        header, seq = fasta.split('\n', 1)
        header = header+' '+gene_id
        header = header.replace(' ', '|')
        fasta = header+'\n'+seq7

        if is_coding:  # Only add sequence to fasta when there's coding DNA.
            fasta_output = fasta_output + fasta

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
