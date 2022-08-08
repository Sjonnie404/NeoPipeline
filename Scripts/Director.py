########################################################################################
#  Script to direct all functions in a streamlined, clear fasion.
#  Made by Shane Pullens, Utrecht University - Theoretical bioinformatics.
#  Version 1.0
########################################################################################

# TODO: Major: clean up files in script for git push & check 'Note'

# Function imports
import Fetch_htseq_counts
import Get_top_counts
import Fetch_transcript_seq
import Get_translations
import Generic_functions

# Normal imports
import time
import pandas as pd
from pathlib import Path
from datetime import datetime

# NOTE think about write file, readfiile etc, etc.

def main():
    """

    :return:
    """
    # TODO: Add nice message formatting, like >>>
    # TODO: Maybe we don't need seperate output files, but all in own output file
    # Predefined variables
    verbose = True
    RNA_only = True
    checkpoints = True  # This variable states that all steps should be saved if something might break.
    project_dir = Path.cwd()
    timestamp = datetime.now().strftime("%d-%b-%Y-h%H-m%M")
    server = 'https://rest.ensembl.org'

    print('Welcome, and thank you for using the Neo Pipeline software v1.0')
    print('Please note: running this pipeline in a separate screen or session is highly advised!')

    # primary_site = input('Please specify the primary site:\t')
    primary_site = 'skin' # TODO remove after debugging
    if primary_site == '':
        print('No primary site specified, using dummy site: "skin"')
        primary_site = 'skin'

    # cancer_project = input('Please specify the project id: \t')
    cancer_project = 'EXCEPTIONAL_RESPONDERS-ER' # TODO remove after debugging
    if cancer_project == '':
        print('No cancer project specified, using dummy project: "EXCEPTIONAL_RESPONDERS-ER"')
        cancer_project = 'EXCEPTIONAL_RESPONDERS-ER'

    # ~~~ Fetching Star count data ~~~
    print('>>> Fetching star count data from GCC database...')
    user_input = {'primary_site': primary_site,
                  'project_id': cancer_project}

    file_name = Fetch_htseq_counts.getCountData(user_input, project_dir)
    new_file_name = Fetch_htseq_counts.extractFiles(file_name, project_dir, user_input, True)

    # ~~~ Getting the top counts ~~~
    print('>>> Defining top Genes, based on count data...')
    input_paths = Path(project_dir / 'Output' / 'Counts' / new_file_name / 'Raw_counts').glob('*')
    Get_top_counts.getReadsMatrix(input_paths, RNA_only, True)  # Note: save matrix should be True for Deseq
    Get_top_counts.runDeseq(Path(project_dir / 'Output' / 'Counts' / new_file_name))

    # ~~~ Fetching transcript sequences ~~~
    print('>>> Fetching transcript sequences from Ensembl datase...')
    print('\t (please note, this step could take a while!)')

    genes_file = open(Path(project_dir / 'Output' / 'Counts' / new_file_name / 'significant_genes.csv'), 'r')
    genes_file = genes_file.readlines()
    gene_list = []
    for line in genes_file:
        # print(line.replace('"','').replace("\n","").split('.')[0])
        gene_id = line.replace('"','').replace("\n","").split('.')[0]
        gene_list.append(gene_id)

    total_genes = len(gene_list)
    print('>>> Amount of significant genes:\t'+str(total_genes))
    # gene_scores = read_data(Path(project_dir / 'Counts' / ))
    backbone_fasta_output = ''
    trans_fasta_output = ''

    for i, gene_id in enumerate(gene_list):
        perc = round((100 * i / total_genes), 2)
        print(f' >> {perc}% - Working on gene {i} out of {total_genes}')
        if verbose:
            print('  > Started process for:\t' + str(gene_id))

        backbone_trans_id, extend_5, extend_3, trans_id_list = Fetch_transcript_seq.get_Ensenble_data(gene_id, verbose=False)
        backbone_fasta_seq = Fetch_transcript_seq.get_backbone_sequence(backbone_trans_id, extend_5, extend_3, gene_id, verbose=False)
        trans_fasta_seqs = Fetch_transcript_seq.get_sequence(trans_id_list, gene_id, verbose=False)

        backbone_fasta_output = backbone_fasta_output +backbone_fasta_seq
        trans_fasta_output = trans_fasta_output +trans_fasta_seqs
        if verbose:
            print('  > Finished process for:\t' + str(gene_id))

    out_path = Path(project_dir / 'Output' / 'Counts' / new_file_name)
    if checkpoints:
        print('>>> Writing sequence files...')
        Generic_functions.write_file(backbone_fasta_output, new_file_name+'_backbone_sequences', out_path, timestamp, True, True)
        Generic_functions.write_file(trans_fasta_output, new_file_name+'_transcripts_sequences', out_path, timestamp, True, True)
        print('>>> Successfully wrote sequence files.')
    exit()

    # ~~~ Translating sequences ~~~
    # NOTE: This might not be needed
    # transcript_file = open(Path(project_dir / "Output" / test_dir / "Fasta" / filename_trans))  # relative path
    # transcript_file = transcript_file.read()
    # dir_path = Path(project_dir / "Output" / test_dir / "Fasta")
    # backbone_file = open(Path(project_dir / "Output" / test_dir / "Fasta" / filename_back))
    # backbone_file = backbone_file.read()
    transcript_file = trans_fasta_output
    backbone_file = backbone_fasta_output

    translated_trans_file = Get_translations.canonical_translation(transcript_file)
    # translated_backbone_file = Get_translations.cryptic_translation(backbone_file)

    print('>>> Writing translated transcription file...')
    write_path = Path(project_dir / 'Output' / 'Fasta')
    Generic_functions.write_file(translated_trans_file, new_file_name+'_translated.fasta', write_path, timestamp, False, True)

    # Note: important
    # TODO: Find cryptic translation library
    # print('>>> Writing cryptic translated backbone file...')
    # write_file(translated_backbone_file, new_file_name+'_translated.fasta', write_path, timestamp, False, True) # Note needs to be converted to backbone file.
    # translated_backbone_file.write()
    print('Succesfully wrote files...')


    return None
# ----------------- Common workflow ----------------------

# X Fetch_htseq_counts.py                       # This fetches the star-count data from the GDC server.
# X Get_top_counts.py                             # Filters out the best genes.
# X Fetch_transcript_seq.py                       # Converts the gene name into canonical & cryptic DNA sequences
# X Get_translations.py                           # Converts the DNA sequences to Proteins (also does cryptic translations
# MHCPAN                                        # converts protein sequences to small peptides and scores these


if __name__ == "__main__":
    main()