########################################################################################
#  Script to direct all functions in a streamlined, clear fasion.
#  Made by Shane Pullens, Utrecht University - Theoretical bioinformatics.
#  Version 1.0
########################################################################################

# Note: don't forget cleanup R script.
# TODO: add sys args for the script
# TODO: Major: clean up files in script for git push & check 'Note'
# TODO: Make sure all brakable things have nice exeptions
# TODO: Other parts could (also) be parallelized using multiple threads (like populating the downloadlist)

# Function imports
import Fetch_htseq_counts
import Get_top_counts
import Fetch_transcript_seq
import Get_translations
import Generic_functions
from Fetch_transcript_seq import threading_transcript_sequences
import MHCpan
from MHCpan import run_mhcPan_threading
from tqdm import tqdm

# Normal imports
import time
import pandas as pd
from pathlib import Path
from datetime import datetime
import os
from queue import Queue
from threading import Thread
from concurrent.futures import ThreadPoolExecutor
from concurrent import futures
from functools import partial


class ORFWorker(Thread):
    def __init__(self, queue):
        Thread.__init__(self)
        self.queue = queue

    def run(self):
        # Place for variables

        while True:
            gene_list, nthread = self.queue.get()
            try:
                Fetch_transcript_seq.threading_transcript_sequences(gene_list, nthread)
            finally:
                self.queue.task_done()

            return None



def main():
    """
    Main function to direct all functions down the road.
    :return: None
    """
    # TODO: Add argument parser
    # Argument parser
    # my_parser = argparse.ArgumentParser(description='Neopipe is a Neo-antigen detection pipeline.')

    # Predefined variables
    verbose = True
    # RNA_only = True
    checkpoints = True  # This variable states that all steps should be saved if something might break.
    n_threads = 60
    og_n_threads = n_threads
    project_dir = Path.cwd()
    timestamp = datetime.now().strftime("%d-%b-%Y-h%H-m%M")
    server = 'https://rest.ensembl.org'

    print('Welcome, and thank you for using the Neo Pipeline software v1.0')
    print('Please note: running this pipeline in a separate screen or session is highly advised!')
    print('Not sure about primary site and/or project ids?, please visit https://portal.gdc.cancer.gov/')
    print('Or leave the inputs blank for dummy data!')
    print("-"*80)
    time.sleep(1)

    # TEMP
    # project_dir = Path.cwd()
    # new_file_name = 'skin_TCGA-SKCM_20220822_132048.170180'
    # # new_file_name_cryptic_translated = new_file_name+'_cryptic_sequences_translated.fasta'
    # # new_file_name_canonical_translated = new_file_name+'_canonical_sequences_translated.fasta'
    # out_path = Path(project_dir / 'Output' / 'Counts' / new_file_name)
    #
    # cryptic_chunk_file_names = ''
    # canonical_chunk_file_names = ''
    # new_file_name_cryptic_translated = 'testing_cryptic_sequences_translated.fasta'
    # new_file_name_canonical_translated = 'testing_canonical_sequences_translated.fasta'
    # cryptic_chunk_file_names, cryptic_n_threads = MHCpan.mhcPan_thread_ripper(out_path,
    #                                                                           new_file_name_cryptic_translated, n_threads)   # 'rips' the target fasta file apart in chunks.
    # canonical_chunk_file_names, canonical_n_threads = MHCpan.mhcPan_thread_ripper(out_path,
    #                                                                               new_file_name_canonical_translated, n_threads)
    # combined_threads = cryptic_n_threads + canonical_n_threads
    # if combined_threads < n_threads:
    #     print('Amount of reads is lower than selected threads, lowering threads to ', str(combined_threads))
    #     n_threads = combined_threads
    # alleles = ['HLA-A01:01', 'HLA-A02:01']
    #
    # MHCpan_cryptic_thread_function = partial(run_mhcPan_threading, alleles, Path(out_path / 'Tmp'),
    #                                  canonical=False, safe_mode='AB')
    # MHCpan_canonical_thread_function = partial(run_mhcPan_threading, alleles, Path(out_path / 'Tmp'),
    #                                       canonical=True, safe_mode='AB')
    # jobs = []
    # with ThreadPoolExecutor(max_workers=n_threads) as executor:
    #     print('Fetching sequences over multiple threads...')
    #     for cryptic_filename in cryptic_chunk_file_names:
    #         # print('Working on filename:\t', filename)
    #         jobs.append(executor.submit(MHCpan_cryptic_thread_function, cryptic_filename))
    #     for canonical_filename in canonical_chunk_file_names:
    #         jobs.append(executor.submit(MHCpan_canonical_thread_function, canonical_filename))
    #
    # MHCpan.mhcPan_thread_assembly(Path(out_path / 'Tmp'), canonical=True)     # Adds all chunk outputs together in a single file.
    # MHCpan.mhcPan_thread_assembly(Path(out_path / 'Tmp'), canonical=False)
    # # MHCpan.remove_tmp_dir(Path(out_path / 'Tmp'))
    # print('Got to here')
    # exit()
    # #### END of TMP

    primary_site = input('Please specify the primary site:\t')
    # primary_site = 'skin' # TODO remove after debugging
    if primary_site == '':
        print('No primary site specified, using dummy site: "skin"')
        primary_site = 'skin'

    cancer_project = input('Please specify the project id: \t')
    if cancer_project == '':
        print('No cancer project specified, using dummy project: "TCGA-SKCM"')
        cancer_project = 'TCGA-SKCM'

    RNA_only = Generic_functions.str2bool(input('RNA only? (T/F)\t'))
    isLenient = Generic_functions.str2bool(input('Lenient? (True/False)\t'))
    info = "RNA only:\t" +str(RNA_only)+ "\n Lenient:\t" +str(isLenient)

    print('Selected variables:')
    print('RNA only:\t', str(RNA_only))
    print('is lenient:\t', str(isLenient))
    # cancer_project = 'EXCEPTIONAL_RESPONDERS-ER' # TODO remove after debugging
    if cancer_project == '':
        print('No cancer project specified, using dummy project: "EXCEPTIONAL_RESPONDERS-ER"')
        cancer_project = 'EXCEPTIONAL_RESPONDERS-ER'

    # ~~~ Fetching Star count data ~~~
    print('>>> Fetching star count data from GCC database...')
    user_input = {'primary_site': primary_site,
                  'project_id': cancer_project}

    file_name = Fetch_htseq_counts.getCountData(user_input, project_dir)
    new_file_name = Fetch_htseq_counts.extractFiles(file_name, project_dir, user_input, True)

    with open(Path(project_dir / 'Output' / 'Counts' / new_file_name / 'Session_info.txt'), 'w') as f:   # Debug, maybe add later as feature
        f.write(info)
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

    ntotal_genes = len(gene_list)
    print('>>> Amount of significant genes:\t'+str(ntotal_genes))
    # gene_scores = read_data(Path(project_dir / 'Counts' / ))

    backbone_fasta_output = ''
    trans_fasta_output = ''
    #TODO: REMOVE:

    if n_threads == 1:
        for i, gene_id in enumerate(gene_list):
            # perc = round((100 * i / total_genes), 2)
            # print(f' >> {perc}% - Working on gene {i} out of {total_genes}')
            # if verbose:
            #     print('  > Started process for:\t' + str(gene_id))

            backbone_trans_id, extend_5, extend_3, trans_id_list = Fetch_transcript_seq.get_Ensenble_data(gene_id, verbose=False)
            backbone_fasta_seq = Fetch_transcript_seq.get_backbone_sequence(backbone_trans_id, extend_5, extend_3, gene_id, verbose=False)
            trans_fasta_seqs = Fetch_transcript_seq.get_sequence(trans_id_list, gene_id, verbose=False)

            backbone_fasta_output = backbone_fasta_output + backbone_fasta_seq
            trans_fasta_output = trans_fasta_output + trans_fasta_seqs

    elif n_threads > 1:
        if n_threads > ntotal_genes:
            print('Workload is smaller than #threads, lowering #threads to ', str(len(gene_list)))
            n_threads = ntotal_genes
            gene_chunk = gene_list
        else:
            n_threads = 10
            gene_list_mod = ntotal_genes % n_threads
            gene_chunk_size = int((ntotal_genes - gene_list_mod) / n_threads)
            gene_chunks = [gene_list[x:x + gene_chunk_size] for x in range(0, ntotal_genes, gene_chunk_size)]
            gene_chunks[-2] = gene_chunks[-2] + gene_chunks[-1]
            del gene_chunks[-1]

            with ThreadPoolExecutor(max_workers=n_threads) as executor:
                jobs = []
                results_done = []

                print(len(gene_chunks))
                print('Fetching sequences over multiple threads...')
                for thread_gene_list in tqdm(gene_chunks):
                    # print(thread_gene_list)
                    time.sleep(3)
                    # print('Starting new thread...')
                    jobs.append(executor.submit(threading_transcript_sequences, thread_gene_list))

                    # executor.map(threading_transcript_sequences(thread_gene_list))
                print('All jobs successfully initialized!')

                backbone_fasta_output = ''
                trans_fasta_output = ''
                for job in futures.as_completed(jobs):
                    backbone_fasta_thread_output, trans_fasta_thread_output = job.result()

                    backbone_fasta_output = backbone_fasta_output + backbone_fasta_thread_output
                    trans_fasta_output = trans_fasta_output + trans_fasta_thread_output

    out_path = Path(project_dir / 'Output' / 'Counts' / new_file_name)
    if checkpoints:
        print('>>> Writing sequence files...')
        Generic_functions.write_file(backbone_fasta_output, new_file_name+'_backbone_sequences',
                                     out_path, timestamp, True, True)
        Generic_functions.write_file(trans_fasta_output, new_file_name+'_transcripts_sequences',
                                     out_path, timestamp, True, True)
        print('>>> Successfully wrote sequence files.')

    print('>>>  Translating sequences...')
    # NOTE: This might not be needed # TODO This should only be added when checkpoint flag is on
    # transcript_file = open(Path(project_dir / "Output" / test_dir / "Fasta" / filename_trans))  # relative path
    # transcript_file = transcript_file.read()
    # dir_path = Path(project_dir / "Output" / test_dir / "Fasta")
    # backbone_file = open(Path(project_dir / "Output" / test_dir / "Fasta" / filename_back))
    # backbone_file = backbone_file.read()
    transcript_file = trans_fasta_output  # debug
    backbone_file = backbone_fasta_output

    translated_canonical_file = Get_translations.canonical_translation(transcript_file)
    translated_cryptic_file = Get_translations.cryptic_translation(backbone_file, isLenient, verbose=False)

    print('>>> Writing translated transcription file...')
    write_path = Path(project_dir / 'Output' / 'Fasta')
    Generic_functions.write_file(translated_canonical_file, new_file_name+'_canonical_sequences_translated.fasta',
                                 out_path, timestamp, False, True)

    print('>>> Writing cryptic translated backbone file...')
    Generic_functions.write_file(translated_cryptic_file, new_file_name+'_cryptic_sequences_translated.fasta',
                                 out_path, timestamp, False, True)
    print('Succesfully wrote translated files...')


    # ~~~ Getting best binding peptides ~~~
    print('Scanning peptides for best binding...')
    new_file_name_cryptic_translated = new_file_name+'_cryptic_sequences_translated.fasta'
    new_file_name_canonical_translated = new_file_name+'_canonical_sequences_translated.fasta'

    alleles = ['HLA-A01:01', 'HLA-A02:01']  # These will be used by mhcPan
    n_threads = og_n_threads
    cryptic_chunk_file_names = ''

    if n_threads == 1:
        print('>>>Single thread selected, discarding multithreading...')
        new_file_name_canonical_translated = new_file_name+'_canonical_sequences_translated.fasta'
        new_file_name_cryptic_translated = new_file_name+'_cryptic_sequences_translated.fasta'
        MHCpan.run_mhcPan(alleles, out_path, new_file_name_canonical_translated,
                          canonical=True, safe_mode='AB')
        MHCpan.run_mhcPan(alleles, out_path, new_file_name_cryptic_translated,
                          canonical=False, safe_mode='AB')

    elif n_threads > 1:
        cryptic_chunk_file_names, cryptic_n_threads = MHCpan.mhcPan_thread_ripper(out_path,
                                                                              new_file_name_cryptic_translated, n_threads)   # 'rips' the target fasta file apart in chunks.
        canonical_chunk_file_names, canonical_n_threads = MHCpan.mhcPan_thread_ripper(out_path,
                                                                                  new_file_name_canonical_translated, n_threads)
        combined_threads = cryptic_n_threads + canonical_n_threads
        if combined_threads < n_threads:
            print('Amount of reads is lower than selected threads, lowering threads to ', str(combined_threads))
            n_threads = combined_threads
        alleles = ['HLA-A01:01', 'HLA-A02:01']

        MHCpan_cryptic_thread_function = partial(run_mhcPan_threading, alleles, Path(out_path / 'Tmp'),
                                                 canonical=False, safe_mode='AB')
        MHCpan_canonical_thread_function = partial(run_mhcPan_threading, alleles, Path(out_path / 'Tmp'),
                                                   canonical=True, safe_mode='AB')
        jobs = []
        with ThreadPoolExecutor(max_workers=n_threads) as executor:
            print('Fetching sequences over multiple threads...')
            for cryptic_filename in cryptic_chunk_file_names:
                # print('Working on filename:\t', filename)
                jobs.append(executor.submit(MHCpan_cryptic_thread_function, cryptic_filename))

            for canonical_filename in canonical_chunk_file_names:
                jobs.append(executor.submit(MHCpan_canonical_thread_function, canonical_filename))

        MHCpan.mhcPan_thread_assembly(Path(out_path / 'Tmp'), canonical=True)     # Adds all chunk outputs together in a single file.
        MHCpan.mhcPan_thread_assembly(Path(out_path / 'Tmp'), canonical=False)
        # MHCpan.remove_tmp_dir(Path(out_path / 'Tmp'))

    print('!!! Successfully finished the Neo-antigen detection pipeline !!!')
    exit()

    return None
# ----------------- Common workflow ----------------------

# X Fetch_htseq_counts.py                       # This fetches the star-count data from the GDC server.
# X Get_top_counts.py                             # Filters out the best genes.
# X Fetch_transcript_seq.py                       # Converts the gene name into canonical & cryptic DNA sequences
# X Get_translations.py                           # Converts the DNA sequences to Proteins (also does cryptic translations
# X MHCPAN                                        # converts protein sequences to small peptides and scores these
# MHCpan analytics & comparison to known peptides

if __name__ == "__main__":
    main()