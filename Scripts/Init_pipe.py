########################################################################################
#  Script to direct all functions in a streamlined, clear fashion.
#  Made by Shane Pullens, Utrecht University - Theoretical bioinformatics.
#  Version 1.0
########################################################################################
# Function imports
import Fetch_htseq_counts
import Get_top_counts
import Fetch_transcript_seq
import Get_translations
import Generic_functions
from Fetch_transcript_seq import threading_transcript_sequences
import Peptide_selection
import MHCpan
from MHCpan import run_mhcPan_threading

# Normal imports
import time
import pandas as pd
from pathlib import Path
from datetime import datetime
import argparse
import shutil
from concurrent.futures import ThreadPoolExecutor
from concurrent import futures
from functools import partial
from tqdm import tqdm


def main():
    """
    Main function to direct all functions in the correct order.
    :return: None
    """

    # Argument parser
    parser = argparse.ArgumentParser(description='Neo-antigen detection pipeline')
    parser.add_argument('-site', '--primary_site', default='skin', type=str, dest='site',
                        help='Define primary cancer site')
    parser.add_argument('-project', '--cancer_project', default='TCGA-SKCM', type=str, dest='project',
                        help='Define cancer project')
    parser.add_argument('-v', '--verbose', action='store_true', dest='verbose',
                        help='Enables verbose mode to show more output')
    parser.add_argument('-t', '--threads', default=4, type=int, dest='thread',
                        help='Defines number of threads to use for parallelization, default is 1')
    parser.add_argument('-check', '--checkpoints', action='store_true', dest='checkpoint',
                        help='Enables save mode to save in-between-step files')
    parser.add_argument('-time', '--add_time', action='store_true', dest='timestamp',
                        help='Adds a timestamp to the output folder')
    parser.add_argument('-o', '--output', type=str, dest='output', help='Defines the specified output folder name')
    parser.add_argument('-rna', '--rna_only', action='store_true', dest='rna',
                        help='Defines the use of RNA-only mode, this discarted all non RNA-related genes')
    parser.add_argument('-stop', '--to_stop', action='store_true', dest='stop',
                        help='Sets the paramater for translation. to_stop = True means a stop codon also needs to be '
                             'found for a complete translation')
    parser.add_argument('-SBthreshold', '--strong_binding_threshold', default=0.5, type=float, dest='SBthreshold',
                        help="Set the strong binding threshold, this defines the max rank a peptide can have to be "
                             "defined as 'Strong Binder'")
    parser.add_argument('-WBthreshold', '--weak_binding_threshold', default=2.0, type=float, dest='WBthreshold',
                        help="Set the weak binding threshold, this defines the max rank a peptide can have to be "
                             "defined as 'Weak Binder'")
    parser.add_argument('-cutoffp', '--peptide_cutoff_percentage', default=1.0, type=float, dest='percentage_cutoff',
                        help='Set the cutoff for peptide selection based on percentage, get overruled when '
                             'absolute cutoff is used')
    parser.add_argument('-cutoffa', '--peptide_cutoff_absolute', default=0, type=int, dest='absolute_cutoff',
                        help='Set the cutoff for peptide selection based on absolute numbers gets '
                             'overrules percentage cutoff')
    parser.add_argument('-PEPinc', '--peptide_inclusive', action='store_true', dest='peptide_inclusive',
                        help='Include the peptides that have the same ER rank, but fall off due to cutoffs')

    # command line parser arguments
    args = parser.parse_args()
    primary_site = args.site
    cancer_project = args.project
    verbose = args.verbose
    n_threads = args.thread
    checkpoints = args.checkpoint  # This variable states that all steps should be saved if something might break.
    add_timestamp = args.timestamp
    output_dir = args.output
    RNA_only = args.rna
    to_stop = args.stop
    SB_threshold = args.SBthreshold
    WB_threshold = args.WBthreshold
    percentage_cutoff = args.percentage_cutoff
    absolute_cutoff = args.absolute_cutoff
    peptide_inclusive = args.peptide_inclusive

    # Predefined variables
    og_n_threads = n_threads
    project_dir = Path.cwd()
    timestamp = datetime.now().strftime("%d-%b-%Y-h%H-m%M")
    server = 'https://rest.ensembl.org'
    cryptic_only = False
    if RNA_only == True:
        cryptic_only = True
    # Welcome message:
    print('Welcome, and thank you for using the Neo Pipeline software v1.0')
    print('Please note: running this pipeline in a separate screen or session is highly advised!')
    print('Not sure about primary site and/or project ids?, please visit https://portal.gdc.cancer.gov/')
    print('Or leave the inputs blank for a toy dataset!')
    print("-"*80)

    print('Selected variables:')
    print('RNA only:\t', str(RNA_only))
    print('To stop:\t', str(to_stop))
    print('Output name:\t', output_dir)
    print()
    session_info = f'This file contains the used parameters for this specific running session of this script.\nUsed arguments: {args}\nlTme this run was initiated: {time.time()}'

# ~~~ Fetching Star count data ~~~
    print('>>> Fetching star count data from GCC database...')
    user_input = {'primary_site': primary_site,
                  'project_id': cancer_project}

    file_name = Fetch_htseq_counts.getCountData(user_input, project_dir)
    new_file_name = Fetch_htseq_counts.extractFiles(file_name, project_dir, user_input, True)

    # Writes used parameters to file.
    with open(Path(project_dir / 'Output' / 'Counts' / new_file_name / 'Session_info.txt'), 'w') as f:
        f.write(session_info)

    # ~~~ Getting the top counts ~~~
    print('>>> Defining top Genes, based on count data...')
    input_paths = Path(project_dir / 'Output' / 'Counts' / new_file_name / 'Raw_counts').glob('*')
    Get_top_counts.getReadsMatrix(input_paths, RNA_only)
    Get_top_counts.runDeseq(Path(project_dir / 'Output' / 'Counts' / new_file_name))

    # ~~~ Fetching transcript sequences ~~~
    print('>>> Fetching transcript sequences from Ensembl datase...')
    # Opens the newly generated significant genes list
    genes_file = open(Path(project_dir / 'Output' / 'Counts' / new_file_name / 'significant_genes.csv'), 'r')
    genes_file = genes_file.readlines()
    gene_list = []

    for line in genes_file:
        gene_id = line.replace('"','').replace("\n","").split('.')[0]
        gene_list.append(gene_id)

    ntotal_genes = len(gene_list)
    print('>>> Amount of significant genes:\t'+str(ntotal_genes))

    backbone_fasta_output = ''
    trans_fasta_output = ''

    # This is ran when a single thread is used
    if n_threads == 1:
        for i, gene_id in enumerate(gene_list):
            backbone_trans_id, extend_5, extend_3, trans_id_list = Fetch_transcript_seq.get_Ensemble_data(gene_id, verbose=False)
            backbone_fasta_seq = Fetch_transcript_seq.get_backbone_sequence(backbone_trans_id, extend_5, extend_3, gene_id, verbose=False)
            trans_fasta_seqs = Fetch_transcript_seq.get_sequence(trans_id_list, gene_id, verbose=False)

            # Add output from all genes together in dedicated fasta file
            backbone_fasta_output = backbone_fasta_output + backbone_fasta_seq
            trans_fasta_output = trans_fasta_output + trans_fasta_seqs

    # This is ran when multiple threads are used
    elif n_threads > 1:
        if n_threads > ntotal_genes:
            print('Workload is smaller than #threads, lowering #threads to ', str(len(gene_list)))
            n_threads = ntotal_genes
            gene_chunk = gene_list
        else:
            n_threads = 10
            # Here, we splice up the gene list and saving it in new files.
            # This is done because netMHCpan only takes an input file, and no input strings.
            gene_list_mod = ntotal_genes % n_threads
            gene_chunk_size = int((ntotal_genes - gene_list_mod) / n_threads)
            gene_chunks = [gene_list[x:x + gene_chunk_size] for x in range(0, ntotal_genes, gene_chunk_size)]
            # Here we add the few remaining genes to the last file.
            gene_chunks[-2] = gene_chunks[-2] + gene_chunks[-1]
            del gene_chunks[-1]

            # Initialize parallelization
            with ThreadPoolExecutor(max_workers=n_threads) as executor:
                jobs = []
                results_done = []
                print('Fetching sequences over multiple threads...')
                for thread_gene_list in tqdm(gene_chunks):
                    # wait period is added to reduce Ensemble server traffic.
                    time.sleep(3)
                    jobs.append(executor.submit(threading_transcript_sequences, thread_gene_list))
                print('All jobs successfully initialized!')

                backbone_fasta_output = ''
                trans_fasta_output = ''
                for job in futures.as_completed(jobs):
                    # Extracts outputs from threads and stitches them together
                    backbone_fasta_thread_output, trans_fasta_thread_output = job.result()
                    backbone_fasta_output = backbone_fasta_output + backbone_fasta_thread_output
                    trans_fasta_output = trans_fasta_output + trans_fasta_thread_output

    out_path = Path(project_dir / 'Output' / 'Counts' / new_file_name)
    # When checkpoint boolean is enabled, save these fasta files.
    if checkpoints:
        print('>>> Writing sequence files...')
        Generic_functions.write_file(backbone_fasta_output, new_file_name+'_backbone_sequences',
                                     out_path, timestamp, True, True)
        Generic_functions.write_file(trans_fasta_output, new_file_name+'_transcripts_sequences',
                                     out_path, timestamp, True, True)
        print('>>> Successfully wrote sequence files.')

    print('>>>  Translating sequences...')
    translated_canonical_file = Get_translations.canonical_translation(trans_fasta_output)
    translated_cryptic_file = Get_translations.cryptic_translation(backbone_fasta_output, to_stop, verbose=True)

    print('>>> Writing translated transcription file...')
    Generic_functions.write_file(translated_canonical_file, str(new_file_name)+'_canonical_sequences_translated.fasta',
                                 out_path, timestamp, False, True)

    print('>>> Writing cryptic translated backbone file...')
    Generic_functions.write_file(translated_cryptic_file, str(new_file_name)+'_cryptic_sequences_translated.fasta',
                                 out_path, timestamp, False, True)
    print('Succesfully wrote translated files...')


    # ~~~ Getting best binding peptides ~~~
    print('Scanning peptides for best binding...')
    new_file_name_cryptic_translated = new_file_name+'_cryptic_sequences_translated.fasta'
    new_file_name_canonical_translated = new_file_name+'_canonical_sequences_translated.fasta'

    alleles = ['HLA-A01:01', 'HLA-A02:01']  # These will be used by mhcPan
    n_threads = og_n_threads

    # variable assignment to suppress 'possible referenced before assignment' warning
    MHCpan_cryp_df = MHCpan_can_df = pd.DataFrame()
    cryp_hla_a01_df = can_hla_a01_df = cryp_hla_a02_df = can_hla_a02_df = pd.DataFrame()
    can_hla_a01_output_peptides = can_hla_a02_output_peptides = pd.DataFrame()

    if n_threads == 1:
        print('>>>Single thread selected, discarding multithreading...')
        new_file_name_canonical_translated = new_file_name+'_canonical_sequences_translated.fasta'
        new_file_name_cryptic_translated = new_file_name+'_cryptic_sequences_translated.fasta'
        MHCpan.run_mhcPan(alleles, out_path, new_file_name_canonical_translated,
                          canonical=True, safe_mode='AB')
        MHCpan.run_mhcPan(alleles, out_path, new_file_name_cryptic_translated,
                          canonical=False, safe_mode='AB')

    # When there's multiple threads, the translated fasta files get ripped apart, with an even distribution of reads
    # This is done because local netMHCpan only takes files as input.
    elif n_threads > 1:
        # Rips files apart for both canonical and cryptic
        cryptic_chunk_file_names, cryptic_n_threads = MHCpan.mhcPan_thread_ripper(out_path,
                                                                              new_file_name_cryptic_translated, n_threads)   # 'rips' the target fasta file apart in chunks.
        canonical_chunk_file_names, canonical_n_threads = MHCpan.mhcPan_thread_ripper(out_path,
                                                                                  new_file_name_canonical_translated, n_threads)
        # Since we can't rip a single read apart, it could happen that we have excess threads, if so discard them.
        combined_threads = cryptic_n_threads + canonical_n_threads
        if combined_threads < n_threads:
            print('Amount of reads is lower than selected threads, lowering threads to ', str(combined_threads))
            n_threads = combined_threads

        alleles = ['HLA-A01:01', 'HLA-A02:01']
        # Define threaded MHCpan functions for both canonical and cryptic
        MHCpan_cryptic_thread_function = partial(run_mhcPan_threading, alleles, Path(out_path / 'Tmp'),
                                                 canonical=False, safe_mode='AB', SB_threshold=SB_threshold, WB_threshold=WB_threshold)
        MHCpan_canonical_thread_function = partial(run_mhcPan_threading, alleles, Path(out_path / 'Tmp'),
                                                   canonical=True, safe_mode='AB', SB_threshold=SB_threshold, WB_threshold=WB_threshold)
        jobs = []
        with ThreadPoolExecutor(max_workers=n_threads) as executor:
            print('Fetching sequences over multiple threads...')
            # Here we start all thread jobs, for predefined cryptic file chunks
            for cryptic_filename in cryptic_chunk_file_names:
                jobs.append(executor.submit(MHCpan_cryptic_thread_function, cryptic_filename))
            # Here we start all thread jobs, for predefined cryptic file chunks
            for canonical_filename in canonical_chunk_file_names:
                jobs.append(executor.submit(MHCpan_canonical_thread_function, canonical_filename))

        print('Actually finished all threads.')
        # When all threads are finished, the multiple MHCpan outputs get stiched back together into a single file
        MHCpan_cryp_df = MHCpan.mhcPan_thread_assembly(Path(out_path / 'Tmp'), canonical=False)
        if not cryptic_only:
            MHCpan_can_df = MHCpan.mhcPan_thread_assembly(Path(out_path / 'Tmp'), canonical=True)     # Adds all chunk outputs together in a single file.
        # MHCpan.remove_tmp_dir(Path(out_path / 'Tmp'))

    if MHCpan_cryp_df.empty:
        print('No binders have been found with the given parameters for cryptic peptides!')
    else:
        # If the dataframe is not empthy, and cryptic only is true, subsample the HLA's from the output file.
        try:
            cryp_hla_a01_df = MHCpan_cryp_df[MHCpan_cryp_df['MHC'] == 'HLA-A*01:01']
            cryp_hla_a02_df = MHCpan_cryp_df[MHCpan_cryp_df['MHC'] == 'HLA-A*02:01']
        except:
            print('Something went wrong with subsampling HLAs from cryptic dataframe!')
            exit()

    if not cryptic_only:
        # When we also want to look for canonical peptides (cryptic_only = False), check if the dataframe is not empty
        # And subsample the HLA's from the output file.
        if MHCpan_can_df.empty:
            print('No binders have been found with the given parameters for canonical peptides!')
        else:
            try:
                can_hla_a01_df = MHCpan_can_df[MHCpan_can_df['MHC'] == 'HLA-A*01:01']
                can_hla_a02_df = MHCpan_can_df[MHCpan_can_df['MHC'] == 'HLA-A*02:01']
            except:
                print('Something went wrong with subsampling HLAs from canonical dataframe!')
                exit()

    print('Starting HLA1...')
    # Here, we first select the top x (defined by user) peptides, than they are compared with known peptides.
    # This will be done for both alleles, and also canonical if defined by the user.
    cryp_hla_a01_top_peptides = Peptide_selection.peptide_filtering(cryp_hla_a01_df, percentage_cutoff, absolute_cutoff, inclusive=peptide_inclusive)
    cryp_hla_a01_output_peptides, cryp_hla1_FP = Peptide_selection.peptide_comparison(cryp_hla_a01_top_peptides, hla='HLA-A01')

    if not cryptic_only:
        can_hla_a01_top_peptides = Peptide_selection.peptide_filtering(can_hla_a01_df, percentage_cutoff, absolute_cutoff, inclusive=peptide_inclusive)
        can_hla_a01_output_peptides, can_hla1_FP = Peptide_selection.peptide_comparison(can_hla_a01_top_peptides, hla='HLA-A01')

    print('Starting HLA2...')
    cryp_hla_a02_top_peptides = Peptide_selection.peptide_filtering(cryp_hla_a02_df, percentage_cutoff, absolute_cutoff, inclusive=peptide_inclusive)
    cryp_hla_a02_output_peptides, cryp_hla2_FP = Peptide_selection.peptide_comparison(cryp_hla_a02_top_peptides, hla='HLA-A02')

    if not cryptic_only:
        can_hla_a02_top_peptides = Peptide_selection.peptide_filtering(can_hla_a02_df, percentage_cutoff, absolute_cutoff, inclusive=peptide_inclusive)
        can_hla_a02_output_peptides, can_hla2_FP = Peptide_selection.peptide_comparison(can_hla_a02_top_peptides, hla='HLA-A02')

    print('Saving best candidates...')
    # Here, all peptide candidates and false positives are saved.
    cryp_hla_a01_output_peptides.to_csv(Path(out_path / 'Candidate_HLA-A01-cryptic_peptides.csv'))
    cryp_hla_a02_output_peptides.to_csv(Path(out_path / 'Candidate_HLA-A02-cryptic_peptides.csv'))
    if not cryptic_only:
        can_hla_a01_output_peptides.to_csv(Path(out_path / 'Candidate_HLA-A01-canonical_peptides.csv'))
        can_hla_a02_output_peptides.to_csv(Path(out_path / 'Candidate_HLA-A02-canonical_peptides.csv'))

    print('saving False positives')
    try:
        cryp_hla1_FP.to_csv(Path(out_path / 'hla1_cryp_FP.csv'))
        cryp_hla2_FP.to_csv(Path(out_path / 'hla2_cryp_FP.csv'))
    except:
        print('Could not save cryptic false positives!')
    if not cryptic_only:
        try:
            can_hla1_FP.to_csv(Path(out_path / 'hla1_can_FP.csv'))
            can_hla2_FP.to_csv(Path(out_path / 'hla2_can_FP.csv'))
        except:
            print('Could not save canonical false positives!')

    if output_dir != '':  # Renames the output folder
        shutil.move(out_path, Path(out_path.parent / output_dir))
    print('!!! Successfully finished the Neo-antigen detection pipeline !!!')
    return None


if __name__ == "__main__":
    main()