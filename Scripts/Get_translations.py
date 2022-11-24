########################################################################################
#  Script translate backbone and transcript fasta non-canonical & canonical respectfully
#  Input: Backbone fasta & Transscript fasta.
#  Output: Fasta files with translated backbone en transcript sequences
#  Made by Shane Pullens, Utrecht University - Theoretical bioinformatics.
########################################################################################

# Imports
import numpy as np
from numpy import ma
from Bio.Seq import Seq
import itertools as iter
from difflib import SequenceMatcher as SeqMatch
import re


def main():
    return None


def cleanup_fasta(seq):
    """
    This function takes sequences that have gaps and weird stops in the sequence and removes them.
    These get joined together again in fasta format.
    :param (str) seq: Sequence that needs to be cleaned
    :return (str): Cleaned sequence in fasta format.
    """
    seq = seq.replace(' ', '').replace('\n', '')
    temp_list = [seq[index: index + 60] for index in range(0, len(seq), 60)]
    cleaned_seq = '\n'.join(temp_list)

    return cleaned_seq


def canonical_translation(file):
    """
    Method translates the DNA sequences to Amino Acid sequences using the Biopython translation module.
    :param file: Fasta file containing DNA reads
    :return: Fasta file containing Amino acid reads
    """
    trans_file = ''
    # Here we split on '>', without removing the '>'
    fastas = file.replace('>', '$$>').split('$$')[1:]
    for fasta in fastas:
        # For each read, the header is isolated, the sequence is checked on missing nucleotides and is fully translated.
        header, seq = fasta.split('\n', 1)
        seq = N_parser(seq.replace('\n', ''))  # Adds trailing N, for sequences that aren't 3 fold
        seq = Seq(seq)  # Generated Seq object for biopython
        # CDS is set to False here, due to Ensemble already removes the stop codon.
        translated_seq = str(seq.translate(cds=False))
        translated_seq = cleanup_fasta(translated_seq)
        translated_fasta = header + '\n' + translated_seq

        trans_file = trans_file + translated_fasta + '\n'  # Adds reads together in fasta format
    return trans_file


def cryptic_translation(file, to_stop=True, verbose=False):
    """
    Method to translate the reads in a cryptic manner (see paper).
    Alternative translation is done based on alternate Open Reading Frames (aORFs)
    This is highly due to change due to it being brand-new
    ORF naming:
    ORF_1:  normal, forward strand
    ORF_2:  forward strand +1
    ORF_3:  forward strand +2
    ORF_4:  normal, reverse strand
    ORF_5:  reverse strand -1
    ORF_6:  reverse strand -2

    :param file: Fasta file containing the reads
    :param to_stop: boolean to make stop codon mandatory for adding cryptic read to list
    :param verbose: boolean to enable verbose mode
    :return: fasta file with translated reads.
    """

    trans_file = ''
    fastas = file.replace('>', '$$>').split('$$')[1:]
    ORFs = ['ORF_1', 'ORF_2', 'ORF_3', 'ORF_4', 'ORF_5', 'ORF_6']
    alt_ORFs = {}

    for fasta in fastas:
        base_header, base_seq = fasta.split('\n', 1)
        translated_ORFs = ''

        # Checks are performed for each ORF
        for i, ORF in enumerate(ORFs):
            if verbose:
                print('> Started checks for:\t', ORF)
            # This checks if the strand should be reversed or not
            if i <= 2:
                seq = base_seq
            elif i >= 3:
                seq = base_seq[::-1]  # Reverses the sequence

            # This checks what reading frame should be used and parses it to altORF_finder
            # altORF_finder searches for alternative ORFs to use.
            if i == 0 or i == 3:
                seq = N_parser(seq.replace('\n', ''))
                alt_ORFs = altORF_finder(seq, to_stop, False)
            elif i == 1 or i == 4:
                seq = N_parser(seq.replace('\n', '')[1:])
                alt_ORFs = altORF_finder(seq, to_stop, False)
            elif i == 2 or i == 5:
                seq = N_parser(seq.replace('\n', '')[2:])
                alt_ORFs = altORF_finder(seq, to_stop, False)

            # Checks the amount of ORFs, if there's more than 1 ORF, do conflate check.
            if len(alt_ORFs.keys()) <= 1:
                a = 1
            elif len(alt_ORFs.keys()) > 1:
                conflated_ORFs = conflate_ORFs(alt_ORFs)
                alt_ORFs = conflated_ORFs

            # For each found ORF, translate the coordinates to Amino acids and safe these
            for alt_ORF_seq, alt_ORF_meta in alt_ORFs.items():
                translated_ORF_seq = str(Seq(alt_ORF_seq).translate(to_stop=False, stop_symbol='N', cds=False, gap='n'))
                ORF_header = base_header + '|' + ORF + '|' + alt_ORF_meta[0] + '|' + alt_ORF_meta[1]
                translated_ORF_seq = cleanup_fasta(translated_ORF_seq)
                translated_ORF = ORF_header + '\n' + translated_ORF_seq
                translated_ORFs = translated_ORFs + '\n' + translated_ORF
        trans_file = trans_file + translated_ORFs

    return trans_file


def conflate_ORFs(ORF_dict, threshold=0.9):
    """
    This function checks if there's overlap in the found ORFs, if found, they are conflated.
    :param ORF_dict:
    :param threshold: The threshold on which overlapping ratio the ORFs should be inflated
    :return: returns a dictionary with the same information, with conflated ORFs
    """
    for a, b in iter.combinations(list(ORF_dict.keys()), 2):
        ratio = SeqMatch(None, a, b, autojunk=False).ratio()  # calculates the overlap ratio
        if ratio >= threshold:
            # Check which ORF has a longer read, discard the shortest
            if len(a) >= len(b):
                try:
                    del ORF_dict[b]
                except:
                    print('Something went wrong when trying to delete item from dict')
            elif len(a) < len(b):
                try:
                    del ORF_dict[a]
                except:
                    print('Something went wrong when trying to delete item from dict')

    return ORF_dict


def N_parser(seq):
    """
    # Checks the length of the sequence and adds extra filler nucleotides 'N', so other libraries won't give warning.
    # Good Bioinformatician practises to keep seq lengths in 3 fold.
    :param seq: nucleotide sequence
    :return: seq: nucleotide with possible added trailing dummy nucleotides 'N'
    """
    if len(seq) % 3 == 1:  # Ads trailing padding to make sure the CDS is correct.
        seq = seq + 'NN'
    elif len(seq) % 3 == 2:
        seq = seq + 'N'
    else:
        seq = seq
    return seq


def altORF_finder(seq, to_stop=True, verbose=False):
    """
    Method to find alternative ORFs in the sequence.
    :param seq: The sequence that is analysed on alternative ORFs
    :param to_stop: Boolean that decides if ORF candidates without a stop codon are discarded
     or get a simulated stop codon at the end of their region (e.g 5UTR, CDS)
    :param verbose: Boolean for verbose output.
    :return: A dict of strings containing the alternate ORFs with different kind of ORFs
    """

    # predefining variables
    found_5UTR_start_codon, found_CDS_start_codon, found_3UTR_start_codon, = False, False, False
    found_5UTR_stop_codon, found_CDS_stop_codon, found_3UTR_stop_codon = False, False, False
    start_5UTR_codons, start_cds_codons, start_3UTR_codons = [], [], []
    stop_5UTR_codons, stop_CDS_codons, stop_3UTR_codons = [], [], []
    ORF_dict = {}

    seq_cds_check = seq  # We need to 2 seq variables because the check uses capitalized letters, and translation uses
    # full capitalized data.
    seq = seq.upper()

    # Need to manually translate the codons, since Biopython doesn't support custom translations tables.
    # The closest alternative translation table also seems to implement other vertebrates

    # ---------- Generates list indexes of found start and stop codons -----------
    start_list = [p for (c, p) in codons(seq, 1) if c == 'TAC'
                  or c == 'CAG'
                  or c == 'CAG']
    start_list.sort()
    start_list = np.array(start_list)

    stop_list = [p for (c, p) in codons(seq, 1) if c == 'TAA'
                 or c == 'TAG'
                 or c == 'TGA']
    stop_list.sort()
    stop_list = np.array(stop_list)

    # ---------- Checks if there's a CDS ----------
    contains_utr = re.search('[a-z]', seq_cds_check)  # Check if the sequence has utrs
    if contains_utr is None:
        has_cds = False
    else:
        seq_cds_check.replace('N', 'n')
        has_cds = [idx for idx in range(len(seq_cds_check)) if seq_cds_check[idx].isupper()]

    if verbose:
        print('~'*20)
        print('Starting ORF detection algorithm...')
        print('To stop codon has been set to\t'+str(to_stop))

    # ---------- Tree when CDS has been detected ----------
    if has_cds:
        if verbose:
            print("CDS detected")

        cds_index = [has_cds[0], has_cds[-1]]  # Define CDS location in sequence
        UTR5 = seq[:cds_index[0]]  # Define UTR locations based on known CDS location
        UTR3 = seq[cds_index[1] + 1:]

        # Check if there's start- and stop codons in UTRs and define masks
        start_5UTR_mask = start_list + 3 < cds_index[0] + 1
        start_CDS_mask = ma.getmask(ma.masked_where((start_list + 3 > cds_index[0] + 1) &
                                                    (start_list + 3 < cds_index[1] + 1), start_list))
        start_3UTR_mask = start_list + 3 > cds_index[1] + 1

        # When we find any hits in the start codons, create a codon list and switch the designated boolean
        if start_5UTR_mask.any():
            start_5UTR_codons = start_list[start_5UTR_mask]
            found_5UTR_start_codon = True

        if start_CDS_mask.any():
            start_cds_codons = start_list[start_CDS_mask]
            found_CDS_start_codon = True

        if start_3UTR_mask.any():
            start_3UTR_codons = start_list[start_3UTR_mask]
            found_3UTR_start_codon = True

        # Check if there's start- and stop codons in UTRs
        stop_5UTR_mask = stop_list + 3 < cds_index[0]
        stop_CDS_mask = ma.getmask(ma.masked_where((stop_list + 3 > cds_index[0] + 1) &
                                                   (stop_list + 3 < cds_index[1] + 1), stop_list))
        stop_3UTR_mask = stop_list > cds_index[1]

        if stop_5UTR_mask.any():
            stop_5UTR_codons = stop_list[stop_5UTR_mask]
            found_5UTR_stop_codon = True

        if stop_CDS_mask.any():
            stop_CDS_codons = stop_list[stop_CDS_mask]
            found_CDS_stop_codon = True

        if stop_3UTR_mask.any():
            stop_3UTR_codons = stop_list[stop_3UTR_mask]
            found_3UTR_stop_codon = True

        # ---------------------------------------------------------------------------
        # Check tree 5UTR start codon
        # Check for uORF & uoORF
        if verbose:
            print('Checking 5UTR for uORF & uoORF...')
        if found_5UTR_start_codon:
            if found_CDS_stop_codon:  # uoORF
                if verbose:
                    print('\tuoORF detected!')
                altORF = seq[start_5UTR_codons[0]:stop_CDS_codons[-1]]
                ORF_dict.update({altORF: ['uoORF', str(start_5UTR_codons[0])+'-'+str(stop_CDS_codons[-1])]})

            elif not found_CDS_stop_codon:
                if found_5UTR_stop_codon:  # uORF
                    if verbose:
                        print('\tuORF detected!')
                    altORF = seq[start_5UTR_codons[0]:stop_5UTR_codons[-1]]
                    ORF_dict.update({altORF: ['uORF', str(start_5UTR_codons[0])+'-'+str(stop_5UTR_codons[-1])]})

                elif not found_5UTR_stop_codon:
                    # print('Found no alternate stop codon in 5UTR, and CDS...')
                    # print('Discarded since it is most likely CDS')
                    a = 1

        elif not found_5UTR_start_codon:
            # print('\tno start codons detected in 5 UTR.')
            a = 1
        # ---------------------------------------------------------------------------
        # Check tree CDS start codon
        # Check for doORF
        if verbose:
            print('Checking CDS for possible doORF...')
        if found_CDS_start_codon:
            if found_3UTR_stop_codon:  # doORF
                if verbose:
                    print('\tdoORF detected!')
                altORF = seq[start_cds_codons[0]:stop_3UTR_codons[-1]]
                ORF_dict.update({altORF: ['doORF', str(start_cds_codons[0])+'-'+str(stop_3UTR_codons[-1])]})

            elif not found_3UTR_stop_codon:  # doORF, simulating
                if to_stop:
                    a = 1
                if not to_stop:
                    if verbose:
                        print('\tdoORF detected! (simulated)')
                    altORF = seq[start_cds_codons[0]:]
                    ORF_dict.update({altORF: ['doORF(sim)', str(start_cds_codons[0])+'-end']})
        # ---------------------------------------------------------------------------
        # Check tree 3UTR start codon
        # Check for dORF
        elif found_3UTR_start_codon:
            if verbose:
                print('\tno start codons detected in CDS.')
                print('Checking 3UTR for possible dORF...')
            if found_3UTR_stop_codon:
                if verbose:
                    print('\tdORF detected!')
                altORF = seq[start_3UTR_codons[0]:stop_3UTR_codons[-1]]
                ORF_dict.update({altORF: ['dORF', str(start_3UTR_codons[0])+'-'+str(stop_3UTR_codons[-1])]})

            elif not found_3UTR_stop_codon:
                if to_stop:
                    # print('Found no stop codon, discarting')
                    a = 1
                if not to_stop:
                    if verbose:
                        print('\tdORF detected! (simulated)')
                    altORF = seq[start_3UTR_codons[0]:]
                    ORF_dict.update({altORF: ['dORF(sim)', str(start_3UTR_codons[0])+'-end']})

    # ---------- Tree when NO CDS has been detected ----------
    elif not has_cds:
        if verbose:
            print("\tNo CDS detected")
            print('Checking sequence for lncORF...')
        if start_list.any():
            if stop_list.any():
                if verbose:
                    print('\tlncORF detected!')
                altORF = seq[start_list[0]:stop_list[-1]]
                ORF_dict.update({altORF: ['lncORF', str(start_list[0])+'-'+str(stop_list[-1])]})

            if not stop_list.any():
                if to_stop:
                    a = 1
                if not to_stop:
                    if verbose:
                        print('\tlncORF detected! (simulated)')
                    altORF = seq[start_list[0]:]
                    ORF_dict.update({altORF: ['lncORF(sim)', str(start_list[0])+'-end']})

        if not start_list.any():
            if verbose:
                print('No start codons detected.')

    return ORF_dict


def codons(seq, frame):
    """Generator function that yields DNA in one-codon blocks
    :param seq:
    :param frame:
    :return: tuple containing (codon, position relative to start)
    """
    start = frame - 1 # redreading frame is 1-based, index for the nucleotide position is 0-based
    while start + 3 <= len(seq):
        yield seq[start:start + 3], start
        start += 3


if __name__ == "__main__":
    main()
