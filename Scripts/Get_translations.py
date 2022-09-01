########################################################################################
#  Script translate backbone and transcript fasta non-canonical & canonical respectfully
#  Input: Backbone fasta & Transscript fasta.
#  Output: Fasta files with translated backbone en transcript sequences
#  Made by Shane Pullens, Utrecht University - Theoretical bioinformatics.
# Note, only need to connect to server when using NetMHCpan, and for deployment
########################################################################################
import os

import numpy as np
from numpy import ma
from Bio.Seq import Seq
from pathlib import Path
from datetime import datetime
import itertools as iter
from difflib import SequenceMatcher as SeqMatch
import re
import string


# TODO: Need to check if there's also backbone files with no known CDS
def main():
    timestamp = datetime.now().strftime("%d-%b-%Y-h%H-m%M")  # Note duplicate (For now) #TODO: switch out this timestamp
    project_dir = Path.cwd()
    filename_trans = '05-Aug-2022-h14-m54skin_EXCEPTIONAL_RESPONDERS-ER_20220805_125442.145250_transcripts_sequences.fasta'
    filename_back = '05-Aug-2022-h14-m54skin_EXCEPTIONAL_RESPONDERS-ER_20220805_125442.145250_backbone_sequences.fasta'
    test_dir = 'testing_breast'  # NOTE: Testdir should be removed when deploying.

    dir_path = Path(project_dir / "Output" / test_dir / "Fasta")
    # transcript_file = open(Path(project_dir / "Output" / test_dir / "Fasta" / filename_trans)).read()  # relative path
    backbone_file = open(Path(project_dir / "Output" / test_dir / "Fasta" / filename_back)).read()

    # translated_trans_file = canonical_translation(transcript_file)
    translated_backbone_file = cryptic_translation(backbone_file, verbose=False)

    print('Writing canonical translated transcription file...')
    # write_file(translated_trans_file, Path(filename_trans).stem + '_translated.fasta', dir_path, timestamp, False, True)

    print('Writing cryptic translated backbone file...')
    write_file(translated_backbone_file, Path(filename_back).stem + '_translated.fasta', dir_path, timestamp, False, True)
    print('Succesfully wrote files...')
    exit()


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


def canonical_translation(file):
    trans_file = ''
    fastas = file.replace('>', '$$>').split('$$')[1:]
    for fasta in fastas:
        header, seq = fasta.split('\n', 1)
        seq = N_parser(seq.replace('\n', ''))
        seq = Seq(seq)
        translated_seq = str(seq.translate(to_stop=True, stop_symbol=None, cds=True))
        translated_seq = cleanup_fasta(translated_seq)
        translated_fasta = header + '\n' + translated_seq

        trans_file = trans_file + translated_fasta + '\n'
    return trans_file


def cryptic_translation(file, verbose=False):
    """
    Module to automatically perfom the cryptic translation.
    V2. [current], finds ORFs and take the whole ORF sequence without gaps
    V3. also implement gaps in the ORFs


    ORF naming:
    ORF_1:  normal forward strand
    ORF_2:  forward strand +1
    ORF_3:  forward strand +2
    ORF_4:  normal reverse strand
    ORF_5:  reverse strand -1
    ORF_6:  reverse strand -2

    :param file:
    :return:
    """

    trans_file = ''
    fastas = file.replace('>', '$$>').split('$$')[1:]
    ORFs = ['ORF_1', 'ORF_2', 'ORF_3', 'ORF_4', 'ORF_5', 'ORF_6']
    alt_ORFs = {}

    for fasta in fastas:
        base_header, base_seq = fasta.split('\n', 1)
        translated_ORFs = ''

        for i, ORF in enumerate(ORFs):
            if verbose:
                print('> Started checks for:\t', ORF)
            # This checks if the strand should be reversed or not
            if i <= 2:
                seq = base_seq
            elif i >= 3:
                seq = base_seq[::-1]  # Reverses the sequence

            # This checks what reading frame should be used
            if i == 0 or i == 3:
                seq = N_parser(seq.replace('\n', ''))
                alt_ORFs = altORF_finder(seq, True, False)
            elif i == 1 or i == 4:
                seq = N_parser(seq.replace('\n', '')[1:])
                alt_ORFs = altORF_finder(seq, True, False)
            elif i == 2 or i == 5:
                seq = N_parser(seq.replace('\n', '')[2:])
                alt_ORFs = altORF_finder(seq, True, False)

            if len(alt_ORFs.keys()) <= 1:
                # print('Empty or single ORF, no need to conflate')
                a = 1
            elif len(alt_ORFs.keys()) > 1:
                # print(str(len(alt_ORFs.keys()))+' ORFs have been found.\nRunning conflation checks')
                conflated_ORFs = conflate_ORFs(alt_ORFs)  # debugging.
                alt_ORFs = conflated_ORFs  # debugging

            for alt_ORF_seq, alt_ORF_meta in alt_ORFs.items():
                translated_ORF_seq = str(Seq(alt_ORF_seq).translate(to_stop=False, stop_symbol='N', cds=False, gap='n'))
                ORF_header = base_header + '|' + ORF + '|' + alt_ORF_meta[0] + '|' + alt_ORF_meta[1]
                translated_ORF_seq = cleanup_fasta(translated_ORF_seq)
                translated_ORF = ORF_header + '\n' + translated_ORF_seq
                translated_ORFs = translated_ORFs + '\n' + translated_ORF
        trans_file = trans_file + translated_ORFs
    # print(trans_file)
    # print('here')

    return trans_file


def conflate_ORFs(ORF_dict, threshold=0.9):
    # nice to have: Choose the distance distance metric (e.g. levenstein)
    #  this could be achieved with the jellyfish library
    """
    This function checks if there's overlap in the found ORFs, if found, they are conflated.
    :param ORF_dict:
    :param threshold: The threshold on which overlapping ratio the ORFs should be inflated
    :return: returns a dictionary with the same information, with conflated ORFs
    """
    for a, b in iter.combinations(list(ORF_dict.keys()), 2):
        ratio = SeqMatch(None, a, b, autojunk=False).ratio()
        if ratio >= threshold:
            # print('Found overlapping ORFs within the specified threshold:\t',threshold)
            # print('Conflating')
            if len(a) >= len(b):
                # print('remove b')
                try:
                    del ORF_dict[b]
                except:
                    print('Something went wrong when trying to delete item from dict')
            elif len(a) < len(b):
                # print('remove a')
                try:
                    del ORF_dict[a]
                except:
                    print('Something went wrong when trying to delete item from dict')

    return ORF_dict


def N_parser(seq):
    # NOTE: print statements should be removed when done debugging
    """
    # Checks the length of the sequence and adds extra filler nucleotides 'N', so other libraries won't give warning.
    # Good Bioinformatician practises to keep seq lengths in 3 fold.
    :param seq: nucleotide sequence
    :return: seq: nucleotide with possible added trailing dummy nucleotides 'N'
    """
    if len(seq) % 3 == 1:  # Ads trailing padding to make sure the CDS is correct.
        # print('1 trailing nucleotide')
        # print(len(seq) % 3)
        # seq = rreplace(seq, '\n', 'NN\n', 1)
        seq = seq + 'NN'
        # seq = seq.replace('\n','N\n', -1)
    elif len(seq) % 3 == 2:
        # print('2 trailing nucleotides')
        # seq = rreplace(seq, '\n', 'N\n', 1)
        seq = seq + 'N'
    else:
        # print('No trailing nucleotides')
        seq = seq
    return seq


def altORF_finder(seq, lenient=False, verbose=False): # TODO: Cleanup print statements
    """
    # Note: could probably just give the 5&3 prime vars to function instead of figuring out.
    Method to find alternative ORFs in the sequence.
    # Nice to have, the paper also removes <16AA , unmaped ORFs & non-AUG ORFs
    :param seq: The sequence that is analysed on alternative ORFs
    :param lenient: Boolean that decides if ORF candidates without a start codon are discarded
     or get a simulated stop codon at the end of their region (e.g 5UTR, CDS)
    :param verbose: Boolean for verbose output.
    :return: A dict of strings containing the alternate ORFs with different kind of ORFs
    """

    # predefining variables
    found_5UTR_start_codon, found_CDS_start_codon, found_3UTR_start_codon, = False, False, False
    found_5UTR_stop_codon, found_CDS_stop_codon, found_3UTR_stop_codon = False, False, False

    start_5UTR_codons, start_cds_codons, start_3UTR_codons = [], [], []
    stop_5UTR_codons, stop_CDS_codons, stop_3UTR_codons = [], [], []
    # ORF_list = []
    ORF_dict = {}



    seq_cds_check = seq
    seq = seq.upper()

    # Need to manually translate the codons, since Biopython doesn't support custom translations tables.
    # The closest alternative translation table also seems to implement other vertebrates

    # codons_list = [seq[i:i + 3] for i in range(0, len(seq), 3)]  # Note: might not be needed.
    # print(codons_list)

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

    # Checks if there's a CDS
    contains_utr = re.search('[a-z]', seq_cds_check)  # Check if the sequence has utrs

    if contains_utr is None:
        has_cds = False
    else:
        seq_cds_check.replace('N', 'n')
        has_cds = [idx for idx in range(len(seq_cds_check)) if seq_cds_check[idx].isupper()]
    if verbose:
        print('~'*20)
        print('Starting ORF detection algorithm...')
        print('Lenient mode is set to\t'+str(lenient))
    if has_cds:
        if verbose:
            print("CDS detected")
        # Define cds location in sequence
        cds_index = [has_cds[0], has_cds[-1]]
        UTR5 = seq[:cds_index[0]] + '-'
        UTR3 = '-' + seq[cds_index[1] + 1:]

        # Check if there's start- and stop codons in UTRs
        start_5UTR_mask = start_list + 3 < cds_index[0] + 1
        start_CDS_mask = ma.getmask(ma.masked_where((start_list + 3 > cds_index[0] + 1) &
                                                    (start_list + 3 < cds_index[1] + 1), start_list))
        start_3UTR_mask = start_list + 3 > cds_index[1] + 1

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
                    print()
                    # print('Found no alternate stop codon in 5UTR, and CDS...')
                    # print('Discarded since it is most likely CDS')

        elif not found_5UTR_start_codon:
            # print('\tno start codons detected in 5 UTR.')
            print()
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
                if lenient:
                    if verbose:
                        print('\tdoORF detected! (simulated)')
                    altORF = seq[start_cds_codons[0]:]
                    ORF_dict.update({altORF: ['doORF(sim)', str(start_cds_codons[0])+'-end']})
                if not lenient:
                    print()
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
                if lenient:
                    if verbose:
                        print('\tdORF detected! (simulated)')
                    altORF = seq[start_3UTR_codons[0]:]
                    ORF_dict.update({altORF: ['dORF(sim)', str(start_3UTR_codons[0])+'-end']})

                if not lenient:
                    # print('Found no stop codon, discarting')
                    print()

    # ---------------------------------------------------------------------------
    # Check tree no CDS
    # Check for lncORF
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
                if lenient:
                    if verbose:
                        print('\tlncORF detected! (simulated)')
                    altORF = seq[start_list[0]:]
                    ORF_dict.update({altORF: ['lncORF(sim)', str(start_list[0])+'-end']})
                if not lenient:
                    print()

        if not start_list.any():
            if verbose:
                print('No start codons detected.')

    return ORF_dict


def codons(seq, frame):
    """Generator function that yields DNA in one-codon blocks

    returns a tuple containing (codon, position relative to start)
    note: reading frame is 1-based, index for the nucleotide position is 0-based
    """
    start = frame - 1
    while start + 3 <= len(seq):
        yield seq[start:start + 3], start
        start += 3


def find_all(sub, a_str):  # NOTE: might not be needed anymore
    # [line[i:i+n] for i in range(0, len(line), n)]
    """
    Generator that finds all the given strings in another string.
    :param sub: target string to find
    :param a_str: string to find substring in.
    :return: generator object, containing the indexes of the found substrings.
    """
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += len(sub)  # use start += 1 to find overlapping matches


def rreplace(s, old, new, occurrence):
    """
    Replaces the last occurrence in string
    :param s: String
    :param old: substring replace
    :param new: replace to put in check
    :param occurrence: Nth occurrences from the right will be replaced
    :return: string with replaced from the right.
    """
    li = s.rsplit(old, occurrence)
    return new.join(li)



# NOTE: This is a duplicate from 'Fetch_transcript_seq.py and should be deleted when a package is created.
def write_file(text, filename, path, timestamp, add_timestamp=True, overwrite_check=False, verbose=False):
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
        filename = timestamp + filename  # +'.fasta' This was only here for testing right?

    # Check if the user already added an extension, skip adding another.
    filename = Path(filename)
    if filename.suffix not in ('.fasta', '.fa'):
        filename = str(filename) + '.fasta'

    absolute_path = Path(path / filename)

    # Check if the user enabled overwrite check. This checks if the user is going to overwrite a file and stops it.
    if overwrite_check:
        if os.path.isfile(absolute_path):
            print('>!! File: ' + str(filename) + ' already exists, aborted writing.')
            return None

    try:
        with open(absolute_path, 'w') as f:
            f.write(text)
        if verbose:
            print('>>> succesfully written ' + filename)
    except:
        print('! An error occurred while trying to write ' + filename)
    return None


# def sequence_analysis(fasta): DEPRECATED FUNCTION
#     """
#     This should be used for the cryptic translation
#     :param fasta:
#     :return:
#     """
#     # DEPRECATED
#     rm_lower = str.maketrans('', '', string.ascii_lowercase)
#     rm_upper = str.maketrans('', '', string.ascii_uppercase)
#
#     # TODO: Check what is the best way to characterize the UTRS, (now going for a added '|' before and after UTR)
#     header, sequence = fasta.split('\n', 1)
#     cds = sequence.translate(rm_lower).replace(' ', '').replace('\n', '')
#
#     # Nice to have: find better solution for splitting and saving the split string.
#     UTR_5, temp_seq = cds.replace('ATG', '1ATG', 1).split('1', 1)
#
#     print(cds)
#     print(len(cds))
#     print('-' * 80)
#     # print(len(cds) % 3)
#     testing = Seq(cds)
#
#     # n = 3
#     # split_strings = [cds[index: index + n] for index in range(0, len(cds), n)]
#     # print(split_strings)
#     # print(split_strings.index('ATG'))
#     # print('-'*80)
#
#     # print('Simple translation:')
#     print(testing.translate())
#     print(len(testing.translate()))
#     print(testing.translate().find('M'))
#     exit()
#     print('cds translation:')
#     print(len(testing.translate(cds=True)))
#     print(testing.translate(cds=True))
#     exit()
#
#     # print(len(UTR_5))
#     # print(len(UTR_5) % 3)
#     print(UTR_5)
#     print('-----')
#     print(cleanup_fasta(UTR_5))
#
#     temp_seq = temp_seq.translate(rm_lower).replace(' ', '').replace('\n', '')
#
#     refurbished_seq = cleanup_fasta(temp_seq)
#     print(len(refurbished_seq))
#     print(len(refurbished_seq) % 3)
#
#     # seq = N_parser(seq)
#     # seq = seq.replace('\n', '',)
#
#     # print('!!!!!!!!')
#     # print(seq)
#     # print(len(seq) % 3)
#     # gene = Seq(seq)
#
#     # print(len(gene))
#     # print(gene.translate())
#     # print(test)
#     # try:
#     #     print(gene.translate(cds=True))
#     #     print('Found gene in:\t' + str(trans_id))
#     #     print('adding to file...')
#     # except:
#     #     print('No gene found in:\t' + str(trans_id))
#
#     # cds = seq[prime_extender_5:-prime_extender_3-1]
#     # seq = prime5+cds+prime3
#     # fasta = header+'\n'+seq
#
#     return None

# Note: This might need to be changed when adding the scripts to a pipeline
if __name__ == "__main__":
    main()
