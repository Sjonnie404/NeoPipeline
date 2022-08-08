########################################################################################
#  Script translate backbone and transcript fasta non-canonical & canonical respectfully
#  Input: Backbone fasta & Transscript fasta.
#  Output: Fasta files with translated backbone en transcript sequences
#  Made by Shane Pullens, Utrecht University - Theoretical bioinformatics.
# Note, only need to connect to server when using NetMHCpan, and for deployment
########################################################################################
import os
from Bio.Seq import Seq
from pathlib import Path
from datetime import datetime
import re
import string

# TODO: Need to check if there's also backbone files with no known CDS
def main():
    timestamp = datetime.now().strftime("%d-%b-%Y-h%H-m%M")  # Note duplicate (For now) #TODO: switch out this timestamp
    project_dir = Path.cwd()
    filename_trans = '05-Aug-2022-h14-m54skin_EXCEPTIONAL_RESPONDERS-ER_20220805_125442.145250_transcripts_sequences.fasta'
    filename_back = '05-Aug-2022-h14-m54skin_EXCEPTIONAL_RESPONDERS-ER_20220805_125442.145250_backbone_sequences.fasta'
    test_dir = 'testing_breast' # NOTE: Testdir should be removed when deploying.

    dir_path = Path(project_dir / "Output" / test_dir / "Fasta")
    transcript_file = open(Path(project_dir / "Output" / test_dir / "Fasta" / filename_trans)).read()  # relative path
    backbone_file = open(Path(project_dir / "Output" / test_dir / "Fasta" / filename_back)).read()

    translated_trans_file = canonical_translation(transcript_file)
    translated_backbone_file = cryptic_translation(backbone_file)

    print('Writing translated transcription file...')
    write_file(translated_trans_file, Path(filename_trans).stem+'_translated.fasta', dir_path, timestamp, False, True)

    print('Writing cryptic translated backbone file...')
    write_file(translated_trans_file, Path(filename_trans).stem+'_translated.fasta', dir_path, timestamp, False, True) # Note needs to be converted to backbone file.
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
        translated_fasta = header+'\n'+translated_seq

        trans_file = trans_file+translated_fasta+'\n'
    return trans_file


def cryptic_translation(file):
    """
    Module to automatically perfom the cryptic translation.

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

    for fasta in fastas:
        base_header, base_seq = fasta.split('\n', 1)
        translated_ORFs = ''

        for i, ORF in enumerate(ORFs):
            print('Started checks for:\t', ORF)
            # This checks if the strand should be reversed or not
            if i <= 2:
                seq = base_seq
            elif i >= 3:
                seq = base_seq[::-1]   # Reverses the sequence

            # This checks what reading frame should be used
            if i == 0 or i == 3:
                seq = N_parser(seq.replace('\n', ''))
                print(base_header)
                altORF_finder(seq)
            elif i == 1 or i == 4:
                seq = N_parser(seq.replace('\n', '')[1:])
            elif i == 2 or i == 5:
                seq = N_parser(seq.replace('\n', '')[2:])

            translated_seq = str(Seq(seq).translate(to_stop=False, stop_symbol='N', cds=False, gap='n'))
            header = base_header+'_'+ORF

            translated_seq = cleanup_fasta(translated_seq)
            translated_ORF = header+'\n'+translated_seq
            translated_ORFs = translated_ORFs + '\n' + translated_ORF

        trans_file = trans_file+translated_ORFs

    exit()
    return trans_file


def old_cryptic_translation(file):
    """
    Module to automatically perfom the cryptic translation.

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
    # TODO: Find a way to simulate translation for non-canonical proteins.
    # NOTE: v1: translate whole sequence (with 6 ORFs) without start or stop.


    # TODO: This should be done after the ORF finder.
    trans_file = ''
    fastas = file.replace('>', '$$>').split('$$')[1:]
    ORFs = ['ORF_1', 'ORF_2', 'ORF_3', 'ORF_4', 'ORF_5', 'ORF_6']

    for fasta in fastas:
        base_header, base_seq = fasta.split('\n', 1)
        print(base_header)
        print('\n')
        print(base_seq)
        exit()
        seq = base_seq
        translated_ORFs = ''

        for i, ORF in enumerate(ORFs):
            # This checks if the strand should be reversed or not
            if i <= 2:
                seq = base_seq
            elif i >= 3:
                seq = base_seq[::-1]   # Reverses the sequence

            # This checks what reading frame should be used
            if i == 0 or i == 3:
                seq = N_parser(seq.replace('\n', ''))
            elif i == 1 or i == 4:
                seq = N_parser(seq.replace('\n', '')[1:])
            elif i == 2 or i == 5:
                seq = N_parser(seq.replace('\n', '')[2:])

            translated_seq = str(Seq(seq).translate(to_stop=False, stop_symbol='N', cds=False, gap='n'))
            header = base_header+'_'+ORF

            translated_seq = cleanup_fasta(translated_seq)
            translated_ORF = header+'\n'+translated_seq
            translated_ORFs = translated_ORFs + '\n' + translated_ORF

        trans_file = trans_file+translated_ORFs

    exit()
    return trans_file

# Note: deprecated by using pathlib module.
# def read_data(path):
#     wd = os.getcwd()
#     wd = wd.rsplit('\\', 1)[0]
#     file = open(wd+path, 'r')
#     file = file.read()
#
#     return file


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
        #seq = rreplace(seq, '\n', 'NN\n', 1)
        seq = seq+'NN'
        # seq = seq.replace('\n','N\n', -1)
    elif len(seq) % 3 == 2:
        # print('2 trailing nucleotides')
        #seq = rreplace(seq, '\n', 'N\n', 1)
        seq = seq+'N'
    else:
        # print('No trailing nucleotides')
        seq = seq
    return seq


def altORF_finder(seq):
    """

    :param seq:
    :return:
    """
    seq = seq.replace('N', 'n')
    res = [idx for idx in range(len(seq)) if seq[idx].isupper()]
    if res:
        print("CDS detected")
        cds_index = [res[0], res[-1]]
        UTR5 = seq[:cds_index[0]]
        UTR3 = seq[cds_index[1]+1:]

        codon_finder

        print(cds_index)

        print(seq)
        print(UTR5)
        print(seq[cds_index[0]:cds_index[1]+1])
        print(UTR3)

    if not res:
        print("There's not a CDS")

    exit()


    # [m.start() for m in re.finditer('(?=tt)', 'ttt')]
    hit_list = [m.start() for m in re.finditer('(?=ATG)', seq)]
    print(hit_list)
    exit()
    print('First found start codon in the CDS')
    print('Last found stop codon in the CDS')
    ORF_list = []



    exit()
    return ORF_list

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



def sequence_analysis(fasta):
    """
    This should be used for the cryptic translation Note: WIP
    :param fasta:
    :return:
    """
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

### Interesting code ######
# >>> from Bio import SeqIO
# >>> record = SeqIO.read("NC_005816.fna", "fasta")
# >>> table = 11
# >>> min_pro_len = 100
# Here is a neat trick using the Seq object’s split method to get a list of all the possible ORF translations in the six reading frames:
#
# >>> for strand, nuc in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
# ...     for frame in range(3):
# ...         length = 3 * ((len(record)-frame) // 3) #Multiple of three
# ...         for pro in nuc[frame:frame+length].translate(table).split("*"):
# ...             if len(pro) >= min_pro_len:
# ...                 print("%s...%s - length %i, strand %i, frame %i" \
# ...                       % (pro[:30], pro[-3:], len(pro), strand, frame))
# GCLMKKSSIVATIITILSGSANAASSQLIP...YRF - length 315, strand 1, frame 0
# KSGELRQTPPASSTLHLRLILQRSGVMMEL...NPE - length 285, strand 1, frame 1
# GLNCSFFSICNWKFIDYINRLFQIIYLCKN...YYH - length 176, strand 1, frame 1
# VKKILYIKALFLCTVIKLRRFIFSVNNMKF...DLP - length 165, strand 1, frame 1
# NQIQGVICSPDSGEFMVTFETVMEIKILHK...GVA - length 355, strand 1, frame 2
# RRKEHVSKKRRPQKRPRRRRFFHRLRPPDE...PTR - length 128, strand 1, frame 2
# TGKQNSCQMSAIWQLRQNTATKTRQNRARI...AIK - length 100, strand 1, frame 2
# QGSGYAFPHASILSGIAMSHFYFLVLHAVK...CSD - length 114, strand -1, frame 0
# IYSTSEHTGEQVMRTLDEVIASRSPESQTR...FHV - length 111, strand -1, frame 0
# WGKLQVIGLSMWMVLFSQRFDDWLNEQEDA...ESK - length 125, strand -1, frame 1
# RGIFMSDTMVVNGSGGVPAFLFSGSTLSSY...LLK - length 361, strand -1, frame 1
# WDVKTVTGVLHHPFHLTFSLCPEGATQSGR...VKR - length 111, strand -1, frame 1
# LSHTVTDFTDQMAQVGLCQCVNVFLDEVTG...KAA - length 107, strand -1, frame 2
# RALTGLSAPGIRSQTSCDRLRELRYVPVSL...PLQ - length 119, strand -1, frame 2
######################################################################################################################


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
        filename = timestamp + filename  #+'.fasta' This was only here for testing right?

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
            print('>>> succesfully written '+filename)
    except:
        print('! An error occurred while trying to write '+filename)
    return None


# Note: This might need to be changed when adding the scripts to a pipeline
if __name__ == "__main__":
    main()