import os
from Bio.Seq import Seq


def main():
    transcript_file = read_data("\Output\\Fasta\\04-Mar-2022-h14-m33_Testing_api_transcripts.fasta")  # relative path
    backbone_file = read_data("\Output\\Fasta\\04-Mar-2022-h14-m33_Testing_api_backbone.fasta")

    translated_trans_file = canonical_translation(transcript_file)
    print(translated_trans_file)
    exit()

    translated_backbone_file = cryptic_translation(backbone_file)

    # transcript_fastas = get_fasta_list(transcript_file)[1:]
    # backbone_fastas = get_fasta_list(backbone_file)[1:]
    #
    # for transcript_fasta in transcript_fastas:
    #     transcript_fasta_trans = get_translation(transcript_fasta)

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
        seq = N_parser(seq)
        seq = Seq(seq.replace('\n', ''))
        translated_seq = str(seq.translate(to_stop=True, stop_symbol=None))
        translated_seq = cleanup_fasta(translated_seq)
        translated_fasta = header+'\n'+translated_seq

        trans_file = trans_file+translated_fasta+'\n'
    return trans_file


def cryptic_translation(file):
    trans_file = ''
    fastas = file.replace('>', '$$>').split('$$')[1:]
    for fasta in fastas:

        header, seq = fasta.split('\n', 1)
        print(seq)
        exit()
        # TODO: Find a way to simulate translation for non-canonical proteins.

        seq = N_parser(seq)
        seq = Seq(seq.replace('\n', ''))

        translated_seq = str(seq.translate(to_stop=True, stop_symbol=None)) # TODO: Find out why Bio still gives warning
        translated_seq = cleanup_fasta(translated_seq)
        translated_fasta = header+'\n'+translated_seq
        # print('-'*80)
        # print(translated_fasta)
        trans_file = trans_file+translated_fasta+'\n'
    print(trans_file)
    exit()
    return trans_file


def read_data(path):
    wd = os.getcwd()
    wd = wd.rsplit('\\', 1)[0]
    file = open(wd+path, 'r')
    file = file.read()

    return file


def N_parser(seq):
    """
    # Checks the length of the sequence and adds extra filler nucleotides 'N', so other libraries won't give warning.
    # Good Bioinformatician practises to keep seq lengths in 3 fold.
    :param seq: nucleotide sequence
    :return: seq: nucleotide with possible added trailing dummy nucleotides 'N'
    """
    if len(seq) % 3 == 1:  # Ads trailing padding to make sure the CDS is correct.
        # print('1 trailing nucleotide')
        # print(len(seq) % 3)
        seq = rreplace(seq, '\n', 'NN\n', 1)
        # seq = seq.replace('\n','N\n', -1)
    elif len(seq) % 3 == 2:
        # print('2 trailing nucleotides')
        # print(len(seq) % 3)
        seq = rreplace(seq, '\n', 'N\n', 1)
    else:
        # print('No trailing nucleotides')
        seq = seq
    return seq


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

# Note: This might need to be changed when adding the scripts to a pipeline
if __name__ == "__main__":
    main()