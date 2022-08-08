from mhctools import NetMHCpan4_EL
import mhctools
from pathlib import Path
from mhctools import NetMHCpan4
# NOTE: This needs a netMHCpan4 installation, eather in /usr/bin/ or referenced in $PATH!


# TODO: Testing if this works the same as the webversion (gold standard)

# Note: test fasta file to perform checks.
# >ENST00000598322.3|ENSG00000269825|ht_seq:_x_|gene_name:_y_|gene_type:_z_
# MLHEKATKKTKEKETRMALPQGCLTFRDVAIEFSLEEWKCLNPAQRALYRAVMLENYRNL
# ESVDSSLKSMTEFSSTGHGNTGEVIHTGTLQRHKSHHIGDFCFPEIKKDIHDFEFQWQEI
# KRNGHEAPMTKIKKLTGSTDRRDQRHAGNKPIKDQLGSSFYSHLPELHEFQTEGKIDNQV
# EKSINNASLVSTSQRISCRLKTHISNRYGKNFLHSSLLTQIQEEHMREKPFQCNECGKAF
# NYSSHLRRHHVTHSGEKQYKCDVCGKVFHQKQYLAWHHRVHTGEKPYKCNECSKTFGHKS
# SLTRHHRLHTGEKPYKCNECGKTFSQTSSLVGHRRRHTGEKPYKCEGCDKVYSCRSQLET
# HRRIHTGEKPYKCKVCDKAFRHNSCLSRHNRVHTEEKPYTCNECGKVFQRDSYLAQHQRV
# HTGEKPYTCNECGKVFNQKAHLACHYRLHTGEKPYKCNECGKTFSQKSSLVGHRRLHTGE
# KPYNCHECGKTFARNSSLLIHKAIHTGEKPYKCNECGKVFNQQSNLAQHQRVHTGEKPYR
# CNECGKTFSHMSSFVYHYRLHSGEKPYKCNECGKTFSHMSSFVCHHRLHTGENPYKCNEC
# GKAFSGQSSLIHHQAIHGIGKLYKCNDSHKVLSNATSIANH



def main():
    # Run NetMHCpan for alleles HLA-A*01:01 and HLA-A*02:01
    predictor = NetMHCpan4_EL(alleles=["A*02:01", "hla-a0101"])
    project_dir = Path.cwd()
    filename = '28-Jun-2022-h14-m14overlapping_genes_transcripts_TEST_translated.fasta'
    test_dir = 'testing_breast' # Note this should be removed for deployment.
    # transcript_file = open(Path(project_dir / "Output" / "Fasta" / filename_trans))  # relative path
    # fasta = open(Path(project_dir / 'Output' / test_dir / 'Fasta' / filename))
    # fasta = fasta.read()

    # fasta_dict = fastaToDict(fasta)
    # protein_sequences = fasta_dict

    # scan the short proteins 1L2Y and 1L3Y for epitopes
    # protein_sequences = {
    #   "1L2Y": "NLYIQWLKDGGPSSGRPPPS",
    #   "1L3Y": "ECDTINCERYNGQVCGGPGRGLCFCGKCRCHPGFEGSACQA"
    # }
    protein_sequences = {
        'ENST00000598322.3': 'MLHEKATKKTKEKETRMALPQGCLTFRDVAIEFSLEEWKCLNPAQRALYRAVMLENYRNLESVDSSLKSMTEFSSTGHGNTGEVIHTGTLQRHKSHHIGDFCFPEIKKDIHDFEFQWQEIKRNGHEAPMTKIKKLTGSTDRRDQRHAGNKPIKDQLGSSFYSHLPELHEFQTEGKIDNQVEKSINNASLVSTSQRISCRLKTHISNRYGKNFLHSSLLTQIQEEHMREKPFQCNECGKAFNYSSHLRRHHVTHSGEKQYKCDVCGKVFHQKQYLAWHHRVHTGEKPYKCNECSKTFGHKSSLTRHHRLHTGEKPYKCNECGKTFSQTSSLVGHRRRHTGEKPYKCEGCDKVYSCRSQLETHRRIHTGEKPYKCKVCDKAFRHNSCLSRHNRVHTEEKPYTCNECGKVFQRDSYLAQHQRVHTGEKPYTCNECGKVFNQKAHLACHYRLHTGEKPYKCNECGKTFSQKSSLVGHRRLHTGEKPYNCHECGKTFARNSSLLIHKAIHTGEKPYKCNECGKVFNQQSNLAQHQRVHTGEKPYRCNECGKTFSHMSSFVYHYRLHSGEKPYKCNECGKTFSHMSSFVCHHRLHTGENPYKCNECGKAFSGQSSLIHHQAIHGIGKLYKCNDSHKVLSNATSIANH'
    }

    # protein_sequences = {
    #   "ENST00000612152": "MLKLYAMFLTLVFLVELVAAIVGFVFRHEIKNSFKNNYEKALKQYNSTGDYRSHAVDKIQNTLHCCGVTDYRDWTDTNYYSEKGFPKSCCKLEDCTPQRDADKVNNELIGIFLAYCLSRAITNNQYEIV",
    #   "ENST00000614008": "MLKLYAMFLTLVFLVELVAAIVGFVFRHEIKNSFKNNYEKALKQYNSTGDYRSHAVDKIQNTLHCCGVTDYRDWTDTNYYSEKGFPKSCCKLEDCTPQRDADKVNNEGCFIKVMTIIESEMGVVAGISFGVACFQDI"
    # }

    print('Predicting binders...')
    binding_predictions = predictor.predict_subsequences(protein_sequences, peptide_lengths=[10]) # because there's so much overlap, only take 9 peptides.
    # print(binding_predictions)

    # flatten binding predictions into a Pandas DataFrame
    output_df = binding_predictions.to_dataframe()


   # epitope collection is sorted by percentile rank
   # of binding predictions
   #  for binding_prediction in binding_predictions:
   #      if binding_prediction.affinity < 100:
   #          print("Strong binder: %s" % (binding_prediction,))
    exit()
    # print(df)
    print('Writing df to file...')
    filename = str(filename)+'_MHCpan.csv'
    output_df.to_csv(Path(project_dir / 'Output' / 'MHCpan' / filename))
    print('Succesfully wrote df to file!')
    # print(binding_predictions)
    return None


def fastaToDict(fasta_string):
    """
    NOTE: a lot of data is duplicate of get_translation-canonical_translations function
    :param fasta:
    :return:
    """
    fastas = fasta_string.replace('>', '$$>').split('$$')[1:]
    fasta_dict = {}
    for fasta in fastas:
        header, seq = fasta.split('\n', 1)
        id = header.split('|',1)[0].split('.')[0].replace('>','')
        seq = seq.replace('\n', '').replace('X','')  # X is used as stopcodon, and not supported by mhcPAN
        fasta_dict[id] = seq
    return fasta_dict


# Note: This might need to be changed when adding the scripts to a pipeline
if __name__ == "__main__":
    main()