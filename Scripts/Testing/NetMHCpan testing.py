from mhctools import NetMHCpan4
import mhctools
import pathlib

# NOTE: This needs a netMHCpan4 installation, eather in /usr/bin/ or referenced in $PATH!

# Run NetMHCpan for alleles HLA-A*01:01 and HLA-A*02:01
predictor = NetMHCpan4(alleles=["A*02:01", "hla-a0101"])
mhctools.Net

# scan the short proteins 1L2Y and 1L3Y for epitopes
# protein_sequences = {
#   "1L2Y": "NLYIQWLKDGGPSSGRPPPS",
#   "1L3Y": "ECDTINCERYNGQVCGGPGRGLCFCGKCRCHPGFEGSACQA"
# }
protein_sequences = {
  "ENST00000612152": "MLKLYAMFLTLVFLVELVAAIVGFVFRHEIKNSFKNNYEKALKQYNSTGDYRSHAVDKIQNTLHCCGVTDYRDWTDTNYYSEKGFPKSCCKLEDCTPQRDADKVNNELIGIFLAYCLSRAITNNQYEIV",
  "ENST00000614008": "MLKLYAMFLTLVFLVELVAAIVGFVFRHEIKNSFKNNYEKALKQYNSTGDYRSHAVDKIQNTLHCCGVTDYRDWTDTNYYSEKGFPKSCCKLEDCTPQRDADKVNNEGCFIKVMTIIESEMGVVAGISFGVACFQDI"
}

binding_predictions = predictor.predict_subsequences(protein_sequences, peptide_lengths=[9, 10])

# flatten binding predictions into a Pandas DataFrame
df = binding_predictions.to_dataframe()

# epitope collection is sorted by percentile rank
# of binding predictions
for binding_prediction in binding_predictions:
    if binding_prediction.affinity < 100:
        print("Strong binder: %s" % (binding_prediction,))

print(df)
print(binding_predictions)
