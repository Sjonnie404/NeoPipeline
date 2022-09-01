# WARNING: Reprecated
# Note: Please do not use, this package is only available in pip, and not in conda, this means dat this package needs to
#  be manually installed in a pip directory, inside the virtual environment.
#  Link to answer: https://stackoverflow.com/questions/41060382/using-pip-to-install-packages-to-anaconda-environment


from PyNetMHCpan import NetMHCpan


h = NetMHCpan.Helper(mhc_class='I',
                     alleles=['HLA-A03:02', 'HLA-B07:02'],
                     n_threads=12)

# If the configuration needs to be updates, a prompt should appear (i.e. you need to set the paths
# to NetMHCpan/NetMHCIIpan).

# If it is the first time you have run the program and a prompt did not appear, you should
# run the following line (as it is incredibly unlikely that the default paths to NetMHCpan and
# NetMHCIIpan are correct for your system)

h.update_config()

# A prompt will appear showing your the current configuration and asking if you would like to update.

# Note that if you only use one of the tools (i.e. NetMHCpan or NetMHCIIpan) you can leave the path to
# the other as it is or change it to something nonsensical. The program will only check the path of the tool
# it needs (according to the mhc_class parameter).

h.add_peptides(['N.ENSDESYSEEEEEMPDSD.E',
                'F.EC[+57]NNSESGYLYC[+57]C[+57]GTC[+57]Y.Y',
                'K.AM[+15]IM[+15]C[+57]QGC[+57]GAFC[+57]HDDC[+57]IG.P'])

# Note that there is no need to remove modifications or previous/next amino acids from your
# search engine output! That is automatically taken care of.

# now run the predictor.

h.make_predictions()

# Once it is done, access them as follows:

h.predictions

# returns a pandas dataframe:

#               Peptide      Allele     Rank      Binder
# 0   ECNNSESGYLYCCGTCY  HLA-A03:02   85.000  Non-binder
# 1  ENSDESYSEEEEEMPDSD  HLA-A03:02  100.000  Non-binder
# 2  AMIMCQGCGAFCHDDCIG  HLA-A03:02  100.000  Non-binder
# 3   ECNNSESGYLYCCGTCY  HLA-B07:02   92.500  Non-binder
# 4  ENSDESYSEEEEEMPDSD  HLA-B07:02  100.000  Non-binder
# 5  AMIMCQGCGAFCHDDCIG  HLA-B07:02  100.000  Non-binder

# Note that the dataframe is in long format. If you have multiple alleles, there will be multiple
# obeservations for each peptide (one per allele).