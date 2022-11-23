
# Cryptic Neo-antigen finder.

This python pipeline tool was developed for my major research project:
#### Automating biomarker identification for immunotherapies: Non-canonical peptides presented on MHC molecules"
 
This tool takes a single `Project-id` and `Primary site` from the National Cancer Insitute's 
GDC data portal as input, and outputs a short list of candidate peptides that are predicted 
to be good binders on MHC molecules.

## Installation
Please note that this pipeline also makes use of the R language and uses a local version of netMHCpan 4.1.
This is essential for the pipeline to run sucessfully.

#### R
R can be downloaded from: https://cran.r-project.org/.
#### netMHCpan 4.1
netMHCpan 4.1 can be downloaded from https://services.healthtech.dtu.dk/service.php?NetMHCpan-4.1 under the 'downloads' section.\
Please note the extra needed files on the website: [data.tar.gz](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/data.tar.gz) & [test.tar.gz](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/test.tar.gz).\
Make sure netMHCpan is added to $PATH [(tutorial)](https://askubuntu.com/a/60219).

#### Pipeline
This project can be installed by cloning the git repository into a new folder.
```bash
  git clone https://github.com/Sjonnie404/NeoPipeline.git neoPipe
  cd neoPipe
  conda create --name <env> --file requirements.txt
  conda activate <env>
```
    
## Usage/Examples

#### Please Note:
The duration of this script can take up multiple hours (dependent on the data and number of cores that are being used).
Running this pipeline in a seperate Bash screen or session is highly advised due to connectivity disruptions.

```bash
conda activate <env>
cd neoPipe
python3 Scripts/Director.py
```

Check the parameter list below or use the `-h` flag for a list of all parameters inside the terminal.

### Parameters
| Options | Explanation |
| ------------- | ------------- |
| -h, --help | show this help message and exit  |
| -site, --primary_site | define primary cancer site  |
| -project, --cancer_project | define cancer project  |
| -v, --verbose | Enables verbose mode to show more output |
| -t, --threads | Defines number of threads to use for parallelization, default is 1 |
| -check, --checkpoints | Enables save mode to save in-between-step files |
| -time, --add_time | Adds a timestamp to the output folder |
| -o, --output | Defines the specified output folder name |
| -rna, --rna_only | Defines the use of RNA-only mode, this discarted all non RNA-related genes |
| -stop, --to_stop | Sets the paramater for translation. to_stop = True means a stop codon also<br> needs to be found for a complete translation |
| -SBthreshold, --strong_binding_threshold | Set the strong binding threshold, this defines the<br> max rank a peptide can have to be defined as 'Strong Binder' |
| -WBthreshold, --weak_binding_threshold | Set the weak binding threshold, this defines the<br> max rank a peptide can have to be defined as 'Weak Binder' |
| -cutoffp, --peptide_cutoff_percentage | Set the cutoff for peptide selection based on percentage,<br> get overruled when absolute cutoff is used |
| -cutoffa, --peptide_cutoff_absolute | Set the cutoff for peptide selection based on absolute numbers<br> gets overrules percentage cutoff |
| -PEPinc, --peptide_inclusive | Include the peptides that have the same ER rank, but fall off due to cutoffs |

## Acknowledgements

 - [Assoc. Prof. Can Kesmir](https://tbb.bio.uu.nl/kesmir/index.html)
 - [Dr. Sabrina Santos Oliveira](https://cellbiology.science.uu.nl/research-groups/sabrina-oliveira-molecular-targeted-therapies/)
 - [Drs. J.K. van Amerongen](https://www.uu.nl/staff/JKvanAmerongen)
 - [S.D. (Shreya) Dharadhar](https://www.uu.nl/staff/SDDharadhar)
 - [Theoretical Biology & Bioinformatics group](https://tbb.bio.uu.nl/)
 - [NetMHCpan developers](https://pubmed.ncbi.nlm.nih.gov/32406916/)
## Authors

- [@Shane Pullens](https://www.github.com/Sjonnie404)


## License

MIT License

Copyright (c), 2022, Shane Ian Pullens

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

