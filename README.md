
# Cryptic Neo-antigen finder.

This python pipeline tool was developed for my major research project:
#### Automating biomarker identification for immunotherapies: Non-canonical peptides presented on MHC molecules"
 
This tool takes a single `Project-id` and `Primary site` from the National Cancer Insitute's 
GDC data portal as input, and outputs a short list of candidate peptides that are predicted 
to be good binders on MHC molecules.

## Installation

This project can be installed by cloning the git repository into a new folder.

```bash
  git clone https://github.com/Sjonnie404/NeoPipeline.git neoPipe
  cd neoPipe
  conda create --name <env> --file requirements.txt
  conda activate <env>
```
    
## Usage/Examples

```bash
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
| -stop, --to_stop | Sets the paramater for translation. to_stop = True means a stop codon also\ needs to be found for a complete translation |
| -SBthreshold, --strong_binding_threshold | Set the strong binding threshold, this defines the\ max rank a peptide can have to be defined as 'Strong Binder' |
| -WBthreshold, --weak_binding_threshold | Set the weak binding threshold, this defines the\ max rank a peptide can have to be defined as 'Weak Binder' |
| -cutoffp, --peptide_cutoff_percentage | Set the cutoff for peptide selection based on percentage,\ get overruled when absolute cutoff is used |
| -cutoffa, --peptide_cutoff_absolute | Set the cutoff for peptide selection based on absolute numbers\ gets overrules percentage cutoff |
| -PEPinc, --peptide_inclusive | Include the peptides that have the same ER rank, but fall off due to cutoffs |

## Acknowledgements

  TODO: Add URLS
 - [UU](https://)
 - [TBB](https://)
 - [Computing source](https://)
 - [NetMHCpan developers](https://)
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

