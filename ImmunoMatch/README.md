ImmunoMatch is a machine learning framework for deciphering the molecular rules governing the pairing of antibody chains. Fine-tuned on an antibody-specific language model ([AntiBERTa2](https://www.google.com/url?sa=t&source=web&rct=j&opi=89978449&url=https://www.biorxiv.org/content/10.1101/2023.12.12.569610v1)), ImmunoMatch learns from paired H and L sequences from single human B cells to distinguish cognate H-L pairs and randomly paired sequences. 

A total of three variants of ImmunoMatch, trained on different subsets of the data, are made available on huggingface:

| Checkpoint name | Trained on |
| --------------- | ---------- |
| [ImmunoMatch](https://huggingface.co/fraternalilab/immunomatch) | A mixture of antibodies with both κ and λ light chains |
| [ImmunoMatch-κ](https://huggingface.co/fraternalilab/immunomatch-kappa) | Antibodies with κ light chains |
| [ImmunoMatch-λ](https://huggingface.co/fraternalilab/immunomatch-lambda) | Antibodies with λ light chains |

Please note that the ImmunoMatch models are provided under a CC-BY-NC-4.0 license.

This package provides functions that can easily run immunomatch on batches of sequences stored in dataframe on your local computer, and attach the calculated pairing scores as a separate column for further analysis.
# Use

Import the package: `from ImmunoMatch import *`

`run_immunomatch_batches(data,hseq_col,lseq_col,ltype_col)`: where data is a dataframe, hseq_col is the name of the column holding the VH sequences, lseq_col is the name of the column holding the VL sequences, ltype_col is the name of the column indicating whether it's IGK or IGL.
