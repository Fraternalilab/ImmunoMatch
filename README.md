# ImmunoMatch

![ImmunoMatch logo](ImmunoMatch_logo.png)

ImmunoMatch is a machine learning framework for deciphering the molecular rules governing the pairing of antibody chains. Fine-tuned on an antibody-specific language model ([AntiBERTA2](https://www.google.com/url?sa=t&source=web&rct=j&opi=89978449&url=https://www.biorxiv.org/content/10.1101/2023.12.12.569610v1)), ImmunoMatch learns from paired H and L sequences from single human B cells to distinguish cognate H-L pairs and randomly paired sequences. 

A total of three variants of ImmunoMatch, trained on different subsets of the data, are made available on huggingface:

| Checkpoint name | Trained on |
| --------------- | ---------- |
| [ImmunoMatch](https://huggingface.co/fraternalilab/immunomatch) | A mixture of antibodies with both κ and λ light chains |
| [ImmunoMatch-κ](https://huggingface.co/fraternalilab/immunomatch-kappa) | Antibodies with κ light chains |
| [ImmunoMatch-λ](https://huggingface.co/fraternalilab/immunomatch-lambda) | Antibodies with λ light chains |

Please note that the ImmunoMatch models are provided under a CC-BY-NC-4.0 license.

### Try it out on Google Colab
`Run_ImmunoMatch.ipynb` contains example code on how to apply any ImmunoMatch model to obtain H-L pairing scores for a given VH-VL sequence pair, or to annotate sequences in batch upon supplying a CSV. You can also try it out on Google Collaboratory:

[![Google Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Fraternalilab/ImmunoMatch/blob/main/Run_ImmunoMatch.ipynb)

### Requirements
There are no specific prerequisites to use ImmunoMatch beyond standard installation of Huggingface libraries on Python. On a clean virtual environment on Google Colab, the installation of these libraries took around 1 minute.

### Figure reproducibility
Folder `figure_code` contains all Python and R code used to generate figure panels in the manuscript.

## Cite

If you have used any of the ImmunoMatch models in your research please cite:

```
@article {Guo2025.02.11.637677,
	author = {Guo, Dongjun and Dunn-Walters, Deborah K and Fraternali, Franca and Ng, Joseph CF},
	title = {ImmunoMatch learns and predicts cognate pairing of heavy and light immunoglobulin chains},
	elocation-id = {2025.02.11.637677},
	year = {2025},
	doi = {10.1101/2025.02.11.637677},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2025/02/15/2025.02.11.637677},
	eprint = {https://www.biorxiv.org/content/early/2025/02/15/2025.02.11.637677.full.pdf},
	journal = {bioRxiv}
}
```


