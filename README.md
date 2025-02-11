# ImmunoMatch

![ImmunoMatch logo](ImmunoMatch_logo.png)

ImmunoMatch is a machine learning framework for deciphering the molecular rules governing the pairing of antibody chains. Fine-tuned on an antibody-specific language model ([AntiBERTA2](https://www.google.com/url?sa=t&source=web&rct=j&opi=89978449&url=https://www.biorxiv.org/content/10.1101/2023.12.12.569610v1)), ImmunoMatch learns from paired H and L sequences from single human B cells to distinguish cognate H-L pairs and randomly paired sequences. 

A total of three variants of ImmunoMatch, trained on different subsets of the data, are made available on huggingface:

| Checkpoint name | Trained on |
| --------------- | ---------- |
| [ImmunoMatch](https://huggingface.co/fraternalilab/immunomatch) | A mixture of antibodies with both κ and λ light chains |
| [ImmunoMatch-κ](https://huggingface.co/fraternalilab/immunomatch-kappa) | Antibodies with κ light chains |
| [ImmunoMatch-λ](https://huggingface.co/fraternalilab/immunomatch-lambda) | Antibodies with λ light chains |


### Try it out on Google Colab
`Run_ImmunoMatch.ipynb` contains example code on how to apply any ImmunoMatch model to obtain H-L pairing scores for a given VH-VL sequence pair, or to annotate sequences in batch upon supplying a CSV. You can also try it out on Google Collaboratory:

[![Google Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Fraternalilab/ImmunoMatch/blob/main/Run_ImmunoMatch.ipynb)

### Figure reproducibility
Folder `figure_code` contains all Python and R code used to generate figure panels in the manuscript.

## Cite

If you have used any of the ImmunoMatch models in your research please cite:

```

```


