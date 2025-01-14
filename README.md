# ImmunoMatch
This is the repo for the antibody heavy and light chain pairing tool - ImmunoMatch

ImmunoMatch is a protein language model finetuned from AntiBERTa2, aiming at investigating the heavy and light chain pairing preferences in antibody. The input sequence to the model should be a pair of sequences of VH and VL domains.

Different variants of ImmunoMatch is available on huggingface, according to the use of interest:

| Checkpoint name | Trained on |
| --------------- | ---------- |
| [ImmunoMatch](https://huggingface.co/fraternalilab/immunomatch) | A mixture of antibodies with both κ and λ light chains |
| [ImmunoMatch-κ](https://huggingface.co/fraternalilab/immunomatch-kappa) | Antibodies with κ light chains |
| [ImmunoMatch-λ](https://huggingface.co/fraternalilab/immunomatch_lambda) | Antibodies with λ light chains |
