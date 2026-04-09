from transformers import (
        RoFormerTokenizer,
        pipeline,
        RoFormerForSequenceClassification,
        Trainer,TrainingArguments
    )

import pandas as pd
import numpy as np

import torch
from datasets import load_dataset

from sklearn.metrics import roc_curve,roc_auc_score

from functools import partial

import warnings
with warnings.catch_warnings():
    warnings.simplefilter('ignore')


# Functions for the calculation of the pairing score for batches of VH-VL

def preprocess_seq(example,hseqcol="input_Hseq",lseqcol="input_Lseq"):
    return {"input_Hseq":" ".join(list(example[hseqcol])), "input_Lseq":" ".join(list(example[lseqcol]))}

def tokenize_function(examples, tokenizer, hseqcol="input_Hseq",lseqcol="input_Lseq",max_length=256,return_tensors="pt"):
    return tokenizer(examples[hseqcol], examples[lseqcol], padding="max_length", truncation=True, max_length=max_length, return_tensors=return_tensors)

def tokenize_the_datasets(df_dir,hseq_col,lseq_col,tokenizer):
    """
    Tokenize the datasets
    args:input:
    df_dir: str, the directory of the dataset
    hseq_col: str, the column name of the heavy chain sequence
    lseq_col: str, the column name of the light chain sequence
    """
    df=pd.read_csv(df_dir)
    if len(df) == 0:
        return df, None
    datasets=load_dataset("csv", data_files={"test":df_dir})
    tokenized_datasets=datasets.map(partial(preprocess_seq,hseqcol=hseq_col,lseqcol=lseq_col))
    tokenized_datasets=tokenized_datasets.map(partial(tokenize_function,tokenizer=tokenizer),batched=True)
    return df, tokenized_datasets

def pairing_scores_batches (df_dir,hseq_col,lseq_col,model_checkpoint):
  """
  Load the model and make the pairing prediction on batches of sequences
  args:input:
  df_dir: the directory of the csv files holding the sequences of pairs of VH and VL sequences
  hseq_col: the column name of the column of VH sequences
  lseq_col: the column name of the column of VL sequences
  model_checkpoint: the chesck point of the version of ImmunoMatch of your interest
  """

  tokenizer = RoFormerTokenizer.from_pretrained(model_checkpoint)
  model=RoFormerForSequenceClassification.from_pretrained(model_checkpoint)

  df,tokenized_datasets = tokenize_the_datasets (df_dir, hseq_col, lseq_col, tokenizer)

  device="cuda" if torch.cuda.is_available() else "cpu"
  batch_size=48
  args = TrainingArguments(
    f"tmp",
    eval_strategy = "no",
    save_strategy = "epoch",
    #learning_rate=2e-5,
    per_device_train_batch_size=batch_size,
    per_device_eval_batch_size=batch_size,
    #num_train_epochs=5,
    #weight_decay=0.01,
    report_to="none"
    #load_best_model_at_end=True,
    #metric_for_best_model="accuracy",
    #push_to_hub=True,
)

  if tokenized_datasets is None:
      pairing_scores = []
  else:
      trainer = Trainer(
                model,
                args,
                processing_class=tokenizer,
                )
      pred_result=trainer.predict(tokenized_datasets["test"])
      pairing_scores=torch.nn.functional.softmax(torch.tensor(pred_result.predictions),dim=1)[:,1].tolist()
      
  df["pairing_scores"]=pairing_scores

  return df

def run_immunomatch_batches(data,hseq_col,lseq_col,ltype_col):
  """
  Run the ImmunoMatch model on batches of sequences
  args:input:
  df_dir: the directory of the csv files holding the sequences of pairs of VH and VL sequences
  hseq_col: the column name of the column of VH sequences
  lseq_col: the column name of the column of VL sequences
  ltype_col: the column name of the column of light chain type
  """
  kappa_data = data.loc[data[ltype_col].apply(lambda x: "IGK" in x)]
  lambda_data = data.loc[data[ltype_col].apply(lambda x: "IGL" in x)]

  k_data_dir = "kappa_data.csv"
  l_data_dir = "lambda_data.csv"
  #kappa_data.to_csv(k_data_dir, index=False)
  #lambda_data.to_csv(l_data_dir, index=False)
  kappa_data.to_csv(k_data_dir)
  lambda_data.to_csv(l_data_dir)

  k_pairing_batch_result = pairing_scores_batches(k_data_dir,hseq_col,lseq_col, "fraternalilab/immunomatch-kappa")    
  l_pairing_batch_result = pairing_scores_batches(l_data_dir,hseq_col,lseq_col, "fraternalilab/immunomatch-lambda")

  pairing_batch_result = pd.concat([k_pairing_batch_result, l_pairing_batch_result]).sort_values(by=["Unnamed: 0"]).drop(columns=["Unnamed: 0"]).reset_index(drop=True)
  return pairing_batch_result


