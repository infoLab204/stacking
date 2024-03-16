import numpy as np
from scipy import stats
import datetime

from sklearn.model_selection import train_test_split
import pandas as pd

import warnings 
warnings.filterwarnings('ignore')

import sys    
import os

geno_data=pd.read_csv(sys.argv[1], sep="\t")  ## genotype data 
pheno_data=pd.read_csv(sys.argv[2], sep="\t")  ## phenotype data

pheno_idx=int(sys.argv[3]) ## select phenotype index
rand_idx=int(sys.argv[4])  ## select randon state number


X_data=geno_data.transpose()
y_target=pheno_data.iloc[:,pheno_idx]

## y normalization
y_target=(y_target-y_target.mean())/y_target.std()

t01 = datetime.datetime.now()
X_train, X_test, y_train, y_test=train_test_split(X_data, y_target, test_size=0.2, random_state=rand_idx)

## Numpy transpose
X_train_n=X_train.values
X_test_n=X_test.values
y_train_n=y_train.values
y_test_n=y_test.values

X_tr_name="X_training_test_full"
y_tr_name="y_training_test_full"

X_tr_output=open(X_tr_name, "w")
y_tr_output=open(y_tr_name, "w")

for i in X_train_n :
    for j in i :
        X_tr_output.write(f"{j}\t")
        
    X_tr_output.write("\n")

for i in X_test_n :
    for j in i :
        X_tr_output.write(f"{j}\t")
        
    X_tr_output.write("\n")
        
         
for i in y_train_n :
    y_tr_output.write(f"{i}\n")
        
for i in y_test_n :
    y_tr_output.write(f"{i}\n")

X_tr_output.close()
y_tr_output.close()

os.system("/opt/R/4.3.0/bin/Rscript base_model.R")
os.system("/opt/R/4.3.0/bin/Rscript stacking_base_model.R")
os.system("python meta_model.py")

t02 = datetime.datetime.now()
print("run time :", t02 - t01)

os.system("rm -rf *.dat")
