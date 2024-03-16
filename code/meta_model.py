from tensorflow.keras.layers import Input, Dense
from tensorflow.keras import optimizers
from tensorflow.keras.models import Model
import numpy as np
import pandas as pd
from scipy import stats
import os
import sys


# load phenotype file
y_pheno=pd.read_csv("y_training_test_full", sep="\n", header=None)

rrBLUP_pred=pd.read_csv("rrBLUP_yHat_training_stacking_predicted.txt", sep="\n",header=None) 
BA_pred=pd.read_csv("BA_yHat_training_stacking_predicted.txt", sep="\n",header=None) 
BB_pred=pd.read_csv("BB_yHat_training_stacking_predicted.txt", sep="\n",header=None) 
BC_pred=pd.read_csv("BC_yHat_training_stacking_predicted.txt", sep="\n",header=None) 
BL_pred=pd.read_csv("BL_yHat_training_stacking_predicted.txt", sep="\n",header=None) 
BRR_pred=pd.read_csv("gBLUP_yHat_training_stacking_predicted.txt", sep="\n",header=None) 

train_n=int(y_pheno.shape[0]*0.8)

# CV stacking : test_predicted.txt
rrBLUP_test_pred=pd.read_csv("rrBLUP_yHat_test_stacking_predicted.txt", sep="\n",header=None) 
BA_test_pred=pd.read_csv("BA_yHat_test_stacking_predicted.txt", sep="\n",header=None) 
BB_test_pred=pd.read_csv("BB_yHat_test_stacking_predicted.txt", sep="\n",header=None) 
BC_test_pred=pd.read_csv("BC_yHat_test_stacking_predicted.txt", sep="\n",header=None) 
BL_test_pred=pd.read_csv("BL_yHat_test_stacking_predicted.txt", sep="\n",header=None) 
BRR_test_pred=pd.read_csv("gBLUP_yHat_test_stacking_predicted.txt", sep="\n",header=None) 

Stack_final_X_train=np.concatenate((rrBLUP_pred, BA_pred, BB_pred, BC_pred, BL_pred, BRR_pred), axis=1)
Stack_final_X_test=np.concatenate((rrBLUP_test_pred, BA_test_pred, BB_test_pred, BC_test_pred, BL_test_pred, BRR_test_pred),axis=1)

y_pheno=y_pheno.values
y_train_full=y_pheno[:train_n]
y_test=y_pheno[train_n:]

Epoch=[400]

for i in Epoch :
    inputs=Input(shape=(6,))
    hidden=Dense(12,activation='sigmoid')(inputs)
    outputs=Dense(1,activation='linear')(hidden)

    E=i
    BS=48

    linear_model=Model(inputs, outputs)

    adam=optimizers.Adam(learning_rate=0.001, beta_1=0.9, beta_2=0.999, epsilon=1e-08)
    linear_model.compile(optimizer=adam, loss='mse')

    hist=linear_model.fit(Stack_final_X_train,y_train_full,epochs=E, batch_size=BS, verbose=1, validation_split=0.2)

    test_pred=linear_model.predict(Stack_final_X_test)
    np.savetxt(f"meta_nn_predicted_test.txt",test_pred, fmt="%.8f")
    train_pred=linear_model.predict(Stack_final_X_train)
    np.savetxt(f"meta_nn_predicted_train.txt",train_pred, fmt="%.8f")

    co1,p1=stats.pearsonr(y_test[:,0], test_pred[:,0])
    co1_1,p1_1=stats.pearsonr(y_train_full[:,0], train_pred[:,0])


    meta_cor=f"meta_NN_correlation.txt"
    outfile=open(meta_cor,"w")
    outfile.write(f"meta test {co1} training {co1_1}\n")
    outfile.close()
