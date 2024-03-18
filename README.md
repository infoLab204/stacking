# stacking : Stacked generalization as a computational method for the genomic selection with Python and R
Sunhee Kim, Sang-Ho Chu, Yong-Jin Park and Chang-Yong Lee

Stacking involves the integration of multiple models, called base models, along with an additional model, called the meta-model. The base models generate predictions for the data, and the meta-model is tasked with learning how to optimally combine these predictions to produce the final predictions. We selected six base models derived from the linear mixed and Bayesian models widely used in GS. To effectively combine the results of the base models, we used a neural network of a multi-layer perceptron as our meta-model.    

The proposed stacking was applied to open-access resources of rice, maize, barley, and mice. Our analysis compared the performance of the stacking model with that of its constituent base models. To compare the performance of the models, we evaluated quantities such as overfitting and mean squared error (MSE) between observed and predicted phenotype values, which served as a measure of the robustness and prediction accuracy of the models. We also performed the hypothesis tests for prediction accuracy between the proposed model and each base model. We have provided the Python and R scripts with datasets for the readers to reproduce the results discussed in the manuscript.

## Install prerequisites : 
* Python : version 3.6 or later
* R : version 4 or later
* R-packages : rrBLUP, BGLR
* Python-packages : scikit-learn, scipy, tensorflow (version 2.2 or later), keras, numpy, pandas 

## Data sets
1.	preprocessed rice: genotype and phenotype data : https://github.com/infoLab204/stacking/blob/main/data/rice.tar.gz
2.	preprocessed barley : genotype and phenotype data : https://github.com/infoLab204/stacking/blob/main/data/barley.tar.gz
3.	preprocessed maize : genotype and phenotype data: https://github.com/infoLab204/stacking/blob/main/data/maize.tar.gz
4.	preprocessed mice : genotype and phenotype data : https://github.com/infoLab204/stacking/blob/main/data/mice.tar.gz

The original data of the species used can be found in the information below.
* Rice :  http://www.ricediversity.org/data/sets/44kgwas/ 
* Barley : https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0164494#sec019
* Mazie : https://www.panzea.org/data
* Mice : R-packages install BGLR 


## Python and R scripts tutorial

1.	Download the script from the github repository 
    * stacking.py            ## gs start function
    * base_model.R          ## estimating each base model  
    * stacking_base_model.R  ## estimating stacking base model
    * meta_model.py         ## estimating meta model used a neural network of a multi-layer perceptron

2.	Python scripts for estimating predictions using base mode and meta model    
    Usage : python stacking.py genotype_file phenotype_file Abbreviation_phenotype_name
  	
    (ex) GS implementation example for barley phenotype SSW    
        python stacking.py barley_genotype.txt barley_phenotype.txt SSW
  	
    |Abbreviation|Full name|
    |---|---|
    |SSW|Standardized seed weight|
    |F2.2|Weight of seeds with the size less than 2.2|
    |F2.5|Weight of seeds with the size between 2.5 and 2.8|
    |F2.8|Weight of seeds with the size greater than 2.8|
    |PC|Protein content|
    |TW|Test weight|
    |EC|Ergosterol content|
    |PY|Protein yield|

    Output
    1. Base mode predicted value : test and training data
        * rrBLUP_yHat_test or train_predicted.txt
        * gBLUP_yHat_test or train_predicted.txt
        * BA_yHat_test or train_predicted.txt
        * BB_yHat_test or train_predicted.txt
        * BC_yHat_test or train_predicted.txt
        * BL_yHat_test or train_predicted.txt
    2. meta model predicted value : test and training data
        * meta_yHat_test or train_predicted.txt
