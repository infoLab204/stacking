library(BGLR)
library(rrBLUP)

SNP <- read.table('X_training_test_full',sep="\t", header=F)
Y0 <- read.table('y_training_test_full',sep="\t", header=F)

n_fold=5
p<-ncol(SNP)
print(p)

SNP <- SNP[,1:(p-1)]

n<-nrow(SNP)
train_n=as.integer(n*0.8)
cat(n, train_n)

SNP0 <- SNP[1:train_n,]
print(nrow(SNP0))
y<-Y0[1:train_n,1]
print(y)

split_n=as.integer(train_n/n_fold)
print(split_n)

rrBLUP_yHat<-c(1:train_n)
BA_yHat<-c(1:train_n)
BB_yHat<-c(1:train_n)
BC_yHat<-c(1:train_n)
BL_yHat<-c(1:train_n)
BRR_yHat<-c(1:train_n)

test_idx=as.integer(seq(1,train_n, by=split_n))

for(i in 1:n_fold) {
    if(i<n_fold) test<-c(test_idx[i]:(test_idx[i+1]-1))
    else test<-c(test_idx[i]:train_n)

    ## rrBLUP
    pheno_train <- y[-test]
    pheno_test <- y[test]
    geno_train <- as.matrix(SNP0[-test,])
    geno_test <- as.matrix(SNP0[test,])

    rrblupmodel <- mixed.solve(y=pheno_train, Z=geno_train, K=NULL, X=NULL)
    marker_effects <- rrblupmodel$u
    rrBLUE <- as.vector(rrblupmodel$beta)

    rrblup_predicted_test <- geno_test %*% marker_effects
    rrblup_predicted_result <-as.vector(rrblup_predicted_test[,1])+rrBLUE

    rrblup_predicted_train <- geno_train %*% marker_effects
    rrblup_predicted_train_result <-as.vector(rrblup_predicted_train[,1])+rrBLUE
    
    rrBLUP_yHat[test] <- rrblup_predicted_result

    write.table(rrBLUP_yHat,file="rrBLUP_yHat_train_stacking_predicted.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")


    yNA <- y
    yNA[test] <-NA

    ## BayesA
    ETA <-list(list(X=SNP0, model="BayesA"))
    saveAt=paste0("BayesA_",i,"_")
    fmBA <-BGLR(y=yNA, ETA=ETA, nIter=5000, burnIn=1000, saveAt=saveAt)
    BA_yHat[test]=fmBA$yHat[test]

    write.table(BA_yHat,file="BA_yHat_training_stacking_predicted.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")

  
    ## BayesB
    ETA <-list(list(X=SNP0, model="BayesB"))
    saveAt=paste0("BayesB_",i,"_")
    fmBB <-BGLR(y=yNA, ETA=ETA, nIter=5000, burnIn=1000, saveAt=saveAt)
    BB_yHat[test]=fmBB$yHat[test]

    write.table(BB_yHat,file="BB_yHat_training_stacking_predicted.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")


    ## BayesC
    ETA <-list(list(X=SNP0, model="BayesC"))
    saveAt=paste0("BayesC_",i,"_")
    fmBC <-BGLR(y=yNA, ETA=ETA, nIter=5000, burnIn=1000, saveAt=saveAt)
    BC_yHat[test]=fmBC$yHat[test]

    write.table(BC_yHat,file="BC_yHat_training_stacking_predicted.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")


    ## Bayes LASSO
    ETA <-list(list(X=SNP0, model="BL"))
    saveAt=paste0("BL_",i,"_")
    fmBL <-BGLR(y=yNA, ETA=ETA, nIter=5000, burnIn=1000, saveAt=saveAt)
    BL_yHat[test]=fmBL$yHat[test]

    write.table(BL_yHat,file="BL_yHat_training_stacking_predicted.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")

    ## Bayes GBLUP
    ETA <-list(list(X=SNP0, model="BRR"))
    saveAt=paste0("BRR_",i,"_")
    fmBRR <-BGLR(y=yNA, ETA=ETA, nIter=5000, burnIn=1000, saveAt=saveAt)
    BRR_yHat[test]=fmBRR$yHat[test]

    write.table(BRR_yHat,file="gBLUP_yHat_training_stacking_predicted.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")

}
print("train learining END")


## rrBLUP test set predict
rrBLUP_yHat_test<-c(1:train_n)
test <-c((train_n+1):n)
y<-c(1:n)
train_idx=seq(1,train_n)
y<-Y0[,1]

pheno_train <- y[-test]
pheno_test <- y[test]
geno_train <- as.matrix(SNP[-test,])
geno_test <- SNP[test,]
geno_test <- as.matrix(geno_test)

rrblupmodel <- mixed.solve(y=pheno_train, Z=geno_train, K=NULL, X=NULL)
marker_effects <- rrblupmodel$u
rrBLUE <- as.vector(rrblupmodel$beta)

rrblup_predicted_test <- geno_test %*% marker_effects
rrblup_predicted_result <-as.vector(rrblup_predicted_test[,1])+rrBLUE

rrBLUP_yHat_test[test] <- rrblup_predicted_result
write.table(rrBLUP_yHat_test[test],file="rrBLUP_yHat_test_stacking_predicted.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")

test <-c((train_n+1):n)
y<-Y0[,1]
yNA <- y
yNA[test] <-NA

## BayesA test set predict
ETA <-list(list(X=SNP, model="BayesA"))
fmBA <-BGLR(y=yNA, ETA=ETA, nIter=5000, burnIn=1000, saveAt="BayesA_")
write.table(fmBA$yHat[test],file="BA_yHat_test_stacking_predicted.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")

## BayesB test set predict
ETA <-list(list(X=SNP, model="BayesB"))
fmBB <-BGLR(y=yNA, ETA=ETA, nIter=5000, burnIn=1000, saveAt="BayesB_")
write.table(fmBB$yHat[test],file="BB_yHat_test_stacking_predicted.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")

## BayesC test set predict
ETA <-list(list(X=SNP, model="BayesC"))
fmBC <-BGLR(y=yNA, ETA=ETA, nIter=5000, burnIn=1000, saveAt="BayesC_")
write.table(fmBC$yHat[test],file="BC_yHat_test_stacking_predicted.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")

## Bayes Lasso test set predict
ETA <-list(list(X=SNP, model="BL"))
fmBL <-BGLR(y=yNA, ETA=ETA, nIter=5000, burnIn=1000, saveAt="BL_")
write.table(fmBL$yHat[test],file="BL_yHat_test_stacking_predicted.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")

## gBLUP test set predict
ETA <-list(list(X=SNP, model="BRR"))
fmBRR <-BGLR(y=yNA, ETA=ETA, nIter=5000, burnIn=1000, saveAt="BRR_")
write.table(fmBRR$yHat[test],file="gBLUP_yHat_test_stacking_predicted.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")
