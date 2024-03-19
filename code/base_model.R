library(BGLR)
library(rrBLUP)

SNP <- read.table('X_training_test_full',sep="\t", header=F)
Y0 <- read.table('y_training_test_full',sep="\t", header=F)

n<-nrow(SNP)
train_n=as.integer(n*0.8)
p<-ncol(SNP)
SNP <- SNP[,1:(p-1)]
test <-c((train_n+1):n)
y<-Y0[,1]
yNA <- y
test <-c((train_n+1):n)
train_idx <-setdiff(1:n, test)
yNA[test] <-NA

test_file_name="y_test_rrBLUP.txt"
train_file_name="y_training_rrBLUP.txt"

## Linear Mixed Model(LMM)
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

rrblup_predicted_train <- geno_train %*% marker_effects
rrblup_predicted_train_result <-as.vector(rrblup_predicted_train[,1])+rrBLUE

write.table(rrblup_predicted_result,file=test_file_name,row.names=F,col.names=F,quote=F, append=F,sep="\n")
write.table(rrblup_predicted_train_result,file=train_file_name, row.names=F,col.names=F,quote=F, append=F,sep="\n")

# BayesA
ETA <-list(list(X=SNP, model="BayesA"))
fmBA <-BGLR(y=yNA, ETA=ETA, nIter=5000, burnIn=1000, saveAt="BA_")
write.table(fmBA$yHat[test],file="y_test_BA.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")
write.table(fmBA$yHat[-test],file="y_training_BA.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")

## BayesB
ETA <-list(list(X=SNP, model="BayesB"))
fmBB <-BGLR(y=yNA, ETA=ETA, nIter=5000, burnIn=1000, saveAt="BB_")
write.table(fmBB$yHat[test],file="y_test_BB.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")
write.table(fmBB$yHat[-test],file="y_training_BB.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")

## BayesC
ETA <-list(list(X=SNP, model="BayesC"))
fmBC <-BGLR(y=yNA, ETA=ETA, nIter=5000, burnIn=1000, saveAt="BC_")
write.table(fmBC$yHat[test],file="y_test_BC.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")
write.table(fmBC$yHat[-test],file="y_train_BC.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")

## Basyes LASSO
ETA <-list(list(X=SNP, model="BL"))
fmBL <-BGLR(y=yNA, ETA=ETA, nIter=5000, burnIn=1000, saveAt="BL_")
write.table(fmBL$yHat[test],file="y_test_BL.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")
write.table(fmBL$yHat[-test],file="y_training_BL.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")

## GBLUP
ETA <-list(list(X=SNP, model="BRR"))
fmBRR <-BGLR(y=yNA, ETA=ETA, nIter=5000, burnIn=1000, saveAt="gBLUP_")
write.table(fmBRR$yHat[test],file="y_test_gBLUP.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")
write.table(fmBRR$yHat[-test],file="y_training_gBLUP.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")
## end
