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

test_file_name="rrBLUP_yHat_test_predicted.txt"
train_file_name="rrBLUP_yHat_training_predicted.txt"

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

rrBLUPco=cor(rrblup_predicted_result, pheno_test)
rrBLUPco1=cor(rrblup_predicted_train_result, pheno_train)

write.table(rrblup_predicted_result,file=test_file_name,row.names=F,col.names=F,quote=F, append=F,sep="\n")
write.table(rrblup_predicted_train_result,file=train_file_name, row.names=F,col.names=F,quote=F, append=F,sep="\n")

# Bayesian methods
# BayesA
ETA <-list(list(X=SNP, model="BayesA"))
fmBA <-BGLR(y=yNA, ETA=ETA, nIter=5000, burnIn=1000, saveAt="BA_")
write.table(fmBA$yHat[test],file="BA_yHat_test_predicted.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")
write.table(fmBA$yHat[-test],file="BA_yHat_training_predicted.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")

BAco<-cor(fmBA$yHat[test],y[test])
BAco1<-cor(fmBA$yHat[-test],y[-test])

## bayesB
ETA <-list(list(X=SNP, model="BayesB"))
fmBB <-BGLR(y=yNA, ETA=ETA, nIter=5000, burnIn=1000, saveAt="BB_")
write.table(fmBB$yHat[test],file="BB_yHat_test_predicted.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")
write.table(fmBB$yHat[-test],file="BB_yHat_training_predicted.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")

BBco<-cor(fmBB$yHat[test],y[test])
BBco1<-cor(fmBB$yHat[-test],y[-test])

## BayesC
ETA <-list(list(X=SNP, model="BayesC"))
fmBC <-BGLR(y=yNA, ETA=ETA, nIter=5000, burnIn=1000, saveAt="BC_")
write.table(fmBC$yHat[test],file="BC_yHat_test_predicted.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")
write.table(fmBC$yHat[-test],file="BC_yHat_train_predicted.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")

BCco<-cor(fmBC$yHat[test],y[test])
BCco1<-cor(fmBC$yHat[-test],y[-test])

## Basyes LASSO
ETA <-list(list(X=SNP, model="BL"))
fmBL <-BGLR(y=yNA, ETA=ETA, nIter=5000, burnIn=1000, saveAt="BL_")
write.table(fmBL$yHat[test],file="BL_yHat_test_predicted.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")
write.table(fmBL$yHat[-test],file="BL_yHat_training_predicted.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")

BLco<-cor(fmBL$yHat[test],y[test])
BLco1<-cor(fmBL$yHat[-test],y[-test])

## GBLUP
ETA <-list(list(X=SNP, model="BRR"))
fmBRR <-BGLR(y=yNA, ETA=ETA, nIter=5000, burnIn=1000, saveAt="gBLUP_")
write.table(fmBRR$yHat[test],file="gBLUP_yHat_test_predicted.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")
write.table(fmBRR$yHat[-test],file="gBLUP_yHat_training_predicted.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")

BRRco<-cor(fmBRR$yHat[test],y[test])
BRRco1<-cor(fmBRR$yHat[-test],y[-test])

cor_save="base_model_correlation.txt"

paste("rrBLUP test", rrBLUPco, "training",rrBLUPco1)
paste("gBLUP test", BRRco,"training", BRRco1)
paste("BA test", BAco,"training", BAco1)
paste("BB test", BBco,"training", BBco1)
paste("BC test", BCco,"training", BCco1)
paste("BL test", BLco,"training", BLco1)

write.table(paste("rrBLUP test", rrBLUPco, "training",rrBLUPco1),file=cor_save,append=FALSE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\n")
write.table(paste("gBLUP test", BRRco,"training", BRRco1),file=cor_save,append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\n")
write.table(paste("BA test", BAco,"training", BAco1),file=cor_save,append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\n")
write.table(paste("BB test", BBco,"training", BBco1),file=cor_save,append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\n")
write.table(paste("BC test", BCco,"training", BCco1),file=cor_save,append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\n")
write.table(paste("BL test", BLco,"training", BLco1),file=cor_save,append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\n")
