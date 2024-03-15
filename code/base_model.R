library(BGLR)
library(rrBLUP)

SNP <- read.table('X_tr_output_full',sep="\t", header=F)
Y0 <- read.table('y_tr_output_full',sep="\t", header=F)

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

test_file_name="rrBLUP_yHat_test_predicted_each.txt"
train_file_name="rrBLUP_yHat_train_predicted_each.txt"

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
fmBA <-BGLR(y=yNA, ETA=ETA, nIter=5000, burnIn=1000, saveAt="BayesA_")
write.table(fmBA$yHat[test],file="BayesA_yHat_test_predicted_each.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")
write.table(fmBA$yHat[-test],file="BayesA_yHat_train_predicted_each.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")

BAco<-cor(fmBA$yHat[test],y[test])
BAco1<-cor(fmBA$yHat[-test],y[-test])

## bayesB
ETA <-list(list(X=SNP, model="BayesB"))
fmBB <-BGLR(y=yNA, ETA=ETA, nIter=5000, burnIn=1000, saveAt="BayesB_")
write.table(fmBB$yHat[test],file="BayesB_yHat_test_predicted_each.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")
write.table(fmBB$yHat[-test],file="BayesB_yHat_train_predicted_each.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")

BBco<-cor(fmBB$yHat[test],y[test])
BBco1<-cor(fmBB$yHat[-test],y[-test])

## BayesC
ETA <-list(list(X=SNP, model="BayesC"))
fmBC <-BGLR(y=yNA, ETA=ETA, nIter=5000, burnIn=1000, saveAt="BayesC_")
write.table(fmBC$yHat[test],file="BayesC_yHat_test_predicted_each.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")
write.table(fmBC$yHat[-test],file="BayesC_yHat_train_predicted_each.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")

BCco<-cor(fmBC$yHat[test],y[test])
BCco1<-cor(fmBC$yHat[-test],y[-test])

## Basyes LASSO
ETA <-list(list(X=SNP, model="BL"))
fmBL <-BGLR(y=yNA, ETA=ETA, nIter=5000, burnIn=1000, saveAt="BL_")
write.table(fmBL$yHat[test],file="BL_yHat_test_predicted_each.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")
write.table(fmBL$yHat[-test],file="BL_yHat_train_predicted_each.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")

BLco<-cor(fmBL$yHat[test],y[test])
BLco1<-cor(fmBL$yHat[-test],y[-test])

## GBLUP
ETA <-list(list(X=SNP, model="BRR"))
fmBRR <-BGLR(y=yNA, ETA=ETA, nIter=5000, burnIn=1000, saveAt="BL_")
write.table(fmBRR$yHat[test],file="BRR_yHat_test_predicted_each.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")
write.table(fmBRR$yHat[-test],file="BRR_yHat_train_predicted_each.txt",row.names=F,col.names=F,quote=F, append=F,sep="\n")

BRRco<-cor(fmBRR$yHat[test],y[test])
BRRco1<-cor(fmBRR$yHat[-test],y[-test])

cor_save="base_model_correlation.txt"

paste("rrBLUP test", rrBLUPco, "rrBLUP train",rrBLUPco1)
paste("BayesA test", BAco,"train", BAco1)
paste("BayesB test", BBco,"train", BBco1)
paste("BayesC test", BCco,"train", BCco1)
paste("Bayes_LASSO test", BLco,"train", BLco1)
paste("Bayes_GBLUP test", BRRco,"train", BRRco1)

write.table(paste("rrBLUP test", rrBLUPco, "train",rrBLUPco1),file=cor_save,append=FALSE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\n")
write.table(paste("BayesA test", BAco,"train", BAco1),file=cor_save,append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\n")
write.table(paste("BayesB test", BBco,"train", BBco1),file=cor_save,append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\n")
write.table(paste("BayesC test", BCco,"train", BCco1),file=cor_save,append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\n")
write.table(paste("Bayes_LASSO test", BLco,"train", BLco1),file=cor_save,append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\n")
write.table(paste("Bayes_GBLUP test", BRRco,"train", BRRco1),file=cor_save,append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\n")
