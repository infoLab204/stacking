## phenotype data load
pheno <- read.table("y_training_test_full")

n <- as.integer(nrow(pheno)*0.8)
y_test <- pheno$V1[(n+1):nrow(pheno)]
N=length(y_test)

## predicted data load
stacking_test=read.table("meta_yHat_test_predicted.txt")
rrBLUP_test=read.table("rrBLUP_yHat_test_predicted.txt")
BA_test=read.table("BA_yHat_test_predicted.txt")
BB_test=read.table("BB_yHat_test_predicted.txt")
BC_test=read.table("BC_yHat_test_predicted.txt")
BL_test=read.table("BL_yHat_test_predicted.txt")
gBLUP_test=read.table("gBLUP_yHat_test_predicted.txt")

stacking_diff <- abs(y_test-stacking_test)
rrBLUP_diff <- abs(y_test-rrBLUP_test)
BA_diff <- abs(y_test-BA_test)
BB_diff <- abs(y_test-BB_test)
BC_diff <- abs(y_test-BC_test)
BL_diff <- abs(y_test-BL_test)
gBLUP_diff <- abs(y_test-gBLUP_test)

rrBLUP_mean <- mean(rrBLUP_diff$V1)
BA_mean <- mean(BA_diff$V1)
BB_mean <- mean(BB_diff$V1)
BC_mean <- mean(BC_diff$V1)
BL_mean <- mean(BL_diff$V1)
gBLUP_mean <- mean(gBLUP_diff$V1)
stacking_mean <- mean(stacking_diff$V1)

rrBLUP_mean_diff <- rrBLUP_mean-stacking_mean
BA_mean_diff <- BA_mean-stacking_mean
BB_mean_diff <- BB_mean-stacking_mean
BC_mean_diff <- BC_mean-stacking_mean
BL_mean_diff <- BL_mean-stacking_mean
gBLUP_mean_diff <- gBLUP_mean-stacking_mean

rrBLUP_var <- var(rrBLUP_diff$V1)
BA_var <- var(BA_diff$V1)
BB_var <- var(BB_diff$V1)
BC_var <- var(BC_diff$V1)
BL_var <- var(BL_diff$V1)
gBLUP_var <- var(gBLUP_diff$V1)
stacking_var <- var(stacking_diff$V1)

rrBLUP_power=sqrt(N/2)*(abs(rrBLUP_mean-stacking_mean)/sqrt((rrBLUP_var+stacking_var)/2))-qnorm(0.975)
BA_power=sqrt(N/2)*(abs(BA_mean-stacking_mean)/sqrt((BA_var+stacking_var)/2))-qnorm(0.975)
BB_power=sqrt(N/2)*(abs(BB_mean-stacking_mean)/sqrt((BB_var+stacking_var)/2))-qnorm(0.975)
BC_power=sqrt(N/2)*(abs(BC_mean-stacking_mean)/sqrt((BC_var+stacking_var)/2))-qnorm(0.975)
BL_power=sqrt(N/2)*(abs(BL_mean-stacking_mean)/sqrt((BL_var+stacking_var)/2))-qnorm(0.975)
gBLUP_power=sqrt(N/2)*(abs(gBLUP_mean-stacking_mean)/sqrt((gBLUP_var+stacking_var)/2))-qnorm(0.975)

power_test <- c(1:6)
power_test[1] <- pnorm(rrBLUP_power)
power_test[2] <- pnorm(gBLUP_power)
power_test[3] <- pnorm(BA_power)
power_test[4] <- pnorm(BB_power)
power_test[5] <- pnorm(BC_power)
power_test[6] <- pnorm(BL_power)

print("power test : ")
paste(power_test[1],power_test[2],power_test[3],power_test[4],power_test[5],power_test[6])
write.table(paste(power_test[1],power_test[2],power_test[3],power_test[4],power_test[5],power_test[6]),file=paste0("power_test.txt"),row.names=F,col.names=F,quote=F, append=F,sep="\t")

common<- 2*((qnorm(0.975)+qnorm(0.80))**2)

rrBLUP_n=common*((rrBLUP_var+stacking_var)/2)/((rrBLUP_mean-stacking_mean)**2)
BA_n=common*((BA_var+stacking_var)/2)/((BA_mean-stacking_mean)**2)
BB_n=common*((BB_var+stacking_var)/2)/((BB_mean-stacking_mean)**2)
BC_n=common*((BC_var+stacking_var)/2)/((BC_mean-stacking_mean)**2)
BL_n=common*((BL_var+stacking_var)/2)/((BL_mean-stacking_mean)**2)
gBLUP_n=common*((gBLUP_var+stacking_var)/2)/((gBLUP_mean-stacking_mean)**2)

print("required sample size : ")

sample_size <- c(1:6)
sample_size[1] <- rrBLUP_n
sample_size[2] <- gBLUP_n
sample_size[3] <- BA_n
sample_size[4] <- BB_n
sample_size[5] <- BC_n
sample_size[6] <- BL_n

paste(sample_size[1],sample_size[2],sample_size[3],sample_size[4],sample_size[5],sample_size[6])
write.table(paste(sample_size[1],sample_size[2],sample_size[3],sample_size[4],sample_size[5],sample_size[6]),file="required_sample_size.txt",row.names=F,col.names=F,quote=F, append=F,sep="\t")

## margin
rrBLUP_delta <- (rrBLUP_mean-stacking_mean)+sqrt((rrBLUP_var+stacking_var)/N)*(qnorm(0.975)+qnorm(0.8))
BA_delta <- (BA_mean-stacking_mean)+sqrt((BA_var+stacking_var)/N)*(qnorm(0.975)+qnorm(0.8))
BB_delta <- (BB_mean-stacking_mean)+sqrt((BB_var+stacking_var)/N)*(qnorm(0.975)+qnorm(0.8))
BC_delta <- (BC_mean-stacking_mean)+sqrt((BC_var+stacking_var)/N)*(qnorm(0.975)+qnorm(0.8))
BL_delta <- (BL_mean-stacking_mean)+sqrt((BL_var+stacking_var)/N)*(qnorm(0.975)+qnorm(0.8))
gBLUP_delta <-(gBLUP_mean-stacking_mean)+sqrt((gBLUP_var+stacking_var)/N)*(qnorm(0.975)+qnorm(0.8))

delta_margin <- c(1:6)
delta_margin[1] <- rrBLUP_delta
delta_margin[2] <- gBLUP_delta
delta_margin[3] <- BA_delta
delta_margin[4] <- BB_delta
delta_margin[5] <- BC_delta
delta_margin[6] <- BL_delta

print("margin : ")
paste0(rrBLUP_delta," ",BA_delta," ", BB_delta," ", BC_delta, " ", BL_delta, " ", gBLUP_delta)
write.table(paste(delta_margin[1],delta_margin[2],delta_margin[3],delta_margin[4],delta_margin[5],delta_margin[6]),file="margin.txt",row.names=F,col.names=F,quote=F, append=F,sep="\t")


rrBLUP <- wilcox.test(rrBLUP_diff$V1,stacking_diff$V1, alternative="greater", conf.int=TRUE, mu=-rrBLUP_delta, paired=TRUE)
BA <- wilcox.test(BA_diff$V1,stacking_diff$V1, alternative="greater", conf.int=TRUE, mu=-BA_delta, paired=TRUE)
BB <- wilcox.test(BB_diff$V1,stacking_diff$V1, alternative="greater", conf.int=TRUE, mu=-BB_delta, paired=TRUE)
BC <- wilcox.test(BC_diff$V1,stacking_diff$V1, alternative="greater", conf.int=TRUE, mu=-BC_delta, paired=TRUE)
BL <- wilcox.test(BL_diff$V1,stacking_diff$V1, alternative="greater", conf.int=TRUE, mu=-BL_delta, paired=TRUE)
gBLUP <- wilcox.test(gBLUP_diff$V1,stacking_diff$V1, alternative="greater", conf.int=TRUE, mu=-gBLUP_delta, paired=TRUE)

pvalue <-c(1:6)
pvalue[1] <- rrBLUP$p.value
pvalue[2] <- BA$p.value
pvalue[3] <- BB$p.value
pvalue[4] <- BC$p.value
pvalue[5] <- BL$p.value
pvalue[6] <- gBLUP$p.value

print("noninferiority test p-value")
paste(pvalue[1],pvalue[2],pvalue[3],pvalue[4],pvalue[5],pvalue[6])
write.table(paste(pvalue[1],pvalue[2],pvalue[3],pvalue[4],pvalue[5],pvalue[6]),file="noninferiority_test_pvalue.txt",row.names=F,col.names=F,quote=F, append=F,sep="\t")

# end of hypothesis_tests.R
