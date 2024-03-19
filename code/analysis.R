y=read.table("y_training_test_full")

n <- as.integer(nrow(y)*0.8) 
y_training <- y$V1[1:n]  ## training
y_test <- y$V1[(n+1):nrow(y)]  ## test

stacking_test=read.table("meta_yHat_test_predicted.txt")
rrBLUP_test=read.table("rrBLUP_yHat_test_predicted.txt")
BA_test=read.table("BA_yHat_test_predicted.txt")
BB_test=read.table("BB_yHat_test_predicted.txt")
BC_test=read.table("BC_yHat_test_predicted.txt")
BL_test=read.table("BL_yHat_test_predicted.txt")
gBLUP_test=read.table("gBLUP_yHat_test_predicted.txt")


## correlation 
correlation=c(1:7)
correlation[1]<-cor(y_test, rrBLUP_test$V1)
correlation[2]<-cor(y_test, gBLUP_test$V1)
correlation[3]<-cor(y_test, BA_test$V1)
correlation[4]<-cor(y_test, BB_test$V1)
correlation[5]<-cor(y_test, BC_test$V1)
correlation[6]<-cor(y_test, BL_test$V1)
correlation[7]<-cor(y_test, stacking_test$V1)
paste("rrBLUP","gBLUP","BA","BB","BC","BL","Stacking")
print(correlation)
write.table(paste("rrBLUP","gBLUP","BA","BB","BC","BL","Stacking"),file="correlation.txt",row.names=F,col.names=F,quote=F, append=F,sep="\t")
write.table(paste(correlation[1],correlation[2],correlation[3],correlation[4],correlation[5],correlation[6], correlation[7]),file="correlation.txt",row.names=F,col.names=F,quote=F, append=T,sep="\t")

## test mse
rrBLUP_sum <- sum((y_test-rrBLUP_test)**2)
gBLUP_sum <- sum((y_test-gBLUP_test)**2)
BA_sum <- sum((y_test-BA_test)**2)
BB_sum <- sum((y_test-BB_test)**2)
BC_sum <- sum((y_test-BC_test)**2)
BL_sum <- sum((y_test-BL_test)**2)
stacking_sum <- sum((y_test-stacking_test)**2)

rrBLUP <- rrBLUP_sum/length(y_test)
gBLUP <- gBLUP_sum/length(y_test)
BA <- BA_sum/length(y_test)
BB <- BB_sum/length(y_test)
BC <- BC_sum/length(y_test)
BL <- BL_sum/length(y_test)
stacking <- stacking_sum/length(y_test)


test_mse=c(1:7)
test_mse[1]<-rrBLUP
test_mse[2]<-gBLUP
test_mse[3]<-BA
test_mse[4]<-BB
test_mse[5]<-BC
test_mse[6]<-BL
test_mse[7]<-stacking

paste(test_mse)
write.table(paste("rrBLUP","gBLUP","BA","BB","BC","BL","Stacking"),file="MSE_test.txt",row.names=F,col.names=F,quote=F, append=F,sep="\t")
write.table(paste(test_mse[1],test_mse[2],test_mse[3],test_mse[4],test_mse[5],test_mse[6], test_mse[7]),file="MSE_test.txt",row.names=F,col.names=F,quote=F, append=T,sep="\t")



## training mse
stacking_training=read.table("meta_yHat_training_predicted.txt")
rrBLUP_training=read.table("rrBLUP_yHat_training_predicted.txt")
BA_training=read.table("BA_yHat_training_predicted.txt")
BB_training=read.table("BB_yHat_training_predicted.txt")
BC_training=read.table("BC_yHat_training_predicted.txt")
BL_training=read.table("BL_yHat_training_predicted.txt")
gBLUP_training=read.table("gBLUP_yHat_training_predicted.txt")

rrBLUP_sum <- sum((y_training-rrBLUP_training)**2)
gBLUP_sum <- sum((y_training-gBLUP_training)**2)
BA_sum <- sum((y_training-BA_training)**2)
BB_sum <- sum((y_training-BB_training)**2)
BC_sum <- sum((y_training-BC_training)**2)
BL_sum <- sum((y_training-BL_training)**2)
stacking_sum <- sum((y_training-stacking_training)**2)

rrBLUP <- rrBLUP_sum/length(y_training)
gBLUP <- gBLUP_sum/length(y_training)
BA <- BA_sum/length(y_training)
BB <- BB_sum/length(y_training)
BC <- BC_sum/length(y_training)
BL <- BL_sum/length(y_training)
stacking <- stacking_sum/length(y_training)


training_mse=c(1:7)
training_mse[1]<-rrBLUP
training_mse[2]<-gBLUP
training_mse[3]<-BA
training_mse[4]<-BB
training_mse[5]<-BC
training_mse[6]<-BL
training_mse[7]<-stacking
paste("rrBLUP","gBLUP","BA","BB","BC","BL","Stacking")
print(training_mse)
write.table(paste("rrBLUP","gBLUP","BA","BB","BC","BL","Stacking"),file="MSE_training.txt",row.names=F,col.names=F,quote=F, append=F,sep="\t")
write.table(paste(training_mse[1],training_mse[2],training_mse[3],training_mse[4],training_mse[5],training_mse[6], training_mse[7]),file="MSE_training.txt",row.names=F,col.names=F,quote=F, append=T,sep="\t")

## overfitting
overfitting=abs(test_mse-training_mse)
print(overfitting)
write.table(paste("rrBLUP","gBLUP","BA","BB","BC","BL","Stacking"),file="overfitting.txt",row.names=F,col.names=F,quote=F, append=F,sep="\t")
write.table(paste(overfitting[1],overfitting[2],overfitting[3],overfitting[4],overfitting[5],overfitting[6], overfitting[7]),file="overfitting.txt",row.names=F,col.names=F,quote=F, append=T,sep="\t")

