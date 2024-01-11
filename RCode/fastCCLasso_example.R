#!/usr/bin/env Rscript
#---------------------------------------------------------------
source("fastCCLasso.R")

method_name <- c("fastCCLasso","SparCC","CCLasso","COAT");
method_num <- length(method_name);
set.seed(1212); 

filein = "example_input.csv";
fileout = "example_output.csv";

OTUdata <- read.csv(filein, header = TRUE,row.names = 1);
OTUdata <- OTUdata + 0.5;
total_count <- rowSums(OTUdata);
xMat <- (1/total_count) * OTUdata;  # compositional data
n <- dim(xMat)[1];
p <- dim(xMat)[2];
cor_lower <- matrix(0,p*(p-1)/2,method_num);
colnames(cor_lower) <- method_name;
for(i in 1:method_num){
  result.i <- Callallmethods(method=method_name[i], xMat=xMat,cv_k=3,Edge_eps=0.05);
  cor_lower[,i] <- result.i$est_lower;
}
vnames <- colnames(OTUdata)
lower.ind <- which(lower.tri(matrix(0,p,p)), arr.ind = TRUE)
CorEst <- data.frame(variable_1=vnames[lower.ind[,1]],
                     variable_2=vnames[lower.ind[,2]],
                     cor_lower)
write.table(CorEst, file=fileout, row.names=F,sep = ",")
