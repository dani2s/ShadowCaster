#! /usr/bin/env Rscript

args<-commandArgs(TRUE)
# Mean of orthologs identity distributions
ortho <- read.csv(args[1], sep = "\t", header = TRUE)
# Estimated by the method of the moments
means = ortho$mean

# Standard deviation of orthologs identity distributions
# Estimated by the method of the moments
stds = ortho$std

# Vector of orthology probability 
p_orth = ortho$Pr_orthologs

# Read in dataset
#CAMBIO IDENT 0 -- 0.5!!
dat = read.csv(args[2], sep = "\t", header = FALSE)

dimDataAxis <- dim.data.frame(dat)[1]
dimDataRows <- dim.data.frame(dat)[2]
L = vector(length = dimDataAxis)
for (j in 1:dimDataAxis){
  iden <- dat[j,2:dimDataRows];
iden <- as.vector(iden);
probs = vector(length = dimDataRows-1)
for (i in 1:length(iden)){
 if (iden[i] > 0.5){
    probs[i] = dnorm(means[i], means[i], stds[i]) * p_orth[i]
  }
  else{
    probs[i] = 1 - p_orth[i]  
  }
}
L[j] <- sum(log(probs))
}
#print(L)
#plot(L)

qHgt = quantile(L, probs = c(0.2))
png(filename = args[3], width = 800, height = 700)
hist(L, main = "alien CDS")
dev.off()

test <- list('id'=dat$V1, 'log'= L)
outResult <- as.data.frame(test)
write.csv(outResult, file = args[4])

