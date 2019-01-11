## Do the DeLong's test on HE p-values against the semi-parallel and the gold standard
library(pROC)
## Load the ROC library

S_list = read.table("snpMat.txt", header=T)
S = matrix(unlist(S_list, use.names=F), ncol = length(S_list), byrow=F)
covariate = read.csv("covariates.csv", header=T, sep=",")
## Read the SNP and covariates data

y = unlist(covariate[2], use.names=F)
rawX0 = covariate[3:5]

normalize <- function(x){(x-min(x))/(max(x)-min(x))}

X0_list = lapply(rawX0, function(x) ifelse(is.na(x), mean(x, na.rm=T), x))
X0 = (matrix(unlist(X0_list), ncol = length(X0_list), byrow=F))
X0 = normalize(X0)
X = cbind(1, X0)
n = nrow(S)
m = ncol(S)
k = ncol(X0)

mod0 = glm(y ~ X-1, family = binomial("logit"))
p = mod0$fitted
w = p * (1 - p)
z = log(p / (1 - p)) + (y - p) / (p * (1 - p))
xtw = t(X * w)
U1 = xtw %*% z
U2 = solve(xtw %*% X, U1)
ztr = z  - X %*% U2
U3 = xtw %*% S
U4 = solve(xtw %*% X, U3)
Str = S - X %*% U4
Str2 = colSums(w * Str^2)
b = crossprod(ztr * w, Str)/Str2
err = sqrt(1/ Str2)
pvalue_parallel = 2 * pnorm(-abs(b / err))
pval_sp = as.numeric(as.matrix(pvalue_parallel))
## Perform the semi-parallel log-regression

result_he = read.table("HE_Pvals.txt",header = F)
pval_he = matrix(unlist(result_he[2]))
## Read the p-value of the HE solution

pval = matrix(pval_he,ncol = 1)
## Working variable on p-values


b1 = NULL
err1 = NULL
pval_gs = NULL
for(i in 1:m){
  mod1 = glm(y ~ S[, i] + X0, family = binomial ("logit"))
  b1[i] = summary(mod1)$coef[2, 1]
  err1[i] = summary(mod1)$coef[2, 2]
  pval_gs[i] = summary(mod1)$coef[2, 4]
}
## Perform the gold-standard log-regression

num_sample = 1499
num_round = 10
num_prog = ncol(pval)
num_comp = num_prog*(num_prog-1)/2
num_sample = 1499

roc_list_sp<-matrix(list(), nrow = 4, ncol = num_round)
comparison_sp<-matrix(list(), nrow = 4, ncol = num_round)
pval_sp<-matrix(c(pval, pval_sp), ncol=num_prog+1)

tp_sp<-array(0,c(num_prog,4,num_round))
prec_sp<-array(0,c(num_prog,4,num_round))
f1_sp<-array(0,c(num_prog,4,num_round))
roc_list_sp<-matrix(list(), nrow = 4, ncol = num_round)

for (i in 2:5) {
  threshold = 10^(-i)
  lgs = 1*(pval_gs<=threshold)
  pval_eval_sp = cbind(pval_sp, lgs)
  for (j in 1:num_round) {
    roc_sp<-list()
    round_sp<-list()
    pval_sample_sp = pval_eval_sp[sample(nrow(pval_sp), num_sample), ]
    ind_gs2 = sum(pval_sample_sp[,num_prog+2])
    while (ind_gs2>=num_sample | ind_gs2<=0) {
      pval_sample_sp = pval_eval_sp[sample(nrow(pval_sp), num_sample), ]
      ind_gs2 = sum(pval_sample_sp[,num_prog+2])  
    }
    for (k in 1:num_prog) {
      tp_sp[k,i-1,j] = length(which(pval_sample_sp[which(pval_sample_sp[,num_prog+1]<=threshold),k]<=threshold))/length(which(pval_sample_sp[,num_prog+1]<=threshold))
      prec_sp[k,i-1,j] = length(which(pval_sample_sp[which(pval_sample_sp[,num_prog+1]<=threshold),k]<=threshold))/length(which(pval_sample_sp[,k]<=threshold))
      f1_sp[k,i-1,j] = 2 * tp_sp[k,i-1,j] *  prec_sp[k,i-1,j] / (prec_sp[k,i-1,j]+tp_sp[k,i-1,j])      
    }
    for (k in 1:(num_prog+1)) {
      roc_sp[[k]]<-roc(pval_sample_sp[,num_prog+2], pval_sample_sp[,k])
    }
    roc_list_sp[[i-1,j]]<-roc_sp
    for (l in 1:num_prog) {
      round_sp[[l]]<-roc.test(roc_sp[[l]], roc_sp[[num_prog+1]])
    }
    comparison_sp[[i-1,j]]<-round_sp
  }
}
## comparison_sp contains all the ROC test data

pstat_sp = matrix(0, nrow = num_round, ncol = num_prog)
pmean_sp = matrix(0, nrow = 4, ncol = num_prog)
pdev_sp = matrix(0, nrow = 4, ncol = num_prog)
pstat_sp = matrix(0, nrow = num_round, ncol = num_prog)
## Prepare for collecting the statistics

for (i in 1:4) {
  for (j in 1:num_round) {
    for (k in 1:num_prog) {
      pstat_sp[j,k] = comparison_sp[[i,j]][[k]]$p.value
    }
  }
  for (k in 1:num_prog) {
    pdev_sp[i,k] = sd(pstat_sp[,k])
  }
  pmean_sp[i,] = colMeans(pstat_sp)
}

## Now we should have results like the following:
## > pmean_sp
##           [,1]
## [1,] 0.4038009
## [2,] 0.5357250
## [3,] 0.6404474
## [4,] 0.8959000
## > pdev_sp
##           [,1]
## [1,] 0.3001436
## [2,] 0.2704072
## [3,] 0.2638097
## [4,] 0.2194620

