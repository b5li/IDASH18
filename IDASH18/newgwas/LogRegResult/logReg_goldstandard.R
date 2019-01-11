# Check the speed of semi-parallel approach
# for logistic regression with covariates
# Karolina Sikorska and Paul Eilers, 2013

#####sample dataset from the paper#########
#set.seed(2013)
#n = 10000
#m = 1000
#k = 1
#S = matrix(2 * runif(n * m), n, m)
#X0 = matrix(rnorm(n * k), n, k)
#X = cbind(1, X0)
#y = rbinom(n, size = 1, prob = c(0.5, 0.5))
#####sample dataset from the paper##########


######import dataset extracted from PGP########
#SNPs
S_list = read.table("snpMat_forEval.txt", header=T)
S = matrix(unlist(S_list, use.names=F), ncol = length(S_list), byrow=F)

covariate = read.csv("covariates.csv", header=T, sep=",")

#disease/control condition
y = unlist(covariate[2], use.names=F)

#covariates w/ missing values
rawX0 = covariate[3:5]

#norm
normalize <- function(x){(x-min(x))/(max(x)-min(x))}


#impute missing values w/ mean
X0_list = lapply(rawX0, function(x) ifelse(is.na(x), mean(x, na.rm=T), x))
X0 = (matrix(unlist(X0_list), ncol = length(X0_list), byrow=F))
X0 = normalize(X0)



#covariates
X = cbind(1, X0)

#number of individuals
n = nrow(S)
#number of SNPs
m = ncol(S)
#number of covariates
k = ncol(X0)
######import dataset extracted from PGP########


# Do the computations
t0 = proc.time()[1]
# mod0 = glm( y ~ X-1, family = binomial("logit")) 
# p = mod0$fitted
# w = p * (1 - p)
# z = log(p / (1 - p)) + (y - p) / (p * (1 - p))
# xtw = t(X * w)
# U1 = xtw %*% z
# U2 = solve(xtw %*% X, U1)
# ztr = z  - X %*% U2
# U3 = xtw %*% S
# U4 = solve(xtw %*% X, U3)
# Str = S - X %*% U4
# Str2 = colSums(w * Str^2)
# b = crossprod(ztr * w, Str)/Str2
# err = sqrt(1/ Str2)
# pvalue_parallel = 2 * pnorm(-abs(b / err))
# pval = as.numeric(as.matrix(pvalue_parallel))

b1 = NULL
err1 = NULL
pval_gs = NULL
for(i in 1:m){
  mod1 = glm(y ~ S[, i] + X0, family = binomial ("logit"))  
  b1[i] = summary(mod1)$coef[2, 1]
  err1[i] = summary(mod1)$coef[2, 2]
  pval_gs[i] = summary(mod1)$coef[2, 4]
}

# Report time
t1 = proc.time()[1] - t0
msip = 1e-06 * n * m / t1
cat(sprintf("Speed: %2.1f Msips\n", msip))

###check results

t2 = proc.time()[1]
# b1 = NULL
# err1 = NULL
# pval1 = NULL
# for(i in 1:m){
#   mod1 = glm(y ~ S[, i] + X0, family = binomial ("logit"))  
#   b1[i] = summary(mod1)$coef[2, 1]
#   err1[i] = summary(mod1)$coef[2, 2]
#   pval1[i] = summary(mod1)$coef[2, 4]
# }

result=read.table("snu2_pval.csv", header = F)
# pval_he = matrix(unlist(result[2],use.names = F), ncol=1)

# result = read.csv("result-DB.txt", header = F)
Zscore = matrix(unlist(result,use.names = F), ncol=1, byrow=F)
pval_he = 2*pnorm(-abs(Zscore))

t3 = proc.time()[1] - t2
msip1 = 1e-06 * n * m / t3
cat(sprintf("Speed: %2.1f Msips\n", msip1))


threshold = 0.00001


#true positive
print("true positive rate")
tp = length(which(pval_he[which(pval_gs<=threshold)]<=threshold))/length(which(pval_gs<=threshold))


#type i error/false positive
print("type i error/false positive rate")
fp = length(which(pval_he[which(pval_gs>threshold)]<=threshold))/(length(which(pval_he[which(pval_gs>threshold)]<=threshold)) + length(which(pval_he[which(pval_gs>threshold)]>threshold)))


#type ii error/false negative
print("type ii error/false negative rate")
fn = length(which(pval_he[which(pval_gs<=threshold)]>threshold))/length(which(pval_gs<=threshold))

#true negative
print("true negative rate")
tn = length(which(pval_he[which(pval_gs>threshold)]>threshold))/length(which(pval_gs>threshold))

precision = length(which(pval_he[which(pval_gs<=threshold)]<=threshold))/length(which(pval_he<=threshold))
f1score = 2 * tp *  precision / (precision+tp)

qp<-qplot(pval_he, pval_gs, main=paste("F1 score",toString(round(f1score,digits = 3)),sep = "="), log="xy", xlab="Secure Semi-Parallel Logistic Regression", ylab="Gold standard")
qp+geom_abline(intercept = 0, slope = 1, color = "blue", size = 2)+geom_vline(xintercept = threshold, color = "red")+geom_hline(yintercept = threshold, color="red")
