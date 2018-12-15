newHEpvals = read.table("HE_Pvals.txt", header=F)
newHEpvals = matrix(unlist(newHEpvals, use.names=F), ncol = length(newHEpvals), byrow=F)
newHEpvals = newHEpvals[,2]

library(MLmetrics)

cutoff = 0.01

label.gold = ifelse(pval1<cutoff,0,1)
label.semilog = ifelse(pval<cutoff,0,1)
label.newHElog = ifelse(newHEpvals<cutoff,0,1)

F1_Score(label.gold,label.newHElog,positive="0")

F1_Score(label.semilog,label.newHElog,positive="0")

