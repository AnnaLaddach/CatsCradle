
source('CradleWare.R')

M = matrix(c(1,5,2,7,3,6),nrow=2)

rownames(M) = c('A','B')
colnames(M) = c('X','Y','Z')

p = sankeyFromMatrix(M)

print(p)
