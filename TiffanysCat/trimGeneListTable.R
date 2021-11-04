
library(HandyPack)
library(stringr)

rm(list=ls())

source('../CradleWare.R')
## ####################################################
## ####################################################

geneTables = c('c2.all_adult.txt',
               'c2.all_combined.txt',
               'c5.go.bp_adult.txt',
               'c5.go.bp_combined.txt',
               'hallmark_adult.txt',
               'hallmark_combined.txt')

for(f in geneTables)
{
    fileIn = paste0('tables/',f)
    fileOut = str_replace(fileIn,'\\.txt','_trimmed.txt')

    trimGeneListTable(fileIn,fileOut)
}
