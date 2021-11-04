

library(HandyPack)
library(stringr)

rm(list=ls())

source('../CradleWare.R')
## ####################################################
## ####################################################
assays = c('RNA','integrated')
resolutions = 1:2
files = c('c2.all_integrated.txt',
          'c5.go.bp_integrated.txt',
          'hallmark_integrated.txt')

for(assay in assays)
{
    for(res in resolutions)
    {
        for(f in files)
        {
            fileIn = paste0(assay,
                            '_resolution_',
                            res,
                            '/',
                            f)
            fileOut = str_replace(fileIn,'\\.txt','_trimmed.txt')

            if(! file.exists(fileIn))
            {
                writeLines(paste('skipping',fileIn))
                next
            }
            
            trimGeneListTable(fileIn,fileOut)
        }
    }
}
