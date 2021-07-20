
library(networkD3)
library(plotly)
library(HandyPack)
library(BioHandy)
library(pheatmap)
library(stringr)

rm(list=ls())

source('CradleWare.R')
source('~/FranzeSingleCell07/SeuratWare02.R')

## ###################################################
## ###################################################

background = c(integrated=2000,
               RNA=11134)

for(res in 1)
{
    for(assay in c('integrated','RNA'))
    {
        ## Get the DE genes:
        geneLists = list()
        DEdir = paste0(assay,'_resolution_',res,'/DE/f')
        files = Sys.glob(paste0(DEdir,'/*'))
        groups = str_replace(files,'.*group_','')
        groups = str_replace(groups,'\\.txt','')
        groups = as.numeric(groups)
        groups = groups[order(groups)]
        
        for(i in 1:length(files))
        {
            geneLists[[getShortName(groups[i])]] =
                Read.Table(files[i])$id
        }
        
        ## Make the overlap matrix:
        N = length(groups)
        M = matrix(0,nrow=N,ncol=N)

        rownames(M) = getShortNames()
        colnames(M) = getShortNames()
        
        for(i in groups)
        {
            for(j in groups)
            {
                a = geneLists[[getShortName(i)]]
                b = geneLists[[getShortName(j)]]
                c = intersect(a,b)
                
                A = length(a)
                B = length(b)
                C = length(c)
                M[getShortName(i),getShortName(j)] =
                    -log10(geneListPValue(A,B,C,background[assay]))
            }
        }

        ## A hack:
        M[! is.finite(M)] = max(M[is.finite(M)]) + 1

        ## Get rid of diagonal for Sankey graph:
        MM = M
        for(i in 1:N)
            MM[i,i] = 0
        
        p = sankeyFromMatrix(MM)
        print(p)

        outDir = 'DEGeneOverlap'
        if(! dir.exists(outDir))
            dir.create(outDir)

        fileName = paste0(outDir,'/',assay,'_DEGeneOverlap.html')
        htmlwidgets::saveWidget(as_widget(p), fileName)

        fileName = str_replace(fileName,'\\.html','.jpg')
        jpeg(fileName,
             height=8,width=8,units='in',res=100)
        ## Take log for heatmap:
        pheatmap(log10(M+1),
                 treeheight_row=0,
                 treeheight_col=0)
        dev.off()
        
    } ## assay
} ## res
