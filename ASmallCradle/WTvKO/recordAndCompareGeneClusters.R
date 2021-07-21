
library(Seurat)
library(HandyPack)
library(BioHandy)

rm(list=ls())
graphics.off()

source('../CradleWare.R')
source('QuadWare.R')

## ###################################################
## ###################################################

res = 1
assays = c('integrated','RNA')

for(assay in assays)
{
    objectQuad = getObjectQuad(assay,res)
    fWT = objectQuad$fWT
    fKO = objectQuad$fKO
    fPrimeWT = objectQuad$fPrimeWT
    fPrimeKO = objectQuad$fPrimeKO

    WTGeneDF = data.frame(gene=rownames(fPrimeWT@meta.data),
                          cluster=fPrimeWT@meta.data$seurat_clusters)
    KOGeneDF = data.frame(gene=rownames(fPrimeKO@meta.data),
                          cluster=fPrimeKO@meta.data$seurat_clusters)    

    outDir = paste0(assay,'_resolution_',res)
    if(! dir.exists(outDir))
        dir.create(outDir)
    
    fileName = paste0(outDir,'/WT_geneClusters.txt')
    Write.Table(WTGeneDF,
                fileName)

    fileName = paste0(outDir,'/KO_geneClusters.txt')
    Write.Table(KOGeneDF,
                fileName)
    
    WTclusters = unique(WTGeneDF$cluster)
    WTclusters = WTclusters[order(WTclusters)]

    KOclusters = unique(KOGeneDF$cluster)
    KOclusters = KOclusters[order(KOclusters)]

    M = matrix(0,nrow=length(WTclusters),ncol=length(KOclusters))
    rownames(M) = as.character(WTclusters)
    colnames(M) = as.character(KOclusters)

    for(i in WTclusters)
    {
        for(j in KOclusters)
        {
            a = WTGeneDF$gene[WTGeneDF$cluster == i]
            b = KOGeneDF$gene[KOGeneDF$cluster == j]
            c = intersect(a,b)

            M[as.character(i),as.character(j)] =
                geneListPValue(length(a),
                               length(b),
                               length(c),
                               nrow(WTGeneDF))
        }
    }
    p = sankeyFromMatrix(M)
    print(p)
}
