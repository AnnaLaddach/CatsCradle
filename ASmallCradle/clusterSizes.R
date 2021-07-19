
library(HandyPack)
library(ggplot2)
library(cowplot)

rm(list=ls())
graphics.off()

## ###################################################
## ###################################################

res = 1
assays = c('integrated','RNA')
figs = list()

for(assay in assays)
{
    fileName = paste0(assay,'_resolution_',res,'/geneClusters.txt')
    clusterDF = Read.Table(fileName)

    clusters = unique(clusterDF$geneCluster)
    clusters = clusters[order(clusters)]

    count = c()
    for(cluster in clusters)
        count[as.character(cluster)] = sum(clusterDF$geneCluster == cluster)

    countDF = data.frame(clusters,count)

    figs[[assay]] = ggplot(countDF,aes(x=clusters,y=count)) +
        geom_col() +
        scale_x_continuous(breaks=clusters) +
        ggtitle(paste('Cluster sizes',assay,'resolution',res))
}

g = plot_grid(figs[[1]],figs[[2]])
print(g)

ggsave(plot=g,
       filename=paste0('figures/clusterSizes_resolution_',res,'.jpg'))
