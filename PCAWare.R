
library(ggplot2)
library(plotly)
library(cowplot)
library(stringr)

## ####################################################
plotPCs = function(df,dataCols,coloring,labeling,outDir)
{
    dim = length(dataCols)

    nameAndMakeDir(outDir)

    plotList = list()
    finger = 1

    for(i in 1:dim)
    {
        for(j in 2:dim)
        {
            if('gene' %in% names(df))
                g = ggplot(df,aes_string(x=dataCols[j],
                                         y=dataCols[i],
                                         color=coloring,
                                         gene='gene')) +
                    geom_point()
            else
                g = ggplot(df,aes_string(x=dataCols[j],
                                         y=dataCols[i],
                                         color=coloring)) +
                    geom_point()
                

            if(i < j)
            {
                p = ggplotly(g,tooltip=labeling)
                fileName = paste0(outDir,
                               '/',
                               dataCols[i],
                               '_vs_',
                               dataCols[j],
                               '.html')
                
                htmlwidgets::saveWidget(as_widget(p),fileName)
            }

            if(max(i,j) > 5)
                next

            if(j <= i)
                plotList[[finger]] = NULL
            else
                plotList[[finger]] = g
            finger = finger + 1
        }
    }

    if(! dim %in% c(3,4,5))
    {
        writeLines('Currently supporting plots for 3, 4 or 5 dimensional data')
        return(NULL)
    }
    
    
    if(dim == 3)
        h = plot_grid(plotList[[1]],plotList[[2]],plotList[[3]],plotList[[4]],
                      ncol=dim-1)
    
    if(dim == 4)
        h = plot_grid(plotList[[1]],plotList[[2]],plotList[[3]],plotList[[4]],
                      plotList[[5]],plotList[[6]],plotList[[7]],plotList[[8]],
                      plotList[[9]],
                      ncol=dim-1)

    if(dim == 5)
        h = plot_grid(plotList[[1]],plotList[[2]],plotList[[3]],plotList[[4]],
                      plotList[[5]],plotList[[6]],plotList[[7]],plotList[[8]],
                      plotList[[9]],plotList[[10]],plotList[[11]],plotList[[12]],
                      plotList[[13]],plotList[[14]],plotList[[15]],plotList[[16]],                      
                      ncol=dim-1)

    fileName = paste0(outDir,
                      '/dimensions1-',
                      length(dataCols),
                      '.jpg')
    ggsave(plot=h,
           filename=fileName,
           width=15,height=15,units='in')
  
                
    return(h)
}

## ####################################################
makePCDF = function(f,dim)
{
    pcs = paste0('PC_',1:dim)
    df = FetchData(f,pcs)
    df$seurat_clusters = f$seurat_clusters

    return(df)
}

## ####################################################
makeGenePCDF = function(f,dim)
{
    df = makePCDF(f,dim)
    df$gene = colnames(f)

    return(df)
}

## ####################################################
makeCentroidDF = function(df)
{
    ## Get the PCs
    idx = str_detect(names(df),'PC_')
    M = data.matrix(df[,idx])
    pcs = names(df)[idx]
    clusters = unique(df$seurat_clusters)

    for(i in 1:length(clusters))
    {
        cluster = clusters[i]
        idx = df$seurat_clusters == cluster
        thisM = M[idx,]

        a = colMeans(thisM)
        A = matrix(a,nrow=1)
        colnames(A) = pcs

        b = data.frame(A)
        b$seurat_clusters = cluster

        if(i == 1)
            centroidDF = b
        else
            centroidDF = rbind(centroidDF,b)
    }
    return(centroidDF)
}

