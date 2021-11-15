
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
    df = df[order(as.numeric(df$seurat_clusters)),]

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
    centroidDF = centroidDF[order(as.numeric(centroidDF$seurat_clusters)),]
    return(centroidDF)
}

## ####################################################
makePCSeparationDFEachVsEach = function(f,dim)
{
    pcDF = makePCDF(f,dim)
    centroidDF = makeCentroidDF(pcDF)
    
    clusters = unique(centroidDF$seurat_clusters)
    clusters = as.character(clusters)
    clusters = as.numeric(clusters)
    clusters = clusters[order(clusters)]

    I = c()
    J = c()
    PC = c()
    separation = c()
    finger = 1
    for(i in 1:(length(clusters)-1))
    {
        for(j in (i+1):length(clusters))
        {
            for(pc in 1:dim)
            {
                I[finger] = clusters[i]
                J[finger] = clusters[j]
                PC[finger] = pc
                separation[finger] = abs(centroidDF[i,pc] - centroidDF[j,pc])
                finger = finger + 1
            }
        }
    }

    separationDF = data.frame(cluster1=I,cluster2=J,
                              PC,separation)

    separationDF$comparison = paste(separationDF$cluster1,
                                    'vs',
                                    separationDF$cluster2)
    comparisonLevels = c()
    for(i in clusters)
        comparisonLevels = c(comparisonLevels,
                             paste(i,'vs',clusters))
    idx = comparisonLevels %in% separationDF$comparison
    comparisonLevels = comparisonLevels[idx]
    separationDF$comparison = factor(separationDF$comparison,
                                     levels=comparisonLevels)
    
    return(separationDF)
}

## ####################################################
makePCSeparationPlotEachVsEach = function(f,dim)
{
    df = makePCSeparationDFEachVsEach(f,dim)


    idx = df$cluster1 == 5 |
        df$cluster2 == 5
    

    g = ggplot(df[idx,],aes(x=PC,y=separation)) +
        geom_line() +
        geom_point() +
        facet_wrap(~comparison,ncol=2)

    return(g)
}
    
## ####################################################
makePCSeparationPlotEachVsEach = function(f,dim,outDir=NULL)
{
    df = makePCSeparationDFEachVsEach(f,dim)

    clusters = unique(c(df$cluster1,df$cluster2))
    plotList = list()

    for(cluster in clusters)
    {
        idx = (df$cluster1 == cluster |
               df$cluster2 == cluster)
    
    xTicks = seq(from=min(clusters),to=max(clusters),b=2)
    g = ggplot(df[idx,],aes(x=PC,y=separation)) +
        geom_line() +
        geom_point() +
        facet_wrap(~comparison,ncol=2) +
        scale_x_continuous(breaks=xTicks)

        plotList[[paste0('cluster_',cluster)]] = g

        if(! is.null(outDir))
        {
            fileName = paste0(outDir,'/cluster_',cluster,'.jpg')
            ggsave(plot=g,
                   filename=fileName)
        }
    }
    return(plotList)
}
    
