## ----setup, include = FALSE, warning = FALSE----------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

## ----message = FALSE, warning=FALSE, eval=FALSE-------------------------------
#  library(CatsCradle)
#  STranspose = transposeSeuratObject(S)

## ----eval = FALSE-------------------------------------------------------------
#  library(plotly)
#  umap = FetchData(STranspose,c('UMAP_1','UMAP_2','seurat_clusters'))
#  umap$gene = colnames(STranspose)
#  plot = ggplot(umap,aes(x=UMAP_1,y=UMAP_2,color=seurat_clusters,label=gene) +
#         geom_point()
#  browseable = ggplotly(plot)
#  print(browseable)
#  htmlwidgets::saveWidget(as_widget(browseable),'genesOnUMAP.html')

## ----eval = FALSE-------------------------------------------------------------
#  library(pheatmap)
#  averageExpMatrix = getAverageExpressionMatrix(S,STranspose)
#  averageExpMatrix = tagRowAndColNames(averageExpMatrix,
#                                       ccTag='cellClusters_',
#                                       gcTag='geneClusters_')
#  pheatmap(averageExpMatrix,
#        treeheight_row=0,
#        treeheight_col=0,
#        fontsize_row=8,
#        fontsize_col=8,
#        cellheight=10,
#        cellwidth=10)

## ----eval = FALSE-------------------------------------------------------------
#  catsCradle =  = sankeyFromMatrix(averageExpMatrix,
#                                disambiguation=c('cells_','genes_'),
#                                plus='cyan',minus='pink',
#                                height=800)
#  print(catsCradle)

## ----eval = FALSE-------------------------------------------------------------
#  g2mGenes = intersect(colnames(STranspose),
#                       hallmark[['HALLMARK_G2M_CHECKPOINT']])
#  stats = getSeuratSubsetClusteringStatistics(STranspose,
#                                        g2mGenes,
#                                        numTrials=1000)

## ----eval = FALSE-------------------------------------------------------------
#  geneSet = intersect(colnames(STranspose),
#                      hallmark[['HALLMARK_INTERFERON_ALPHA_RESPONSE']])
#  geometricallyNearbyGenes = getNearbyGenes(STranspose,geneSet,radius=0.2,metric='umap')
#  theGeometricGenesThemselves = names(geometricallyNearbyGenes)
#  combinatoriallyNearbyGenes = getNearbyGenes(STranspose,geneSet,radius=1,metric='NN')
#  theCombinatoricGenesThemselves = names(combinatoriallyNearbyGenes)
#  

## -----------------------------------------------------------------------------
library(CatsCradle)
annotatedGenes = annotateGenesByGeneSet(hallmark)
names(annotatedGenes[['Myc']])

## -----------------------------------------------------------------------------
 Myc = annotateGeneAsVector('Myc',hallmark)
 MycNormalised = annotateGeneAsVector('Myc',hallmark,TRUE)

## -----------------------------------------------------------------------------
predicted = predictAnnotation('Myc',hallmark,STranspose,radius=.5)
predicted$Myc

