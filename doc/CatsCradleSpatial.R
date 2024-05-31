## ----setup, include = FALSE, warning = FALSE----------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    fig.dim = c(6,6),
    comment = "#>"
)

## -----------------------------------------------------------------------------
library(Seurat)
library(CatsCradle)
ImageDimPlot(smallXenium, cols = "polychrome", size = 1)

## -----------------------------------------------------------------------------
centroids = GetTissueCoordinates(smallXenium)
rownames(centroids) = centroids$cell
delaunayNeighbours = computeNeighboursDelaunay(centroids)
head(delaunayNeighbours)

## -----------------------------------------------------------------------------
idx = (delaunayNeighbours$nodeA == 16307 |
       delaunayNeighbours$nodeB == 16307)
nbhd = unique(c(delaunayNeighbours$nodeA[idx],
                delaunayNeighbours$nodeB[idx]))
nbhd		

## -----------------------------------------------------------------------------
extendedNeighboursList = getExtendedNBHDs(delaunayNeighbours, 4)
extendedNeighbours = collapseExtendedNBHDs(extendedNeighboursList, 4)

## -----------------------------------------------------------------------------
idx = (extendedNeighbours$nodeA == 16307 |
       extendedNeighbours$nodeB == 16307)
nbhd = unique(c(extendedNeighbours$nodeA[idx],
                extendedNeighbours$nodeB[idx]))
length(nbhd)		

## -----------------------------------------------------------------------------
euclideanNeighbours = computeNeighboursEuclidean(centroids,
threshold=20)

## -----------------------------------------------------------------------------
agg = aggregateSeuratGeneExpression(smallXenium,extendedNeighbours,
                                    verbose=FALSE)
smallXenium$aggregateNBHDClusters = agg@active.ident
ImageDimPlot(smallXenium,group.by='aggregateNBHDClusters',cols='polychrome')

## -----------------------------------------------------------------------------
NBHDByCTMatrix = computeNBHDByCTMatrix(delaunayNeighbours, clusters)

## -----------------------------------------------------------------------------
NBHDByCTMatrixExtended = 
  computeNBHDByCTMatrix(extendedNeighbours, clusters)

## -----------------------------------------------------------------------------
cellTypesPerCellTypeMatrix = 
  computeCellTypesPerCellTypeMatrix(NBHDByCTMatrix,clusters)

## -----------------------------------------------------------------------------
colours = DiscretePalette(length(levels(clusters)), palette = "polychrome")
names(colours) = levels(clusters)

cellTypesPerCellTypeGraphFromCellMatrix(cellTypesPerCellTypeMatrix, 
                                    minWeight = 0.05, colours = colours)

## -----------------------------------------------------------------------------
library(pheatmap)
pheatmap(cellTypesPerCellTypeMatrix)

## -----------------------------------------------------------------------------
cellTypesPerCellTypeMatrixExtended = computeCellTypesPerCellTypeMatrix(NBHDByCTMatrixExtended,clusters)

cellTypesPerCellTypeGraphFromCellMatrix(cellTypesPerCellTypeMatrixExtended,
minWeight = 0.05, colours = colours)

## -----------------------------------------------------------------------------
cellTypesPerCellTypePValues = computeNeighbourEnrichment(delaunayNeighbours, 
                                                     clusters, verbose = F)

## -----------------------------------------------------------------------------
cellTypesPerCellTypePValuesNegLog = -log10(cellTypesPerCellTypePValues)
pheatmap(cellTypesPerCellTypePValuesNegLog)

## -----------------------------------------------------------------------------
NBHDByCTSeurat = computeNBHDVsCTSeurat(NBHDByCTMatrix,verbose=FALSE)

## -----------------------------------------------------------------------------
NBHDByCTSeurat$cellType = clusters

## -----------------------------------------------------------------------------
DimPlot(NBHDByCTSeurat, group.by = c("cellType"), cols = "polychrome", reduction = "umap")
DimPlot(NBHDByCTSeurat, group.by = c("neighbourhood_clusters"), cols = "polychrome", reduction = "umap")

## -----------------------------------------------------------------------------
smallXenium$NBHDCluster = NBHDByCTSeurat@active.ident
ImageDimPlot(smallXenium, group.by = "NBHDCluster", size = 1, cols = "polychrome")

## -----------------------------------------------------------------------------
NBHDByCTSeuratExtended = computeNBHDVsCTSeurat(NBHDByCTMatrixExtended,
                                               verbose=FALSE)

## -----------------------------------------------------------------------------
NBHDByCTSeuratExtended$cellType = clusters

## -----------------------------------------------------------------------------
DimPlot(NBHDByCTSeuratExtended, group.by = c("cellType"), cols = "polychrome", reduction = "umap")
DimPlot(NBHDByCTSeuratExtended, group.by = c("neighbourhood_clusters"), cols = "polychrome", reduction = "umap")

## -----------------------------------------------------------------------------
smallXenium$NBHDClusterExtended= 
  NBHDByCTSeuratExtended@active.ident
ImageDimPlot(smallXenium, group.by = c("NBHDClusterExtended"), 
             size = 1, cols = "polychrome")

## -----------------------------------------------------------------------------
CTByNBHDCluster = table(NBHDByCTSeurat$cellType,NBHDByCTSeurat@active.ident)
CTByNBHDCluster = CTByNBHDCluster/rowSums(CTByNBHDCluster)

rownames(CTByNBHDCluster) = paste0("CellType",rownames(CTByNBHDCluster))
colnames(CTByNBHDCluster) = paste0("NBHDCluster",colnames(CTByNBHDCluster))

pheatmap(CTByNBHDCluster,
      fontsize_row=8,
      fontsize_col=8,
      cellheight=10,
      cellwidth=10)

sankeyFromMatrix(CTByNBHDCluster)

## -----------------------------------------------------------------------------
CTByNBHDClusterExtended = table(NBHDByCTSeuratExtended$cellType,NBHDByCTSeuratExtended@active.ident)
CTByNBHDClusterExtended = CTByNBHDClusterExtended/rowSums(CTByNBHDClusterExtended)

rownames(CTByNBHDClusterExtended) = paste0("CellType",rownames(CTByNBHDClusterExtended))
colnames(CTByNBHDClusterExtended) = paste0("NBHDCluster",colnames(CTByNBHDClusterExtended))

pheatmap(CTByNBHDClusterExtended,
      fontsize_row=8,
      fontsize_col=8,
      cellheight=10,
      cellwidth=10)

sankeyFromMatrix(CTByNBHDClusterExtended)

## -----------------------------------------------------------------------------
CTByNBHDSeurat = 
  computeNBHDVsCTSeurat(t(NBHDByCTMatrix), npcs = 10, 
                        transpose = T, resolution = 1, n.neighbors = 5,
			verbose=FALSE)

CTByNBHDSeurat$cellType = colnames(CTByNBHDSeurat)

DimPlot(CTByNBHDSeurat, group.by = "cellType", cols = "polychrome", 
        reduction = "umap", label = T)

## -----------------------------------------------------------------------------
CTByNBHDSeurat= computeGraphEmbedding(CTByNBHDSeurat)

DimPlot(CTByNBHDSeurat,group.by = "cellType", cols = "alphabet", reduction = "graph", label = T)

## -----------------------------------------------------------------------------
ImageDimPlot(smallXenium, cols = "polychrome", size = 1)

## -----------------------------------------------------------------------------
pca = Embeddings(CTByNBHDSeurat, reduction = "pca")
res = pheatmap(pca)

## -----------------------------------------------------------------------------
CTClust = cutree(res$tree_row, k = 11)
CTByNBHDSeurat$neighbourhood_clusters = factor(CTClust)

## -----------------------------------------------------------------------------
CTComposition = table(CTByNBHDSeurat$cellType, CTByNBHDSeurat$neighbourhood_clusters)
pheatmap(CTComposition)

## -----------------------------------------------------------------------------
averageExpMatrix = getAverageExpressionMatrix(NBHDByCTSeurat,
                           CTByNBHDSeurat,
			    clusteringName='neighbourhood_clusters')
averageExpMatrix = tagRowAndColNames(averageExpMatrix,
                                     ccTag='neighbourhoodClusters_',
                                     gcTag='cellTypeClusters_')
pheatmap(averageExpMatrix,
      fontsize_row=8,
      fontsize_col=8,
      cellheight=10,
      cellwidth=10)


sankeyFromMatrix(averageExpMatrix)

## -----------------------------------------------------------------------------
moransI = runMoransI(smallXenium, delaunayNeighbours, assay = "SCT", 
                     layer = "data", nSim = 100, verbose = FALSE)

## -----------------------------------------------------------------------------
head(moransI)

## -----------------------------------------------------------------------------
tail(moransI)

## -----------------------------------------------------------------------------
ImageFeaturePlot(smallXenium, "Nwd2")

## -----------------------------------------------------------------------------
ImageFeaturePlot(smallXenium, "Trbc2")

## -----------------------------------------------------------------------------
ligandReceptorResults = performLigandReceptorAnalysis(smallXenium, delaunayNeighbours, 
                                                "mouse", clusters,verbose=FALSE)

## -----------------------------------------------------------------------------
head(ligandReceptorResults$interactionsOnEdges[,1:10],10)

## -----------------------------------------------------------------------------
head(ligandReceptorResults$totalInteractionsByCluster[,1:10],10)

## -----------------------------------------------------------------------------
head(ligandReceptorResults$meanInteractionsByCluster[,1:10],10)

## -----------------------------------------------------------------------------
head(ligandReceptorResults$simResults[,1:10],10)

## -----------------------------------------------------------------------------
head(ligandReceptorResults$pValues,10)

## -----------------------------------------------------------------------------
ligRecMatrix = makeLRInteractionHeatmap(ligandReceptorResults, clusters, colours = colours, labelClusterPairs = F)

## -----------------------------------------------------------------------------
cellTypePerCellTypeLigRecMatrix = makeSummedLRInteractionHeatmap(ligandReceptorResults, clusters, "total")

## -----------------------------------------------------------------------------
cellTypePerCellTypeLigRecMatrix = makeSummedLRInteractionHeatmap(ligandReceptorResults, clusters, "mean", logScale = T)

## -----------------------------------------------------------------------------
hist(cellTypePerCellTypeLigRecMatrix)

## -----------------------------------------------------------------------------
cellTypesPerCellTypeGraphFromCellMatrix(cellTypePerCellTypeLigRecMatrix, 
                                    minWeight = 0.4, colours = colours)

## -----------------------------------------------------------------------------
scaleFactor = 3
cellTypesPerCellTypeGraphFromCellMatrix(cellTypePerCellTypeLigRecMatrix/scaleFactor, 
                                    minWeight = 0.4/scaleFactor, colours = colours)

## -----------------------------------------------------------------------------
edgeSeurat = computeEdgeSeurat(ligandReceptorResults, centroids)

## -----------------------------------------------------------------------------
ImageFeaturePlot(edgeSeurat, features = "Pdyn-Npy2r")

## -----------------------------------------------------------------------------
edgeNeighbours = computeEdgeGraph(delaunayNeighbours)

## -----------------------------------------------------------------------------
moransILigandReceptor = runMoransI(edgeSeurat, edgeNeighbours, assay = "RNA", 
                     layer = "counts", nSim = 100)

## -----------------------------------------------------------------------------
head(moransILigRec)

## -----------------------------------------------------------------------------
tail(moransILigRec)

## -----------------------------------------------------------------------------
ImageFeaturePlot(edgeSeurat, "Penk-Htr1f")

## -----------------------------------------------------------------------------
ImageFeaturePlot(edgeSeurat, "Sst-Gpr17")

