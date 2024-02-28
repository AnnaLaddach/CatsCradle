
## ####################################################
#' S
#'
#' A Seurat object of 2000 genes by 1445 cells.
#'
#' @format A Seurat object
#'
#' \describe{
#' A Seurat object of cells.  It includes a UMAP of the
#' cells and annotated clustering into cell types.
#' }  
#' @source This is subset from the data associated with
#'     https://www.nature.com/articles/s41586-021-04006-z
"S"

## ####################################################
#' STranspose
#'
#' A Seurat object of 1445 cells by 2000 genes
#'
#' @format A Seurat object
#' \describe{
#' This is the transpose of S. It includes a UMAP of the
#' genes and and their clustering into gene types.
#' }
#' @source Produced from S by transposeSeuratObject()
"STranspose"

## ####################################################
#' averageExpMatrix
#'
#' The average expression matrix for cell clusters of S
#' and gene clusters of STranspose.
#'
#' @format A 30 x 13 matrix
#' \describe{
#' This gives the average gene expression for each of the
#' gene clusters in each of the cell clusters.
#' }
#' @source Produced from S and STranspose by
#' getAverageExpression()
"averageExpMatrix"

## ####################################################
#' hallmark
#'
#' The mouse hallmark gene sets.
#'
#' @format A named list
#' \describe{
#' The mouse hallmark gene sets.
#' }
#' @source https://www.gsea-msigdb.org/gsea/msigdb/mouse/collections.jsp
"hallmark"

## ####################################################
#' shorterHallmark
#'
#' An abbreviation of mouse hallmark gene sets.
#'
#' @format A named list
#' \describe{
#' The alphabetically first ten gene hallmark gene sets.
#' }
#' @source Abbreviated from hallmark
"shorterHallmark"

## ####################################################
#' clusterDF
#'
#' A data frame of genes and their gene cluster membership
#'
#' @format A data frame with columns gene and geneCluster
#' \describe{
#' This gives the cluster membership of the genes of STranspose.
#' }
#' @source Extracted from STranspose
"clusterDF"

## ####################################################
#' NN
#'
#' The nearest neighbor graph of STranspose
#'
#' @format A data frame with columns nodeA, nodeB and weight
#' \describe{
#' This gives the weigthed edges of the nearest neighbor graph
#' of the genes in STranspose.
#' }
#' @source Extracted from STranspose by getNearestNeighborListsSeurat
"NN"

## ####################################################
#' centroids
#'
#' Cell centroids from xenium spatial data
#'
#' @format a dataframe
#' where rownames are cellnames and columns contain x 
#' and y coordinates respectively.
#' \describe{
#' This gives the x and y coordinates for cell centroids from xenium  
#' mouse brain spatial data.
#' }
#' @source tiny subset from
#' https://www.10xgenomics.com/resources/datasets/fresh-frozen-mouse-brain-for-xenium-explorer-demo-1-standard
"centroids"

## ####################################################
#' clusters
#'
#' Clusters from xenium spatial data
#'
#' @format a named vector of cluster where names are each cell and
#' clusters are a factor.
#' \describe{
#' This contains cluster annotations for xenium mouse brain data extracted from
#' a Seurat analysis. 
#' }
#' @source clusters from Seurat analysis of 
#' https://www.10xgenomics.com/resources/datasets/fresh-frozen-mouse-brain-for-xenium-explorer-demo-1-standard
"clusters"

## ####################################################
#' humanLRN
#'
#' A data frame giving 12019 human ligand receptor pairs
#'
#' @format a data frame with two columns, 'from' and 'to'
#' \describe{A data frame with two columns, 'from' and 'to'.
#' Each row represents a human ligand - receptor pair.
#' }
#' @source This is derived from the nichenetr human ligand -
#' receptor network.
"humanLRN"

## ####################################################
#' mouseLRN
#'
#' A data frame giving 11592 mouse ligand receptor pairs
#'
#' @format a data frame with two columns, 'from' and 'to'
#' \describe{A data frame with two columns, 'from' and 'to'.
#' Each row represents a mouse ligand - receptor pair.
#' }
#' @source This is derived from the nichenetr mouse ligand -
#' receptor network.
"mouseLRN"

## ####################################################
#' samllXenium
#'
#' A spatial Seurat object of 4261 cells and 248 genes
#'
#' @format A Seurat object
#'
#' \describe{
#' A spatial Seurat object subset from the Xenium object used in
#' https://satijalab.org/seurat/articles/seurat5_spatial_vignette_2
#' Other data objects whose names begin with small are derived from
#' this object using CatsCradle.
#' }
#' 
#' @source This is subset from the Xenium spatial Seurat object
#' https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs.zip
#' to include a small region of the field of view surrounding the dentate gyrus. 
"smallXenium"

## ####################################################
#' smallCentroids
#'
#' A data.frame of the centroids of the cells in the
#' smallXenium Seurat object
#'
#' @format A data frame with columns x, y, and cell.  Rownames
#' also give the names of the cells.
#'
#' \describe{
#' This data frame gives the centroids of the cells in smallXenium.
#' }
#'
#' @source Extracted from smallXenium with GetTissueCoordinates().
#' The names of the cells are appended as rownames.
"smallCentroids"

## ####################################################
#' smallDelaunayTriangulation
#'
#' A data frame of the Delaunay triangulation of smallCentroids,
#' the centroids of sammXenium
#'
#' @format A data frame with character columns nodeA and nodeB.
#'
#' \describe{
#' The entries in this data frame are cell names from smallXenium.
#' Each row represents an undirected edge in the Delaunay
#' triangulation of smallCentroids.
#' }
#'
#' @source This is computed from smallCentroids using
#' computeNeighboursDelaunay()
"smallDelaunayTriangulation"

## ####################################################
#' smallNbhds
#'
#' A named list of characters in which each name is the
#' name of a cell from smallXenium and its vector gives
#' the names of the surrounding cells.
#' 
#' @format A named list.  Each entry in this lists gives
#' the names of the cells surrounding a given cell.
#'
#' \describe{
#' This list gives the combinatorial balls of radius 1
#' around the cells in smallDelaunayTriangulation
#'}
#' 
#' @source Computed from smallDelaunayTriangulation using
#' computeNeighbourhoods() with addSelf=TRUE
"smallNbhds"

## ####################################################
#' smallCombNbhds
#'
#' A named list giving the combinatorial neighbourhood
#' of radius 2 around each cell
#'
#' @format A named list.  Each entry in this lists gives
#' the names of the cells in the combinatorial ball of
#' radius 2 surrounding a given cell.
#'
#' \describe{
#' A named list.  This gives the cells in the combinatorial
#' ball of radius 2 around each cell.
#' }
#'
#' @source This is computed from smallDelaunayTriangulation using
#' findCombinatorialNeighbourhoods() with radius 2
"smallCombNbhds"

## ####################################################
#' smallNbhdMatrix
#'
#' A neighbourhoods by cell types matrix giving the count
#' of cells of each type in each neighbourhood
#'
#' @format A matrix whose rows correspond to neighbourhoods
#' (and therefore, to cells) and whose columns correspond to
#' cell types.
#'
#' \describe{
#' This matrix gives the counts for each of the cell types
#' in each neighbourhoods.  Since it is taken from the subset
#' smallXenium rather from the full Xenium object, not all
#' cell types appear.
#' }
#'
#' @source This is computed from smallNbhds and
#' smallXenium$seurat_clusters using computeNeighbourhoodByCTMatrix().
"smallNbhdMatrix"

## ####################################################
#' smallNbhdObj
#'
#' A Seurat object computed from smallNbhdMatrix. Think of
#' neighbourhoods "expressing" cell types.
#'
#' @format A Seurat object consisting of 4261 samples (the
#' neighbourhoods and 24 features (the cell types).
#'
#' \describe{
#' This is a Seurat object created by taking smallNbhdMatrix
#' as the counts.
#' }
#'
#' @source Created from smallNbhdMatrix by
#' computeNeighbourhoodByCTSeurat()
"smallNbhdObj"

## ####################################################
#' smallCellTypesPerCellTypeMatrix
#'
#' For each cell type, this matrix shows the fraction
#' of the neighbourhoods of that cell type composed of
#' each cell type.
#'
#' @format A matrix whose rows and columns correspond to
#' cell types.
#'
#' \describe{
#' Each row of this matrix corresponds to a cell type.  On
#' that row we see the proportions of all neighbourhoods
#' surrounding cells of that cell type as regards the cell types
#' they contain.  In particular, each row sums to 1.  Since
#' this comes from a subset of the full Xenium object, not all
#' cell types appear.
#' }
#'
#' @source This is created from smallNbhdMatrix and the seurat_clusters
#' of smallXenium using cellTypesPerCellTypeMatrix()
"smallCellTypesPerCellTypeMatrix"









