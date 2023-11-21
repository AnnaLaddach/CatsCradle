
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







