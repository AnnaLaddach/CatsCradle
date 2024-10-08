
## ####################################################
#' exSeuratObj
#'
#' A Seurat object of 2000 genes by 540 cells.
#'
#' @format A Seurat object
#'
#' \describe{
#' A Seurat object of cells.  It includes a UMAP of the
#' cells and annotated clustering into cell types. It has
#' been severely reduced in size to accommodate Bioconductor
#' size restrictions.
#' }  
#' @source This is subset from the data associated with
#'     https://www.nature.com/articles/s41586-021-04006-z
"exSeuratObj"

## ####################################################
#' humanLRN
#'
#' A data frame giving 12019 human ligand receptor pairs
#'
#' @format a data frame with two columns, 'from' and 'to'
#' \describe{A data frame with two columns, 'from' and 'to'.
#' Each row represents a human ligand - receptor pair.
#' }
#' @source This is taken from the nichenetr package,
#' url = {https://www.nature.com/articles/s41592-019-0667-5}.
#' Specifically we use the human ligand - receptor network.
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
#' @source This is taken from the nichenetr package,
#' url = {https://www.nature.com/articles/s41592-019-0667-5}.
#' Specifically, we use the mouse ligand - receptor network.
"mouseLRN"

## ####################################################
#' smallXenium
#'
#' A spatial Seurat object of 4261 cells and 248 genes
#'
#' @format A Seurat object
#'
#' \describe{
#' A spatial Seurat object subset from the Xenium object used in
#' https://satijalab.org/seurat/articles/seurat5_spatial_vignette_2.
#' }
#' 
#' @source This is subset from the Xenium spatial Seurat object
#' https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs.zip
#' to include a small region of the field of view surrounding the dentate gyrus. 
"smallXenium"

## ####################################################
#' moransI
#'
#' A data fame containing Moran's I and related pvalues.
#'
#' @format A data fame containing Moran's I and related pvalues.
#'
#' \describe{
#' Moran's I values calculated for the genes in smallXenium (using the SCT 
#' assay). Pvalues derived using 100 permutations.
#' }
#'
#' @source Created from smallXenium and delaunayNeighbours by using
#' runMoransI()
"moransI"


## ####################################################
#' ligandReceptorResults
#'
#' The result of performLigandReceptorAnalysis(smallXenium, delaunayNeighbours, 
#' "mouse", clusters,verbose=FALSE) 
#'
#' @format A list of data frames.
#'
#' \describe{
#' A list containing:
#' interactionsOnEdges - a data frame whose first two columns give
#' the neighbouring cells and next two columns give their corresponding 
#' clusters. Each of the remaining columns is a logical
#' corresponding to a ligand-receptor pair telling whether the ligand
#' is expressed in the first cell and the receptor is expressed in the
#' second cell.
#' totalInteractionsByCluster - a dataframe where the first column gives a 
#' directed (sender-receiver) pair of clusters. The second column gives the 
#' total number of edges between those clusters. The remaining columns give the 
#' total numbers of edges on which particular ligand receptor interactions are 
#' present.
#' meanInteractionsByCluster - a dataframe where the first column gives a 
#' directed (sender-receiver) pair of clusters. The second column gives the 
#' total number of edges between those clusters. The remaining columns give the 
#' total numbers of edges on which particular ligand receptor interactions are 
#' present (for that cluster pair) divided by the total number of edges between 
#' those clusters.
#' simResults - a dataframe where the rownames are sender-receiver cluster pairs 
#' and column names are ligand receptor pairs. Values give the number of 
#' simulations for which observed values are greater than simulated values.
#' pValues - a dataframe where the rownames are sender-receiver cluster pairs 
#' and column names are ligand receptor pairs. Entries are uppertail pvalues 
#' describing whether a particular ligand receptor interaction is observed more 
#' frequently between 2 clusters than expected.
#' }
#'
#' @source Created from smallXenium and delaunayNeighbours by using
#' performLigandReceptorAnalysis(()
"ligandReceptorResults"

## ####################################################
#' moransILigandReceptor
#' 
#' Moran's I for the ligand receptor pairs
#'
#' @format A data frame showing the spatial autocorrelation of the
#' 28 ligand receptor pairs
#'
#' \describe{
#' A data frame with rownames giving the  28 ligand-receptor pairs and columns moransI
#' and pValues
#' }
#'
#' @source Computed using the function runMoransI on the object edgeSeurat and
#' neighbours edgeNeighbours = computeEdgeGraph(delaunayNeighbours) with 100 trials.
#' For more informations see the CatsCradleSpatial vignette.
"moransILigandReceptor"

## ####################################################
#' seuratGenes
#'
#' A vector of genes used for subsetting exSeuratObj
#'
#' @format A vector of genes
#'
#' \describe{
#' A vector of the top 100 most variable genes in exSeuratObj
#' used to subset this object to give toy examples.
#' }
#'
#' @source Computed by retrieving the data layer from exSeuratObj and
#' subsetting to the 100 genes with the highest standard deviation.
"seuratGenes"

## ####################################################
#' seuratCells
#'
#' A vector of cells used for subsetting exSeuratObj
#'
#' @format A vector of cells
#'
#' \describe{
#' A vector of cells consisting of half the cells from
#' each seurat_cluster in exSeuratObj used to subset this
#' object to give toy examples.
#' }
#'
#' @source Computed by retrieving half the cells from each
#' cluster in exSeuratObj
"seuratCells"

## ####################################################
#' xeniumCells
#'
#' A vector of cells used for subsetting exSeuratObj
#'
#' @format A vector of cells
#'
#' \describe{
#' A vector of cells consisting of approximately one
#' quarter of the cells in smallXenium used to subset
#' this object to give toy examples.
#' }
#'
#' @source We extracted a rectangle whose width and
#' height were one half the width and height of smallXenium
#' and which was centered in the field of view of
#' smallXenium
"xeniumCells"






