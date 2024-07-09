
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
#' S_sce
#'
#' A SingleCellExperiment object of 2000 genes by 1445 cells.
#'
#' @format A SingleCellExperiment object
#'
#' \describe{
#' A SingleCellExperiment version of S.
#' }  
#' @source This was made from S using as.SingleCellExperiment
"S_sce"

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
#' STranspose_sce
#'
#' A SingleCellExperiment object of 1445 cells by 2000 genes
#'
#' @format A SingleCellExperiment object
#' \describe{
#' This is a SingleCellExperiment version of STranspose.
#' }
#' @source Produced from STranspose using as.SingleCellExperiment
#' together with ancillary code for copying nearest neighbour
#' graphs
"STranspose_sce"

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
#' @source subset of tiny subset from
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
#' @source clusters from Seurat analysis of subset of
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
#' delaunayNeighbours
#'
#' A data frame of the Delaunay triangulation of centroids of cells in
#' smallXenium
#'
#' @format A data frame with character columns nodeA and nodeB.
#'
#' \describe{
#' The entries in this data frame are cell names from smallXenium.
#' Each row represents an undirected edge in the Delaunay
#' triangulation of centroids.
#' }
#'
#' @source This is computed from centroids using
#' computeNeighboursDelaunay()
"delaunayNeighbours"


## ####################################################
#' euclideanNeighbours
#'
#' A data frame of nearest neighbours calculated 
#' using euclidean distance as a cutoff.
#'
#' @format A data frame with character columns nodeA and nodeB.
#'
#' \describe{
#' The entries in this data frame are cell names from smallXenium.
#' Each row represents an undirected edge in the nearest neighbour 
#' graph using a cutoff euclidean distance of 20 \\mu m.
#' }
#'
#' @source This is computed from centroids using
#' computeNeighboursEuclidean()
"euclideanNeighbours"


## ####################################################
#' extendedNeighboursList
#'
#' A list of nth degree neighbour graphs calculated from delaunayNeighbours.
#'
#' @format A list of dataframes with character columns nodeA and nodeB.
#'
#' \describe{
#' A named list of neighbour graphs, where each graph contains edges 
#' connecting vertices of degree n. Each graph (list entry) is named according to degree n.
#' }
#'
#' @source This is computed from delaunayNeighbours using
#' getExtendedNBHDs()
"extendedNeighboursList"

## ####################################################
#' extendedNeighbours
#'
#' A neighbour graph where nodes are connected to all nodes up to 
#' and including degree 4 in the original graph (delaunayNeighbours)
#'
#' @format A data frame with character columns nodeA and nodeB.
#'
#' \describe{
#' The entries in this data frame are cell names from smallXenium.
#' Each row represents an undirected edge in the extended neighbour 
#' graph. Nodes are connected to all nodes up to 
#' and including degree 4 in the original graph (delaunayNeighbours).
#' }
#'
#' @source This is computed from extendedNeighboursList using
#' collapseExtendedNBHDs() 
"extendedNeighbours"

## ####################################################
#' edgeNeighbours
#'
#' A neighbour graph where nodes represent interactions between cells. 
#'
#' @format A data frame with character columns nodeA and nodeB.
#'
#' \describe{
#' A spatial graph where edges in the original delaunayNeighbours become nodes 
#' and A-B edges (in the original graph) become connected to
#' all A- edges and all B- edges. 
#' }
#'
#' @source This is computed from delaunayNeighbours using
#' computeEdgeGraph()
"edgeNeighbours"


## ####################################################
#' NBHDByCTMatrix
#'
#' A neighbourhoods by cell types matrix giving the count
#' of cells of each type in each neighbourhood.
#'
#' @format A matrix whose rows correspond to neighbourhoods
#' (and therefore, to cells) and whose columns correspond to
#' cell types.
#'
#' \describe{
#' This matrix gives the counts for each of the cell types
#' in each neighbourhood.
#' }
#'
#' @source This is computed from delaunayNeighbours and
#' clusters using computeNBHDByCTMatrix().
"NBHDByCTMatrix"


## ####################################################
#' NBHDByCTMatrixExtended
#'
#' A neighbourhoods by cell types matrix giving the count
#' of cells of each type in each extended neighbourhood.
#'
#' @format A matrix whose rows correspond to extended neighbourhoods
#' (and therefore, to cells) and whose columns correspond to
#' cell types.
#'
#' \describe{
#' This matrix gives the counts for each of the cell types
#' in each extended neighbourhood.
#' }
#'
#' @source This is computed from delaunayNeighbours and
#' clusters using computeNBHDByCTMatrix().
"NBHDByCTMatrixExtended"


## ####################################################
#' NBHDByCTSeurat
#'
#' A Seurat object computed from NBHDByCTMatrix. Think of
#' neighbourhoods "expressing" cell types.
#'
#' @format A Seurat object consisting of 4261 samples (the
#' neighbourhoods and 24 features (the cell types).
#'
#' \describe{
#' This is a Seurat object created by taking NBHDByCTMatrix
#' as the counts.
#' }
#'
#' @source Created from NBHDByCTMatrix by
#' computeNBHDVsCTSeurat()
"NBHDByCTSeurat"


## ####################################################
#' NBHDByCTSeuratExtended
#'
#' A Seurat object computed from NBHDByCTMatrixExtended. Think of
#' neighbourhoods "expressing" cell types.
#'
#' @format A Seurat object consisting of 4261 samples (the
#' neighbourhoods) and 24 features (the cell types).
#'
#' \describe{
#' This is a Seurat object created by taking NBHDByCTMatrixExtended
#' as the counts.
#' }
#'
#' @source Created from NBHDByCTMatrixExtended by
#' computeNBHDVsCTSeurat()
"NBHDByCTSeuratExtended"

## ####################################################
#' cellTypesPerCellTypeMatrix
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
#' they contain.  In particular, each row sums to 1.
#' }
#'
#' @source This is created from NBHDByCTMatrix and the clusters
#' using cellTypesPerCellTypeMatrix()
"cellTypesPerCellTypeMatrix"

## ####################################################
#' cellTypesPerCellTypeMatrixExtended
#'
#' For each cell type, this matrix shows the fraction
#' of the neighbourhoods of that cell type composed of
#' each cell type. This uses the extended neighbourhoods
#' of combinatorial radius 4
#'
#' @format A matrix whose rows and columns correspond to
#' cell types.
#'
#' \describe{
#' Each row of this matrix corresponds to a cell type.  On
#' that row we see the proportions of all neighbourhoods
#' surrounding cells of that cell type as regards the cell types
#' they contain.  In particular, each row sums to 1.  This uses
#' the extended neighbourhoods of combinatorial radius 4.
#' }
#'
#' @source This is created from NBHDByCTMatrixExtended and the
#' clusters using cellTypesPerCellTypeMatrix()
"cellTypesPerCellTypeMatrixExtended"


## ####################################################
#' cellTypesPerCellTypePValues
#' A symmetric matrix containing P values describing whether cell types are more 
#' frequently neighbours than expected by chance. 
#' 
#' @format A matrix whose rows and columns correspond to
#' cell types.
#'
#' \describe{
#' Rows and columns of this matrix correspond to a cell types. Matrix give p 
#' values describing whether cell types are more frequently neighbours than 
#' expected by chance.
#' }
#'
#' @source This is created from delaunayNeighbours and the clusters
#' using computeNeighbourEnrichment()
"cellTypesPerCellTypePValues"


## ####################################################
#' colours
#'
#' A character vector of colours where names are clusters.
#'
#' @format A character vector of colours where names are clusters.
#'
#' \describe{
#' A character vector of colours from the polychrome palette where names 
#' are clusters.
#' }
#'
#' @source Created using Seurat DiscretePalette() 
"colours"


## ####################################################
#' CTByNBHDSeurat
#'
#' A Seurat object computed from the transpose of NBHDByCTMatrix. Think of
#' cell types "expressing" (being found in) neighbourhoods.
#'
#' @format A Seurat object consisting of 24 samples (the cell types) and 
#' 4261 features (the neighbourhoods).
#'
#' \describe{
#' This is a Seurat object created by taking t(NBHDByCTMatrix)
#' as the counts.
#' }
#'
#' @source Created from t(NBHDByCTMatrix) by
#' computeNBHDVsCTSeurat()
"CTByNBHDSeurat"



## ####################################################
#' CTByNBHDSeuratExtended
#'
#' A Seurat object computed from the transpose of NBHDByCTMatrixExtended.
#' Think of cell types "expressing" (being found in) neighbourhoods.
#'
#' @format A Seurat object consisting of 16 samples (cell types) and 
#' 4261 features (the neighbourhoods).
#'
#' \describe{
#' This is a Seurat object created by taking t(NBHDByCTMatrixExtended)
#' as the counts.
#' }
#'
#' @source Created from t(NBHDByCTMatrixExtended) by
#' computeNBHDVsCTSeurat()
"CTByNBHDSeuratExtended"


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
#' edgeSeurat
#'
#' A Seurat object computed from ligandReceptorResults using 
#' computeEdgeSeurat()
#'
#' @format A Seurat object consisting of 25518 samples (edges between cells
#' and 28 features (ligand-receptor pairs).
#'
#' \describe{
#' This is a Seurat object where 
#' each point represents an edge between cells, and spatial coordinates are the 
#' centroids of edges between cells. The "expression matrix" is the 
#' binarised presence/absence of an interaction (ligand receptor pair) on an edge. 
#' }
#'
#' @source Created from ligandReceptorResults by
#' computeEdgeSeurat()
"edgeSeurat"

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


