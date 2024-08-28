
## ####################################################
## This functions focus on the analysis of spatial
## neighbourhoods.
## ####################################################

## ####################################################
#' This function computes a matrix where neighbourhoods are rows and
#' cell types are columns. The values in the matrix indicate the
#' number of cells of a given type within a neighbourhood.
#'
#' @param  spatialGraph - a spatial graph in neighbour list format.
#' @param cellTypes - named vector of cell types where names are each cell and
#' cell types are a factor
#' @return a matrix of neighbourhoods by cell types
#' @export
#' @examples
#' getExample = make.getExample()
#' clusters = getExample('clusters')
#' delaunayNeighbours = getExample('delaunayNeighbours')
#' NBHDByCTMatrix = computeNBHDByCTMatrix(delaunayNeighbours,clusters)
computeNBHDByCTMatrix = function(spatialGraph, cellTypes){
  
  spatialGraphBA = spatialGraph[,c(2,1)]
  names(spatialGraphBA) = c("nodeA","nodeB")
  spatialGraph = rbind(spatialGraph,spatialGraphBA)
  spatialGraph[,2] = cellTypes[spatialGraph[,2]]
  NBHDByCTMatrix = table(spatialGraph[,1],spatialGraph[,2])
  NBHDByCTMatrix  = as.data.frame.matrix(NBHDByCTMatrix)
  NBHDByCTMatrix  =   NBHDByCTMatrix[names(cellTypes),]
  
  return(NBHDByCTMatrix)
}

## ####################################################
#' This function creates a seurat object using a neighbourhood by cell type 
#' matrix
#' 
#' @param dataMatrix - a matrix of neighbourhoods by cell types or its 
#' transpose.
#' @param resolution - resolution for clustering (default 0.1).
#' @param npcs - number of pcs used for PCA, defaults to 10.
#' @param n.neighbors - number of neighbors used by UMAP, defaults to 30.
#' @param transpose - defaults to FALSE.
#' @param verbose - defaults to TRUE, used to limit trace if FALSE
#' @param returnType - Will return a SingleCellExperiment if this is either
#' of SCE, SingleCellExperiment or their lower-case equivalents.  Otherwise,
#' returns a Seurat object
#' @return a seurat object based on a neighbourhood by cell type matrix or its 
#' transpose, containing clusters and UMAP. This can also be a
#' SingleCellExperiment depending on the parameter returnType.
#' @import SeuratObject
#' @export
#' @examples
#' NBHDByCTMatrix = make.getExample()('NBHDByCTMatrix')
#' NBHDByCTSeurat = computeNBHDVsCTObject(NBHDByCTMatrix)
#' NBHDByCTSingleCell_sce = computeNBHDVsCTObject(NBHDByCTMatrix,returnType='SCE')
computeNBHDVsCTObject= function(dataMatrix, resolution = 0.1, 
                                npcs = 10, n.neighbors = 30L, 
                                transpose = FALSE,
                                verbose=TRUE,
                                returnType='Seurat'){
    dataMatrix = t(dataMatrix)
    NBHDSeurat = CreateSeuratObject(dataMatrix)
    NBHDSeurat[['RNA']]$data = NBHDSeurat[['RNA']]$counts
    NBHDSeurat = ScaleData(NBHDSeurat,verbose=verbose)
    NBHDSeurat = RunPCA(NBHDSeurat, assay = "RNA", 
                        features = rownames(NBHDSeurat), 
                        npcs = npcs,
                        verbose=verbose)
    if (transpose){
      NBHDSeurat = RunUMAP(NBHDSeurat,assay='RNA',
                           dims = seq_len(npcs), n.neighbors = n.neighbors,
                           verbose=verbose)
  
    } else{
      NBHDSeurat = RunUMAP(NBHDSeurat,assay='RNA',
                           features=rownames(NBHDSeurat), 
                           n.neighbors = n.neighbors,
                           verbose=verbose)
    }
    if (transpose){
        NBHDSeurat = FindNeighbors(NBHDSeurat, dims = seq_len(npcs),
                                   verbose=verbose)
    } else{
        NBHDSeurat = FindNeighbors(NBHDSeurat, 
                                   features=rownames(NBHDSeurat),
                                   verbose=verbose)
    }
    NBHDSeurat = FindClusters(NBHDSeurat,
                              resolution=resolution)

    ## Rename seurat_clusters to neighbourhood_clusters
    idx = names(NBHDSeurat@meta.data) == 'seurat_clusters'
    names(NBHDSeurat@meta.data)[idx] = 'neighbourhood_clusters'
    
    return(returnAs(NBHDSeurat,returnType))
}

## ####################################################
#' This function adds a force directed graph embedding to a seurat object
#' 
#' @param seuratObj - a seurat object of SingleCellExperiment to be
#' turned into a Seurat object
#' @param graph - which graph to extract.  Defaults to
#' paste0(f@active.assay,'_snn')
#' @param returnType - Will return a SingleCellExperiment if this is either
#' of SCE, SingleCellExperiment or their lower-case equivalents.  Otherwise,
#' returns a Seurat object
#' @return a seurat object with a "graph" dimensionality reduction. Can also
#' be a SingleCellExperiment depending on parameter returnType.
#' @importFrom igraph graph_from_adjacency_matrix layout_with_fr
#' @importFrom igraph V E V<- E<-
#' @export
#' @examples
#' NBHDByCTSeurat = make.getExample()('NBHDByCTSeurat')
#' objWithEmbedding = computeGraphEmbedding(NBHDByCTSeurat)
computeGraphEmbedding = function(seuratObj, graph=defaultGraph(seuratObj),
                                 returnType='Seurat'){
    seuratObj = acceptor(seuratObj)
    graph = seuratObj@graphs[[graph]]
    igraphGraph = igraph::graph_from_adjacency_matrix(graph)
    graphLayout = igraph::layout_with_fr(igraphGraph)
    colnames(graphLayout) = c("graph_1","graph_2")
    rownames(graphLayout) = colnames(seuratObj)
    graphDimReduc = CreateDimReducObject(embeddings = graphLayout,   
                                         key = "graph",   assay = "RNA")
    seuratObj[["graph"]] = graphDimReduc
    return(returnAs(seuratObj,returnType))
}


## ####################################################
#' For each cell type, this function looks at the neighbourhoods
#' around cells of that type and discovers the fractions of those
#' cells of each type.
#'
#' @param nbhdByCellType - A matrix whose rows are neighbourhoods
#' each denoted by the cell at their center, whose columns are
#' cell types, and whose entries are counts.
#' @param cellTypes - named vector of cell types where names are each cell and
#' cell types are a factor
#' @return A square matrix whose rownames and colnames are the
#' seurat_clusters as character strings.  Each row corresponds
#' to neighbourhoods around all cells of that type and the entries
#' give the fractions of those neighbourhoods occupied by cells
#' of each type.
#' @export
#' @examples
#' getExample = make.getExample()
#' NBHDByCTMatrix = getExample('NBHDByCTMatrix')
#' clusters = getExample('clusters')
#' cellTypesPerCellType = computeCellTypesPerCellTypeMatrix(NBHDByCTMatrix,clusters)
computeCellTypesPerCellTypeMatrix = function(nbhdByCellType,cellTypes)
{
  MM = aggregate(nbhdByCellType, list(cellTypes), sum)
  rownames(MM) = MM$Group.1
  MM = MM[,seq(from=2,to=ncol(MM))]
  MM = MM/rowSums(MM)
  MM = as.matrix(MM)
  return(MM)
}


## ####################################################
#' This function converts a matrix as found by
#' cellTypesPerCellTypeMatrix into a directed igraph
#' whose vertices correspond to seurat_clusters and whose
#' edge correspond to occupancy fraction.
#'
#' @param M - a matrix as found by cellTypesPerCellTypeMatrix. Note,
#' however, that this matrix may need to be reduced to a square matrix
#' as the matrix produced from a subset object may be missing certain
#' cell types as rows.
#' @param colours - a named vector of colours used to colour the
#' vertices of the graph.  The names are the seurat_clusters
#' as character strings.
#' @param selfEdges - a logical which determines whether to include
#' self edges.  Defaults to FALSE
#' @param minWeight - Allows one to exclude edges of low weight.
#' Defaults to 0, thus including all edges.
#' @param edgeWeighting - a parameter used to thicken the edges
#' in the display.  Defaults to 20.
#' @param edgeCurved - a parameter to set curvature of the edges.
#' Defaults to 0.2
#' @param arrowSize - a parameter to set arrow size. Defaults to 4.
#' @param arrowWidth - a parameter to set arrow width. Defaults to 4.
#' @param plotGraph - a logical which determines whether to
#' plot the graph.  Defaults to TRUE.
#' @return This returns a directed igraph whose vertices are
#' the cell types and whose arrows indicate "ownership" of
#' cells of the target type by neighbourhoods of cells of the
#' source type.  Layout is done witht the FR algorithm and
#' coordinates are found in the coords attribute of G.  If colours
#' were supplied these are found in color attribute of V(G).  Edge
#' weights and widths are found in the weight and width attributes
#' of E(G).
#' @export
#' @examples
#' getExample = make.getExample()
#' cellTypesPerCellTypeMatrix = getExample('cellTypesPerCellTypeMatrix')
#' colours = getExample('colours')
#' G = cellTypesPerCellTypeGraphFromCellMatrix(cellTypesPerCellTypeMatrix, 
#'                                    minWeight = 0.05, colours = colours)
cellTypesPerCellTypeGraphFromCellMatrix = function(M,
                                               colours=NULL,
                                               selfEdges=FALSE,
                                               minWeight=0,
                                               edgeWeighting=20,
                                               edgeCurved=0.2,
                                               arrowSize=4,
                                               arrowWidth=4,
                                               plotGraph=TRUE)
{
    idx = M >= minWeight
    M[!idx] = 0
    
    G = graph_from_adjacency_matrix(M,
                                    mode='directed',
                                    weighted=TRUE,
                                    diag=selfEdges)

    if(! is.null(colours))    
        V(G)$color = colours[names(V(G))]
    G$coords = layout_with_fr(G)
    E(G)$width = edgeWeighting * E(G)$weight

    if(plotGraph)
    {
        if(! is.null(colours))
            plot(G,
                 layout=G$coords,
                 vertex.color=V(G)$color,
                 edge.width=E(G)$width,
                 edge.curved = edgeCurved,
                 arrow.size = arrowSize,
                 arrow.width = arrowWidth)
        else
            plot(G,
                 layout=G$coords,
                 edge.width=E(G)$width,
                 edge.curved = edgeCurved,
                 arrow.size = arrowSize,
                 arrow.width = arrowWidth)
    }
    return(G)    
}

## ####################################################
#' This function takes a neighbourhood-by-cell type
#' matrix and produces a directed igraph showing the
#' fractions of cells of each type in the neighbourhoods
#' around cells of each type.
#'
#' @param nbhdByCellType - A matrix whose rows are neighbourhoods
#' each denoted by the cell at their center, whose columns are
#' cell types, and whose entries are counts.
#' @param clusters - a named vector whose names are the cells
#' and whose entries are their seurat_clusters.
#' @param colours - a named vector of colours used to colour the
#' vertices of the graph.  The names are the seurat_clusters
#' as character strings.
#' @param selfEdges - a logical which determines whether to include
#' self edges.  Defaults to FALSE
#' @param minWeight - Allows one to exclude edges of low weight.
#' Defaults to 0, thus including all edges.
#' @param edgeWeighting - a parameter used to thicken the edges
#' in the display.  Defaults to 20.
#' @param edgeCurved - a parameter to set curvature of the edges.
#' Defaults to 0.2
#' @param arrowSize - a parameter to set arrow size. Defaults to 4.
#' @param arrowWidth - a parameter to set arrow width. Defaults to 4.
#' @param plotGraph - a logical which determines whether to
#' plot the graph.  Defaults to TRUE.
#' @return This returns a directed igraph whose vertices are
#' the cell types and whose arrows indicate "ownership" of
#' cells of the target type by neighbourhoods of cells of the
#' source type.  Layout is done witht the FR algorithm and
#' coordinates are found in the coords attribute of G.  If colours
#' were supplied these are found in the color attribute of V(G).
#' Edge weights and widths are found in the weight and width
#' attributes of E(G).
#' @export
cellTypesPerCellTypeGraphFromNbhdMatrix = function(nbhdByCellType,
                                                   clusters,
                                                   colours=NULL,
                                                   selfEdges=FALSE,
                                                   minWeight=0,
                                                   edgeWeighting=20,
                                                   edgeCurved=0.2,
                                                   arrowSize=4,
                                                   arrowWidth=4,
                                                   plotGraph=TRUE)
{
    M = computeCellTypesPerCellTypeMatrix(nbhdByCellType,clusters)
    G = cellTypesPerCellTypeGraphFromCellMatrix(M,
                                            colours=colours,
                                            selfEdges=selfEdges,
                                            minWeight=minWeight,
                                            edgeWeighting=edgeWeighting,
                                            plotGraph=plotGraph,
                                            edgeCurved = edgeCurved,
                                            arrowSize = arrowSize,
                                            arrowWidth = arrowWidth)

    
    return(G)
}


## ####################################################
#' This function calculates P values for whether cell types are more frequently 
#' neighbours than expected by chance. It does this by comparison to randomised
#' neighbour graphs where edges are randomised but the degree of each node is 
#' preserved. 
#'
#' @param spatialGraph - a spatial graph in neighbour list format.
#' @param cellTypes - named vector of cell types where names are each cell and
#' cell types are a factor.
#' @param nSim - the number of randomised graphs to create for pvalue 
#' calculation.
#' @param maxTries - the maximum number of tries to remove self edges during 
#' graph randomisation. If self edges are remeining this will be reported.
#' @param verbose - whether to print trace.  Defaults to TRUE
#' @return A square matrix containing upper tail p values describing whether two 
#' cell types are more frequently found together than expected by chance.
#' @importFrom abind abind
#' @export
#' @examples
#' getExample = make.getExample()
#' delaunayNeighbours = getExample('delaunayNeighbours')
#' clusters = getExample('clusters')
#' cellTypesPerCellTypePValues = computeNeighbourEnrichment(delaunayNeighbours, 
#'                                         clusters, nSim = 10, verbose = FALSE)
computeNeighbourEnrichment = function(spatialGraph, cellTypes, nSim = 1000,
                                      maxTries = 1000,
                                      verbose=TRUE){
  results = list()
  spatialGraphOrig = spatialGraph
  NBHDByCTmatrix = computeNBHDByCTMatrix(spatialGraphOrig,cellTypes) 
  cellTypeMatrix = computeCellTypesPerCellTypeMatrix(NBHDByCTmatrix, cellTypes) 
  for (i in seq_len(nSim)){ 
    spatialGraph = spatialGraphOrig
    simGraph = randomiseGraph(spatialGraph, maxTries = maxTries)
    NBHDByCTmatrix = computeNBHDByCTMatrix(simGraph,cellTypes)
    results[[i]] = computeCellTypesPerCellTypeMatrix(NBHDByCTmatrix,cellTypes)
    if (i %% 10 == 0 & verbose){
      writeLines(as.character(i))
    }
  }
  results = lapply(results, function(x, y) y > x, y = cellTypeMatrix)
  results = abind(results, along = 3L)
  results = rowSums(results, dims = 2)
  results = abs((results - nSim)/nSim) 
  results = pmax(results,(1/nSim))
  return(results)
}


## ####################################################
#' This function takes a Seurat object and a list of
#' neighbourhoods and creates a Seurat object where the
#' columns are the neighbourhoods, the rows are are the
#' genes and the values are gene expression totals for
#' the cells in each neighbourhood
#'
#' @param f - a Seurat object with layer counts or a SingleCellExperiment
#' to be turned into a Seurat object
#' @param neighbourhoods - Neighbourhoods as given by a
#' collapsed expanded edge graph, as produced by
#' collapseNeighbourhoods. In particular, each cell should
#' appear as nodeA.
#' @param verbose - used to control trace, defaults to TRUE
#' @param returnType - Will return a SingleCellExperiment if this is either
#' of SCE, SingleCellExperiment or their lower-case equivalents.  Otherwise,
#' returns a Seurat object or SingleCellExperiment, depending on the
#' parameter returnType.
#' @return a Seurat object giving total gene expression
#' in each neighbourhood or SingleCellExperiment
#' @export
#' @examples
#' getExample = make.getExample()
#' smallXenium = getExample('smallXenium')
#' extendedNeighbours = getExample('extendedNeighbours')
#' agg = aggregateGeneExpression(smallXenium,extendedNeighbours,verbose=FALSE)
aggregateGeneExpression = function(f,neighbourhoods,verbose=TRUE,
                                         returnType='Seurat')
{
    f = acceptor(f)
    
    cells = colnames(f)
    nbhds = unique(c(neighbourhoods$nodeA,
                     neighbourhoods$nodeB))

    ## Sanity check:
    stopifnot(identical(cells[order(cells)],
                        nbhds[order(nbhds)]))

    genes = rownames(f)
    counts = FetchData(f,genes,layer='counts')
    counts = t(counts)
    counts = data.matrix(counts)

    nbhdList =nbhdsAsEdgesToNbhdsAsList(cells,
                                        neighbourhoods)
    
    C = aggregateFeatureMatrix(counts, nbhdList, rowSums)

    nbhdObj = CreateSeuratObject(counts=C)
    nbhdObj = NormalizeData(nbhdObj,verbose=verbose)
    nbhdObj = ScaleData(nbhdObj,verbose=verbose)
    nbhdObj = FindVariableFeatures(nbhdObj,verbose=verbose)
    nbhdObj = RunPCA(nbhdObj,verbose=verbose)
    nbhdObj = RunUMAP(nbhdObj,dims=seq_len(20),verbose=verbose)
    nbhdObj = FindNeighbors(nbhdObj,verbose=verbose)
    nbhdObj = FindClusters(nbhdObj,verbose=verbose)

    ## Rename the clusters:
    idx = names(nbhdObj@meta.data) == 'seurat_clusters'
    names(nbhdObj@meta.data)[idx] = 'aggregation_clusters'
    
    return(returnAs(nbhdObj,returnType))
}

## ####################################################
#' This function takes a matrix where rows are features and columns are cells,
#' and a neighbourhood list, and creates an matrix where columns are the 
#' neighbourhoods, the rows are are the features and the values are aggregated 
#' expression values for cells in each neighbourhood.
#'
#' @param M - a matrix where column names are cells and row names are 
#' features.
#' @param nbhdList - a named list with memberships of the neighbourhoods
#' of cells
#' @param aggregateFunction - a function to aggregate expression (e.g. rowSums,
#' rowMeans)
#' @return a matrix giving aggregated gene expression for a cell's neighbourhood.
#' @export
aggregateFeatureMatrix = function(M, nbhdList, aggregateFunction)
{
    cells = colnames(M)
    res = lapply(cells, function(x, M, nbhdList)
                        aggregateFunction(M[,nbhdList[[x]], drop = FALSE]), 
                 M = M, nbhdList = nbhdList)
    
    aggrM = do.call(cbind, res)
    colnames(aggrM ) = cells

    return(aggrM)
}

## ####################################################
#' This function takes a matrix where rows are features and columns are cells,
#' and a neighbourhood list, and computes Moran's I. 
#'
#' @param M - a matrix where column names are cells and row names are features.
#' @param nbhdList - a named list with memberships of the neighbourhoods
#' of cells
#' @return a matrix giving aggregated gene expression for a cell's neighbourhood.
#' @export
computeMoransI = function(M,nbhdList){
  aggrM = aggregateFeatureMatrix(M,nbhdList, rowMeans)
  means = rowMeans(M)
  zi = M - means
  zj = aggrM - means
  moransI = rowSums(zi*zj)/rowSums((zi^2))
  return(moransI)
}

## ####################################################
#' This function takes a matrix where rows are features and columns are cells,
#' and a neighbourhood list, and computes Moran's I. 
#'
#' @param obj - a Seurat object
#' @param spatialGraph - a data frame of neighbouring
#' cell pairs. 
#' @param assay - assay to pull data from, defaults to RNA.
#' @param layer - layer to pull data from, defaults to data.
#' @param nSim - number of simulations to perform for p value calculation. 
#' Defaults to 100.
#' @param verbose - whether to print trace, defaults to TRUE
#' @import SingleCellExperiment
#' @import SpatialExperiment
#' @importFrom S4Vectors SelfHits
#' @return a dataframe containing Moran's I and p values for each feature.
#' @export
#' @examples
#' getExample = make.getExample()
#' smallXenium = getExample('smallXenium')
#' delaunayNeighbours = getExample('delaunayNeighbours')
#' moransI = runMoransI(smallXenium, delaunayNeighbours, assay = "SCT", 
#' layer = "data", nSim = 10, verbose = FALSE)

runMoransI = function(obj, spatialGraph, assay = "RNA", layer = "data",
                      nSim = 100, verbose = TRUE){
    obj = acceptor(obj)
    
    spatialGraph = symmetriseNN(spatialGraph)
    
    M = as.matrix(LayerData(obj, assay = assay, layer = layer))
    nbhdList = nbhdsAsEdgesToNbhdsAsList(colnames(M),
                                         spatialGraph)
    
    moransI = computeMoransI(M,nbhdList)
    results = list()
    for (i in seq_len(nSim)){
        permuted = permuteMatrix(M)
        results[[i]] = computeMoransI(permuted,nbhdList)
        if (i %% 10 == 0 & verbose){
            writeLines(as.character(i))
        }
    }
    
    results = do.call(cbind, results)
    simResults = rowSums(moransI > results) 
    pValues = abs((simResults - nSim)/nSim) 
    pValues = pmax(pValues, (1/nSim))
    results = cbind(moransI, pValues)
    results = results[order(moransI, decreasing = TRUE),]
    results = data.frame(results)
    return(results)
}
