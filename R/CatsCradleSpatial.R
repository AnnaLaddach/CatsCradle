## ####################################################
#' This function computes a spatial graph where 
#' neighbors are identified based on Delaunay triangulation.
#' 
#' @param centroids - a dataframe containing centroids 
#' where rownames are cellnames and the first two columns
#' contain x and y coordinates respectively.
#' @return a graph in neighbour format, i.e., a data frame with
#'     columns nodeA and nodeB.
#' @import geometry
#' @import Rfast
#' @export
#' @examples
#' delaunayNeighbours = computeNeighboursDelaunay(centroids)
computeNeighboursDelaunay = function(centroids){

    ## This is confounded when centroids has extra columns:
    centroids = centroids[,1:2]
  
    ##get cell names
    cellNames = rownames(centroids)
    
    ##compute delaunay triangulation
    triangles = delaunayn(centroids)
    
    ##extract neighbours
    results = rbind(triangles[,c(1,2)],triangles[,c(1,3)],triangles[,c(2,3)])
    
    ##sort rows
    results = rowSort(results)
    
    ## remove duplicates
    results = unique(results)
    
    ##convert indices to cell names
    results[,1] = cellNames[as.numeric(results[,1])]
    results[,2] = cellNames[as.numeric(results[,2])]

    ## ## This shouldn't be, but it is:
    ## idx = results[,1] == results[,2]
    ## results = results[!idx,]
    
    ## Convert to data.frame and name as nodeA, nodeB:
    results = as.data.frame(results)
    names(results) = c('nodeA','nodeB')
    
    return(results)
}


## ####################################################
#' This function computes a spatial graph where 
#' neighbors are identified based on euclidean distance and a 
#' user defined threshold. 
#' 
#' @param centroids - a dataframe containing centroids 
#' where rownames are cellnames and columns contain x 
#' and y coordinates respectively.
#' @param threshold - a distance cut off to compute neighbours.
#' @return a graph in neighbour format, i.e., a data frame with
#'     columns nodeA and nodeB.
#' @import rdist 
#' @import reshape2 
#' @import Rfast
#' @export
#' @examples
#' euclideanNeighbours = computeNeighboursEuclidean(centroids,20)
computeNeighboursEuclidean = function(centroids, threshold){
  centroids = centroids[,c(1,2)]
  colnames(centroids) = c("x","y")
  maxX = max(centroids[,"x"])
  minX = min(centroids[,"x"])
  maxY = max(centroids[,"y"])
  minY = min(centroids[,"y"])
  XStep = (maxX - minX)/10
  YStep = (maxY - minY)/10
  results = list()
  k = 1
  
  ##calculate distances for cells in overlapping windows  
  for (i in 0:9){
    for (j in 0:9){
      x1 = i * XStep + minX 
      x2 = x1 + 2*XStep
      y1 = j * YStep + minY 
      y2 = y1 + 2*YStep
      selected = centroids[(centroids[,"x"] >= x1) & (centroids[,"x"] <= x2) & (centroids[,"y"] >= y1) & (centroids[,"y"] <= y2), ]
      distances = pdist(selected) 
      distances = reshape2::melt(distances)
      
      ##retain edges within distance threshold
      distances = distances[distances$value < threshold,]
      distances = distances[,c(1,2)]
      
      ##sort rows
      distances =  rowSort(as.matrix(distances))
      distances[,1] = rownames(selected)[distances[,1]]
      distances[,2] = rownames(selected)[as.numeric(distances[,2])]
      
      results[[k]] = distances
      k = k+1
    }
  }
  
  results = do.call(rbind, results)
  
  ##remove duplicate edges
  results = unique(results)
  colnames(results) = c("nodeA", "nodeB")
  results = results[results[,"nodeA"] != results[,"nodeB"],]
  results = as.data.frame(results)
  return(results)
}

## ####################################################
## Not exported:
extractCells = function(NN)
{
    cells = unique(c(NN$nodeA,NN$nodeB))
    cells = cells[order(cells)]
    return(cells)
}


## ####################################################
#' neighbourhoodDiameter
#'
#' This function takes a list of neighbourhoods and and the
#' centroids of the cells and finds their diameters, i.e.,
#' for each neighbourhood, the maximum distance between.
#'
#' @param neighbourhoods - a list of neighbourhoods as
#' returned by nbhdsAsEdgesToNbhdsAsList
#' @param centroids - the centroids of the cells
#' @import pracma
#' @return a named numeric.  The names are the names
#' of the list neighbourhoods and the values are the
#' maximum distance within each neighbourhood
#' @export
#' @examples
#' cells = unique(c(delaunayNeighbours$nodeA,delaunayNeighbours$nodeB))
#' nbhds = nbhdsAsEdgesToNbhdsAsList(cells,delaunayNeighbours)
#' diameters = neighbourhoodDiameter(nbhds[1:100],centroids)
neighbourhoodDiameter = function(neighbourhoods,centroids)
{
  rownames(centroids) = centroids$cell
  nbhds = names(neighbourhoods)
  diameters = c()
  for(nbhd in nbhds)
  {
    theseCells = neighbourhoods[[nbhd]]
    N = length(theseCells)
    M = matrix(0,nrow=N,ncol=2)
    M[1:N,1:2] = as.matrix(centroids[theseCells,1:2])
    D = distmat(M,M)
    
    diameters[nbhd] = max(D)
  }
  
  return(diameters)
}


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
#' NBHDByCTMatrix = computeNBHDByCTMatrix(delaunayNeighbours,
#'                                          clusters)
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
#' @import Seurat
#' @import SeuratObject
#' @export
#' @examples
#' NBHDByCTSeurat = computeNBHDVsCTSeurat(NBHDByCTMatrix)
#' NBHDByCTSeurat_sce = computeNBHDVsCTSeurat(NBHDByCTMatrix,returnType='SCE')
computeNBHDVsCTSeurat= function(dataMatrix, resolution = 0.1, 
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
                           dims = 1:npcs, n.neighbors = n.neighbors,
                           verbose=verbose)
  
    } else{
      NBHDSeurat = RunUMAP(NBHDSeurat,assay='RNA',
                           features=rownames(NBHDSeurat), 
                           n.neighbors = n.neighbors,
                           verbose=verbose)
    }
    if (transpose){
        NBHDSeurat = FindNeighbors(NBHDSeurat, dims = 1:npcs,
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
#' @import Seurat  
#' @import igraph
#' @export
#' @examples
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
#' This function takes a nearest neighbour graph and a radius
#' and calculates nth degree neighbour graphs where max(n) == radius
#'
#' @param spatialGraph - a nearest neighbour graph
#' @param n - the maximum degree to calculate a neighbour graph with edges 
#' connecting vertices of degree n for.
#' @return A named list of neighbour graphs, where each graph contains edges 
#' connecting vertices of degree n. Each graph is named according to degree n.
#' @import data.table
#' @export
#' @examples
#' extendedNeighboursList = getExtendedNBHDs(delaunayNeighbours, 4)
getExtendedNBHDs = function(spatialGraph, n){
  spatialGraph = data.table(spatialGraph)
  
  spatialGraphR = spatialGraph
  names(spatialGraphR) = c("nodeB","nodeA")
  spatialGraph = rbind(spatialGraph,spatialGraphR[,c(2,1)])
  neighbours = list()
  neighbours[[1]] = spatialGraph
  
  for (i in (2:n)){
    writeLines(paste('radius',i))
    graph = merge(neighbours[[i-1]], neighbours[[1]], by.x = "nodeB", 
                  by.y = "nodeA", allow.cartesian = TRUE)
    graph = graph[,c("nodeA","nodeB.y")]
    names(graph) = c("nodeA","nodeB")
    graph = unique(graph)
    orig = c(paste0(neighbours[[i-1]]$nodeB,"_",neighbours[[i-1]]$nodeA))
    if (i > 2){
      orig = c(orig,paste0(neighbours[[i-2]]$nodeB,"_",neighbours[[i-2]]$nodeA))
    }
    graph = graph[graph$nodeA != graph$nodeB,]
    new = paste0(graph$nodeA,"_",graph$nodeB)
    graph = graph[!(new %in% orig),]
    neighbours[[i]] = graph
  }
  return(neighbours)
}



## ####################################################
#' This function takes an expanded neighbourhood list and collapses it to a 
#' nearest neighbourhood graph where all neighbours of degree <= n in the 
#' original graph are considered first neighbours.
#'
#' @param extendedNeighboursList - the results of getExtendedNBHDs()
#' @param n - the maximum degree to connect neighbours. Defaults to the maximum 
#' degree neighbourhoods were expanded to in the results of getExtendedNBHDs().
#' @return a graph in neighbour format, i.e., a data frame with
#'     columns nodeA and nodeB, where nodes that were originally of degree <= n 
#'     are connected.
#' @export
#' @examples
#' extendedNeighbours = collapseExtendedNBHDs(extendedNeighboursList, 4)
collapseExtendedNBHDs = function(extendedNeighboursList, n = length(extendedNeighboursList)){
  
  collapsedGraph = as.data.frame(do.call(rbind,extendedNeighboursList[1:n]))

  return(collapsedGraph)
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
#' @import dplyr
#' @export
#' @examples
#' cellTypesPerCellType = computeCellTypesPerCellTypeMatrix(NBHDByCTMatrix,
#'                                                      clusters)
computeCellTypesPerCellTypeMatrix = function(nbhdByCellType,cellTypes)
{
  MM = aggregate(nbhdByCellType, list(cellTypes), sum)
  rownames(MM) = MM$Group.1
  MM = MM[,2:ncol(MM)]
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
#' coordinates are found in G$coords.  If colours were supplied
#' these are found in V(G)$color.  Edge weights and widths are
#' found in E(G)$weight and E(G)$width.
#' @export
#' @examples
#' G = cellTypesPerCellTypeGraphFromCellMatrix(cellTypesPerCellTypeMatrix, 
#' minWeight = 0.05, colours = colours)
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
#' coordinates are found in G$coords.  If colours were supplied
#' these are found in V(G)$color.  Edge weights and widths are
#' found in E(G)$weight and E(G)$width. 
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
    M = cellTypesPerCellTypeMatrix(nbhdByCellType,clusters)
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
#' This function performs degree-preserving randomisation of neighbour graphs. 
#'
#' @param spatialGraph - a spatial graph in neighbour list format.
#' @param maxTries - the maximum number of tries to remove self edges during 
#' graph randomisation. If self edges are remeining this will be reported.
#' @return A randomised graph where degree from the original graph is preserved.
randomiseGraph = function(spatialGraph, maxTries = 1000){
  n = nrow(spatialGraph)
  toFlip = sample(1:n, size = round(n/2))
  simGraph = spatialGraph[toFlip,c(2,1)]
  names(simGraph) =  c("nodeA","nodeB")
  simGraph = rbind(simGraph, spatialGraph[!(1:n %in% toFlip),])
  simGraph[,2] = sample(simGraph[,2])
  selfEdge = which(simGraph$nodeA == simGraph$nodeB)
  i = 1
  while((length(selfEdge) > 0) & (i <= maxTries)){
    nodeB = simGraph$nodeB
    simGraph$nodeB[selfEdge] = nodeB[(selfEdge + 1) %% n]
    simGraph$nodeB[(selfEdge + 1) %% n] = nodeB[selfEdge]
    selfEdge = which(simGraph$nodeA == simGraph$nodeB)
    i = i + 1
  }
  nSelf = length(selfEdge)
  if (nSelf > 0){
    print(paste0(nSelf, " self edges not removed"))
  }
  return(simGraph)
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
#' @import abind
#' @export
#' @examples
#' cellTypesPerCellTypePValues = computeNeighbourEnrichment(delaunayNeighbours, 
#' clusters, nSim = 10, verbose = FALSE)
computeNeighbourEnrichment = function(spatialGraph, cellTypes, nSim = 1000,
                                      maxTries = 1000,
                                      verbose=TRUE){
  results = list()
  spatialGraphOrig = spatialGraph
  NBHDByCTmatrix = computeNBHDByCTMatrix(spatialGraphOrig,cellTypes) 
  cellTypeMatrix = computeCellTypesPerCellTypeMatrix(NBHDByCTmatrix, cellTypes) 
  for (i in 1:nSim){ 
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
#' This function retrieves the Nichenetr ligand-
#' receptor network for mouse or human.
#'
#' @param species - either 'human' or 'mouse'
#' @return This returns a data frame whose first two
#' columns are from and to, i.e., ligand and receptor.
#' These are derived from the nichenetr ligand receptor
#' networks.
#' @export
#' @examples
#' lrn = getLigandReceptorNetwork('human')
getLigandReceptorNetwork = function(species)
{
    stopifnot(species %in% c('human','mouse'))
    
    if(species == 'human')
        return(humanLRN)

    if(species == 'mouse')
        return(mouseLRN)
}

## ####################################################
#' This functions takes an Seurat object, its species
#' and a ligand receptor network and subsets the ligand
#' receptor network to those pairs that occur in the
#' panel
#'
#' @param obj - a Seurat object or SingleCellExperiment to
#' be converted to a Seurat object
#' @param species - either 'human' or 'mouse'
#' @param lrn - a ligand-receptor network, i.e., a
#' data frame with columns from and to.  By default, it
#' retrieves the nichenetr ligand receptor network
#' @return This returns a data frame with columns ligand and
#' receptor
#' @export
#' @examples 
#' lrPairs = getLigandReceptorPairsInPanel(smallXenium, "mouse")
getLigandReceptorPairsInPanel = function(obj,species,
                                         lrn = getLigandReceptorNetwork(species))
{
    stopifnot(species %in% c('mouse','human'))

    obj = acceptor(obj)
  
  ## The panel
  panel = rownames(obj)
  panel = str_replace(panel,"-",".")
  lrn$from = str_replace(lrn$from,"-",".")
  lrn$to = str_replace(lrn$to,"-",".")
  
  ## Ligands and receptors:
  pairs = paste(lrn$from,lrn$to,sep='-')
  ligands = unique(lrn$from)
  receptors = unique(lrn$to)
  ligandsFound = intersect(ligands,panel)
  receptorsFound = intersect(receptors,panel)
  
  panelPairs = c()
  for(a in ligandsFound)
    for(b in receptorsFound)
      panelPairs = c(panelPairs,paste(a,b,sep='-'))
  idx = panelPairs %in% pairs
  pairsFound = panelPairs[idx]
  
  a = str_split(pairsFound,'-')
  pairsFoundDF = data.frame(ligand=unlist(lapply(a,function(x) return(x[1]))),
                            receptor=unlist(lapply(a,function(x) return(x[2]))))
  
  return(pairsFoundDF)
}


## ####################################################
#' This function takes a binarised expression matrix, a set of ligand receptor
#' pairs and a set of edges denoting neighbouring cells and
#' annotates these with the ligand receptor interactions taking
#' place on those edges in each direction.
#'
#' @param M - a binarised expression matrix where rows are genes and columns
#' are cells.
#' @param pairDF - a data frame giving the ligand-receptor pairs
#' @param spatialGraph - a data frame of neighbouring
#' cell pairs.  Note that each row is a directed edge (A,B) so
#' that this data frame should have both the edge (A,B) and the
#' edge (B,A)
#' @return This returns a data frame whose first two columns give
#' the neighbouring cells.  Each of the remaining columns is a logical
#' corresponding to a ligand-receptor pair telling whether the ligand
#' is expressed in the first cell and the receptor is expressed in the
#' second cell.
getInteractionsOnEdges = function(M,pairDF,spatialGraph)
{
  ## Find the interactions on the edges:
  edges = spatialGraph
  
  for(i in 1:nrow(pairDF))
  {
    tag = paste(pairDF$ligand[i],pairDF$receptor[i],sep='-')
    edges[,tag] = (M[pairDF$ligand[i],edges$nodeA] &
                     M[pairDF$receptor[i],edges$nodeB])
  }
  
  return(edges)
}

## ####################################################
#' This function permutes the rows of a matrix.
#'
#' @param M - a binarised expression matrix where rows are genes and columns
#` are cells.
#' @return This returns a matrix in which the values have been permuted within
#' rows.
permuteMatrix = function(M){
  n = ncol(M)
  for (i in 1:nrow(M)){
    M[i,] = M[i,sample(n)]
  }
  return(M)
}


## ####################################################
#' This functions retrieves an expression matrix from a seurat object and 
#' binarises it.
#'
#' @param obj - a Seurat object or SingleCellExperiment to be
#' turned into a Seurat object
#' @param cutoff - a cutoff for binarisation. Defaults to 0.
#' @param layer - layer to fetch data from. Defaults to count.
#' @return A binarised expression matrix where rows are genes and columns are 
#' cells.
getBinarisedMatrix = function(obj, cutoff = 0, layer = 'count'){
    obj = acceptor(obj)
    M = FetchData(obj,rownames(obj),layer='count')
    M = data.matrix(t(M))
    cutoff = 0
    M = M > cutoff
    rownames(M) = str_replace(rownames(M),"-",".")
    return(M)
}

## ####################################################
#' This function takes a listing of the neighbouring
#' cells together with the presence or absence of each
#' ligand-receptor pair on each edge and produces a count
#' showing for each cell, how many neighbours it has with
#' that interaction either as source or as target
#'
#' @param edges - A data frame of neighbouring cells
#' together with their interactions as produced by
#' getInteractionsOnEdges()
#' @param sourceOrTarget - a character, either 'source' or
#' 'target' telling which direction of interaction to count
#' @return This returns a data frame with one row for each
#' cell and a column giving the name of that cell and the
#' other columns giving the counts of interactions that it
#' has with its neighbours.
#' @export
countLRInteractionsPerCell = function(edges,sourceOrTarget)
{
    stopifnot(sourceOrTarget %in% c('source','target'))
    
    if(sourceOrTarget == 'source')
        by = factor(edges$nodeA)
    if(sourceOrTarget == 'target')
        by = factor(edges$nodeB)

    edges = edges[,3:ncol(edges)]
    interactionCountDF = aggregate(edges,by=list(by),FUN=sum) 

    interactionCountDF$Group.1 =
        str_replace(interactionCountDF$Group.1,'Group.','')
    rownames(interactionCountDF) = interactionCountDF$Group.1
    names(interactionCountDF)[1] = 'cell'
     
    return(interactionCountDF)
}


## ####################################################
#' This takes a data frame of interaction counts as found
#' by countLRInteractionsPerCell(), the underlying Seurat object
#' and the neighbourhood Seurat object and annotates the counts
#' with the cell type and the neighbourhood type corresponding
#' to the cells of the interaction counts.
#'
#' @param interactionCounts - as found by countLRInteractionsPerCell()
#' @param obj - a Seurat object, or SingleCellExperiment to be turned
#' into a Seurat object
#' @param nbhdObj - a neighbourhood x cell type Seurat object or a
#' SingleCellExperiment to be turned into a Seurat object
#' @return This returns the interaction counts annotated with the
#' cell type and neighbourhood type of each cell.
#' @export
annotateLRInteractionCounts = function(interactionCounts,obj,nbhdObj)
{
    obj = acceptor(obj)
    nbhdObj = acceptor(nbhdObj)
    
    ## Get nbhd and cell types:
    annotated = data.frame(cell=colnames(obj),
                           cellType=obj$seurat_clusters)
    rownames(annotated) = annotated$cell

    ## Append nbhd type:
    annotated$nbhdType = nbhdObj$seurat_clusters[annotated$cell]

    ## Bung in the interaction counts:
    pairs = names(interactionCounts)[2:ncol(interactionCounts)]

    annotated[,pairs] = 0
    annotated[rownames(interactionCounts),pairs] =
        interactionCounts[rownames(interactionCounts),pairs]

    return(annotated)
}





## ####################################################
#' Given a seurat object, a spatial graph, clusters and species this function 
#' identifies ligand-receptor interactions between neighbouring cells, 
#' identifies ligand-receptor interactions within and between clusters and 
#' calculates whether these are observed more frequently than expected by 
#' chance.
#'
#' @param obj - a Seurat object
#' @param spatialGraph - a data frame of neighbouring
#' cell pairs. 
#' @param clusters - named vector of clusters where names are each cell and
#' clusters are a factor
#' @param species - either 'human' or 'mouse'
#' @param nSim - number of simulations to perform for p value calculation.
#' @param lrn - a ligand-receptor network, i.e., a
#' data frame with columns from and to.  By default, it
#' retrieves the nichenetr ligand receptor network
#' @param verbose - whether to print trace, defaults to TRUE
#' @return A list containing:
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
#' @export
#' @examples
#' performLigandReceptorAnalysis(smallXenium, delaunayNeighbours, 
#'                                       "mouse", clusters, nSim = 10,
#'                                        verbose=FALSE)

performLigandReceptorAnalysis = function(obj, spatialGraph, species, clusters,
                                      nSim = 1000, 
                                      lrn = getLigandReceptorNetwork(species),
                                      verbose=TRUE){
  
  #symmetrise spatial graph
  spatialGraphBA = spatialGraph[,c(2,1)]
  names(spatialGraphBA) = c("nodeA","nodeB")
  spatialGraph = rbind(spatialGraph,spatialGraphBA)
  spatialGraph = unique(spatialGraph)
  
  #get ligand receptor pairs
  lrPairs = getLigandReceptorPairsInPanel(obj, species)
  
  #get binarised expression matrix for ligand receptor pairs
  M = getBinarisedMatrix(obj)
  M = M[rownames(M) %in% c(lrPairs$ligand,lrPairs$receptor),]
  
  #get ligand receptor interactions between cells
  interactionsOnEdges = getInteractionsOnEdges(M,lrPairs,spatialGraph)
  
  #annotate with clusters
  interactionsOnEdges = cbind(clusters[interactionsOnEdges$nodeA],
                              clusters[interactionsOnEdges$nodeB], 
                              interactionsOnEdges)
  names(interactionsOnEdges)[1:2] = c("clusterA","clusterB") 
  
  #get sum of interactions within and between clusters
  pair = 
    paste0(interactionsOnEdges$clusterA,  "-", interactionsOnEdges$clusterB) 
  totalInteractionsByCluster = 
    aggregate(interactionsOnEdges[,5:ncol(interactionsOnEdges)], list(pair), 
              sum)
  rownames(totalInteractionsByCluster) = totalInteractionsByCluster$Group.1
  totalInteractionsByCluster = 
    totalInteractionsByCluster[,2:ncol(totalInteractionsByCluster)]
  
  #get total edges per cluster pair
  totalEdges = table(pair)
  
  totalInteractionsByCluster = cbind(totalEdges[rownames(totalInteractionsByCluster)],
                                totalInteractionsByCluster)
  colnames(totalInteractionsByCluster)[1:2] = c("clusterPair", "totalEdges")
  
  meanInteractionsByCluster = totalInteractionsByCluster[,3:ncol(totalInteractionsByCluster)]/totalInteractionsByCluster$totalEdges
  meanInteractionsByCluster = cbind(totalInteractionsByCluster[,1:2], meanInteractionsByCluster)
  
  
  
  #perform simulations
  results = list()
  for (i in 1:nSim){
    permuted = permuteMatrix(M)
    sim = getInteractionsOnEdges(permuted,lrPairs,spatialGraph)
    sim = aggregate(sim[,3:ncol(sim)], list(pair), sum)
    rownames(sim) = sim$Group.1
    results[[i]] = sim[,2:ncol(sim)]
    if (i %% 10 == 0 & verbose){
      writeLines(as.character(i))
    }
  }
  
  #calculate summary statistics for simulation results
  results = lapply(results, function(x, y) y > x, y = 
                     totalInteractionsByCluster[,3:ncol(totalInteractionsByCluster)])
  results = abind(results, along = 3L)
  simResults = rowSums(results, dims = 2)
  rownames(simResults) = rownames(totalInteractionsByCluster)
  pValues = abs((simResults - nSim)/nSim) 
  pValues = pmax(pValues, (1/nSim))
  return(list("interactionsOnEdges" = interactionsOnEdges, 
              "totalInteractionsByCluster" = totalInteractionsByCluster,
              "meanInteractionsByCluster" = meanInteractionsByCluster,
              "simResults" = as.data.frame(simResults),
              "pValues" = as.data.frame(pValues)))
}

    
## ####################################################
#' This function takes ligandReceptorResults and plots a heatmap of -log10(pvalues).
#'
#' @param ligandReceptorResults - as returned by performLigandReceptorAnalysis()
#' @param clusters - named vector of cell types where names are each cell and
#' clusters are a factor
#' @param colours - a named list of colours where names are clusters. If not 
#' specified the default pheatmap colour scheme will be used.
#' @param  pValCutoffClusterPair - a cutoff for showing interactions between two
#' clusters. A cluster pair must have at least one ligand-receptor interaction
#' pvalue <  pValCutoffClusterPair. Defaults to 0.05.
#' @param  pValCutoffLigRec  - a cutoff for showing interactions between a 
#' ligand and receptor. At least one cluster pair must have 
#' pvalue <  pValCutoffLigRec for ligand-receptor pair. Defaults to 0.05.
#' @param  labelClusterPairs - show labels for cluster pairs. Defaults to TRUE.
#' @import Rfast
#' @import pheatmap
#' @return matrix of -log10(pvalues) that underlies the heatmap.
#' @export
#' @examples
#' ligRecMatrix = makeLRInteractionHeatmap(ligandReceptorResults, 
#' clusters, colours = colours, labelClusterPairs = FALSE)
makeLRInteractionHeatmap = function(ligandReceptorResults,
                                  clusters,
                                  colours = c(),
                                  pValCutoffClusterPair = 0.05, 
                                  pValCutoffLigRec = 0.05,
                                  labelClusterPairs = TRUE)
{
  pValues = as.matrix(ligandReceptorResults$pValues)
  selectedPValues = pValues[rowMins(pValues, value = TRUE) < pValCutoffClusterPair,
                                             colMins(pValues, value = TRUE) < pValCutoffLigRec]
  negLog10PValues = -log10(selectedPValues)
  rowAnno = str_split_fixed(rownames(selectedPValues), pattern = "-", 2)
  rowAnno = as.data.frame(rowAnno)
  names(rowAnno) = c("sender","receiver")
  rowAnno$sender = factor(rowAnno$sender, levels = levels(clusters))
  rowAnno$receiver = factor(rowAnno$receiver, levels = levels(clusters))
  rownames(rowAnno) = rownames(selectedPValues)
  rowAnno = rowAnno[,c("receiver","sender")]
  if (length(colours) > 0){
  pheatmap(negLog10PValues, annotation_row = rowAnno, annotation_colors = list("sender" = colours, 
                                                                    "receiver" = colours),
           show_rownames = labelClusterPairs)
  } else{
    pheatmap(negLog10PValues, annotation_row = rowAnno, show_rownames = labelClusterPairs)
  }
  return(negLog10PValues)
}


## ####################################################
#' This function takes ligandReceptorResults and plots a heatmap of the total 
#' number of ligand receptor interactions between clusters.
#'
#' @param ligandReceptorResults - as returned by performLigandReceptorAnalysis()
#' @param clusters - named vector of cell types where names are each cell and
#' clusters are a factor
#' @param type - "total" or "mean" to plot raw total interactions or mean interactions per edge.
#' @param  logScale - plot heatmap using log scale (defaults to T)
#' @import pheatmap
#' @return matrix of total ligand receptor interactions that underlies the heatmap.
#' @export
#' @examples 
#' cellTypePerCellTypeLigRecMatrix = 
#' makeSummedLRInteractionHeatmap(ligandReceptorResults, clusters, "mean")
makeSummedLRInteractionHeatmap = function(ligandReceptorResults, clusters, type, logScale = TRUE){ 
  if (type == "total"){
    interactionsByCluster = ligandReceptorResults$totalInteractionsByCluster
  } 
  if (type == "mean"){
    interactionsByCluster = ligandReceptorResults$meanInteractionsByCluster
  } 
  summedInteractionsByCluster = rowSums(interactionsByCluster[,3:ncol(interactionsByCluster)])
  pair = str_split_fixed(names(summedInteractionsByCluster), pattern = "-", 2)
  summedInteractionsByCluster = as.data.frame(cbind(pair,summedInteractionsByCluster))
  colnames(summedInteractionsByCluster) = c("Sender", "Receiver", "nInteractions")
  summedInteractionsByCluster$nInteractions = as.numeric(summedInteractionsByCluster$nInteractions)
  clusterNames = levels(clusters)
  nClusters = length(clusterNames)
  summedInteractionsByClusterMatrix = matrix(0, ncol = nClusters, nrow = nClusters)
  for (i in 1:nClusters){
    for (j in 1:nClusters){
      value = summedInteractionsByCluster$nInteractions[(summedInteractionsByCluster$Sender == clusterNames[i]) & (summedInteractionsByCluster$Receiver == clusterNames[j])]
      if (length(value) > 0){
        summedInteractionsByClusterMatrix[i,j] = value
      } 
    }  
  }
  
  colnames(summedInteractionsByClusterMatrix) = clusterNames
  rownames(summedInteractionsByClusterMatrix) = clusterNames
  if (logScale){
    pheatmap(log(summedInteractionsByClusterMatrix +1))
    return(log(summedInteractionsByClusterMatrix +1))
  } else{
    pheatmap((summedInteractionsByClusterMatrix))
    return(summedInteractionsByClusterMatrix)
  }
}

## ####################################################
#' This function takes interactionResults and creates a seurat object where 
#' each point represents an edge between cells, and spatial coordinates are the 
#' centroids of edges between cells. The "expression matrix" is the 
#' binarised presence/absence of an interaction (ligand receptor pair) on an edge. 
#'
#' @param ligandReceptorResults - as returned by performLigandReceptorResultsAnalysis()
#' @param centroids - a dataframe containing centroids 
#' where rownames are cellnames and the first two columns
#' contain x and y coordinates respectively.
#' @param npcs - number of pcs used for PCA, defaults to 10
#' @param returnType Determines whether to return a Seurat object or a
#' SpatialExperiment.  Will do the later if this is set to either SCE,
#' SingleCellExperiment or lower case versions of either.
#' @import Seurat
#' @return This returns a seurat object where 
#' each point represents an edge between cells, and spatial coordinates are the 
#' centroids of edges between cells. The "expression matrix" is the 
#' binarised presence/absence of an interaction (ligand receptor pair) on an edge.
#' Depending on the parameter returnType, this can alternatively be returned as
#' a SpatialExperiment.
#' @export
#' @examples
#' edgeSeurat = computeEdgeSeurat(ligandReceptorResults, centroids)
computeEdgeSeurat = function(ligandReceptorResults, centroids, npcs = 10,
                             returnType='Seurat'){
  interactionsOnEdges = ligandReceptorResults$interactionsOnEdges
  rownames(interactionsOnEdges) = paste0(interactionsOnEdges$nodeA, "-", interactionsOnEdges$nodeB)
  interactionsOnEdgesMat = as.matrix(interactionsOnEdges[,5:ncol(interactionsOnEdges)])
  interactionsOnEdgesMat= 1 * interactionsOnEdgesMat
  edgeSeurat = CreateSeuratObject(t(interactionsOnEdgesMat), meta = interactionsOnEdges[,1:4])
  edgeCoords = as.data.frame(cbind(centroids[interactionsOnEdges$nodeA, 1:2], 
                                  centroids[interactionsOnEdges$nodeB, 1:2]))
  
  edgeCoords$edgeX = 0.6 * edgeCoords[,1] + 0.4 * edgeCoords[,3]
  edgeCoords$edgeY = 0.6 * edgeCoords[,2] + 0.4 * edgeCoords[,4] 
  
  
  edgeCentroidDF = data.frame(
    x = edgeCoords$edgeX,
    y = edgeCoords$edgeY,
    cell = colnames(edgeSeurat),
    stringsAsFactors = FALSE
  )
  
  centroidData <- list(
    "centroids" = CreateCentroids(edgeCentroidDF)
  )
  coords = CreateFOV(
    coords = centroidData,
    type = c("centroids"),
    assay = "RNA"
  )
  
  edgeSeurat[["global"]] = coords
  return(returnAs(edgeSeurat,returnType,spatial=TRUE))
}



## ####################################################
#' nbhdsAsEdgesToNbhdsAsList
#'
#' This function takes a set of neighbourhoods given
#' by edges and turns it into a named list giving the
#' memberships of each neighbourhood
#'
#' @param cells - The cells whose neighbourhoods to
#' extract. 
#' @param neighbourhoods - neighbourhoods given as a
#' data frame with columns nodeA and nodeB, for example
#' the output of collapseNeighbourhoods
#' @return a named list with memberships of the neighbourhoods
#' of cells
#' @export
#' @examples
#' cells = unique(c(delaunayNeighbours$nodeA,delaunayNeighbours$nodeB))
#' nbhdsList = nbhdsAsEdgesToNbhdsAsList(cells,delaunayNeighbours)
nbhdsAsEdgesToNbhdsAsList = function(cells,
                                     neighbourhoods)
{
    nbhdList = list()
    for(cell in cells)
    {
        nbhdList[[cell]] = c(cell,
                             neighbourhoods$nodeB[neighbourhoods$nodeA == cell])
    }

    return(nbhdList)
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
#' @import Seurat
#' @export
#' @examples 
#' agg = aggregateSeuratGeneExpression(smallXenium,extendedNeighbours,
#' verbose=FALSE)
aggregateSeuratGeneExpression = function(f,neighbourhoods,verbose=TRUE,
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
    
    C = aggregateFeatureMatrix(counts, nbhdList, rowsums)

    nbhdObj = CreateSeuratObject(counts=C)
    nbhdObj = NormalizeData(nbhdObj,verbose=verbose)
    nbhdObj = ScaleData(nbhdObj,verbose=verbose)
    nbhdObj = FindVariableFeatures(nbhdObj,verbose=verbose)
    nbhdObj = RunPCA(nbhdObj,verbose=verbose)
    nbhdObj = RunUMAP(nbhdObj,dims=1:20,verbose=verbose)
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
    res = lapply(cells, function(x, M, nbhdList) aggregateFunction(M[,nbhdList[[x]], drop = FALSE]), 
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
## #' @param aggregateFunction - a function to aggregate expression (e.g. rowSums,
#' rowMeans)

#' @return a matrix giving aggregated gene expression for a cell's neighbourhood.
#' @export
computeMoransI = function(M,nbhdList){
  aggrM = aggregateFeatureMatrix(M,nbhdList, rowMeans)
  means = rowmeans(M)
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
#' @import Seurat
#' @import SingleCellExperiment
#' @import SpatialExperiment
#' @importFrom S4Vectors SelfHits
#' @return a dataframe containing Moran's I and p values for each feature.
#' @export
#' @examples
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
    for (i in 1:nSim){
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


## ####################################################
#' This function takes a spatial graph and computes a new spatial graph where
#' edges become nodes and A-B edges (in the original graph) become connected to
#' all A- edges and all B- edges. 
#' 
#' @param spatialGraph - a data frame of neighbouring edge pairs. 
#' @param selfEdges - a logical determining whether to include self edges. 
#' Defaults to False.
#' @return a graph in neighbour format  where edges in the original graph 
#' become nodes and A-B edges (in the original graph) become connected to
#' all A- edges and all B- edges. 
#' @import data.table
#' @export
#' @examples
#' edgeNeighbours = computeEdgeGraph(delaunayNeighbours)

computeEdgeGraph = function(spatialGraph, selfEdges = FALSE){
  spatialGraph = data.table(spatialGraph)
  spatialGraph$edge = paste0(spatialGraph$nodeA, "-", spatialGraph$nodeB)
  spatialGraphEdgesA = merge(spatialGraph, spatialGraph[,c(1,3)], by = "nodeA", allow.cartesian = TRUE)
  spatialGraphEdgesB = merge(spatialGraph, spatialGraph[,c(2,3)], by = "nodeB", allow.cartesian = TRUE)
  
  spatialGraphEdges = rbind(spatialGraphEdgesA[,c(3,4)],spatialGraphEdgesB[,c(3,4)])
  names(spatialGraphEdges) = c("nodeA","nodeB")
  if (! selfEdges){
    spatialGraphEdges = spatialGraphEdges[spatialGraphEdges$nodeA != spatialGraphEdges$nodeB,]
  }
  
  return(spatialGraphEdges)
}


## ####################################################
## The following functions are used internally to manage
## conversions between Seurat objects on the one hand and
## SingleCellExperiments and SpatialExperiments on the
## other.
## ####################################################
acceptor = function(obj)
{
    ## We only convert if we need to:
    if(isa(obj,'SingleCellExperiment'))
        return(SCEtoSeurat(obj))

    return(obj)
}

## ####################################################
returnAs = function(obj,returnType,spatial=FALSE)
{
    if(returnType %in% c('SCE',
                         'sce',
                         'SingleCellExperiment',
                         'singlecellexperiment'))
        return(SeuratToSCE(obj,spatial))

    return(obj)
}

## ####################################################
SCEtoSeurat = function(sce)
{
    f = as.Seurat(sce)

    numCells = ncol(sce)
    for(n in colPairNames(sce))
    {
        M = matrix(0,nrow=numCells,ncol=numCells)
        rownames(M) = colnames(sce)
        colnames(M) = colnames(sce)

        df = colPair(sce,n)
        df = as.data.frame(df)
        M[df$from,df$to] = df$value
        M[df$to,df$from] = df$value
        graph = as.Graph(M)
        f@graphs[[n]] = graph
    }

    return(f)
}

## ####################################################
SeuratToSCE = function(f,spatial)
{
    sce = as.SingleCellExperiment(f)
    for(n in names(f@graphs))
    {
        NN = getNearestNeighborListsSeurat(f,n)
        numPairs = nrow(NN)

        a = 1:ncol(sce)
        names(a) = colnames(sce)
        cell1 = a[NN$nodeA]
        cell2 = a[NN$nodeB]

        names(cell1) = NULL
        names(cell2) = NULL
        
        colPair(sce,n) = SelfHits(cell1,
                                  cell2,
                                  nnode=ncol(sce),
                                  value=NN$weight)
    }

    if(spatial)
    {
        sce = SpatialExperiment(sce)

        ## Copy in the centroids:
        centroids = GetTissueCoordinates(f)
        rownames(centroids) = centroids$cell
        centroids = centroids[,1:2]
        centroids = data.matrix(centroids)

        spatialCoords(sce) = centroids
    }
    return(sce)
}
