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
#' smallTriangulation = computeNeighboursDelaunay(smallCentroids)
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
## #' @examples

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
#' This function extracts neighbourhoods from a spatial graph in neighbour list
#' format
#' 
#' @param  spatialGraph - a spatial graph in neighbour list format.
#' @param cellNames - defaults to the names in the spatial graph
#' @param addSelf - logical determining whether to include the cell at
#' the center of the neighbourhood.  Defaults to FALSE
#' @return a named list of neighbourhoods where a neighborhood is a set of a 
#' cell's nearest neighbours. 
#' @export
#' @examples
#' smallNbhds = computeNeighbourhoods(smallDelaunayTriangulation,addSelf=TRUE)
computeNeighbourhoods = function(spatialGraph,
                                 cellNames=extractCells(spatialGraph),
                                 addSelf = F){
  neighbourhoods = list()
  for (cell in cellNames){
    neighbours = 
      c(spatialGraph[,2][cell == spatialGraph[,1]],
        spatialGraph[,1][cell == spatialGraph[,2]])
    neighbourhoods[[cell]] = neighbours
    if (addSelf){
      neighbourhoods[[cell]] = c(neighbourhoods[[cell]], cell)
    }
  }
 return(neighbourhoods) 
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
#' smallNbhdMatrix = computeNeighbourhoodByCTMatrix(smallDelaunayTriangulation,
#'                                          smallXenium$seurat_clusters)
computeNeighbourhoodByCTMatrix = function(spatialGraph, cellTypes){
  
  spatialGraphBA = spatialGraph[,c(2,1)]
  names(spatialGraphBA) = c("nodeA","nodeB")
  spatialGraph = rbind(spatialGraph,spatialGraphBA)
  spatialGraph[,2] = cellTypes[spatialGraph[,2]]
  neighbourhoodbyCTmatrix = table(spatialGraph[,1],spatialGraph[,2])
  neighbourhoodbyCTmatrix = as.data.frame.matrix(neighbourhoodbyCTmatrix)
  neighbourhoodbyCTmatrix = neighbourhoodbyCTmatrix[names(cellTypes),]
  
  return(neighbourhoodbyCTmatrix)
}


## ####################################################
#' This function creates a seurat object using a neighbourhood by cell type 
#' matrix
#' 
#' @param neighbourhoodByCT - a matrix of neighbourhoods by cell types
#' @param resolution - resolution for clustering (default 0.1)
#' @param npcs - number of pcs used for PCA, defaults to 10
#' @param n.neighbors - number of neighbors used by UMAP, defaults to 30
#' @param compute_UMAP - defaults to TRUE
#' @param transpose - defaults to FALSE
#' @return a seurat object based on a neighbourhood by cell type matrix,
#' containing clusters and UMAP.
#' @import Seurat
#' @export
#' @examples
#' smallNbhdObj = computeNeighbourhoodByCTSeurat(smallNbhdMatrix)
computeNeighbourhoodByCTSeurat= function(neighbourhoodByCT, resolution = 0.1, 
                                         npcs = 10, n.neighbors = 30L, 
                                         compute_UMAP = T, transpose = F){
    dataMatrix = t(neighbourhoodByCT)
    neighbourhoodSeurat = CreateSeuratObject(dataMatrix)
    neighbourhoodSeurat[['RNA']]$data = dataMatrix
    neighbourhoodSeurat = ScaleData(neighbourhoodSeurat)
    neighbourhoodSeurat = RunPCA(neighbourhoodSeurat, assay = "RNA", 
                                 features = rownames(neighbourhoodSeurat), 
                                 npcs = npcs)
    if (transpose){
        neighbourhoodSeurat = RunUMAP(neighbourhoodSeurat,assay='RNA',
                                      dims = 1:10, n.neighbors = n.neighbors)
  
    } else{
        neighbourhoodSeurat = RunUMAP(neighbourhoodSeurat,assay='RNA',
                                      features=rownames(neighbourhoodSeurat), 
                                      n.neighbors = n.neighbors)
    }
    if (transpose){
        neighbourhoodSeurat = FindNeighbors(neighbourhoodSeurat, dims = 1:npcs)
    } else{
        neighbourhoodSeurat = FindNeighbors(neighbourhoodSeurat, 
                                            features=rownames(neighbourhoodSeurat))
    }
    neighbourhoodSeurat = FindClusters(neighbourhoodSeurat,
                                       resolution=resolution)

    ## Rename seurat_clusters to neighbourhood_clusters
    idx = names(neighbourhoodSeurat@meta.data) == 'seurat_clusters'
    names(neighbourhoodSeurat@meta.data)[idx] = 'neighbourhood_clusters'
    
    return(neighbourhoodSeurat)
}

## ####################################################
#' This function adds a force directed graph embedding to a seurat object
#' 
#' @param seuratObj - a seurat object. 
#' @param graph - which graph to extract.  Defaults to
#' paste0(f@active.assay,'_snn')
#' @return a seurat object with a "graph" dimensionality reduction.
#' @import Seurat  
#' @import igraph
#' @export
#' @examples
#' objWithEmbedding = computeGraphEmbedding(smallNbhdObj)
computeGraphEmbedding = function(seuratObj, graph=defaultGraph(seuratObj)){
  graph = seuratObj@graphs[[graph]]
  igraphGraph = igraph::graph_from_adjacency_matrix(graph)
  graphLayout = igraph::layout_with_fr(igraphGraph)
  colnames(graphLayout) = c("graph_1","graph_2")
  rownames(graphLayout) = colnames(seuratObj)
  graphDimReduc = CreateDimReducObject(embeddings = graphLayout,   
                                       key = "graph",   assay = "RNA")
  seuratObj[["graph"]] = graphDimReduc
  return(seuratObj)
}



## #' ## ####################################################
## #' #' This function takes a nearest neighbour graph and a radius
## #' #' and finds the combinatorial ball around each cell
## #' #'
## #' #' @param NN - a nearest neighbour graph
## #' #' @param radius - combinatorial radius
## #' #' @return This returns a named list.  In each case, the name is the
## #' #'     cell at the center of the ball.  Each item in the list is a
## #' #'     named numeric.  The values give combinatorial distance from the
## #' #'     center and the names give the cells.
## #' #' @export
## #' #' @examples
## #' #' smallCombNbhds = findCombinatorialNeighbourhoods(smallDelaunayTriangulation,2)
## #' findCombinatorialNeighbourhoods = function(NN,radius)
## #' {
## #'     combinatorialNbhdsImpl = function(edges,radius)
## #'     {
## #'         ## ####################################################
## #'         ## The base case:
## #'         if(radius == 0)
## #'         {
## #'             ball = list()
## #'             vertices = names(edges)
## #'             for(v in vertices)
## #'             {
## #'                 a = 0
## #'                 names(a) = v
## #'                 ball[[v]] = a
## #'             }
## #'             return(ball)
## #'         }
## #'         ## ####################################################
## #'         ## Recurse:
## #'         else
## #'         {
## #'             smaller = combinatorialNbhdsImpl(edges,radius-1)
## #'             for(v in names(smaller))
## #'             {
## #'                 frontier = smaller[[v]][smaller[[v]]==radius-1]
## #'                 proposed = c()
## #'                 for(frontierV in names(frontier))
## #'                     proposed = c(proposed,edges[[frontierV]])
## #'                 proposed = unique(proposed)
## #'                 idx = proposed %in% names(smaller[[v]])
## #'                 proposed = proposed[!idx]
## #'                 theNew = rep(radius,length(proposed))
## #'                 names(theNew) = proposed
## #'                 smaller[[v]] = c(smaller[[v]],theNew)
## #'             }
## #'             return(smaller)
## #'         }
## #'     } ## End of impl
## #' 
## #'     ## ####################################################    
## #'     ## Get the edges:
## #'     vertices = unique(c(NN$nodeA,NN$nodeB))
## #'     edges = list()
## #'     for(v in vertices)
## #'     {
## #'         idx = NN$nodeA == v
## #'         neighboursB = unique(NN$nodeB[idx])
## #'         idx = NN$nodeB == v
## #'         neighboursA = unique(NN$nodeA[idx])
## #'         
## #'         edges[[v]] = unique(c(neighboursA,neighboursB))
## #'     }
## #' 
## #'     return(combinatorialNbhdsImpl(edges,radius))
## #' }
## #' 
## #' 
## #' ## ####################################################
## #' #' This function reduces combinatorial balls to an
## #' #' extended nearest neighbour graph
## #' #'
## #' #' @param balls - A named list as produced by the function
## #' #'     findCombinatorialNeighbourhoods
## #' #' @return This returns a nearest neighbour data frame where the cells
## #' #'     in each combinatorial ball are now considered the center's
## #' #'     neighbours.
## #' #' @export
## #' #' @examples
## #' #' reduced = reduceCombinatorialBalls(smallCombNbhds)
## #' reduceCombinatorialBalls = function(balls)
## #' {
## #'     nodeA = c()
## #'     nodeB = c()
## #'     for(cell in names(balls))
## #'     {
## #'         newNeighbours = names(balls[[cell]])
## #'         N = length(newNeighbours)
## #'         nodeA = c(nodeA,rep(cell,N))
## #'         nodeB = c(nodeB,newNeighbours)
## #'     }
## #'     NN = data.frame(nodeA,nodeB)
## #'     idx = NN$nodeA == NN$nodeB
## #'     NN = NN[!idx,]
## #' 
## #'     return(NN)
## #' }

## ####################################################
#' This function takes a nearest neighbour graph and a radius
#' and calculates nth neighbour graphs where max(n) == radius
#'
#' @param spatialGraph - a nearest neighbour graph
#' @param n - the maximum degree to calculate a neighbour graph with edges 
#' connecting vertices of degree n for.
#' @return A named list of neighbour graphs, where each graph contains edges 
#' connecting vertices of degree n. Each graph is named according to degree n.
#' @import data.table
#' @export
#' @examples
#' smallNbhds = expandNeighbourhoods(smallDelaunayTriangulation,4)
expandNeighbourhoods  = function(spatialGraph, n){
  spatialGraph = data.table(spatialGraph)
  
  spatialGraphR = spatialGraph
  names(spatialGraphR) = c("nodeB","nodeA")
  spatialGraph = rbind(spatialGraph,spatialGraphR[,c(2,1)])
  neighbours = list()
  neighbours[[1]] = spatialGraph
  
  for (i in (2:n)){
    writeLines(paste('radius',i))
    graph = merge(neighbours[[i-1]], neighbours[[1]], by.x = "nodeB", 
                  by.y = "nodeA", allow.cartesian = T)
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
#' This function takes expanded neighbourhoods and collapses them to a nearest 
#' neighbourhood graph where all neighbours of degree <= n in the original graph 
#' are considered first neighbours.
#'
#' @param expandedNeighbours - the results of expandNeighbourhoods()
#' @param n - the maximum degree to connect neighbours. Defaults to the maximum 
#' degree neighbourhoods were expanded to in the results of expandNeighbourhoods().
#' @return a graph in neighbour format, i.e., a data frame with
#'     columns nodeA and nodeB, where nodes that were originally of degree <= n 
#'     are connected.
#' @export
collapseNeighbourhoods = function(expandedNeighbours, n = length(expandedNeighbours)){
  
  collapsedGraph = as.data.frame(do.call(rbind,expandedNeighbours[1:n]))

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
#' smallCellTypesPerCellType = computeCellTypesPerCellTypeMatrix(smallNbhdMatrix,
#'                                                      smallXenium$seurat_clusters)
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
#' @param colors - a named vector of colors used to color the
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
#' coordinates are found in G$coords.  If colors were supplied
#' these are found in V(G)$color.  Edge weights and widths are
#' found in E(G)$weight and E(G)$width.
#' @export
#' @examples
#' idx = colSums(smallCellTypesPerCellTypeMatrix) > 0
#' M = smallCellTypesPerCellTypeMatrix[,idx]
#' G = cellTypesPerCellTypeGraphFromMatrix(M,plotGraph=FALSE)
cellTypesPerCellTypeGraphFromMatrix = function(M,
                                               colors=NULL,
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

    if(! is.null(colors))    
        V(G)$color = colors[names(V(G))]
    G$coords = layout_with_fr(G)
    E(G)$width = edgeWeighting * E(G)$weight

    if(plotGraph)
    {
        if(! is.null(colors))
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
#' @param colors - a named vector of colors used to color the
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
#' coordinates are found in G$coords.  If colors were supplied
#' these are found in V(G)$color.  Edge weights and widths are
#' found in E(G)$weight and E(G)$width. 
#' @export
cellTypesPerCellTypeGraph = function(nbhdByCellType,
                                     clusters,
                                     colors=NULL,
                                     selfEdges=FALSE,
                                     minWeight=0,
                                     edgeWeighting=20,
                                     edgeCurved=0.2,
                                     arrowSize=4,
                                     arrowWidth=4,
                                     plotGraph=TRUE)
{
    M = cellTypesPerCellTypeMatrix(nbhdByCellType,clusters)
    G = cellTypesPerCellTypeGraphFromMatrix(M,
                                            colors=colors,
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
#' @param verbose - whether to print trace.  Defaults to TRUE
#' @return A square matrix containing upper tail p values describing whether two 
#' cell types are more frequently found together than expected by chance.
#' @import abind
#' @export
## #' @examples
computeNeighbourEnrichment = function(spatialGraph, cellTypes, nSim = 1000,
                                      verbose=TRUE){
  results = list()
  spatialGraphOrig = spatialGraph
  neighbourhoodByCT = computeNeighbourhoodByCTMatrix(spatialGraphOrig,cellTypes) 
  cellTypeMatrix = computeCellTypesPerCellTypeMatrix(neighbourhoodByCT, cellTypes) 
  for (i in 1:nSim){ 
    spatialGraph = spatialGraphOrig
    spatialGraphBA = spatialGraph[,c(2,1)]
    names(spatialGraphBA) = c("nodeA","nodeB")
    spatialGraph = rbind(spatialGraph,spatialGraphBA)
    spatialGraph = unique(spatialGraph)
    spatialGraph[,2] = sample(spatialGraph[,2])
    spatialGraph[,2] = cellTypes[spatialGraph[,2]]
    neighbourhoodbyCTmatrix = table(spatialGraph[,1],spatialGraph[,2])
    neighbourhoodbyCTmatrix = as.data.frame.matrix(neighbourhoodbyCTmatrix)
    neighbourhoodbyCTmatrix = neighbourhoodbyCTmatrix[names(cellTypes),]
    results[[i]] = computeCellTypesPerCellTypeMatrix(neighbourhoodbyCTmatrix,cellTypes)
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
#' @param obj - a Seurat object
#' @param species - either 'human' or 'mouse'
#' @param lrn - a ligand-receptor network, i.e., a
#' data frame with columns from and to.  By default, it
#' retrieves the nichenetr ligand receptor network
#' @return This returns a data frame with columns ligand and
#' receptor
#' @export
getLigandReceptorPairsInPanel = function(obj,species,
                                         lrn = getLigandReceptorNetwork(species))
{
    stopifnot(species %in% c('mouse','human'))

    ## The panel
    panel = rownames(obj)

    ## Ligands and receptors:
    pairs = paste(lrn$from,lrn$to,sep='_')
    ligands = unique(lrn$from)
    receptors = unique(lrn$to)
    ligandsFound = intersect(ligands,panel)
    receptorsFound = intersect(receptors,panel)

    panelPairs = c()
    for(a in ligandsFound)
        for(b in receptorsFound)
            panelPairs = c(panelPairs,paste(a,b,sep='_'))
    idx = panelPairs %in% pairs
    pairsFound = panelPairs[idx]

    a = str_split(pairsFound,'_')
    pairsFoundDF = data.frame(ligand=unlist(lapply(a,function(x) return(x[1]))),
                              receptor=unlist(lapply(a,function(x) return(x[2]))))

    return(pairsFoundDF)
}

## #' ## ####################################################
## #' #' This function takes a Seurat object, a set of ligand receptor
## #' #' pairs and a set of edges denoting neighbouring cells and
## #' #' annotates these with the ligand receptor interactions taking
## #' #' place on those edges in each direction.
## #' #'
## #' #' @param obj - a Seurat object
## #' #' @param pairDF - a data frame giving the ligand-receptor pairs
## #' #' @param delaunayNeighbours - a data frame of neighbouring
## #' #' cell pairs.  Note that each row is a directed edge (A,B) so
## #' #' that this data frame should have both the edge (A,B) and the
## #' #' edge (B,A)
## #' #' @return This returns a data frame whose first two columns give
## #' #' the neighbouring cells.  Each of the remaining columns is a logical
## #' #' corresponding to a ligand-receptor pair telling whether the ligand
## #' #' is expressed in the first cell and the receptor is expressed in the
## #' #' second cell.
## #' #' @export
## #' getInteractionsOnEdges = function(obj,pairDF,delaunayNeighbours)
## #' {
## #'     ## Discretize expression:
## #'     M = FetchData(obj,rownames(obj),layer='count')
## #'     M = data.matrix(t(M))
## #'     cutoff = 0
## #'     M = M > cutoff
## #' 
## #'     ## Find the interactions on the edges:
## #'     edges = delaunayNeighbours
## #' 
## #'     for(i in 1:nrow(pairDF))
## #'     {
## #'         tag = paste(pairDF$ligand[i],pairDF$receptor[i],sep='_')
## #'         edges[,tag] = (M[pairDF$ligand[i],edges$nodeA] &
## #'                    M[pairDF$receptor[i],edges$nodeB])
## #'     }
## #' 
## #'     return(edges)
## #' }



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
#' @export
getInteractionsOnEdges = function(M,pairDF,spatialGraph)
{
  ## Find the interactions on the edges:
  edges = spatialGraph
  
  for(i in 1:nrow(pairDF))
  {
    tag = paste(pairDF$ligand[i],pairDF$receptor[i],sep='_')
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
#' @param obj - a Seurat object
#' @param cutoff - a cutoff for binarisation. Defaults to 0.
#' @return A binarised expression matrix where rows are genes and columns are 
#' cells.
getBinarisedMatrix = function(obj, cutoff = 0){
  M = FetchData(obj,rownames(obj),layer='count')
  M = data.matrix(t(M))
  cutoff = 0
  M = M > cutoff
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
countInteractionsPerCell = function(edges,sourceOrTarget)
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
#' by countInteractionsPerCell(), the underlying Seurat object
#' and the neighbourhood Seurat object and annotates the counts
#' with the cell type and the neighbourhood type corresponding
#' to the cells of the interaction counts.
#'
#' @param interactionCounts - as found by countInteractionsPerCell()
#' @param obj - a Seurat object
#' @param nbhdObj - a neighbourhood x cell type Seurat object
#' @return This returns the interaction counts annotated with the
#' cell type and neighbourhood type of each cell.
#' @export
annotateInteractionCounts = function(interactionCounts,obj,nbhdObj)
{
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
#' @return This returns a data frame whose first two columns give
#' the neighbouring cells.  Each of the remaining columns is a logical
#' corresponding to a ligand-receptor pair telling whether the ligand
#' is expressed in the first cell and the receptor is expressed in the
#' second cell.
#' @export

performInteractionAnalysis = function(obj, spatialGraph, species, clusters,
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
    paste0(interactionsOnEdges$clusterA,  "_", interactionsOnEdges$clusterB) 
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
              "simResults" = simResults,
              "pValues" = pValues))
}

    
## ####################################################
#' This function takes interactionResults and plots a heatmap of -log10(pvalues).
#'
#' @param interactionResults - as returned by performInteractionAnalysis()
#' @param clusters - named vector of cell types where names are each cell and
#' clusters are a factor
#' @param colors - a named list of colors where names are clusters. If not 
#' specified the default pheatmap color scheme will be used.
#' @param  pValCutoffClusterPair - a cutoff for showing interactions between two
#' clusters. A cluster pair must have at least one ligand-receptor interaction
#' pvalue <  pValCutoffClusterPair. Defaults to 0.05.
#' @param  pValCutoffLigRec  - a cutoff for showing interactions between a 
#' ligand and receptor. At least one cluster pair must have 
#' pvalue <  pValCutoffLigRec for ligand-receptor pair. Defaults to 0.05.
#' @param  labelClusterPairs - show labels for cluster pairs. Defaults to TRUE.
#' @import Rfast
#' @import pheatmap
#' @export
makeInteractionHeatmap = function(interactionResults,
                                  clusters,
                                  colors = c(),
                                  pValCutoffClusterPair = 0.05, 
                                  pValCutoffLigRec = 0.05,
                                  labelClusterPairs = T)
{
  selectedPValues = interactionResults$pValues[rowMins(interactionResults$pValues, value = T) < pValCutoffClusterPair,
                                             colMins(interactionResults$pValues, value = T) < pValCutoffLigRec]
  negLog10PValues = -log10(selectedPValues)
  rowAnno = str_split_fixed(rownames(selectedPValues), pattern = "_", 2)
  rowAnno = as.data.frame(rowAnno)
  names(rowAnno) = c("sender","receiver")
  rowAnno$sender = factor(rowAnno$sender, levels = levels(clusters))
  rowAnno$receiver = factor(rowAnno$receiver, levels = levels(clusters))
  rownames(rowAnno) = rownames(selectedPValues)
  if (length(colors) > 0){
  pheatmap(negLog10PValues, annotation_row = rowAnno, annotation_colors = list("sender" = colors, 
                                                                    "receiver" = colors),
           show_rownames = labelClusterPairs)
  } else{
    pheatmap(negLog10PValues, annotation_row = rowAnno, show_rownames = labelClusterPairs)
  }
}


## ####################################################
#' This function takes interactionResults and plots a heatmap of the total 
#' number of ligand receptor interactions between clusters.
#'
#' @param interactionResults - as returned by performInteractionAnalysis()
#' @param clusters - named vector of cell types where names are each cell and
#' clusters are a factor
#' @param type - "total" or "mean" to plot raw total interactions or mean interactions per edge.
#' @param  logScale - plot heatmap using log scale (defaults to T)
#' @import pheatmap
#' @export
makeSummedInteractionHeatmap = function(interactionResults, clusters, type, logScale = T){ 
  if (type == "total"){
    interactionsByCluster = interactionResults$totalInteractionsByCluster
  } 
  if (type == "mean"){
    interactionsByCluster = interactionResults$meanInteractionsByCluster
  } 
  summedInteractionsByCluster = rowSums(interactionsByCluster[,3:ncol(interactionsByCluster)])
  pair = str_split_fixed(names(summedInteractionsByCluster), pattern = "_", 2)
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
  } else{
    pheatmap((summedInteractionsByClusterMatrix))
  }
}

## ####################################################
#' This function takes interactionResults and creates a seurat object where 
#' each point represents an edge between cells, and spatial coordinates are the 
#' centroids of edges between cells. The "expression matrix" is the 
#' binarised presence/absence of an interaction (ligand receptor pair) on an edge. 
#'
#' @param interactionResults - as returned by performInteractionAnalysis()
#' @param centroids - a dataframe containing centroids 
#' where rownames are cellnames and the first two columns
#' contain x and y coordinates respectively.
#' @param npcs - number of pcs used for PCA, defaults to 10
#' @import Seurat
#' @return This returns a seurat object where 
#' each point represents an edge between cells, and spatial coordinates are the 
#' centroids of edges between cells. The "expression matrix" is the 
#' binarised presence/absence of an interaction (ligand receptor pair) on an edge. 
#' @export

computeEdgeSeurat = function(interactionResults, centroids, npcs = 10){
  interactionsOnEdges = interactionResults$interactionsOnEdges
  rownames(interactionsOnEdges) = paste0(interactionsOnEdges$nodeA, "_", interactionsOnEdges$nodeB)
  interactionsOnEdgesMat = as.matrix(interactionsOnEdges[,5:ncol(interactionsOnEdges)])
  interactionsOnEdgesMat= 1 * interactionsOnEdgesMat
  edgeSeurat = CreateSeuratObject(t(interactionsOnEdgesMat), meta = interactionsOnEdges[,1:5])
  edgeCoords = cbind(centroids[interactionsOnEdges$nodeA, 1:2], centroids[interactionsOnEdges$nodeB, 1:2])
  
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
  return(edgeSeurat)
}



## ####################################################
#' nbhdsAsEdgesToNbhdsAsList
#'
#' This function takes a set of neighbourhoods given
#' by edges and turns it into a named list giving the
#' memberships of each neighbourhood
#'
#' @param cells - The cells whose neighbourhoods to
#' extract.  Defaults to neighbourhoods$nodeA
#' @param neighbourhoods - neighbourhoods given as a
#' data frame with columns nodeA and nodeB, for example
#' the output of collapseNeighbourhoods
#' @return a named list with memberships of the neighbourhoods
#' of cells
#' @export
nbhdsAsEdgesToNbhdsAsList = function(cells=unique(neighbourhoods$nodeA),
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
#' @param f - a Seurat object with layer counts
#' @param neighbourhoods - Neighbourhoods as given by a
#' collapsed expanded edge graph, as produced by
#' collapseNeighbourhoods. In particular, each cell should
#' appear as nodeA.
#' @import Seurat
#' @return a Seurat object giving total gene expression
#' in each neighbourhood.
#' @export
aggregateGeneExpression = function(f,neighbourhoods)
{
    cells = colnames(f)
    nbhds = unique(neighbourhoods$nodeA)

    ## Sanity check:
    stopifnot(identical(cells[order(cells)],
                        nbhds[order(nbhds)]))

    genes = rownames(f)
    counts = FetchData(f,genes,layer='counts')
    counts = t(counts)
    counts = data.matrix(counts)

    nbhdList =nbhdsAsEdgesToNbhdsAsList(cells,
                                        neighbourhoods)
    
    C = aggregateFeature(counts, nbhdList, rowsums)

    nbhdObj = CreateSeuratObject(counts=C)
    nbhdObj = NormalizeData(nbhdObj)
    nbhdObj = ScaleData(nbhdObj)
    nbhdObj = FindVariableFeatures(nbhdObj)
    nbhdObj = RunPCA(nbhdObj)
    nbhdObj = RunUMAP(nbhdObj,dims=1:20)
    nbhdObj = FindNeighbors(nbhdObj)
    nbhdObj = FindClusters(nbhdObj)
    
    return(nbhdObj)
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
aggregateFeature = function(M, nbhdList, aggregateFunction)
{
    cells = colnames(M)
    res = lapply(cells, function(x, M, nbhdList) aggregateFunction(M[,nbhdList[[x]], drop = F]), 
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
#' @param aggregateFunction - a function to aggregate expression (e.g. rowSums,
#' rowMeans)

#' @return a matrix giving aggregated gene expression for a cell's neighbourhood.
#' @export
computeMoransI = function(M,nbhdList){
  aggrM = aggregateFeature(M,nbhdList, rowMeans)
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
#' @return a dataframe containing Moran's I and p values for each feature.
#' @export

runMoransI = function(obj, spatialGraph, assay = "RNA", layer = "data",
                      nSim = 100, verbose = T){
  spatialGraph = symmetriseNN(spatialGraph)
  
  M = as.matrix(LayerData(obj, assay = assay, layer = layer))
  nbhdList = nbhdsAsEdgesToNbhdsAsList(colnames(M),
                                       neighbourhoods)
  
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
  results = results[order(moransI, decreasing = T),]
  results = data.frame(results)
  return(results)
}
