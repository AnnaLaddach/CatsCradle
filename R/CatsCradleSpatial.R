## ####################################################
#' This function computes a spatial graph where 
#' neighbors are identified based on Delaunay triangulation.
#' 
#' @param centroids - a dataframe containing centroids 
#' where rownames are cellnames and columns contain x 
#' and y coordinates respectively.
#' @return a graph in neighbour format, i.e., a data frame with
#'     columns nodeA and nodeB.
#' @import interp
#' @export
## #' @examples
computeNeighboursDelaunay = function(centroids){
    
    ##compute delaunay triangulation
    delaunay = tri.mesh(centroids[,1], centroids[,2])
    
    ##extract triangles from results dataframe
    triangles = delaunay$trlist
    
    results = matrix(nrow = 0, ncol = 2)
    cellNames = rownames(centroids)
    
    ##convert triangle indices to cell names
    triangles[,1] = cellNames[as.numeric(triangles[,1])]
    triangles[,2] = cellNames[as.numeric(triangles[,2])]
    triangles[,3] = cellNames[as.numeric(triangles[,3])]
    
    ##create neighbour list
    for (i in 1:nrow(triangles)){
        results = rbind(results,sort(c(triangles[i,1],triangles[i,2])))
        results = rbind(results,sort(c(triangles[i,1],triangles[i,3])))
        results = rbind(results,sort(c(triangles[i,2],triangles[i,3])))
    }
    
    ##remove duplicate edges
    results = unique(results)
    
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
#' @export
## #' @examples
computeNeighboursEuclidean = function(coords, threshold){
    colnames(centroids) = c("x","y")
    results = matrix(nrow = 0, ncol = 2)
    
    for (i in 1:nrow(coords)){
        cellX = centroids[i,1]
        cellY = centroids[i,2]
        cellName = rownames(centroids)[i]
        
        ##only calculate distances for cells that could be within distance threshold.
        possibleNeighbours = centroids[(centroids[,"x"] < (cellX + threshold)) &
                                       (centroids[,"x"] > (cellX - threshold)) &
                                       (centroids[,"y"] < (cellY + threshold)) &
                                       (centroids[,"y"] > (cellY - threshold)),,
                                       drop = FALSE]
        
        neighbours = rownames(possibleNeighbours)[cdist(coords[i,,drop=FALSE]
                                                       ,possibleNeighbours,)[1,] < threshold]
        
        for (neighbour in neighbours){
            if (cellName != neighbour){
                results = rbind(results,sort(c(cellName,neighbour)))
            }
        }
    }
    ##remove duplicate edges
    results = unique(results)
    
    ## Convert to data.frame and name as nodeA, nodeB:
    results = as.data.frame(results)
    names(results) = c('nodeA','nodeB')
    
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
#' This function extracts neighbourhoods from a spatial graph in neighbour list
#' format
#' 
#' @param  spatialGraph - a spatial graph in neighbour list format.
#' @return a named list of neighbourhoods where a neighborhood is a set of a 
#' cell's nearest neighbours. 
#' @export
## #' @examples
#' 
#' 
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
#' @param neighbourhoods- a named list of neighbourhoods where a neighbourhood 
#' is a set of cells. 
#' @param cellTypes - named vector of cell types where names are each cell and
#' cell types are a factor
#' @return a matrix of neighbourhoods by cell types
#' @export
## #' @examples

computeNeighbourhoodByCTMatrix = function(neighbourhoods, cellTypes){
  
    cellNames = names(cellTypes)
    neighbourhoodByCT = matrix(nrow = 0, ncol = length(levels(cellTypes)))
    
    for (neighbourhood in neighbourhoods){
        neighbourhoodByCT = rbind(neighbourhoodByCT,table(cellTypes[neighbourhood]))
    }
    rownames(neighbourhoodByCT) = names(neighbourhoods)
    
    return(neighbourhoodByCT)
}
  

## #' ## ####################################################
## #' #' This function computes a matrix where cells are rows and 
## #' #' cell types are columns. The values in the matrix indicate the  
## #' #' number of spatial neighbours of a given cell type a particular cell has. 
## #' #' 
## #' #' @param  spatialGraph - a spatial graph in neighbour list format.
## #' #' @param cellTypes - named vector of cell types where names are each cell and
## #' #' cell types are a factor.
## #' #' @return a matrix of cells by cell types.
## #' #' @export
## #' #' @examples
## #' 
## #' computeCellByCTMatrix = function(spatialGraph, cellTypes){
## #'   
## #'   cellNames = names(cellTypes)
## #'   cellByCt = matrix(nrow = 0, ncol = length(levels(cellTypes)))
## #' 
## #'   for (cell in cellNames){
## #'     neighbours = 
## #'       c(spatialGraph[,2][cell == spatialGraph[,1]],
## #'         spatialGraph[,1][cell == spatialGraph[,2]])
## #'     cellByCt = rbind(cellByCt,table(xenium.obj@active.ident[neighbours]))
## #'   }
## #'   
## #'   return(cellByCt)
## #' }


## ####################################################
#' This function creates a seurat object using a neighbourhood by cell type 
#' matrix
#' 
#' @param neighbourhoodByCT - a matrix of neighbourhoods by cell types
#' @param resolution - resolution for clustering (default 0.1)
#' @return a seurat object based on a neighbourhood by cell type matrix,
#' containing clusters and UMAP.
#' @import Seurat
#' @export
## #' @examples

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
                                      n.neighbors = n.neighbors)
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
## #' @examples

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



## ####################################################
#' This function takes a nearest neighbour graph and a radius
#' and finds the combinatorial ball around each cell
#'
#' @param NN - a nearest neighbour graph
#' @param radius - combinatorial radius
#' @return This returns a named list.  In each case, the name is the
#'     cell at the center of the ball.  Each item in the list is a
#'     named numeric.  The values give combinatorial distance from the
#'     center and the names give the cells.
#' @export
findCombinatorialNeighbourhoods = function(NN,radius)
{
    combinatorialNbhdsImpl = function(edges,radius)
    {
        ## ####################################################
        ## The base case:
        if(radius == 0)
        {
            ball = list()
            vertices = names(edges)
            for(v in vertices)
            {
                a = 0
                names(a) = v
                ball[[v]] = a
            }
            return(ball)
        }
        ## ####################################################
        ## Recurse:
        else
        {
            smaller = combinatorialNbhdsImpl(edges,radius-1)
            for(v in names(smaller))
            {
                frontier = smaller[[v]][smaller[[v]]==radius-1]
                proposed = c()
                for(frontierV in names(frontier))
                    proposed = c(proposed,edges[[frontierV]])
                proposed = unique(proposed)
                idx = proposed %in% names(smaller[[v]])
                proposed = proposed[!idx]
                theNew = rep(radius,length(proposed))
                names(theNew) = proposed
                smaller[[v]] = c(smaller[[v]],theNew)
            }
            return(smaller)
        }
    } ## End of impl

    ## ####################################################    
    ## Get the edges:
    vertices = unique(c(NN$nodeA,NN$nodeB))
    edges = list()
    for(v in vertices)
    {
        idx = NN$nodeA == v
        neighboursB = unique(NN$nodeB[idx])
        idx = NN$nodeB == v
        neighboursA = unique(NN$nodeA[idx])
        
        edges[[v]] = unique(c(neighboursA,neighboursB))
    }

    return(combinatorialNbhdsImpl(edges,radius))
}


## ####################################################
#' This function reduces combinatorial balls to an
#' extended nearest neighbour graph
#'
#' @param balls - A named list as produced by the function
#'     findCombinatorialNeighbourhoods
#' @return This returns a nearest neighbour data frame where the cells
#'     in each combinatorial ball are now considered the center's
#'     neighbours.
#' @export
reduceCombinatorialBalls = function(balls)
{
    nodeA = c()
    nodeB = c()
    for(cell in names(balls))
    {
        newNeighbours = names(balls[[cell]])
        N = length(newNeighbours)
        nodeA = c(nodeA,rep(cell,N))
        nodeB = c(nodeB,newNeighbours)
    }
    NN = data.frame(nodeA,nodeB)

    return(NN)
}

## ####################################################
#' For each cell type, this function looks at the neighbourhoods
#' around cells of that type and discovers the fractions of those
#' cells of each type.
#'
#' @param nbhdByCellType - A matrix whose rows are neighbourhoods
#' each denoted by the cell at their center, whose columns are
#' cell types, and whose entries are counts.
#' @param seurat_clusters - a named vector whose names are the cells
#' and whose entries are their seurat_clusters.
#' @return A square matrix whose rownames and colnames are the
#' seurat_clusters as character strings.  Each row corresponds
#' to neighbourhoods around all cells of that type and the entries
#' give the fractions of those neighbourhoods occupied by cells
#' of each type.
#' @export
cellTypesPerCellTypeMatrix = function(nbhdByCellType,seurat_clusters)
{
    clusters = unique(seurat_clusters)
    clusters = clusters[order(clusters)]
    N = length(clusters)
    M = matrix(0,nrow=N,ncol=N)
    Clusters = as.character(clusters)
    rownames(M) = Clusters
    colnames(M) = Clusters

    for(cell in rownames(nbhdByCellType))
    {
        type = as.character(seurat_clusters[cell])
        M[type,] = M[type,] + nbhdByCellType[cell,]
    }

    rowTotals = rowSums(M)
    MM = M
    for(i in 1:nrow(M))
        MM[i,] = MM[i,]  / rowTotals[i]

    return(MM)
}

## ####################################################
#' This function converts a matrix as found by
#' cellTypesPerCellTypeMatrix into a directed igraph
#' whose vertices correspond to seurat_clusters and whose
#' edge correspond to occupancy fraction.
#'
#' @param M - a matrix as found by cellTypesPerCellTypeMatrix
#' @param colors - a named vector of colors used to color the
#' vertices of the graph.  The names are the seurat_clusters
#' as character strings.
#' @param selfEdges - a logical which determines whether to include
#' self edges.  Defaults to FALSE
#' @param minWeight - Allows one to exclude edges of low weight.
#' Defaults to 0, thus including all edges.
#' @param edgeWeighting - a parameter used to thicken the edges
#' in the display.  Defaults to 20.
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
cellTypesPerCellTypeGraphFromMatrix = function(M,
                                               colors=NULL,
                                               selfEdges=FALSE,
                                               minWeight=0,
                                               edgeWeighting=20,
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
                 edge.width=E(G)$width)
        else
            plot(G,
                 layout=G$coords,
                 edge.width=E(G)$width)
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
#' @param seurat_clusters - a named vector whose names are the cells
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
                                     seurat_clusters,
                                     colors=NULL,
                                     selfEdges=FALSE,
                                     minWeight=0,
                                     edgeWeighting=20,
                                     plotGraph=TRUE)
{
    M = cellTypesPerCellTypeMatrix(nbhdByCellType,seurat_clusters)
    G = cellTypesPerCellTypeGraphFromMatrix(M,
                                            colors=colors,
                                            selfEdges=selfEdges,
                                            minWeight=minWeight,
                                            edgeWeighting=edgeWeighting,
                                            plotGraph=plotGraph)

    return(G)
}
