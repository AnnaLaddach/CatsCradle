## ####################################################
#' This function computes a spatial graph where 
#' neighbors are identified based on Delaunay triangulation.
#' 
#' @param centroids - a dataframe containing centroids 
#' where rownames are cellnames and columns contain x 
#' and y coordinates respectively.
#' @return a graph in neighbour format.
#' @import interp
#' @export
#' @examples
computeNeighboursDelaunay = function(centroids){
  
  #compute delaunay triangulation
  delaunay = tri.mesh(centroids[,1], centroids[,2])
  
  #extract triangles from results dataframe
  triangles = delaunay$trlist
  
  results = matrix(nrow = 0, ncol = 2)
  cellNames = rownames(centroids)
  
  #convert triangle indices to cell names
  triangles[,1] = cellNames[as.numeric(triangles[,1])]
  triangles[,2] = cellNames[as.numeric(triangles[,2])]
  triangles[,3] = cellNames[as.numeric(triangles[,3])]
  
  #create neighbour list
  for (i in 1:nrow(triangles)){
    results = rbind(results,sort(c(triangles[i,1],triangles[i,2])))
    results = rbind(results,sort(c(triangles[i,1],triangles[i,3])))
    results = rbind(results,sort(c(triangles[i,2],triangles[i,3])))
  }
  
  #remove duplicate edges
  results = unique(results)
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
#' @return a graph in neighbour format.
#' @import rdist
#' @export
#' @examples
computeNeighboursEuclidean = function(coords, threshold){
  colnames(centroids) = c("x","y")
  results = matrix(nrow = 0, ncol = 2)
  
  for (i in 1:nrow(coords)){
    cellX = centroids[i,1]
    cellY = centroids[i,2]
    cellName = rownames(centroids)[i]
    
    #only calculate distances for cells that could be within distance threshold.
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
  #remove duplicate edges
  results = unique(results)
  
  return(results)
}



## ####################################################
#' This function extracts neighbourhoods from a spatial graph in neighbour list
#' format
#' 
#' @param  spatialGraph - a spatial graph in neighbour list format.
#' @return a named list of neighbourhoods where a neighborhood is a set of a 
#' cell's nearest neighbours. 
#' @export
#' @examples
#' 
#' 
computeNeighbourhoods = function(spatialGraph, cellNames, addSelf = F){
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
#' @examples

computeNeighbourhoodByCTMatrix = function(neighbourhoods, cellTypes){
  
  cellNames = names(cellTypes)
  neighbourhoodByCT = matrix(nrow = 0, ncol = length(levels(cellTypes)))
  
  for (neighbourhood in neighbourhoods){
    neighbourhoodByCT = rbind(neighbourhoodByCT,table(cellTypes[neighbourhood]))
  }
  
  return(neighbourhoodByCT)
}
  

#' ## ####################################################
#' #' This function computes a matrix where cells are rows and 
#' #' cell types are columns. The values in the matrix indicate the  
#' #' number of spatial neighbours of a given cell type a particular cell has. 
#' #' 
#' #' @param  spatialGraph - a spatial graph in neighbour list format.
#' #' @param cellTypes - named vector of cell types where names are each cell and
#' #' cell types are a factor.
#' #' @return a matrix of cells by cell types.
#' #' @export
#' #' @examples
#' 
#' computeCellByCTMatrix = function(spatialGraph, cellTypes){
#'   
#'   cellNames = names(cellTypes)
#'   cellByCt = matrix(nrow = 0, ncol = length(levels(cellTypes)))
#' 
#'   for (cell in cellNames){
#'     neighbours = 
#'       c(spatialGraph[,2][cell == spatialGraph[,1]],
#'         spatialGraph[,1][cell == spatialGraph[,2]])
#'     cellByCt = rbind(cellByCt,table(xenium.obj@active.ident[neighbours]))
#'   }
#'   
#'   return(cellByCt)
#' }


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
#' @examples

computeNeighbourhoodByCTSeurat= function(neighbourhoodByCT, resolution = 0.1, 
                                         npcs = 10, n.neighbors = 30L, 
                                         compute_UMAP = T, transpose = F){
  neighbourhoodSeurat = CreateSeuratObject(t(neighbourhoodByCT))
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
#' @examples

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