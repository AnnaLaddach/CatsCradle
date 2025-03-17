
## ####################################################
## These functions focus on creating and manipulating
## neighbourhoods.
## ####################################################

## ####################################################
#' This function computes a spatial graph where 
#' neighbors are identified based on Delaunay triangulation.
#' 
#' @param centroids - a dataframe containing centroids 
#' where rownames are cellnames and the first two columns
#' contain x and y coordinates respectively.
#' @return a graph in neighbour format, i.e., a data frame with
#'     columns nodeA and nodeB.
#' @importFrom geometry delaunayn
#' @importFrom Rfast rowSort
#' @export
#' @examples
#' centroids = make.getExample()('centroids')
#' delaunayNeighbours = computeNeighboursDelaunay(centroids)
computeNeighboursDelaunay = function(centroids){
    
    ## This is confounded when centroids has extra columns:
    centroids = centroids[,seq_len(2)]
    
    ##get cell names
    cellNames = rownames(centroids)
    
    ##compute delaunay triangulation
    triangles = delaunayn(centroids)
    
    ##extract neighbours
    results = rbind(triangles[,c(1,2)],
                    triangles[,c(1,3)],
                    triangles[,c(2,3)])
    
    ##sort rows
    results = rowSort(results)
    
    ## remove duplicates
    results = unique(results)
    
    ##convert indices to cell names
    results[,1] = cellNames[as.numeric(results[,1])]
    results[,2] = cellNames[as.numeric(results[,2])]
    
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
#' @importFrom rdist pdist
#' @export
#' @examples
#' centroids = make.getExample()('centroids')
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
    for (i in seq(from=0,to=9)){
        for (j in seq(from=0,to=9)){
            x1 = i * XStep + minX 
            x2 = x1 + 2*XStep
            y1 = j * YStep + minY 
            y2 = y1 + 2*YStep
            selected = centroids[(centroids[,"x"] >= x1) &
                                 (centroids[,"x"] <= x2) &
                                 (centroids[,"y"] >= y1) &
                                 (centroids[,"y"] <= y2), ]
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
#' neighbourhoodDiameter
#'
#' This function takes a list of neighbourhoods and and the
#' centroids of the cells and finds their diameters, i.e.,
#' for each neighbourhood, the maximum distance between.
#'
#' @param neighbourhoods - a list of neighbourhoods as
#' returned by nbhdsAsEdgesToNbhdsAsList
#' @param centroids - the centroids of the cells
#' @return a named numeric.  The names are the names
#' of the list neighbourhoods and the values are the
#' maximum distance within each neighbourhood
#' @export
#' @examples
#' getExample = make.getExample()
#' centroids = getExample('centroids')
#' delaunayNeighbours = getExample('delaunayNeighbours')
#' cells = unique(c(delaunayNeighbours[,'nodeA'],delaunayNeighbours[,'nodeB']))
#' nbhds = nbhdsAsEdgesToNbhdsAsList(cells,delaunayNeighbours)
#' diameters = neighbourhoodDiameter(nbhds[seq_len(100)],centroids)
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
        M[seq_len(N),seq_len(2)] = as.matrix(centroids[theseCells,seq_len(2)])
        D = distmat(M,M)
        
        diameters[nbhd] = max(D)
    }
    
    return(diameters)
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
#' @importFrom data.table data.table merge.data.table
#' @export
#' @examples
#' delaunayNeighbours = make.getExample()('delaunayNeighbours')
#' extendedNeighboursList = getExtendedNBHDs(delaunayNeighbours, 4)
getExtendedNBHDs = function(spatialGraph, n){
    spatialGraph = data.table(spatialGraph)
    
    spatialGraphR = spatialGraph
    names(spatialGraphR) = c("nodeB","nodeA")
    spatialGraph = rbind(spatialGraph,spatialGraphR[,c(2,1)])
    neighbours = list()
    neighbours[[1]] = spatialGraph
    
    for (i in seq(from=2,to=n)){
        writeLines(paste('radius',i))
        graph = merge.data.table(neighbours[[i-1]], neighbours[[1]],
                                 by.x = "nodeB", 
                                 by.y = "nodeA",
                                 allow.cartesian = TRUE)
        graph = graph[,c("nodeA","nodeB.y")]
        names(graph) = c("nodeA","nodeB")
        graph = unique(graph)
        orig = c(paste0(neighbours[[i-1]]$nodeB,
                        "_",
                        neighbours[[i-1]]$nodeA))
        if (i > 2){
            orig = c(orig,paste0(neighbours[[i-2]]$nodeB,
                                 "_",
                                 neighbours[[i-2]]$nodeA))
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
#' extendedNeighboursList = make.getExample()('extendedNeighboursList',toy=TRUE)
#' extendedNeighbours = collapseExtendedNBHDs(extendedNeighboursList, 4)
collapseExtendedNBHDs = function(extendedNeighboursList,
                                 n = length(extendedNeighboursList)){ 
    collapsedGraph = as.data.frame(do.call(rbind,
                                           extendedNeighboursList[seq_len(n)]))
    collapsedGraph = desymmetriseNN(collapsedGraph)
    
    return(collapsedGraph)
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
#' @param self - include cell in its neighbourhood, defaults to FALSE
#' @return a named list with memberships of the neighbourhoods
#' of cells
#' @export
#' @examples
#' delaunayNeighbours = make.getExample()('delaunayNeighbours')
#' cells = unique(c(delaunayNeighbours[,'nodeA'],delaunayNeighbours[,'nodeB']))
#' nbhdsList = nbhdsAsEdgesToNbhdsAsList(cells,delaunayNeighbours)
nbhdsAsEdgesToNbhdsAsList = function(cells,
                                     neighbourhoods, self = FALSE)
{
    nbhdList = list()
    if (self) {
      for(cell in cells)
      {
          nbhdList[[cell]] = c(cell,
                               neighbourhoods$nodeB[neighbourhoods$nodeA == cell])
      }
    } else {
      for(cell in cells)
      {
        nbhdList[[cell]] = c(neighbourhoods$nodeB[neighbourhoods$nodeA == cell])
      }
    }
    return(nbhdList)
}

