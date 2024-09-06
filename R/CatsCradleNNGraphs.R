
## ####################################################
## These functions focus on manipulation of nearest
## neighbour graphs.
## ####################################################


## ####################################################
#' This function extracts a shared nearest neighbor network
#' from a Seurat object
#'
#' @param f - a Seurat object or SingleCellExperiment to
#' be converted to a Seurat object
#' @param graph - which graph to extract.  Defaults to
#' paste0(f@active.assay,'_snn')
#' @return - This returns dataframe of neighbors:
#' nodeA - node names for node A 
#' nodeB - node names for node B
#' weight - edge weight
#' @export
#' @examples
#' STranspose = make.getExample()('STranspose',toy=TRUE)
#' NN = getNearestNeighbourLists(STranspose)
getNearestNeighbourLists = function(f, graph=defaultGraph(f)){

    f = acceptor(f)
  
    ## convert to TsparseMatrix and extract relevant information
    graph = as(f@graphs[[graph]], "TsparseMatrix") 
    neighborListDf = data.frame("nodeA" = graph@Dimnames[[1]][graph@i+1],
                                "nodeB" =  graph@Dimnames[[2]][graph@j+1], 
                                "weight" = graph@x)
    
    ## remove self-loops
    neighborListDf = 
        neighborListDf[neighborListDf$nodeA != neighborListDf$nodeB,]
    return(neighborListDf)
}


## ####################################################
#' This function takes the data frame of neighbor genes
#' and reduces it so that each undirected edge is
#' represented by only one directed edge.  This ensures
#' that randomisation does not magically split undirected
#' edges into two edges.
#'
#' @param NN - a dataframe containing the neighborlist
#' @return - a neighborListDF with only one directed edge per
#' undirected edge.
#' @export
#' @examples
#' NN = make.getExample()('NN',toy=TRUE)
#' print(dim(NN))
#' NNN = desymmetriseNN(NN)
#' print(dim(NNN))
desymmetriseNN = function(NN)
{
    ## Use this to order the genes on each edge
    orderGenes = function(i)
    {
        if(order(c(NN$nodeA[i],NN$nodeB[i]))[1] == 2)
            return(TRUE)
        return(FALSE)
    }

    ## Order the genes on each edge:
    idx = unlist(lapply(seq_len(nrow(NN)),orderGenes))
    NN[idx,seq_len(2)] = NN[idx,seq(from=2,to=1)]

    ## Delete the duplicates:
    tag = paste(NN$nodeA,NN$nodeB,sep='___')
    NN = NN[!duplicated(tag),]

    return(NN)
}

## ####################################################
#' This function generates random indices for node B
#'
#' @param neighborListDf - a dataframe containing the neighborlist
#' @param n - the number of times to randomise indices
#' @param useWeights - whether to preserve edgeweights.
#' @return - a matrix with randomised indices for node B
#' @export
#' @examples
#' NN = make.getExample()('NN')
#' NN = desymmetriseNN(NN)
#' randomIndices = randomiseNodeIndices(NN,10,TRUE)
randomiseNodeIndices = function(neighborListDf, n = 100, useWeights = FALSE){
    NN = desymmetriseNN(neighborListDf)
    if(!identical(NN,neighborListDf))
        stop(paste0('randomiseNodeIndices is meant to be used',
                    'with a desymmetrised neighborList'))
  
  #determine number of edges and create empty matrix for randomised indices
  nEdges = nrow(neighborListDf)
  indices = seq_len(nEdges)
  randomIndices = matrix(, nrow = nEdges, ncol = n)
  
  #check if weights are to be used
  if (useWeights){
    
    #determine unique weights
    weights = unique(neighborListDf$weight)
    for (i in seq_len(n)){
      
      #randomise indices within each weight category
      for (weight in weights){
        selected = neighborListDf$weight == weight
        randomIndices[selected,i] =
          sample(indices[selected])
      }
    }
  }
  
  #otherwise ignore weights and randomise indices
  else {
    for (i in seq_len(n)){
      randomIndices[,i] = sample(indices)
    }
  }
  return(randomIndices)
}


## ####################################################
#' Tests whether a nearest neighbor graph is symmetric
#'
#' The nearest neighbor relationship is not inherently
#' symmetric.  This tests whether the nearest neighbor graph
#' retrieved from a Seurat object is.
#'
#' @param NN - a nearest neighbor graph.  This is in the form
#' of a data frame  as returned by getNearestNeighbourLists.
#' Its coloumns include nodeA and nodeB.
#' @return TRUE or FALSE
#' @export
#' @examples
#'  NN = make.getExample()('NN',toy=TRUE)
#' symmetryTest = symmetryCheckNN(NN)
symmetryCheckNN = function(NN)
{
    tagA = paste(NN$nodeA,NN$nodeB)
    tagA = tagA[order(tagA)]
    
    tagB = paste(NN$nodeB,NN$nodeA)
    tagB = tagB[order(tagB)]
    
    return(identical(tagA,tagB))
}


## ####################################################
#' This symmetrises a nearest neighbors graph.
#'
#' This first checks to see if the NN graph is symmetric
#' and if not symmetrises it.
#'
#' @param NN - a nearest neighbors graph as returned
#' by getNearestNeighbourLists
#'
#' @return a nearest neighbors graph
#' 
#' @export
#' @examples
#' NN = make.getExample()('NN',toy=TRUE)
#' NNStar = symmetriseNN(NN)
symmetriseNN = function(NN)
{
    ## Is it already symmetric?
    if(symmetryCheckNN(NN))
        return(NN)

    ## If not, symmetrise:
    NN2 = NN
    NN2[,seq_len(2)] = NN2[,seq(from=2,to=1)]
    NN = rbind(NN,NN2)

    tag = paste(NN$nodeA,NN$nodeB)
    idx = ! duplicated(tag)

    NN = NN[idx,]

    return(NN)
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
  toFlip = sample(seq_len(n), size = round(n/2))
  simGraph = spatialGraph[toFlip,c(2,1)]
  names(simGraph) =  c("nodeA","nodeB")
  simGraph = rbind(simGraph, spatialGraph[!(seq_len(n) %in% toFlip),])
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

