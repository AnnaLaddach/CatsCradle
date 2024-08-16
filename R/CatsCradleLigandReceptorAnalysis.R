
## ####################################################
## These functions focus of analysing ligand-receptor
## interactions between neighbouring cells is spatail
## data,
## ####################################################

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

    getExample = make.getExample()
    
    if(species == 'human')
    {
        humanLRN = getExample('humanLRN')
        return(humanLRN)
    }

    if(species == 'mouse')
    {
        mouseLRN = getExample('mouseLRN')
        return(mouseLRN)
    }
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
#' smallXenium = make.getExample()('smallXenium')
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
  
  for(i in seq_len(nrow(pairDF)))
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
  for (i in seq_len(nrow(M))){
    M[i,] = M[i,sample(n)]
  }
  return(M)
}

## ####################################################
#' This functions retrieves an expression matrix from a
#' seurat object or SingleCellExperiment and binarises it.
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

    edges = edges[,seq(from=3,to=ncol(edges))]
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
    pairs = names(interactionCounts)[seq(from=2,to=ncol(interactionCounts))]

    annotated[,pairs] = 0
    annotated[rownames(interactionCounts),pairs] =
        interactionCounts[rownames(interactionCounts),pairs]

    return(annotated)
}

## ####################################################
#' Given a seurat object, a spatial graph, clusters and
#' species this function identifies ligand-receptor 
#' interactions between neighbouring cells, identifies
#' ligand-receptor interactions within and between clusters
#' and calculates whether these are observed more frequently 
#' than expected by chance.
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
#' getExample = make.getExample()
#' smallXenium = getExample('smallXenium')
#' delaunayNeighbours = getExample('delaunayNeighbours')
#' clusters = getExample('clusters')
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
  names(interactionsOnEdges)[seq_len(2)] = c("clusterA","clusterB") 
  
  #get sum of interactions within and between clusters
  pair = 
    paste0(interactionsOnEdges$clusterA,  "-", interactionsOnEdges$clusterB) 
  totalInteractionsByCluster = 
    aggregate(interactionsOnEdges[,seq(from=5,to=ncol(interactionsOnEdges))], list(pair), 
              sum)
  rownames(totalInteractionsByCluster) = totalInteractionsByCluster$Group.1
  totalInteractionsByCluster = 
    totalInteractionsByCluster[,seq(from=2,to=ncol(totalInteractionsByCluster))]
  
  #get total edges per cluster pair
  totalEdges = table(pair)
  
  totalInteractionsByCluster = cbind(totalEdges[rownames(totalInteractionsByCluster)],
                                totalInteractionsByCluster)
  colnames(totalInteractionsByCluster)[seq_len(2)] = c("clusterPair", "totalEdges")
  
    meanInteractionsByCluster = totalInteractionsByCluster[,seq(from=3,to=ncol(totalInteractionsByCluster))] /
        totalInteractionsByCluster$totalEdges
  meanInteractionsByCluster = cbind(totalInteractionsByCluster[,seq_len(2)], meanInteractionsByCluster)
  
  
  
  #perform simulations
  results = list()
  for (i in seq_len(nSim)){
    permuted = permuteMatrix(M)
    sim = getInteractionsOnEdges(permuted,lrPairs,spatialGraph)
    sim = aggregate(sim[,seq(from=3,to=ncol(sim))], list(pair), sum)
    rownames(sim) = sim$Group.1
    results[[i]] = sim[,seq(from=2,to=ncol(sim))]
    if (i %% 10 == 0 & verbose){
      writeLines(as.character(i))
    }
  }
  
  #calculate summary statistics for simulation results
  results = lapply(results, function(x, y) y > x, y = 
                     totalInteractionsByCluster[,seq(from=3,to=ncol(totalInteractionsByCluster))])
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
#' @importFrom Rfast rowMins colMins
#' @import pheatmap
#' @return matrix of -log10(pvalues) that underlies the heatmap.
#' @export
#' @examples
#' getExample = make.getExample()
#' clusters = getExample('clusters')
#' colours = getExample('colours')
#' ligandReceptorResults = getExample('ligandReceptorResults')
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
      pheatmap(negLog10PValues, annotation_row = rowAnno,
               annotation_colors = list("sender" = colours, "receiver" = colours),
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
#' @param type - "total" or "mean" to plot raw total interactions or mean
#' interactions per edge.
#' @param  logScale - plot heatmap using log scale (defaults to TRUE)
#' @import pheatmap
#' @return matrix of total ligand receptor interactions that underlies t
#' he heatmap.
#' @export
#' @examples
#' getExample = make.getExample()
#' clusters = getExample('clusters')
#' ligandReceptorResults = getExample('ligandReceptorResults')
#' cellTypePerCellTypeLigRecMatrix = 
#' makeSummedLRInteractionHeatmap(ligandReceptorResults, clusters, "mean")
makeSummedLRInteractionHeatmap = function(ligandReceptorResults, clusters, type, logScale = TRUE){ 
  if (type == "total"){
    interactionsByCluster = ligandReceptorResults$totalInteractionsByCluster
  } 
  if (type == "mean"){
    interactionsByCluster = ligandReceptorResults$meanInteractionsByCluster
  } 
  summedInteractionsByCluster = rowSums(interactionsByCluster[,seq(from=3,to=ncol(interactionsByCluster))])
  pair = str_split_fixed(names(summedInteractionsByCluster), pattern = "-", 2)
  summedInteractionsByCluster = as.data.frame(cbind(pair,summedInteractionsByCluster))
  colnames(summedInteractionsByCluster) = c("Sender", "Receiver", "nInteractions")
  summedInteractionsByCluster$nInteractions = as.numeric(summedInteractionsByCluster$nInteractions)
  clusterNames = levels(clusters)
  nClusters = length(clusterNames)
  summedInteractionsByClusterMatrix = matrix(0, ncol = nClusters, nrow = nClusters)
  for (i in seq_len(nClusters)){
    for (j in seq_len(nClusters)){
        value = summedInteractionsByCluster$nInteractions[(summedInteractionsByCluster$Sender == clusterNames[i]) &
                                                          (summedInteractionsByCluster$Receiver == clusterNames[j])]
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
#' @return This returns a seurat object where 
#' each point represents an edge between cells, and spatial coordinates are the 
#' centroids of edges between cells. The "expression matrix" is the 
#' binarised presence/absence of an interaction (ligand receptor pair) on an edge.
#' Depending on the parameter returnType, this can alternatively be returned as
#' a SpatialExperiment.
#' @export
#' @examples
#' getExample = make.getExample()
#' centroids = getExample('centroids')
#' ligandReceptorResults = getExample('ligandReceptorResults')
#' edgeSeurat = computeEdgeObject(ligandReceptorResults, centroids)
computeEdgeObject = function(ligandReceptorResults, centroids, npcs = 10,
                             returnType='Seurat'){
  interactionsOnEdges = ligandReceptorResults$interactionsOnEdges
  rownames(interactionsOnEdges) = paste0(interactionsOnEdges$nodeA, "-", interactionsOnEdges$nodeB)
  interactionsOnEdgesMat = as.matrix(interactionsOnEdges[,seq(from=5,to=ncol(interactionsOnEdges))])
  interactionsOnEdgesMat= 1 * interactionsOnEdgesMat
  edgeSeurat = CreateSeuratObject(t(interactionsOnEdgesMat), meta.data = interactionsOnEdges[,seq_len(4)])
  edgeCoords = as.data.frame(cbind(centroids[interactionsOnEdges$nodeA, seq_len(2)], 
                                  centroids[interactionsOnEdges$nodeB, seq_len(2)]))
  
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
#' @importFrom Matrix sparseMatrix
#' @export
#' @examples
#' delaunayNeighbours = make.getExample()('delaunayNeighbours')
#' edgeNeighbours = computeEdgeGraph(delaunayNeighbours)

computeEdgeGraph = function(spatialGraph, selfEdges = FALSE){
  spatialGraph = data.table(spatialGraph)
  spatialGraph$edge = paste0(spatialGraph$nodeA, "-", spatialGraph$nodeB)
  spatialGraphEdgesA = merge.data.table(spatialGraph,
                             spatialGraph[,c(1,3)], by = "nodeA", allow.cartesian = TRUE)
  spatialGraphEdgesB = merge.data.table(spatialGraph,
                                        spatialGraph[,c(2,3)], by = "nodeB", allow.cartesian = TRUE)
  
  spatialGraphEdges = rbind(spatialGraphEdgesA[,c(3,4)],spatialGraphEdgesB[,c(3,4)])
  names(spatialGraphEdges) = c("nodeA","nodeB")
  if (! selfEdges){
    spatialGraphEdges = spatialGraphEdges[spatialGraphEdges$nodeA != spatialGraphEdges$nodeB,]
  }
  
  return(spatialGraphEdges)
}



 


