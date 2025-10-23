
## ####################################################
## These functions focus of analysing ligand-receptor
## interactions between neighbouring cells is spatial
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
  
  pairsFoundDF = lrn[(lrn$from %in% panel) & (lrn$to %in% panel),]
  colnames(pairsFoundDF) = c("ligand","receptor")
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
#' This function permutes the columns of a sparse dcG matrix.
#'
#' @param M - a binarised expression matrix in sparse format where rows are cells and columns
#` are genes.
#' @return This returns a matrix in which the values have been permuted within
#' columns.
permuteColumns = function(M){
  nRow = M@Dim[1]
  rowInd = vector(mode="integer", length=length(M@i))
  j = 1
  for (i in diff(M@p)){
    if (i != 0){
      value = sort(sample(0:(nRow-1),i,replace = FALSE))
      rowInd[j:(j+i-1)] = value
      j = j + i
    }
  }
  M@i = rowInd
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
#' @return A binarised sparse expression matrix where rows are genes and 
#' columns are cells.
getBinarisedMatrix = function(obj, cutoff = 0, layer = 'counts'){
  obj = acceptor(obj)
  M = GetAssayData(obj, layer = "counts")
  M = M > 0
  rownames(M) = str_replace(rownames(M),"-",".")
  M = t(M)
  M = as(M, "dMatrix")
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

#' Given a seurat object, a spatial graph, clusters and
#' species this function identifies ligand-receptor 
#' interactions between neighbouring cells, identifies
#' ligand-receptor interactions within and between clusters
#' and calculates whether these are observed more frequently 
#' than expected by chance. If the "analytical" method is selected, an upper 
#' tail p-value for observing a given number of A-B edges positive for a given 
#' interaction is calculated using a binomial test (R pbinom function) where:
#'  q = number of A-B edges positive for an interaction
#'  size = total number of A-B edges   
#'  prob = pL*pR
#' Where pL is the probability of a cell expressing a specific ligand (number 
#' of cells positive for a ligand/total cells), 
#' and pR is the probability of a cell expressing a specific receptor (number 
#' of cells positive for a receptor/total cells). 
#' If conditional = True  p-values will be calculated given the proportion of 
#' cells that express ligands and receptors in the specific clusters (pL = 
#' number of cells in cluster A positive for a ligand/number of cells in cluster
#' A, pR = number of cells in cluster B positive for a receptor/number of cells
#' in cluster B).
#'   
#' We recommend to use the analytical method, which has a much faster runtime 
#' than the permutation-based method, however for legacy purposes and user 
#' flexibility we retain the permutation-based method.
#' @param obj - a Seurat object
#' @param spatialGraph - a data frame of neighbouring
#' cell pairs. 
#' @param clusters - named vector of clusters where names are each cell and
#' clusters are a factor
#' @param species - either 'human' or 'mouse'
#' @param method - method for computing p-values. Defaults to "analytical". 
#' If "permutation" is selected p-values are calculated by comparison to 
#' randomised graphs (note this is slower than the analytical approach).
#' @param conditional - if method is "analytical" and conditional is true, 
#' p-values will be calculated given the proportion of cells that express 
#' ligands and receptors in the specific clusters. Otherwise global proportions 
#' of ligand and receptor expression are used. Defaults to FALSE.
#' @param minEdgesPos - the minimum edges that need to be positive for a 
#' ligand-receptor interaction between two clusters for a p-value to be 
#' calculated. Only taken into consideration when the analytical method is 
#' selected.   
#' @param nSim - number of simulations to perform for pvalue calculation.
#' @param lrn - a ligand-receptor network, i.e., a
#' data frame with columns from and to.  By default, it
#' retrieves the nichenetr ligand receptor network
#' @param verbose - whether to print trace, defaults to TRUE
#' @return A list containing:
#' interactionsOnEdges - a sparse matrix where the rownames give pairs of 
#' neighbouring cells and column names give ligand-receptor pairs. 
#' Entries are TRUE if the ligand is expressed in the first cell and the receptor is expressed in 
#' the second cell and FALSE if not.
#' interactionsOnEdgesMeta - a dataframe where the first two columns are the
#' cells that comprise the edges in interactionsOnEdges, and the next two
#' columns are their clusters.
#' totalInteractionsByCluster - a dataframe where the rownames are 
#' sender-receiver cluster pairs and column names are ligand receptor pairs. 
#' Entries are total numbers of edges on which particular ligand receptor 
#' interactions are present.
#' meanInteractionsByCluster - a dataframe where the rownames are 
#' sender-receiver cluster pairs and column names are ligand receptor pairs. 
#' Entries are total numbers of edges on which particular ligand receptor i
#' nteractions are present (for that cluster pair) divided by the total number 
#' of edges between those clusters.
#' simResults - a dataframe where the rownames are sender-receiver cluster 
#' pairs and column names are ligand receptor pairs. 
#' Values give the number of simulations for which observed values are greater 
#' than simulated values. Only returned if method = "permutation".
#' pValues - a dataframe where the rownames are sender-receiver cluster pairs 
#' and column names are ligand receptor pairs. Entries are uppertail p-values 
#' describing whether a particular ligand receptor interaction is observed more 
#' frequently between 2 clusters than expected.
#' totalEdges - a vector of total edges between cluster pairs.
#' @importFrom stringr str_split_1
#' @export
#' @examples
#' getExample = make.getExample()
#' smallXenium = getExample('smallXenium')
#' delaunayNeighbours = getExample('delaunayNeighbours')
#' clusters = getExample('clusters')
#' performLigandReceptorAnalysis(smallXenium, delaunayNeighbours, 
#'                                       "mouse", clusters)

performLigandReceptorAnalysis = function(obj, spatialGraph, species, clusters,
                                      method = "analytical",
                                      conditional = FALSE,
                                      minEdgesPos = 10,
                                      nSim = 1000, 
                                      lrn = getLigandReceptorNetwork(species),
                                      verbose=TRUE){
  
  stopifnot(method %in% c('analytical','permutation')) 
  if (method == "analytical"){
    results = performLigandReceptorAnalysisAnalytical(obj, spatialGraph, 
                                            species, clusters, 
                                            conditional = conditional,
                                            lrn = lrn,
                                            minEdgesPos = minEdgesPos)
  }
  if (method == "permutation"){
    results = performLigandReceptorAnalysisPermutation (obj, spatialGraph, 
                                              species, clusters, nSim = 1000, 
                                              lrn = lrn,  minEdgesPos = minEdgesPos,
                                              verbose=verbose)
  }
  return(results)
}


## ####################################################
#' Given a seurat object, a spatial graph, clusters and
#' species this function identifies ligand-receptor 
#' interactions between neighbouring cells, identifies
#' ligand-receptor interactions within and between clusters
#' and calculates whether these are observed more frequently 
#' than expected by chance using a permutation-based approach. 
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
#' @param minEdgesPos - the minimum edges that need to be positive for a 
#' ligand-receptor interaction between two clusters for a p-value to be 
#' calculated. Only taken into consideration when the analytical method is 
#' selected.   
#' @param verbose - whether to print trace, defaults to TRUE
#' @return A list containing:
#' interactionsOnEdges - a sparse matrix where the rownames give pairs of 
#' neighbouring cells and column names give ligand-receptor pairs. 
#' Entries are TRUE if the ligand is expressed in the first cell and the receptor is expressed in 
#' the second cell and FALSE if not.
#' interactionsOnEdgesMeta - a dataframe where the first two columns are the
#' cells that comprise the edges in interactionsOnEdges, and the next two
#' columns are their clusters.
#' totalInteractionsByCluster - a dataframe where the rownames are 
#' sender-receiver cluster pairs and column names are ligand receptor pairs. 
#' Entries are total numbers of edges on which particular ligand receptor 
#' interactions are present.
#' meanInteractionsByCluster - a dataframe where the rownames are 
#' sender-receiver cluster pairs and column names are ligand receptor pairs. 
#' Entries are total numbers of edges on which particular ligand receptor i
#' nteractions are present (for that cluster pair) divided by the total number 
#' of edges between those clusters.
#' simResults - a dataframe where the rownames are sender-receiver cluster 
#' pairs and column names are ligand receptor pairs. 
#' Values give the number of simulations for which observed values are greater 
#' than simulated values. 
#' pValues - a dataframe where the rownames are sender-receiver cluster pairs 
#' and column names are ligand receptor pairs. Entries are uppertail p-values 
#' describing whether a particular ligand receptor interaction is observed more 
#' frequently between 2 clusters than expected.
#' totalEdges - a vector of total edges between cluster pairs.
#' @examples
#' getExample = make.getExample()
#' smallXenium = getExample('smallXenium')
#' delaunayNeighbours = getExample('delaunayNeighbours')
#' clusters = getExample('clusters')
#' performLigandReceptorAnalysis(smallXenium, delaunayNeighbours, 
#'                                       "mouse", clusters,  minEdgesPos = 10, nSim = 10,
#'                                        verbose=FALSE)

performLigandReceptorAnalysisPermutation = function(obj, spatialGraph, species, 
                                              clusters,
                                              nSim = 1000, 
                                              lrn = getLigandReceptorNetwork(species),
                                              minEdgesPos = 10,
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
  M = M[,colnames(M) %in% c(lrPairs$ligand,lrPairs$receptor)]
  
  #get ligand receptor interactions between cells
  interactionsOnEdges = M[spatialGraph$nodeA,lrPairs$ligand] &
    M[spatialGraph$nodeB,lrPairs$receptor]
  interactionsOnEdges = t(interactionsOnEdges)
  
  ligRecNames = paste(lrPairs$ligand,lrPairs$receptor,sep='_')
  cellNames = paste(spatialGraph$nodeA, spatialGraph$nodeB, sep = "_")
  
  pair = 
    paste0(clusters[spatialGraph$nodeA],  "-", clusters[spatialGraph$nodeB])
  interactionsOnEdgeslgT = as(interactionsOnEdges, "TsparseMatrix")
  interactionsOnEdgesAnno = cbind(pair[(interactionsOnEdgeslgT@j+1)],ligRecNames[(interactionsOnEdgeslgT@i+1)])
  totalInteractionsByCluster = table(interactionsOnEdgesAnno[,1],interactionsOnEdgesAnno[,2])
  totalInteractionsByCluster = as.data.frame.matrix(totalInteractionsByCluster)
  
  #get total edges per cluster pair
  totalEdges = table(pair)
  
  meanInteractionsByCluster = totalInteractionsByCluster/
    totalEdges[rownames(totalInteractionsByCluster)]
  
  interactionsOnEdges = t(interactionsOnEdges)
  colnames(interactionsOnEdges) = ligRecNames
  rownames(interactionsOnEdges) = cellNames
  
  #perform simulations
  results = list()
  clusterPairs = rownames(totalInteractionsByCluster)
  ligRecPairs = colnames(totalInteractionsByCluster)
  
  for (i in seq_len(nSim)){
    permuted = permuteColumns(M)
    sim = permuted[spatialGraph$nodeA,lrPairs$ligand] &
      permuted[spatialGraph$nodeB,lrPairs$receptor]
    sim = t(sim)
    sim = as(sim, "TsparseMatrix")
    sim = cbind(pair[(sim@j+1)],ligRecNames[(sim@i+1)])
    sim = table(sim[,1],sim[,2])
    sim = as.data.frame.matrix(sim)
    
    missingRows = clusterPairs[!(clusterPairs %in% rownames(sim))]
    missingColumns = ligRecPairs[!(ligRecPairs %in% colnames(sim))]
    
    for (row in missingRows){
      sim[row,] = rep(0, ncol(sim))
    }
    for (col in missingColumns){
      sim[,col] = rep(0, nrow(sim))
    }
    
    sim = sim[clusterPairs,ligRecPairs]
    results[[i]] = sim
    
    if (i %% 10 == 0 & verbose){
      writeLines(as.character(i))
    }
  }
  interactionsOnEdgesMetaData = cbind(spatialGraph, clusters[spatialGraph$nodeA], clusters[spatialGraph$nodeB])
  names(interactionsOnEdgesMetaData)[3:4] = c("cellTypeA","cellTypeB")
  if (nSim > 0){
    #calculate summary statistics for simulation results
    results = lapply(results, function(x, y) y > x, y = 
                       totalInteractionsByCluster)
    results = abind(results, along = 3L)
    simResults = rowSums(results, dims = 2)
    rownames(simResults) = rownames(totalInteractionsByCluster)
    pValues = abs((simResults - nSim)/nSim) 
    pValues = pmax(pValues, (1/nSim))
    pValues[totalInteractionsByCluster < minEdgesPos] = NA
    return(list("interactionsOnEdges" = interactionsOnEdges, 
                "interactionsOnEdgesMeta" = interactionsOnEdgesMetaData,
                "totalInteractionsByCluster" = totalInteractionsByCluster,
                "meanInteractionsByCluster" = meanInteractionsByCluster,
                "simResults" = as.data.frame(simResults),
                "pValues" = as.data.frame(pValues),
                "totalEdges" = totalEdges))
  } else {
    return(list("interactionsOnEdges" = interactionsOnEdges,
                "interactionsOnEdgesMeta" = interactionsOnEdgesMetaData,
                "totalInteractionsByCluster" = totalInteractionsByCluster,
                "meanInteractionsByCluster" = meanInteractionsByCluster,
                "totalEdges" = totalEdges))
  }
}

## ####################################################
#' Given a seurat object, a spatial graph, clusters and
#' species this function identifies ligand-receptor 
#' interactions between neighbouring cells, identifies
#' ligand-receptor interactions within and between clusters
#' and calculates whether these are observed more frequently 
#' than expected by chance using an analytical approach.
#'
#' @param obj - a Seurat object
#' @param spatialGraph - a data frame of neighbouring
#' cell pairs. 
#' @param clusters - named vector of clusters where names are each cell and
#' clusters are a factor
#' @param species - either 'human' or 'mouse'
#' @param conditional - if method is "analytical" and conditional is true, 
#' p-values will be calculated given the proportion of cells that express 
#' ligands and receptors in the specific clusters. Otherwise global proportions 
#' of ligand and receptor expression are used. Defaults to FALSE.
#' @param lrn - a ligand-receptor network, i.e., a
#' data frame with columns from and to.  By default, it
#' retrieves the nichenetr ligand receptor network
#' @param minEdgesPos - the minimum edges that need to be positive for a 
#' ligand-receptor interaction between two clusters for a p-value to be 
#' calculated. 
#' @importFrom stats pbinom phyper
#' @return A list containing:
#' interactionsOnEdges - a sparse matrix where the rownames give pairs of 
#' neighbouring cells and column names give ligand-receptor pairs. 
#' Entries are TRUE if the ligand is expressed in the first cell and the receptor is expressed in 
#' the second cell and FALSE if not.
#' interactionsOnEdgesMeta - a dataframe where the first two columns are the
#' cells that comprise the edges in interactionsOnEdges, and the next two
#' columns are their clusters.
#' totalInteractionsByCluster - a dataframe where the rownames are 
#' sender-receiver cluster pairs and column names are ligand receptor pairs. 
#' Entries are total numbers of edges on which particular ligand receptor 
#' interactions are present.
#' meanInteractionsByCluster - a dataframe where the rownames are 
#' sender-receiver cluster pairs and column names are ligand receptor pairs. 
#' Entries are total numbers of edges on which particular ligand receptor 
#' interactions are present (for that cluster pair) divided by the total number 
#' of edges between those clusters.
#' pValues - a dataframe where the rownames are sender-receiver cluster pairs 
#' and column names are ligand receptor pairs. Entries are uppertail p-values 
#' describing whether a particular ligand receptor interaction is observed more 
#' frequently between 2 clusters than expected.
#' totalEdges - a vector of total edges between cluster pairs.
performLigandReceptorAnalysisAnalytical = function(obj, spatialGraph, species, clusters,
                                                    conditional = FALSE,
                                                    lrn = getLigandReceptorNetwork(species),
                                                    minEdgesPos = 10)
  {
  
  #symmetrise spatial graph
  spatialGraph = spatialGraph[,c(1,2)]
  spatialGraphBA = spatialGraph[,c(2,1)]
  names(spatialGraphBA) = c("nodeA","nodeB")
  spatialGraph = rbind(spatialGraph,spatialGraphBA)
  spatialGraph = unique(spatialGraph)
  
  #get ligand receptor pairs
  lrPairs = getLigandReceptorPairsInPanel(obj, species)
  
  #get binarised expression matrix for ligand receptor pairs
  M = getBinarisedMatrix(obj)
  M = M[,colnames(M) %in% c(lrPairs$ligand,lrPairs$receptor)]
  
  #get ligand receptor interactions between cells
  interactionsOnEdges = M[spatialGraph$nodeA,lrPairs$ligand] &
    M[spatialGraph$nodeB,lrPairs$receptor]
  interactionsOnEdges = t(interactionsOnEdges)
  
  ligRecNames = paste(lrPairs$ligand,lrPairs$receptor,sep='_')
  cellNames = paste(spatialGraph$nodeA, spatialGraph$nodeB, sep = "_")
  
  #get sum of interactions within and between clusters
  pair = 
    paste0(clusters[spatialGraph$nodeA],  "-", clusters[spatialGraph$nodeB])
  interactionsOnEdgeslgT = as(interactionsOnEdges, "TsparseMatrix")
  interactionsOnEdgesAnno = cbind(pair[(interactionsOnEdgeslgT@j+1)],ligRecNames[(interactionsOnEdgeslgT@i+1)])
  totalInteractionsByCluster = table(interactionsOnEdgesAnno[,1],interactionsOnEdgesAnno[,2])
  totalInteractionsByCluster = as.data.frame.matrix(totalInteractionsByCluster)
  
  #get total edges per cluster pair
  totalEdges = table(pair)
  
  meanInteractionsByCluster = totalInteractionsByCluster/
    totalEdges[rownames(totalInteractionsByCluster)]
  
  interactionsOnEdges = t(interactionsOnEdges)
  colnames(interactionsOnEdges) = ligRecNames
  rownames(interactionsOnEdges) = cellNames
  
  
  geneTotals = colSums(M)
  clusterTotals = table(clusters)
  totalCells = length(clusters)
  
  #calculate stats
  results = list()
  clusterPairs = rownames(totalInteractionsByCluster)
  ligRecPairs = colnames(totalInteractionsByCluster)
  
  genes = colnames(M)
  MlgT = as(M, "TsparseMatrix")
  MAnno = cbind(as.character(clusters)[(MlgT@i+1)],genes[(MlgT@j+1)])
  clusterMatrix = table(MAnno[,1],MAnno[,2])
  clusterMatrix = as.data.frame.matrix(clusterMatrix)
  
  pValues = matrix(NA, nrow = nrow(totalInteractionsByCluster), ncol = ncol(totalInteractionsByCluster))
  rownames(pValues) = rownames(totalInteractionsByCluster)
  colnames(pValues) = colnames(totalInteractionsByCluster)
  for (CTPair in rownames(totalInteractionsByCluster)){
    CTs = str_split_1(CTPair, "-")
    for (genePair in colnames(totalInteractionsByCluster)){
      genes = str_split_1(genePair, "_")
      if (conditional){
        pG1 = clusterMatrix[CTs[1],genes[1]]/clusterTotals[CTs[1]]
        pG2 = clusterMatrix[CTs[2],genes[2]]/clusterTotals[CTs[2]]
      } else {
      pG1 = geneTotals[genes[1]]/totalCells
      pG2 = geneTotals[genes[2]]/totalCells
      }
      nEdgePos =   totalInteractionsByCluster[CTPair,genePair] 
      if(nEdgePos >= minEdgesPos){
        pValues[CTPair,genePair] = pbinom(nEdgePos,totalEdges[CTPair],pG1*pG2, lower.tail = F)
      }
    }
  } 
  
  interactionsOnEdgesMetaData = cbind(spatialGraph, clusters[spatialGraph$nodeA], clusters[spatialGraph$nodeB])
  names(interactionsOnEdgesMetaData)[3:4] = c("cellTypeA","cellTypeB")
  return(list("interactionsOnEdges" = interactionsOnEdges, 
              "interactionsOnEdgesMeta" = interactionsOnEdgesMetaData, 
              "totalInteractionsByCluster" = totalInteractionsByCluster,
              "meanInteractionsByCluster" = meanInteractionsByCluster,
              "pValues" = as.data.frame(pValues),
              "totalEdges" = totalEdges))
}



## ####################################################
#' This function takes ligandReceptorResults and plots a heatmap of -log10(pvalues).
#' If the minimum p-value is 0 a pseudocount of 0.001 will be added before log 
#' transformation.
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
#' ligandReceptorResults = getExample('ligandReceptorResults')
#' cellTypePerCellTypeLigRecMatrix = 
#' makeSummedLRInteractionHeatmap(ligandReceptorResults, clusters, "mean")
makeLRInteractionHeatmap = function(ligandReceptorResults,
                                  clusters,
                                  colours = c(),
                                  pValCutoffClusterPair = 0.05, 
                                  pValCutoffLigRec = 0.05,
                                  labelClusterPairs = TRUE)
{ 
  pValues = as.matrix(ligandReceptorResults$pValues)
  if (min(pValues, na.rm = T) == 0){
    pValues = pValues + 0.001 
  }
  pValuesMod = pValues
  pValuesMod[is.na(pValuesMod)] = 1
  
  selectedPValues = pValues[rowMins(pValuesMod, value = TRUE) < pValCutoffClusterPair,
                            colMins(pValuesMod, value = TRUE) < pValCutoffLigRec]
  selectedPValuesMod = pValuesMod[rowMins(pValuesMod, value = TRUE) < pValCutoffClusterPair,
                                  colMins(pValuesMod, value = TRUE) < pValCutoffLigRec]

  negLog10PValues = -log10(selectedPValues + 0.001)
  negLog10PValuesMod = -log10(selectedPValuesMod + 0.001)

  rowAnno = str_split_fixed(rownames(selectedPValues), pattern = "-", 2)
  rowAnno = as.data.frame(rowAnno)
  names(rowAnno) = c("sender","receiver")
  rowAnno$sender = factor(rowAnno$sender, levels = levels(clusters))
  rowAnno$receiver = factor(rowAnno$receiver, levels = levels(clusters))
  rownames(rowAnno) = rownames(selectedPValues)
  rowAnno = rowAnno[,c("receiver","sender")]
  p1 = pheatmap(negLog10PValuesMod, annotation_row = rowAnno, show_rownames = labelClusterPairs)
  if (length(colours) > 0){
    p2 = pheatmap(negLog10PValues, annotation_row = rowAnno,  annotation_colors = list("sender" = colours, "receiver" = colours),
          show_rownames = labelClusterPairs, cluster_rows = F,cluster_cols = F, na_col = "grey")
    p2$tree_row = p1$tree_row
    p2$tree_col = p1$tree_col
    print(p2)
  } else{
    p2 = pheatmap(negLog10PValues, annotation_row = rowAnno, show_rownames = labelClusterPairs, cluster_rows = F,cluster_cols = F, na_col = "grey")
    p2$tree_row = p1$tree_row
    p2$tree_col = p1$tree_col
    print(p2)
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
  summedInteractionsByCluster = rowSums(interactionsByCluster)
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
#' This is a utility function for converting entries in ligandReceptorResults 
#' to long format. 
#'
#' @param data - item from ligandReceptorResults
#' @param name - name to give column of returned data
#' @importFrom reshape2 melt
#' @return dataframe with item from ligandReceptorResults in long format
formatData = function(data, name){
  data$clusterPair = rownames(data)
  data = melt(data, id.vars = c("clusterPair"), variable.name = "interaction")
  names(data)[3] = name
  return(data)
}


## ####################################################
#' This is a utility function for converting ligandReceptor cluster-level 
#' results to long format and calculates adjusted p-values.
#'
#' @param ligandReceptorResults - ligandReceptorReults calculated using 
#' performLigandReceptorAnalysis()
#' @importFrom stringr str_split_fixed
#' @return ligand receptor results in long format

convertToLong = function(ligandReceptorResults){
  meanLong = formatData(ligandReceptorResults$meanInteractionsByCluster, "mean")
  totalLong = formatData(ligandReceptorResults$totalInteractionsByCluster, "total")
  pvalLong = formatData(ligandReceptorResults$pValues, "pValue")
  resultsLong = cbind(totalLong, meanLong$mean, pvalLong$pValue)
  resultsLong = cbind(resultsLong, str_split_fixed(resultsLong$clusterPair, pattern = "-", 2))
  names(resultsLong)[4:7] = c("mean","pValue","sender","receiver")
  resultsLong$negLog10PValue = -log10(resultsLong$pValue + 0.001)
  resultsLong$padj = p.adjust(resultsLong$pValue, method = "fdr")
  return(resultsLong)
}


## ####################################################
#' This is a function to create a dotplot using the ligand receptor results
#'
#' @param ligandReceptorResults - ligandReceptorResults calculated using 
#' performLigandReceptorAnalysis().
#' @param senderClusters - sender clusters to plot (defaults to all).
#' @param receiverClusters - receiver clusters to plot (defaults to all).
#' @param padjCutoff - only plot results with p-adj < padjCutoff (defaults to 
#' 0.05).
#' @param pvalCutoff - only plot results with p-value < pvalCutoff (defaults to 
#' False in which case padjCutoff is used).
#' @param pvalCutoff - only plot results with p-value < pvalCutoff (defaults to 
#' False in which case padjCutoff is used).
#' @param splitBy - split plots by "sender" or "receiver" cell types (defaults 
#' to sender).
#' @importFrom stringr str_split_fixed
#' @import ggplot2
#' @return matrix of total ligand receptor interactions that underlies the heatmap.
#' @export
#' @examples
#' getExample = make.getExample()
#' centroids = getExample('centroids')
#' ligandReceptorResults = getExample('ligandReceptorResults')
#' p = plotLRDotplot(ligandReceptorResults)
plotLRDotplot = function(ligandReceptorResults, senderClusters = unique(ligandReceptorResults$interactionsOnEdgesMeta$cellTypeA),
                         receiverClusters = unique(ligandReceptorResults$interactionsOnEdgesMeta$cellTypeB),  padjCutoff = 0.05,pvalCutoff = F, splitBy = "sender"){
  ligandReceptorResultsLong = convertToLong(ligandReceptorResults)
  resultsLongSelected = ligandReceptorResultsLong[(ligandReceptorResultsLong$sender %in% senderClusters)
                                            & (ligandReceptorResultsLong$receiver %in%  receiverClusters),]
  resultsLongSelected = resultsLongSelected[!(is.na(resultsLongSelected$padj)),]
  if (pvalCutoff){
    resultsLongSelected = resultsLongSelected[resultsLongSelected$pValue < pvalCutoff,]
  } else {resultsLongSelected = resultsLongSelected[resultsLongSelected$padj < padjCutoff,]}
  if (splitBy == "receiver"){
  p = ggplot(resultsLongSelected, aes(x=interaction, y=sender)) + geom_point(aes(size = mean,color = negLog10PValue)) +  facet_wrap(~receiver) +
    theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_colour_gradient(low = "#ffcccb", high = "#8B0000") +
    guides(color=guide_legend(title="-Log10(p-val + 0.001)"), size=guide_legend(title="Mean")) + ylab("Sender") + xlab("") + ggtitle("Receiver")
  }
  if (splitBy == "sender"){
    p = ggplot(resultsLongSelected, aes(x=interaction, y=receiver)) + geom_point(aes(size = mean,color = negLog10PValue)) +  facet_wrap(~sender) +
    theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_colour_gradient(low = "#ffcccb", high = "#8B0000") +
    guides(color=guide_legend(title="-Log10(p-val + 0.001)"), size=guide_legend(title="Mean")) + ylab("Receiver") + xlab("") + ggtitle("Sender")
  }
  print(p)
  return(p)
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
  
  #rownames(interactionsOnEdges) = paste0(interactionsOnEdges$nodeA, "-", interactionsOnEdges$nodeB)
  #interactionsOnEdgesMat = as.matrix(interactionsOnEdges[,seq(from=5,to=ncol(interactionsOnEdges))])
  interactionsOnEdges= 1 * interactionsOnEdges
  edgeSeurat = CreateSeuratObject(t(interactionsOnEdges), meta.data = ligandReceptorResults$interactionsOnEdgesMeta)
  edgeCoords = as.data.frame(cbind(centroids[ligandReceptorResults$interactionsOnEdgesMeta$nodeA, seq_len(2)], 
                                  centroids[ligandReceptorResults$interactionsOnEdgesMeta$nodeB, seq_len(2)]))
  
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



 


