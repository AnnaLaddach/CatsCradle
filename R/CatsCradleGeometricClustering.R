
## ####################################################
## These functions focus on detecting geometric clustering
## of points, e.g., gene sets on a gene UMAP.
## ####################################################


## ####################################################
#' This function computes a p-value for the geometric
#' clustering of a gene set (in UMAP or PCA reduction)
#' based on the median distance from its complement to
#' the set.
#'
#' @param fPrime - a transposed Seurat object, i.e. a
#' Seurat object of genes or SingleCellExperiment to
#' be converted to a Seurat object
#' @param geneSubset - a subset of the genes which can
#' be given as a character vector as a logical vector
#' @param numTrials - the number of random trials to be
#' carried out for randomised testing. Defaults to 1000.
#' @param reduction - can be 'UMAP' or 'PCA', defaults
#' to 'UMAP'
#' @param numPCs - number of PCs to use if reduction is
#' 'PCA'
#' @return A p-value reporting how often a random subset
#' of the same size is sufficiently clustered to produce
#' an equally large distance from its complement.
#' @importFrom pracma distmat Norm ceil std
#' @importFrom stats median aggregate fisher.test kmeans p.adjust sd
#' @importFrom graphics hist
#' @importFrom methods as
#' @importFrom utils data
#' @export
#' @examples
#' getExample = make.getExample()
#' STranspose = getExample('STranspose')
#' hallmark = getExample('hallmark',toy=TRUE)
#' geneSubset = intersect(colnames(STranspose),hallmark[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]])
#' p = getObjectSubsetClusteringPValue(STranspose,geneSubset,100)
getObjectSubsetClusteringPValue = function(fPrime,
                                           geneSubset,
                                           numTrials=1000,
                                           reduction='UMAP',
                                           numPCs=10)
{
    fPrime = acceptor(fPrime)
    
    ## Test for non-empty subset:
    if(isa(geneSubset,'logical'))
        numInSubset = sum(geneSubset)
    if(isa(geneSubset,'character'))
        numInSubset = length(intersect(geneSubset,
                                       colnames(fPrime)))
    stopifnot(numInSubset > 0)
    
    stats = getObjectSubsetClusteringStatistics(fPrime,
                                                geneSubset,
                                                numTrials,
                                                reduction,
                                                numPCs)
    
    return(stats$pValue)
}

## ####################################################
#' This function computes statistics for the geometric
#' clustering of a gene set (in UMAP or PCA reduction)
#' based on the median distance from its complement to
#' the set.
#'
#' @param fPrime - a transposed Seurat object, i.e. a
#' Seurat object of genes or SingleCellExperiment to
#' be converted to a Seurat object
#' @param geneSubset - a subset of the genes which can
#' be given as a character vector or as a logical vector
#' @param numTrials - the number of random trials to be
#' carried out for randomised testing. Defaults to 1000.
#' @param reduction - can be 'UMAP' or 'PCA', defaults
#' to 'UMAP'
#' @param numPCs - number of PCs to use if reduction is
#' 'PCA'
#' @return A list of statistics resulting from the
#' testing of randomised subsets of the same size as the
#' given gene subset.  These include subsetDistance, the 
#' actual median complement distance; randomSubsetDistance,
#' the median complement distances for randomised subsets;
#' pValue, computed by comparing the real and randomised
#' distances; and zScore, the z-distance of the actual
#' median distance from the mean of the randomised distances.
#' @export
#' @examples
#' getExample = make.getExample()
#' STranspose = getExample('STranspose',toy=TRUE)
#' hallmark = getExample('hallmark')
#' geneSubset = intersect(colnames(STranspose),hallmark[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]])
#' stats = getObjectSubsetClusteringStatistics(STranspose,geneSubset,100)
getObjectSubsetClusteringStatistics = function(fPrime,
                                               geneSubset,
                                               numTrials=1000,
                                               reduction='UMAP',
                                               numPCs=10)
{
    fPrime = acceptor(fPrime)
    
    if(reduction == 'UMAP')
        S = fetchUMAP(fPrime)

    if(reduction == 'PCA')
    {
        pcs = paste0('PC_',seq_len(numPCs))
        S = FetchData(fPrime,pcs)
    }

    ## The coordinates of all the genes in the
    ## chosen reduction:
    S = as.matrix(S)
    rownames(S) = colnames(fPrime)

    answer = runGeometricClusteringTrials(S,
                                          geneSubset,
                                          numTrials)

    return(answer)
}


## ####################################################
#' This runs random trials to determine the statistical
#' significance of the clustering of a set of points
#' within a larger set.
#' 
#' This function takes a matrix whose rows are geometric
#' coordinates and a subset of these points either given
#' as a character vector which is a subset of the rownames
#' or as a logical vector.  It returns statistics on the
#' mean distance of the complement to the subset.
#'
#' @param S - a set of points given as a matrix. The rows
#' are the coordinates of these points
#' @param geneSubset - this is either a subset of the rownames of
#' S or a logical whose length is nrow(S)
#' @param numTrials - the number or random trials to perform
#' @return This returns a list. subsetDistance gives the
#' median complement distance for the actual set,
#' randomSubsetDistance gives the complement distances for
#' the numTrials random sets, pValue gives a p-value based
#' on the rank of the actual distance among the random
#' distances and zScore gives its z-score.
#' @export
#' @examples
#' library(Seurat)
#' getExample = make.getExample()
#' STranspose = getExample('STranspose',toy=TRUE)
#' hallmark = getExample('hallmark')
#' S = data.matrix(FetchData(STranspose,c('umap_1','umap_2')))
#' geneSubset = rownames(S) %in% hallmark[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]]
#' geneClustering = runGeometricClusteringTrials(S,geneSubset,100)
runGeometricClusteringTrials = function(S,
                                        geneSubset,
                                        numTrials)
{
    if(inherits(geneSubset,'character'))
        geneSubset = rownames(S) %in% geneSubset
    
    answer = list()
    answer$subsetDistance = medianComplementDistance(S,geneSubset)

    randomSubsetDistance = c()
    for(i in seq_len(numTrials))
    {
        randomSubset = sample(rownames(S),sum(geneSubset))
        randomSubset = rownames(S) %in% randomSubset
        randomSubsetDistance[i] = medianComplementDistance(S,randomSubset)
    }
    answer$randomSubsetDistance = randomSubsetDistance

    count = sum(randomSubsetDistance > answer$subsetDistance)
    count = max(count,1)

    answer$pValue = count / numTrials

    mu = mean(randomSubsetDistance)
    std = sd(randomSubsetDistance)
    answer$zScore = (answer$subsetDistance - mu) / std

     return(answer)
}


## ####################################################
## Used internally
medianComplementDistance = function(S,geneSubset)
{
    if(inherits(geneSubset,'character'))
        geneSubset = colnames(S) %in% geneSubset
    
    ## The complement:
    A = S[!geneSubset,]
    
    ## The subset:
    B = S[geneSubset,]
    
    D = distmat(A,B)
    
    rowMin = c()
    for(i in seq_len(nrow(D)))
        rowMin[i] = min(D[i,])
    d = median(rowMin)
    
    return(d)
}


## ####################################################
#' This finds the directed Hausdorf distance from A to B
#'
#' @param A - an m x d matrix representing m points in
#' dimension d
#' @param B - an n x d matrix representing n points in
#' dimension d
#' @return This returns the distance of the furthest point
#' in A from its nearest point in B.
#' @export
#' @examples
#' A = matrix(seq_len(8),ncol=2)
#' B = matrix(seq(from=3,to=16),ncol=2)
#' d_hausdorf = directedHausdorfDistance(A,B)
directedHausdorfDistance = function(A,B)
{
    D = distmat(A,B)
    
    rowMin = c()
    for(i in seq_len(nrow(D)))
        rowMin[i] = min(D[i,])
    d = max(rowMin)
    
    return(d)
}

## ####################################################
directedMedianDistance =function(A,B)
{
    D = distmat(A,B)
    
    rowMin = c()
    for(i in seq_len(nrow(D)))
        rowMin[i] = min(D[i,])
    d = median(rowMin)
    
    return(d)
}

## ####################################################
#' This takes a set S of n points in dimension d given
#' by an n x d matrix and a subset A given by a logical
#' and returns the median distance from the complement to
#' the given subset.
#'
#' @param S - an n x d matrix representing a set of n points
#' in dimension d
#' @param idx - a logical of length n representing a subset of
#' S.  This should not be the empty set or all of S.
#' @return This returns the median distance from the complement
#' to the subset
#' @export
#' @examples
#' S = matrix(seq_len(12),ncol=2)
#' idx = c(rep(FALSE,3),rep(TRUE,3))
#' compDist = medianComplementDistance(S,idx)
medianComplementDistance = function(S,idx)
{
    ## The complement:
    A = S[!idx,]
    
    ## The subset:
    B = S[idx,]
    
    return(directedMedianDistance(A,B))
}

## ####################################################
#' This takes a set S of n points in dimension d and a subset
#' A and computes a p-value for the co-localization of the subset
#' by comparing the median complement distance for the given set
#' to values of the median complement distance computed for random
#' subsets of the same size.
#'
#' @param  S - an n x d matrix representing a set of n points
#' in dimension d
#' @param idx - a logical of length n representing a subset of
#' S.  This should not be the empty set or all of S.
#' @param numTrials - the number of random trials to perform,
#' defaults to 1000
#' @param returnTrials - whether to report the real and random median
#' complement distances.
#' @return By default this reports a p-value.  If returnTrials is set,
#' this returns a list giving the p-value, the actual complement distance
#' and the random complement distances.
#' @export
#' @examples
#' library(Seurat)
#' getExample = make.getExample()
#' STranspose = getExample('STranspose',toy=TRUE)
#' hallmark = getExample('hallmark')
#' S = data.matrix(FetchData(STranspose,c('umap_1','umap_2')))
#' idx = colnames(STranspose) %in% hallmark[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]]
#' mcpv = medianComplementPValue(S,idx,numTrials=100)
medianComplementPValue = function(S,idx,numTrials=1000,returnTrials=FALSE)
{
    actual = medianComplementDistance(S,idx)
    
    numSuper = nrow(S)
    numSub = sum(idx)
    
    random = c()

    for(i in seq_len(numTrials))
    {
        IDX = rep(FALSE,numSuper)
        r = sample(numSuper,numSub)
        IDX[r] = TRUE
        random[i] = medianComplementDistance(S,IDX)
    }
    n = sum(random >= actual)
    n = max(1,n)
    
    answer = list()
    answer$pValue = n / numTrials    
    answer$medianDist = actual
    answer$random = random

    if(!returnTrials)
        return(answer$pValue)
    
    return(answer)
}

 
