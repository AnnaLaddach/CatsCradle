
## ####################################################
## Functions for quality control of neighbouring cells
## ####################################################

## ####################################################
#' This function annotates edges with their distance and
#' the types of cells they connect
#'
#' @param edges - A data frame with columns nodeA and nodeB giving the
#'     cells of each edge
#' @param clusters - the clusters of each cell
#' @param centroids - the centroids of each cell
#' @return a data frame giving the edges (as nodeA and nodeB), their
#'     lengths and the cell type pair.
#' @export
#' @examples
#' getExample = make.getExample()
#' centroids = getExample('centroids')
#' clusters = getExample('clusters')
#' delaunayNeighbours = getExample('delaunayNeighbours')
#' annEdges = edgeLengthsAndCellTypePairs(delaunayNeighbours,
#'                    clusters,centroids)
edgeLengthsAndCellTypePairs = function(edges,clusters,centroids)
{
    centr = data.matrix(centroids[,seq_len(2)])
    delta = centr[edges$nodeA,] - centr[edges$nodeB,]

    getLength = function(i)
    {
        return(Norm(delta[i,]))
    }

    theRun = seq_len(nrow(delta))
    edges$length = unlist(lapply(theRun,getLength))

    getClusterPair = function(i)
    {
        thePair = c(clusters[edges$nodeA[i]],
                    clusters[edges$nodeB[i]])
        thePair = thePair[order(thePair)]
        tag = paste(thePair,collapse='_')

        return(tag)
    }

    edges$cellTypePair = unlist(lapply(theRun,getClusterPair))
    
    return(edges)                      
}


## ####################################################
#' This finds proposed cutoffs for edge lengths by clustering
#' the lengths of the edges for each cell type pair using k-means
#' clustering with k  = 2
#'
#' @param annEdges - a data frame with columns nodeA, nodeB, length
#'     and cellTypePair as produced by edgeLengthsAndCellTypePairs.
#' @return This returns a data frame with columns cellTypePair and
#'     cutoff. 
#' @export
#' @examples
#' getExample = make.getExample()
#' centroids = getExample('centroids')
#' clusters = getExample('clusters')
#' delaunayNeighbours = getExample('delaunayNeighbours')
#' annEdges =
#'     edgeLengthsAndCellTypePairs(delaunayNeighbours,clusters,centroids)
#' cutoffDF = edgeCutoffsByClustering(annEdges)
edgeCutoffsByClustering = function(annEdges)
{
    cellTypePairs = unique(annEdges$cellTypePair)
    cutoff = c()
    for(ctp in cellTypePairs)
    {
        idx = annEdges$cellTypePair == ctp
        lengths = annEdges$length[idx]
        centers = c(min(lengths),max(lengths))
        if(length(lengths) <= 2)
        {
            cutoff[ctp] = mean(lengths)
            next
        }
        km = kmeans(lengths,2)

        A = lengths[km$cluster == 1]
        B = lengths[km$cluster == 2]

        if(max(A) < min(B))
            cutoff[ctp] = (max(A) + min(B)) / 2
        else
            cutoff[ctp] = (max(B) + min(A)) / 2
    }
    df = data.frame(cellTypePair=names(cutoff),cutoff)

    return(df)
}

## ####################################################
#' This finds edge cutoffs by percentile
#'
#'@param annEdges - a data frame with columns nodeA, nodeB, length
#'     and cellTypePair as produced by edgeLengthsAndCellTypePairs.
#' @param percentileCutoff - a numeric
#' @return This returns a data frame with columns cellTypePair and
#'     cutoff.
#' @export
#' @examples
#' getExample = make.getExample()
#' centroids = getExample('centroids')
#' clusters = getExample('clusters')
#' delaunayNeighbours = getExample('delaunayNeighbours')
#' annEdges =
#'     edgeLengthsAndCellTypePairs(delaunayNeighbours,clusters,centroids)
#' cutoffDF = edgeCutoffsByPercentile(annEdges,percentileCutoff=95)
edgeCutoffsByPercentile = function(annEdges,
                                   percentileCutoff)
{
    cellTypePairs = unique(annEdges$cellTypePair)

    cutoff = c()
    for(ctp in cellTypePairs)
    {
        idx = annEdges$cellTypePair == ctp
        lengths = annEdges$length[idx]
        lengths = lengths[order(lengths)]
        num = length(lengths)
        n = ceil(percentileCutoff * num / 100)
        n = max(1,n)
        n = min(n,num)
        cutoff[ctp] = lengths[n]
    }
    cutoffDF = data.frame(cellTypePair=names(cutoff),
                          cutoff)

    return(cutoffDF)
}


## ####################################################
#' This finds edge cutoffs by z-score
#'
#'@param annEdges - a data frame with columns nodeA, nodeB, length
#'     and cellTypePair as produced by edgeLengthsAndCellTypePairs.
#' @param zCutoff - a numeric
#' @return This returns a data frame with columns cellTypePair and
#'     cutoff.
#' @export
#' @examples
#' getExample = make.getExample()
#' centroids = getExample('centroids')
#' clusters = getExample('clusters')
#' delaunayNeighbours = getExample('delaunayNeighbours') 
#' annEdges =
#'     edgeLengthsAndCellTypePairs(delaunayNeighbours,clusters,centroids)
#' cutoffDF = edgeCutoffsByZScore(annEdges,zCutoff=1.5)
edgeCutoffsByZScore = function(annEdges,zCutoff)
{
    cellTypePairs = unique(annEdges$cellTypePair)

    cutoff = c()
    for(ctp in cellTypePairs)
    {
        idx = annEdges$cellTypePair == ctp
        lengths = annEdges$length[idx]

        cutoff[ctp] = mean(lengths) + zCutoff * std(lengths)
    }
    cutoffDF = data.frame(cellTypePair=names(cutoff),
                          cutoff)

    return(cutoffDF)       
}


## ####################################################
#' This finds proposed cutoffs for edge lengths by computing
#' the histogram of edge lengths for each cell type pair and
#' then using the watershed algorithm to find the hump of the
#' histogram containing the median.
#'
#' @param annEdges - a data frame with columns nodeA, nodeB, length
#'     and cellTypePair as produced by edgeLengthsAndCellTypePairs.
#' @param tolerance - the tolerance parameter for the watershed
#'     algorithm.
#' @param nbins - the number of bins for the histogram
#' @return This returns a data frame with columns cellTypePair and
#'     cutoff.
#' @importFrom EBImage watershed
#' @export
#' @examples
#' getExample = make.getExample()
#' centroids = getExample('centroids')
#' clusters = getExample('clusters')
#' delaunayNeighbours = getExample('delaunayNeighbours') 
#' annEdges =
#'     edgeLengthsAndCellTypePairs(delaunayNeighbours,clusters,centroids)
#' cutoffDF = edgeCutoffsByWatershed(annEdges)
edgeCutoffsByWatershed = function(annEdges,nbins=15,tolerance=10)
{
    cellTypePairs = unique(annEdges$cellTypePair)
    cutoff = c()

    for(ctp in cellTypePairs)
    {
        lengths = annEdges$length[annEdges$cellTypePair == ctp]
        
        a = hist(lengths,nbins,plot=FALSE)
        d = a$density
        d = matrix(d,nrow=1)
        res = watershed(d,tolerance)

        med = median(lengths)
        ## Where median occurs:
        idx = a$mids >= med

        ## Maybe there are no bins beyond the median:
        if(sum(idx) == 0)
        {
            delta = a$mids[2] - a$mids[1]
            cutoff[ctp] = a$mids[length(a$mids)] + delta
            next
        }
        clusteredAs = res[min(which(idx),na.rm=TRUE)]
        upTo = max(which(res == clusteredAs),na.rm=TRUE)

        if(length(a$mid) == upTo)
            cutoff[ctp] = a$mid[upTo]
        else
            cutoff[ctp] = a$mid[upTo+1]
    }
    cutoffDF = data.frame(cellTypePair=names(cutoff),
                          cutoff)

    return(cutoffDF)
}


## ####################################################
#' edgeLengthPlot
#'
#' This plots histograms of the edge lengths broken out
#' by the cell types of the cells they connect.  It
#' optionally plots a cutoff for each pair of types.
#'
#' @param annEdges - A data frame as produced by
#'     edgeLengthsAndCellTypePairs
#' @param cutoffDF - A data frame with columns cellTypePair and
#'     cutoff. This defaults to NULL in which case no cutoffs will be
#'     plotted. 
#' @param whichPairs - Which cellTypePairs to plot.  If this is NULL,
#'     we plot all pairs.  If this is a numeric, we plot only pairs
#'     that have at least this many edges.  If this is a character
#'     vector, we plot the pairs in this list.
#' @param xLim - limits the extent of the plots. Defaults to 100.  Can
#'     be set to NULL.
#' @param legend - Show legend, defaults to FALSE
#' @return This returns a ggplot object
#' @export
#' @examples
#' getExample = make.getExample()
#' centroids = getExample('centroids')
#' clusters = getExample('clusters')
#' delaunayNeighbours = getExample('delaunayNeighbours') 
#' annEdges = edgeLengthsAndCellTypePairs(delaunayNeighbours,
#'                    clusters,centroids)
#' cutoffDF = edgeCutoffsByPercentile(annEdges,95)
#' g = edgeLengthPlot(annEdges,cutoffDF,whichPairs=60)
edgeLengthPlot = function(annEdges,
                          cutoffDF,
                          whichPairs,
                          xLim=100,
                          legend=FALSE)
{
    ## ####################################################
    ## What to plot:

    ## All:
    if(is.null(whichPairs))
    {
        useThesPairs = unique(annEdges$cellTypePair)
    }

    ## Sufficiently large:
    if(isa(whichPairs,'numeric'))
    {
        a = table(annEdges$cellTypePair)
        b = as.numeric(a)
        names(b) = names(a)
        useThesePairs = names(b[b >= whichPairs])
    }

    ## A pre-chosen set:
    if(isa(whichPairs,'character'))
    {
        useThesePairs = whichPairs
    }

    ## Make the plots:
    idx = annEdges$cellTypePair %in% useThesePairs
    plotDF = annEdges[idx,]

    g = ggplot(plotDF,aes(x=length,fill=cellTypePair)) +
        geom_density() +
        facet_wrap(~cellTypePair)

    ## Do we want a legend:
    if(! legend)
        g = g + theme(legend.position='none')

    ## Do we want an xlim?
    if(! is.null(xlim))
        g = g + xlim(0,xLim)

    ## Do we want cutoffs?
    ## No. We're done:
    if(is.null(cutoffDF))
        return(g)

    ## Yes.  Plot cutoffs:
    idx = cutoffDF$cellTypePair %in% useThesePairs
    g = g +
        geom_vline(data=cutoffDF[idx,],aes(xintercept=cutoff)) +
        facet_wrap(~cellTypePair)

    return(g)
}


## ####################################################
#' This subsets edges by our chosen critera
#'
#' @param annEdges - a data frame with columns nodeA, nodeB, length
#'     and cellTypePair as produced by edgeLengthsAndCellTypePairs.
#' @param cutoffSpec - This can be either a numeric value which will
#'     be applied across all edges as an upper limit or a data frame
#'     with columns cellTypePair and cutoff as produced by any of the
#'     edgeCutoffsBy functions
#' @return This returns a subset of the annotated edges
#' @export
#' @examples
#' getExample = make.getExample()
#' centroids = getExample('centroids')
#' clusters = getExample('clusters')
#' delaunayNeighbours = getExample('delaunayNeighbours') 
#' annEdges =
#'     edgeLengthsAndCellTypePairs(delaunayNeighbours,clusters,centroids)
#' tolerance = 5
#' nbins = 15
#' cutoffDFWater = edgeCutoffsByWatershed(annEdges,
#'                                       tolerance=tolerance,
#'                                       nbins=nbins)
#' culledEdges = cullEdges(annEdges,cutoffDFWater)
cullEdges = function(annEdges,cutoffSpec)
{
    if(isa(cutoffSpec,'numeric'))
    {
        idx = annEdges$length <= cutoffSpec
        annEdges = annEdges[idx,]

        return(annEdges)
    }

    ## Otherwise cutoffSpec is a data frame:
    for(i in seq_len(nrow(cutoffSpec)))
    {
        idx = annEdges$cellTypePair == cutoffSpec$cellTypePair[i] &
            annEdges$length <= cutoffSpec$cutoff[i]

        a = annEdges[idx,]

        if(i == 1)
            culledEdges = a
        else
            culledEdges = rbind(culledEdges,a)
    }

    return(culledEdges)
}


