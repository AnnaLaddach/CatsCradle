

## ####################################################
## These functions are used internally and are not
## exported.
## ####################################################

## ####################################################
## Extracts a p-value for gene set overlap from their
## sizes and the size of the background.
geneListPValue = function(A,B,C,background=25000)
{
  M = matrix(0,nrow=2,ncol=2)
  
  M[1,1] = background - A
  M[1,2] = A
  M[2,1] = B - C
  M[2,2] = C
  
  f = fisher.test(M,alternative='greater')
  
  return(f$p.value)
}


## ####################################################
## Extracts a default graph from a Seurat object.
defaultGraph = function(f)
{
    graph = paste0(f@active.assay,'_snn')
    return(graph)
}

## ####################################################
## Used to fetch umap agnostically as to upper or lower case:
fetchUMAP = function(f)
{
    umapIsCalled = names(f@reductions$umap)
    if(is.null(umapIsCalled))
        umapIsCalled = names(f@reductions$UMAP)
    df = FetchData(f,umapIsCalled)

    return(df)
}

## ####################################################
## This extracts the cells from a nearest neighbour graph.
extractCells = function(NN)
{
    cells = unique(c(NN$nodeA,NN$nodeB))
    cells = cells[order(cells)]
    return(cells)
}


## ####################################################
## The following functions are used internally to manage
## conversions between Seurat objects on the one hand and
## SingleCellExperiments and SpatialExperiments on the
## other.  These are not necessarily lossless conversions.
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
    if('counts' %in% names(assays(sce))) {
        f = as.Seurat(sce)
    } else {
        f = as.Seurat(sce,counts=NULL) }
    
    numCells = ncol(sce)
    for(n in colPairNames(sce))
    {
        df = colPair(sce,n)
        df = as.data.frame(df)

        ## Factorize:
        df[,1] = factor(df[,1])
        df[,2] = factor(df[,2])

        ## Make sparse matrix:
        M = sparseMatrix(i=as.numeric(df[,1]),
                         j=as.numeric(df[,2]),
                         x = as.numeric(df[,3]),
                         dimnames = list(levels(df[,1]), levels(df[,2])))

        ## Make graph and install:
        graph = as.Graph(M)
        f@graphs[[n]] = graph
    }

    return(f)
}

## ####################################################
SeuratToSCE = function(f,spatial=FALSE)
{
    sce = as.SingleCellExperiment(f)
    
    for(n in names(f@graphs))
    {
        NN = getNearestNeighbourLists(f,n)
        numPairs = nrow(NN)

        a = seq_len(ncol(sce))
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

        for(n in names(f@meta.data))
            colData(sce)[,n] = f@meta.data[,n]

        ## Copy in the centroids:
        centroids = GetTissueCoordinates(f)
        rownames(centroids) = centroids$cell
        centroids = centroids[,seq_len(2)]
        centroids = data.matrix(centroids)

        spatialCoords(sce) = centroids
    }
    return(sce)
}
