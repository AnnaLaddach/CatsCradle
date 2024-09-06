
## ####################################################
## These functions transpose a single cell object and
## use this to cluster and annotate genes.
## ####################################################



## ###################################################
#' Create the transpose of a Seurat object
#'
#' This takes a Seurat object f and creates a new Seurat object whose
#' expression matrix is the transpose of that of f.  This can
#' also be a SingleCellExperiment which will be converted to a
#' Seurat object
#'
#' @param f - a Seurat object
#' @param active.assay - the assay to use.  Defaults to 'RNA'
#' @param npcs - number of principal components, defaults to 30
#' @param dims - dimensions to use for umap and nearest neighbors,
#'     defaults to 1:20
#' @param res - the clustering resolution, defaults to 1
#' @param returnType - Will return a SingleCellExperiment if this is either
#' of SCE, SingleCellExperiment or their lower-case equivalents.  Otherwise,
#' returns a Seurat object
#' @param verbose - Controls whether to display trace from the Seurat
#'     functions. Defaults to FALSE
#' @return A Seurat object or SingleCellExperiment
#' @export
#' @import Seurat
#' @examples
#' exSeuratObj = make.getExample()('exSeuratObj',toy=TRUE)
#' STranspose = transposeObject(exSeuratObj)
#' STransposeAsSCE = transposeObject(exSeuratObj,returnType='SCE')
transposeObject = function(f,active.assay='RNA',
                           npcs=30,dims=seq_len(20),res=1,
                           returnType='Seurat',
                           verbose=FALSE)
 {
    f = acceptor(f)

    data = f[[active.assay]]$data
    data = t(data)
    rownames(data) = str_replace_all(rownames(data),'_','-')

    fPrime = CreateSeuratObject(data,assay=active.assay)
    fPrime[[active.assay]]$data = data
    ## fPrime[[active.assay]]$counts = NULL

    fPrime = ScaleData(fPrime,verbose=verbose)
    fPrime = FindVariableFeatures(fPrime,assay=active.assay,verbose=verbose)
    fPrime = RunPCA(fPrime,npcs=npcs,verbose=verbose)
    fPrime = RunUMAP(fPrime,reduction='pca',dims=dims,verbose=verbose)
    fPrime = FindNeighbors(fPrime,verbose=verbose)
    fPrime = FindClusters(fPrime,resolution=res,verbose=verbose)

    return(returnAs(fPrime,returnType))
 }

## ####################################################
#' This compares the gene clusters to other gene sets
#' e.g., GO, Hallmark, and determines the p-value for
#' their overlaps when compared to a set of background
#' genes.
#'
#' @param geneSets - a named list of gene sets
#' @param clusterDF - a data frame giving the cluster
#' membership of each gene with columns gene and geneCluster
#' @param backgroundGenes - a character vector of genes
#' @param adjust - a logical deciding whether to adjust
#' p values.  Defaults to FALSE.
#' @return a matrix of p-values rows correspond to the gene
#' sets and the columns correspond the the CatsCradle gene
#' clusters
#' @export
#' @examples
#' getExample = make.getExample()
#' STranspose = getExample('STranspose',toy=TRUE)
#' clusterDF = data.frame(gene=colnames(STranspose),
#'                        geneCluster=STranspose$seurat_clusters)
#' hallmark = getExample('hallmark',toy=TRUE)
#' geneSet = intersect(hallmark[[1]],colnames(STranspose))
#' pvalueMatrix = geneSetsVsGeneClustersPValueMatrix(geneSet,
#'                                               clusterDF,
#'                                               colnames(STranspose))
geneSetsVsGeneClustersPValueMatrix = function(geneSets,
                                              clusterDF,
                                              backgroundGenes,
                                              adjust=FALSE)
{
    background = length(backgroundGenes)
    clusters = as.character(unique(clusterDF$geneCluster))
    clusters = as.numeric(unique(clusterDF$geneCluster))
    clusters = clusters[order(clusters)]
    NClusters = length(clusters)
    NGeneSets = length(geneSets)
    
    M = matrix(0,nrow=NGeneSets,ncol=NClusters)
    rownames(M) = names(geneSets)
    colnames(M) = as.character(clusters)
    
    for(i in seq_len(NGeneSets))
    {
        for(j in seq_len(NClusters))
        {
            cluster = clusters[j]
            idx = clusterDF$geneCluster == cluster
            clusterGenes = clusterDF$gene[idx]
            
            A = length(geneSets[[i]])
            B = length(clusterGenes)
            C = length(intersect(geneSets[[i]],clusterGenes))
            
            ## Maybe there's a better way to do this:
            M[i,j] = geneListPValue(A,B,C,background)
        }
    }
    
    if(adjust)
    {
        N = p.adjust(M,method='fdr')
        N = matrix(N,nrow=nrow(M))
        rownames(N) = rownames(M)
        colnames(N) = colnames(M)
        
        return(N)
    } else {
        return(M)
    }
}    

## ####################################################
#' This function gets the neighbors of a given gene
#' using either the gene Seurat object or its nearest
#' neighbor graph returned from getNearestNeighbourLists
#'
#' @param gene - the gene in question
#' @param NN - either the gene Seurat object or its nearest
#' neighbor graph as found by getNearestNeighbourLists.
#' This can also be a SingleCellExperiment which will be converted
#' to a Seurat object
#' @return the neighboring genes
#' @export
#' @examples
#' library(Seurat)
#' getExample = make.getExample()
#' STranspose = getExample('STranspose',toy=TRUE)
#' NN = getExample('NN',toy=TRUE)
#' neighbors = getGeneNeighbors("Ccl6",STranspose)
#' neighborsAgain = getGeneNeighbors("Ccl6",NN)
getGeneNeighbors = function(gene,NN)
{
    ## Somewhat of a fudge when NN is a NN graph:
    NN = acceptor(NN)
    
    if(inherits(NN,'Seurat'))
        NN = getNearestNeighbourLists(NN)

    NN = symmetriseNN(NN)

    idx = NN$nodeA == gene

    return(NN$nodeB[idx])
}


## ####################################################
#' Discovers the combinatorial ball of a given radius
#' around a fixed set of genes in the nearest neighbor
#' graph of a Seurat object.
#'
#' @param NN - a nearest neighbors graph
#' @param origin - a gene or list of genes
#' @param radius - the radius of the combinatorial ball
#' to be found.
#'
#' @return This returns a data frame whose columns are the
#' gene name, the radius from the origin at which it is found
#'
#' @export
#' @examples
#' getExample = make.getExample()
#' NN = getExample('NN',toy=TRUE)
#' STranspose = getExample('STranspose',toy=TRUE)
#' spheres = combinatorialSpheres(NN,'Ccl6',3)
#' hallmark = getExample('hallmark',toy=TRUE)
#' geneSet = intersect(hallmark[[1]],colnames(STranspose))
#' sphereAroundSet = combinatorialSpheres(NN,geneSet,1)
combinatorialSpheres = function(NN,origin,radius)
{
    NN = symmetriseNN(NN)
    
    nodes = unique(NN$nodeA)
    ball = data.frame(nodes,
                      radius=-1)
    rownames(ball) = nodes

    ## Initialise:
    ball[origin,'radius'] = 0

    numGenes = length(nodes)
    numFound = sum(ball$radius > -.5)
    ## Iterate:
    r = 0
    while(numFound < numGenes)
    {
        r = r + 1
        if(is.numeric(radius) &
           r > radius)
            break
           
        writeLines(paste('radius',r))
        ## Which vertices can we sprout from:
        fromThese = rownames(ball[ball$radius==r-1,])

        ## Candidate nodes:
        idx = NN$nodeA %in% fromThese
        toThese = NN$nodeB[idx]
        toThese = unique(toThese)

        ## Which ones are new:
        idx = rownames(ball) %in% toThese &
            ball$radius == -1

        ball$radius[idx] = r

        newFound = sum(ball$radius > -.5)
        writeLines(paste('   ',newFound))

        ## Can only happen in a disconnected graph
        ## or with radius too large:
        if(newFound == numFound)
            break
    }

    ## Trim:
    idx = ball$radius != -1
    ball = ball[idx,]

    ## Order:
    ball = ball[order(rownames(ball)),]
    ball = ball[order(ball$radius),]

    return(ball)
}


## ####################################################
#' This function annotates genes with terms
#'
#' This essentially inverts a list of gene sets.  It takes
#' a list (e.g., Hallmark or GO) where each list item is a
#' name of a gene set and gives the genes in that set and
#' returns a list where each item is a gene and gives the
#' gene sets that gene is in.
#'
#' @param geneSets - a list of gene sets, e.g., as produced by
#' readGmt
#' @return - A list where names are genes and values are lists of terms
#' @export
#' @examples
#' hallmark = make.getExample()('hallmark')
#' annotatedGenes = annotateGenesByGeneSet(hallmark)
annotateGenesByGeneSet = function(geneSets)
{
    genes = unique(unlist(geneSets))
    genes = genes[order(genes)]

    numSets = length(geneSets)
    numGenes = length(genes)

    M = matrix(NA,nrow=numGenes,ncol=numSets)
    rownames(M) = genes
    colnames(M) = names(geneSets)
    

    for(geneSet in names(geneSets))
        for(gene in geneSets[[geneSet]])
            M[gene,geneSet] = geneSet

    genesAnno = list()
    for(gene in genes)
    {
        genesAnno[[gene]] = M[gene,]
        genesAnno[[gene]] = genesAnno[[gene]][!is.na(genesAnno[[gene]])]
    }
    return(genesAnno)   
}


## ####################################################
#' This function returns a numeric indicating which gene
#' sets it does and does not belong to.  This vector can
#' be normalised to account for the sizes of the sets.
#'
#' @param gene - the gene to annotate
#' @param geneSets - a list of gene sets
#' @param normalise - whether to normalise by set size
#' @return a numeric
#' @export
#' @examples
#' hallmark = make.getExample()('hallmark')
#' Myc = annotateGeneAsVector('Myc',hallmark)
#' MycNormalised = annotateGeneAsVector('Myc',hallmark,TRUE)
annotateGeneAsVector = function(gene,geneSets,normalise=FALSE)
{
    N = length(geneSets)
    annotation = rep(0,N)
    names(annotation) = names(geneSets)
    for(set in names(geneSets))
    {
        if(! gene %in% geneSets[[set]])
            next
        if(normalise)
            annotation[set] = 1 / length(geneSets[[set]])
        else
            annotation[set] = 1
    }
    return(annotation)
}


## ####################################################
#' This function makes annotation predictions for a
#' set of genes based on gene sets (e.g., hallmark)
#' and a CatsCradle object by considering the annotations of
#' its neighboring genes.
#'
#' @param genes - a character vector of genes
#' @param geneSets - a set of annotations, e.g., hallmark
#' or GO
#' @param fPrime - a Seurat object of genes SingleCellExperiment 
#' to be converted to a Seurat object
#' @param radius - radius for prediction neighborhood
#' @param metric - reduction or NN, defaults to umap
#' @param numPCs - used only if reduction is pca, defaults to NULL
#' @param normaliseByGeneSet - determines whether vector annotations
#'     are normalised by gene set size.  Defaults to TRUE
#' @param normaliseByDistance - determines whether neighbor
#'     contributions are normalised by edge weight.  Defaults to
#'     TRUE.
#' @param normaliseToUnitVector - determines whether to normalise
#'     returned values to unit length.  Defaults to TRUE
#' @return This returns a list of prediction vectors, one vector
#' for each gene in genes, each vector corresponding to the sets
#' in geneSets
#' @export
#' @examples
#' getExample = make.getExample()
#' STranspose = getExample('STranspose',toy=TRUE)
#' STranspose_sce = getExample('STranspose_sce',toy=TRUE)
#' hallmark = getExample('hallmark',toy=TRUE)
#' set.seed(100)
#' genes = sample(colnames(STranspose),5)
#' predictions = predictAnnotation(genes,hallmark,STranspose,radius=.5)
#' predictions_sce = predictAnnotation(genes,hallmark,STranspose_sce,radius=.5)
predictAnnotation = function(genes,
                             geneSets,
                             fPrime,
                             radius,
                             metric='umap',
                             numPCs=NULL,
                             normaliseByGeneSet=TRUE,
                             normaliseByDistance=TRUE,
                             normaliseToUnitVector=TRUE)
{
    fPrime = acceptor(fPrime)
    
    genesAnno = annotateGenesByGeneSet(geneSets)

    predictions = list()
    for(gene in genes)
    {
        ## This returns only the gene sets that show up:
        p = predictGeneAnnotationImpl(gene,fPrime,genesAnno,
                                      radius,metric,numPCs,
                                      normaliseByDistance)

        ## The 'wrapping' here is to produce a vector using
        ## all the gene sets:
        wrapped = rep(0,length(geneSets))
        names(wrapped) = names(geneSets)

        for(n in names(p))
        {
            if(normaliseByGeneSet)
                wrapped[n] = p[n] / length(geneSets[[n]])
            else
                wrapped[n] = p[n]
        }

        if(normaliseToUnitVector & sum(wrapped) != 0)
            wrapped = wrapped / Norm(wrapped)

        predictions[[gene]] = wrapped
    }

    return(predictions)
}


## ####################################################
#' This function is the implementation for predicting
#' the functions of a gene based on the functions of its
#' neighbours.
#'
#' @param gene - gene to annotate
#' @param fPrime - a Seurat object of genes or
#' SingleCellExperiment to be converted to a Seurat object
#' @param genesAnno - genes annotated with gene sets
#' @param radius - radius of neighbours to consider
#' @param metric - which metric to use to discover
#' neighbours, can be one of 'umap', 'tsne', 'pca', 'NN',
#' defaults to umap
#' @param numPCs - used only if metric is pca. Defaults to NULL
#' @param normaliseByDistance - choose whether to normalise
#' contributions of neighbors by their distance, defaults to
#' TRUE
#' @return This returns a named list.  The names are
#' the anotations that apply to the neighbour genes, the values
#' are the relative wieghts of the contributions.
#' @export 
#' @examples
#' getExample = make.getExample()
#' STranspose = getExample('STranspose',toy=TRUE)
#' hallmark = getExample('hallmark',toy=TRUE)
#' genesAnno = annotateGenesByGeneSet(hallmark)
#' predictions = predictGeneAnnotationImpl('Myc',STranspose,genesAnno,
#' radius=.5,metric='umap')
predictGeneAnnotationImpl = function(gene,fPrime,genesAnno,
                                     radius,metric,numPCs=NULL,
                                     normaliseByDistance=TRUE)
{
    fPrime = acceptor(fPrime)
    
    predictedTerms = list()

    ## Use weights only in this case:
    weights = (metric == 'NN' & normaliseByDistance==TRUE)
    nearby = getNearbyGenes(fPrime,gene,radius,metric,numPCs,weights)

    ## Neighbour genes are names.
    ## Distances are values.
    for(neighbour in names(nearby))
    {
        if (!(neighbour %in% names(genesAnno)))
            next
        
        ## extract terms for neighbour
        terms = genesAnno[[neighbour]]

        ## add these to predicted terms with appropriate weights
        for (term in terms)
        {
            if(term %in% names(predictedTerms))
                predictedTerms[[term]] =
                    predictedTerms[[term]] + 1 / nearby[neighbour]
            else
                predictedTerms[[term]] = 1 / nearby[neighbour]
        }
    }
    termNames = names(predictedTerms)
    predictedTerms = as.numeric(predictedTerms)
    names(predictedTerms) = termNames
    predictedTerms =  
        predictedTerms[order(-predictedTerms)]
    
    return(predictedTerms)
}


## ####################################################
#' This function predicts the functions of all genes based
#' on the functions of their neighbours.
#'
#' @param geneSets - a set of gene sets, e.g., hallmark
#' @param fPrime - a transposed Seurat object (generated with 
#' transposeObject()) or SingleCellExperiment to
#' be converted to a Seurat object
#' @param radius - radius of the region to use for prediction
#' @param metric - reduction or NN, defaults to umap
#' @param normaliseByGeneSet - normalise by size of each gene set,
#' defaults to TRUE
#' @param normaliseByDistance - attenutate neighbour contributions
#' based on distance, defaults to TRUE
#' @param normaliseToUnitVector - return results as unit
#' vectors, defaults to TRUE
#' @return - A list where names are genes and values are vectors
#' of gene annotations whose entries correspond to the geneSets
#' @export
#' @examples
#' getExample = make.getExample()
#' STranspose = getExample('STranspose',toy=TRUE)
#' hallmark = getExample('hallmark',toy=TRUE)
#' predictions = predictAnnotationAllGenes(hallmark,STranspose,radius=.5)
predictAnnotationAllGenes = function(geneSets,
                                     fPrime,
                                     radius,
                                     metric='umap',
                                     normaliseByGeneSet=TRUE,
                                     normaliseByDistance=TRUE,
                                     normaliseToUnitVector=TRUE)
{

    fPrime = acceptor(fPrime)
    genes = colnames(fPrime)
    
    return(predictAnnotation(genes,
                             geneSets,
                             fPrime,
                             radius,
                             metric,
                             normaliseByGeneSet,
                             normaliseByDistance,
                             normaliseToUnitVector))
}


## ####################################################
#' Nearby genes
#'
#' This finds the genes near a give subset using either
#' a dimensional reduction or the nearest neighbor graph
#'
#' @param fPrime - a Seurat object of genes or
#' SingleCellExperiment to be converted to a Seurat object
#' @param geneSet - set of genes
#' @param metric - the metric to use, one of umap, tsne,
#' pca or nearest neighbor
#' @param radius - the distance around the given set
#' @param numPCs - used only if the metric is pca
#' @param weights - whether to use edge weights in the NN case
#' @return This returns a named vector whose values are distance
#' from geneSet and whose names are the nearby genes.
#' @export
#' @examples
#' getExample = make.getExample()
#' STranspose = getExample('STranspose',toy=TRUE)
#' hallmark = getExample('hallmark',toy=TRUE)
#' geneSet = intersect(colnames(STranspose),hallmark[[1]])
#' geometricallyNearby = getNearbyGenes(STranspose,geneSet,radius=0.2,metric='umap')
#' combinatoriallyNearby = getNearbyGenes(STranspose,geneSet,radius=1,metric='NN')
#' weightedNearby = getNearbyGenes(STranspose,'Myc',radius=1,metric='NN',weights=TRUE)
getNearbyGenes = function(fPrime,geneSet,radius,metric='umap',
                       numPCs=NULL,weights=FALSE)
{
    stopifnot(metric %in% c('umap','tsne','pca','NN'))

    fPrime = acceptor(fPrime)

    ## ####################################################
    ## Combinatorial with weights:
    if(weights)
    {
        stopifnot(length(geneSet) == 1 &
                  metric == 'NN' &
                  radius == 1)
        return(nearbyGenesWeighted(fPrime,geneSet))
    }

    ## ####################################################
    ## Combinatorial:
    if(metric == 'NN')
    {
        NN = getNearestNeighbourLists(fPrime)
        spheres = combinatorialSpheres(NN,geneSet,radius)
        idx = ! spheres$nodes %in% geneSet
        spheres = spheres[idx,]

        
        nearby = spheres$radius
        names(nearby) = spheres$nodes
        
        return(nearby)
    }

    ## ####################################################
    ## Metric:
    ## Get the coords:
    if(metric == 'umap')
        a = fetchUMAP(fPrime)
    if(metric == 'tsne')
        a = FetchData(fPrime,c('tSNE_1','tSNE_2'))
    if(metric == 'pca')
    {
        pcs = paste0('PC_',seq_len(numPCs))
        a = FetchData(fPrime,pcs)
    }
    S = data.matrix(a)
    rownames(S) = colnames(fPrime)

    idx = rownames(S) %in% geneSet
    X = S[idx,]
    SHat = S[!idx,]

    D = distmat(SHat,X)
    rowMin = c()
    for(i in seq_len(nrow(D)))
        rowMin[i] = min(D[i,])
    idx = rowMin <= radius

    nearby = rowMin[idx]
    names(nearby) = rownames(SHat)[idx]

    return(nearby)
}

## ####################################################
nearbyGenesWeighted = function(fPrime,gene)
{
    NN = getNearestNeighbourLists(fPrime)
    NN = symmetriseNN(NN)
    idx = NN$nodeA == gene
    neighbours = NN$nodeB[idx]
    distance = 1 / NN$weight
    names(distance) = neighbours

    return(distance)
}


## ####################################################
#' This orders the gene set p-values (or -log10 p-values) and
#' applies a cutoff (if given) to show only the significant
#' gene sets for each gene cluster
#'
#' @param M - A matrix of gene set p-values (or their logs)
#' to be ordered by their significance
#' @param ascending - Direction in which to order the columns.
#' Defaults to TRUE, so that p-values will be ordered according
#' to decreasing significance, should be set to FALSE if ordering
#' -log p-value
#' @param cutoff - if non-null this is used to extract only
#' significant cases
#' @param nameTag - can be used to modify the names of the list.
#' @return This returns a list of whose entries are data frames,
#' one for each gene cluster, each giving the significant gene sets
#' for that cluster and their significance.
#' @export
orderGeneSetPValues = function(M,ascending=TRUE,cutoff=NULL,nameTag='')
{
    answer = list()
    
    for(i in seq_len(ncol(M)))
    {
        tag = paste0(nameTag,colnames(M)[i])

        ## Get the two column df:
        a = data.frame(geneSet=rownames(M),
                       value=M[,i],
                       stringsAsFactors=FALSE)
        rownames(a) = NULL

        ## Re-order the two-col df:
        if(ascending)
            a = a[order(a$value),]
        else
            a = a[order(-a$value),]
        
        ## Apply the cutoff:
        if(!is.null(cutoff))
        {
            if(ascending)
                idx = a$value <= cutoff
            else
                idx = a$value >= cutoff

            a = a[idx,]
        }
        answer[[tag]] = a
    }
    return(answer)
}

## ####################################################
#' This gets the clusters in their cannonical order
#'
#' This deals with skullduggery in which seurat_clusters
#' has been converted from a factor to a character or a
#' numeric. 
#'
#' @param f - a Seurat object with meta.data column seurat_clusters
#'  or SingleCellExperiment to be turned into a Seurat object
#' @return A vector of these unique values in order
#' @export
#' @examples
#' STranspose = make.getExample()('STranspose',toy=TRUE)
#' geneClusters = getClusterOrder(STranspose)
getClusterOrder = function(f)
{
    f = acceptor(f)
    a = unique(f$seurat_clusters)
    a = a[order(as.numeric(a))]

    return(a)
}

## ####################################################
#' This function reads in gene sets in .gmt format
#'
#' @param gmtFile - a .gmt file containing gene sets, e.g., Hallmark of GO
#' @param addDescr - include gene set description (2nd column in .gmt file) in 
#' gene set name  
#' @return - A named list of gene sets
#' @export
readGmt = function(gmtFile, addDescr = FALSE){
    lines = readLines(gmtFile)
    geneSets = list()
    for (line in lines){
        info = strsplit(line, "\t")[[1]]
        if (addDescr){
            name = paste(info[1],info[2])
        } else {
            name = info[1]
        }
        geneSets[name] = list(info[seq(from=3,
                                       to=length(info))])
    } 
    return(geneSets)
}

## ####################################################
#' This function strips out non-gene information from
#' the beginning of GO sets, etc.
#'
#' @param geneSet - a list of gene sets
#' @return a named list of gene sets
#' @export
stripGeneSet = function(geneSet)
{
    names = character(length(geneSet))
    for(i in seq_len(length(geneSet)))
    {
        names[i] = geneSet[[i]][1]
        geneSet[[i]] = geneSet[[i]][seq(from=3,to=length(geneSet[[i]]))]
    }
    
    return(geneSet)
}
