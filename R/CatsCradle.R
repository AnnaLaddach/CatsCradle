


## ###################################################
#' Create the transpose of a Seurat object
#'
#' This takes a Seurat object f and creates a new Seurat object whose
#' expression matrix is the transpose of that of f.
#'
#' @param f - a Seurat object
#' @param active.assay - the assay to use.  Defaults to the
#'     active.assay
#' @param npcs - number of principal components, defaults to 30
#' @param dims - dimensions to use for umap and nearest neighbors,
#'     defaults to 1:20
#' @param res - the clustering resolution, defaults to 1
#' @return A Seurat object
#' @export
#' @import Seurat
#' @examples
#' STranspose = transposeSeuratObject(S)
transposeSeuratObject = function(f,active.assay=f@active.assay,
                                 npcs=30,dims=1:20,res=1)
{
  f@active.assay = active.assay
  M = f@assays[[active.assay]]@data
  M = as.matrix(M)
  MPrime = t(M)
  rownames(MPrime) = colnames(f)
  colnames(MPrime) = rownames(f)
  
  fPrime = CreateSeuratObject(MPrime,assay=active.assay)
  
  fPrime = FindVariableFeatures(fPrime)
  fPrime = ScaleData(fPrime)
  fPrime = RunPCA(fPrime,npcs=npcs)
  fPrime = RunUMAP(fPrime,reduction='pca',dims=dims)
  fPrime = FindNeighbors(fPrime)
  fPrime = FindClusters(fPrime,resolution=res)
  
  return(fPrime)
}

## ###################################################
#' This computes average expression of each gene cluster in
#' each cell cluster and returns the result as a matrix
#'
#' @param f - The Seurat object of cells
#' @param fPrime - The Seurat object of genes
#' @return A matrix of the average expression
#' @export
#' @examples
#' M = getAverageExpressionMatrix(S,STranspose)
getAverageExpressionMatrix = function(f,fPrime)
{
  f$seurat_clusters = as.character(f$seurat_clusters)
  cellCluster = unique(f$seurat_clusters)
  cellCluster = cellCluster[order(as.numeric(cellCluster))]
  
  fPrime$seurat_clusters = as.character(fPrime$seurat_clusters)
  geneCluster = unique(fPrime$seurat_clusters)
  geneCluster = geneCluster[order(as.numeric(geneCluster))]
  
  ## Get assay data:
  X = GetAssayData(f,slot='scale')

  ## Seems X can be smaller:
  f = f[rownames(X),]
  fPrime = fPrime[,rownames(X)]
   
  M = matrix(0,nrow=length(cellCluster),ncol=length(geneCluster))
  rownames(M) = cellCluster
  colnames(M) = geneCluster
  
  for(i in cellCluster)
  {
    for(j in geneCluster)
    {
      idxI = f$seurat_clusters == i
      idxJ = fPrime$seurat_clusters == j
      
      M[i,j] = sum(X[idxJ,idxI]) / (sum(idxI) * sum(idxJ))
    }
  }
  return(M)
}

## ####################################################
#' This gussies up the rownames and colnames of M
#'
#' @param M - a matrix, typically the average expression matrix
#' @param ccTag - a prefix for the row (cell cluster) names
#' @param gcTag - a prefix for the column (gene cluster) names
#' @export
#' @examples
#' averageExpMatrix = tagRowAndColNames(averageExpMatrix)
tagRowAndColNames = function(M,ccTag='CC_',gcTag='GC_')
{
  rownames(M) = paste0(ccTag,rownames(M))
  colnames(M) = paste0(gcTag,colnames(M))
  
  return(M)
}

## ####################################################
#' This converts an average gene expression matrix to a
#' data frame.
#'
#' @param M - An average gene expression matrix.
#' @return A data frame with columns cellCluster, geneCluster
#' and average expression
#' @import reshape2
#' @export
#' @examples
#' averageExpDF = getAverageExpressionDF(averageExpMatrix)
getAverageExpressionDF = function(M)
{
    df = melt(M)
    names(df) = c('cellCluster','geneCluster','averageExpression')

    return(df)
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
  for(i in 1:length(geneSet))
  {
    names[i] = geneSet[[i]][1]
    geneSet[[i]] = geneSet[[i]][3:length(geneSet[[i]])]
  }
  
  return(geneSet)
}

## ####################################################
## Used internally:
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
#' @return a matrix of p-values
#' @export
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
    
    for(i in 1:NGeneSets)
    {
        for(j in 1:NClusters)
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
    
    for(i in 1:ncol(M))
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
#' This makes a sankey graph from a matrix of average
#' expression.  Our "Cat's Cradle".
#'
#' @param M - a matrix of gene expression
#' @param disambiguation - used to distinguish between
#' the row names and the column names if these overlap
#' @param fontSize - defaults to 20
#' @param minus - color to use for links with negative
#' values
#' @param plus - color for positive values
#' @param height - height in pixels, defaults to 1200
#' @param width - width in pixels, defaults to 900
#' @return A sankey graph
#' @export
#' @import networkD3
#' @import stringr
#' @examples
#' set.seed(100)
#' M = matrix(runif(12)-.3,nrow=3)
#' rownames(M) = as.character(1:3)
#' colnames(M) = as.character(1:4)
#' S = sankeyFromMatrix(M)
sankeyFromMatrix = function(M,disambiguation=c('R_','C_'),
                            fontSize=20,minus='red',plus='blue',
                            height=1200,width=900)
{
  ## Maybe the matrix doesn't have row and column names:
  if(is.null(rownames(M)))
    rownames(M) = paste0(disambiguation[1],1:nrow(M))
  if(is.null(colnames(M)))
    colnames(M) = paste0(disambiguation[1],1:ncol(M))
  
  ## Create the links DF:
  from = rownames(M)
  to = colnames(M)
  
  if(length(intersect(from,to)) > 0)
  {
    from = paste0(disambiguation[1],from)
    to = paste0(disambiguation[2],to)
  }
  
  source = rep(from,each=length(to))
  target = rep(to,length(from))
  
  value = as.numeric(t(M))
  
  links = data.frame(source,target,value,
                     stringsAsFactors=FALSE)
  
  ## Color the links DF:
  idx = links$value >= 0
  links$group = ''
  links$group[idx] = 'plus'
  links$group[!idx] = 'minus'
  
  links$value = abs(links$value)
  
  ## Create the nodes DF:
  nodes = unique(c(links$source,links$target))
  nodes = data.frame(name=nodes,
                     stringsAsFactors=FALSE)
  
  links$IDsource = match(links$source,nodes$name) - 1
  links$IDtarget = match(links$target,nodes$name) - 1
  
  linkColor = 'd3.scaleOrdinal() .domain(["minus","plus"]) .range(["X", "Y"])'
  linkColor = str_replace(linkColor,'X',minus)
  linkColor = str_replace(linkColor,'Y',plus)  
  
  
  p = sankeyNetwork(Links = links, Nodes = nodes, Source = "IDsource", Target = "IDtarget", 
                    Value = "value", NodeID = "name", 
                    colourScale=linkColor,
                    LinkGroup="group",
                    fontSize=fontSize,
                    height=height,
                    width=width)
  
  return(p)
}

## ####################################################
#' This produces a matrix giving the average expression of gene
#' clusters in cells.  By default, it uses all cells and all gene
#' clusters.
#'
#' @param f - the cell Seurat object
#' @param fPrime - the genes Seurat object
#' @param cells - the cells to compute this for
#' @param geneClusters - the geneClusters to compute average
#' expression for
#' @return A matrix where the rows correspond to cells, the columns
#' correspond to geneClusters and the entries give average expression
#' for each cluster in each cell
#' @export
#' @examples
#' clusterExpression = getGeneClusterAveragesPerCell(S,STranspose)
getGeneClusterAveragesPerCell = function(f,
                                         fPrime,
                                         cells=colnames(f),
                                         geneClusters=getClusterOrder(fPrime))
{
    M = matrix(0,nrow=length(cells),ncol=length(geneClusters))
    rownames(M) = cells
    colnames(M) = paste0('cluster_',geneClusters)

    for(i in 1:length(geneClusters))
    {
        cluster = geneClusters[i]
        idx = fPrime$seurat_clusters == cluster
        theseGenes = colnames(fPrime)[idx]
        expression = data.matrix(FetchData(f,theseGenes))
        expression = expression[rownames(expression) %in% cells,]

        M[,i] = rowMeans(expression)
    }

    return(M)
}

## ####################################################
#' This gets the clusters in their cannonical order
#'
#' This deals with skullduggery in which seurat_clusters
#' has been converted to a character or a numeric. 
#'
#' @param f - a Seurat object with meta.data column seurat_clusters
#' @return A vector of these unique values in order
#' @export
getClusterOrder = function(f)
{
    a = unique(f$seurat_clusters)
    a = a[order(as.numeric(a))]

    return(a)
}

## ####################################################
## Default graph
## Used internally
defaultGraph = function(f)
{
    graph = paste0(f@active.assay,'_snn')
    return(graph)
}

## ####################################################
#' This function extracts a shared nearest neighbor network
#' from a Seurat object
#'
#' @param f - a Seurat object
#' @param graph - which graph to extract.  Defaults to
#' paste0(f@active.assay,'_snn')
#' @return - This returns dataframe of neighbors:
#' nodeA - node names for node A 
#' nodeB - node names for node B
#' weight - edge weight
#' @export
#' @examples
#' NN = getNearestNeighborListsSeurat(STranspose)
getNearestNeighborListsSeurat = function(f, graph=defaultGraph(f)){
  
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
#' This function gets the neighbors of a given gene
#' using either the gene Seurat object or its nearest
#' neighbor graph returned from getNearestNeighborListsSeurat
#'
#' @param gene - the gene in question
#' @param NN - either the gene Seurat object or its nearest
#' neighbor graph as found by getNearestNeighborListsSeurat
#' @return the neighboring genes
#' @export
#' @examples
#' neighbors = getGeneNeighbors("Ccl6",STranspose)
#' neighborsAgain = getGeneNeighbors("Ccl6",NN)
getGeneNeighbors = function(gene,NN)
{
    if(class(NN) == 'Seurat')
        NN = getNearestNeighborListsSeurat(NN)

    NN = symmetrise(NN)

    idx = NN$nodeA == gene

    return(NN$nodeB[idx])
}

## ####################################################
#' This function takes the data frame of neighbor genes
#' and reduces it so that each undirected edge is
#' represented by only one directed edge.  This ensures
#' that randomisation does not magically split undirected
#' edges into two edges.
#'
#' @param neighborListDf - a dataframe containing the neighborlist
#' @return - a neighborListDF with only one directed edge per
#' undirected edge.
#' @export
#' @examples
#' print(dim(NN))
#' NNN = desymmetrise(NN)
#' print(dim(NNN))
desymmetriseNN = function(NN)
{
    orderGenes = function(i)
    {
        if(order(c(NN$nodeA[i],NN$nodeB[i]))[1] == 2)
            return(TRUE)
        return(FALSE)
    }

    idx = unlist(lapply(1:nrow(NN),orderGenes))
    NN[idx,1:2] = NN[idx,2:1]
    
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
#' randomIndices = randomiseNodeIndices(NN,10,TRUE)
randomiseNodeIndices = function(neighborListDf, n = 100, useWeights = F){
    NN = desymmetrize(neighborListDf)
    if(!identical(NN,neighborListDf))
        stop(paste0('randomiseNodeIndices is meant to be used',
                    'with a desymmetrised neighborList'))
  
  #determine number of edges and create empty matrix for randomised indices
  nEdges = nrow(neighborListDf)
  indices = 1:nEdges
  randomIndices = matrix(, nrow = nEdges, ncol = n)
  
  #check if weights are to be used
  if (useWeights){
    
    #determine unique weights
    weights = unique(neighborListDf$weight)
    for (i in 1:n){
      
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
    for (i in 1:n){
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
#' of a data frame  as returned by getNearestNeighborListsSeurat.
#' Its coloumns include nodeA and nodeB.
#'
#' @return TRUE or FALSE
#' @export
#' @examples
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
#' by getNearestNeighborListsSeurat
#'
#' @return a nearest neighbors graph
#' 
#' @export
#' @examples
#' NNStar = symmetriseNN(NN)
symmetriseNN = function(NN)
{
    ## Is it already symmetric?
    if(symmetryCheckNN(NN))
        return(NN)

    ## If not, symmetrise:
    NN2 = NN[,c(2,1,3)]
    NN = rbind(NN,NN2)

    tag = paste(NN$nodeA,NN$nodeB)
    idx = ! duplicated(tag)

    NN = NN[idx,]

    return(NN)
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
#' spheres = combinatorialSpheres(NN,'Ccl6',3)
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
#' This function reads in gene sets in .gmt format
#'
#' @param gmtFile - a .gmt file containing gene sets, e.g., Hallmark of GO
#' @param addDescr - include gene set description (2nd column in .gmt file) in 
#' gene set name  
#' @return - A named list of gene sets
#' @export
readGmt = function(gmtFile, addDescr = F){
  lines = readLines(gmtFile)
  geneSets = list()
  for (line in lines){
    info = strsplit(line, "\t")[[1]]
    if (addDescr){
      name = paste(info[1],info[2])
    } else {
      name = info[1]
    }
    geneSets[name] = list(info[3:length(info)])
  } 
  return(geneSets)
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
#' annotatedGenes = annotateGenes(hallmark)
annotateGenes = function(geneSets)
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
#' @param fPrime - a Seurat object of genes
#' @param normaliseByGeneSet - determines whether vector annotations
#'     are normalised by gene set size.  Defaults to TRUE
#' @param normaliseByEdgeWeight - determines whether neighbor
#'     contributions are normalised by edge weight.  Defaults to
#'     TRUE.
#' @param normaliseToUnitVector - determines whether to normalise
#'     returned values to unit length.  Defaults to TRUE
#' @export
#' @examples
#' set.seed(100)
#' genes = sample(colnames(STranspose),5)
#' predictions = predictTerms(genes,hallmark,STranspose)
#' 
predictTerms = function(genes,
                        geneSets,
                        fPrime,
                        normaliseByGeneSet=TRUE,
                        normaliseByEdgeWeight=TRUE,
                        normaliseToUnitVector=TRUE)
{
    genesAnno = annotateGenes(geneSets)

    predictions = list()
    for(gene in genes)
    {
        p = predictTermsImpl(gene,fPrime,genesAnno,normaliseByEdgeWeight)

        wrapped = rep(0,length(geneSets))
        names(wrapped) = names(geneSets)

        for(n in names(p))
        {
            if(normaliseByGeneSet)
                wrapped[n] = p[n] / length(geneSets[[n]])
            else
                wrapped[n] = p[n]
        }

        if(normaliseToUnitVector)
            wrapped = wrapped / Norm(wrapped)

        predictions[[gene]] = wrapped
    }

    return(predictions)
}

## ####################################################
#' This function is the implementation for predicting the
#' functions of a gene based on the functions of its
#' neighbours.
#'
#' @param gene - gene to annotate
#' @param fPrime - a transposed Seurat object (generated with 
#' transposeSeuratObject())
#' @param genesAnno - genes annotated with gene sets
#' @param normaliseByEdgeWeight - choose whether to normalise 
#' (divide scores by total weight of edges)
#' @param graph - which graph to use.  Defaults to the
#' active.assay followed by _snn
#' @return - A lists of terms where values are scores.
#' Scores are calculated according to the weights of the SNN graph.
#'
#' @export
#' @examples
#' genesAnno = annotateGenes(hallmark)
#' predictions = predictTermsImpl('Myc',STranspose,genesAnno)
predictTermsImpl = function(gene, fPrime, genesAnno,
                            normaliseByEdgeWeight = TRUE,
                            graph=defaultGraph(fPrime)){
  
    predictedTerms = list()
    
    ## determine neighbors
    theGraph = fPrime@graphs[[graph]]
    genes = rownames(theGraph)
    neighbors = genes[theGraph[gene,] > 0]
    neighbors = neighbors[neighbors != gene]
  
    ## calculate total weight of edges to neighbors
    total = sum(fPrime@graphs[[graph]][gene,])
  
    ## iterate through neighbors
    for (neighbor in neighbors){
        
        ## determine weight of connecting edge and normalise if T
        weight = fPrime@graphs[[graph]][gene,neighbor]
        if (normaliseByEdgeWeight){
            weight = weight/total
        }
        if (!(neighbor %in% names(genesAnno))){
            next
        }
    
        ## extract terms for neighbor
        terms = genesAnno[[neighbor]]
    
        ## add these to predicted terms with appropriate weights
        for (term in terms){
            if (term %in% names(predictedTerms)){
                predictedTerms[[term]] = 
                    predictedTerms[[term]] + weight
            } else {
                predictedTerms[[term]] = weight 
            }
        }
        termNames = names(predictedTerms)
        predictedTerms = as.numeric(predictedTerms)
        names(predictedTerms) = termNames
        predictedTerms =  
            predictedTerms[order(-predictedTerms)]
    }
    return(predictedTerms)
}


## ####################################################
#' This function predicts the functions of all genes based
#' on the functions of their neighbours.
#'
#' @param fPrime - a transposed Seurat object (generated with 
#' transposeSeuratObject())
#' @param genesAnno - genes annotated with gene sets
#' @param genes - the genes to be predicted, defaults to all
#' @param normalise - whether to normalise to gene set size, defauls
#' to TRUE
#' @return - A list where names are genes and values are lists of terms.
#' The values of the lists of terms are calculated according to the weights
#' of the edges connecting the neighbors.
#' @export
#' @examples
#' set.seed(100)
#' genes = sample(colnames(STranspose),50)
#' genesAnno = annotateGenes(hallmark)
#' predictions = predictTermsAllGenes(STranspose,genesAnno,genes)
predictTermsAllGenes = function(fPrime,
                                genesAnno,
                                genes=colnames(fPrime),
                                normalise = T){
  
  genesPredictedTerms = list()
  i = 1
  
  #iterate through genes
  for (gene in genes){
    if (i %% 100 == 0){
        writeLines(paste(i, "genes processed"))
    }
    i = i + 1
    genesPredictedTerms[[gene]] = predictTerms(gene, fPrime, genesAnno, 
                                               normalise)
  }
  return(genesPredictedTerms)
}

## ####################################################
#' This function computes a p-value for the geometric
#' clustering of a gene set (in UMAP or PCA reduction)
#' based on the median distance from its complement to
#' the set.
#'
#' @param fPrime - a transposed Seurat object, i.e. a
#' Seurat object of genes
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
#' @import pracma
#' @export
#' @examples
#' geneSubset = intersect(colnames(STranspose),hallmark[[1]])
#' p = getSubsetClusteringPValue(STranspose,geneSubset,100)
getSubsetClusteringPValue = function(fPrime,
                                     geneSubset,
                                     numTrials=1000,
                                     reduction='UMAP',
                                     numPCs=10)
{
    stats = getSubsetClusteringStatistics(fPrime,
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
#' Seurat object of genes
#' @param geneSubset - a subset of the genes which can
#' be given as a character vector as a logical vector
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
#' geneSubset = intersect(colnames(STranspose),hallmark[[1]])
#' stats = getSubsetClusteringStatistics(STranspose,geneSubset,100)
getSubsetClusteringStatistics = function(fPrime,
                                         geneSubset,
                                         numTrials=1000,
                                         reduction='UMAP',
                                         numPCs=10)
{
    if(reduction == 'UMAP')
        S = FetchData(fPrime,c('UMAP_1','UMAP_2'))

    if(reduction == 'PCA')
    {
        pcs = paste0('PC_',1:numPCs)
        S = FetchData(fPrime,pcs)
    }

    ## The coordinates of all the genes in the
    ## chosen reduction:
    S = as.matrix(S)
    rownames(S) = colnames(fPrime)

    answer = runClusteringTrials(S,
                                 geneSubset,
                                 numTrials)

    return(answer)
}

## ####################################################
## Used internally
runClusteringTrials = function(S,
                               geneSubset,
                               numTrials)
{
    if(class(geneSubset) == 'character')
        geneSubset = rownames(S) %in% geneSubset
    
    answer = list()
    answer$subsetDistance = medianComplementDistance(S,geneSubset)

    randomSubsetDistance = c()
    for(i in 1:numTrials)
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
    if(class(geneSubset) == 'character')
        geneSubset = colnames(S) %in% geneSubset
    
    ## The complement:
    A = S[!geneSubset,]
    
    ## The subset:
    B = S[geneSubset,]
    
    D = distmat(A,B)
    
    rowMin = c()
    for(i in 1:nrow(D))
        rowMin[i] = min(D[i,])
    d = median(rowMin)
    
    return(d)
}
