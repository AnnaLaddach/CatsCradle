

## ####################################################
## These functions study how expression of the gene
## clusters vary across cells and cell clusters.
## ####################################################



## ###################################################
#' This computes average expression of each gene cluster in
#' each cell cluster and returns the result as a matrix
#'
#' @param f - The Seurat object of cells, or SingleCellExperiment
#' to be turned into a Seurat object
#' @param fPrime - The Seurat object of genes,  or SingleCellExperiment
#' to be turned into a Seurat object
#' @param clusteringName In many cases, this will be the cell
#'     clustering, i.e., seurat_clusters, which is the default, but
#'     for neighbourhood Seurat objects, this can be
#'     neighbourhood_clusters.
#' @param layer - layer to use for expression values
#' @return A matrix of the average expression where the rows
#' correspond to cell clusters and the columns correspond to
#' gene clusters.
#' @export
#' @importFrom stringr str_replace str_replace_all
#' @importFrom stringr str_split str_split_fixed
#' @examples
#' getExample = make.getExample()
#' STranspose = getExample('STranspose',toy=TRUE)
#' exSeuratObj = getExample('exSeuratObj',toy=TRUE)
#' M = getAverageExpressionMatrix(exSeuratObj,STranspose,layer='data')
getAverageExpressionMatrix = function(f,fPrime,
                                      clusteringName='seurat_clusters',
                                      layer='scale.data')
{
    f = acceptor(f)
    fPrime = acceptor(fPrime)
    
    f$clustering = as.character(f@meta.data[,clusteringName])
    cellCluster = unique(f$clustering)
    cellCluster = cellCluster[order(as.numeric(cellCluster))]
    
    fPrime$clustering = as.character(fPrime@meta.data[,clusteringName])
    geneCluster = unique(fPrime$clustering)
    geneCluster = geneCluster[order(as.numeric(geneCluster))]
    
    ## Get assay data:
    X = GetAssayData(f,layer=layer)
    
    ## Seems X can be smaller:
    selected = rownames(X)
    f = f[selected,]
    selectedUnderscore = str_replace(selected,'-','_')
    if (sum(selectedUnderscore %in% colnames(fPrime)) > 
        sum(selected %in% colnames(fPrime))){
      fPrime = fPrime[,selectedUnderscore]
    } else {
      fPrime = fPrime[,selected]
    }
   
    M = matrix(0,nrow=length(cellCluster),ncol=length(geneCluster))
    rownames(M) = cellCluster
    colnames(M) = geneCluster
  
    for(i in cellCluster)
    {
        for(j in geneCluster)
        {
            idxI = f$clustering == i
            idxJ = fPrime$clustering == j
            
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
#' @return The same matrix with fancier row and col names
#' @export
#' @examples
#' getExample = make.getExample()
#' averageExpMatrix = getExample('averageExpMatrix',toy=TRUE)
#' averageExpMatrix = tagRowAndColNames(averageExpMatrix,'cellCluster_','geneCluster_')
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
#' @importFrom reshape2 melt
#' @export
#' @examples
#' getExample = make.getExample()
#' averageExpMatrix = getExample('averageExpMatrix',toy=TRUE)
#' averageExpDF = getAverageExpressionDF(averageExpMatrix)
getAverageExpressionDF = function(M)
{
    df = melt(M)
    names(df) = c('cellCluster','geneCluster','averageExpression')

    return(df)
}


## ####################################################
#' This makes a sankey graph from a matrix of average
#' expression.  Our "Cat's Cradle".
#'
#' @param M - a matrix of gene expression
#' @param disambiguation - used to distinguish between
#' the row names and the column names if these overlap
#' @param fontSize - defaults to 20
#' @param minus - colour to use for links with negative
#' values
#' @param plus - colour for positive values
#' @param height - height in pixels, defaults to 1200
#' @param width - width in pixels, defaults to 900
#' @return A sankey graph
#' @export
#' @importFrom networkD3 sankeyNetwork
#' @examples
#' set.seed(100)
#' M = matrix(runif(12)-.3,nrow=3)
#' rownames(M) = as.character(seq_len(3))
#' colnames(M) = as.character(seq_len(4))
#' sankey = sankeyFromMatrix(M)
sankeyFromMatrix = function(M,disambiguation=c('R_','C_'),
                            fontSize=20,minus='red',plus='blue',
                            height=1200,width=900)
{
  ## Maybe the matrix doesn't have row and column names:
  if(is.null(rownames(M)))
    rownames(M) = paste0(disambiguation[1],seq_len(nrow(M)))
  if(is.null(colnames(M)))
    colnames(M) = paste0(disambiguation[1],seq_len(ncol(M)))
  
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
  
  ## Colour the links DF:
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
  
  linkColour = 'd3.scaleOrdinal() .domain(["minus","plus"]) .range(["X", "Y"])'
  linkColour = str_replace(linkColour,'X',minus)
  linkColour = str_replace(linkColour,'Y',plus)  
  
  
    p = sankeyNetwork(Links = links, Nodes = nodes,
                      Source = "IDsource", Target = "IDtarget", 
                      Value = "value", NodeID = "name", 
                      colourScale=linkColour,
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
#' @param f - the cell Seurat object or SingleCellExperiment to be
#' turned into a Seurat object
#' @param fPrime - the genes Seurat object or SingleCellExperiment to be
#' turned into a Seurat object
#' @param cells - the cells to compute this for
#' @param geneClusters - the geneClusters to compute average
#' expression for
#' @param layer - the data layer to use, defaults to 'data'
#' @return A matrix where the rows correspond to cells, the columns
#' correspond to geneClusters and the entries give average expression
#' for each cluster in each cell
#' @export
#' @examples
#' getExample = make.getExample()
#' exSeuratObj = getExample('exSeuratObj',toy=TRUE)
#' STranspose = getExample('STranspose',toy=TRUE)
#' clusterExpression = getGeneClusterAveragesPerCell(exSeuratObj,STranspose)
getGeneClusterAveragesPerCell = function(f,
                                         fPrime,
                                         cells=colnames(f),
                                         geneClusters=getClusterOrder(fPrime),
                                         layer='data')
{
    f = acceptor(f)
    fPrime = acceptor(fPrime)
    
    M = matrix(0,nrow=length(cells),ncol=length(geneClusters))
    rownames(M) = cells
    colnames(M) = paste0('cluster_',geneClusters)

    for(i in seq_len(length(geneClusters)))
    {
        cluster = geneClusters[i]
        idx = fPrime$seurat_clusters == cluster
        theseGenes = colnames(fPrime)[idx]
        expression = data.matrix(FetchData(f,theseGenes,layer=layer))
        expression = expression[rownames(expression) %in% cells,]

        M[,i] = rowMeans(expression)
    }

    return(M)
}


## ####################################################
#' Mean gene cluster on cell umap
#'
#' This function paints gene expression for a
#' given gene cluster on cell umap.
#'
#' @param f - a Seurat object of cells or SingleCellExperiment to
#' be converted to a Seurat object
#' @param fPrime - the corresponding Seurat object of genes
#' SingleCellExperiment to be converted to a Seurat object
#' @param geneCluster - a gene cluster of fPrime
#' @return This returns a ggplot object
#' @importFrom ggplot2 ggplot geom_point geom_density theme xlim ggtitle
#' @importFrom ggplot2 geom_vline facet_wrap aes
#' @export
#' @examples
#' getExample = make.getExample()
#' exSeuratObj = getExample('exSeuratObj',toy=TRUE)
#' STranspose = getExample('STranspose',toy=TRUE)
#' g = meanGeneClusterOnCellUMAP(exSeuratObj,STranspose,geneCluster=0)
meanGeneClusterOnCellUMAP = function(f,fPrime,geneCluster)
{
    f = acceptor(f)
    fPrime = acceptor(fPrime)
    
    plotDF = fetchUMAP(f)
    idx = fPrime$seurat_clusters == geneCluster
    genes = colnames(fPrime)[idx]
    genes = intersect(genes,rownames(f))
    expression = FetchData(f,genes)
    plotDF$expression = rowMeans(expression)

    g = ggplot(plotDF,aes(x=UMAP_1,y=UMAP_2,color=expression)) +
        geom_point()

    return(g)
}

## ####################################################
#' This gets z-scores for the values of features
#'
#' @param f - a Seurat object of cells or SingleCellExperiment to
#' be converted to a Seurat object
#' @param featurs - a set of features to retrieve z-scores for,
#' defaults to rownames(f)
#' @param layer - the data layer to retrieve
#' @return This returns a data frame with a column for each
#' feature and a row for each cell
#' @export
#' @examples
#' getExample = make.getExample()
#' exSeuratObj = getExample('exSeuratObj',toy=TRUE)
#' df = getFeatureZScores(exSeuratObj)
getFeatureZScores = function(f,features=rownames(f),layer='data')
{
    f = acceptor(f)
    
    df = FetchData(f,features,layer=layer)
    
    for(i in 1:ncol(df))
        df[,i] = (df[,i] - mean(df[,i])) / sd(df[,i])

    return(df)
}

## ####################################################
#' This finds the mean z-score for features in subsets
#' of cells e.g., in each of the seurat_clusters
#'
#' @param f - a Seurat object of cells or SingleCellExperiment to
#' be converted to a Seurat object
#' @param features - a set of features of f
#' @param clusterBy - the name of the column of f@meta.data to
#' be used to subset the cells
#' @param layer - the data layer to be used for z-scores
#' @return This returns a data frame each of whose columns
#' corresponds to a value of the clusterBy data.  In the case
#' where the clusterBy data is a factor or numeric, it prepends
#' cluster_ to the column name.
#' @export
#' @examples
#' getExample = make.getExample()
#' exSeuratObj = getExample('exSeuratObj',toy=TRUE)
#' STranspose = getExample('STranspose',toy=TRUE)
#' df = meanZPerCluster(exSeuratObj,features=colnames(STranspose),
#'                      clusterBy='shortName')
meanZPerCluster = function(f,
                           features,
                           clusterBy='seurat_clusters',
                           layer='data')
{
    f = acceptor(f)
    
    clusters = unique(f@meta.data[,clusterBy])
    clusters = clusters[order(clusters)]
    if(isa(clusters,'numeric') |
       isa(clusters,'factor'))
    {
        Clusters = paste0('cluster_',clusters)
    } else {
        Clusters = as.character(clusters)
    }

    z = getFeatureZScores(f,features=features,layer=layer)

    df = list()
    for(i in 1:length(clusters))
    {
        cluster = clusters[i]
        Cluster = Clusters[i]

        cells = colnames(f)[f@meta.data[,clusterBy] == cluster]
        meanZ = colMeans(z[cells,])

        df[[Cluster]] = meanZ
    }
    df = as.data.frame(df)
    rownames(df) = features

    return(df)
}

## ####################################################
#' This collects together mean z-score data together with
#' UMAP coordinates from the gene seurat object for plotting.
#'
#' @param f - a Seurat object of cells or SingleCellExperiment to
#' be converted to a Seurat object
#' @param fPrime - the corresponding Seurat object of genes
#' SingleCellExperiment to be converted to a Seurat object
#' @param clusterBy - the name of the column of f@meta.data to
#' be used to subset the cells
#' @param layer - the data layer to be used for z-scores
#' @return This returns a data frame with the UMAP coordinates
#' of the gene Seurat object and the average z-score for each
#' gene within each of the cell clusters defined by the clusterBy
#' column of the meta.data of f.
#' @export
#' @examples
#' getExample = make.getExample()
#' exSeuratObj = getExample('exSeuratObj',toy=TRUE)
#' STranspose = getExample('STranspose',toy=TRUE)
#' df = meanZPerClusterOnUMAP(exSeuratObj,STranspose,clusterBy='shortName')
meanZPerClusterOnUMAP = function(f,
                                 fPrime,
                                 clusterBy='seurat_clusters',
                                 layer='data')
{
    f = acceptor(f)
    fPrime = acceptor(fPrime)
    
    umap = FetchData(fPrime,c('umap_1','umap_2'))
    umap = cbind(gene=colnames(fPrime),umap)
    rownames(umap) = umap$gene

    meanZ = meanZPerCluster(f,
                            features=colnames(fPrime),
                            clusterBy=clusterBy,
                            layer=layer)
    umap = cbind(umap,meanZ)

    return(umap)
}


