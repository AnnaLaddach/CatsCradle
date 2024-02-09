## ####################################################
#' This function computes a spatial graph where 
#' neighbors are identified based on Delaunay triangulation.
#' 
#' @param centroids - a dataframe containing centroids 
#' where rownames are cellnames and columns contain x 
#' and y coordinates respectively.
#' @return a graph in neighbour format, i.e., a data frame with
#'     columns nodeA and nodeB.
#' @import interp
#' @export
## #' @examples
computeNeighboursDelaunay = function(centroids){
    
    ##compute delaunay triangulation
    delaunay = tri.mesh(centroids[,1], centroids[,2])
    
    ##extract triangles from results dataframe
    triangles = delaunay$trlist
    
    results = matrix(nrow = 0, ncol = 2)
    cellNames = rownames(centroids)
    
    ##convert triangle indices to cell names
    triangles[,1] = cellNames[as.numeric(triangles[,1])]
    triangles[,2] = cellNames[as.numeric(triangles[,2])]
    triangles[,3] = cellNames[as.numeric(triangles[,3])]
    
    ##create neighbour list
    for (i in 1:nrow(triangles)){
        results = rbind(results,sort(c(triangles[i,1],triangles[i,2])))
        results = rbind(results,sort(c(triangles[i,1],triangles[i,3])))
        results = rbind(results,sort(c(triangles[i,2],triangles[i,3])))
    }
    
    ##remove duplicate edges
    results = unique(results)
    
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
#' @import rdist
#' @export
## #' @examples
computeNeighboursEuclidean = function(centroids, threshold){
    colnames(centroids) = c("x","y")
    results = matrix(nrow = 0, ncol = 2)
    
    for (i in 1:nrow(centroids)){
        cellX = centroids[i,1]
        cellY = centroids[i,2]
        cellName = rownames(centroids)[i]
        
        ##only calculate distances for cells that could be within distance threshold.
        possibleNeighbours = centroids[(centroids[,"x"] < (cellX + threshold)) &
                                       (centroids[,"x"] > (cellX - threshold)) &
                                       (centroids[,"y"] < (cellY + threshold)) &
                                       (centroids[,"y"] > (cellY - threshold)),,
                                       drop = FALSE]
        
        neighbours = rownames(possibleNeighbours)[cdist(centroids[i,,drop=FALSE]
                                                       ,possibleNeighbours,)[1,] < threshold]
        
        for (neighbour in neighbours){
            if (cellName != neighbour){
                results = rbind(results,sort(c(cellName,neighbour)))
            }
        }
    }
    ##remove duplicate edges
    results = unique(results)
    
    ## Convert to data.frame and name as nodeA, nodeB:
    results = as.data.frame(results)
    names(results) = c('nodeA','nodeB')
    
    return(results)
}

## ####################################################
## Not exported:
extractCells = function(NN)
{
    cells = unique(c(NN$nodeA,NN$nodeB))
    cells = cells[order(cells)]
    return(cells)
}


## ####################################################
#' This function extracts neighbourhoods from a spatial graph in neighbour list
#' format
#' 
#' @param  spatialGraph - a spatial graph in neighbour list format.
#' @param cellNames - defaults to the names in the spatial graph
#' @param addSelf - logical determining whether to include the cell at
#' the center of the neighbourhood.  Defaults to FALSE
#' @return a named list of neighbourhoods where a neighborhood is a set of a 
#' cell's nearest neighbours. 
#' @export
## #' @examples
#' 
#' 
computeNeighbourhoods = function(spatialGraph,
                                 cellNames=extractCells(spatialGraph),
                                 addSelf = F){
  neighbourhoods = list()
  for (cell in cellNames){
    neighbours = 
      c(spatialGraph[,2][cell == spatialGraph[,1]],
        spatialGraph[,1][cell == spatialGraph[,2]])
    neighbourhoods[[cell]] = neighbours
    if (addSelf){
      neighbourhoods[[cell]] = c(neighbourhoods[[cell]], cell)
    }
  }
 return(neighbourhoods) 
}
  

## ####################################################
#' This function computes a matrix where neighbourhoods are rows and 
#' cell types are columns. The values in the matrix indicate the  
#' number of cells of a given type within a neighbourhood. 
#' 
#' @param neighbourhoods - a named list of neighbourhoods where a neighbourhood 
#' is a set of cells. 
#' @param cellTypes - named vector of cell types where names are each cell and
#' cell types are a factor
#' @return a matrix of neighbourhoods by cell types
#' @export
## #' @examples

computeNeighbourhoodByCTMatrix = function(neighbourhoods, cellTypes){
  
    cellNames = names(cellTypes)
    neighbourhoodByCT = matrix(nrow = 0, ncol = length(levels(cellTypes)))
    
    for (neighbourhood in neighbourhoods){
        neighbourhoodByCT = rbind(neighbourhoodByCT,table(cellTypes[neighbourhood]))
    }
    rownames(neighbourhoodByCT) = names(neighbourhoods)
    
    return(neighbourhoodByCT)
}
  

## #' ## ####################################################
## #' #' This function computes a matrix where cells are rows and 
## #' #' cell types are columns. The values in the matrix indicate the  
## #' #' number of spatial neighbours of a given cell type a particular cell has. 
## #' #' 
## #' #' @param  spatialGraph - a spatial graph in neighbour list format.
## #' #' @param cellTypes - named vector of cell types where names are each cell and
## #' #' cell types are a factor.
## #' #' @return a matrix of cells by cell types.
## #' #' @export
## #' #' @examples
## #' 
## #' computeCellByCTMatrix = function(spatialGraph, cellTypes){
## #'   
## #'   cellNames = names(cellTypes)
## #'   cellByCt = matrix(nrow = 0, ncol = length(levels(cellTypes)))
## #' 
## #'   for (cell in cellNames){
## #'     neighbours = 
## #'       c(spatialGraph[,2][cell == spatialGraph[,1]],
## #'         spatialGraph[,1][cell == spatialGraph[,2]])
## #'     cellByCt = rbind(cellByCt,table(xenium.obj@active.ident[neighbours]))
## #'   }
## #'   
## #'   return(cellByCt)
## #' }


## ####################################################
#' This function creates a seurat object using a neighbourhood by cell type 
#' matrix
#' 
#' @param neighbourhoodByCT - a matrix of neighbourhoods by cell types
#' @param resolution - resolution for clustering (default 0.1)
#' @param npcs - number of pcs used for PCA, defaults to 10
#' @param n.neighbors - number of neighbors used by UMAP, defaults to 30
#' @param compute_UMAP - defaults to TRUE
#' @param transpose - defaults to FALSE
#' @return a seurat object based on a neighbourhood by cell type matrix,
#' containing clusters and UMAP.
#' @import Seurat
#' @export
## #' @examples

computeNeighbourhoodByCTSeurat= function(neighbourhoodByCT, resolution = 0.1, 
                                         npcs = 10, n.neighbors = 30L, 
                                         compute_UMAP = T, transpose = F){
    dataMatrix = t(neighbourhoodByCT)
    neighbourhoodSeurat = CreateSeuratObject(dataMatrix)
    neighbourhoodSeurat[['RNA']]$data = dataMatrix
    neighbourhoodSeurat = ScaleData(neighbourhoodSeurat)
    neighbourhoodSeurat = RunPCA(neighbourhoodSeurat, assay = "RNA", 
                                 features = rownames(neighbourhoodSeurat), 
                                 npcs = npcs)
    if (transpose){
        neighbourhoodSeurat = RunUMAP(neighbourhoodSeurat,assay='RNA',
                                      n.neighbors = n.neighbors)
    } else{
        neighbourhoodSeurat = RunUMAP(neighbourhoodSeurat,assay='RNA',
                                      features=rownames(neighbourhoodSeurat), 
                                      n.neighbors = n.neighbors)
    }
    if (transpose){
        neighbourhoodSeurat = FindNeighbors(neighbourhoodSeurat, dims = 1:npcs)
    } else{
        neighbourhoodSeurat = FindNeighbors(neighbourhoodSeurat, 
                                            features=rownames(neighbourhoodSeurat))
    }
    neighbourhoodSeurat = FindClusters(neighbourhoodSeurat,
                                       resolution=resolution)
    return(neighbourhoodSeurat)
}

## ####################################################
#' This function adds a force directed graph embedding to a seurat object
#' 
#' @param seuratObj - a seurat object. 
#' @param graph - which graph to extract.  Defaults to
#' paste0(f@active.assay,'_snn')
#' @return a seurat object with a "graph" dimensionality reduction.
#' @import Seurat  
#' @import igraph
#' @export
## #' @examples

computeGraphEmbedding = function(seuratObj, graph=defaultGraph(seuratObj)){
  graph = seuratObj@graphs[[graph]]
  igraphGraph = igraph::graph_from_adjacency_matrix(graph)
  graphLayout = igraph::layout_with_fr(igraphGraph)
  colnames(graphLayout) = c("graph_1","graph_2")
  rownames(graphLayout) = colnames(seuratObj)
  graphDimReduc = CreateDimReducObject(embeddings = graphLayout,   
                                       key = "graph",   assay = "RNA")
  seuratObj[["graph"]] = graphDimReduc
  return(seuratObj)
}



## ####################################################
#' This function takes a nearest neighbour graph and a radius
#' and finds the combinatorial ball around each cell
#'
#' @param NN - a nearest neighbour graph
#' @param radius - combinatorial radius
#' @return This returns a named list.  In each case, the name is the
#'     cell at the center of the ball.  Each item in the list is a
#'     named numeric.  The values give combinatorial distance from the
#'     center and the names give the cells.
#' @export
findCombinatorialNeighbourhoods = function(NN,radius)
{
    combinatorialNbhdsImpl = function(edges,radius)
    {
        ## ####################################################
        ## The base case:
        if(radius == 0)
        {
            ball = list()
            vertices = names(edges)
            for(v in vertices)
            {
                a = 0
                names(a) = v
                ball[[v]] = a
            }
            return(ball)
        }
        ## ####################################################
        ## Recurse:
        else
        {
            smaller = combinatorialNbhdsImpl(edges,radius-1)
            for(v in names(smaller))
            {
                frontier = smaller[[v]][smaller[[v]]==radius-1]
                proposed = c()
                for(frontierV in names(frontier))
                    proposed = c(proposed,edges[[frontierV]])
                proposed = unique(proposed)
                idx = proposed %in% names(smaller[[v]])
                proposed = proposed[!idx]
                theNew = rep(radius,length(proposed))
                names(theNew) = proposed
                smaller[[v]] = c(smaller[[v]],theNew)
            }
            return(smaller)
        }
    } ## End of impl

    ## ####################################################    
    ## Get the edges:
    vertices = unique(c(NN$nodeA,NN$nodeB))
    edges = list()
    for(v in vertices)
    {
        idx = NN$nodeA == v
        neighboursB = unique(NN$nodeB[idx])
        idx = NN$nodeB == v
        neighboursA = unique(NN$nodeA[idx])
        
        edges[[v]] = unique(c(neighboursA,neighboursB))
    }

    return(combinatorialNbhdsImpl(edges,radius))
}


## ####################################################
#' This function reduces combinatorial balls to an
#' extended nearest neighbour graph
#'
#' @param balls - A named list as produced by the function
#'     findCombinatorialNeighbourhoods
#' @return This returns a nearest neighbour data frame where the cells
#'     in each combinatorial ball are now considered the center's
#'     neighbours.
#' @export
reduceCombinatorialBalls = function(balls)
{
    nodeA = c()
    nodeB = c()
    for(cell in names(balls))
    {
        newNeighbours = names(balls[[cell]])
        N = length(newNeighbours)
        nodeA = c(nodeA,rep(cell,N))
        nodeB = c(nodeB,newNeighbours)
    }
    NN = data.frame(nodeA,nodeB)

    return(NN)
}

## ####################################################
#' For each cell type, this function looks at the neighbourhoods
#' around cells of that type and discovers the fractions of those
#' cells of each type.
#'
#' @param nbhdByCellType - A matrix whose rows are neighbourhoods
#' each denoted by the cell at their center, whose columns are
#' cell types, and whose entries are counts.
#' @param seurat_clusters - a named vector whose names are the cells
#' and whose entries are their seurat_clusters.
#' @return A square matrix whose rownames and colnames are the
#' seurat_clusters as character strings.  Each row corresponds
#' to neighbourhoods around all cells of that type and the entries
#' give the fractions of those neighbourhoods occupied by cells
#' of each type.
#' @export
cellTypesPerCellTypeMatrix = function(nbhdByCellType,seurat_clusters)
{
    clusters = unique(seurat_clusters)
    clusters = clusters[order(clusters)]
    N = length(clusters)
    M = matrix(0,nrow=N,ncol=N)
    Clusters = as.character(clusters)
    rownames(M) = Clusters
    colnames(M) = Clusters

    for(cell in rownames(nbhdByCellType))
    {
        type = as.character(seurat_clusters[cell])
        M[type,] = M[type,] + nbhdByCellType[cell,]
    }

    rowTotals = rowSums(M)
    MM = M
    for(i in 1:nrow(M))
        MM[i,] = MM[i,]  / rowTotals[i]

    return(MM)
}

## ####################################################
#' This function converts a matrix as found by
#' cellTypesPerCellTypeMatrix into a directed igraph
#' whose vertices correspond to seurat_clusters and whose
#' edge correspond to occupancy fraction.
#'
#' @param M - a matrix as found by cellTypesPerCellTypeMatrix
#' @param colors - a named vector of colors used to color the
#' vertices of the graph.  The names are the seurat_clusters
#' as character strings.
#' @param selfEdges - a logical which determines whether to include
#' self edges.  Defaults to FALSE
#' @param minWeight - Allows one to exclude edges of low weight.
#' Defaults to 0, thus including all edges.
#' @param edgeWeighting - a parameter used to thicken the edges
#' in the display.  Defaults to 20.
#' @param plotGraph - a logical which determines whether to
#' plot the graph.  Defaults to TRUE.
#' @return This returns a directed igraph whose vertices are
#' the cell types and whose arrows indicate "ownership" of
#' cells of the target type by neighbourhoods of cells of the
#' source type.  Layout is done witht the FR algorithm and
#' coordinates are found in G$coords.  If colors were supplied
#' these are found in V(G)$color.  Edge weights and widths are
#' found in E(G)$weight and E(G)$width.
#' @export
cellTypesPerCellTypeGraphFromMatrix = function(M,
                                               colors=NULL,
                                               selfEdges=FALSE,
                                               minWeight=0,
                                               edgeWeighting=20,
                                               plotGraph=TRUE)
{
    idx = M >= minWeight
    M[!idx] = 0
    
    G = graph_from_adjacency_matrix(M,
                                    mode='directed',
                                    weighted=TRUE,
                                    diag=selfEdges)

    if(! is.null(colors))    
        V(G)$color = colors[names(V(G))]
    G$coords = layout_with_fr(G)
    E(G)$width = edgeWeighting * E(G)$weight

    if(plotGraph)
    {
        if(! is.null(colors))
            plot(G,
                 layout=G$coords,
                 vertex.color=V(G)$color,
                 edge.width=E(G)$width)
        else
            plot(G,
                 layout=G$coords,
                 edge.width=E(G)$width)
    }
    return(G)    
}

## ####################################################
#' This function takes a neighbourhood-by-cell type
#' matrix and produces a directed igraph showing the
#' fractions of cells of each type in the neighbourhoods
#' around cells of each type.
#'
#' @param nbhdByCellType - A matrix whose rows are neighbourhoods
#' each denoted by the cell at their center, whose columns are
#' cell types, and whose entries are counts.
#' @param seurat_clusters - a named vector whose names are the cells
#' and whose entries are their seurat_clusters.
#' @param colors - a named vector of colors used to color the
#' vertices of the graph.  The names are the seurat_clusters
#' as character strings.
#' @param selfEdges - a logical which determines whether to include
#' self edges.  Defaults to FALSE
#' @param minWeight - Allows one to exclude edges of low weight.
#' Defaults to 0, thus including all edges.
#' @param edgeWeighting - a parameter used to thicken the edges
#' in the display.  Defaults to 20.
#' @param plotGraph - a logical which determines whether to
#' plot the graph.  Defaults to TRUE.
#' @return This returns a directed igraph whose vertices are
#' the cell types and whose arrows indicate "ownership" of
#' cells of the target type by neighbourhoods of cells of the
#' source type.  Layout is done witht the FR algorithm and
#' coordinates are found in G$coords.  If colors were supplied
#' these are found in V(G)$color.  Edge weights and widths are
#' found in E(G)$weight and E(G)$width. 
#' @export
cellTypesPerCellTypeGraph = function(nbhdByCellType,
                                     seurat_clusters,
                                     colors=NULL,
                                     selfEdges=FALSE,
                                     minWeight=0,
                                     edgeWeighting=20,
                                     plotGraph=TRUE)
{
    M = cellTypesPerCellTypeMatrix(nbhdByCellType,seurat_clusters)
    G = cellTypesPerCellTypeGraphFromMatrix(M,
                                            colors=colors,
                                            selfEdges=selfEdges,
                                            minWeight=minWeight,
                                            edgeWeighting=edgeWeighting,
                                            plotGraph=plotGraph)

    
    return(G)
}

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
    
    if(species == 'human')
        return(humanLRN)

    if(species == 'mouse')
        return(mouseLRN)
}

## ####################################################
#' This functions takes an Seurat object, its species
#' and a ligand receptor network and subsets the ligand
#' receptor network to those pairs that occur in the
#' panel
#'
#' @param obj - a Seurat object
#' @param species - either 'human' or 'mouse'
#' @param lrn - a ligand-receptor network, i.e., a
#' data frame with columns from and to.  By default, it
#' retrieves the nichenetr ligand receptor network
#' @return This returns a data frame with columns ligand and
#' receptor
#' @export
getLigandReceptorPairsInPanel = function(obj,species,
                                         lrn = getLigandReceptorNetwork(species))
{
    stopifnot(species %in% c('mouse','human'))

    ## The panel
    panel = rownames(obj)

    ## Ligands and receptors:
    pairs = paste(lrn$from,lrn$to,sep='_')
    ligands = unique(lrn$from)
    receptors = unique(lrn$to)
    ligandsFound = intersect(ligands,panel)
    receptorsFound = intersect(receptors,panel)

    panelPairs = c()
    for(a in ligandsFound)
        for(b in receptorsFound)
            panelPairs = c(panelPairs,paste(a,b,sep='_'))
    idx = panelPairs %in% pairs
    pairsFound = panelPairs[idx]

    a = str_split(pairsFound,'_')
    pairsFoundDF = data.frame(ligand=unlist(lapply(a,function(x) return(x[1]))),
                              receptor=unlist(lapply(a,function(x) return(x[2]))))

    return(pairsFoundDF)
}

## ####################################################
#' This function takes a Seurat object, a set of ligand receptor
#' pairs and a set of edges denoting neighbouring cells and
#' annotates these with the ligand receptor interactions taking
#' place on those edges in each direction.
#'
#' @param obj - a Seurat object
#' @param pairDF - a data frame giving the ligand-receptor pairs
#' @param delaunayNeighbours - a data frame of neighbouring
#' cell pairs.  Note that each row is a directed edge (A,B) so
#' that this data frame should have both the edge (A,B) and the
#' edge (B,A)
#' @return This returns a data frame whose first two columns give
#' the neighbouring cells.  Each of the remaining columns is a logical
#' corresponding to a ligand-receptor pair telling whether the ligand
#' is expressed in the first cell and the receptor is expressed in the
#' second cell.
#' @export
getInteractionsOnEdges = function(obj,pairDF,delaunayNeighbours)
{
    ## Discretize expression:
    M = FetchData(obj,rownames(obj),layer='count')
    M = data.matrix(t(M))
    cutoff = 0
    M = M > cutoff

    ## Find the interactions on the edges:
    edges = delaunayNeighbours

    for(i in 1:nrow(pairDF))
    {
        tag = paste(pairDF$ligand[i],pairDF$receptor[i],sep='_')
        edges[,tag] = (M[pairDF$ligand[i],edges$nodeA] &
                   M[pairDF$receptor[i],edges$nodeB])
    }

    return(edges)
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
countInteractionsPerCell = function(edges,sourceOrTarget)
{
    stopifnot(sourceOrTarget %in% c('source','target'))
    
    if(sourceOrTarget == 'source')
        by = factor(edges$nodeA)
    if(sourceOrTarget == 'target')
        by = factor(edges$nodeB)

    edges = edges[,3:ncol(edges)]
    interactionCountDF = aggregate(edges,by=list(by),FUN=sum) 

    interactionCountDF$Group.1 =
        str_replace(interactionCountDF$Group.1,'Group.','')
    rownames(interactionCountDF) = interactionCountDF$Group.1
    names(interactionCountDF)[1] = 'cell'
     
    return(interactionCountDF)
}


## ####################################################
#' This takes a data frame of interaction counts as found
#' by countInteractionsPerCell(), the underlying Seurat object
#' and the neighbourhood Seurat object and annotates the counts
#' with the cell type and the neighbourhood type corresponding
#' to the cells of the interaction counts.
#'
#' @param interactionCounts - as found by countInteractionsPerCell()
#' @param obj - a Seurat object
#' @param nbhdObj - a neighbourhood x cell type Seurat object
#' @return This returns the interaction counts annotated with the
#' cell type and neighbourhood type of each cell.
#' @export
annotateInteractionCounts = function(interactionCounts,obj,nbhdObj)
{
    ## Get nbhd and cell types:
    annotated = data.frame(cell=colnames(obj),
                           cellType=obj$seurat_clusters)
    rownames(annotated) = annotated$cell

    ## Append nbhd type:
    annotated$nbhdType = nbhdObj$seurat_clusters[annotated$cell]

    ## Bung in the interaction counts:
    pairs = names(interactionCounts)[2:ncol(interactionCounts)]

    annotated[,pairs] = 0
    annotated[rownames(interactionCounts),pairs] =
        interactionCounts[rownames(interactionCounts),pairs]

    return(annotated)
}
    

