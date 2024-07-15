
## ####################################################
#' This function makes the function whichretrieves and
#' makes example data objects.
#'
#' @return This returns the function which retrieves
#' and makes example data objects. The latter saves any
#' object it has found for quicker return.  Using the value 'list'
#' causes it to return the list of all objects found so far.
#' @export
#' @examples
#' getExample = make.getExample()
#' ## Provided:
#' smallXenium = getExample('smallXenium')
#' ## Computed:
#' delaunayNeighbours = getExample('delaunayNeighbours')
make.getExample = function()
{
    
    alreadyFound = list(S=S,
                        smallXenium=smallXenium,
                        ligandReceptorResults=ligandReceptorResults,
                        moransI=moransI,
                        moransILigandReceptor=moransILigandReceptor,
                        hallmark=hallmark,
                        humanLRN=humanLRN,
                        mouseLRN=mouseLRN)
    
    
    getExample = function(whichOne)
    {
        if(whichOne %in% names(alreadyFound))
            return(alreadyFound[[whichOne]])
        
        if(whichOne == 'list')
            return(alreadyFound)
        
        if(! whichOne %in% exampleObjects())
            stop(paste(whichOne,
                       'is not an available example'))
        
        ## ####################################################
        ## Otherwise, compute:
        
        if(whichOne == 'STranspose')
        {
            S = getExample('S')
            answer = transposeSeuratObject(S)
        }
        
        if(whichOne == 'S_sce')
        {
            S = getExample('S')
            answer = SeuratToSCE(S)
        }
        
        if(whichOne == 'STranspose_sce')
        {
            STranspose = getExample('STranspose')
            answer = SeuratToSCE(STranspose)
        }
        
        if(whichOne == 'delaunayNeighbours')
        {
            centroids = getExample('centroids')
            delaunayNeighbours = computeNeighboursDelaunay(centroids)
            answer = delaunayNeighbours
        }
        
        if(whichOne == 'NBHDByCTSeuratExtended')
        {
            NBHDByCTMatrixExtended = getExample('NBHDByCTMatrixExtended')
            answer = computeNBHDVsCTSeurat(NBHDByCTMatrixExtended,
                                           verbose=FALSE)
        }
        
        if(whichOne == 'CTByNBHDSeuratExtended')
        {
            ## Do we use this anywhere?
        }
        
        if(whichOne == 'edgeSeurat')
        {
            ligandReceptorResults = getExample('ligandReceptorResults')
            centroids = getExample('centroids')
            answer = computeEdgeSeurat(ligandReceptorResults, centroids)
        }
        
        if(whichOne == 'centroids')
        {
            smallXenium = getExample('smallXenium')
            centroids = GetTissueCoordinates(smallXenium)
            rownames(centroids) = centroids$cell
            answer = centroids
        }
        
        if(whichOne == 'clusters')
        {
            smallXenium = getExample('smallXenium')
            answer = smallXenium@active.ident
        }
        
        if(whichOne == 'extendedNeighbours')
        {
            extendedNeighboursList = getExample('extendedNeighboursList')
            answer = collapseExtendedNBHDs(extendedNeighboursList, 4)
        }
        
        if(whichOne == 'extendedNeighboursList')
        {
            delaunayNeighbours = getExample('delaunayNeighbours')
            answer = getExtendedNBHDs(delaunayNeighbours,4)
        }
        
        if(whichOne == 'NBHDByCTSeurat')
        {
            NBHDByCTMatrix = getExample('NBHDByCTMatrix')
            answer = computeNBHDVsCTSeurat(NBHDByCTMatrix,verbose=FALSE)
        }
        
        if(whichOne == 'NN')
        {
            STranspose = getExample('STranspose')
            answer = getNearestNeighborListsSeurat(STranspose)
        }
        
        if(whichOne == 'CTByNBHDSeurat')
        {
            NBHDByCTMatrix = getExample('NBHDByCTMatrix')
            answer =
                computeNBHDVsCTSeurat(t(NBHDByCTMatrix), npcs = 10, 
                                      transpose = TRUE, resolution = 1, n.neighbors = 5,
                                      verbose=FALSE)
        }
        
        if(whichOne == 'edgeNeighbours')
        {
            delaunayNeighbours = getExample('delaunayNeighbours')
            answer =  computeEdgeGraph(delaunayNeighbours)
        }
        
        if(whichOne == 'colours')
        {
            clusters = getExample('clusters')
            clusterNames = levels(clusters) 
            colours = c('#5A5156',
                        '#E4E1E3',
                        '#F6222E',
                        '#FE00FA',
                        '#16FF32',
                        '#3283FE',
                        '#FEAF16',
                        '#B00068',
                        '#1CFFCE',
                        '#90AD1C',
                        '#2ED9FF',
                        '#DEA0FD',
                        '#AA0DFE',
                        '#F8A19F',
                        '#325A9B',
                        '#C4451C')
            names(colours) = clusterNames
            answer = colours
        }
        
        if(whichOne == 'euclideanNeighbours')
        {
            centroids = getExample('centroids')
            answer =  computeNeighboursEuclidean(centroids,threshold=20)
        }
        
        if(whichOne == 'NBHDByCTMatrixExtended')
        {
            extendedNeighbours = getExample('extendedNeighbours')
            clusters = getExample('clusters')
            answer =  computeNBHDByCTMatrix(extendedNeighbours, clusters)
        }
        
        if(whichOne == 'NBHDByCTMatrix')
        {
            delaunayNeighbours = getExample('delaunayNeighbours')
            clusters = getExample('clusters')
            answer = computeNBHDByCTMatrix(delaunayNeighbours, clusters) 
        }
        
        if(whichOne == 'shorterHallmark')
        {
            hallmark = getExample('hallmark')
            answer = hallmark[1:10]
        }
        
        if(whichOne == 'clusterDF')
        {
            ## I don't think we use this.
        }
        
        if(whichOne == 'cellTypesPerCellTypeMatrixExtended')
        {
            NBHDByCTMatrixExtended = getExample('NBHDByCTMatrixExtended')
            clusters = getExample('clusters')
            answer = computeCellTypesPerCellTypeMatrix(NBHDByCTMatrixExtended,clusters)
        }
        
        if(whichOne == 'cellTypesPerCellTypeMatrix')
        {
            NBHDByCTMatrix = getExample('NBHDByCTMatrix')
            clusters = getExample('clusters')
            answer = computeCellTypesPerCellTypeMatrix(NBHDByCTMatrix,clusters)
        }
        
        if(whichOne == 'cellTypesPerCellTypePValues')
        {
            delaunayNeighbours = getExample('delaunayNeighbours')
            clusters = getExample('clusters')
            answer = computeNeighbourEnrichment(delaunayNeighbours, 
                                                clusters, verbose = FALSE)
        }

        if(whichOne == 'averageExpMatrix')
        {
            STranspose = getExample('STranspose')
            answer = getAverageExpressionMatrix(S,STranspose,layer='data')
        }
                
        ## Save and return:
        alreadyFound[[whichOne]] <<- answer
        
        return(answer)
    }

    return(getExample)
}


## ####################################################
#' This returns the names of available example objects.
#'
#' @return A character vector of the names of available
#' example data objects
#' @export
#' @examples
#' availableObjects = exampleObjects()
exampleObjects = function()
{
    objects = c('S',
                'STranspose',
                'STranspose_sce',
                'S_sce',
                'NBHDByCTSeuratExtended',
                'edgeSeurat',
                'smallXenium',
                'extendedNeighbours',
                'extendedNeighboursList',
                'NBHDByCTSeurat',
                'NN',
                'CTByNBHDSeurat',
                'edgeNeighbours',
                'ligandReceptorResults',
                'centroids',
                'delaunayNeighbours',
                'mouseLRN',
                'euclideanNeighbours',
                'NBHDByCTMatrixExtended',
                'hallmark',
                'humanLRN',
                'NBHDByCTMatrix',
                'clusters',
                'shorterHallmark',
                'averageExpMatrix',
                'moransI',
                'cellTypesPerCellTypeMatrixExtended',
                'cellTypesPerCellTypeMatrix',
                'cellTypesPerCellTypePValues',
                'moransILigandReceptor',
                'colours')

    return(objects)
}




