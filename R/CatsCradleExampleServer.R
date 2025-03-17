
## ####################################################
#' This function makes the function whichretrieves and
#' makes example data objects.
#'
#' @return This returns the function which retrieves
#' and makes example data objects. The latter saves any
#' object it has found for quicker return.  Using the value 'list'
#' causes it to return the list of all objects found so far.
#' @importFrom SummarizedExperiment assays
#' @importFrom msigdbr msigdbr
#' @export
#' @examples
#' getExample = make.getExample()
#' ## Provided:
#' smallXenium = getExample('smallXenium')
#' ## Computed:
#' delaunayNeighbours = getExample('delaunayNeighbours')
make.getExample = function()
{
    
    alreadyFound = list()
    toysFound = list()

    
    
    getExample = function(whichOne,toy=FALSE)
    {
        packageVars = c('exSeuratObj',
                        'smallXenium',
                        'ligandReceptorResults',
                        'moransI',
                        'moransILigandReceptor',
                        'humanLRN',
                        'mouseLRN')

        if(whichOne %in% names(alreadyFound) & !toy)
        {
            answer = alreadyFound[[whichOne]]
            return(answer)
        }

        if(whichOne %in% names(toysFound) & toy)
        {
            answer = toysFound[[whichOne]]
            return(answer)
        }
   
        
        if(whichOne %in% packageVars)
        {
            if(whichOne == 'exSeuratObj')
            {
                data(exSeuratObj,envir=environment())
                answer = exSeuratObj
                if(toy)
                {
                    data(seuratGenes,envir=environment())
                    data(seuratCells,envir=environment())
                    answer = answer[seuratGenes,seuratCells]
                }
            }

            if(whichOne == 'smallXenium')
            {
                data(smallXenium,envir=environment())
                answer = smallXenium
                if(toy)
                {
                    data(xeniumCells,envir=environment())
                    answer = answer[,xeniumCells]
                }
            }

            if(whichOne == 'ligandReceptorResults')
            {
                data(ligandReceptorResults,envir=environment())
                answer = ligandReceptorResults
            }

            if(whichOne == 'moransI')
            {
                data(moransI,envir=environment())
                answer = moransI
            }

            if(whichOne == 'moransILigandReceptor')
            {
                data(moransILigandReceptor,envir=environment())
                answer = moransILigandReceptor
            }

            if(whichOne == 'humanLRN')
            {
                data(humanLRN,envir=environment())
                answer = humanLRN
            }

            if(whichOne == 'mouseLRN')
            {
                data(mouseLRN,envir=environment())
                answer = mouseLRN
            }
        }
            

             
        if(whichOne == 'list')
            return(names(alreadyFound))
        
        if(! whichOne %in% exampleObjects())
            stop(paste(whichOne,
                       'is not an available example'))
        
        ## ####################################################
        ## Otherwise, compute:

        
        
        if(whichOne == 'STranspose')
        {
            exSeuratObj = getExample('exSeuratObj',toy)
            answer = transposeObject(exSeuratObj)
        }
        
        if(whichOne == 'S_sce')
        {
            exSeuratObj = getExample('exSeuratObj',toy)
            answer = SeuratToSCE(exSeuratObj)
        }
        
        if(whichOne == 'STranspose_sce')
        {
            STranspose = getExample('STranspose',toy)
            answer = SeuratToSCE(STranspose)
        }
        
        if(whichOne == 'delaunayNeighbours')
        {
            centroids = getExample('centroids',toy)
            delaunayNeighbours = computeNeighboursDelaunay(centroids)
            answer = delaunayNeighbours
        }
        
        if(whichOne == 'NBHDByCTSeuratExtended')
        {
            NBHDByCTMatrixExtended = getExample('NBHDByCTMatrixExtended',toy)
            answer = computeNBHDVsCTObject(NBHDByCTMatrixExtended,
                                           verbose=FALSE)
        }

        if(whichOne == 'NBHDByCTSingleCellExtended_sce')
        {
            NBHDByCTSeuratExtended = getExample('NBHDByCTSeuratExtended',toy)
            answer = SeuratToSCE(NBHDByCTSeuratExtended)
        }
        
        if(whichOne == 'CTByNBHDSeuratExtended')
        {
            ## Do we use this anywhere?
        }
        
        if(whichOne == 'edgeSeurat')
        {
            ligandReceptorResults = getExample('ligandReceptorResults',toy)
            centroids = getExample('centroids',toy)
            answer = computeEdgeObject(ligandReceptorResults, centroids)
        }
        
        if(whichOne == 'centroids')
        {
            smallXenium = getExample('smallXenium',toy)
            centroids = GetTissueCoordinates(smallXenium)
            rownames(centroids) = centroids$cell
            answer = centroids
        }
        
        if(whichOne == 'clusters')
        {
            smallXenium = getExample('smallXenium',toy)
            answer = smallXenium@active.ident
        }
        
        if(whichOne == 'extendedNeighbours')
        {
            extendedNeighboursList = getExample('extendedNeighboursList',toy)
            answer = collapseExtendedNBHDs(extendedNeighboursList, 4)
        }
        
        if(whichOne == 'extendedNeighboursList')
        {
            delaunayNeighbours = getExample('delaunayNeighbours',toy)
            answer = getExtendedNBHDs(delaunayNeighbours,4)
        }
        
        if(whichOne == 'NBHDByCTSeurat')
        {
            NBHDByCTMatrix = getExample('NBHDByCTMatrix',toy)
            answer = computeNBHDVsCTObject(NBHDByCTMatrix,verbose=FALSE)
        }
        
        if(whichOne == 'NN')
        {
            STranspose = getExample('STranspose',toy)
            answer = getNearestNeighbourLists(STranspose)
        }
        
        if(whichOne == 'CTByNBHDSeurat')
        {
            NBHDByCTMatrix = getExample('NBHDByCTMatrix',toy)
            answer =
                computeNBHDVsCTObject(t(NBHDByCTMatrix), npcs = 10, 
                                      transpose = TRUE, resolution = 1,
                                      n.neighbors = 5,
                                      verbose=FALSE)
        }
        
        if(whichOne == 'edgeNeighbours')
        {
            delaunayNeighbours = getExample('delaunayNeighbours',toy)
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

        if(whichOne == 'hallmark')
        {
            h = msigdbr(species = "mouse", category = "H")
            N = unique(h$gs_name)
            
            answer = list()
            for(ell in N)
                answer[[ell]] = h$gene_symbol[h$gs_name == ell]
        }
        
        if(whichOne == 'euclideanNeighbours')
        {
            centroids = getExample('centroids',toy)
            answer =  computeNeighboursEuclidean(centroids,threshold=20)
        }
        
        if(whichOne == 'NBHDByCTMatrixExtended')
        {
            extendedNeighbours = getExample('extendedNeighbours',toy)
            clusters = getExample('clusters',toy)
            answer =  computeNBHDByCTMatrix(extendedNeighbours, clusters)
        }
        
        if(whichOne == 'NBHDByCTMatrix')
        {
            delaunayNeighbours = getExample('delaunayNeighbours',toy)
            clusters = getExample('clusters',toy)
            answer = computeNBHDByCTMatrix(delaunayNeighbours, clusters) 
        }
        
        if(whichOne == 'shorterHallmark')
        {
            hallmark = getExample('hallmark',toy)
            answer = hallmark[seq_len(10)]
        }
        
        if(whichOne == 'clusterDF')
        {
            ## I don't think we use this.
        }
        
        if(whichOne == 'cellTypesPerCellTypeMatrixExtended')
        {
            NBHDByCTMatrixExtended = getExample('NBHDByCTMatrixExtended',toy)
            clusters = getExample('clusters')
            answer = computeCellTypesPerCellTypeMatrix(NBHDByCTMatrixExtended,
                                                       clusters)
        }
        
        if(whichOne == 'cellTypesPerCellTypeMatrix')
        {
            NBHDByCTMatrix = getExample('NBHDByCTMatrix',toy)
            clusters = getExample('clusters',toy)
            answer = computeCellTypesPerCellTypeMatrix(NBHDByCTMatrix,clusters)
        }
        
        if(whichOne == 'cellTypesPerCellTypePValues')
        {
            delaunayNeighbours = getExample('delaunayNeighbours',toy)
            clusters = getExample('clusters',toy)
            answer = computeNeighbourEnrichment(delaunayNeighbours, 
                                                clusters, verbose = FALSE)
        }

        if(whichOne == 'averageExpMatrix')
        {
            exSeuratObj = getExample('exSeuratObj',toy)
            STranspose = getExample('STranspose',toy)
            answer = getAverageExpressionMatrix(exSeuratObj,STranspose,layer='data')
        }

        if(whichOne == 'X_sce')
        {
            smallXenium = getExample('smallXenium',toy)
            answer = SeuratToSCE(smallXenium,spatial=TRUE)
        }
                
        ## Save and return:
        if(! toy)
            alreadyFound[[whichOne]] <<- answer
        if(toy)
            toysFound[[whichOne]] <<- answer

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
    objects = c('exSeuratObj',
                'STranspose',
                'STranspose_sce',
                'S_sce',
                'NBHDByCTSeuratExtended',
                'NBHDByCTSingleCellExtended_sce',                
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
                'colours',
                'X_sce')

    return(objects)
}






