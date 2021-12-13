
library(dplyr)
library(networkD3)
library(stringr)

## ####################################################
sankeyFromMatrix = function(M,disambiguation=c('R_','C_'),
                            fontSize=20,minus='red',plus='blue')
{
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
    idx = links$value > 0
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
                      colourScale=linkColor, LinkGroup="group")

    return(p)
}

set.seed(5)
m = 10
n = 15
M = matrix(runif(m*n),nrow=m)
M = M - .5

p = sankeyFromMatrix(M,minus='pink',plus='steelblue')
print(p)

