
library(networkD3)
library(dplyr)

source = c('A','B','B','C')
target = c('X','X','Y','Y')
value = c(5,-4,6,-3)

links = data.frame(source,target,value,
                   stringsAsFactors=FALSE)

nodes = unique(c(links$source,links$target))
nodes = data.frame(name=nodes,
                   stringsAsFactors=FALSE)

links$IDsource = match(links$source,nodes$name) - 1
links$IDtarget = match(links$target,nodes$name) - 1

idx = links$value > 0
links$group = ''
links$group[idx] = 'blue'
links$group[!idx] = 'red'

links$value = abs(links$value)

red = 'red'
blue = 'black'
linkColor = 'd3.scaleOrdinal() .domain(["red","blue"]) .range(["red", "green"])'

p = sankeyNetwork(Links = links, Nodes = nodes, Source = "IDsource", Target = "IDtarget", 
                  Value = "value", NodeID = "name", 
                  colourScale=linkColor, LinkGroup="group")

print(p)
