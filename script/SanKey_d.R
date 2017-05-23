  ### this script target: show the long Align with SanKeys diagram ###
  
  library(networkD3)
  library(RColorBrewer)
  
  mytest.dataset <- read.delim("~/Work/test.data.txt",header = T)
  mynames <- read.delim("~/Work/test.data.name.txt",header = T)
  sankeyNetwork(Links = mytest.dataset, 
                Nodes = mynames,
                colourScale = JS("d3.scaleOrdinal(d3.schemeCategory20c);"),
                #colourScale = JS('d3.scale.ordinal(range(["#7d3945","#e0677b", "#244457"]));'),
                NodeGroup = 'group',
                LinkGroup = 'ltype',
                Source = 'Source',
                Target = 'Target', 
                Value = 'value', 
                NodeID = 'name',sinksRight = T,
                #units = 'TWh', 
                fontSize = 12, nodeWidth = 10)