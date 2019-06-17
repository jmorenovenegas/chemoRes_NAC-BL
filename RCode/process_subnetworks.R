## Process the module subnetworks for linkcomm
get_gml <- function(subnetworkfile){
  interactome = strsplit(subnetworkfile, split = "_")[[1]][1]
  module = strsplit(subnetworkfile, split = "_")[[1]][4]
  df = read.csv(file = file.path('Results/expanded_modules_networks', subnetworkfile))
  df = df[,2:3]
  names(df) = c("node1", "node2")
  network = graph_from_data_frame(df,directed = F)
  network = igraph::simplify(network)
  cl <- clusters(network) # get largest connected component
  lcc = induced.subgraph(network, V(network)[which(cl$membership == which.max(cl$csize))])
  lcc = igraph::simplify(lcc)
  write_graph(lcc, file = paste('Results/linkcomm', paste(interactome, "_subnetwork_connected.gml",sep = ""), sep = "/"), format = "gml")
}
