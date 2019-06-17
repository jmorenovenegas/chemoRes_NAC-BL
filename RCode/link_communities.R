# Link Communities of a graph
# Note. Better run in Rstudio, so you can tune the number and size of linkcomms

link_comms <- function(networkfile){
  #### Load the network ###
  subnet = strsplit(networkfile, split = "\\.")[[1]][1]
  network = read_graph(file.path('Results/expanded_modules_networks', networkfile), format = "gml") 
  network.df <- get.data.frame(network,what = "edges")# linkcomm needs an edgelist
  names(network.df) <- c("node1", "node2")
  
  # The weighted network produces 4 communities, two of them tiny and the other two are massive (ca. 250 nodes). 
  # This might be due because of the ctkernel weight. I´m using the unweighted network as well
  network.df = network.df %>% dplyr::select(node1, node2)
  # the unweighted network has 11 communities (for metabolic ineractome derived).
  
  #### Get Link Communities
  
  # We get lc, an object of class linkcomm, which is a list containing the following components:
  # * numbers. An integer vector with the number of edges, nodes, and communities.
  # * hclust. An object of class hclust, which contains information about the hierarchical  clustering of links.
  # * pdmax. A numerical value indicating the height of the dendrogram at which the partition density is maximised.
  # * pdens. A numerical matrix with 2 columns; the first is the heights at which clusters  appear and the second is the partition density.
  # * nodeclusters. A data frame consisting of 2 columns; the first contains node names, and the second contains single community IDs for each node. All communities and their nodes are represented, but not necessarily all nodes.
  # * clusters. A list of integer vectors containing the link IDs that belong to each community. Community IDs are the numerical position of the communities in the list.
  # * edges. A data frame with 3 columns; the first two contain nodes that interact with each other, and the third is an integer vector of community IDs indicating community membership for each link.
  # * numclusters. A named integer vector. Names are node names and integer values are the number of communities to which each node belongs.
  # * clustsizes. A named integer vector. Names are community IDs and integer values indicate the number of nodes that belong in each community.
  # * igraph. An object of class igraph. The network is represented here as an igraph object.
  # * edgelist. A character matrix with 2 columns containing the nodes that interact with each other.
  # * directed. Logical indicating whether the network is directed.
  # * bipartite. Logical indicating whether the network is bi-partite.
  
  
  lc <- getLinkCommunities(network.df, directed = F ) # network is undirected
  cat(paste("This network has ", lc$numbers[3], "communities. \n"))
  
  # For nested comms run `getAllNestedComm(lc)`
  
  
  # We can get many communities and some completely nested. The metacommunities allow us to merge very similar linkcomms
  
  cr <- getClusterRelatedness(lc, hcmethod = "ward", cutat = 0.99)
  
  cutDendrogramAt(cr, cutat = 1.5) # plot the clustering cut at 1.5 and the number of communities produced
  
  mc <- newLinkCommsAt(lc, cutat = 0.997) #17 communities with Ghiassian
  getClusterRelatedness(mc, hcmethod = "ward")
  
  cat(paste("When merged similar communities we got", mc$numbers[3], "communities. \n"))
  
  
  # save the linkcomm object
  save(lc, file = 'Results/linkcomm', paste(subnet,"_lc.RData", sep = "/"))
  save(mc, file = 'Results/linkcomm', paste(subnet,"_mc.RData", sep = "/"))
  
  plot(mc, type = "graph")
  
  commsizes = mc$clustsizes
  binsize <- diff(range(commsizes))/40
  ggplot(NULL, aes(x= commsizes)) + 
    geom_histogram(binwidth=binsize,fill="white", colour="black") +
    theme_linedraw() +
    scale_x_continuous(breaks=seq(0,300,by = 50)) +
    labs(x = "community size", y = "number of communities")
  
  
  # Add the link communities as edge attribute and save the graph for further analysis.
  k <- as_data_frame(mc$igraph)
  names(k) <- c("node1", "node2")
  l <- mc$edges
  names(l) <- c("node1", "node2", "linkcomm")
  network.mc <- graph_from_data_frame(merge(k,l),directed = F)
  write_graph(network.mc, file = 'Results/linkcomm', paste(subnet,"_mc.gml", sep = "/"), format = "gml")
  
  mcs = as.numeric(unique((mc$edges)$cluster))
  for(i in mcs){
    mc_enrich(linkcomms = mc, community = i)
  }
}

#### Enrichment of the link communities. Note working with the merged communities
mc_enrich <- function(linkcomms, community){
  genes_df = mc$edges %>% dplyr::filter(cluster==community) %>% dplyr::select(node1, node2)
  genes = unique(c(genes_df$node1, genes_df$node2))
  enrich_rs <- enrich_cp(genes, comparison = paste(subnet,"linkcomm",community, sep = "_"))
  enrich_summary <- enrich_rs$summary %>% arrange(qvalue)
  enrich_summary <- convert_enriched_ids(enrich_summary, entrezsymbol = entrezsymbol) %>% filter(qvalue <=0.05) %>% arrange(-Count)
  write.csv(enrich_summary, file = file.path('Results/linkcomm',paste(subnet,"_enrich-linkcomm_",community,".csv", sep = "")),quote = F,row.names = F)                                       
}
