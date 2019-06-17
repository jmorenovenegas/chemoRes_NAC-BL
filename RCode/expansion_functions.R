# pybin : Path to python source
# network_path: Path to network edge list file (interactome)
# seed_genes_path: Path to seed genes list file
# n: Number of DIAMOnD iterations
# alpha: Weight of the seeds, default value is set to 1
# diamond_path : DIAMOnD path
# outfilename : Output file name

run_DIAMOnD <- function(pybin, diamond_path, network_path, seed_genes_path, n, alpha = 1, outfilename){
  system(paste(pybin, diamond_path, network_path, seed_genes_path, n, alpha, outfilename, sep = " "))
}

expand_module <- function(pybin, diamond_path, module, n=300, network_path){
  interactome <- unlist(strsplit(network_path, "/|\\."))[2]
  
  aux.df <- read.csv(file = paste("Results/modules_enrichments/",module,"_module_genes.csv", sep=""))
  genes <- as.numeric(aux.df$entrezID)
  seed_genes_path <- paste("Results/",module,"_seed_genes.txt", sep="")
  diamond_res_outfilename <- "Results/DIAMOnD_res.csv"
  
  # Write temporary module seed genes file
  write.table(genes ,file = seed_genes_path, row.names = F, col.names = F)
  
  run_DIAMOnD(pybin, diamond_path, network_path, seed_genes_path, n, outfilename = diamond_res_outfilename)
  
  aux.df <- read.csv(file = diamond_res_outfilename, header = T, sep = "\t")
  new_genes <- aux.df$DIAMOnD_node
  expanded_module_genes <- c(genes, new_genes)
  
  # Delete temporary module seed genes file
  system(paste("rm", seed_genes_path, sep=" "))
  # Delete temporary DIAMOnD result file
  system(paste("rm", diamond_res_outfilename, sep=" "))
  res.filename <- paste("Results/expanded_modules/",interactome,"_expanded_",module,"_module.txt", sep = "")
  
  cat(paste("Expanded ",module, " module genes in ",res.filename, sep = ""))
  # Write module expanded genes file
  write.table(expanded_module_genes, file = res.filename, row.names = F, col.names = F)
}

expanded_module_enrichment <- function(genes_file_path){
  outfilename <- paste("Results/expanded_modules_enrichments/",unlist(strsplit(genes_file_path, '/|\\.'))[3], "_enrichment.csv",sep = "")
  genes <- readLines(con = genes_file_path)
  
  cat(paste(genes_file_path, "enrichment...", '\n'))
  
  enrich_rs <- enrich_cp(genes, comparison = "ChemoRes")
  enrich_summary <- enrich_rs$summary %>% arrange(qvalue)
  enrich_summary <- convert_enriched_ids(enrich_summary, entrezsymbol = entrezsymbol) %>% filter(qvalue <=0.01) %>% arrange(-Count)
  
  write.csv(enrich_summary, file = outfilename)
  
  cat(paste(genes_file_path, "enrichment DONE", '\n'))
  cat(paste("Results in", outfilename, '\n'))
  return("OK")
}

entrez2symbol <- function(entrez, dict){
  indx <- match(entrez, dict$ENTREZID)
  if(!is.na(indx)){
    return(dict$SYMBOL[indx])
  }else{
    return(NA)
  }
}

network_building <- function(genes_file_path){
  expanded_module <- unlist(strsplit(genes_file_path, '/|\\.'))[3]
  outfilename <- paste("Results/expanded_modules_networks/",expanded_module, "_network.csv",sep = "")
  interactome_name <- paste(unlist(strsplit(expanded_module, '_'))[1], 
                            unlist(strsplit(expanded_module, '_'))[2],
                            sep='_')
  
  interactome_file_path <- paste("Interactomes/",interactome_name,".tsv",sep="")
  
  genes <- as.numeric(readLines(con = genes_file_path))
  dict <- bitr(genes, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
  
  if(interactome_file_path == 'Interactomes/STRINGdb_interactome.tsv'){
    interactome <- read.table(file = interactome_file_path, header = F)[1:3]
    colnames(interactome) <- c('gene_id1', 'gene_id2','score')
  }else{
    interactome <- read.table(file = interactome_file_path, header = F)[1:2]
    colnames(interactome) <- c('gene_id1', 'gene_id2')
  }
  
  network <- interactome %>% 
    dplyr::filter(gene_id1 %in% genes & gene_id2 %in% genes)
  
  gene_symbol1 <- unlist(lapply(network$gene_id1, entrez2symbol, dict))
  gene_symbol2 <- unlist(lapply(network$gene_id2, entrez2symbol, dict))
  
  network <- network %>%
    dplyr::mutate(gene_symbol1 = gene_symbol1, 
                  gene_symbol2 = gene_symbol2)
  write.csv(network, file = outfilename, row.names = F)
}

