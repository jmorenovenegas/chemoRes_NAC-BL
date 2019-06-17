#########################################################################################
# STRINGdb interactome processing
#########################################################################################

# We download the whole interactome of Homo sapiens from https://string-db.org/cgi/download.pl?sessionId=ED2PAnsg3O9r&species_text=Homo+sapiens.
# We need to process the data frame to apply a score filter and to map the ensembl protein ids to gene ids.

# Read file
string_interactome <- read.table(file = "Data/stringdb_9606.protein.links.v11.0.txt", header = TRUE)
string_interactome <- string_interactome %>% filter(combined_score > 700) %>% 
  tidyr::separate(col = protein1, into = c("organism1","N1"), sep = '\\.', remove = T) %>% 
  tidyr::separate(col = protein2, into = c("organism2","N2"), sep = '\\.', remove = T) %>%
  dplyr::select(N1, N2, combined_score)

ensembl = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Dictionary ensembl protein id to gene id construction
interactome_ids <- unique(c(string_interactome$N1,string_interactome$N2))
dict <- getBM(attributes = c("ensembl_peptide_id","entrezgene"), 
              filters = "ensembl_peptide_id",
              values = dict,
              mart = ensembl)

# Mapping ensembl protein ids to entrez gene ids
ENSP_to_geneID <- function(ENSP_id, dict){
  if(ENSP_id %in% dict$ensembl_peptide_id){
    indx <- match(ENSP_id, dict$ensembl_peptide_id)
    return(dict$entrezgene[indx])
  }else{
    return(NA)
  }
}

string_interactome <- na.omit(string_interactome %>% 
                                dplyr::mutate(gene_id1 = unlist(lapply(string_interactome$N1, ENSP_to_geneID, dict))
                                              , gene_id2 = unlist(lapply(string_interactome$N2, ENSP_to_geneID, dict))) %>%
                                dplyr::select(gene_id1, gene_id2, combined_score))

# Save processed interactome in a file
write.table(string_interactome, row.names = F, col.names = F, file = "Interactomes/processed_STRINGdb_Hs_interactome.tsv")
rm(string_interactome)
rm(dict)
rm(interactome_ids)
rm(ensembl)

#########################################################################################
# iRef interactome processing
#########################################################################################

irefindex_tab <- iRefR::get_irefindex(tax_id = "9606", iref_version = "13.0", data_folder = "Data")
#load(file = "Data/9606.mitab.08122013.RData")

# creating an id conversion table to translate icROGs to entrez gene ids
id_conversion_table = create_id_conversion_table(irefindex_tab,
                                                 data_folder = "Data",
                                                 output_filename = "id_conversion_table_9606", 
                                                 IDs_to_include = "entrezgene/locuslink")

#load(file="Data/id_conversion_table_9606.RData")

# Mapping protein iROGs to gene ids
iROGid2geneid <- function(irogid, dict){
  if(irogid %in% dict$irogid){
    indx <- match(irogid, dict$irogid)
    return(dict$id_value[indx])
  }else{
    return(NA)
  }
}

# Apply the previous function, omit rows with missing values and keep gene ids and confidence.
iRef_interactome <- na.omit(irefindex_tab %>% 
                              dplyr::select(irogida, irogidb, confidence) %>%
                              dplyr::mutate(gene_id1 = as.numeric(unlist(lapply(irogida, 
                                                                                iROGid2geneid, 
                                                                                as.data.frame(id_conversion_table)))),
                                            gene_id2 = as.numeric(unlist(lapply(irogidb, 
                                                                                iROGid2geneid, 
                                                                                as.data.frame(id_conversion_table))))) %>% 
                              dplyr::select(gene_id1, gene_id2, confidence)) %>% 
                              dplyr::distinct();
# Filter the interactions on basis on confidence
iRef_interactome <- iRef_interactome %>% 
  tidyr::separate(confidence, c('hpr','lpr','np'), sep='\\|') %>%
  dplyr::select(-hpr, -lpr) %>%
  tidyr::separate(np, c('np', 'np_score'), sep = ':') %>%
  dplyr::select(-np);

# Save processed interactome in a file.
write.table(iRef_interactome, file = "Interactomes/iRef_interactome.tsv", row.names = F, col.names = F)
rm(iRef_interactome)
rm(irefindex_tab)
rm(id_conversion_table)
#########################################################################################
# Metabolism-centered interactome processing
#########################################################################################

# We download the interactome from https://bioinfo.uth.edu/ccmGDB/download/2072genes_PathwayCommons.txt?csrt=5379671582968607161. 
# This interactome has Pathway Commons format (see http://www.pathwaycommons.org/pc/sif_interaction_rules.do). 
# We have to process it to keep only the INTERACTS_WITH relationships.

met.interactome <- na.omit(read.table(file = "Data/2072genes_PathwayCommons.txt", header = F, stringsAsFactors = F, fill = T, col.names = c('??','entity1','relationship_type','entity2'))) %>% dplyr::filter(relationship_type == 'interacts-with') %>% dplyr::select(entity1, entity2)

symbol_ids <- unique(c(met.interactome$entity1, met.interactome$entity2))
dict <- bitr(symbol_ids, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

symbol_to_entrez <- function(symbol_id, dict){
  indx <- match(symbol_id, dict$SYMBOL)
  if(!is.na(indx)){
    return(dict$ENTREZID[indx])
  }else{
    return(NA)
  }
}

met.interactome <- na.omit(met.interactome %>% 
                             mutate(N1 = as.numeric(unlist(lapply(entity1, symbol_to_entrez, dict))), 
                                    N2 = as.numeric(unlist(lapply(entity2, symbol_to_entrez, dict))) %>% 
                                      dplyr::select(N1, N2)));
                           
# Save processed interactome in a file.
write.table(met.interactome, file = "Interactomes/metabolism-centered_interactome.tsv", row.names = F, col.names = F);
rm(symbol_ids);
rm(dict);
rm(met.interactome)