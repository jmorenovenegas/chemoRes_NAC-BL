## Core Script for enrichment analysis on GO BP, KEGG and Reactome
# Input: Dataframe with genes associated with a phenotype via gene co-expression networks
# Output: table with enrichment results

# Based on `https://github.com/hbc/clark-nanostring`, `https://rawgit.com/hbc/clark-nanostring/master/tumor-set/tumor-set.html#t-test`, `http://bioconductor.org/packages/devel/bioc/vignettes/NanoStringQCPro/inst/doc/vignetteNanoStringQCPro.pdf`

##### Libraries ######
orgdb = "org.Hs.eg.db"
biomart_dataset = "hsapiens_gene_ensembl"
keggname = "hsa"
reactomename = "human"
loadpkg("dplyr")
loadpkg("clusterProfiler")
loadpkg("org.Hs.eg.db")
loadpkg("biomaRt")
loadpkg("pathview")
loadpkg("ReactomePA")


mart = biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
entrezsymbol = biomaRt::getBM(attributes = c("entrezgene", "hgnc_symbol"), mart = mart)
entrezsymbol$entrezgene = as.character(entrezsymbol$entrezgene)
rm(mart)

summarize_cp = function(res, comparison) {
  summaries = data.frame()
  for (ont in names(res)) {
    ontsum = as.data.frame(res[[ont]])
    ontsum$ont = ont
    summaries = rbind(summaries, ontsum)
  }
  summaries$comparison = comparison
  return(summaries)
}

enrich_cp <- function(genes, comparison){
  
  mf = enrichGO(genes, OrgDb = orgdb, ont = "MF",
                qvalueCutoff = 1, pvalueCutoff = 1)
  cc = enrichGO(genes,  OrgDb = orgdb, ont = "CC",
                qvalueCutoff = 1, pvalueCutoff = 1)
  bp = enrichGO(genes,  OrgDb = orgdb, ont = "BP",
                qvalueCutoff = 1, pvalueCutoff = 1)
  kg = enrichKEGG(gene = genes, organism = keggname, pvalueCutoff = 1,
                  qvalueCutoff = 1, pAdjustMethod = "BH")
  re = enrichPathway(gene=genes, organism = reactomename, pvalueCutoff=1, qvalueCutoff = 1, pAdjustMethod = "BH")
  all = list(mf = mf, cc = cc, bp = bp, kg = kg, re = re)
  all[["summary"]] = summarize_cp(all, comparison)
  return(all)
}

convert_enriched_ids = function(res, entrezsymbol) {
  res = res %>% mutate(geneID = strsplit(as.character(geneID), "/")) %>% tidyr::unnest(geneID) %>% 
    left_join(entrezsymbol, by = c(geneID = "entrezgene")) %>% group_by(ID, 
                                                                        Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, Count, ont, 
                                                                        comparison) %>% summarise(geneID = paste(geneID, collapse = "/"), symbol = paste(hgnc_symbol, 
                                                                                                                                                         collapse = "/"))
  return(res)
}
