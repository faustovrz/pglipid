
get_entrez_from_pathway <- function(pathway_name, sp = "zma"){

  if(length(pathway_name) != 1){
    stop("Pathway name must be character(1)")
  }
  sp_pathway <- KEGGREST::keggList("pathway",sp)
  is_pathway <- grepl(pathway_name,sp_pathway)
  if(!any(is_pathway)){
    stop("Pathway search string not found!")
  }
  path_code <- names(sp_pathway[is_pathway])
  if(sum(is_pathway) > 1){
    stop("Ambiguous name. Found multiple pathways:\n", paste(path_code, collapse =" "))
  }
  pathway <- KEGGREST::keggGet(path_code)
  record_lines <- length(pathway[[1]]$GENE)
  entrez <- pathway[[1]]$GENE[c(TRUE,FALSE)] # selecting every other line
  return(
    list(name    = pathway[[1]]$PATHWAY_MAP,
         entry   = pathway[[1]]$ENTRY,
         entrez  = entrez %>% as.integer())
  )
}


xref <- read.table(config$ref$xref, sep = "\t", header = TRUE, na.strings = "")


get_cyc_xrefs <- function(entrez_ids){
  data.frame (
    Entrez =  entrez_ids
  ) %>%
    dplyr::left_join(xref) %>%
    dplyr::select( Entrez, v4_gene_model) %>%
    dplyr::left_join(corncyc_pathway) %>%
    dplyr::select(Entrez, v4_gene_model, Gene.id, Gene.name) %>%
    unique()
}


expand_cyc <- function(cyc_ids, pathway_cover){
    data.frame(
    Gene.name = union(
      cyc_ids,
      pathway_cover %>%
        dplyr::filter(n_test >1) %>%  # remove incidental pathways
        dplyr::select(Pathway.id) %>%
        dplyr::inner_join(corncyc_pathway) %>%
        dplyr::select(Gene.name) %>%
        dplyr::filter(Gene.name != "unknown") %>%
        dplyr::arrange(Gene.name) %>%
        dplyr::filter(!is.na(Gene.name)) %>%
        dplyr::pull(Gene.name)
    )
  )
}

get_unassigned_cyc <- function(cyc_ids){
  data.frame( Gene.name = cyc_ids) %>%
    dplyr::inner_join(genes_col, by = c(Gene.name = "NAME")) %>%
    dplyr::select(Gene.id = UNIQUE.ID, Gene.name) %>%
    dplyr::inner_join(orphan_enz) %>%
    dplyr::select(Gene.id,Gene.name) %>% unique()
}
