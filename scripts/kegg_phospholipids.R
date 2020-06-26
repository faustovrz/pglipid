source("corncyc.R")
library(dplyr)
library(KEGGREST)


pathway_names <-  c("Glycerophospholipid metabolism", "Glycerolipid metabolism")

candidates  <- character()

hits <- list()


for( kegg_name in pathway_names  ){
  kegg_result <- get_entrez_from_pathway(kegg_name)

  pathway_genes <- get_cyc_xrefs(kegg_result$entrez)

  hits[[kegg_result$entry]] <- pathway_genes$v4_gene_model

  pathway_cover <- corncyc_classify(
    pathway_genes %>%
    dplyr::filter(!is.na(Gene.name)) %>%  # Exclude genes that are not assigned to pathways
    dplyr::pull(Gene.name) %>%
    unique()
  )


  expanded_cyc <- expand_cyc(pathway_genes$Gene.name, pathway_cover)

  get_unassigned_cyc(expanded_cyc$Gene.name)


#  corncyc_classify(expanded_cyc$Gene.name) %>% print (n=50)

  write.csv(expanded_cyc,
            file = file.path(
                      config$output$dir,
                      paste0(kegg_result$entry,
                             "_expanded_cyc",
                             ".csv"
                      )
            ),
            quote = FALSE,
            row.names = FALSE)

  candidates <- union( candidates,expanded_cyc$Gene.name)
}

candidates <- candidates[!is.na(candidates)]
length(candidates)

# pathways to remove
remove <- c( "PWY-6446",  "PWY-5143", "PWY66-21",
             "PWY-6457",  "PWY-2501", "PWY-6333",
            "PWY1F-467", "PWY0-1313", "PWY-3941")

# pathways to add
add <- "PWY-6825"


curated_candidates <- corncyc_classify(candidates) %>%
  dplyr::filter(cover == 1) %>%
  dplyr::filter(!Pathway.id %in% remove) %>%  # manually remove pathways
  dplyr::select(Pathway.id) %>%               # manually add pathways
  dplyr::add_row(Pathway.id = add) %>%
  dplyr::inner_join(corncyc_pathway) %>%
  dplyr::select(Gene.name) %>%
  dplyr::filter(Gene.name != "unknown") %>%
  dplyr::arrange(Gene.name) %>%
  dplyr::pull(Gene.name) %>%
  unique()

curated_candidates

curated_candidates <- curated_candidates[!is.na(curated_candidates)]
corncyc_classify(curated_candidates) %>% print(n = 60)

gene_models <- union(union(hits$zma00561, hits$zma00561),  drop_transcript_suffix(curated_candidates))
gene_models <- gene_models[!is.na(gene_models)]
 corncyc_pathway

get_unassigned_cyc(
  data.frame( v4_gene_model  = gene_models) %>%
    dplyr::inner_join(corncyc_pathway) %>%
    dplyr::select(Gene.id, Gene.name) %>%
    unique() %>%
    dplyr::pull(Gene.name)
)



get_unassigned_cyc(expanded_cyc$Gene.name)


# Mostly aldehyde dehydrogenases in these orphan enzymes

data.frame( v4_gene_model  = gene_models) %>%
  dplyr::inner_join(corncyc_pathway) %>%
  dplyr::select(Gene.id, Gene.name) %>%
  dplyr::pull(Gene.name) %>%
  unique() %>%
  corncyc_classify() %>% print(n = 60)

gene_models %>% length()


write.csv( data.frame(gene = gene_models),
          file = file.path(
            config$output$dir,"pglipid_candidates.csv"),
          quote = FALSE,
          row.names = FALSE)
