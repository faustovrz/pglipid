# The configuration is stored at "sysdata.rda"

# install directory from R studio: /R/library/path/pglipid
internal_data <- file.path("R","sysdata.rda")
if(file.exists(internal_data)){
  load(internal_data)
}

# install directory from R studio: /source/path/pglipid
internal_data <- file.path("..","R","sysdata.rda")
if(file.exists(internal_data)){
  load(internal_data)
}

print("pathways_col not found!")

make_pathways_col <- function()(
  if(!exists("pathways_col")) {
    print("pathways_col not found!")
    return(
      read_col(
        file.path(config$corncyc$dir,
                  config$corncyc$pathways_col)
      )
    )
  }else{
    return(pathways_col)
  }

)

print("pathways_col not found 2!")

pathways_col <-  make_pathways_col()


make_genes_col <- function(ref)(
  if(!exists("genes_col")) {
    return(
      read_col(
        file.path(config$corncyc$dir,
                  config$corncyc$genes_col)
      )
    )
  }else{
    return(genes_col)
  }

)

genes_col <- make_genes_col()

genes_col$v4_gene_model <- drop_transcript_suffix(genes_col$NAME)



make_xref <- function()(
  if(!exists("xref")){
    return(
    read.table(ref$AGPv4$xref$file,
                       sep = "\t", header = TRUE, na.strings = "")
    )
  }

)

make_corncyc_pathway <- function(overwrite = FALSE){

  corncyc_pathway <-   pathways_col %>%
    # pivot GENE.ID
    dplyr::select(UNIQUE.ID, NAME, starts_with("GENE.ID")) %>%
    tidyr::pivot_longer(
      cols = starts_with("GENE.ID"),
      values_to = "Gene.id",
      values_drop_na = TRUE
    ) %>%
    dplyr::select(-name) %>%
    dplyr::rename(Pathway.id = UNIQUE.ID, Pathway.name = NAME) %>%
    dplyr::left_join(
      genes_col %>%
        dplyr::select(Gene.id = UNIQUE.ID, Gene.name = NAME)
    ) %>%
    dplyr::mutate( v4_gene_model = drop_transcript_suffix(Gene.name))

  return(corncyc_pathway)
}


# genes_dat <- read_dat("genes_dat")

make_corncyc_gene_synonym <- function(){
  # pin1 has no associated transcript!!!!!
  #
  # genes_dat[["GDQC-114725"]]
  #
  #                                                             COMMON-NAME
  #                                                                  "pin1"
  #                                                             ACCESSION-1
  #                                                                   "NIL"
  #                                                             ACCESSION-2
  #                                                                   "NIL"
  #                                                                 DBLINKS
  # "(ENSEMBL-PROTEIN \"ARRAY(0x17998d8)\" NIL |zhang| 3701270650 NIL NIL)"
  #                                                                 DBLINKS
  #        "(MAIZEGDB \"ARRAY(0x27fd8e8)\" NIL |zhang| 3701270582 NIL NIL)"
  #
  # But I do know that it is  Zm00001d044812
  # This is an error in the database curation
  # I do not know if it is worth trying to curate the table or the dat file
  # probably in version 9 this is fixed

  # with(genes_col, genes_col[UNIQUE.ID == "GDQC-114725",])
  # with(genes_col, genes_col[grepl("Zm00001d044812", NAME),])
  # sum(grepl("Zm00001d044812",genes_col$NAME))

  corncyc_gene_synonym <- genes_col %>%
    # pivot synonym
    dplyr::select( UNIQUE.ID,NAME, dplyr::starts_with("SYNONYM")) %>%
    tidyr::pivot_longer(
      cols = dplyr::starts_with("SYNONYM"),
      values_to = "Synonym",
      values_drop_na = TRUE) %>%
    dplyr::select(-name) %>%
    dplyr::rename(Gene.id = UNIQUE.ID) %>%
    dplyr::arrange( Synonym)

    return(corncyc_gene_synonym)
}





# corncyc_pathway %>%
#   dplyr::group_by(Pathway.id) %>%
#   dplyr::summarize(n = length(Gene.id)) %>%
#   dplyr::arrange(-n)



make_enz_rxn_path  <- function(){
  enzymes_col <- read_col(
    file.path(config$corncyc$dir,
              config$corncyc$enzymes_col)
  )

  cyc_enz_sub <- enzymes_col %>%
    dplyr::select(UNIQUE.ID,Subunit = SUBUNIT.COMPOSITION) %>%
    dplyr::mutate(Subunit = gsub("\\d\\*|,\\d\\*","", Subunit)) %>%
    tidyr::separate_rows(Subunit, sep = ",",  convert = TRUE)

  enz_rxn_path <- enzymes_col %>%
    dplyr::select( ENZRXN = UNIQUE.ID, dplyr::starts_with("PATHWAYS")) %>%
    tidyr::pivot_longer(
      cols = dplyr::starts_with("PATHWAYS"),
      values_to = "Pathway.id",
      values_drop_na = TRUE) %>%
    dplyr::select(-name) %>%
    dplyr::left_join(
      enzymes_col %>%
        dplyr::select(ENZRXN = UNIQUE.ID, Protein.id = SUBUNIT.COMPOSITION) %>%
        dplyr::mutate(Protein.id = gsub("\\d\\*|,\\d\\*","", Protein.id)) %>%
        tidyr::separate_rows(Protein.id, sep = ",",  convert = TRUE)
    ) %>%
    dplyr::mutate(
      Gene.id = gsub("-MONOMER","",Protein.id)
    )

 return(enz_rxn_path)
}



make_enz_rxn <- function(overwrite = FALSE){
  proteins_dat <- read_dat(
    file.path(config$corncyc$dir,
              config$corncyc$proteins_dat)
  )

  enz_rxn <-NULL

  enz_rxn <- lapply(proteins_dat, function(x){

    catalyzes <- x[names(x) == "CATALYZES"]
    if(length(catalyzes) == 0){
      ENZRXN <- NA
    } else{
      ENZRXN <- catalyzes
    }

    genes <- x[names(x) == "GENE"]
    if(length(genes) == 0){
      GENE <- NA
    } else{
      GENE <- genes
    }

    data.frame(
      Protein.id = x["UNIQUE-ID"],
      Gene.id    = GENE,
      ENZRXN     = ENZRXN
    )
  }
  ) %>% dplyr::bind_rows()

  return(enz_rxn )
}




make_orphan_enz <- function(enz_rxn = NULL, enz_rxn_path = NULL){
  if(is.null(enz_rxn) | is.null(enz_rxn)){
    stop("You must provide tables")
  }
  orphan_enz <- enz_rxn %>%
    dplyr::left_join(enz_rxn_path) %>%
    dplyr::group_by(Protein.id) %>%
    dplyr::filter(all(is.na(Pathway.id))) %>%
    unique()

  return(orphan_enz)
}


# reactions_dat <- read_dat("reactions_dat")

# reactions_dat <-NULL



make_map_id <- function(version = "AGPv4"){

  transcript <- get_gene_annot( version = version,
                                organellar = FALSE) %>%
    subset(type == "mRNA")

  print("Mapping transcript ids (Zm...) to corncyc Zm..._TXXX.1) ")

  if(is.null(transcript)){
    stop(c("Transcript is ",NULL))
  }
  # Without rownames as.dataframe() will fail

  names(transcript) <- transcript$transcript_id

  transcript$gene_id <- drop_transcript_suffix(transcript$transcript_id)

  map_cDNA_id <- transcript %>%
    as.data.frame() %>%
    dplyr::group_by(gene_id) %>%
    dplyr::select(gene_id,transcript_id)

  cyc_transcript <- genes_col%>%
    dplyr::select(Gene.name = NAME) %>%
    dplyr::mutate(
      gene_id = drop_transcript_suffix(Gene.name),
      transcript_id = gsub("\\..*","",Gene.name)) %>%
    dplyr::filter(grepl("^Zm",.[,"Gene.name"]))

  id_map <- map_cDNA_id  %>%
    dplyr::left_join(cyc_transcript) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::arrange(Gene.name, transcript_id) %>%
    dplyr::slice(1)

  id_map$merge <- paste0(id_map$transcript_id,".1")
  id_map$merge[!is.na(id_map$Gene.name)] <- id_map$Gene.name[!is.na(id_map$Gene.name)]

  return(id_map)
}




