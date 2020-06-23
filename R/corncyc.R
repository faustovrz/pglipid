library(Biostrings)
library(dplyr)

config.yaml <- file.path('','Volumes', 'GoogleDrive','My\ Drive','repos','pglipid','extdata','config.yaml')

file.exists(config.yaml)
configr::eval.config.sections(config.yaml)

config <- configr::read.config(file = config.yaml)

corncyc <-  configr::eval.config(
  config = "corncyc",
  file = config.yaml)

cyc_file <- function(x){
  file.path(corncyc$dir,corncyc[x])
}


read_col <- function(input = NULL) {
  if ( cyc_file(input) %>% file.exists()) {
    input <- cyc_file(input)
  }
  read.table(
    input,
    sep = "\t",
    header = TRUE,
    quote = "",
    fill = TRUE,
    na.strings = "")
}


read_dat <- function(input = NULL){

  if ( file.exists(cyc_file(input))) {
    input <- cyc_file(input)
  }

  dat <- readLines(input, encoding = "UTF-8") %>%
    gsub("^ +", "", .,) %>%    # remove leading spaces
    gsub(" +$", "", .,)        # remove trailing spaces

  dat <- dat[!grepl("^#|^/$", dat)]  # remove pound comments,
                                     # and single slash lines
  dat <- paste(
    iconv(dat,
          from="UTF-8",              # changing the damn encoding
          to = "UTF-8"),             # from UTF8 to UTF8 because of BOM?
    collapse = "\n"
  )  %>%
    gsub("\\.\n/?=[^/]",             # Multiple COMMENT values start with /
         "\\.\nCOMMENT - ",
         ., perl = TRUE) %>%
    strsplit("^//\n|\n//\n") %>%     # Chunks split by \n//\n
    unlist()                         # str_split returns a list

  dat <- lapply(dat, FUN = function(x) strsplit(x,"\n") %>% unlist())
  names(dat) <- gsub("UNIQUE-ID - ", "",
                     sapply(dat, "[[", 1))
  dat2 <-list()

  for( name in names(dat)){
    key_val  <- strsplit(dat[[name]], " - ")
    val <-  sapply(key_val, "[", 2)
    names(val) <- sapply(key_val,"[", 1)
    dat2[[name]] <- val
  }

  return(dat2)
}

genes_col <- read_col("genes_col")


pathways_col <- read_col("pathways_col")


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
  )


# genes_dat <- read_dat("genes_dat")

corncyc_gene_synonym <- genes_col %>%
  # pivot synonym
  dplyr::select( UNIQUE.ID,NAME, dplyr::starts_with("SYNONYM")) %>%
  tidyr::pivot_longer(
    cols = dplyr::starts_with("SYNONYM"),
    values_to = "Synonym",
    values_drop_na = TRUE) %>%
  dplyr::select(-name) %>%
  dplyr::rename(Gene.id = UNIQUE.ID) %>%
  dplyr::arrange( Synonym) %>%
  print(n = 200)

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


corncyc_pathway %>%
  dplyr::group_by(Pathway.id) %>%
  dplyr::summarize(n = length(Gene.id)) %>%
  dplyr::arrange(-n)

enzymes_col <- read_col("enzymes_col")

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



proteins_dat <- read_dat("proteins_dat")

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

rownames(enz_rxn) <- NULL

orphan_enz <- enz_rxn %>%
  dplyr::left_join(enz_rxn_path) %>%
  dplyr::group_by(Protein.id) %>%
  dplyr::filter(all(is.na(Pathway.id))) %>%
  unique()

# reactions_dat <- read_dat("reactions_dat")

# reactions_dat <-NULL

pathway_n <- corncyc_pathway %>%
  dplyr::pull(Gene.id) %>%
  unique %>% length()


gene_n <- nrow(genes_col)
# gene_n <- length(names(corncyc_seq))

unassined_n <-  gene_n - pathway_n

colnames(corncyc_pathway)

add_transcript_version<-function(x){
  # Well, there are no other versions other than 1
  # I checked:
  # Zm <- corncyc_pathway$Gene.name[grep("Zm",corncyc_pathway$Gene.name)]
  # Zm.ver <- gsub(".*\\.","",Zm) %>% as.integer()
  # table(Zm.ver)
  # Zm.ver
  # 1
  # 9120
  paste(x,1, sep = '.' )

}

drop_transcript_suffix <- function(x) {
  gsub("_[TP]\\d\\d\\d.*", "", x)
}


test_count <- function(test_genes) {
  gene_ids <- test_genes %>%
              as.data.frame()
  colnames(gene_ids)[1] <- "Gene.name"
  test_count <-  gene_ids %>%
    dplyr::left_join(corncyc_pathway) %>%
    dplyr::group_by(Pathway.id,Pathway.name) %>%
    dplyr::summarise(n_test = length(Gene.name)) %>%
    dplyr::arrange(-n_test)

    test_count[is.na(test_count)] <-"unassigned"
    test_count
    }

cyc_count <- corncyc_pathway %>%
  dplyr::group_by(Pathway.id,Pathway.name) %>%
  dplyr::filter(Gene.name != "unknown")  %>%
  dplyr::summarise(n = length(Gene.name))



cyc_test <- function(test_genes) {
  cyc_test <- cyc_count %>%
    dplyr::right_join(test_count(test_genes)) %>%
    dplyr::mutate(cover = n_test/n) %>%
    dplyr::arrange(-cover) %>%
    as.data.frame()
  # cyc_test <- within(cyc_test,{
  #   n[Pathway.name == "unassigned"] <- unassined_n
  #   cover <- n_test/n
  #   })
}


corncyc_classify <- function(test_genes, bg = NULL){

  sum_test <- length(test_genes)
  cyc_test <- cyc_test(test_genes)
  cyc_genes <-  corncyc_pathway %>%
    dplyr::filter(Pathway.id %in% cyc_test$Pathway.id) %>%
    dplyr::filter(Gene.name != "unknown") %>%
    dplyr::select(Gene.name) %>%
    unique() %>% t() %>% as.vector()
  sum_cyc <- length(cyc_genes)

  if (is.null(bg)){
    bg <- union(test_genes,cyc_genes)
  }

  if (! all(test_genes %in% bg)){
    stop(paste("Test genes absent from background:",
               test_genes[!test_genes %in% bg])
         )
  }

  sum_total <- length(bg)

  fisher <- apply(cyc_test[3:5],1, FUN = function(x){
    if(any(is.na(x))){ return(c(NA,NA)) }
    else {
      m <- c(x["n_test"]            , x["n"] - x["n_test"],
             sum_test - x["n_test"] , sum_total - sum_test +  x["n_test"] - x["n"] ) %>%
        matrix( ncol =2, nrow = 2) %>% t()
      f <- fisher.test(m,
                       alternative="greater")
      c(f$p.value, f$estimate)}
}) %>% t()
colnames(fisher) <- c("p_value", "odds_ratio")

cyc_test <- cbind(cyc_test,fisher)


cyc_test %>%as_tibble() %>%
  dplyr::mutate(FDR = p.adjust(cyc_test$p_value)) %>%
  dplyr::select(Pathway.id:cover,odds_ratio, p_value,FDR) %>%
  dplyr::arrange(FDR,p_value)
}





# legacy functions when I was working with PMN12.5 corncyc.fasta.
# source:  ftp://ftp.plantcyc.org/Pathways/BLAST_sets/PMN_BLAST_archives/PMN12.5/corncyc.fasta
# Now I work with the ensembl cdna fasta file

corncyc_seq <- Biostrings::readAAStringSet(
  filepath = config$pathwayc$fasta
)

corncyc_genes <- data.frame(
  #defline = gsub(" \\| \\S+: | \\| ","\t",names(corncyc_seq), perl =TRUE)
  defline =  names(corncyc_seq)
) %>%
  tidyr::separate(defline,
                  c("fasta.id",
                    "transcript",
                    "description",
                    "Species",
                    "Gene.id"),
                  sep = " \\| \\S+: | \\| ") %>%
  dplyr::mutate(Gene.name = gsub("_\\S\\d{3}\\.\\d+","",transcript, perl =TRUE)) %>%
  dplyr::select(fasta.id,Gene.id,Gene.name, everything())

row.names(corncyc_genes) <- corncyc_genes$fasta.id

# colnames(corncyc_genes)
# corncyc_genes[-5]

pathway_n <- corncyc_pathway %>%
  dplyr::select(Gene.id) %>%
  dplyr::filter(Gene.id != "unknown") %>%
  unique() %>%
  nrow()
