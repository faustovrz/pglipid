library(Biostrings)
library(dplyr)

config.yaml <- file.path('extdata','config.yaml')

configr::eval.config.sections(config.yaml)

config <- configr::read.config(file = config.yaml)

read_cyc_col <- function(x){
  read.table(
    x,
    sep = "\t",
    header = TRUE,
    quote = "",
    fill = TRUE,
    na.strings = "")
}

path_cyc_wide <- read_cyc_col(config$corncyc$pathway)


corncyc_pathway <- cbind(
  # pivot GENE.ID
  path_cyc_wide %>%
  dplyr::select(UNIQUE.ID, NAME, starts_with("GENE.ID")) %>%
  tidyr::pivot_longer(
    cols = starts_with("GENE.ID"),
    values_to = "Gene.id",
    values_drop_na = TRUE
  ) %>%
  dplyr::select(-name),

  # pivot GENE.NAME
  path_cyc_wide %>%
    dplyr::select(UNIQUE.ID, starts_with("GENE.NAME")) %>%
    tidyr::pivot_longer(
      cols = starts_with("GENE.NAME"),
      values_to = "Gene.name",
      values_drop_na = TRUE) %>%
    dplyr::select(-UNIQUE.ID,-name)
  ) %>%
  dplyr::rename(Pathway.id = UNIQUE.ID, Pathway.name = NAME)

corncyc_genes <- read_cyc_col(config$corncyc$genes) %>%
  dplyr::rename(Gene.id = UNIQUE.ID, Gene.name = NAME)

corncyc_gene_synonym <- corncyc_genes %>%
  dplyr::select( Gene.name, dplyr::starts_with("SYNONYM")) %>%
  tidyr::pivot_longer(
    cols = dplyr::starts_with("SYNONYM"),
    values_to = "Synonym",
    values_drop_na = TRUE) %>%
  dplyr::select(-name) %>%
  dplyr::arrange(Synonym) %>%
  print(n = 200)


corncyc_enzymes <- read_cyc_col(config$corncyc$enzymes)



corncyc_pathway %>%
  dplyr::group_by(Pathway.id) %>%
  dplyr::summarize(n = length(Gene.id)) %>%
  dplyr::arrange(-n)


corncyc_seq <- Biostrings::readAAStringSet(
  filepath = config$pathwayc$fasta
  )
corncyc_pathway[50:55,]

cyc_enz_sub <- corncyc_enzymes %>%
  dplyr::select(UNIQUE.ID,Subunit = SUBUNIT.COMPOSITION) %>%
  dplyr::mutate(Subunit = gsub("\\d\\*|,\\d\\*","", Subunit)) %>%
  tidyr::separate_rows(Subunit, sep = ",",  convert = TRUE)

enz_rxn_path <- corncyc_enzymes %>%
  dplyr::select( ENZRXN = UNIQUE.ID, dplyr::starts_with("PATHWAYS")) %>%
  tidyr::pivot_longer(
    cols = dplyr::starts_with("PATHWAYS"),
    values_to = "Pathway.id",
    values_drop_na = TRUE) %>%
  dplyr::select(-name) %>%
  dplyr::left_join(
    corncyc_enzymes %>%
      dplyr::select(ENZRXN = UNIQUE.ID, Protein.id = SUBUNIT.COMPOSITION) %>%
      dplyr::mutate(Protein.id = gsub("\\d\\*|,\\d\\*","", Protein.id)) %>%
      tidyr::separate_rows(Protein.id, sep = ",",  convert = TRUE)
  ) %>%
  dplyr::mutate(
    Gene.id = gsub("-MONOMER","",Protein.id)
  )

invalid_utf8_ <- function(x){

  !is.na(x) & is.na(iconv(x, "UTF-8", "UTF-8"))

}

detect_invalid_utf8 <- function(string, seperator){

  stringSplit <- unlist(strsplit(string, seperator))

  invalidIndex <- unlist(lapply(stringSplit, invalid_utf8_))

  data.frame(
    word = stringSplit[invalidIndex],
    stringIndex = which(invalidIndex == TRUE)
  )

}

gsub("red","blue", s)
read_dat <- function(x){
  dat <- NULL
  dat <- readLines(x, encoding = "UTF-8") %>%
          gsub("^ +", "", .,) %>%    # remove leading spaces
          gsub(" +$", "", .,)        # remove trailing spaces

  dat <- dat[!grepl("^#|^/$", dat)]  # remove pound comments,
                                     # and single slash lines
  dat <- paste(
    iconv(dat, from="UTF-8",         # changing the damn encoding
          to = "UTF-8"),             # from UTF8 to UTF8 because BOM?
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

  for( p in names(dat)){
    key_val  <- strsplit(dat[[p]], " - ")
    val <-  sapply(key_val, "[", 2)
    names(val) <- sapply(key_val,"[", 1)
    dat2[[p]] <- val
  }

  return(dat2)
}

protein_dat <- read_dat(config$corncyc$proteins_dat)

enz_rxn <-NULL

enz_rxn <- lapply(protein_dat, function(x){

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

enz_rxn %>%
  dplyr::left_join(enz_rxn_path)

reactions_dat <- read_dat(config$corncyc$reactions_dat)


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

gene_n <- nrow(corncyc_genes)

gene_n <- length(names(corncyc_seq))

unassined_n <-  gene_n - pathway_n




test_count <- function(test_genes) {
  gene_ids <- as.data.frame(test_genes)
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

