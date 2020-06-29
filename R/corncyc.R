# The configuration is stored at "sysdata.rda"

internal_data <- file.path("R","sysdata.rda")
if(file.exists(internal_data)){
  load(internal_data)
}

internal_data <- file.path("..","R","sysdata.rda")
if(file.exists(internal_data)){
  load(internal_data)
}

#' @title Read Biocyc col report tables
#' @description The col files are tab delimited files representing info in data files
#' @details
#' @param input Biocyc col file
#'
#' @export
#' @return a dataframe. Note that many column names are serialized because
#'                      they are different values for a common key in dat files
#'
#' @examples
#' \dontrun{
#' read_col("proteins.dat")
#' }
#'
read_col <- function(input = NULL) {
  read.table(
    input,
    sep = "\t",
    header = TRUE,
    quote = "",
    fill = TRUE,
    na.strings = "")
}




genes_col <- read_col(
  file.path(
    config$corncyc$dir,
    config$corncyc$genes_col
  )
)



cyc_test <- function(test_genes) {

  cyc_count <- pglipid::corncyc_pathway %>%
    dplyr::group_by(Pathway.id,Pathway.name) %>%
    dplyr::filter(Gene.name != "unknown")  %>%
    dplyr::summarise(n = length(Gene.name))

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


test_count <- function(test_genes) {
  gene_ids <- test_genes %>%
    as.data.frame()
  colnames(gene_ids)[1] <- "Gene.name"
  test_count <-  gene_ids %>%
    dplyr::left_join(pglipid::corncyc_pathway) %>%
    dplyr::group_by(Pathway.id, Pathway.name) %>%
    dplyr::summarise(n_test = length(Gene.name)) %>%
    dplyr::arrange(-n_test)

  test_count[is.na(test_count)] <-"unassigned"
  test_count
}


corncyc_classify <- function(test_genes, bg = NULL){

  pathway_n <- pglipid::corncyc_pathway %>%
    dplyr::pull(Gene.id) %>%
    unique %>% length()

  gene_n <- nrow(genes_col)
  # gene_n <- length(names(corncyc_seq))

  unassined_n <-  gene_n - pathway_n

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
  if(is.numeric(bg) & length(bg) == 1){
    sum_total <- bg
  }else{
    if (! all(test_genes %in% bg)){
      stop(paste("Test genes absent from background:",
                 test_genes[!test_genes %in% bg])
      )
    }
    sum_total <- length(bg)
  }


  fisher <- apply(cyc_test[3:5],1, FUN = function(x){
    if(any(is.na(x))){ return(c(NA,NA)) }
    else {
      m <- c(x["n_test"]            , x["n"] - x["n_test"],
             sum_test - x["n_test"] , sum_total - sum_test +  x["n_test"] - x["n"] ) %>%
        matrix( ncol =2, nrow = 2) %>% t()
      f <- fisher.test(m, alternative="greater")
      c(f$p.value, f$estimate)}
  }) %>% t()
  colnames(fisher) <- c("p_value", "odds_ratio")

  cyc_test <- cbind(cyc_test,fisher)


  cyc_test %>%as_tibble() %>%
    dplyr::mutate(FDR = p.adjust(cyc_test$p_value)) %>%
    dplyr::select(Pathway.id:cover,odds_ratio, p_value,FDR) %>%
    dplyr::arrange(FDR,p_value)
}



get_cyc_id <- function(gene_id){
  data.frame(gene_id = gene_id) %>%
    dplyr::left_join(pglipid::id_map) %>%
    dplyr::pull(merge)
}
