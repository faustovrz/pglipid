config.yaml <- file.path('config.yaml')

system.file('extdata','config.yaml', package = "pglipid", mustWork= TRUE)

file.exists(config.yaml)
configr::eval.config.sections(config.yaml)

config <- configr::read.config(file = config.yaml)

corncyc <-  configr::eval.config(
  config = "corncyc",
  file = config.yaml)

ref <-  configr::eval.config(
  config = "ref",
  file = config.yaml)

input <- configr::eval.config(
  config = "input",
  file = config.yaml)

output <- configr::eval.config(
  config = "output",
  file = config.yaml)

#' @title Make cyc file paths from configuaration keys
#' @description  Make cyc file paths from configuaration keys
#' @details
#' @param input Biocyc col file
#'
#' @export
cyc_file <- function(x){
  file.path(corncyc$dir,corncyc[x])
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

load(file.path('data','corncyc_pathway.rda'))

pathway_n <- corncyc_pathway %>%
  dplyr::pull(Gene.id) %>%
  unique %>% length()

genes_col <- read_col("genes_col")
gene_n <- nrow(genes_col)
# gene_n <- length(names(corncyc_seq))

unassined_n <-  gene_n - pathway_n

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
    dplyr::left_join(id_map) %>%
    dplyr::pull(merge)
}
