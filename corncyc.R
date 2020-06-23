library(Biostrings)
library(dplyr)

corncyc_pathway <- read.delim("/ref/zea/corncyc/corncyc_pathways.20180702", header = TRUE)
corncyc_seq <- Biostrings::readAAStringSet(filepath = "/ref/zea/corncyc/corncycPMN13.fasta")

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

