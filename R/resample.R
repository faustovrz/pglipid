library(dplyr)
library(GenomicRanges)
library(regioneR) 
library(ggplot2)
library(ggpubr)

# Sets of functions for resampling and plotting gene set statistics
# with respect to genic or genomic background

resample <- function(SNP,             # GRanges object
                     gene_annotation, # GRanges object
                     stat = "P",      # Dataframe column to resample
                     sampling = "nr", # "nr": Non redundant  "r" :redundant 
                     bg = "genomic",  # background "genomic" or "genic"
                     test = NULL,     # test gene ids from gene_annotation
                     reps = 10000,    # number of replicates
                     name = NULL){    # replicate name or identifier
  
  if(is.null(test)) {
    stop("No gene test subset provided")
  }
  
  # Find genic SNPs for the background null distribution
  
  genic <- findOverlaps(SNP,gene_annotation) %>% queryHits() %>% unique()
  
  test_genes <- gene_annotation[ names(gene_annotation) %in% test]
  
  # some test genes are not in this genome annotation
  not_annotated <- test[!test %in% names(gene_annotation)]
  
  # Find SNPs in test genes
  
  test_genes_olap <- findOverlaps(SNP, test_genes)
  
  hit_SNP_idx <- queryHits(test_genes_olap)
  hit_gene_idx <- subjectHits(test_genes_olap)
  
  test_genes_with_SNP <- hit_gene_idx %>% unique()
  in_test_genes <- hit_SNP_idx %>% unique()
  
  # hits <- data.frame( 
  #   chr  = seqnames(SNP)[queryHits(test_genes_olap)],
  #   pos   =  start(SNP)[queryHits(test_genes_olap)],
  #   gene = names(test_genes)[subjectHits(test_genes_olap)]) %>% 
  #   
  #   dplyr::inner_join(as.data.frame(SNP), by=c(chr = "seqnames", pos = "start")) %>% 
  #   dplyr::arrange(chr, pos)
  
  if (sampling == "nr"){
    # Non redundant sampling 
    # A particular SNP appears a single time in the test set
    test_stat <-  mcols(SNP)[in_test_genes,stat]
    names(test_stat) <- names(SNP[in_test_genes,])
  }else if(sampling == "r"){
    # Redundant sampling 
    # A particular SNP appears as many times as it overlaps genes
    # test_stat <-  hits[,stat]
    test_stat <-  as.data.frame(SNP)[hit_SNP_idx,stat]
    # names(test_stat) <-  hits[,"gene"]
    names(test_stat) <-  as.data.frame(SNP)[hit_SNP_idx,"gene"]
  }
  
  
  if (length(bg) == 1){
    if ( bg == "genomic" ){
      ## Genomic background, i.e. all SNPs genic and intergenic
      bg <- mcols(SNP)[,stat]
      bg_type ="genomic"
    }else if(bg == "genic"){
      bg <- mcols(SNP)[genic,stat]
      bg_type ="genic"
    }
  }
  
  n <- length(test_stat)
  
  # Null distribution
  nd <- replicate(reps, mean(sample(bg, n, replace = FALSE)))
  
  test_mean <- mean(test_stat)
  
  pvalue <- 1 - ecdf(nd)(test_mean) 

  list( name = name,
        not_annotated = not_annotated,
        sampling = sampling,
        test_genes_with_SNP = test_genes_with_SNP,
        in_test_genes =  in_test_genes,
        n = n,
        test_genes = names(test_genes),
        nd = nd,                         # Null distribution of the statistic
        nd_mean = mean(nd),              # parameters
        bg_n = length(bg),               # for building htitle for plot
        bg_type = bg_type,               # or table
        test_mean = test_mean,           #
        pvalue  = pvalue)                #
  
}



negLogP <- function(x){expression( -"Log"[10] ( italic(P) ))}




pvalue_string <- function(x, n = 10000){
  if(x == 0){
    paste("p <", formatC(1/n, format = "G", digits = 2))
  } else{
    paste("p =", formatC(x, format = "G", digits = 2))
  }
}

hist_title <- function(x, verbose = TRUE){
  if (verbose == TRUE){
      with(x, paste( pvalue_string(pvalue), "\n",
                "bg:",bg_n, bg_type, "SNPs\n",
                "test:", n, sampling, "SNPs",
                "from", length(test_genes_with_SNP), "genes")
      )
    } else {
      with(x, pvalue_string(pvalue))
    }
}

plot_resample <-  function(resample,
                           verbose = TRUE,
                           title = NULL, ...){
  if (is.null(title)){
  title <- hist_title(resample, verbose = verbose)
  }
  data.frame(nd = resample$nd)  %>% 
    ggplot2::ggplot(aes(x = nd)) +
    ggplot2::ggtitle(title) +
    ggplot2::ylab("Frequency") +
    ggplot2::geom_histogram( col = "black",
                             bins = 50,
                             fill = "black") +
    ggplot2::geom_vline(xintercept = resample$nd_mean, 
                        col = "darkgreen", lwd = 1) +
    ggplot2::geom_vline(xintercept =  resample$test_mean,
                        col = "red", lwd =1 ) +
    ggpubr::theme_classic2(base_size = 12) +
    ggplot2::theme(
      plot.title = element_text(
        face  = "bold",
        hjust = 0.5)
    )
}

