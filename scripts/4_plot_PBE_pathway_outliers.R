source("/Volumes/GoogleDrive/My\ Drive/repos/pglipid/R/corncyc.R")

genes <- subset(annot,
                type == "gene" &
                  biotype != "transposable_element" &
                  biotype != "pseudogene" &
                  seqnames %in% 1:10) + 10000 # + 10Kb upstream and downstream


pops <- c("US","MH","GH","AN")

union_set <- vector()
to_intersect <- list()
to_intersect_bg <- list()
union_bg <- vector()

for (pop in pops) {

  bg_file <- paste0(pop,"_PBE_bg_genes.csv")
  bg <- read.csv(bg_file,col.names = "gene", header = TRUE)
  bg$Gene.name <- get_cyc_id(bg$gene)

  to_intersect_bg[[pop]] <- bg$gene

  outlier_file <- paste0(pop,"_PBE_outlier_genes.csv")
  test_set <- read.csv(outlier_file,
                       col.names = "gene",
                       header = TRUE)  %>%
    dplyr::mutate(Gene.game = get_cyc_id(gene) ) %>%
    dplyr::filter(Gene.game %in%  bg$Gene.name )




  cyc_summary_file <- paste0(pop,"_PBE_outlier_corncyc_summary.csv")
  corncyc_classify(test_set$Gene.game, bg = bg$Gene.name) %>% print(n =20)
  write.csv( corncyc_classify(test_set$Gene.game, bg = bg$Gene.name),
             file = cyc_summary_file)


  union_set <-union(union_set,test_set$gene)
  to_intersect[[pop]] <- test_set$gene

  union_bg <- union(union_bg, bg$gene)
}



intersect_set <- Reduce(intersect, to_intersect )
intersect_set<- intersect_set[intersect_set %in% union_bg]

intersect_bg <- Reduce(intersect, to_intersect_bg )



write.csv( intersect_set,
           file = "intersect_PBE_outlier_genes.csv")

write.csv( corncyc_classify(get_cyc_id(intersect_set), bg = get_cyc_id(intersect_bg)),
           file = "intersect_PBE_outlier_corncyc_summary.csv")

corncyc_classify(get_cyc_id(intersect_set), bg = get_cyc_id(intersect_bg)) %>% dplyr::arrange(-n)



pbe_files <- c("SW_US.allPBE.txt",  "MexHigh.allPBE.txt",
               "GuaHigh.allPBE.txt", "Andes.allPBE.txt")

  colnames(pbe)[1:2]<-c("CHR","BP")
  colnames(pbe)[ncol(pbe)] <- "PBE"

  chain_file <- "/ref/zea/chain_files/AGPv3_to_B73_RefGen_v4.chain"
  ch <-import.chain(chain_file)

  SNP <- makeGRangesFromDataFrame(pbe,
                                  keep.extra.columns=TRUE,
                                  ignore.strand=TRUE,
                                  seqinfo = seqinfo(AGPv4),
                                  seqnames.field="CHR",
                                  start.field="BP",
                                  end.field="BP") %>%
    liftOver(ch) %>% unlist()

  SNP$BP <- start(SNP)
  SNP$CHR <- seqnames(SNP)


  test_genes <- genes[ names(genes) %in% testgeneID]

  # some test genes are not in this genome annotation
  not_annotated <- testgeneID[!testgeneID %in% names(genes)]

  # Find PBE stats for SNPs overlaping test genes

  test_genes_olap <- findOverlaps(SNP, test_genes)

  test_genes_with_SNP <- subjectHits(test_genes_olap) %>% unique()
  in_test_genes <- queryHits(test_genes_olap) %>% unique()

  test_gene_hits <- rbind(test_gene_hits,
                          data.frame(
                            pop = pop,
                            CHR  = SNP$CHR[queryHits(test_genes_olap)],
                            BP   =  SNP$BP[queryHits(test_genes_olap)],
                            gene = names(test_genes)[subjectHits(test_genes_olap)]) %>%
                            dplyr::inner_join(as.data.frame(SNP), by=c("CHR","BP")) %>%
                            arrange(CHR,BP))

  # Find test genes overlapping with top 5% SNPs (outliers) in this population
  top <- 0.05
  thresh <- 1 - top

  outlier_SNP <- SNP[ SNP$PBE > quantile(SNP$PBE, thresh)]
  outlier_olap <- findOverlaps(outlier_SNP, test_genes)

  with_outlier_SNP <- subjectHits(outlier_olap) %>% unique()
  write.csv( data_frame(gene = with_outlier_SNP),
             file = paste0(pop,"_PBE_outlier_genes.csv"),
             row.names = FALSE,
             quote = FALSE)

  outlier_hits <- rbind(outlier_hits,
                        data.frame( pop = pop,
                                    gene    = names(test_genes)[with_outlier_SNP],
                                    outlier = 1))

}

length(test_genes)
length(not_annotated)



test_gene_hits
outlier_hits

# 3. Build the outlier tables as Li Wang


test_genes_PBE <- geneset %>%
  select(gene) %>%
  dplyr::left_join(outlier_hits) %>%
  tidyr::pivot_wider(names_from = pop,
                     values_from = outlier,
                     values_fill = 0) %>%
  dplyr::select(-`NA`) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(sum = rowSums(.[2:5])) %>%
  dplyr::left_join(
    test_gene_hits %>%
      dplyr::group_by(gene) %>%
      dplyr::summarise(max_PBE = max(PBE),
                       max_PBE_pop = pop[which.max(PBE)])) %>%
  dplyr::arrange(-sum, -max_PBE)




write.csv(test_genes_PBE,
          file = "glycerolipid_pathways_PBE_overlap_outlier.csv",
          row.names = FALSE)


shapes <- c("\u25B2","\u25CF","\u25A0","\u25BC")


test_genes_PBE %>%
  dplyr::left_join(test_gene_hits) %>%
  dplyr::group_by(pop,gene) %>%
  dplyr::summarise(mean_PBE = mean(PBE), max_PBE = max(PBE),sum = max(sum)) %>%
  ggplot2::ggplot(aes(x=factor(sum), y=max_PBE, fill=factor(pop))) +
  ggplot2::geom_text(aes(color = factor(pop), label = shapes[factor(pop)]),
                     position = position_jitterdodge()) +
  ggplot2::xlab("# Population Selected") +
  ggplot2::ylab("Max PBE") +
  ggpubr::theme_classic2()


#-------------------------------------------------------------------------------
# Upset plot from the outliers

set_m <- test_genes_PBE[,c("gene",pops,"sum")] %>%
  dplyr::filter( sum >0) %>%
  dplyr::select(gene:AN) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()


# n <- length(names(genes))
n <- length(testgeneID)
Result <- SuperExactTest::supertest(SET_input(set_m) ,n=n, degree = 2:4)

pvalues <- Result$P.value

pdf(file = "Figure2_A_SET_style.pdf")
plot(Result, Layout="landscape", sort.by="size", keep=FALSE)
dev.off()


# Blue plot

pdf(file = "Figure2_A.pdf")
plot_blue(set_m, Result$P.value)
dev.off()





