library(dplyr)
library(GenomicRanges)
library(regioneR)
library(liftOver)
# 1. Initialize configuration ##################################################
config.yaml <- file.path('extdata','config.yaml')

file.exists(config.yaml)
configr::eval.config.sections(config.yaml)

config <- configr::read.config(file = config.yaml)

ref <-  configr::eval.config(
  config = "ref",
  file = config.yaml)

input <- configr::eval.config(
  config = "input",
  file = config.yaml)

output <- configr::eval.config(
  config = "output",
  file = config.yaml)

analysis <- configr::eval.config(
  config = "analysis",
  file = config.yaml)



################################################################################


genes <- (
  subset(
  get_gene_annot( version = output$PBE$B73_version, organellar = FALSE),
                  type    == "gene" &
                  biotype != "transposable_element" &
                  biotype != "pseudogene" # &
  ) + 10000 ) %>%        # + 10Kb upstream and downstream
  GenomicRanges::trim()


# 1. Get genes corresponding to the 5% outlier SNPs  by population

pbe_files <-  input$PBE$files
pops <- names(input$PBE$files)

names(pbe_files) <- pops

chain_file <-  file.path(ref$chain_file$dir,
                         ref$chain_file$AGPv3_AGPv4)

ch <-rtracklayer::import.chain(chain_file)

pbe_outlier <- data.frame( gene = genes$gene_id)
gene_hits <- data.frame()
outlier_hits <- data.frame()

for (pop in pops) {

  pbe_file <- file.path(input$PBE$dir, pbe_files[pop])

  SNP <- GenomicRanges::makeGRangesFromDataFrame(
      read.delim(pbe_file, header = TRUE),
      keep.extra.columns = TRUE,
      ignore.strand = TRUE,
      seqinfo = get_chr_GR("AGPv3") %>% seqinfo(),
      seqnames.field = "CHROM",
      start.field = "POS",
      end.field = "POS"
      ) %>%
    liftOver(ch) %>%
    unlist()

  seqlevels(SNP) <- as.character(1:10)
  seqinfo(SNP) <-  seqinfo(genes)

  # Find genes overlapping with top 5% SNPs (outliers) in this population
  top <-  analysis$PBE$outlier$top
  thresh <- 1 - top

  outlier_SNP <- SNP[ SNP$PBE0 > quantile(SNP$PBE0, thresh)]
  outlier_olap <- findOverlaps(outlier_SNP, genes)
  with_outlier_SNP <- subjectHits(outlier_olap) %>% unique()

  write.csv( data.frame(
            gene      = genes$gene_id[with_outlier_SNP]),
            file      = file.path(
                         output$dir,
                         paste0(pop, output$PBE$outlier$suffix)
                        ),
            row.names = FALSE,
            quote     = FALSE)

  pbe_outlier[pop] <- 0
  pbe_outlier[with_outlier_SNP,pop] <- 1

  genic_olap <- findOverlaps(SNP,genes)
  with_PBE_SNP <- subjectHits(genic_olap) %>% unique()
  bg_genes <- genes$gene_id[with_PBE_SNP]

  write.csv( data.frame(gene = bg_genes),
             file = file.path(
                      output$dir,
                      paste0(pop, output$PBE$outlier$bg_suffix)
                    ),
             row.names = FALSE,
             quote = FALSE)

  gene_hits <- rbind(
     gene_hits,
     data.frame(
       pop  = pop,
       chr  = seqnames(SNP)[queryHits(genic_olap)],
       pos  = start(SNP)[queryHits(genic_olap)],
       gene = genes$gene_id[subjectHits(genic_olap)]
       ) %>%
       dplyr::inner_join(
         as.data.frame(SNP),
         by=c(chr= "seqnames", pos = "start")
         ) %>%
       arrange(chr,pos)
    )

  outlier_hits <- rbind(
    outlier_hits,
    data.frame( pop     = pop,
                gene    = genes$gene_id[with_outlier_SNP],
                outlier = 1)
    )

}


# Get PBE statistics for genes in these populations----------------------------
# this is really ugly but it  gets the results

outlier_hits$gene %>% unique()

pbe_outlier <- data.frame(gene = genes$gene_id) %>%
  dplyr::left_join(outlier_hits) %>%
  tidyr::pivot_wider(names_from = pop,
                     values_from = outlier,
                     values_fill = 0) %>%
  dplyr::select(-`NA`) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(sum = rowSums(.[2:5])) %>%
  dplyr::left_join(
    gene_hits %>%
      dplyr::group_by(gene) %>%
      dplyr::summarise(max_PBE = max(PBE0),
                       max_PBE_pop = pop[which.max(PBE0)])) %>%
  dplyr::arrange(-sum, -max_PBE) #%>%
  # dplyr::left_join(
  #   as.data.frame(genes) %>%
  #    dplyr::select(gene =gene_id, description)
  #   )

write.csv(pbe_outlier,
          file = file.path(output$dir,
                           output$PBE$outlier$summary),
          quote = FALSE, row.names = FALSE)



quartz()
pbe_outlier[pbe_outlier$sum >0,]
nrow(pbe_outlier)
pbe_outlier%>%
  dplyr::left_join(outlier_hits) %>%
  dplyr::group_by(pop,gene) %>%
  # dplyr::summarise(mean_PBE = mean(PBE), max_PBE = max(PBE),sum = max(sum)) %>%
  ggplot2::ggplot(  ggplot2::aes(x=factor(sum),
                                 y=max_PBE,
                                 color=factor(pop),
                                 shape =max_PBE_pop
                                )
                    ) +
  ggplot2::geom_point( position = ggplot2::position_jitterdodge()) +
  ggplot2::xlab("# Population Selected") +
  ggplot2::ylab("Max PBE") +
  ggpubr::theme_classic2()


#-------------------------------------------------------------------------------
# Upset plot from the outliers

set_m <- pbe_outlier[,c("gene",pops,"sum")] %>%
  dplyr::filter( sum >0) %>%
  dplyr::select(gene:AN) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()

source("R/upset_blue.R")
# n <- length(names(genes))
n <- length(genes)
Result <- SuperExactTest::supertest(SET_input(set_m) ,n=n, degree = 2:4)
str(Result)
pvalues <- Result$P.value
# correct P values to e-320 when p = 0

pdf(file = file.path(
              output$dir,
              output$PBE$outlier$barplot
            )
    )
plot(Result, Layout="landscape", sort.by="size", keep=FALSE)
dev.off()


# Blue plot

# correct P values to e-320 when p = 0
pdf(file = file.path(
    output$dir,
    output$PBE$outlier$bluebar
  )
)

plot_blue(set_m, Result$P.value)
dev.off()



