source("corncyc.R")
library(dplyr)
library(KEGGREST)
library(Biostrings)



# KEGG "Glycerophospholipid metabolism" ----------------------------------------

glpm <- "Glycerophospholipid metabolism"
map <- names(keggFind("pathway","Glycerophospholipid metabolism"))
glpm <- gsub("path:map", "zma:", map)

zma <- keggList("pathway","zma")

glpm_id <- names(zma[grepl("Glycerophospholipid metabolism",zma)])
glpm <- keggGet(glpm_id)
glpm_ngene <- length(glpm[[1]]$GENE)
glpm_entrez <- glpm[[1]]$GENE[seq(1,glpm_ngene,by =2)]

glpm_genes <- paste0("zma:",glpm_entrez)

# get the Sequences
# You can only get 10 at a time!!!!!!!!!!!!!
#
# seqs <- Biostrings::DNAStringSet()
# seqlist <- lapply(chunks,
#   function(x){
#     seqs <- append(seqs, keggGet(x, "ntseq"))
#   }
# )

# seqs <- Biostrings::DNAStringSet()
#
# for(x in names(seqlist)){
#          seqs <- append(seqs,seqlist[[x]])
# }

# Biostrings::writeXStringSet(seqs, filepath = "zma00564.fasta")





# Corncyc classfication 1 ----- KEGG "Glycerophospholipid metabolism" ---------
#
#
# Make blast vs corncycPMN13 to get gene identifiers
#
# blastcmd <- 'blastx -query zma00564.fasta  -db /ref/zea/corncyc/corncycPMN13  -subject_besthit  -num_alignments 1 -outfmt 6 > blast'
# system(blastcmd)
#

# Make blast vs AGPv4 cdna to get gene identifiers
#
#blastcmd <- 'blastn -query /Users/fvrodriguez/Projects/NCSU/05_PLA1/zma00564.fasta  -db /ref/zea/Zea_mays.B73_RefGen_v4.cdna.all  -subject_besthit  -num_alignments 1 -outfmt 6 >blast'
#system(blastcmd)

# add colnames to blast table?
# blast <- read.delim("blast", header = FALSE)

# zma00564 <- data.frame (
#   Gene.name = gsub("_.*","",
#                    blast$V2[blast$V3>95],
#                    perl = TRUE) %>%
#     unique()
# )

# zma00564 <- data.frame (
#   Gene.name = blast$V2[blast$V3>95] %>%  unique()  %>% add_transcript_version()
# )
# nrow(zma00564)

# ID mapping to  AGPv4

xref <- read.table(config$ref$xref, sep = "\t", header = TRUE)

zma00564 <- data.frame (
  Entrez = glpm_entrez %>% as.integer()
)


corncyc_pathway$v4_gene_model <- drop_transcript_suffix(corncyc_pathway$Gene.name)

zma00564 <- zma00564 %>%
  dplyr::left_join(xref) %>%
  dplyr::select( Entrez, v4_gene_model) %>%
  dplyr::left_join(corncyc_pathway) %>%
  dplyr::select(Entrez, v4_gene_model,Gene.id, Gene.name) %>%
  unique()

pathway_cover <- corncyc_classify(
  zma00564 %>%
  dplyr::filter(!is.na(Gene.name)) %>%
  dplyr::pull(Gene.name) %>%
  unique()
  ) %>% print (n=30)


zma00564_expanded_cyc <- data.frame(
  Gene.name = union(
    zma00564$Gene.name,
    pathway_cover %>%
      dplyr::filter(n_test >1) %>%  # remove incidental pathways
      dplyr::select(Pathway.id) %>%
      dplyr::inner_join(corncyc_pathway) %>%
      dplyr::select(Gene.name) %>%
      dplyr::filter(Gene.name != "unknown") %>%
      dplyr::arrange(Gene.name) %>%
      dplyr::pull(Gene.name)
  )
)

zma00564_expanded_cyc %>%
  dplyr::inner_join(genes_col, by = c(Gene.name = "NAME")) %>%
  dplyr::select(Gene.id = UNIQUE.ID, Gene.name) %>%
  dplyr::inner_join(orphan_enz) %>%
  dplyr::select(Gene.id,Gene.name) %>% unique()

nrow(zma00564_expanded_cyc)

corncyc_classify(zma00564_expanded_cyc$Gene.name) %>% print (n=50)

write.csv(zma00564_expanded_cyc,
          file = "zma00564_expanded_cyc.csv",
          quote = FALSE,
          row.names = FALSE)


source("corncyc.R")
library(KEGGREST)
library(Biostrings)

# KEGG "Glycerolipid metabolis" ----------------------------------------

pathway_name<- "Glycerolipid metabolism"
path_map <- names(keggFind("pathway",pathway_name))
code <- gsub("path:map", "zma:", map)
zma_code <- gsub("zma:","zma", code)

zma <- keggList("pathway","zma")

path_code <- names(zma[grepl(pathway_name,zma)])

pathway <- keggGet(path_code)
line_n <- length(pathway[[1]]$GENE)
entrez <- pathway[[1]]$GENE[seq(1,line_n,by =2)]

pathway_genes <- paste0("zma:",entrez)
chunks <- split(pathway_genes, ceiling(seq_along(pathway_genes)/10))

# get the Sequences
# You can only get 10 at a time!!!!!!!!!!!!!

# seqlist <- lapply(chunks,
#   function(x){
#      keggGet(x, "ntseq")
#   }
# )
#
# seqs <- Biostrings::DNAStringSet()
#
# for(x in names(seqlist)){
#          seqs <- append(seqs,seqlist[[x]])
# }
#
# pathway_fasta <- paste0(zma_code,".fasta")
# Biostrings::writeXStringSet(seqs, filepath = pathway_fasta)

# Corncyc classfication 2 ----- KEGG "Glycerolipid metabolism" ---------
#
#
# Make blast vs AGPv4 to get gene identifiers
#

# blastcmd <- 'blastn -query /Users/fvrodriguez/Projects/NCSU/05_PLA1/zma00561.fasta  -db /ref/zea/Zea_mays.B73_RefGen_v4.cdna.all  -subject_besthit  -num_alignments 1 -outfmt 6 >blast'
# system(blastcmd)

# add colnames to blast table?
# blast <- read.delim("blast", header = FALSE)
# zma00561 <- data.frame (
#   Gene.name = blast$V2[blast$V3>95] %>%  unique()  %>% add_transcript_version()
# )
# nrow(zma00561)


# ID mapping to  AGPv4

xref <- read.table(config$ref$xref, sep = "\t", header = TRUE)

zma00561 <- data.frame (
  Entrez = entrez %>% as.integer()
)


# corncyc_pathway$v4_gene_model <- drop_transcript_suffix(corncyc_pathway$Gene.name)

zma00561 <- zma00561 %>%
  dplyr::left_join(xref) %>%
  dplyr::select( Entrez, v4_gene_model) %>%
  dplyr::left_join(corncyc_pathway) %>%
  dplyr::select(Entrez, v4_gene_model,Gene.id, Gene.name) %>%
  unique()

pathway_cover <- corncyc_classify(
  zma00561 %>%
    dplyr::filter(!is.na(Gene.name)) %>%
    dplyr::pull(Gene.name) %>%
    unique()
) %>% print (n=30)


zma00561_expanded_cyc <- data.frame(
  Gene.name = union(
    zma00561$Gene.name,
    pathway_cover %>%
      dplyr::filter(n_test >1) %>%  # remove incidental pathways
      dplyr::select(Pathway.id) %>%
      dplyr::inner_join(corncyc_pathway) %>%
      dplyr::select(Gene.name) %>%
      dplyr::filter(Gene.name != "unknown") %>%
      dplyr::arrange(Gene.name) %>%
      dplyr::pull(Gene.name)
  )
)

zma00561_expanded_cyc %>%
  dplyr::inner_join(genes_col, by = c(Gene.name = "NAME")) %>%
  dplyr::select(Gene.id = UNIQUE.ID, Gene.name) %>%
  dplyr::inner_join(orphan_enz) %>%
  dplyr::select(Gene.id,Gene.name) %>% unique()

nrow(zma00561_expanded_cyc)

corncyc_classify(zma00561_expanded_cyc$Gene.name) %>% print (n=50)

write.csv(zma00564_expanded_cyc,
          file = "zma00564_expanded_cyc.csv",
          quote = FALSE,
          row.names = FALSE)




pathway_cover <- corncyc_classify(zma00561$Gene.name) %>% print(n = 30)



zma00561_expanded_cyc <- data.frame(
  Gene.name = union(
    zma00561$Gene.name,
    pathway_cover %>%
      dplyr::filter(n_test >1) %>%  # remove incidental pathways
      dplyr::select(Pathway.id) %>%
      dplyr::inner_join(corncyc_pathway) %>%
      dplyr::select(Gene.name) %>%
      dplyr::filter(Gene.name != "unknown") %>%
      dplyr::arrange(Gene.name) %>%
      dplyr::pull(Gene.name)
  )
)

corncyc_classify(zma00561_expanded_cyc$Gene.name) %>% print(n =60)


write.csv(zma00561_expanded_cyc,
          file = "zma0056_expanded_cyc.csv",
          quote = FALSE,
          row.names = FALSE)

zma00561_union_zma00564 <-
 data.frame(
   Gene.name = union(
      zma00561_expanded_cyc$Gene.name,
      zma00564_expanded_cyc$Gene.name
      )
  )

nrow(zma00561_union_zma00564)

corncyc_classify(
  zma00561_union_zma00564$Gene.name
) %>%
  print(n = 80)  %>%
  write.csv("zma00561_union_zma00564_expanded_corncyc_cover.csv")

# pathways to remove
remove <- c("PWY-6446","PWY-5143", "PWY66-21","PWY-6457","PWY-2501","PWY-6333","PWY1F-467","PWY0-1313", "PWY-3941")

# pathways to add
add <- c("PWY-6825")

nrow(zma00561_union_zma00564)

zma00561_union_zma00564 %>%
  dplyr::inner_join(genes_col, by = c(Gene.name = "NAME")) %>%
  dplyr::select(Gene.id = UNIQUE.ID, Gene.name) %>%
  dplyr::inner_join(orphan_enz) %>%
  dplyr::select(Gene.id,Gene.name) %>% unique()

write.csv(zma00561_union_zma00564,
          file = "zma00561_union_zma00564_expanded_cyc.csv",
          quote = FALSE,
          row.names = FALSE)


final_set <- corncyc_classify(
  zma00561_union_zma00564$Gene.name
  ) %>%
  dplyr::filter(cover == 1) %>%
  dplyr::filter(!Pathway.id %in% remove) %>%  # remove incidental pathways
  dplyr::select(Pathway.id) %>%
  dplyr::add_row(Pathway.id = add) %>%
  dplyr::inner_join(corncyc_pathway) %>%
  dplyr::select(Gene.name) %>%
  dplyr::filter(Gene.name != "unknown") %>%
  dplyr::arrange(Gene.name) %>%
  dplyr::pull(Gene.name) %>%
  unique()

candidates <- union(final_set,union(zma00561$Gene.name, zma00564$Gene.name))
candidates <- candidates[!is.na(candidates)]
corncyc_classify(candidates) %>% print(n = 60)



