library(dplyr)
config_file <- file.path("..", "inst","config.yaml")
config_R    <- file.path("..", "R", "config.R")
read_col_R  <- file.path("..", "R", "read_col.R")
read_dat_R  <- file.path("..", "R", "read_dat.R")
gene_name_R <- file.path("..", "R", "gene_name.R")
genome_R    <- file.path("..", "R", "genome.R")
make_data_R <- file.path("..", "R", "make_data.R")



source(config_R)
source(read_col_R)
source(read_dat_R)
source(gene_name_R)
source(genome_R)

configure(config_file = config_file, overwrite = TRUE)

source(make_data_R)

corncyc_pathway <- make_corncyc_pathway()
gene_synonym    <- make_corncyc_gene_synonym()
enz_rxn_path    <- make_enz_rxn_path()
enz_rxn         <- make_enz_rxn()
orphan_enz      <- make_orphan_enz( enz_rxn, enz_rxn_path)
map_id          <- make_map_id()

usethis::use_data(
  corncyc_pathway,
  gene_synonym,
  enz_rxn_path,
  enz_rxn,
  orphan_enz,
  map_id,
  internal = FALSE,
  overwrite = TRUE
)

