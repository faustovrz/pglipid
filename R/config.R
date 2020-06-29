configure <- function(
  config_file = system.file('config.yaml', package = "pglipid", mustWork= TRUE),
  overwrite = FALSE
){
  config  <- configr::read.config(file = config_file)
  corncyc <- configr::eval.config(file = config_file, config = "corncyc")
  ref     <- configr::eval.config(file = config_file, config = "ref")
  input   <- configr::eval.config(file = config_file, config = "input")
  output  <- configr::eval.config(file = config_file, config = "output")

  usethis::use_data(config, corncyc, ref, input, output,
                    internal = TRUE,
                    overwrite = overwrite)
}


write_config <- function(x,
                         file = system.file('config.yaml', package = "pglipid", mustWork= TRUE)
){
  configr::write.config(x, write.type = "yaml")
}


