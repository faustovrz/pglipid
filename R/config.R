get_config<- function(
  config_file = system.file('config.yaml',
                            package = "pglipid",
                            mustWork= TRUE)
){
  config  <- configr::read.config(file = config_file)
  return(config)
}


write_config <- function(x,
                file = system.file('config.yaml',
                                   package = "pglipid",
                                   mustWork= TRUE)
){
  configr::write.config(x, write.type = "yaml")
}


