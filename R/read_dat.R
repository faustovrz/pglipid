#' @title Read Biocyc dat files
#' @description The dat files are delimited by double slashes // each line being a key-value pair
#' @details
#' @param input Biocyc dat file
#'
#' @export
#' @return a list of named vaectors representing key:value pairs
#'
#' @examples
#' \dontrun{
#' read_dat("proteins.dat")
#' }
#

read_dat <- function(input = NULL){
  if ( file.exists(cyc_file(input))) {
    input <- cyc_file(input)
  }

  dat <- readLines(input, encoding = "UTF-8") %>%
    gsub("^ +", "", .,) %>%          # remove leading spaces
    gsub(" +$", "", .,)              # remove trailing spaces

  dat <- dat[!grepl("^#|^/$", dat)]  # remove pound comments,
  # and single slash lines
  dat <- paste(
    iconv(dat,
          from="UTF-8",              # changing the damn encoding
          to = "UTF-8"),             # from UTF8 to UTF8 because of BOM?
    collapse = "\n"
  )  %>%
    gsub("\\.\n/?=[^/]",             # Multiple COMMENT values start with /
         "\\.\nCOMMENT - ",
         ., perl = TRUE) %>%
    strsplit("^//\n|\n//\n") %>%     # Chunks split by \n//\n
    unlist()                         # str_split returns a list

  dat <- lapply(dat, FUN = function(x) strsplit(x,"\n") %>% unlist())
  names(dat) <- gsub("UNIQUE-ID - ", "",
                     sapply(dat, "[[", 1))
  out <- list()

  for( name in names(dat)){
    key_val  <- strsplit(dat[[name]], " - ")  #  " - " delimits key from value
    val <-  sapply(key_val, "[", 2)           #  retrieve values
    names(val) <- sapply(key_val,"[", 1)      #  retrieve keys
    out[[name]] <- val
  }

  return(out)
}
