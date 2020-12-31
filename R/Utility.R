# `ExtractField` function was implemented in `Seurat` 2.0
ExtractField <- function(string, field = 1, delim = "_") {
  fields <- as.numeric(
    x = unlist(x = strsplit(x = as.character(x = field), split = ","))
  )
  if (length(x = fields) == 1) {
    return(strsplit(x = string, split = delim)[[1]][field])
  }
  return(paste(strsplit(x = string, split = delim)[[1]][fields], collapse = delim))
}

CreateResultsFolder <- function(name) {
  assertthat::is.string(name)
  if (name == 'R' || name == 'Data') {
    dir.create(file.path('Results'))
  } else {
    dir.create(file.path(name))
  }
}
