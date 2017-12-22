# read in raw data

get_samples <- function(file) {
  dat <- read_excel(file)
  return(dat)
}

get_sites <- function(file) {
  dat <- read.csv(file)
  return(dat)
}
  