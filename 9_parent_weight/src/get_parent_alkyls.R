get_parents <- function() {
  parent_compounds <- c('Benz(a)anthracene', 'Chrysene', 'Pyrene', 'Fluoranthene', 
                      'Fluorene', 'Naphthalene', 'Anthracene', 'Phenanthrene')
  return(parent_compounds)
}

get_alkyls <- function() {
  cs <- c('C1', 'C2', 'C3', 'C4')
  alkylated_compounds <- c(paste0(cs, '-Chrysenes'), paste0(cs[1:3], '-Fluoranthenes/Pyrenes'), 
                           paste0(cs[1:3], '-Fluorenes'), paste0(cs, '-Naphthalenes'),
                           paste0(cs, '-Phenanthrenes/Anthracenes'))
  return(alkylated_compounds)
}


