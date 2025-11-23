
## These functions override library function to silence all package messages.
## I keep them at the b
silence = function(x) suppressPackageStartupMessages(suppressWarnings({invisible({capture.output(x)})}))
library = function(...) silence({base::library(...)})

library('data.table')
library('reshape2')
library('dqrng')
library('compositions')

## Rarefy to even library size for some downstream operations
## This is necessary to assess diversity.
rarefy = function(x, n)
  ## Rarefy function for large sample-sizes. In essence, this function is just the
  ## fastest interpretation of the following operations:
  ## 1. Convert from wide to long
  ## 2. Expand the long population matrix into a data.table with each row
  ##    representing an individual
  ## 3. Subsample rows without replacement.
  ## 4. Tabulate the subsampled data.
  ## 5. Convert from long back to wide, add zeroes back into the matrix & return.
  ## Params:
  ## x = a count matrix, where columns are samples, to rarefy
  ## n = the threshold to rarefy the counts.
  ## Returns:
  ## A rarefied matrix where columns are samples.
{
  paste('subsampling to', n, 'otus...') %>% message()
  
  ## Wide to long.
  mdt = reshape2::melt(x[, colSums(x) >= n]) %>% as.data.table()
  
  ## Drop zero-value OTUs
  mdt = mdt[value > 0,]
  
  ## Subsample without replacement
  subsam = mdt[,rep(Var1, value) 
               %>% dqsample(size = n, replace = FALSE) 
               %>% table() 
               %>% as.data.table()
               , by = Var2]
  
  ## Change colnames
  colnames(subsam) = c('sample', 'otu', 'subsamn')
  
  ## Long back to wide
  outdt = data.table::dcast(subsam
                            , formula = otu ~ sample
                            , value.var = 'subsamn')
  
  ## Cast as matrix, removing first column (the ids).
  outmat = outdt[,-1] %>% as.matrix()
  rownames(outmat) = outdt[,otu]
  
  ## Clean up NAs and remove zero-count rows.
  outmat[is.na(outmat)] = 0
  outmat = outmat[rowSums(outmat) > 0,]
  
  paste(ncol(x) - ncol(outmat), 'samples were removed.') %>% message()
  paste(nrow(x) - nrow(outmat), 'otus were removed.') %>% message()
  
  return(outmat)
}

agglomerate = function(counts, taxonomy, rank)
  ## Agglomerate function for large sample-sizes.
  ## 1. Convert from wide to long
  ## 2. Bind to taxonomy
  ## 3. Sum counts within selected rank
  ## 4. Convert from long back to wide, add zeroes back into the matrix & return.
  ## Params:
  ## counts = a count matrix, where columns are samples, to agglomerate.
  ## taxonomy = a matrix or data frame where rows are OTUs, representing taxonomy
  ## rank = One of the columns of the taxonomy matrix (a string) to agglomerate on
  ## Returns:
  ## A count matrix where columns are samples and rows are the rank of interest.
{
  mcounts = reshape2::melt(counts) %>% as.data.table() %>% setkey(Var1)
  mtax = taxonomy[mcounts,][,.(value = sum(value)), by = .(Var2, R = get(rank))]
  
  dtax = data.table::dcast(mtax, R ~ Var2, value.var = 'value')
  aggmatrix = dtax[,-1] %>% as.matrix()
  rownames(aggmatrix) = dtax[,R]
  
  return(aggmatrix)
}

## Prevalence filter.
## Params:
## counts = a count matrix where samples are columns.
## min_prevalence = minimum prevalence threshold (as a decimal).
## Returns:
# A filtered count matrix with rows in fewer samples than the threshold of
# min_prevalence removed.
prevalence = function(counts, min_prevalence)
  counts[rowSums(counts > 0) > (ncol(counts)*min_prevalence),];

## Abundance filter.
## Params:
## counts = a count matrix where samples are columns.
## min_abundance = minimum abundance threshold (as a decimal).
## Returns:
## A filtered count matrix with samples at a lower abundance in any sample than
## the threshold removed.
abundance = function(counts, min_abundance)
  counts[(rowSums(apply(counts, 2, function(x) x/sum(x)) >= min_abundance) > 0) %>% which(),];

## CLR transform.
## Params:
## counts = a count matrix where samples are rows.
## pseudocount = pseudocount to be added.
## Returns:
## CLR-transformed matrix, with column names reapplied.
clrtrans = function(x, pseudocount)
{
  tr = compositions::clr(x + pseudocount) %>% as.matrix()
  colnames(tr) = colnames(x)
  return(tr)
}


## data.table to a matrix
## Params: x = a data table; typefun = cast entries as this type
## Returns: a matrix of type where the rownames are the first column of the dt.
dt2mat = function(x, typefun = NULL)
{
  if(is.null(typefun))
  {
    ty = x[1,2] %>% unlist() %>% class()
    typefun = paste0('as.', ty) %>% get()
  }
  
  mat = x[,-1] %>% apply(2, typefun) %>% as.matrix()
  rownames(mat) = x[,1] %>% unlist() %>% as.character()
  return(mat)
}

## data.table to a data.frame
## Params: x = a data table
## Returns: a data.frame where the rownames are the first column of the dt.
dt2df = function(x)
{
  df = x[,-1] %>% as.data.frame()
  rownames(df) = x[,1] %>% unlist() %>% as.character()
  return(df)
}