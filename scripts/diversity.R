library('data.table')
library('magrittr')
library('DirichletMultinomial')
library('BiocGenerics')

## Make a matrix compositional across either columns (margin = 2) or rows (margin = 1)
make.compositional = function(matrix, margin = 2) apply(matrix, margin, function(x) x/sum(x));

make.uncompositional = function(matrix)
{
  ## Convert back to integer counts by multiplying by the log10-exponent
  ## inverse to the exponent of the minimum value (i.e.: if the min value
  ## is 1e-6, multiply the entire dataset by 1e6)
  mult_fac = 10^-(log10(matrix[matrix!=0]) %>% min())
  int_matrix = (matrix * mult_fac)
  
  ## Get rownames
  rownames(int_matrix) = rownames(matrix)
  colnames(int_matrix) = colnames(matrix)
  
  return(int_matrix)
}

## Run effective microbial richness per Reitmer et al. 2021
## https://www.nature.com/articles/s43705-021-00033-z
## min.abund is the minimum threshold abundance (default is 0.25%)
## matrix is the matrix of observations
## optional arguments are passed to make.compositonal
EMR = function(matrix, min.abund = 0.25 * 1/100, ...)
{
  
  .emr = function(x) sum(x > min.abund)
  compositional = make.compositional(matrix, ...)
  
  emr_out = apply(compositional, 2, .emr)
  
  return(emr_out)
}

## Run DMMs.
## Arguments are:
## matrix = a community matrix with samples-as-columns.
## nmaxdmms = (optional) a number of dmms to maximally assign
## (defaults to sqrt(1/2 sample_count, minimum 4)
## Returns a list containing:
## A named vector indicating group associations according to best-laplace
## The index of the fit object with the best-laplace
## The returned fit object.
generate.dmms = function(matrix, nmaxdmms = NULL)
{
  int_matrix = copy(matrix)
  int_matrix = make.uncompositional(int_matrix)
  int_matrix = int_matrix %>% apply(2, as.integer)
  
  ## Transpose; dmm expects species-as-columns
  int_matrix = t(int_matrix)
  colnames(int_matrix) = rownames(matrix)
  
  ## Default to sqrt (1/2); this is a good guess based on my experience
  if(is.null(nmaxdmms)) 
  {
    nmaxdmms = ceiling(sqrt(ncol(matrix) / 2 ));
    nmaxdmms = max(nmaxdmms, 4)
  }
  dmm_fit = lapply(2:nmaxdmms, dmn, count = int_matrix, verbose = FALSE)
  
  ## Return laplace
  laplace_scores = dmm_fit %>% vapply(laplace, numeric(1))
  min_laplace = which.min(laplace_scores)
  
  ## Get best fit
  best_fit = dmm_fit[[min_laplace]]
  ndmms = min_laplace + 1
  best_fit_probs = mixture(best_fit)
  best_fit_quanta = 1 - round(1 - best_fit_probs)
  
  colnames(best_fit_quanta) = 1:ndmms %>% paste0('DMM',.)

  dmm_melt = reshape2::melt(best_fit_quanta)
  dmm_melt2 = dmm_melt[dmm_melt$value != 0,]
  
  dmms = setNames(dmm_melt2[,'Var2'], dmm_melt2[,'Var1'])
  
  return(list(dmms, min_laplace, dmm_fit))
}
