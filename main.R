library('rmarkdown')
source('scripts/preprocessing.R')

## List of reports to knit
file_list = c(
  "import-and-preprocess.Rmd", ## Cleanse data
  "diversity-calc.Rmd",        ## Run some prelim ecological calculations
  "beta-diversity.Rmd",        ## Covariate analysis via beta diversity
  "differential-abundance.Rmd" ## Differential abundance analysis
)

## Output directory
output_dir = "reports"
if (!dir.exists(output_dir)) dir.create(output_dir);

## Knit files.
silence({
lapply(file_list, function(rmd_file)
{
	output_file = file.path(output_dir, tools::file_path_sans_ext(basename(rmd_file)))
	render(rmd_file, output_file = paste0(output_file, ".html"))
})})