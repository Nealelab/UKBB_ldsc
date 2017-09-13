###
# run rmarkdown reports for each phenotype
###
# setwd("/Users/raymondwalters/Documents/Code/github/ukbb_ldsc_pages/")

# load packages
rm(list=ls())
require("rmarkdown")

# load data
load("ukbb_h2part.RData")

# loop
for(phen in dat$phenotype[c(1:5,sample(1:nrow(dat),5))]){
	render('h2_part_template.Rmd',
		   output_file=paste0('h2_summary_',phen,'.html'),
		   output_dir='./docs')
}

