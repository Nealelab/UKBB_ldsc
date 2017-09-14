###
# run rmarkdown reports for each phenotype
###
setwd("/Users/raymondwalters/Documents/Code/github/UKBB_ldsc/")

# load packages
rm(list=ls())
require("rmarkdown")

# load univar data
load("./results/ukbb_h2univar.RData")
isUnivar <- T

# univar browser
# render('./site/h2_browser.Rmd',output_file='h2_univar_browser.html',output_dir = './docs')

# load partitioned data
load("./results/ukbb_h2part.RData")
isUnivar <- F

# partitioned browser
# render('./site/h2_browser.Rmd',output_file='h2_browser.html',output_dir = './docs')

# front pages
render('./site/welcome.Rmd',output_file='index.html',output_dir = './docs')
render('./site/details.Rmd',output_file='details.html',output_dir = './docs')

# loop pheno pages
# for(phen in dat$phenotype[c(1:5,sample(1:nrow(dat),5))]){
#	render('./site/h2_part_template.Rmd',
#		   output_file=paste0('h2_summary_',phen,'.html'),
#		   output_dir='./docs')
# }



