###
# run rmarkdown reports for each phenotype
###
setwd("/Users/raymondwalters/Documents/Code/github/UKBB_ldsc/site/")

# load packages
rm(list=ls())
# devtools::install_github("rstudio/rmarkdown",ref="bcdb1a28b231576a98dea398986a868da804e609")
require("rmarkdown")

# build
render_site()

#########

# for yaml for per-pheno sites
# load("../results/ukbb_h2part.RData")

# dat$isBinary <- !is.na(dat$N_case)
# dat$Neff <- dat$N
# dat$Neff[dat$isBinary] <- round( (4/((1/dat$N_case)+(1/dat$N_control)))[dat$isBinary], 2)
# dat <- dat[dat$Neff > 200,]

# d2 <- dat$phenotype
# cat(paste0(" h2_summary_",d2,".html:\n    src: \"h2_part_template.Rmd\"\n    params:\n      pheno: \"",d2,"\"\n      datfile: \"../results/ukbb_h2part.RData\"\n"), file="tmp.yaml")


#####################

### old

# load univar data
# load("./results/ukbb_h2univar.RData")
# isUnivar <- T

# univar browser
# render('./site/h2_browser.Rmd',output_file='h2_univar_browser.html',output_dir = './docs')

# load partitioned data
# load("./results/ukbb_h2part.RData")
# isUnivar <- F

# partitioned browser
# render('./site/h2_browser.Rmd',output_file='h2_browser.html',output_dir = './docs')

# front pages
# render('./site/welcome.Rmd',output_file='index.html',output_dir = './docs')
# render('./site/details.Rmd',output_file='details.html',output_dir = './docs')

# loop pheno pages
# for(phen in dat$phenotype[c(1:5,sample(1:nrow(dat),5))]){
#	dat_sub <- dat[which(dat$phenotype==phen),]
#	render('./site/h2_part_template.Rmd',
#		   output_file=paste0('h2_summary_',phen,'.html'),
#		   output_dir='./docs')
#}
