# /usr/bin/env Rscript

# load packages
rm(list=ls())
# devtools::install_github("rstudio/rmarkdown",ref="bcdb1a28b231576a98dea398986a868da804e609")
require("rmarkdown")

# setup
setwd("/Users/raymondwalters/Documents/Code/github/UKBB_ldsc/site_build/")
process <- FALSE
viz <- FALSE
pages <- FALSE
testing <- TRUE
h2 <- read.delim("../results/round2_final/ukb31063_h2_topline.02Oct2019.tsv.gz")

# build yaml and get files
system("cp ../site_source/yml_source/site_header.yml ./_site.yml")
system("cat ../site_source/yml_source/site_core.yml >> _site.yml")
system("cp ../site_source/rmd_source/index.Rmd ./")
system("cp ../site_source/rmd_source/details.Rmd ./")
system("cp ../site_source/rmd_source/h2_browser.Rmd ./")
system("cp ../site_source/analytics_header.html ./")
system("cp -r ../site_source/icon ./")

if(viz){
  system("cat ../site_source/yml_source/site_viz.yml >> _site.yml")
  system("cp ../site_source/rmd_source/plots_home.Rmd ./")
  system("cp ../site_source/rmd_source/viz_annot.Rmd ./")
  system("cp ../site_source/rmd_source/viz_h2.Rmd ./")
  system("cp ../site_source/rmd_source/viz_qq.Rmd ./")
  system("cp ../site_source/rmd_source/viz_sampsize.Rmd ./")
  system("cp ../site_source/rmd_source/viz_univar.Rmd ./")
}

if(process){
  system("cp ../UKBB_ldsc_scripts/ukbb_h2_select_topline.Rmd ./")
  system("cp ../UKBB_ldsc_scripts/ukbb_h2_process_confidence.Rmd ./")
  system("cp ../UKBB_ldsc_scripts/ukbb_h2_process_significance.Rmd ./")
  system("cat ../site_source/yml_source/site_processing.yml >> _site.yml")
  system("cp -r ../external_images ./images")
}

if(pages){
  if(testing){
    set.seed(5)
    h2 <- h2[sample(1:nrow(h2),10),]

    # cat(paste0(" h2_summary_",d2,".html:\n    src: \"h2_part_template.Rmd\"\n    params:\n      pheno: \"",d2,"\"\n      datfile: \"../results/ukbb_h2part.RData\"\n"), file="tmp.yaml")
  }else{
    set.seed(5)
  }
}
  
# create site
# system("ln -s ./_site.yml ./rmd_source/")
render_site()

# cleanup temp files
system("rm *Rmd")
system("rm analytics_header.html")
system("rm -r icon")
if(process){
  system("rm -r images")
}

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
