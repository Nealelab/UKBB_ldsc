# /usr/bin/env Rscript

# load packages
rm(list=ls())
# devtools::install_github("rstudio/rmarkdown",ref="bcdb1a28b231576a98dea398986a868da804e609")
require("rmarkdown")

# setup
setwd("/Users/raymondwalters/Documents/Code/github/UKBB_ldsc/site_build/")
process <- FALSE
round1 <- FALSE
viz <- TRUE
pages <- FALSE
testing <- TRUE
h2 <- read.delim("../results/round2_final/ukb31063_h2_topline.02Oct2019.tsv.gz")

# build yaml and get files
system("cp ../site_source/yml_source/site_header.yml ./_site.yml")
system("cat ../site_source/yml_source/site_core.yml >> _site.yml")
system("cp ../site_source/rmd_source/index.Rmd ./")
system("cp ../site_source/rmd_source/downloads.Rmd ./")
system("cp ../site_source/rmd_source/methods.Rmd ./")
system("cp ../site_source/rmd_source/h2_credits.Rmd ./")
system("cp ../site_source/rmd_source/h2_browser.Rmd ./")
system("cp ../site_source/analytics_header.html ./")
system("cp ../site_source/sandflat/sandflat.min.css ./")
system("cp ../site_source/dt_style.css ./")
system("cp ../site_source/rmd_source/_code_highlight_fix.Rmd ./")
system("cp ../site_source/rmd_source/_toc_fix.Rmd ./")
system("cp -r ../site_source/icon ./")
system("cp -r ../reference/ukb_ord_codings_warn.txt ./")

if(round1){
  
  system("cat ../site_source/yml_source/site_round1_core.yml >> _site.yml")
  system("cp ../site_source/rmd_source/round1_details.Rmd ./")
  system("cp ../site_source/rmd_source/round1_h2_browser.Rmd ./")
  
  if(viz){
    system("cat ../site_source/yml_source/site_round1_viz.yml >> _site.yml")
    system("cp ../site_source/rmd_source/round1_plots_home.Rmd ./")
    system("cp ../site_source/rmd_source/round1_viz_annot.Rmd ./")
    system("cp ../site_source/rmd_source/round1_viz_h2.Rmd ./")
    system("cp ../site_source/rmd_source/round1_viz_qq.Rmd ./")
    system("cp ../site_source/rmd_source/round1_viz_sampsize.Rmd ./")
    system("cp ../site_source/rmd_source/round1_viz_univar.Rmd ./")
  }
}


if(viz){
  system("cat ../site_source/yml_source/site_viz.yml >> _site.yml")
  system("cp ../site_source/rmd_source/viz_annot.Rmd ./")
  system("cp ../site_source/rmd_source/viz_h2.Rmd ./")
  system("cp ../site_source/rmd_source/viz_sampsize.Rmd ./")
  system("cp ../site_source/rmd_source/viz_R1vR2.Rmd ./")
}

if(process){
  system("cp ../UKBB_ldsc_scripts/ukbb_h2_select_topline.Rmd ./")
  system("cp ../UKBB_ldsc_scripts/ukbb_h2_process_confidence.Rmd ./")
  system("cp ../UKBB_ldsc_scripts/ukbb_h2_process_significance.Rmd ./")
  system("cat ../site_source/yml_source/site_processing.yml >> _site.yml")
  system("cp -r ../external_images ./images")
}else{
  system("cp ../site_source/process_backup/*html ./")
}

if(pages){
  if(testing){
    set.seed(5)
    h2 <- h2[c(sample(1:nrow(h2),10),137,2942,732),] # manual exclusions are extreme high/low pvals
  }
  cat(paste0("\n  h2_summary_",h2$phenotype,".html:\n    src: \"h2_part_template.Rmd\"\n    params:\n      pheno: \"",h2$phenotype,"\"\n      dat_topline: \"../results/round2_final/ukb31063_h2_topline.02Oct2019.tsv.gz\"\n      dat_part: \"../results/round2_final/ukb31063_h2_z4.02Oct2019.tsv.gz\""), file="_site.yml", append=TRUE)
  system("cp ../site_source/rmd_source/h2_part_template.Rmd ./")
}
  
# create site
# system("ln -s ./_site.yml ./rmd_source/")
render_site()

# cleanup temp files
system("rm *Rmd")
system("rm analytics_header.html")
system("rm -r icon")
system("rm sandflat.min.css")
system("rm dt_style.css")
system("rm ukb_ord_codings_warn.txt")
if(process){
  system("rm -r images")
}else{
  system("rm confidence.html")
  system("rm significance.html")
  system("rm select_topline.html")
}

# eof
