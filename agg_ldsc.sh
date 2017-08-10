#! /bin/bash
###
# download and organize ldsc h2 results
###

# configure environ to get gsutil working 
# (custom dotkit on Broad)
source /broad/software/scripts/useuse
reuse cloud

# stricter error catching while download/parse
set -e

#####
# UKBB 1859
#####

# download all
# gsutil cp gs://ukbb-gwas-results/ldsc_results/ukbb1859_h2_results.batch_* ./

# combine
# (remove previous version if necessary)
# if [ -f ukbb1859_h2_results.txt.gz ]; then
#	rm ukbb1859_h2_results.txt.gz
# fi
# for ii in `seq 1 10`; do
#	zcat ukbb1859_h2_results.batch_${ii}.txt.gz | awk -v i=$ii 'NR>1 || i==1' | gzip -c >> ukbb1859_h2_results.txt.gz
# 	rm ukbb1859_h2_results.batch_${ii}.txt.gz
# done

# upload to cloud
# gsutil cp ukbb1859_h2_results.txt.gz gs://ukbb-gwas-results/ldsc_results/ukbb1859_h2_results.txt.gz


#####
# UKBB 1189 ICD codes
#####

# download all
gsutil cp gs://ukbb-gwas-results/ldsc_results/ukbb1189_icd10_h2_results.batch_* ./

# combine
# (remove previous version if necessary)
if [ -f ukbb1189_icd10_h2_results.txt.gz ]; then
      rm ukbb1189_icd10_h2_results.txt.gz
fi
for ii in `seq 1 4`; do
	zcat ukbb1189_icd10_h2_results.batch_${ii}.txt.gz | awk -v i=$ii 'NR>1 || i==1' | gzip -c >> ukbb1189_icd10_h2_results.txt.gz
	rm ukbb1189_icd10_h2_results.batch_${ii}.txt.gz
done

# upload to cloud
gsutil cp ukbb1189_icd10_h2_results.txt.gz gs://ukbb-gwas-results/ldsc_results/ukbb1189_icd10_h2_results.txt.gz

##################

# Combine all apps
zcat ukbb1859_h2_results.txt.gz | gzip -c > ukbb_all_h2_results.txt.gz
zcat ukbb1189_icd10_h2_results.txt.gz | tail -n +2 | gzip -c >> ukbb_all_h2_results.txt.gz

# to the cloud
gsutil cp ukbb_all_h2_results.txt.gz gs://ukbb-gwas-results/ldsc_results/ukbb_all_h2_results.txt.gz

# make Rdata
set +e
use R-3.4
./make_rdata.R ukbb_all_h2_results.txt.gz dat.RData > make_rdata.log

# eof

