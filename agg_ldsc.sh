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

# download all
gsutil cp gs://ukbb-gwas-results/ldsc_results/ukbb1859_h2_results.batch_* ./

# combine ukbb1859
# (remove previous version if necessary)
if [ -f ukbb1859_h2_results.txt.gz ]; then
	rm ukbb1859_h2_results.txt.gz
fi
for ii in `seq 1 10`; do
	zcat ukbb1859_h2_results.batch_${ii}.txt.gz | awk -v i=$ii 'NR>1 || i==1' | gzip -c >> ukbb1859_h2_results.txt.gz
	rm ukbb1859_h2_results.batch_${ii}.txt.gz
done

# upload to cloud
gsutil cp ukbb1859_h2_results.txt.gz gs://ukbb-gwas-results/ldsc_results/ukbb1859_h2_results.txt.gz

# make Rdata
set +e
use R-3.4
./make_rdata.R ukbb1859_h2_results.txt.gz dat.RData > make_rdata.log

# eof

