#! /bin/bash
###
# download and organize ldsc h2 results
###

# configure environ to get gsutil working 
# (custom dotkit on Broad)
# source /broad/software/scripts/useuse
# reuse cloud

# stricter error catching while download/parse
set -e

# #####
# # UKBB 18597
# #####
#
# # download all
# gsutil cp gs://ukbb-gwas-results/ldsc_results/ukbb1859_h2_results.batch_* ./
#
# # combine
# # (remove previous version if necessary)
# if [ -f ukbb1859_h2_results.txt.gz ]; then
# 	rm ukbb1859_h2_results.txt.gz
# fi
# for ii in `seq 1 10`; do
# 	zcat ukbb1859_h2_results.batch_${ii}.txt.gz | awk -v i=$ii 'NR>1 || i==1' | gzip -c >> ukbb1859_h2_results.txt.gz
# 	rm ukbb1859_h2_results.batch_${ii}.txt.gz
# done
#
# # upload to cloud
# gsutil cp ukbb1859_h2_results.txt.gz gs://ukbb-gwas-results/ldsc_results/ukbb1859_h2_results.txt.gz
#
# # get cloud file and corresponding phenotype file
# # gsutil cp gs://ukbb-gwas-results/ldsc_results/ukbb1859_h2_results.txt.gz ukbb1859_h2_results.txt.gz
# gsutil cp gs://ukbb-gwas-results/ukb1859/ukb1859_phenosummary_final.tsv ukb1859_phenosummary_final.tsv
#
# # add labels
# awk -v FS='\t' 'NR==FNR && NR>1{
# 					samp[$1]=99;
# 					desc[$1]=$2;
# 					nn[$1]=$3;
# 					nca[$1]=$5;
# 					nco[$1]=$6;
# 					if($5=="NA"){
# 						prev[$1]="NA"
# 					}else{
# 						prev[$1]=$5/($5+$6)
# 					}
# 				}NR!=FNR{
# 					OFS="\t";
# 					if(FNR==1){
# 						$1 = "prevelence"
# 						print "phenotype","description","N","N_case","N_control",$0;
# 					}else{
# 						key = $1;
# 						if(samp[key]==99){
# 							$1 = prev[key];
# 							print key,desc[key],nn[key],nca[key],nco[key],$0;
# 						}else{
# 							$1 = "NO_PHENO";
# 							print key,"NO_PHENO","NO_PHENO","NO_PHENO","NO_PHENO",$0;
# 						}
# 					}
# 				}' ukb1859_phenosummary_final.tsv <(zless ukbb1859_h2_results.txt.gz) | gzip > ukbb18597_h2univar_results.labeled.tsv.gz
#
#
# #####
# # UKBB 11898 ICD codes
# #####
#
# # download all
# gsutil cp gs://ukbb-gwas-results/ldsc_results/ukbb1189_icd10_h2_results.batch_* ./
#
# # combine
# # (remove previous version if necessary)
# if [ -f ukbb1189_icd10_h2_results.txt.gz ]; then
#       rm ukbb1189_icd10_h2_results.txt.gz
# fi
# for ii in `seq 1 4`; do
# 	zless ukbb1189_icd10_h2_results.batch_${ii}.txt.gz | awk -v i=$ii 'NR>1 || i==1' | gzip -c >> ukbb1189_icd10_h2_results.txt.gz
# 	rm ukbb1189_icd10_h2_results.batch_${ii}.txt.gz
# done
# zless ukbb1189_icd10_h2_results.batch_fix.txt.gz | awk 'NR>1' | gzip -c >> ukbb1189_icd10_h2_results.txt.gz
# rm ukbb1189_icd10_h2_results.batch_fix.txt.gz
#
# # upload to cloud
# gsutil cp ukbb1189_icd10_h2_results.txt.gz gs://ukbb-gwas-results/ldsc_results/ukbb1189_icd10_h2_results.txt.gz
#
# # # get cloud file and corresponding phenotype file
# # gsutil cp gs://ukbb-gwas-results/ldsc_results/ukbb1189_icd10_h2_results.txt.gz ukbb1189_icd10_h2_results.txt.gz
# gsutil cp gs://ukbb-gwas-results/ukb1189/ukb1189_icd10_phenosummary_final.tsv ukb1189_icd10_phenosummary_final.tsv
#
# # add labels
# awk -v FS='\t' 'NR==FNR && NR>1{
# 					sub(/^41202_/, "", $1);
# 					samp[$1]=99;
# 					desc[$1]=$2;
# 					nn[$1]=$3;
# 					nca[$1]=$5;
# 					nco[$1]=$6;
# 					if($5=="NA"){
# 						prev[$1]="NA"
# 					}else{
# 						prev[$1]=$5/($5+$6)
# 					}
# 				}NR!=FNR{
# 					OFS="\t";
# 					if(FNR==1){
# 						$1 = "prevelence"
# 						print "phenotype","description","N","N_case","N_control",$0;
# 					}else{
# 						key = $1;
# 						if(samp[key]==99){
# 							$1 = prev[key];
# 							print key,desc[key],nn[key],nca[key],nco[key],$0;
# 						}else{
# 							$1 = "NO_PHENO";
# 							print key,"NO_PHENO","NO_PHENO","NO_PHENO","NO_PHENO",$0;
# 						}
# 					}
# 				}' ukb1189_icd10_phenosummary_final.tsv <(zless ukbb1189_icd10_h2_results.txt.gz) | gzip > ukbb11898_icd10_h2univar_results.labeled.tsv.gz
#
#
#
#
#####
# UKBB 18597 partitioned
#####

# download all
gsutil cp gs://ukbb-gwas-results/ldsc_results/batches/ukbb18597_h2part_results.batch_* ./

# combine
# (remove previous version if necessary)
if [ -f ukbb18597_h2part_results.txt.gz ]; then
	rm ukbb18597_h2part_results.txt.gz
fi
for ii in `seq 1 10`; do
	zless ukbb18597_h2part_results.batch_${ii}.txt.gz | awk -v i=$ii 'NR>1 || i==1' | gzip -c >> ukbb18597_h2part_results.txt.gz
	rm ukbb18597_h2part_results.batch_${ii}.txt.gz
done

# # upload to cloud
# gsutil cp ukbb18597_h2part_results.txt.gz gs://ukbb-gwas-results/ldsc_results/ukbb18597_h2part_results.txt.gz
#
# # get cloud file and corresponding phenotype file
# # gsutil cp gs://ukbb-gwas-results/ldsc_results/ukbb18597_h2part_results.txt.gz ukbb18597_h2part_results.txt.gz
gsutil cp gs://ukbb-gwas-results/ukb1859/ukb1859_phenosummary_final.tsv ukb1859_phenosummary_final.tsv
#
# add labels
awk -v FS='\t' 'NR==FNR && NR>1{
					samp[$1]=99;
					desc[$1]=$2;
					nn[$1]=$3;
					nca[$1]=$5;
					nco[$1]=$6;
					if($5=="NA"){
						prev[$1]="NA"
					}else{
						prev[$1]=$5/($5+$6)
					}
				}NR!=FNR{
					OFS="\t";
					if(FNR==1){
						$1 = "prevelence"
						print "phenotype","description","N","N_case","N_control",$0;
					}else{
						key = $1;
						if(samp[key]==99){
							$1 = prev[key];
							print key,desc[key],nn[key],nca[key],nco[key],$0;
						}else{
							$1 = "NO_PHENO";
							print key,"NO_PHENO","NO_PHENO","NO_PHENO","NO_PHENO",$0;
						}
					}
				}' ukb1859_phenosummary_final.tsv <(zless ukbb18597_h2part_results.txt.gz) | gzip > ukbb18597_h2part_results.labeled.tsv.gz




#####
# UKBB 11898 ICD codes partitioned
#####

# download all
gsutil cp gs://ukbb-gwas-results/ldsc_results/batches/ukbb11898_icd10_h2part_results.batch_* ./
#
# combine
# (remove previous version if necessary)
if [ -f ukbb11898_icd10_h2part_results.txt.gz ]; then
      rm ukbb11898_icd10_h2part_results.txt.gz
fi
for ii in `seq 1 8`; do
	zless ukbb11898_icd10_h2part_results.batch_${ii}.txt.gz | awk -v i=$ii 'NR>1 || i==1' | gzip -c >> ukbb11898_icd10_h2part_results.txt.gz
	rm ukbb11898_icd10_h2part_results.batch_${ii}.txt.gz
done
zless ukbb11898_icd10_h2part_results.batch_fix.txt.gz | awk 'NR>1' | gzip -c >> ukbb11898_icd10_h2part_results.txt.gz
rm ukbb11898_icd10_h2part_results.batch_fix.txt.gz
#
# # upload to cloud
# gsutil cp ukbb11898_icd10_h2part_results.txt.gz gs://ukbb-gwas-results/ldsc_results/ukbb11898_icd10_h2part_results.txt.gz
#
# # get cloud file and corresponding phenotype file
# # gsutil cp gs://ukbb-gwas-results/ldsc_results/ukbb1189_icd10_h2_results.txt.gz ukbb1189_icd10_h2_results.txt.gz
gsutil cp gs://ukbb-gwas-results/ukb1189/ukb1189_icd10_phenosummary_final.tsv ukb1189_icd10_phenosummary_final.tsv
#
# add labels
awk -v FS='\t' 'NR==FNR && NR>1{
					sub(/^41202_/, "", $1);
					samp[$1]=99;
					desc[$1]=$2;
					nn[$1]=$3;
					nca[$1]=$5;
					nco[$1]=$6;
					if($5=="NA"){
						prev[$1]="NA"
					}else{
						prev[$1]=$5/($5+$6)
					}
				}NR!=FNR{
					OFS="\t";
					if(FNR==1){
						$1 = "prevelence"
						print "phenotype","description","N","N_case","N_control",$0;
					}else{
						key = $1;
						if(samp[key]==99){
							$1 = prev[key];
							print key,desc[key],nn[key],nca[key],nco[key],$0;
						}else{
							$1 = "NO_PHENO";
							print key,"NO_PHENO","NO_PHENO","NO_PHENO","NO_PHENO",$0;
						}
					}
				}' ukb1189_icd10_phenosummary_final.tsv <(zless ukbb11898_icd10_h2part_results.txt.gz) | gzip > ukbb11898_icd10_h2part_results.labeled.tsv.gz



#####
# UKBB 11898 ICD codes partitioned
#####

# download all
gsutil cp gs://ukbb-gwas-results/ldsc_results/batches/ukbb11898_h2part_results.batch_* ./

# combine
# (remove previous version if necessary)
if [ -f ukbb11898_h2part_results.txt.gz ]; then
      rm ukbb11898_h2part_results.txt.gz
fi
for ii in `seq 1 3`; do
	zless ukbb11898_h2part_results.batch_${ii}.txt.gz | awk -v i=$ii 'NR>1 || i==1' | gzip -c >> ukbb11898_h2part_results.txt.gz
	rm ukbb11898_h2part_results.batch_${ii}.txt.gz
done

# upload to cloud
gsutil cp ukbb11898_h2part_results.txt.gz gs://ukbb-gwas-results/ldsc_results/batches/ukbb11898_h2part_results.txt.gz
#
# get cloud file and corresponding phenotype file
# gsutil cp gs://ukbb-gwas-results/ldsc_results/ukbb1189_icd10_h2_results.txt.gz ukbb1189_icd10_h2_results.txt.gz
gsutil cp gs://ukbb-gwas-results/ukb1189/ukb1189_phenosummary_final.tsv ukb1189_phenosummary_final.tsv

# add labels
awk -v FS='\t' 'NR==FNR && NR>1{
					sub(/^41202_/, "", $1);
					samp[$1]=99;
					desc[$1]=$2;
					nn[$1]=$3;
					nca[$1]=$5;
					nco[$1]=$6;
					if($5=="NA"){
						prev[$1]="NA"
					}else{
						prev[$1]=$5/($5+$6)
					}
				}NR!=FNR{
					OFS="\t";
					if(FNR==1){
						$1 = "prevelence"
						print "phenotype","description","N","N_case","N_control",$0;
					}else{
						key = $1;
						if(samp[key]==99){
							$1 = prev[key];
							print key,desc[key],nn[key],nca[key],nco[key],$0;
						}else{
							$1 = "NO_PHENO";
							print key,"NO_PHENO","NO_PHENO","NO_PHENO","NO_PHENO",$0;
						}
					}
				}' ukb1189_phenosummary_final.tsv <(zless ukbb11898_h2part_results.txt.gz) | gzip > ukbb11898_h2part_results.labeled.tsv.gz


##################


# Combine all apps univar h2
# zless ukbb18597_h2univar_results.labeled.tsv.gz | gzip -c > ukbb_all_h2univar_results.txt.gz
# zless ukbb11898_icd10_h2univar_results.labeled.tsv.gz | tail -n +2 | gzip -c >> ukbb_all_h2_results.txt.gz

# combine all apps partitioned h2
zless ukbb18597_h2part_results.labeled.tsv.gz | gzip -c > ukbb_all_h2part_results.txt.gz.tmp
zless ukbb11898_icd10_h2part_results.labeled.tsv.gz | tail -n +2 | gzip -c >> ukbb_all_h2part_results.txt.gz.tmp
zless ukbb11898_h2part_results.labeled.tsv.gz | tail -n +2 | gzip -c >> ukbb_all_h2part_results.txt.gz.tmp

# removing duplicated phenos
zless ukbb_all_h2part_results.txt.gz.tmp | awk -v FS='\t' -v OFS='\t' '{if(samp[$1]!=99){samp[$1]=99; print $0}}' | gzip > ukbb_all_h2part_results.txt.gz


# to the cloud
# gsutil cp ukbb_all_h2univar_results.txt.gz gs://ukbb-gwas-results/ldsc_results/ukbb_all_h2univar_results.txt.gz
# gsutil cp ukbb_all_h2part_results.txt.gz gs://ukbb-gwas-results/ldsc_results/ukbb_all_h2part_results.txt.gz


# make Rdata
set +e
use R-3.4
# ./make_rdata.R ukbb_all_h2univar_results.txt.gz ukbb_h2univar.RData > make_rdata_h2univar.log
./make_rdata.R ukbb_all_h2part_results.txt.gz ukbb_h2part.RData > make_rdata_h2part.log

# cleanup
rm *phenosummary_final.tsv
rm *.labeled.tsv.gz
rm ukbb11898_icd10_h2part_results.txt.gz
# rm ukbb1189_icd10_h2_results.txt.gz
rm ukbb18597_h2part_results.txt.gz
# rm ukbb1859_h2_results.txt.gz
rm ukbb11898_h2part_results.txt.gz

# eof

