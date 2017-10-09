#! /usr/bin/env python

###########################################

#####
# Dropping problematic SNPs from ldsc sumstats files
#
# Note: Only converting for phens with h2_z > 2
# 
# Need to address:
# - SNPs with alleles not matching previous HM3 --merge-alleles file
# - SNPs with strand ambiguous alleles
#
# Curated files for problematic SNPs in UKBB:
mismatch_file = 'gs://ukbb-ldsc-sumstats/miss.list'
ambiguous_file = 'gs://ukbb-ldsc-sumstats/ambiguous_removals.list'
#
# Directory with all ldsc sumstat files
in_dir = 'gs://ukbb-ldsc-sumstats/ldsc'
#
# Directory for output
out_dir = 'gs://ukbb-ldsc-sumstats/ldsc_for_rg'
#
# LDSR h2 results file
# Used to track index of phenotypes
# and get phens with h2_z > 2
h2_file = 'gs://ukbb-gwas-results/ldsc_results/ukbb_all_h2part_results.txt.gz'
#
# Number of parallel processes to run
num_proc=8
#
#####


####
# Package setup
####

import pandas as pd
import csv
import gzip
import datetime
import subprocess
from multiprocessing import Process, Pool
from hail import *

# start hail
hc = HailContext()

# Reinstall google cloud engine (broken by switching to anaconda)
subprocess.call(['/home/anaconda2/bin/pip','install','google-compute-engine'])


####
# Get list of SNPs to drop
####

with hadoop_read(mismatch_file) as f1:
    list1 = pd.read_table(f1)
    drop_snps = list1['SNP'].tolist()

with hadoop_read(ambiguous_file) as f2:
    list2 = pd.read_table(f2,header=None, names = ['SNP'])
    snps2 = list2['SNP'].tolist()

drop_snps.extend(snps2)


####
# Get list of phenotypes
####

with hadoop_read(h2_file) as f3:
    h2 = pd.read_table(f3)

phens = h2.loc[h2['h2_z'] > 2,'phenotype'].tolist()


####
# Function for rewriting file
####

def fix_sumstat(ph, in_dir=in_dir, out_dir=out_dir):
    
    # input filenames
    gsfile_in = in_dir +'/' + str(ph) + '.ldsc.tsv.gz'
    local_in = 'file:///mnt/data/tmp/' + str(ph) + '.ldsc.tsv.gz'
    in_path = '/mnt/data/tmp/' + str(ph) + '.ldsc.tsv.gz'
    
    # output filenames
    gsfile_out = out_dir +'/' + str(ph) + '.ukbb.sumstats.gz'
    local_out = 'file:///mnt/data/tmp/' + str(ph) + '.ukbb.sumstats.gz'
    out_path = '/mnt/data/tmp/' + str(ph) + '.ukbb.sumstats.gz'
    
    # download input
    subprocess.call(['gsutil','cp', gsfile_in, in_path])
    
    # stream through file
    with gzip.open(in_path, 'rb') as fin, gzip.open(out_path, 'wb') as fout:
        fread = csv.reader(fin, delimiter='\t')
        fwrite = csv.writer(fout, delimiter='\t')
        for row in fread:
            if row[0].strip() not in drop_snps:
                fwrite.writerow(row)
        fin.close()
        fout.close()
    
    # upload to cloud
    subprocess.call(['gsutil','cp', out_path, gsfile_out])
    
    # clean up
    subprocess.call(['rm', in_path, out_path])    
    
    # log and exit
    print('Finished '+str(ph) + ' at {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))
    return 0


####
# Parallelize conversion
####

print('Starting '+str(len(phens)) + ' conversions at {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))

pool = Pool(num_proc)
results = pool.imap_unordered(fix_sumstat, phens)
pool.close()
pool.join()

print('DONE!!! Finished '+str(len(phens)) + ' conversions at {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))


