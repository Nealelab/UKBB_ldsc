#! /usr/bin/env python

###########################################


########
#
# Genetic correlation analysis
# 1 phenotype vs. all significant UKBB
#
####
#
# Setup:
#
# local working dir
wd = '/home/mtag'
#
# h2 results file (to get list of UKBB sig. h2)
h2_file = 'gs://ukbb-gwas-results/ldsc_results/ukbb_all_h2part_results.txt.gz'
zthresh = 4
# 
# directory of UKBB sumstat files
gs_sumstat_dir = 'gs://ukbb-ldsc-sumstats/ldsc_for_rg'
#
# other sumstat
sumstat = 'gs://ukbb-rwalters-rg/sumstats/***pheno***.sumstats.gz'
target_name = 'pheno_name'
out_bucket = 'gs://ukbb-rwalters-rg/results'
#
# parallelization
chunk_size = 10
num_proc = 6
#
#
#
########


###
# load packages
print "Loading packages..."
###

from hail import *
import numpy as np
import pandas as pd
import os
import sys
import subprocess
import itertools
import datetime
from io import StringIO
from argparse import Namespace
from scipy import stats
from multiprocessing import Process, Pool
from functools import partial


###
# Install packages (MTAG, fix google-compute-engine)
print "Installing requirements..."
print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))
###

# Reinstall google cloud engine (broken by switching to anaconda)
subprocess.call(['pip','install','google-compute-engine'])

# Install joblib (for MTAG)
subprocess.call(['pip','install','joblib'])

# for ldsc
subprocess.call(['pip','install','bitarray'])

# Download MTAG for ref panel
if not os.path.isdir('/home/mtag/mtag-master'):
    subprocess.call(['wget', '--quiet', '-P', '/home/mtag/', 'https://github.com/omeed-maghzian/mtag/archive/master.zip'])
    subprocess.call(['unzip', '-q', '/home/mtag/master.zip', '-d', '/home/mtag'])
    subprocess.call(['rm','/home/mtag/master.zip'])

# Download ldsc
if not os.path.isdir('/home/mtag/ldsc-master'):
    subprocess.call(['wget', '--quiet', '-P', '/home/mtag/', 'https://github.com/bulik/ldsc/archive/master.zip'])
    subprocess.call(['unzip', '-q', '/home/mtag/master.zip', '-d', '/home/mtag'])    
    subprocess.call(['rm','/home/mtag/master.zip'])
    
# load MTAG and ldsc
sys.path.insert(0, '/home/mtag/ldsc-master/')
sys.path.insert(0, '/home/mtag/mtag-master/')
import ldscore.ldscore as ld
import ldscore.sumstats as sumstats
from mtag import Logger_to_Logging

# reference path
ld_ref_panel = '/home/mtag/mtag-master/ld_ref_panel/eur_w_ld_chr/'


#####
# Prep sumstat files
#####

print "Preparing files..."
print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))

### Get list of UKBB sumstats

# get from cloud
loc_h2_file = 'h2.txt.gz'
loc_h2_path = wd + '/' + loc_h2_file
subprocess.call(['gsutil','cp',h2_file,loc_h2_path])

# read
h2 = pd.read_table(loc_h2_path)

# filter to h2 z-score threshold
phens = h2.loc[h2['h2_z'] > zthresh,'phenotype'].tolist()

# split to chunks
ph_chunks = [phens[i:(i+chunk_size)] for i in xrange(0, len(phens), chunk_size)]

### Load other sumstat file

loc_target_ss = wd + '/' + 'target.sumstats.gz'
subprocess.call(['gsutil','cp',sumstat,loc_target_ss])


#####
# core function to run rg
#   on list of ukbb phenotypes
print "Preparing rg function..."
print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))
#####

def ldsc_rg_target(ph_list, **kwargs):
    """
    Assumes keyword args for:
    wd
    gs_sumstat_dir
    ld_ref_panel
    target_name
    """
    
    # log
    print "Starting phenotypes: "
    print ph_list
    
    # download sumstats for phens
    for ph in ph_list:
        gs_ss_path = gs_sumstat_dir + '/' + str(ph) + '.ukbb.sumstats.gz'
        loc_ss_path = wd + '/' + str(ph) + '.ukbb.sumstats.gz'
        subprocess.call(['gsutil','cp',gs_ss_path,loc_ss_path])    

    # list of files
    ukb_loc_list = ','.join([wd+'/'+str(x)+".ukbb.sumstats.gz" for x in ph_list])
    rg_file_list = ','.join([loc_target_ss,ukb_loc_list])

    # list of names
    ukb_name_list = [str(x)+".ukbb" for x in ph_list]
    rg_name_list = [target_name] + ukb_name_list

    # dummy output name
    rg_out = wd + '/' + 'rg.summary'

    # args for ldsc
    args_ldsc_rg =  Namespace(out=rg_out, 
                              bfile=None,
                              l2=None,
                              extract=None,
                              keep=None,
                              ld_wind_snps=None,
                              ld_wind_kb=None,
                              ld_wind_cm=None,
                              print_snps=None,
                              annot=None,
                              thin_annot=False,
                              cts_bin=None,
                              cts_break=None,
                              cts_names=None,
                              per_allele=False,
                              pq_exp=None,
                              no_print_annot=False,
                              maf=None,
                              h2=None,
                              rg=rg_file_list,
                              ref_ld=None,
                              ref_ld_chr=ld_ref_panel,
                              w_ld=None,
                              w_ld_chr=ld_ref_panel,
                              overlap_annot=False,
                              no_intercept=False,
                              intercept_h2=None, 
                              intercept_gencov=None,
                              M=None,
                              two_step=None,
                              chisq_max=None,
                              print_cov=False,
                              print_delete_vals=False,
                              chunk_size=50,
                              pickle=False,
                              invert_anyway=False,
                              yes_really=False,
                              n_blocks=200,
                              not_M_5_50=False,
                              return_silly_things=False,
                              no_check_alleles=False,
                              print_coefficients=False,
                              samp_prev=None,
                              pop_prev=None,
                              frqfile=None,
                              h2_cts=None,
                              frqfile_chr=None,
                              print_all_cts=False)

    # run rg
    rg_out = sumstats.estimate_rg(args_ldsc_rg, Logger_to_Logging())

    # format output
    rg_tab_txt = sumstats._get_rg_table(rg_name_list, rg_out, args_ldsc_rg)
    rg_df = pd.read_csv(StringIO(rg_tab_txt), delim_whitespace=True)
    
    return rg_df

ldsc_rg_map = partial(ldsc_rg_target,
                      wd=wd,
                      gs_sumstat_dir=gs_sumstat_dir,
                      ld_ref_panel=ld_ref_panel,
                      target_name=target_name)


# dispatch rg
print "Starting ldsc..."
print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))
pool = Pool(num_proc)
results = pool.imap_unordered(ldsc_rg_map, ph_chunks)
pool.close()
pool.join()

# collect results
print "Processing results..."
print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))
dat = pd.concat(results)


####
# write results to file
print "Saving..."
print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))
####
rg_outname = 'rg.'+target_name+'.txt.gz'
rg_local = wd+'/'+rg_outname
rg_cloud = out_bucket+'/'+rg_outname

# write local
dat.to_csv(rg_local, sep='\t', compression='gzip', index=False, na_rep="NA")

# move to cloud
subprocess.call(['gsutil','cp','file://'+rg_local,rg_cloud])

print "################"
print "Finished!!"
print "Results output to:"
print str(rg_cloud)
print "################"

# eof
