#! /usr/bin/env python

###########################################


# Settings:
wd = '/home/mtag/' # root working directory
ld_ref_panel = '/home/mtag/mtag-master/ld_ref_panel/baselineLD_v1.1/baselineLD.' # local path
ld_w_panel = '/home/mtag/mtag-master/ld_ref_panel/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.' # local path
ld_frq_panel = '/home/mtag/mtag-master/ld_ref_panel/1000G_Phase3_frq/1000G.EUR.QC.'
phen_summary = 'gs://ukbb-gwas-results/ukb1189/ukb1189_phenosummary_final.tsv' # in cloud
num_phens = 366
ss_bucket = 'gs://ukbb-gwas-results/ldsc' # 'gs://ukbb_association/ldsc/sumstats' # bucket with sumstats.gz files
out_bucket = 'gs://ukbb-gwas-results/ldsc_results/batches' # ouput google bucket location
num_proc = 6 # number of processes to run


# setup parallelization
import argparse
parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--parsplit', type=int, required=True, help="number of parallel batches to split phenotypes into")
parser.add_argument('--paridx', type=int, required=True, help="which of the phenotype batches to run")

args = parser.parse_args()

idx = xrange(args.paridx-1, num_phens, args.parsplit)
h2_fname = 'ukbb11898_h2part_results.batch_'+str(args.paridx)+'.txt.gz'

###########################################

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
subprocess.call(['/home/anaconda2/bin/pip','install','google-compute-engine'])

# Download MTAG and load it's version of ldsc (mtag allows call from inside python)
if not os.path.isdir('/home/mtag/mtag-master'):
    subprocess.call(['wget', '--quiet', '-P', '/home/mtag/', 'https://github.com/omeed-maghzian/mtag/archive/master.zip'])
    subprocess.call(['unzip', '-q', '/home/mtag/master.zip', '-d', '/home/mtag'])

# load MTAG
sys.path.insert(0, '/home/mtag/mtag-master/')
sys.path.insert(0, '/home/mtag/mtag-master/ldsc_mod/')
import ldsc_mod.ldscore as ldsc
from mtag import Logger_to_Logging

# download partitioned LD scores
subprocess.call(['wget', '--quiet', '-P', '/home/mtag/mtag-master/ld_ref_panel/', 'https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baselineLD_v1.1_ldscores.tgz'])
subprocess.call(['wget', '--quiet', '-P', '/home/mtag/mtag-master/ld_ref_panel/', 'https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_frq.tgz'])
subprocess.call(['wget', '--quiet', '-P', '/home/mtag/mtag-master/ld_ref_panel/', 'https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_weights_hm3_no_MHC.tgz'])
subprocess.call(['tar', '-zxvf', '/home/mtag/mtag-master/ld_ref_panel/1000G_Phase3_baselineLD_v1.1_ldscores.tgz','-C','/home/mtag/mtag-master/ld_ref_panel/'])
subprocess.call(['tar', '-zxvf', '/home/mtag/mtag-master/ld_ref_panel/1000G_Phase3_frq.tgz','-C','/home/mtag/mtag-master/ld_ref_panel/'])
subprocess.call(['tar', '-zxvf', '/home/mtag/mtag-master/ld_ref_panel/1000G_Phase3_weights_hm3_no_MHC.tgz','-C','/home/mtag/mtag-master/ld_ref_panel/'])

###
# Define ldsc calling/handling functions
print "Preparing ldsc functions..."
print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))
###


####
# Define handling of ldsc h2 output
####

def process_h2_part(h2_results, phname, outfile, ncas=None, ncon=None):
    
    # extract ldsc h2 results
    h2obs, h2obs_se = h2_results.tot, h2_results.tot_se
    inter, inter_se = h2_results.intercept, h2_results.intercept_se
    lam = h2_results.lambda_gc
    mchi = h2_results.mean_chisq
    ratio, ratio_se = h2_results.ratio, h2_results.ratio_se
    
    # get p-values
    h2z = h2obs/h2obs_se
    h2p = stats.norm.sf(h2z)
    intz = (inter-1.0)/inter_se
    intp = stats.norm.sf(intz)
    
    
    # get liability-scale values
    if pd.notnull(ncas) and pd.notnull(ncon):
        cc = ldsc.regressions.h2_obs_to_liab(1, P=(ncas/(ncas+ncon)), K=(ncas/(ncas+ncon)))
        h2liab = abs(cc) * h2obs
        h2liab_se = abs(cc) * h2obs_se
    else:
        h2liab = h2obs
        h2liab_se = h2obs_se

    # get global results
    dat1 = {'phenotype' : phname,
            'mean_chi2' : mchi,
            'lambdaGC' : lam,
            'intercept' : inter,
            'intercept_se' : inter_se,
            'intercept_z' : intz,
            'intercept_p' : intp,
            'ratio' : ratio,
            'ratio_se' : ratio_se,
            'h2_observed' : h2obs,
            'h2_observed_se' : h2obs_se,
            'h2_liability' : h2liab,
            'h2_liability_se' : h2liab_se,
            'h2_z' : h2z,
            'h2_p' : h2p
            }
    names = ['phenotype',
             'mean_chi2',
             'lambdaGC',
             'intercept',
             'intercept_se',
             'intercept_z',
             'intercept_p',
             'ratio',
             'ratio_se',
             'h2_observed',
             'h2_observed_se',
             'h2_liability',
             'h2_liability_se',
             'h2_z',
             'h2_p']
    
    
    # get category results and flatten
    cat_df = pd.read_table(outfile)

    vals = [phname]
    for row in cat_df.itertuples():
        cat = str(row[1]).replace("L2_0","")
        # record column names
        names.extend([cat+'::Prop_SNPs',
                      cat+'::Prop_h2',
                      cat+'::Prop_h2_se',
                      cat+'::Enrichment',
                      cat+'::Enrichment_se',
                      cat+'::Enrichment_p',
                      cat+'::Coefficient',
                      cat+'::Coefficient_se',
                      cat+'::Coefficient_z',
                      cat+'::Coefficient_p'])
        # record results
        dat1.update({
            cat+'::Prop_SNPs' : row[2],
            cat+'::Prop_h2' : row[3],
            cat+'::Prop_h2_se' : row[4],
            cat+'::Enrichment' : row[5],
            cat+'::Enrichment_se' : row[6],
            cat+'::Enrichment_p' : row[7],
            cat+'::Coefficient' : row[8],
            cat+'::Coefficient_se' : row[9],
            cat+'::Coefficient_z' : row[10],
            cat+'::Coefficient_p' : 2*stats.norm.sf(abs(row[10]))
        })

    df = pd.DataFrame(data=dat1, index=pd.Series([phname]))
                
    return df[names]




####
# Define core function to run ldsc
####

def ldsc_h2_part(args, **kwargs):
    """
    Runs LD score to estimate h2 for the named UKBB phenotype
    
    Args is a list with elements:
    - args[0] = phenotype name
    - args[1] = N_cases
    - args[2] = N_controls
    
    # keyword args is for global settings:
    - wd            (working directory)
    - ld_ref_panel  (local path, supplied to --ref-ld-chr)
    - ld_w_panel    (local path, supplied to --w-ld-chr)
    - ld_frq_panel  (local path, supplied to --frqfile-chr)
    - ss_bucket     (cloud bucket containing sumstats files)
    
    Using this structure for the sake of multiprocessing.pool.map()
    """
    
    # handle args
    phname = str(args[0])
    ncas = float(args[1])
    ncon = float(args[2])
    
    # define names
    
    ss_name = str(phname)+'.ldsc.tsv.gz'
    sspath_local = wd+'/'+ss_name
    sspath_cloud = ss_bucket+'/'+ss_name
    h2_out = 'h2.part.ukbb.'+str(phname)
    
    # download sumstats file
    subprocess.call(['gsutil','cp',sspath_cloud,sspath_local])
    
    # run ldsc
    args_h2 =  Namespace(out=h2_out, 
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
                         maf=0.05,
                         h2=sspath_local,
                         rg=None,
                         ref_ld=None,
                         ref_ld_chr=ld_ref_panel,
                         w_ld=None,
                         w_ld_chr=ld_w_panel,
                         overlap_annot=True,
                         no_intercept=False,
                         intercept_h2=None,
                         intercept_gencov=None,
                         M=None,
                         two_step=None,
                         chisq_max=9999,
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
                         print_coefficients=True,
                         samp_prev=None,
                         pop_prev=None,
                         frqfile=None,
                         h2_cts=None,
                         frqfile_chr=ld_frq_panel,
                         print_all_cts=False,
                         sumstats_frames=None,
                         rg_mat=False)
                     
    print "Launching ldsc for "+str(phname)
    h2_results = ldsc.sumstats.estimate_h2(args_h2, Logger_to_Logging())
    print "Completed ldsc for "+str(phname)
    
    # cleanup sumstats file
    subprocess.call(['rm',sspath_local])
    
    return process_h2_part(h2_results, phname, h2_out+'.results', float(ncas), float(ncon))




def strip_icd(text):
    if text.startswith("41202_"):
        return text[6:]
    else:
        return text


if __name__ == "__main__":

    print "Starting Hail..."
    print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))
    hc = HailContext()
    

    ####
    # Load phenotype info
    print "Reading phenotype summary..."
    print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))
    ####
    
    with hadoop_read(phen_summary) as f:
        phens = pd.read_table(f)

    # get subset of phenotypes to run
    ph_sub = phens.iloc[idx]
    ph_list = list(ph_sub.index.values)
    
    # strip ICD code prefix (not used in ldsc sumstats file names)
    ph_list = [strip_icd(x) for x in ph_list]

    print "Phenotypes to run:"
    print ph_list

    ####
    # Run phenotypes
    ####
    
    # zip arguments to map
    iter_args = itertools.izip(ph_list, ph_sub['N.cases'], ph_sub['N.controls'])

    # bake in globals
    ldsc_h2_map = partial(ldsc_h2_part, wd=wd, ld_ref_panel=ld_ref_panel, ld_w_panel=ld_w_panel, ld_frq_panel=ld_frq_panel, ss_bucket=ss_bucket)

    # dispatch
    print "Starting ldsc..."
    print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))
    pool = Pool(num_proc)
    results = pool.imap_unordered(ldsc_h2_map, iter_args)
    pool.close()
    pool.join()


    ####
    # Load output to dataframe
    print "Processing results..."
    print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))
    ####
    dat = pd.concat(results)
    # dat = pd.DataFrame(index=ph_list, columns=col_ord)
    # for res in results:
    #     dat.update(pd.DataFrame(data=res, index=pd.Series(res[0][0]), columns=col_ord))


    ####
    # write results to file
    print "Saving..."
    print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))
    ####
    h2_local = wd+'/'+h2_fname
    h2_cloud = out_bucket+'/'+h2_fname

    # write local
    dat.to_csv(h2_local, sep='\t', compression='gzip', index=False, na_rep="NA")

    # move to cloud
    hadoop_copy('file://'+h2_local, h2_cloud)
    
    
    print "################"
    print "Finished!!"
    print "Results output to:"
    print str(h2_cloud)
    print "################"
# eof
