#! /usr/bin/env python

print "Preparing packages, etc..."

# load packages
from hail import *
import subprocess
import os
hc = HailContext()


# Reinstall google cloud engine (broken by switching to anaconda)
subprocess.call(['/home/anaconda2/bin/pip','install','google-compute-engine'])


# utility to check if file exists on cloud
def gs_missing(filepath):
    
    # strip tailing '/'
    if str(filepath).endswith('/'):
        tmp = str(filepath)[:-1]
        filepath = tmp
    
    # check for file
    stat = subprocess.call(['gsutil','-q','stat',str(filepath)])
    
    # check for folder
    if stat:
        stat = subprocess.call(['gsutil','-q','stat',str(filepath)+'/'])
    
    return stat


###
# Load hapmap3 sites vcf to vds
#   hm3 file is from Broad resource bundle:
#   wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/hapmap_3.3_hg19_pop_stratified_af.vcf.gz
###

if gs_missing('gs://ukbb_association/ldsc/ld_ref_panel/hapmap_3.3_hg19_pop_stratified_af.vcf.gz'):

    print "Downloading HapMap3 vcf..."
    
    # downlaod from Broad
    subprocess.call(['wget',
                    '--quiet', 
                    '-P', 
                    '/home/hapmap/',
                    'ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/hapmap_3.3_hg19_pop_stratified_af.vcf.gz'])
    
    # upload to cloud
    subprocess.call(['gsutil', 
                     'cp', 
                     '/home/hapmap/hapmap_3.3_hg19_pop_stratified_af.vcf.gz',
                     'gs://ukbb_association/ldsc/ld_ref_panel/hapmap_3.3_hg19_pop_stratified_af.vcf.gz'])
    

if gs_missing('gs://ukbb_association/ldsc/ld_ref_panel/hm3.r3.hg19.vds'): 
    
    print "Creating vds of HapMap3 sites..."
    
    (hc
        .import_vcf('gs://ukbb_association/ldsc/ld_ref_panel/hapmap_3.3_hg19_pop_stratified_af.vcf.gz', force=True)
        .write('gs://ukbb_association/ldsc/ld_ref_panel/hm3.r3.hg19.vds', overwrite=True)
    )



###
# Create filtered keytable with autosomal, biallelic HM3 snps
# with MAF > .01, INFO > 0.9 in full UKBB data
# (based on mfi file)
###
if gs_missing('gs://ukbb_association/ldsc/ld_ref_panel/hm3.r3.hg19.auto_bi_af.kt') or gs_missing('gs://ukbb_association/ldsc/ld_ref_panel/hm3.r3.hg19.auto_bi_af.ukbb_full_qcpos.kt'):

    ###
    # Create filtered keytable with only autosomal, biallelic snps
    # - also parse to keep allele freqs by poplutation
    print "Creating keytable of autosomal, biallelic HM3 SNPs with population allele freqs..."
    ###
    
    # read full HM3
    hm3_vds = hc.read('gs://ukbb_association/ldsc/ld_ref_panel/hm3.r3.hg19.vds')
    
    # filter to autosomal, biallelic snps and extract allelic freqs
    hm3_keep = (hm3_vds.filter_variants_expr('''
                                                v.isBiallelic() && 
                                                v.isAutosomal() &&
                                                v.altAllele().isSNP()
                                            ''')
                        .annotate_variants_expr('''va.ASW_AF = 1.0 - va.info.ASW[0].split("=")[1].toFloat(), 
                                                   va.CEU_AF = 1.0 - va.info.CEU[0].split("=")[1].toFloat(), 
                                                   va.CHB_AF = 1.0 - va.info.CHB[0].split("=")[1].toFloat(),
                                                   va.CHD_AF = 1.0 - va.info.CHD[0].split("=")[1].toFloat(),
                                                   va.CHS_AF = 1.0 - va.info.CHS[0].split("=")[1].toFloat(),
                                                   va.CLM_AF = 1.0 - va.info.CLM[0].split("=")[1].toFloat(),
                                                   va.FIN_AF = 1.0 - va.info.FIN[0].split("=")[1].toFloat(),
                                                   va.GBR_AF = 1.0 - va.info.GBR[0].split("=")[1].toFloat(),
                                                   va.GIH_AF = 1.0 - va.info.GIH[0].split("=")[1].toFloat(),
                                                   va.IBS_AF = 1.0 - va.info.IBS[0].split("=")[1].toFloat(),
                                                   va.JPT_AF = 1.0 - va.info.JPT[0].split("=")[1].toFloat(),
                                                   va.LWK_AF = 1.0 - va.info.LWK[0].split("=")[1].toFloat(),
                                                   va.MKK_AF = 1.0 - va.info.MKK[0].split("=")[1].toFloat(),
                                                   va.MXL_AF = 1.0 - va.info.MXL[0].split("=")[1].toFloat(),
                                                   va.PUR_AF = 1.0 - va.info.PUR[0].split("=")[1].toFloat(),
                                                   va.TSI_AF = 1.0 - va.info.TSI[0].split("=")[1].toFloat(),
                                                   va.YRI_AF = 1.0 - va.info.YRI[0].split("=")[1].toFloat()
                                                '''))

    # save autosomal, biallelic HM3 snps (w/o QC) keytable
    if gs_missing('gs://ukbb_association/ldsc/ld_ref_panel/hm3.r3.hg19.auto_bi_af.kt'):
        hm3_keep.variants_table().write('gs://ukbb_association/ldsc/ld_ref_panel/hm3.r3.hg19.auto_bi_af.kt', overwrite=True)
    
    
    # get subset passing QC in UKBB full set
    if gs_missing('gs://ukbb_association/ldsc/ld_ref_panel/hm3.r3.hg19.auto_bi_af.ukbb_full_qcpos.kt'):
        
        print "Creating keytable of HM3 SNPs passing UKBB Full Cohort QC..."
        
        qcpos_kt = (
            hc.import_table(
                'gs://ukbb_association/ukb_mfi_v2.tsv.bgz', 
                impute=True, 
                no_header=True)
              .rename({
                  'f0': 'CHR',
                  'f1': 'RS',
                  'f2': 'BP',
                  'f3': 'REF',
                  'f4': 'ALT',
                  'f5': 'UKBB_MAF',
                  'f6': 'UKBB_INFO'})
              .filter('UKBB_MAF > 0.01 && UKBB_INFO > 0.9')
              .annotate('v = Variant("chr"+str(CHR), BP, REF, ALT)')
              .key_by('v')
        )
    
        # get hm3 biallelic
        hm3_vds = hc.read('gs://ukbb_association/ldsc/ld_ref_panel/hm3.r3.hg19.vds')

        hm3_pass = (hm3_keep.filter_variants_table(qcpos_kt, keep=True)
                            .annotate_variants_table(qcpos_kt,
                                                     expr = '''va.UKBB_MAF = table.UKBB_MAF, 
                                                     va.UKBB_INFO = table.UKBB_INFO''')
                    )

        hm3_pass.variants_table().write('gs://ukbb_association/ldsc/ld_ref_panel/hm3.r3.hg19.auto_bi_af.ukbb_full_qcpos.kt', overwrite=True)



###
# Create filtered keytable with autosomal, biallelic HM3 snps
# with MAF > .01, INFO > 0.9 and passing QC in UKBB GWAS analysis set (N=337,202)
#
# Also provides both UKBB and HM3 formatting of variant
###

if gs_missing('gs://ukbb_association/ldsc/ld_ref_panel/hm3.r3.hg19.auto_bi_af.ukbb_gwas_qcpos.kt'):
    
    print "Creating keytable of HM3 SNPs passing UKBB GWAS QC..."

    # get full hm3 list (autosomal, biallelic snps)
    # and drop extra info
    hm3_snps = (
                hc.read_table('gs://ukbb_association/ldsc/ld_ref_panel/hm3.r3.hg19.auto_bi_af.kt')
                  .annotate('''v_hm3 = v,
                               HM3_rsid = va.rsid
                            ''')
                  .select(['v_hm3','HM3_rsid'])
               )

    # get subset passing QC in UKBB analyses set (N=337,202)
    # with with MAF > .01, INFO > 0.9
    qc2_kt = (
                hc.read_table('gs://ukbb_association/gwas_variants.kt')
                  .filter('info > 0.9 && AF > 0.01 && AF < 0.99')
                  .annotate('''v_hm3 = Variant("chr"+v.contig.replace("^0",""), v.start, v.ref, v.alt),
                               v_ukb = v,
                               UKBB_rsid = rsid
                            ''')
                  .select(['v_hm3','v_ukb','UKBB_rsid'])
             )
    
    # merge
    hm3_save = hm3_snps.key_by('v_hm3').join(qc2_kt.key_by('v_hm3'))

    # save
    hm3_save.write('gs://ukbb_association/ldsc/ld_ref_panel/hm3.r3.hg19.auto_bi_af.ukbb_gwas_qcpos.kt', overwrite=True)


print '''
##############
Finished!!

HapMap3 vcf (from ftp.broadinstitute.org/bundle/hg19/):
gs://ukbb_association/ldsc/ld_ref_panel/hapmap_3.3_hg19_pop_stratified_af.vcf.gz

HapMap3 full vds: 
gs://ukbb_association/ldsc/ld_ref_panel/hm3.r3.hg19.vds

HapMap3 autosomal, biallelic SNPs with freqs keytable:
gs://ukbb_association/ldsc/ld_ref_panel/hm3.r3.hg19.auto_bi_af.kt

HapMap3 SNPs passing UKBB full cohort QC:
gs://ukbb_association/ldsc/ld_ref_panel/hm3.r3.hg19.auto_bi_af.ukbb_full_qcpos.kt

HapMap3 SNPs passing UKBB GWAS QC:
gs://ukbb_association/ldsc/ld_ref_panel/hm3.r3.hg19.auto_bi_af.ukbb_gwas_qcpos.kt
##############
'''

# eof