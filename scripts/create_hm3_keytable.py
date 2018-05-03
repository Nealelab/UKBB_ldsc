#! /usr/bin/env python

print("Preparing packages, etc...")

# Note: not reusing the gs://ukbb_association copies of 
# HM3 reference files due to convserion to Hail 0.2

####
# setup
####

# re-run everything?
force_all = False

### input

# where to save files
BUCKET = 'gs://ukb31063-mega-gwas/ldsc/ld_ref_panel/'

# source of HapMap3 sites file (path only)
HM3_FTP = 'ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/'

# name of original HapMap3 sites file
HM3_NAME = 'hapmap_3.3_b37_pop_stratified_af.vcf.gz'

# UKBB variant QC in the full sample (for INFO scores)
GWAS_mfi = 'gs://ukb31063-mega-gwas/qc/ukb31063.imputed_v3.mfi.ht'

# UKBB variant QC in the GWAS samples (for MAF)
GWAS_qc = 'gs://ukb31063-mega-gwas/qc/ukb31063.imputed_v3.variant_qc.autosomes.ht'

# SNPs passing QC to be included in the UKBB GWAS
GWAS_snps = 'gs://ukb31063-mega-gwas/qc/ukb31063.gwas_variants.autosomes.ht'


### output

# name to save HapMap3 sites as matrixtable
HM3_vds = 'hm3.r3.b37.mt'

# name to save autosomal, biallelic SNPs (also has parsed population allele freqs)
HM3_bi_af = 'hm3.r3.b37.auto_bi_af.ht'

# name to save HapMap3 SNPs passing ldsc QC in UKBB GWAS sample
# - will save both Hail table (.ht) and flat file (.tsv)
HM3_qcpos_stem = 'hm3.r3.b37.auto_bi_af.ukbb_gwas_qcpos'



# load packages
import hail as hail
import subprocess
hail.init()


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
# Load hapmap3 sites vcf to hail table
#   hm3 file is from Broad resource bundle:
#   wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/hapmap_3.3_b37_pop_stratified_af.vcf.gz
###

if gs_missing(BUCKET + HM3_NAME) or force_all:

    print("Downloading HapMap3 vcf...")
    
    # download from Broad
    subprocess.call(['wget',
                    '--quiet', 
                    '-P', 
                    '/home/hapmap/',
                    HM3_FTP + HM3_NAME])
    
    # upload to cloud
    subprocess.call(['gsutil', 
                     'cp', 
                     '/home/hapmap/' + HM3_NAME,
                     BUCKET + HM3_NAME])
    

if gs_missing(BUCKET + HM3_vds) or force_all: 
    
    print("Creating Hail table of HapMap3 sites...")
        
    (hail
        .import_vcf(BUCKET + HM3_NAME, force_bgz=True, min_partitions=500)
        .write(BUCKET + HM3_vds, overwrite=True)
    )



###
# Create filtered keytable with autosomal, biallelic HM3 snps
# with MAF > .01, INFO > 0.9 and passing QC in UKBB GWAS analysis set
###
if gs_missing(BUCKET + HM3_bi_af) or force_all:

    ###
    # Create filtered keytable with only autosomal, biallelic snps
    # - also parse to keep allele freqs by population
    print("Creating keytable of autosomal, biallelic HM3 SNPs with population allele freqs...")
    ###
    
    # read full HM3 sites
    hm3_site = hail.read_matrix_table(BUCKET + HM3_vds).rows()
    
    # filter to autosomal, biallelic snps and extract allelic freqs
    hm3_keep = hm3_site.filter(hm3_site.locus.in_autosome() & 
                               (hail.len(hm3_site.alleles) == 2) & 
                               hail.is_snp(hm3_site.alleles[0],hm3_site.alleles[1]))
    
    hm3_keep = hm3_keep.annotate(ASW_AF = 1.0 - hail.float(hm3_keep.info['ASW'][0].split("=")[1]), 
                                 CEU_AF = 1.0 - hail.float(hm3_keep.info['CEU'][0].split("=")[1]),
                                 CHB_AF = 1.0 - hail.float(hm3_keep.info['CHB'][0].split("=")[1]),
                                 CHD_AF = 1.0 - hail.float(hm3_keep.info['CHD'][0].split("=")[1]),
                                 CHS_AF = 1.0 - hail.float(hm3_keep.info['CHS'][0].split("=")[1]),
                                 CLM_AF = 1.0 - hail.float(hm3_keep.info['CLM'][0].split("=")[1]),
                                 FIN_AF = 1.0 - hail.float(hm3_keep.info['FIN'][0].split("=")[1]),
                                 GBR_AF = 1.0 - hail.float(hm3_keep.info['GBR'][0].split("=")[1]),
                                 GIH_AF = 1.0 - hail.float(hm3_keep.info['GIH'][0].split("=")[1]),
                                 IBS_AF = 1.0 - hail.float(hm3_keep.info['IBS'][0].split("=")[1]),
                                 JPT_AF = 1.0 - hail.float(hm3_keep.info['JPT'][0].split("=")[1]),
                                 LWK_AF = 1.0 - hail.float(hm3_keep.info['LWK'][0].split("=")[1]),
                                 MKK_AF = 1.0 - hail.float(hm3_keep.info['MKK'][0].split("=")[1]),
                                 MXL_AF = 1.0 - hail.float(hm3_keep.info['MXL'][0].split("=")[1]),
                                 PUR_AF = 1.0 - hail.float(hm3_keep.info['PUR'][0].split("=")[1]),
                                 TSI_AF = 1.0 - hail.float(hm3_keep.info['TSI'][0].split("=")[1]),
                                 YRI_AF = 1.0 - hail.float(hm3_keep.info['YRI'][0].split("=")[1]))


    # save keytable
    hm3_keep.write(BUCKET + HM3_bi_af, overwrite=True)


###
# Create filtered keytable with autosomal, biallelic HM3 snps
# with MAF > .01, INFO > 0.9 and passing QC in UKBB GWAS analysis set
#
# Also provides both UKBB and HM3 formatting of variant
###

if gs_missing(BUCKET + HM3_qcpos_stem + '.ht') or force_all:
    
    print("Creating Hail table of HM3 SNPs passing UKBB GWAS QC...")


    # get list of SNPs to be used in GWAS
    ukb_snps = hail.read_table(GWAS_snps).repartition(500, shuffle=False)
    
    # get full hm3 list (autosomal, biallelic snps)
    # and drop extra info
    hm3_snps = hail.read_table(BUCKET + HM3_bi_af).repartition(500, shuffle=False)
    hm3_snps = hm3_snps.select(hm3_snps.locus, hm3_snps.alleles, hm3_snps.rsid)    
    
    # merge UKB passing SNPs to HM3 sites
    # - matching on locus + sorted alleles + rsid
    #   (are 1-2k SNPs with matching locus and different rsid)
    hm3_snps = hm3_snps.annotate(sort_al=hail.sorted(hm3_snps.alleles)).key_by('locus','sort_al','rsid')
    ukb_snps = ukb_snps.annotate(sort_al=hail.sorted(ukb_snps.alleles)).key_by('locus','sort_al','rsid')
    ukb_hm3 = ukb_snps.join(hm3_snps).key_by('locus','alleles')
    ukb_hm3 = ukb_hm3.drop(ukb_hm3.alleles_1)
    
    # merge in info scores (from full UKBB sample)
    ukb_mfi = hail.read_table(GWAS_mfi).repartition(500, shuffle=False)
    ukb_qc2 = ukb_hm3.join(ukb_mfi.select(ukb_mfi.locus, ukb_mfi.alleles, ukb_mfi.info))
    
    # merge in MAF from the UKBB GWAS sample
    ukb_qc = hail.read_table(GWAS_qc).repartition(500, shuffle=False)
    ukb_qc2 = ukb_qc2.join(ukb_qc.select(ukb_qc.locus, ukb_qc.alleles, ukb_qc.variant_qc.AF))
    
    # filter to MAF > .01, info > 0.9, not strand-ambi
    ukb_qc2 = ukb_qc2.filter(
                    (ukb_qc2.info > 0.9) & 
                    (ukb_qc2.AF > 0.01) & 
                    (ukb_qc2.AF < 0.99) & 
                    (~hail.is_strand_ambiguous(ukb_qc2.alleles[0], ukb_qc2.alleles[1]))
                )

    # save passing
    ukb_qc2 = ukb_qc2.drop(ukb_qc2.sort_al)
    ukb_qc2.write(BUCKET + HM3_qcpos_stem + '.ht', overwrite=True)



if gs_missing(BUCKET + HM3_qcpos_stem + '.tsv.bgz') or force_all:

    # format tsv version
    ukb_qc2 = hail.read_table(BUCKET + HM3_qcpos_stem + '.ht')
    ukb_qc2 = ukb_qc2.annotate(A1=ukb_qc2.alleles[0],A2=ukb_qc2.alleles[1])
    ukb_qc2 = ukb_qc2.select(ukb_qc2.locus,ukb_qc2.A1,ukb_qc2.A2,ukb_qc2.rsid,ukb_qc2.varid,ukb_qc2.variant,ukb_qc2.AF,INFO=ukb_qc2.info)
    
    # save
    ukb_qc2.export(BUCKET + HM3_qcpos_stem + '.tsv.bgz')



####
#
# Check success and exit
#
####
    
if not gs_missing(BUCKET + HM3_qcpos_stem + '.ht'):

    end_message='''
    ##############
    Finished!!
    
    HapMap3 source:
    {}
    
    HapMap3 raw vcf:
    {}
    
    HapMap3 as Hail MatrixTable: 
    {}
    
    HapMap3 autosomal, biallelic SNPs with freqs keytable:
    {}
    
    HapMap3 SNPs passing UKBB GWAS QC:
    {}
    {}
    ##############
    '''.format( HM3_FTP + HM3_NAME,
                BUCKET + HM3_NAME,
                BUCKET + HM3_vds,
                BUCKET + HM3_bi_af,
                BUCKET + HM3_qcpos_stem + '.ht',
                BUCKET + HM3_qcpos_stem + '.tsv.bgz'
            )
    
    print(end_message)

else:
    
    print("!!! Error: Failed to write? !!!")

# eof
