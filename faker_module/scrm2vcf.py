import allel
import numpy as np
import re 
import sys
import os
import pandas as pd
argv={x.split('=')[0]: x.split('=')[1] for x in sys.argv[1:]}
locus_datafile=argv['locus_datafile']
locus_data = pd.read_csv(locus_datafile,sep='\t')
samples_names =  [ 'A{}'.format(x) for x in range(locus_data.loc[0,'size_popA']) ] + [ 'B{}'.format(x) for x in range(locus_data.loc[0,'size_popB'])]

def get_vcfGT(locus_data,contig_nb):
    pos_list = []
    haplotype_list = []
    for snp in locus_data:
        snp = snp.split(' ')
        pos_list.append(np.floor(float(snp[0])))
        haplotype_list.append(np.array(snp[2:],dtype=int))
    nsnp = len(pos_list)
    snp_info = pd.DataFrame({'#CHROM':np.full(nsnp,'contig_{}'.format(contig_nb)),'POS':np.array(pos_list,dtype=int),
            'ID':np.full(nsnp,'fake'),'REF':np.full(nsnp,'A'),'ALT':np.full(nsnp,'T'),
            'QUAL':np.full(nsnp,'.'),'FILTER':np.full(nsnp,'PASS'),'INFO':np.full(nsnp,'.'),
            'FORMAT':np.full(nsnp,'GT')})
    gt = pd.DataFrame(haplotype_list)
    gt.columns = samples_names
    return pd.concat([snp_info,gt],axis=1)

contig_nb = -1
in_locus = False
locus_str = []
vcf_list = []
header ='position'
for line in sys.stdin:
    line = line.strip()
    if in_locus == False and not re.search(header,line):
        continue
    elif in_locus == False and re.search(header,line):
        in_locus = True
        locus_str = []
        contig_nb +=1
        continue
    elif (in_locus == True) and re.search('scrm',line):
        vcf_list.append(get_vcfGT(locus_str,contig_nb))
        in_locus = False
    else:
        locus_str.append(line)

vcf_list.append(get_vcfGT(locus_str,contig_nb))
vcf=pd.concat(vcf_list)
vcf.to_csv('simulation.vcf',sep='\t',header=True,index=False)
