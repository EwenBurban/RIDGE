import allel
import sys
import pandas as pd
import numpy as np
argv={x.split('=')[0]: x.split('=')[1] for x in sys.argv[1:]}
data_file = argv['data']
contig_file = argv['contig_file']
popfile = argv['popfile']
nameA = argv['nameA']
nameB = argv['nameB']
win_size = int(argv['window_size'])
min_sites = 1 #TODO add a warning about this
output=argv['output']
print('Warning: window containing less than {} sites are rejected'.format(min_sites))
### load data ####
contig_data = pd.read_csv(contig_file,sep='\t')
popfile = pd.read_csv(popfile)
data = allel.read_vcf(data_file)
gt = allel.GenotypeArray(data['calldata/GT'])
popA_samples = popfile[nameA]
popB_samples = popfile[nameB]
popA_index=list(np.where(np.in1d(data['samples'],popA_samples)==True)[0])
popB_index=list(np.where(np.in1d(data['samples'],popB_samples)==True)[0])
acA = gt.count_alleles(subpop=popA_index)
acB = gt.count_alleles(subpop=popB_index)

### process data ###
contig_list = list(contig_data['contig_name'])
chr_stat_list = []
for contig in contig_list:
    sel_snp = np.where(data['variants/CHROM'] == contig)
    contig_length = contig_data[contig_data['contig_name'] == contig]['contig_length'].to_list()[0]
    pos = data['variants/POS'][sel_snp]
    sub_acA = acA[sel_snp]
    sub_acB = acB[sel_snp]
    window_list = allel.windowed_diversity( pos,sub_acA,size=win_size,start=1,stop=contig_length)[1]

    sxA_arr = sxB_arr= sf_arr = ss_arr= nsites = thetaA= thetaB = TajDA = TajDB = Fst = piA = piB = dxy = da = dataset = np.empty(0)
    for window in window_list:
        sel_snp_sfs = np.where(np.logical_and(pos>=window[0],pos<=window[1])==True)[0]
        sfs_nsites = len(sel_snp_sfs)
        if sfs_nsites > min_sites:
            piA_tmp = allel.sequence_diversity(pos,sub_acA,start=window[0],stop=window[1])
            piB_tmp= allel.sequence_diversity(pos,sub_acB,start=window[0],stop=window[1])
            dxy_tmp = allel.sequence_divergence(pos,sub_acA,sub_acB,start=window[0],stop=window[1])
            da_tmp = dxy_tmp - (piA_tmp + piB_tmp)/2
            TajDA_tmp = allel.tajima_d(sub_acA,pos,start=window[0],stop=window[1])
            TajDB_tmp = allel.tajima_d(sub_acB,pos,start=window[0],stop=window[1])
            Fst_tmp = allel.average_hudson_fst(sub_acA[sel_snp_sfs,:],sub_acB[sel_snp_sfs,:],sfs_nsites)[0]
            thetaA_tmp = allel.watterson_theta(pos,acA,start=window[0],stop=window[1])
            thetaB_tmp =  allel.watterson_theta(pos,acB,start=window[0],stop=window[1])
            sfs = allel.joint_sfs(acA[sel_snp_sfs,1],acB[sel_snp_sfs,1])
            sxA = np.sum(sfs[1:-1,(0,-1)])/sfs_nsites
            sxB = np.sum(sfs[(0,-1),1:-1])/sfs_nsites
            sf = (sfs[-1,0] + sfs[0,-1])/sfs_nsites
            ss = np.sum(sfs[1:-1,1:-1])/sfs_nsites
            dataset_tmp = '{}:{}-{}'.format(contig,window[0],window[1])

            nsites = np.append(nsites,sfs_nsites)
            sxA_arr = np.append(sxA_arr,sxA)
            sxB_arr = np.append(sxB_arr,sxB)
            sf_arr = np.append(sf_arr,sf)
            ss_arr = np.append(ss_arr,ss)
            thetaA = np.append(thetaA,thetaA_tmp)
            thetaB = np.append(thetaB,thetaB_tmp)
            piA = np.append(piA,piA_tmp)
            piB = np.append(piB,piB_tmp)
            TajDA = np.append(TajDA,TajDA_tmp)
            TajDB = np.append(TajDB,TajDB_tmp)
            dxy = np.append(dxy,dxy_tmp)
            da = np.append(da,da_tmp)
            Fst = np.append(Fst,Fst_tmp)
            dataset = np.append(dataset,dataset_tmp)
        

    #### format all arrays into a dataframe
    contig_stat= pd.DataFrame({'dataset':dataset,'bialsite_avg':nsites,'piA_avg':piA,'piB_avg':piB,'divAB_avg':dxy,
            'netDivAB_avg':da,'thetaA_avg':thetaA,
            'thetaB_avg':thetaB,'DtajA_avg':TajDA,'DtajB_avg':TajDB,'sxA_avg':sxA_arr,
            'sxB_avg':sxB_arr,'sf_avg':sf_arr,'ss_avg':ss_arr,'FST_avg':Fst})
    chr_stat_list.append(contig_stat)
    print('conting {} done'.format(contig))

locus_stat = pd.concat(chr_stat_list,axis=0)
locus_stat.set_index('dataset',inplace=True)
locus_stat.to_csv(output,sep='\t',header=True,index_label='dataset',float_format='%.5f',mode='w',na_rep='0')
