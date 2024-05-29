import allel
import sys
import pandas as pd
import numpy as np
import time

argv={x.split('=')[0]: x.split('=')[1] for x in sys.argv[1:]}
data_file = argv['data']
bed_file = argv['bed_file']
popfile = argv['popfile']
nameA = argv['nameA']
nameB = argv['nameB']
mask_file = argv['mask_file']
win_size = int(argv['window_size'])
min_sites = 3 
locus_write = eval(argv['locus_write'])
global_write = eval(argv['global_write'])
output_dir=argv['output_dir']
ploidy=int(argv['ploidy'])
print('Warning: window containing less than {} sites are rejected'.format(min_sites))
### load data ####
bed = pd.read_csv(bed_file,sep='\t')
popfile = pd.read_csv(popfile,dtype=str)
if mask_file != 'None':
    mask_file = pd.read_csv(mask_file,sep='\t')


data = allel.read_vcf(data_file)
gt = allel.GenotypeArray(data['calldata/GT'])
popA_samples = popfile[nameA]
popB_samples = popfile[nameB]
popA_index=list(np.where(np.in1d(data['samples'],popA_samples)==True)[0])
popB_index=list(np.where(np.in1d(data['samples'],popB_samples)==True)[0])
acA = gt.count_alleles(subpop=popA_index)
acB = gt.count_alleles(subpop=popB_index)

def get_outlier(vec,method='supp'):
    vec=np.arcsin(np.sqrt(vec))
    Q1=np.quantile(vec,0.25)
    Q3=np.quantile(vec,0.75)
    IQR=Q3-Q1
    if method == 'supp' : 
        upper_outlier_born=Q3+1.5*IQR
        outlier=np.where(vec>upper_outlier_born)
    elif method == 'inf' : 
        under_outlier_born=Q1-1.5*IQR
        outlier=np.where(vec<under_outlier_born)
    return len(outlier[0])/len(vec)


### process data ###
chr_stat_list = []
for contig in set(bed['chr']):
    sel_snp = np.where(data['variants/CHROM'] == contig)
    pos = data['variants/POS'][sel_snp]
    window_list  =np.array(bed[bed['chr']==contig][['start','end']])
    contig_accessibility=mask_file.loc[(mask_file['chr'] == contig),'accessibility'].to_numpy() 
    contig_pos_access=mask_file.loc[(mask_file['chr'] == contig),'pos'].to_numpy()
    contig_accessibility=np.array([contig_accessibility[i] for i in np.argsort(contig_pos_access)])


    sxA_arr = sxB_arr= sf_arr = ss_arr= nsites = thetaA= thetaB = TajDA = TajDB = Fst = piA = piB = dxy = da = dataset =np.empty(0)
    for window in window_list:
        sel_snp_sfs = sel_snp[0][np.where(np.logical_and(pos>=window[0],pos<=window[1])==True)[0]]
        sfs_nsites = len(sel_snp_sfs)
        sub_pos=pos[np.where(np.logical_and(pos>=window[0],pos<=window[1])==True)[0]]
        if isinstance(mask_file,pd.DataFrame): 
            accessibility=contig_accessibility[int(window[0]):(int(window[1])+1)]
        else : 
            accessibility=np.full(int(window[1]-window[0]),1)

        if sfs_nsites >= min_sites:
            relative_sub_pos=sub_pos-(int(window[0])-1)
            tpos,ac=allel.util.mask_inaccessible(accessibility,relative_sub_pos,acA[sel_snp_sfs])
            piA_tmp = allel.sequence_diversity(relative_sub_pos,acA[sel_snp_sfs],start=1,stop=int(window[1])-int(window[0])+1,is_accessible=accessibility)
            piA_tmp_bis = allel.sequence_diversity(relative_sub_pos,acA[sel_snp_sfs],start=1,stop=int(window[1])-int(window[0])+1)
            piB_tmp= allel.sequence_diversity(relative_sub_pos,acB[sel_snp_sfs],start=1,stop=int(window[1])-int(window[0])+1,is_accessible=accessibility)
            thetaA_tmp = allel.watterson_theta(relative_sub_pos,acA[sel_snp_sfs],start=1,stop=int(window[1])-int(window[0])+1,is_accessible=accessibility)
            thetaB_tmp =  allel.watterson_theta(relative_sub_pos,acB[sel_snp_sfs],start=1,stop=int(window[1])-int(window[0])+1,is_accessible=accessibility)

            TajDA_tmp = allel.tajima_d(acA[sel_snp_sfs],sub_pos,start=window[0],stop=window[1])
            TajDB_tmp = allel.tajima_d(acB[sel_snp_sfs],sub_pos,start=window[0],stop=window[1])

            dxy_tmp = allel.sequence_divergence(relative_sub_pos,acA[sel_snp_sfs],acB[sel_snp_sfs],start=1,stop=int(window[1])-int(window[0])+1,is_accessible=accessibility)
            da_tmp = dxy_tmp - (piA_tmp + piB_tmp)/2
            a,b,c=allel.weir_cockerham_fst(gt[sel_snp_sfs],[popA_index,popB_index])
            Fst_tmp = np.nansum(a)/(np.nansum(a)+np.nansum(b)+np.nansum(c))

            sfs = allel.joint_sfs(acA[sel_snp_sfs,1],acB[sel_snp_sfs,1],len(popA_index)*ploidy,len(popB_index)*ploidy)
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
if locus_write == True:
    locus_stat.to_csv('{}/ABCstat_locus.txt'.format(output_dir),sep='\t',header=True,index_label='dataset',float_format='%.5f',mode='w',na_rep='0')

if global_write == True:
    avg = locus_stat.apply(np.nanmean,axis=0)
    med = locus_stat.apply(np.nanmedian,axis=0)
    std = locus_stat.apply(np.nanstd,axis=0)
    std.index = [x.replace('_avg','_std') for x in std.keys()]
    med.index = [x.replace('_avg','_median') for x in med.keys()]
    pearson = {'pearson_r_pi':np.corrcoef(locus_stat['piA_avg'].astype(float),locus_stat['piB_avg'].astype(float)).min(),
            'pearson_r_theta':np.corrcoef(locus_stat['thetaA_avg'].astype(float),locus_stat['thetaB_avg'].astype(float)).min(),
            'pearson_r_divAB_netdivAB':np.corrcoef(locus_stat['divAB_avg'].astype(float),locus_stat['netDivAB_avg'].astype(float)).min(),
            'pearson_r_divAB_FST':np.corrcoef(locus_stat['divAB_avg'].astype(float),locus_stat['FST_avg'].astype(float)).min(),
            'pearson_r_netdivAB_FST':np.corrcoef(locus_stat['netDivAB_avg'].astype(float),locus_stat['FST_avg'].astype(float)).min()
            }
    ss_sf = {'ss_sf':np.count_nonzero(np.logical_and(locus_stat['ss_avg']>0,locus_stat['sf_avg']>0)==True),
            'ss_noSf':np.count_nonzero(np.logical_and(locus_stat['ss_avg']>0,locus_stat['sf_avg']==0)==True),
            'noSs_sf':np.count_nonzero(np.logical_and(locus_stat['ss_avg']==0,locus_stat['sf_avg']>0)==True),
            'noSs_noSf':np.count_nonzero(np.logical_and(locus_stat['ss_avg']==0,locus_stat['sf_avg']==0)==True)
            }
    outlier={'fst_outlier':get_outlier(locus_stat['FST_avg']),'divAB_outlier':get_outlier(locus_stat['divAB_avg']),'netDivAB_outlier':get_outlier(locus_stat['netDivAB_avg']),
            'sf_outlier': get_outlier(locus_stat['sf_avg']),'piA_outlier':get_outlier(locus_stat['piA_avg'],method='inf'),'piB_outlier':get_outlier(locus_stat['piB_avg'],method='inf'),'ss_outlier':get_outlier(locus_stat['ss_avg'],method='inf')}
    glob_stat = pd.DataFrame(pd.concat([avg,std,med,pd.Series(pearson),pd.Series(ss_sf),pd.Series(outlier)])).T

    glob_stat.to_csv('{}/ABCstat_global.txt'.format(output_dir),sep='\t',header=True,index_label='dataset',float_format='%.5f',mode='w',na_rep='0')
