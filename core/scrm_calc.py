import allel
import numpy as np
import re 
import sys
import os
import pandas as pd
# TODO:switch to sgkit  <22-04-22, yourname> #
argv={x.split('=')[0]: x.split('=')[1] for x in sys.argv[1:]}
locus_write=eval(argv['locus_write'])
global_write=eval(argv['global_write'])
locus_datafile=argv['locus_datafile']
num_dataset=int(argv['num_dataset'])
min_sites = 1 #TODO add a warning about this
locus_data = pd.read_csv(locus_datafile,sep='\t')
locus_size_it = iter(locus_data['locus_length'])
subpop_it = ([range(0,locus_data.loc[i,'size_popA']),range(locus_data.loc[i,'size_popA'],(locus_data.loc[i,'size_popA'] + locus_data.loc[i,'size_popB']))] for i in range(locus_data.shape[0]))
### control ###
log = open('error.log','w')
loc_num = iter(range(locus_data.shape[0]))
#######
def get_GT(locus_data):
    pos_list = []
    haplotype_list = []
    for snp in locus_data:
        snp = snp.split(' ')
        pos_list.append(np.floor(float(snp[0])))
        haplotype_list.append(np.array(snp[2:],dtype=int))
    return {'pos':np.array(pos_list,dtype=int) , 'GT':haplotype_list}

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

def get_abcstat(gt,locus_length,subpop,seed):
    # when biased is true, all statistics will be calculated infering missing position as 
    # invariant sites. If biased == False, stats are calculated only on given SNP.
        l_n = next(loc_num)
        nsites = len(gt['pos'])
        if nsites > min_sites:
            pos = gt['pos']
            h = allel.HaplotypeArray(gt['GT'])
            acA  = h.count_alleles(subpop=subpop[0])
            acB  = h.count_alleles(subpop=subpop[1])
            piA = allel.sequence_diversity(pos,acA,start=1,stop=locus_length)
            piB = allel.sequence_diversity(pos,acB,start=1,stop=locus_length)
            dxy = allel.sequence_divergence(pos,acA,acB,start=1,stop=locus_length)
            #w=np.tile(pos,(1,2)).reshape(len(pos),2,order='F')
            da = dxy - (piA + piB)/2
            thetaA = allel.watterson_theta(pos,acA,start=1,stop=locus_length)
            thetaB = allel.watterson_theta(pos,acB,start=1,stop=locus_length)
            TajDA = allel.tajima_d(acA,pos,start=1,stop=locus_length)
            TajDB = allel.tajima_d(acB,pos,start=1,stop=locus_length)
            sfs = allel.joint_sfs(acA[:,1],acB[:,1],len(subpop[0]),len(subpop[1]))
            sxA = np.sum(sfs[1:-1,(0,-1)])/nsites
            sxB = np.sum(sfs[(0,-1),1:-1])/nsites
            sf = (sfs[-1,0] + sfs[0,-1])/nsites
            ss = np.sum(sfs[1:-1,1:-1])/nsites
            num,den=allel.hudson_fst(acA,acB)
            Fst = np.nansum(num)/np.nansum(den)
            return {'bialsite_avg':nsites,'piA_avg':piA,'piB_avg':piB,'divAB_avg':dxy,
                    'netDivAB_avg':da,'thetaA_avg':thetaA,
                    'thetaB_avg':thetaB,'DtajA_avg':TajDA,'DtajB_avg':TajDB,'sxA_avg':sxA,
                    'sxB_avg':sxB,'sf_avg':sf,'ss_avg':ss,'FST_avg':Fst,'seed':seed}
        else :
            log.write(str(l_n) + '\n')

locus_nb = 0
in_locus = False
in_head_locus = False
locus_str = []
stat = []
header ='position'
for line in sys.stdin:
    line = line.strip()

    if in_head_locus == False and re.search('scrm',line) and in_locus==False :
        in_head_locus =  True
        continue
    elif in_head_locus == True and in_locus==False :
        seed = line 
        in_head_locus = False
        continue
    elif in_locus == False and in_head_locus==False  and not re.search(header,line):
        continue
    elif in_locus == False and re.search(header,line):
        in_locus = True
        locus_str = []
        continue
    elif in_locus == True and re.search('scrm',line):
        gt = get_GT(locus_str)
        stat.append(get_abcstat(gt,next(locus_size_it),next(subpop_it),seed))
        in_locus = False
        in_head_locus = True
    else:
        locus_str.append(line)
gt = get_GT(locus_str)
stat.append( get_abcstat(gt,next(locus_size_it),next(subpop_it),seed))
index = np.arange(len(stat))
index = np.delete(index,np.where(np.array(stat)==None)[0])
stat = list(filter(None,stat))
locus_stat = pd.DataFrame(stat)
locus_stat.replace(np.inf,0)

log.close()

if locus_write:
    locus_stat.set_index(index + num_dataset*locus_stat.shape[0],inplace=True)
    if os.path.isfile('ABCstat_locus.txt') == False:
        locus_stat.to_csv('ABCstat_locus.txt',sep='\t',header=True,index_label='dataset',float_format='%.5f',mode='a',na_rep='0')
    else:
        locus_stat.to_csv('ABCstat_locus.txt',sep='\t',header=None,index_label='dataset',float_format='%.5f',mode='a',na_rep='0')

if global_write:
    locus_stat.drop('seed',axis=1,inplace=True)
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
        'sf_outlier': get_outlier(locus_stat['sf_avg']),'piA_outlier':get_outlier(locus_stat['piA_avg'],method='inf'),'piB_outlier':get_outlier(locus_stat['piB_avg'],method='inf')}
    glob_stat = pd.DataFrame(pd.concat([avg,std,med,pd.Series(pearson),pd.Series(ss_sf),pd.Series(outlier)])).T
    glob_stat.set_index(glob_stat.index + num_dataset,inplace=True)
    if os.path.isfile('ABCstat_global.txt') == False:
        glob_stat.to_csv('ABCstat_global.txt',sep='\t',header=True,index_label='dataset',float_format='%.10f',mode='a',na_rep='0')
    else:
        glob_stat.to_csv('ABCstat_global.txt',sep='\t',header=False,index_label='dataset',float_format='%.10f',mode='a',na_rep='0')


