import sys
import os
import pandas as pd
import itertools as it
import numpy as np
import re
## Goal of the script
# It produce command line in the file exec.sh, when executed will generate coalescent simulation through scrm tool and then the result of simulation will be converted by the
# script named scrm_calc.py into summary statistics stored in ABCstat_global.txt [if {global_write}==True] or/and ABCstat_locus.txt [if {locus_write}==True]
# This script use parameters specified in {priorfile}
argv={x.split('=')[0]: x.split('=')[1] for x in sys.argv[1:]}
nMultilocus=int(argv['nMultilocus'])
config_yaml = argv['config_yaml']
locus_datafile=argv['locus_datafile']
locus_write=eval(argv['locus_write'])
global_write=eval(argv['global_write'])
binpath=argv['binpath']
priorfile=argv['priorfile']
Nref=float(argv['Nref'])
#### read locus_datafile
locus_data = pd.read_csv(locus_datafile,sep='\t')
nLoci = locus_data.shape[0]

mscommand = "scrm {totpopsize} 1 -t {theta} -r {rho} {locus_length} -l 100r -I 2 {size_popA} {size_popB} {M_current} -m 1 2 {M12_current} -m 2 1 {M21_current} -n 1 {N1} -n 2 {N2}  -ema {Tsc} 0 0 0 0 -ema {Tam} 0 {M12_ancestral} {M21_ancestral} 0 -ej {Tsplit} 2 1 -eN {Tsplit} {Na} -ema {Tsplit} 0 0 0 0  -seed {seed}--print-model -transpose-segsites -SC abs"

#### read priorfile and correct for aberation and remove useless prior
#glob_prior = pd.read_csv(priorfile,sep='\t',index_col='dataset') 
glob_prior = pd.read_csv(priorfile,sep='\t') 
glob_prior.reset_index(inplace=True,drop=True)
if glob_prior.shape[0] != nMultilocus :
    glob_prior = glob_prior.loc[np.random.choice(range(glob_prior.shape[0]),nMultilocus),:]
    glob_prior.reset_index(inplace=True,drop=True)
glob_prior[glob_prior < 0 ] = 1e-5

glob_prior['Tsc'] = glob_prior.apply(lambda x: x['Tsc'] if x['Tsc'] < x['Tsplit'] else x['Tsplit'],axis=1)
glob_prior['Tam'] = glob_prior.apply(lambda x: x['Tam'] if x['Tam'] < x['Tsplit'] else x['Tsplit'],axis=1)
## for locus dataset, Pbarrier is fixed to 0.5 to generate a dataset with an half under barrrier status and the other hald under non-barrier status
## in order feed well the random forest with the same amount of information for both classes
for migration in ['M_current','M_ancestral'] : 
    glob_prior['Pbarrier' + migration ]=0.5
########################### Generate locus level parameters
def beta_dis(X,a,b): # Create rescale beta distribution
    scalar = np.random.beta(a,b)
    rescale = a / (a +b)
    _ =  X * scalar/rescale 
    return(_)


def build_locusDf(param,locus_df,nLoci):# This function apply the genomic mode defined before to create heterogeneity among locus
    locus_sim = pd.DataFrame([param for x in range(nLoci)]) # repeat the param line nLoci times into a DF
    locus_sim.reset_index(inplace=True,drop=True)
    locus_sim = pd.concat([locus_sim,locus_df],axis=1)
    N = ['Na','N1','N2']
    locus_sim[N] = locus_sim[N].apply(lambda x: beta_dis(x,param['shape_N_a'],param['shape_N_b']),axis=1)
    migration= 'M_current'
    locus_sim[migration] = locus_sim[migration].multiply(np.random.choice([0,1],nLoci,p= [param['Pbarrier'+ migration ],1-param['Pbarrier'+ migration ]]),axis=0)
    locus_sim['M12_current'] = locus_sim['N1'] *4* Nref * locus_sim[migration]
    locus_sim['M21_current']= locus_sim['N2'] *4* Nref * locus_sim[migration]
    migration= 'M_ancestral'
    locus_sim[migration] = locus_sim[migration].multiply(np.random.choice([0,1],nLoci,p= [param['Pbarrier'+ migration ],1-param['Pbarrier'+ migration ]]),axis=0)
    locus_sim['M12_ancestral'] = locus_sim['N1'] *4* Nref * locus_sim[migration]
    locus_sim['M21_ancestral']= locus_sim['N2'] *4* Nref * locus_sim[migration]
    return locus_sim

# For each multilocus dataset, transform it in a dataframe containing prior for each locus in the multilocus dataset
locus_param_df = [build_locusDf(glob_prior.loc[x,:],locus_data,nLoci) for x in range(nMultilocus)]
locus_param_df= pd.concat(locus_param_df,axis=0)
prob_vec_barrier=np.where(locus_param_df['M_current']==0,1,0)
prob_vec_barrier=prob_vec_barrier/np.sum(prob_vec_barrier)
prob_vec_non_barrier=np.where(locus_param_df['M_current']!=0,1,0)
prob_vec_non_barrier=prob_vec_non_barrier/np.sum(prob_vec_non_barrier)
barrier_param_df = locus_param_df.iloc[np.random.choice(range(locus_param_df.shape[0]),int(locus_param_df.shape[0]/2),p=prob_vec_barrier,replace=True),:]
non_barrier_param_df = locus_param_df.iloc[np.random.choice(range(locus_param_df.shape[0]),int(locus_param_df.shape[0]/2),p=prob_vec_non_barrier,replace=True),:]
locus_param_df_full=pd.concat([barrier_param_df,non_barrier_param_df],axis=0)
locus_param_df_full['seed']=np.random.randint(0,high = 1e18, size = locus_param_df_full.shape[0])# add unique seed number to each locus simulation. It’s used as a tag to avoid confusion and mismatch during reference table
print(locus_param_df_full)
################### write the priorfiles and  the ms commands 
# write priorfile.txt which contains global simulation parameters
#if global_write == True:
 #   with open('priorfile.txt','w') as o:
  #      o.write(glob_prior.to_csv(sep="\t",header=True,index_label='dataset',float_format='%.10f'))

if locus_write == True:
    with open('priorfile_locus.txt','w') as lo:
        locus_param_df_full.reset_index(drop=True,inplace=True)
        lo.write(locus_param_df_full.to_csv(sep="\t",header=True,index_label='dataset',float_format='%.10f'))

def short(x,digits=10): # In order to limit storage cost, each numerical value is shorten to 10 digits
    if type(x) == np.dtype('float64'):
        y = round(x,digits)
        if y <= 0: # in case of negative value, set to 1e-5. Because negative value should not happen. 
            y = 10 ** -digits
        return y
    else:
        return x

num_dataset=0
locus_param_df_full=locus_param_df_full.apply(lambda x: [str(short(y)) for y in x ],axis=0)
locus_param_df_full['tag']=np.repeat(np.arange(nMultilocus),nLoci)
split_locus_param=dict(tuple(locus_param_df_full.groupby('tag')))
with open('exec.sh','w') as o:
    for sim in split_locus_param:
        sim= pd.DataFrame(split_locus_param[sim])
        sim.reset_index(drop=True,inplace=True)
        output_command = '{ '
        for locus in range(nLoci):
            mscommand_formated =  mscommand.format(**sim.loc[locus,:]) 
            output_command = output_command +  mscommand_formated + " ;"
        output_command += '}' + ' | python3 {binpath}/scrm_calc.py locus_datafile={locus_datafile} num_dataset={num_dataset} locus_write=True global_write=False \n'.format(binpath=binpath,locus_datafile=locus_datafile,num_dataset=num_dataset) 
        num_dataset += 1
        o.write(output_command)
