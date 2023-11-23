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
#### read locus_datafile
locus_data = pd.read_csv(locus_datafile,sep='\t')
nLoci = locus_data.shape[0]

mscommand = "scrm {totpopsize} 1 -t {theta} -r {rho} {locus_length} -l 100r -I 2 {size_popA} {size_popB} {M_current} -n 1 {N1} -n 2 {N2}  -ema {Tsc} 0 0 0 0 -ema {Tam} 0 {M_ancestral} {M_ancestral} 0 -ej {Tsplit} 2 1 -eN {Tsplit} {Na} -ema {Tsplit} 0 0 0 0  -seed {seed}--print-model -transpose-segsites -SC abs"

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
    locus_sim['seed'] = np.random.randint(0,high = 1e18, size = nLoci) # add unique seed number to each locus simulation. It’s used as a tag to avoid confusion and mismatch during reference table
    N = ['Na','N1','N2']
    locus_sim[N] = locus_sim[N].apply(lambda x: beta_dis(x,param['shape_N_a'],param['shape_N_b']),axis=1)
    migration= 'M_current'
    locus_sim[migration] = locus_sim[migration].multiply(np.random.choice([0,1],nLoci,p= [param['Pbarrier'+ migration ],1-param['Pbarrier'+ migration ]]),axis=0)
    migration= 'M_ancestral'
    locus_sim[migration] = locus_sim[migration].multiply(np.random.choice([0,1],nLoci,p= [param['Pbarrier'+ migration ],1-param['Pbarrier'+ migration ]]),axis=0)
    return locus_sim

# For each multilocus dataset, transform it in a dataframe containing prior for each locus in the multilocus dataset
locus_param_df = [build_locusDf(glob_prior.loc[x,:],locus_data,nLoci) for x in range(nMultilocus)]


################### write the priorfiles and  the ms commands 
# write priorfile.txt which contains global simulation parameters
if global_write == True:
    with open('priorfile.txt','w') as o:
        o.write(glob_prior.to_csv(sep="\t",header=True,index_label='dataset',float_format='%.5f'))

if locus_write == True:
    with open('priorfile_locus.txt','w') as lo:
        locus_param_df_full = pd.concat(locus_param_df,axis=0)
        locus_param_df_full.reset_index(drop=True,inplace=True)
        lo.write(locus_param_df_full.to_csv(sep="\t",header=True,index_label='dataset',float_format='%.5f'))

def short(x,digits=5): # In order to limit storage cost, each numerical value is shorten to 5 digits
    if type(x) == np.dtype('float64'):
        y = round(x,digits)
        if y <= 0: # in case of negative value, set to 1e-5. Because negative value should not happen. 
            y = 10 ** -digits
        return y
    else:
        return x

num_dataset=0
with open('exec.sh','w') as o:
    for sim in locus_param_df:
        sim = sim.apply(lambda x: [str(short(y)) for y in x ],axis=0)
        output_command = '{ '
        for locus in range(nLoci):
            mscommand_formated =  mscommand.format(**sim.loc[locus,:]) 
            output_command = output_command +  mscommand_formated + " ;"
        output_command += '}' + ' | python3 {binpath}/scrm_calc.py locus_datafile={locus_datafile} num_dataset={num_dataset} locus_write=True global_write=False \n'.format(binpath=binpath,locus_datafile=locus_datafile,num_dataset=num_dataset) 
        num_dataset += 1
        o.write(output_command)
