import sys
import os
import pandas as pd
import itertools as it
import numpy as np
import re
argv={x.split('=')[0]: x.split('=')[1] for x in sys.argv[1:]}
nMultilocus=int(argv['nMultilocus'])
model=argv['model'] # add a model check
config_yaml = argv['config_yaml']
locus_datafile=argv['locus_datafile']
locus_write=eval(argv['locus_write'])
global_write=eval(argv['global_write'])
binpath=argv['binpath']
#### read locus_datafile
locus_data = pd.read_csv(locus_datafile,sep='\t')
nLoci = locus_data.shape[0]
mscommand = ""
if "SI" in model:
    mscommand = "scrm {totpopsize} 1 -t {theta} -r {rho} {locus_length} -l 100r -I 2 {size_popA} {size_popB} 0 -n 1 {N1} -n 2 {N2}  -ej {Tsplit} 2 1 -eN {Tsplit} {Na} -seed {seed} --print-model -transpose-segsites -SC abs"
if "AM" in model:
    mscommand = "scrm {totpopsize} 1 -t {theta} -r {rho} {locus_length} -l 100r -I 2 {size_popA} {size_popB} 0 -n 1 {N1} -n 2 {N2} -ema {Tam} 0 {M_ancestral} {M_ancestral} 0 -ej {Tsplit} 2 1 -eN {Tsplit} {Na} -ema {Tsplit} 0 0 0 0 -seed {seed} --print-model -transpose-segsites -SC abs"
if "SC" in model:
    mscommand = "scrm {totpopsize} 1 -t {theta} -r {rho} {locus_length} -l 100r -I 2 {size_popA} {size_popB} {M_current} -n 1 {N1} -n 2 {N2}  -ema {Tsc} 0 0 0 0 -ej {Tsplit} 2 1 -eN {Tsplit} {Na} -seed {seed}--print-model -transpose-segsites -SC abs"
if "IM" in model:
    mscommand = "scrm {totpopsize} 1 -t {theta} -r {rho} {locus_length} -l 100r -I 2 {size_popA} {size_popB} {M_current} -n 1 {N1} -n 2 {N2}  -ej {Tsplit} 2 1 -eN {Tsplit} {Na} -seed {seed} --print-model -transpose-segsites -SC abs"


shape_bound = [0.1, 10]
N_bound = [0, 0] # number of diploid individuals in the population
T_bound = [0, 0] # number of generations
M_bound = [0,0] # number of migrants per generation
P_bound = [0,0.5]
config_yaml = open(config_yaml, 'r')
for i in config_yaml:
    i = i.strip().split(':')
    if(i[0] == 'N_min'):
        N_bound[0] = float(i[1])
    if(i[0] == 'N_max'):
        N_bound[1] = float(i[1])
    if(i[0] == 'Tsplit_min'):
        T_bound[0] = float(i[1])
    if(i[0] == 'Tsplit_max'):
        T_bound[1] = float(i[1])
    if(i[0] == 'M_min'):
        M_bound[0] = float(i[1])
    if(i[0] == 'M_max'):
        M_bound[1] = float(i[1])
    if(i[0] == 'Pbarrier_max'):
        P_bound[1] = float(i[1])
    if(i[0] == 'ploidy'):
        ploidy = int(i[1])
    if(i[0] == 'Nref'):
        Nref = float(i[1])
config_yaml.close()

import numpy as np

def loguniform(low=0, high=1, size=1):
    return np.power(10, np.random.uniform(np.log10(low), np.log10(high), size))

################## convert parameter values in coalescent units
#Nref = (N_bound[1]+N_bound[0])/2.0
N_bound[0] /= Nref
N_bound[1] /= Nref
print(N_bound)
T_bound[0] /= (2*ploidy*Nref)
T_bound[1] /= (2*ploidy*Nref)
min_Tsc = 0.05
max_Tsc = 0.3
min_Tam = 0.5
###### build global priors for nMultilocus datasets ######

glob_prior = pd.DataFrame({'Tsplit': np.random.uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus),
        'Na': np.random.uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus) ,
        'N1': np.random.uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus) ,
        'N2': np.random.uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus) })
if 'SC' in model:
    glob_prior['Tsc'] = glob_prior['Tsplit'].apply(lambda x: np.random.uniform(low = min_Tsc * x, high = max_Tsc * x))
if 'AM' in model:
    migration = 'M_ancestral'
    glob_prior['Tam'] = glob_prior['Tsplit'].apply(lambda x: np.random.uniform(low = min_Tam*x, high =x))
    glob_prior[migration] = loguniform(low=M_bound[0],high=M_bound[1],size = nMultilocus) 
if 'SC' in model or 'IM' in model:
    migration = 'M_current'
    glob_prior[migration] = loguniform(low=M_bound[0],high=M_bound[1],size = nMultilocus)
if '2M' in model:
    glob_prior['shape_' + migration + '_a'] = np.random.uniform(low=shape_bound[0],high=shape_bound[1],size=nMultilocus)   
    glob_prior['shape_' + migration + '_b'] = np.random.uniform(low=shape_bound[0],high=shape_bound[1],size=nMultilocus)   
    glob_prior['Pbarrier' + migration] = np.random.uniform(low=P_bound[0],high=P_bound[1],size=nMultilocus)
if '2N' in model:
    glob_prior['shape_N_a'] = np.random.uniform(low=shape_bound[0],high=shape_bound[1],size=nMultilocus)   
    glob_prior['shape_N_b'] = np.random.uniform(low=shape_bound[0],high=shape_bound[1],size=nMultilocus)   
if '3M' in model:
    glob_prior['Pbarrier' + migration] = np.random.uniform(low=P_bound[0],high=P_bound[1],size=nMultilocus)

########################### Generate locus level parameters

def beta_dis(X,a,b): # In case of modeBarrier == beta, apply the a and b values to a vector X of values
    scalar = np.random.beta(a,b)
    rescale = a / (a +b)
    _ =  X * scalar/rescale 
    return(_)

def build_locusDf(param,locus_df,nLoci):
    locus_sim = pd.DataFrame([param for x in range(nLoci)]) # repeat the param line nLoci times into a DF
    locus_sim.reset_index(inplace=True,drop=True)
    locus_sim = pd.concat([locus_sim,locus_df],axis=1)
    locus_sim['seed'] = np.random.randint(0,high = 1e18, size = nLoci)
    if 'SI' in model:
        migration = 'null'
    else:
        pat  = re.compile('M')
        migration =  list(filter(pat.match,list(param.keys())))[0]
    if 'shape_N_a' in param :
        N = ['Na','N1','N2']
        locus_sim[N] = locus_sim[N].apply(lambda x: beta_dis(x,param['shape_N_a'],param['shape_N_b']),axis=1)
        locus_sim[N] = np.clip(locus_sim[N],9e-4,1e9)
    if 'Pbarrier'+ migration  in param:
        if 'shape_' + migration + '_a' in param : 
            locus_sim[migration] = locus_sim[migration].apply(lambda x: beta_dis(x,param['shape_' + migration + '_a'],param['shape_' + migration + '_b']))
            locus_sim[migration] = np.clip(locus_sim[migration],M_bound[0],M_bound[1])
        locus_sim[migration] = locus_sim[migration].multiply(np.random.choice([0,1],nLoci,p= [param['Pbarrier'+ migration ],1-param['Pbarrier'+ migration ]]),axis=0)
    return locus_sim

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

def short(x,digits=5):
    if type(x) == np.dtype('float64'):
        y = round(x,digits)
        if y <= 0:
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
        output_command += '}' + ' | python3 {binpath}/scrm_calc.py locus_datafile={locus_datafile} num_dataset={num_dataset} locus_write=False global_write=True\n'.format(binpath=binpath,locus_datafile=locus_datafile,num_dataset=num_dataset) 
        num_dataset += 1
        o.write(output_command)
