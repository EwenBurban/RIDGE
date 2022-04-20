import sys
import os
from numpy.random import uniform
from numpy.random import binomial
from numpy.random import beta
from random import shuffle
from random import randint
from numpy import multiply
from pandas import DataFrame
from numpy import repeat
from pandas import concat
import itertools as it
import numpy as np
import pandas as pd
argv={x.split('=')[0]: x.split('=')[1] for x in sys.argv[1:]}
nMultilocus=int(argv['nMultilocus'])
model=argv['model']
bpfile=argv['bpfile']
locus_write=eval(argv['locus_write'])
global_write=eval(argv['global_write'])
posteriors = argv['posteriors']
binpath=argv['binpath']
gene_flow_barrier_proportion = False # relique of old script, allow the script to test various values of nBarriers
barrier_simmetry= False # relique of old script, force simulation to test the case of non-symmetric migration rate.  

test_bpfile = os.path.isfile(bpfile)

if test_bpfile == False:
    sys.exit('\n\tERROR in submit_simulations_2pop.py : the file {0}/bpfile is not found\n'.format(path))
else:
    infile = open(bpfile, 'r')
    tmp = infile.readline()
    tmp = infile.readline().strip().split('\t')
    nlocus = len(tmp)
    size_popA = int(infile.readline().strip().split('\t')[0])
    size_popB = int(infile.readline().strip().split('\t')[0])
    totpopsize= size_popA + size_popB
    infile.close()

first_part_command = 'scrm {totpopsize} 1 '.format(totpopsize=totpopsize)
mscommand = ""
if "SI" in model:
    mscommand = "-t {} -r {} {} -I 2 {} {} 0 -n 1 {} -n 2 {} -en {} 1 {} -en {} 2 {} -ej {} 2 1 -eN {} {}"
if "AM" in model:
    mscommand = "-t {} -r {} {} -I 2 {} {} 0 -n 1 {} -n 2 {} -en {} 1 {} -en {} 2 {} -ema {} 0 {} {} 0 -ej {} 2 1 -eN {} {}"
if "SC" in model:
    mscommand = "-t {} -r {} {} -I 2 {} {} 0 -m 1 2 {} -m 2 1 {} -n 1 {} -n 2 {} -en {} 1 {} -en {} 2 {} -eM {} 0 -ej {} 2 1 -eN {} {}"
if "IM" in model:
    mscommand = "-t {} -r {} {} -I 2 {} {} 0 -n 1 {} -n 2 {} -en {} 1 {} -en {} 2 {} -m 1 2 {} -m 2 1 {} -ej {} 2 1 -eN {} {}"
if "PAN" in model:
    mscommand = "-t {} -r {} {} -eN 0 {} -eN {} {}"

if mscommand == "":
    print("You specified a wrong model: SI_x, AM_x, AM_x or SC_x\n")
    sys.exit()


shape_bound = [0.1, 5]


################ read bpfile  ##################
infile = open(bpfile, "r")
tmp = infile.readline()
col_names = ['L','nsamA','nsamB','theta','rho'] # the list order must be keeped due to the bpfile build method
locus_dic = { col: infile.readline().strip().split("\t") for col in col_names }
locus_df = DataFrame(locus_dic)
locus_df = locus_df.astype(float)
infile.close()
# sum of nsamA + nsamB
locus_df.insert(column='nsam_tot',value=locus_df[['nsamA','nsamB']].astype(int).sum(axis=1),loc=0)
# number of loci
nLoci = locus_df.shape[0]



###### function needed for parameter generation ####
def asint(l):
    _ = [int(x) for x in l]
    return(_)

def asfloat(l):
    _ = [float(x) for x in l]
    return(_)

def flatter(nested_list):
    flat_list = []
    for sublist in nested_list:
            for item in sublist:
                flat_list.append(item)
    return flat_list

def short(x,digits=5):
    if type(x) == np.dtype('float64'):
        return round(x,digits)
    else:
        return x



def beta_dis(X,a,b): # In case of modeBarrier == beta, apply the a and b values to a vector X of values
    scalar = beta(a,b)
    rescale = a / (a +b)
    _ = X * scalar/rescale
    return(_)

def build_locusDf(param,locus_df):
    nLoci = locus_df.shape[0]
    locus_sim = pd.DataFrame([param for x in range(nLoci)]) # repeat the param line nloci times into a DF
    locus_sim.reset_index(inplace=True,drop=True)
    locus_sim = pd.concat([locus_sim,locus_df],axis=1)
    if 'shape_N_a' in param :
        locus_sim[['Na','N1','N2']] = locus_sim[['Na','N1','N2']].apply(lambda x: beta_dis(x,param['shape_N_a'],param['shape_N_b']),axis=1)
    locus_sim[['founders1','founders2']] = locus_sim[['founders1','founders2']].multiply(locus_sim['Na'],axis=0)
    if 'nBarriersM12' in param:
        locus_sim['M12'] = locus_sim['M12'].multiply(np.random.choice([0,1],nLoci,p= [param['nBarriersM12']/nLoci,1-(param['nBarriersM12']/nLoci)]),axis=0)
        locus_sim['M21'] = locus_sim['M21'].multiply(np.random.choice([0,1],nLoci,p= [param['nBarriersM21']/nLoci,1-(param['nBarriersM21']/nLoci)]),axis=0)
        if 'shape_M12_a' in param : 
            locus_sim['M12'] = locus_sim['M12'].apply(lambda x: beta_dis(x,param['shape_M12_a'],param['shape_M12_b']))
            locus_sim['M21'] = locus_sim['M21'].apply(lambda x: beta_dis(x,param['shape_M21_a'],param['shape_M21_b']))
    tmp=locus_sim[['Na','N1','N2','founders1','founders2']].copy()
    tmp[tmp<1e-4]=1e-4
    locus_sim.update(tmp)

    return locus_sim

param_df = pd.read_csv(posteriors,sep='\t') 
if param_df.shape[0] != nMultilocus :
    param_df = param_df.loc[np.random.choice(range(param_df.shape[0]),nMultilocus),:]
    param_df.reset_index(inplace=True,drop=True)
param_df[param_df < 0 ] = 1e-4

if 'nBarriersM21' in param_df.columns:
    tmp = param_df[['nBarriersM12','nBarriersM21']].copy()
    tmp[tmp > nLoci] = nLoci
    param_df.update(tmp)
if 'SC' not in model and 'Tsc' in param_df.columns :
    del param_df['Tsc']
if 'AM' not in model and 'Tam' in param_df.columns :
    del param_df['Tam']

########################### Generate locus level parameters

locus_param_df = [build_locusDf(param_df.loc[x,:],locus_df) for x in range(nMultilocus)]

## define the ms argument order depending of the model
base_prior = ['theta','rho','L','nsamA','nsamB'] # this part is common to all models
if 'IM' in model:
    locus_prior_names = base_prior + ['N1','N2','Tdem1','founders1','Tdem2','founders2','M12','M21','Tsplit','Tsplit','Na']
elif 'SC' in model:
    locus_prior_names = base_prior + ['M12','M21','N1','N2','Tdem1','founders1','Tdem2','founders2','Tsc','Tsplit','Tsplit','Na']
elif 'AM' in model:
    locus_prior_names = base_prior + ['N1','N2','Tdem1','founders1','Tdem2','founders2','Tam','M12','M21','Tsplit','Tsplit','Na']
elif 'SI' in model:
    locus_prior_names = base_prior + ['N1','N2','Tdem1','founders1','Tdem2','founders2','Tsplit','Tsplit','Na']


################### write the priorfiles and pipe the ms argument through print()
# write priorfile.txt which contains global simulation parameters
if global_write == True:
    with open('priorfile.txt','w') as o:
        o.write(param_df.to_csv(sep="\t",header=True,index_label='dataset',float_format='%.5f'))

if locus_write == True:
    with open('priorfile_locus.txt','w') as lo:
        lo.write(locus_param_df[0].to_csv(sep="\t",header=True,index_label='dataset',float_format='%.5f')) # only the header of the first df is use
        for df in locus_param_df[1:] :
            lo.write(df.to_csv(sep="\t",header=False,index_label='dataset',float_format='%.5f'))

num_dataset=0
o=open('exec.sh','w')
for sim in locus_param_df:
    output_command = '{ '
    sim.round(5)
    for locus in range(nLoci):
        locus_prior_str = [str(short(i)) for i in sim.loc[locus,locus_prior_names]]
        mscommand_formated =  mscommand.format(*locus_prior_str) 
        output_command = output_command +  first_part_command + mscommand_formated + " ;"
    output_command += '}' + ' | pypy {binpath}/mscalc_new.py bpfile={bpfile} num_dataset={num_dataset} \n'.format(binpath=binpath,bpfile=bpfile,num_dataset=num_dataset) 
    o.write(output_command)
    num_dataset += 1
o.close()
