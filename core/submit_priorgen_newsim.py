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

argv={x.split('=')[0]: x.split('=')[1] for x in sys.argv[1:]}
nMultilocus=int(argv['nMultilocus'])
model=argv['model'] # add a model check
config_yaml = argv['config_yaml']
bpfile=argv['bpfile']
locus_write=eval(argv['locus_write'])
global_write=eval(argv['global_write'])
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
N_bound = [0, 0] # number of diploid individuals in the population
T_bound = [0, 0] # number of generations
M_bound = [0, 0] # 4.N.m , so the number of diploid migrant copies is 2.N.m
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
    if(i[0] == 'modeBarrier'):
        # is equal to "beta" or "bimodal"
        modeBarrier = i[1].replace(" ", "")
config_yaml.close()


################## convert parameter values in coalescent units
#Nref = (N_bound[1]+N_bound[0])/2.0
Nref = (0+N_bound[1])/2.0
N_bound[0] /= Nref
N_bound[1] /= Nref
T_bound[0] /= (4*Nref)
T_bound[1] /= (4*Nref)

min_Tsc = T_bound[0]
max_Tsc = 0.2
min_Tam = 0.5

################ read bpfile  ##################
infile = open(bpfile, "r")
tmp = infile.readline()
col_names = ['L','nsamA','nsamB','theta','rho'] # the list order must be keeped due to the bpfile build method
locus_dic = { col: infile.readline().strip().split("\t") for col in col_names }
locus_df = DataFrame(locus_dic)
locus_df=locus_df.astype(float)
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

def N(): # Generate 3 random global population size vector :  N1 , N2, Na
    _ = [uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus) for i in range(3)]
    return(_)

def founders(): # Generate 2 founders values vector
    _=[uniform(low = 0.01, high = 1, size = nMultilocus) for i in range(2)]
    return(_)

def M(): # Generate Migration rate vector ( M12 and M21 ) , migration can be symetrical or asymetrical
    M12= list(uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus))
    if barrier_simmetry == False:
        M21 =  list(uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus))
    else:
        M21 = M12[:]
    return([M12,M21])

def Tdemo(): # Generate random value vector for Demographic times  
    Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)
    Tdem1 = [ uniform(low = 0, high = Tsplit[i], size = 1)[0] for i in range(nMultilocus) ]
    Tdem2 = [ uniform(low = 0, high = Tsplit[i], size = 1)[0] for i in range(nMultilocus) ]
    Time_list = [Tsplit,Tdem1,Tdem2]
    if 'SC' in model:
        Tsc = [ uniform(low = min_Tsc, high = max_Tsc * Tsplit[i], size = 1)[0] for i in range(nMultilocus) ]
        Time_list.append(Tsc)
    if 'AM' in model:
        Tam = [ uniform(low = min_Tam*Tsplit[i], high = Tsplit[i], size = 1)[0] for i in range(nMultilocus) ]
        Time_list.append(Tam)
    return(Time_list) # the param order is Tsplit, Tdem , Tdem2 , Tsc or Tam


def shape_beta(): # In case of modeBarrier = beta, generate randomly a and b values vector of the beta distribution
    shape_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
    shape_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
    return([shape_a, shape_b])

def beta_dis(X,a,b): # In case of modeBarrier == beta, apply the a and b values to a vector X of values
    scalar = beta(a,b,size = nLoci)
    rescale = a / (a +b)
    _ = np.array([ X * i/rescale for i in scalar ])
    return(_)

def barriers(): #  generate a vector of nBarriers values, randomly or already determined by gene_flow_barrier_proportion value
    _ =  np.random.randint(int(0.9*nLoci-1),size=nMultilocus) 
    return(_)

def produceBarriers(nLoci, nBarriers): 
    # produces a vector of 0 (non barrier) or 1 (barrier), of size equal to the number of loci
    barriers = [0]*nBarriers + [1]*(nLoci-nBarriers)
    shuffle(barriers)
    return(barriers)



###################### build a DataFrame containing the necessary global parameter of the model

param_name = ['N1','N2','Na','founders1','founders2','Tsplit','Tdem1','Tdem2']
param_values = [N(),founders(),Tdemo()]

if 'SC' in model:
    param_name.append('Tsc')
elif 'AM' in model:
    param_name.append('Tam')

if 'SI' not in model:
    param_name += ['M12','M21']
    param_values += [M()]

if '2N' in model:
    param_name += ['shape_N_a','shape_N_b']
    param_values += [shape_beta()]

if '2M' in model and modeBarrier=='beta':
    param_name += ['shape_M12_a','shape_M12_b','shape_M21_a','shape_M21_b','nBarriersM12','nBarriersM21']
    shape_M12 = shape_beta()
    nBarriersM12 = barriers()
    if barrier_simmetry == True:
        shape_M21 = shape_M12
        nBarriersM21 = nBarriersM12
    else:
        shape_M21 = shape_beta()
        nBarriersM21 = barriers()
    param_values += [shape_M12,shape_M21,[nBarriersM12,nBarriersM21]]

if '2M' in model and modeBarrier=='bimodal':
    param_name += ['nBarriersM12','nBarriersM21']
    nBarriersM12 = barriers()
    if barrier_simmetry == True:
        nBarriersM21 = nBarriersM12
    else:
        nBarriersM21 = barriers()
    param_values += [nBarriersM12,nBarriersM21]

param_values = flatter(param_values)

if len(param_values) != len(param_name) :
    with open('error_priorgen.txt','w') as o:
        o.write('the number of param_values is different from the number of param'+'\n')
        o.write(str(len(param_values)) + '\n')
        o.write(str(len(param_name))+ '\n')
    sys.exit()

param_dic = dict(zip(param_name,param_values)) # build the dictionnary of simulation parameters
param_df = DataFrame(param_dic) # Create global parameters DataFrame


########################### Generate locus level parameters

## define the ms argument order depending of the model
base_prior = ['theta','rho','L','nsamA','nsamB'] # this part is common to all models
if 'IM' in model:
    locus_prior_names = base_prior + ['N1_vec','N2_vec','Tdem1','founders1xNa','Tdem2','founders2xNa','M12_vec','M21_vec','Tsplit','Tsplit','Na_vec']
elif 'SC' in model:
    locus_prior_names = base_prior + ['M12_vec','M21_vec','N1_vec','N2_vec','Tdem1','founders1xNa','Tdem2','founders2xNa','Tsc','Tsplit','Tsplit','Na_vec']
elif 'AM' in model:
    locus_prior_names = base_prior + ['N1_vec','N2_vec','Tdem1','founders1xNa','Tdem2','founders2xNa','Tam','M12_vec','M21_vec','Tsplit','Tsplit','Na_vec']
elif 'SI' in model:
    locus_prior_names = base_prior + ['N1_vec','N2_vec','Tdem1','founders1xNa','Tdem2','founders2xNa','Tsplit','Tsplit','Na_vec']

max_param_loop = 3 # the max number of time of calling entirely a generator
sim_param = {a:it.chain(*[(x for x in param_df[a]) for i in range(max_param_loop)]) for a in param_df.columns}# generate a iterator for each column and store it inside the list

df_list = [DataFrame() for i in range(nMultilocus)] # generate a list of dataframe for each nMulitlocus simulation. Each rows is locus specifiq parameters

## Generate N individuals vectors depending of the model
if 'shape_N_a' in param_df.columns:
    for x in  ['Na','N1','N2']:
        for he_df in df_list : 
            he_df[x+'_vec'] = beta_dis(X=next(sim_param[x]),a=next(sim_param['shape_N_a']),b=next(sim_param['shape_N_b'])) 
else : 
    for x in  ['Na','N1','N2']:
        for he_df in df_list : 
            he_df[x+'_vec'] =repeat(next(sim_param[x]),nLoci)

## Generate Migration vector depending of the model
if 'shape_M12_a' in param_df.columns and 'nBarriersM12' in param_df.columns :
    for x in ['M12','M21']:
        for he_df in df_list : 
            he_df[x+'_vec'] = beta_dis(X=next(sim_param[x]),a=next(sim_param[('shape_'+x+'_a')]),b=next(sim_param[('shape_'+x+'_b')])) * produceBarriers(nLoci=nLoci,nBarriers=next(sim_param[('nBarriers'+x)])) 
elif 'nBarriersM12' in param_df.columns: 
    for x in ['M12','M21']:
        for he_df in df_list : 
            he_df[x+'_vec'] =produceBarriers(nLoci=nLoci,nBarriers=next(sim_param[('nBarriers'+x)])) * next(sim_param[x])
elif 'SI' not in model :
    for x in ['M12','M21']:
        for he_df in df_list : 
             he_df[x+'_vec']=repeat(next(sim_param[x]),nLoci)

if barrier_simmetry == True:
    for he_df in df_list:
        he_df['M21_vec'] = he_df['M12_vec']

## Generate Na x founders vectors:
for x in ['founders1','founders2']:
    for he_df in df_list :
        he_df[x+'xNa'] = he_df['Na_vec'] * next(sim_param[x])
## Add the rest of the parameters in the he_df
for x in [p for p in  locus_prior_names if p not in list(df_list[0].columns)+ list(locus_df.columns)]:
    for he_df in df_list :
        he_df[x] = repeat(next(sim_param[x]),nLoci)
##concat the basis iformation with he_df
param_dflist = []
for he_df  in df_list:
    tmp=he_df[['Na_vec','N1_vec','N2_vec','founders1xNa','founders2xNa']].copy()
    tmp[tmp<1e-4]=1e-4
    he_df.update(tmp)
    param_dflist.append(concat([he_df,locus_df],axis=1)) # create the complete dataframe for each sim

################### write the priorfiles and pipe the ms argument through print()
# write priorfile.txt which contains global simulation parameters
if global_write == True:
    with open('priorfile.txt','w') as o:
        o.write(param_df.to_csv(sep="\t",header=True,index_label='dataset',float_format='%.5f'))

if locus_write == True:
    with open('priorfile_locus.txt','w') as lo:
        lo.write(param_dflist[0].to_csv(sep="\t",header=True,index_label='dataset',float_format='%.5f')) # only the header of the first df is use
        for df in param_dflist[1:] :
            lo.write(df.to_csv(sep="\t",header=False,index_label='dataset',float_format='%.5f'))

num_dataset=0
with open('exec.sh','w') as o:
    for sim in param_dflist:
        sim.round(5)
        output_command = '{ '
        for locus in range(nLoci):
            locus_prior_str = [str(short(i)) for i in sim.loc[locus,locus_prior_names]]
            mscommand_formated =  mscommand.format(*locus_prior_str) 
            output_command = output_command +  first_part_command + mscommand_formated + " ;"
        output_command += '}' + ' | pypy {binpath}/mscalc_new.py bpfile={bpfile} num_dataset={num_dataset} \n'.format(binpath=binpath,bpfile=bpfile,num_dataset=num_dataset) 
        num_dataset += 1
        o.write(output_command)
