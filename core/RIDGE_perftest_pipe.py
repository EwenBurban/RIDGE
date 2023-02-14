# links to the codes
import os
binpath = config['binpath']
core_path = binpath + '/core'
# light or normal mode
lightMode = config['lightMode']

# general property
nmultilocus = 1000 # number of multilocus simulations per iteration (500)
nPosterior_locus = 1000
split_size=int(nmultilocus/4)
split_size_locus=int(nmultilocus/50)
# data generation
# model comparison
if lightMode==False:
    nCPU_R = 10 # number of CPUs for the model comp for the model forest R functions (8)
    ntree = 1000 # number of tree for the random forest (RF) model comparison (1000)
    nIterations_model_comp = 10 # number of subdirectories for the simulations used in the RF model comparison
    nIterations_estim = 8 # number of subdirectories for the simulations used in the nnet param estimates (500)
    n_estimation=5
else:
    nCPU_R = 2 # number of CPUs for the model comp for the model forest R functions (8)
    ntree = 1000 # number of tree for the random forest (RF) model comparison (1000)
    nIterations_model_comp = 5 # number of subdirectories for the simulations used in the RF model comparison
    nIterations_estim = 4 # number of subdirectories for the simulations used in the nnet param estimates (500)
    n_estimation=2
Nthread_shape_it=8
ITERATIONS_MODEL_COMP = range(nIterations_model_comp)
MODELS_COMP = ['SC_1M_1N', 'SC_1M_2N', 'SC_3M_1N', 'SC_3M_2N', 'AM_1M_1N', 'AM_1M_2N', 'AM_3M_1N', 'AM_3M_2N', 'IM_1M_1N', 'IM_1M_2N', 'IM_3M_1N', 'IM_3M_2N', 'SI_1N', 'SI_2N']

# first and second estimation of parameters
nIterations_gof = 10 # number of subdirectories for the simulations used in the nnet param estimates (10)
ITERATIONS_ESTIMATES = range(nIterations_estim)


# informations from the config.yaml file
nameA = config['nameA']
nameB = config['nameB']
useSFS = config['useSFS']
Pbarrier_max=config['Pbarrier_max']
timeStamp = config['timeStamp']
list_dir=os.listdir(timeStamp+'/simdata')
wdir=[timeStamp+'/simdata/'+x for x in list_dir]
nMin = config['nMin']
Nref = (0 + config['N_max'])/2.0 # Nref is the mid point of the prior
window_size = config['window_size']
config_yaml = timeStamp + '/'  +config['config_yaml']
popfile= timeStamp + '/'  +  config['popfile']
print(wdir)
############# singularity parametrisaiton #########
container_path =binpath + '/container' 
Sc='singularity exec --bind {0},{1} {2}'.format(binpath,timeStamp,container_path)
def_container='python3-pd.sif'
############### locus model ################
locus_model=['IM_3M_2N'] # list of model used to test the barrier
# IM test for current and ancestral barriers
# SC test for recent current barriers
# AM test for ancestral barriers
############### wildcards management ############## 
# Prevent Snakemake to use regex behavior on wildcards
wildcard_constraints:
    wdir='|'.join(wdir),
    timeStamp=timeStamp,
    binpath=binpath,
    model='|'.join(MODELS_COMP),
    locus_model = '|'.join(locus_model),
    ii='|'.join([str(x) for x in ITERATIONS_ESTIMATES]),
    j='|'.join([str(x) for x in range(100)])

############## End of Pipeline & Targets ############
rule targets: # edit at the end 
    input:
        ABCstat_global = expand("{wdir}/ABCstat_global.txt",wdir=wdir),
        ABCstat_locus = expand("{wdir}/ABCstat_locus.txt",wdir=wdir),
        sim =expand("{timeStamp}/modelComp/{model}_{i}/ABCstat_global.txt",timeStamp=timeStamp,model=MODELS_COMP,i=ITERATIONS_MODEL_COMP),
##      gof_prior = expand("{wdir}/gof_prior.txt",wdir=wdir),
##      posterior = expand("{wdir}/posterior.txt",wdir=wdir),
##      mw = expand("{wdir}/model_weight.txt",wdir=wdir),
##      sim_post_glob = expand('{wdir}/sim_posterior/ABCstat_global.txt',wdir=wdir),
##      sim_post_loc = expand('{wdir}/sim_posterior/ABCstat_locus.txt',wdir=wdir),
##      gof_posterior = expand("{wdir}/gof_posterior.txt",wdir=wdir),
##      barrier_assignation = expand('{wdir}/Pbarrier.txt',wdir=wdir),
#        table = expand('{wdir}/roc_table.txt',wdir=wdir),
#        fig = expand('{wdir}/roc.pdf',wdir=wdir)
    shell:
            """
            echo done
            """
############################ generating global data ############
rule gen:
    output:
        "{wdir}/locus_datafile",
        "{wdir}/locus_datafile_locussp"
    shell: 
        """
        cp {timeStamp}/locus_datafile* {wildcards.wdir}/
        """




rule simulationsModelComp:
    params:
        nmultilocus={nmultilocus}
    input:
        locus_datafile = expand("{timeStamp}/locus_datafile",timeStamp=timeStamp)
    output:
        "{timeStamp}/modelComp/{model}_{i}/priorfile.txt",
        "{timeStamp}/modelComp/{model}_{i}/ABCstat_global.txt",
    threads: 4
    resources:
        mem_mb=10000
    shell:
        """
        cd {timeStamp}/modelComp/{wildcards.model}_{wildcards.i}
        {Sc}/python.sif python3 {core_path}/submit_priorgen_newsim.py model={wildcards.model} \
            locus_datafile={input.locus_datafile} locus_write=False global_write=True \
            config_yaml={config_yaml} binpath={core_path} nMultilocus={nmultilocus} 
        split -l {split_size} exec.sh sub_
        list=(sub_*)
        chmod a+x sub_*
        for ll in ${{list[@]}}; do echo ./$ll >> tmp_exec.sh; done
        {Sc}/scrm_py.sif parallel -a tmp_exec.sh -j 4 --halt now,fail=1
        rm exec.sh sub_* tmp_exec.sh
        sed -i '2,${{/dataset/d}}' ABC*.txt
        """




############################### Global Parameter estimation  #####  
# First it list the X pertinent models written in models_filter files
# For each model it will proceed to five round of estimation/posterior optimisation
# At the end it will estimate the global parameters of each model and do the ponderated mean on all models values


rule estimation_posterior_and_model_weight:
    input:
        obs= expand("{timeStamp}/simdata/{list_dir}/ABCstat_global.txt",timeStamp=timeStamp,list_dir=list_dir),
        sim =expand("{timeStamp}/modelComp/{model}_{i}/ABCstat_global.txt",timeStamp=timeStamp,model=MODELS_COMP,i=ITERATIONS_MODEL_COMP),
    output:
        expand("{timeStamp}/{list_dir}/posterior.txt",timeStamp=timeStamp,list_dir=list_dir),
        expand("{timeStamp}/{list_dir}/model_weight.txt",timeStamp=timeStamp,list_dir=list_dir)
    shell:
        """
        {Sc}/R.sif Rscript {core_path}/estimate_posterior_and_mw.R  \
                obs_dir={timeStamp}/simdata sim_dir={timeStamp}/modelComp ncores={nCPU_R} ntree={ntree} nPosterior={nPosterior_locus} mode=multi \
                Pbarrier_max={Pbarrier_max}
        """


### gof #####


rule gof_prior:
    input: 
        obs= expand("{timeStamp}/simdata/{list_dir}/ABCstat_global.txt",timeStamp=timeStamp,list_dir=list_dir),
        sim =expand("{timeStamp}/modelComp/{model}_{i}/ABCstat_global.txt",timeStamp=timeStamp,model=MODELS_COMP,i=ITERATIONS_MODEL_COMP),
    output:
        expand("{timeStamp}/simdata/{list_dir}/gof_prior.txt",timeStamp=timeStamp,list_dir=list_dir)
    shell:
        """
        {Sc}/R.sif Rscript obs_dir={timeStamp}/simdata sim_dir={timeStamp}/modelComp mode=multi output=gof_prior.txt nb_replicate={nPosterior_locus}
        """



rule sim_posterior :
    input:
        posteriors = "{wdir}/posterior.txt",
        gof_prior=expand("{timeStamp}/simdata/{list_dir}/gof_prior.txt",timeStamp=timeStamp,list_dir=list_dir)
    output:
        "{wdir}/sim_posterior/priorfile.txt",
        "{wdir}/sim_posterior/priorfile_locus.txt",
        "{wdir}/sim_posterior/ABCstat_global.txt",
        "{wdir}/sim_posterior/ABCstat_locus.txt"
    threads: 4
    resources:
        mem_mb=10000
    shell: 
        """
        mkdir -p {wildcards.wdir}/sim_posterior/
        cd {wildcards.wdir}/sim_posterior/
        {Sc}/python.sif python3 {core_path}/submit_priorgen_gof.py \
            locus_datafile={timeStamp}/locus_datafile locus_write=True global_write=True \
            priorfile={input.posteriors} binpath={core_path} nMultilocus={nmultilocus} 
        split -l {split_size} exec.sh sub_
        chmod a+x sub_*
        list=(sub_*)
        for ll in ${{list[@]}}; do echo ./$ll >> tmp_exec.sh; done
        {Sc}/scrm_py.sif parallel -a tmp_exec.sh -j 4 --halt now,fail=1
        rm exec.sh tmp_exec.sh sub* 
        sed -i '2,${{/dataset/d}}' ABC*.txt
        """

rule gof_posterior:
    input:
        obs= "{wdir}/ABCstat_global.txt",
        posteriors = "{wdir}/posterior.txt",
        sim="{wdir}/sim_posterior/ABCstat_global.txt"
    output:
        "{wdir}/gof_posterior.txt"
    shell:
        """
        {Sc}/R.sif Rscript obs_dir={input.obs} sim_dir={input.sim} mode=single output=gof_posterior.txt nb_replicate={nPosterior_locus}
        """
######################### Detect barrier locus ########################
rule estimate_locus_param:
    input:
        "{wdir}/sim_posterior/ABCstat_locus.txt",
        "{wdir}/sim_posterior/priorfile_locus.txt"
    output:
        '{wdir}/Pbarrier.txt'
    threads: nCPU_R
    shell:
        """
        {Sc}/R.sif Rscript {core_path}/estimate_locus_param_rf.R timeStamp={wildcards.wdir} ncores={nCPU_R}\
                ntree={ntree} obs_dir={wildcards.wdir}/ \
                sim_dir={wildcards.wdir}/sim_posterior/ \
                output={output}  param=M_current 
        """

########################### Calculate theoritical ROC AUC of barrier inference ####################

rule simulation_locus_roc: 
    input:
        pbarrier='{wdir}/Pbarrier.txt',
        post = "{wdir}/posterior.txt",
        locus_datafile = "{wdir}/locus_datafile_locussp"
    output:
        "{wdir}/sim_roc/ABCstat_locus.txt",
        "{wdir}/sim_roc/priorfile_locus.txt"
    threads: 50
    shell:
        """
        mkdir -p {wildcards.wdir}/sim_roc
        cd {wildcards.wdir}/sim_roc
        {Sc}/python.sif python3 {core_path}/submit_priorgen_gof.py  \
            locus_datafile={input.locus_datafile} \
            config_yaml={config_yaml} binpath={core_path} nMultilocus=100 priorfile={input.post} locus_write=True global_write=False 
        split exec.sh -l{split_size_locus} sub_exec
        list_sub=(sub_exec*)
        for sub in ${{list_sub[@]}}; do echo sh $sub >> tmp_exec.sh; done
        {Sc}/scrm_py.sif parallel -a tmp_exec.sh -j 50
        rm tmp_exec.sh sub_exec* exec.sh
        sed -i '2,${{/dataset/d}}' ABCstat_locus.txt
        if [[ -e ABCjsfs_locus.txt ]] ;then sed -i '2,${{/dataset/d}}' ABCjsfs_locus.txt ;fi
        """

rule estimate_roc:
    input:
        locus_estimation = '{wdir}/Pbarrier.txt',
        abc = "{wdir}/sim_roc/ABCstat_locus.txt",
        prior = "{wdir}/sim_roc/priorfile_locus.txt"
    group: 'locus'
    output:
        table = '{wdir}/roc_table.txt',
        fig = '{wdir}/roc.pdf'
    shell:
        """
            {Sc}/R.sif Rscript {core_path}/roc_estimation.R ntree=1000 ncores=8 timeStamp={wildcards.wdir} sim_dir=sim_roc\
                    output_dir={wildcards.wdir} mig_mod=M_current locus_estimation={input.locus_estimation}
        """
