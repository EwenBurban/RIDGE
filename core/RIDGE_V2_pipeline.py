# links to the codes
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
    nCPU_R = 8 # number of CPUs for the model comp for the model forest R functions (8)
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


# informations from the config.yaml file
nameA = config['nameA']
nameB = config['nameB']
useSFS = config['useSFS']
timeStamp = config['timeStamp']
nMin = config['nMin']
Nref = config['Nref'] # Nref is the mid point of the prior
window_size = config['window_size']
config_yaml = timeStamp + '/'  +config['config_yaml']
popfile= timeStamp + '/'  +  config['popfile']
mode= config['mode']
############# singularity parametrisaiton #########
container_path =binpath + '/container' 
Sc='singularity exec --bind {0},{1} {2}'.format(binpath,timeStamp,container_path)
def_container='python3-pd.sif'
############### locus model ################
locus_model=['IM_3M_2N'] # list of model used to test the barrier
# IM test for current and ancestral barriers
# SC test for recent current barriers
# AM test for ancestral barriers
locus_param2estimate = ['M_current']
############### wildcards management ############## 
# Prevent Snakemake to use regex behavior on wildcards
wildcard_constraints:
    timeStamp=timeStamp,
    binpath=binpath,
    model='|'.join(MODELS_COMP),
    locus_model = '|'.join(locus_model),
    param_locus = '|'.join(locus_param2estimate),
    j='|'.join([str(x) for x in range(100)])

############## End of Pipeline & Targets ############
rule targets: # edit at the end 
    input:
        bpfile = expand("{timeStamp}/locus_datafile", timeStamp=timeStamp),
        ABCstat_global = expand("{timeStamp}/ABCstat_global.txt",timeStamp=timeStamp),
        ABCstat_locus = expand("{timeStamp}/ABCstat_locus.txt",timeStamp=timeStamp),
        sim =expand("{timeStamp}/modelComp/{model}_{i}/ABCstat_global.txt",timeStamp=timeStamp,model=MODELS_COMP,i=ITERATIONS_MODEL_COMP),
        param_simulation = expand("{timeStamp}/sim_posterior_gof/ABCstat_global.txt",timeStamp=timeStamp),
        param_estimation = expand("{timeStamp}posterior.txt",timeStamp=timeStamp),
        locus_estimation = expand('{timeStamp}/locus_param_estimation/{locus_model}_{param_locus}_estimate_1Pd.txt',timeStamp=timeStamp,
        locus_model=locus_model,param_locus=locus_param2estimate),
        table = expand('{timeStamp}/roc_curve/roc_table.txt',timeStamp=timeStamp),
        fig = expand('{timeStamp}/roc_curve/roc.pdf',timeStamp=timeStamp)
    shell:
            """
            echo 'done'
            """
############################ generating global data ############

rule simulationsModelComp:
    params:
        nmultilocus={nmultilocus}
    input:
        locus_datafile = "{timeStamp}/locus_datafile",
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
        obs= expand("{timeStamp}/ABCstat_global.txt",timeStamp=timeStamp),
        sim =expand("{timeStamp}/modelComp/{model}_{i}/ABCstat_global.txt",timeStamp=timeStamp,model=MODELS_COMP,i=ITERATIONS_MODEL_COMP),
    output:
        "{timeStamp}/posterior.txt",
        "{timeStamp}/model_weight"
    shell:
        """
        {Sc}/R.sif Rscript {core_path}/estimate_posterior_and_mw.R  \
                obs_dir={timeStamp} sim_dir={timeStamp}/modelComp ncores={nCPU_R} ntree={ntree} 
        """
### gof #####


rule gof_prior:
    input: 
        obs= expand("{timeStamp}/ABCstat_global.txt",timeStamp=timeStamp),
        sim =expand("{timeStamp}/modelComp/{model}_{i}/ABCstat_global.txt",timeStamp=timeStamp,model=MODELS_COMP,i=ITERATIONS_MODEL_COMP),
    output:
        "{timeStamp}/gof_prior.txt"
    shell:
        """
        {Sc}/R.sif Rscript obs_dir={timeStamp}/ sim_dir={timeStamp}/modelComp model=single output={output}
        """

rule sim_posterior_gof :
    input:
        posteriors = "{timeStamp}/posterior.txt"
    params:
        iteration=2
    output:
        "{timeStamp}/sim_posterior/priorfile.txt",
        "{timeStamp}/sim_posterior/ABCstat_global.txt"
    threads: 4
    resources:
        mem_mb=10000
    shell: 
        """
        mkdir -p {timeStamp}/sim_posterior_gof/
        cd {timeStamp}/sim_posterior_gof/
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
        obs= expand("{timeStamp}/ABCstat_global.txt",timeStamp=timeStamp),
        sim=expand("{timeStamp}/sim_posterior/ABCstat_global.txt",timeStamp=timeStamp)
    output:
        "{timeStamp}/gof_posterior.txt"
    shell:
        """
        {Sc}/R.sif Rscript obs_dir={timeStamp}/ sim_dir={timeStamp}/modelComp model=single output={output}
        """

########################### Generate locus data ####################

rule simulation_locus_roc: 
    input:
        post = "{timeStamp}/locus_posteriors_mw.txt",
        locus_datafile = "{timeStamp}/locus_datafile_locussp"
    output:
        "{timeStamp}/locus_sim_param/{locus_model}_roc/ABCstat_locus.txt",
        "{timeStamp}/locus_sim_param/{locus_model}_roc/priorfile_locus.txt"
    threads: 50
    shell:
        """
        mkdir -p {timeStamp}/locus_sim_param/{wildcards.locus_model}_roc
        cd {timeStamp}/locus_sim_param/{wildcards.locus_model}_roc
        {Sc}/python.sif python3 {core_path}/submit_priorgen_locus.py model={wildcards.locus_model} \
            locus_datafile={input.locus_datafile} \
            config_yaml={config_yaml} binpath={core_path} nMultilocus={nmultilocus} priorfile={input.post} locus_write=True global_write=True 
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
        locus_estimation = expand('{timeStamp}/locus_param_estimation/{locus_model}_{param_locus}_estimate_1Pd.txt',timeStamp=timeStamp,locus_model=locus_model,param_locus=locus_param2estimate),
        abc = expand("{timeStamp}/locus_sim_param/{locus_model}_roc/ABCstat_locus.txt",timeStamp=timeStamp,locus_model=locus_model),
        prior = expand("{timeStamp}/locus_sim_param/{locus_model}_roc/priorfile_locus.txt",timeStamp=timeStamp,locus_model=locus_model)
    params:
        locus_model=locus_model
    output:
        table = '{timeStamp}/roc_curve/roc_table.txt',
        fig = '{timeStamp}/roc_curve/roc.pdf'
    shell:
        """
            {Sc}/R.sif Rscript {core_path}/roc_estimation.R ntree=1000 ncores=8 timeStamp={timeStamp} sim_dir=locus_sim_param/{params.locus_model}_roc\
                    output_dir={timeStamp}/roc_curve mig_mod=M_current locus_estimation={input.locus_estimation}
        """
