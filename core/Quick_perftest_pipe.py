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
MODELS_COMP = ['SC_1M_1N', 'SC_1M_2N', 'SC_2M_1N', 'SC_2M_2N', 'AM_1M_1N', 'AM_1M_2N', 'AM_2M_1N', 'AM_2M_2N', 'IM_1M_1N', 'IM_1M_2N', 'IM_2M_1N', 'IM_2M_2N', 'SI_1N', 'SI_2N']
MODELS_COMP = ['SC_1M_1N', 'SC_1M_2N', 'SC_3M_1N', 'SC_3M_2N', 'AM_1M_1N', 'AM_1M_2N', 'AM_3M_1N', 'AM_3M_2N', 'IM_1M_1N', 'IM_1M_2N', 'IM_3M_1N', 'IM_3M_2N', 'SI_1N', 'SI_2N']

# first and second estimation of parameters
nIterations_gof = 10 # number of subdirectories for the simulations used in the nnet param estimates (10)
ITERATIONS_ESTIMATES = range(nIterations_estim)


# informations from the config.yaml file
nameA = config['nameA']
nameB = config['nameB']
useSFS = config['useSFS']
timeStamp = config['timeStamp']
wdir = timeStamp + '/simdata/' + config['simdata_dir']
nMin = config['nMin']
Nref = (0 + config['N_max'])/2.0 # Nref is the mid point of the prior
window_size = config['window_size']
config_yaml = timeStamp + '/'  +config['config_yaml']
popfile= timeStamp + '/'  +  config['popfile']
############# singularity parametrisaiton #########
container_path =binpath + '/container' 
Sc='singularity exec --bind {0},{1} {2}'.format(binpath,timeStamp,container_path)
def_container='python3-pd.sif'
############### Get used model ################
import os
model_filter_file='{wdir}/model_filter_2.txt'.format(wdir=wdir)
if os.path.isfile(model_filter_file):
    with open(model_filter_file,'r') as f:
        used_model = f.readline().strip().replace('"','').split('\t')
else :
    used_model=MODELS_COMP
############### locus model ################
locus_model=['IM_2M_2N'] # list of model used to test the barrier
# IM test for current and ancestral barriers
# SC test for recent current barriers
# AM test for ancestral barriers
locus_param2estimate = ['M_current']
############### wildcards management ############## 
# Prevent Snakemake to use regex behavior on wildcards
wildcard_constraints:
    wdir=wdir,
    timeStamp=timeStamp,
    binpath=binpath,
    model='|'.join(MODELS_COMP),
    used_model = '|'.join(used_model),
    locus_model = '|'.join(locus_model),
    param_locus = '|'.join(locus_param2estimate),
    ii='|'.join([str(x) for x in ITERATIONS_ESTIMATES]),
    j='|'.join([str(x) for x in range(100)])

############## End of Pipeline & Targets ############
rule targets: # edit at the end 
    input:
        bpfile = expand("{wdir}/locus_datafile", wdir=wdir),
        ABCstat_global = expand("{wdir}/ABCstat_global.txt",wdir=wdir),
        ABCstat_locus = expand("{wdir}/ABCstat_locus.txt",wdir=wdir),
        sim =expand("{timeStamp}/modelComp/{model}_{i}/ABCstat_global.txt",timeStamp=timeStamp,model=MODELS_COMP,i=ITERATIONS_MODEL_COMP),
        model_filter =  expand("{wdir}/models_filter.txt",wdir=wdir),
        param_estimation_1 = expand("{wdir}/estimation_param_1/posterior_{used_model}.txt",wdir=wdir,used_model=used_model)
       # param_simulation = expand("{wdir}/estimation_param_{n_est}/{used_model}_{ii}/ABCstat_global.txt",wdir=wdir,n_est=n_estimation,used_model=used_model,ii=ITERATIONS_ESTIMATES),
       # param_estimation = expand("{wdir}/estimation_param_{n_est}/posterior_{used_model}.txt",wdir=wdir,n_est=n_estimation,used_model=used_model),
       # model_weight= expand("{wdir}/models_weight.txt",wdir=wdir),
       # mw_post = expand('{wdir}/locus_posteriors_mw.txt',wdir=wdir),
       # locus_estimation = expand('{wdir}/locus_param_estimation/{locus_model}_{param_locus}_estimate_1Pd.txt',wdir=wdir,
       # locus_model=locus_model,param_locus=locus_param2estimate)
#        table = expand('{wdir}/roc_curve/roc_table.txt',wdir=wdir),
#        fig = expand('{wdir}/roc_curve/roc.pdf',wdir=wdir)
    shell:
            """
#            rm -r {wdir}/locus_sim_param
            echo done
            """
############################ generating global data ############

rule simulationsModelComp:
    params:
        nmultilocus={nmultilocus}
    input:
        locus_datafile = expand("{wdir}/locus_datafile",wdir=wdir)
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



############## Model filtering  ######################
rule model_filtering:
    input:
        obs=expand("{wdir}/ABCstat_global.txt",wdir=wdir),
        sim =expand("{timeStamp}/modelComp/{model}_{i}/ABCstat_global.txt",timeStamp=timeStamp,model=MODELS_COMP,i=ITERATIONS_MODEL_COMP)
    output:
        expand("{wdir}/models_filter.txt",wdir=wdir)
    threads: nCPU_R 
    shell:
        """
        {Sc}/R.sif Rscript {core_path}/model_averaging.R timeStamp={timeStamp} ncores={nCPU_R}\
        ntree={ntree} output_name=models_filter.txt obs_dir={timeStamp}/simdata/ \
        sim_dir={timeStamp}/modelComp obs_pattern=ABCstat_global.txt sim_pattern=ABCstat_global.txt  useSFS={useSFS} mode=multi
        """
############################### Global Parameter estimation  #####  
# First it list the X pertinent models written in models_filter files
# For each model it will proceed to five round of estimation/posterior optimisation
# At the end it will estimate the global parameters of each model and do the ponderated mean on all models values


rule estimation_param_model_1:
    params:
        iteration=1
    input:
        obs= expand("{wdir}/ABCstat_global.txt",wdir=wdir),
        sim =expand("{timeStamp}/modelComp/{model}_{i}/ABCstat_global.txt",timeStamp=timeStamp,model=MODELS_COMP,i=ITERATIONS_MODEL_COMP),
        model_filter=expand("{wdir}/models_filter.txt",wdir=wdir)
    output:
        "{wdir}/estimation_param_1/posterior_{used_model}.txt"
    shell:
        """
        {Sc}/R.sif Rscript {core_path}/point_estimate.R Nref={Nref} nameA={nameA} nameB={nameB} nMin={nMin} \
                sub_dir_sim=modelComp nSubdir={nIterations_model_comp} ntree={ntree} ncores={nCPU_R}\
                useSFS={useSFS} bestModel={wildcards.used_model} timeStamp={timeStamp} \
                nPosterior={nPosterior_locus} binpath={core_path} \
                path2observation={wdir}/ \
                output_dir={wdir}/estimation_param_{params.iteration}
        """
rule param_sim_2 :
    input:
        posteriors = "{wdir}/estimation_param_1/posterior_{used_model}.txt"
    params:
        iteration=2
    output:
        "{wdir}/estimation_param_2/{used_model}_{ii}/priorfile.txt",
        "{wdir}/estimation_param_2/{used_model}_{ii}/ABCstat_global.txt"
    threads: 4
    resources:
        mem_mb=10000
    shell: 
        """
        mkdir -p {wdir}/estimation_param_{params.iteration}/{wildcards.used_model}_{wildcards.ii}
        cd {wdir}/estimation_param_{params.iteration}/{wildcards.used_model}_{wildcards.ii}
        {Sc}/python.sif python3 {core_path}/submit_priorgen_postgen.py model={wildcards.used_model} \
            locus_datafile={wdir}/locus_datafile locus_write=False global_write=True \
            priorfile={input.posteriors} binpath={core_path} nMultilocus={nmultilocus} 
        split -l {split_size} exec.sh sub_
        chmod a+x sub_*
        list=(sub_*)
        for ll in ${{list[@]}}; do echo ./$ll >> tmp_exec.sh; done
        {Sc}/scrm_py.sif parallel -a tmp_exec.sh -j 4 --halt now,fail=1
        rm exec.sh tmp_exec.sh sub* 
        sed -i '2,${{/dataset/d}}' ABC*.txt
        """

rule estimation_param_model_2 :
    input:
        obs= expand("{wdir}/ABCstat_global.txt",wdir=wdir),
        sim =expand("{{wdir}}/estimation_param_2/{{used_model}}_{ii}/ABCstat_global.txt",ii=ITERATIONS_ESTIMATES),
        model_filter=expand("{wdir}/models_filter.txt",wdir=wdir)
    threads: nCPU_R
    params:
        iteration=2
    output:
        posteriors = "{wdir}/estimation_param_2/posterior_{used_model}.txt"
    shell:
        """
        {Sc}/R.sif Rscript {core_path}/estimate_posteriors.R Nref={Nref} nameA={nameA} nameB={nameB} nMin={nMin} \
                sub_dir_sim=estimation_param_{params.iteration} nSubdir={nIterations_estim} ntree={ntree} ncores={nCPU_R}\
                useSFS={useSFS} bestModel={wildcards.used_model} timeStamp={wdir} \
                nPosterior={nPosterior_locus} binpath={core_path} \
                path2observation={wdir}/ \
                output_dir={wdir}/estimation_param_{params.iteration}
        """


use rule param_sim_2 as param_sim_3 with :
    input:
        posteriors = "{wdir}/estimation_param_2/posterior_{used_model}.txt"
    params:
        iteration=3
    output:
        "{wdir}/estimation_param_3/{used_model}_{ii}/priorfile.txt",
        "{wdir}/estimation_param_3/{used_model}_{ii}/ABCstat_global.txt"

use rule estimation_param_model_2 as estimation_param_model_3 with : 
    input:
        sim =expand("{{wdir}}/estimation_param_3/{{used_model}}_{ii}/ABCstat_global.txt",ii=ITERATIONS_ESTIMATES)
    params:
        iteration=3
    output:
        posteriors = "{wdir}/estimation_param_3/posterior_{used_model}.txt"

use rule param_sim_2 as param_sim_4 with:
    input:
        posteriors = "{wdir}/estimation_param_3/posterior_{used_model}.txt"
    params:
        iteration=4
    output:
        "{wdir}/estimation_param_4/{used_model}_{ii}/priorfile.txt",
        "{wdir}/estimation_param_4/{used_model}_{ii}/ABCstat_global.txt"

use rule estimation_param_model_2 as estimation_param_model_4 with:
    input:
        sim =expand("{{wdir}}/estimation_param_4/{{used_model}}_{ii}/ABCstat_global.txt",ii=ITERATIONS_ESTIMATES)
    params:
        iteration=4
    output:
        posteriors = "{wdir}/estimation_param_4/posterior_{used_model}.txt"

use rule param_sim_2 as param_sim_5 with:
    input:
        posteriors = "{wdir}/estimation_param_4/posterior_{used_model}.txt"
    params:
        iteration=5
    output:
        "{wdir}/estimation_param_5/{used_model}_{ii}/priorfile.txt",
        "{wdir}/estimation_param_5/{used_model}_{ii}/ABCstat_global.txt"

use rule estimation_param_model_1 as estimation_param_model_5 with: 
    input:
        sim =expand("{{wdir}}/estimation_param_5/{{used_model}}_{ii}/ABCstat_global.txt",ii=ITERATIONS_ESTIMATES)
    params:
        iteration=5
    output:
        posteriors = "{wdir}/estimation_param_5/posterior_{used_model}.txt"

############## Model weighting  ######################


rule model_averaging: 
    input:
        obs="{wdir}/ABCstat_global.txt",
        sim =expand("{{wdir}}/estimation_param_{n_est}/{used_model}_{ii}/ABCstat_global.txt",used_model=used_model,ii=ITERATIONS_ESTIMATES,n_est=n_estimation),
        posteriors=expand("{{wdir}}/estimation_param_{n_est}/posterior_{used_model}.txt",used_model=used_model,n_est=n_estimation)
    output:
        "{wdir}/models_weight.txt"
    threads: nCPU_R 
    shell:
        """
        {Sc}/R.sif Rscript {core_path}/model_averaging.R timeStamp={wdir} ncores={nCPU_R} \
                ntree={ntree} output_name=models_weight.txt obs_dir={wdir} \
                sim_dir={wdir}/estimation_param_{n_estimation} \
                obs_pattern=ABCstat_global.txt sim_pattern=ABCstat_global.txt  useSFS={useSFS} mode=single
        """
########### Global parameter averaging ############

############### posterior averaging ###############

rule posterior_averaging:
    input:
        model_weight = "{wdir}/models_weight.txt",
        posteriors=expand('{{wdir}}/estimation_param_{n_estimation}/posterior_{used_model}.txt',
            used_model=used_model,n_estimation=n_estimation)
    output:
         '{wdir}/locus_posteriors_mw.txt'
    shell:
        """
        {Sc}/R.sif Rscript {core_path}/posterior_averaging.R \
                master_post_dir={wdir}/estimation_param_{n_estimation} \
                nPosteriors_locus={nPosterior_locus} weight_file={input.model_weight} output={output} 
        """

########################### Generate locus data ####################

rule simulation_locus: 
    input:
        post = "{wdir}/locus_posteriors_mw.txt",
        locus_datafile = "{wdir}/locus_datafile_locussp"
    output:
        "{wdir}/locus_sim_param/{locus_model}/ABCstat_locus.txt",
        "{wdir}/locus_sim_param/{locus_model}/priorfile_locus.txt"
    threads: 50
    shell:
        """
        mkdir -p {wdir}/locus_sim_param/{wildcards.locus_model}
        cd {wdir}/locus_sim_param/{wildcards.locus_model}
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

rule estimate_locus_param:
    input:
        "{wdir}/locus_sim_param/{locus_model}/ABCstat_locus.txt",
        "{wdir}/locus_sim_param/{locus_model}/priorfile_locus.txt"
    output:
        '{wdir}/locus_param_estimation/{locus_model}_{param_locus}_estimate_1Pd.txt'
    threads: nCPU_R
    shell:
        """
        mkdir -p {wdir}/locus_param_estimation/
        {Sc}/R.sif Rscript {core_path}/estimate_locus_param_rf.R timeStamp={wdir} ncores={nCPU_R}\
                ntree={ntree} obs_dir={wdir}/ \
                sim_dir={wdir}/locus_sim_param/{wildcards.locus_model} \
                output={output} model={wildcards.locus_model} param={wildcards.param_locus}  
        """

