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
        gof_prior = expand("{timeStamp}/gof_prior.txt",timeStamp=timeStamp),
        posterior = expand("{timeStamp}/posterior.txt",timeStamp=timeStamp),
        mw = expand("{timeStamp}/model_weight.txt",timeStamp=timeStamp),
        sim_post_glob = expand('{timeStamp}/sim_posterior/ABCstat_global.txt',timeStamp=timeStamp),
        gof_posterior = expand("{timeStamp}/gof_posterior.txt",timeStamp=timeStamp),
        sim_locus = expand('{timeStamp}/sim_locus/ABCstat_locus.txt',timeStamp=timeStamp),
        barrier_assignation = expand('{timeStamp}/Pbarrier.txt',timeStamp=timeStamp),
        report = expand('{timeStamp}/report_barrier_detection.txt',timeStamp=timeStamp)
    shell:
            """
            echo 'done'
            """
#### generate subsample #######

rule generate_bed_file:
    input:
        contig_file = contig_file
    params:
        nloci_per_chr_global=100,
        nloci_per_chr_full_locus=-1 # here -1 mean there is no sampling, all locus are analysed. 
    output:
        bed_global = '{timeStamp}/bed_global_dataset.txt',
        bed_full_locus = '{timeStamp}/bed_full_locus_dataset.txt'
    shell:
        """
            {Sc}/R.sif Rscript {core_path}/generate_bed_sample.R contig_file={contig_file}\
                    window_size={window_size} nLoci_per_chr={params.nloci_per_chr_global}\
                    output={output.bed_global}
            {Sc}/R.sif Rscript {core_path}/generate_bed_sample.R contig_file={contig_file}\
                    window_size={window_size} nLoci_per_chr={params.nloci_per_chr_full_locus}\
                    output={output.bed_full_locus}
        """
rule get_rec_rate:
    input:
        rec_rate_map = rec_rate_map,
        contig_file = contig_file
    output:
        '{timeStamp}/bed_rec_rate.txt'
    shell:
        """
            {Sc}/R.sif Rscript {core_path}/get_rec_rate.R contig_file={contig_file}\
                    window_size={window_size} rho_map={rec_rate_map} output={output}
        """

rule generate_locus_datafile:
    input:
        bed_rec_rate = '{timeStamp}/bed_rec_rate.txt',
        contig_file = contig_file,
        popfile = popfile,
        bed_global = '{timeStamp}/bed_global_dataset.txt'
    output:
        datafile = '{timeStamp}/locus_datafile'
    shell:
        """
            {Sc}/R.sif Rscript {core_path}/generate_locus_datafile.R \
                    rho_map={input.bed_rec_rate} bedfile={input.bed_global} mu={mu}\
                    Nref={Nref} window_size={window_size} nameA={nameA} nameB={nameB}\
                    popfile={popfile} output={output.datafile} ploidy={ploidy}
        """

rule generate_abc_stat: 
    input:
        vcf_file = expand('{timeStamp}/haplotyped.vcf',timeStamp=timeStamp),
        popfile = popfile,
        bed_global = '{timeStamp}/bed_global_dataset.txt',
        bed_locus = '{timeStamp}/bed_full_locus_dataset.txt'
    output:
        '{timeStamp}/ABCstat_global.txt',
        '{timeStamp}/ABCstat_locus.txt'
    shell:
        """
            {Sc}/scrm_py.sif python3 {core_path}/vcf2abc.py  data={input.vcf_file} bed_file={input.bed_global}\
                    popfile={popfile} nameA={nameA} nameB={nameB} window_size={window_size}\
                    locus_write="False" global_write="True" output_dir={timeStamp}
            {Sc}/scrm_py.sif python3 {core_path}/vcf2abc.py  data={input.vcf_file} bed_file={input.bed_locus}\
                    popfile={popfile} nameA={nameA} nameB={nameB} window_size={window_size}\
                    locus_write="True" global_write="False" output_dir={timeStamp}
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
        "{timeStamp}/model_weight.txt"
    shell:
        """
        {Sc}/R.sif Rscript {core_path}/estimate_posterior_and_mw.R  \
                obs_dir={timeStamp} sim_dir={timeStamp}/modelComp ncores={nCPU_R} ntree={ntree} mode=single nPosterior={nPosterior_locus}
        """
### gof #####

rule gof_prior:
    input: 
        obs= expand("{timeStamp}/ABCstat_global.txt",timeStamp=timeStamp),
        sim =expand("{timeStamp}/modelComp/{model}_{i}/ABCstat_global.txt",timeStamp=timeStamp,model=MODELS_COMP,i=ITERATIONS_MODEL_COMP),
    output:
        expand("{timeStamp}/gof_prior.txt",timeStamp=timeStamp)
    shell:
        """
        {Sc}/R.sif Rscript {core_path}/gof_estimate.R obs_dir={timeStamp} \
                sim_dir={timeStamp}/modelComp mode=single output=gof_prior.txt \
                nb_replicate={nPosterior_locus} type=prior
        """

rule sim_posterior :
    input:
        posteriors = "{timeStamp}/posterior.txt",
        gof_prior=expand("{timeStamp}/gof_prior.txt",timeStamp=timeStamp)
    output:
        "{timeStamp}/sim_posterior/priorfile.txt",
        "{timeStamp}/sim_posterior/ABCstat_global.txt",
    threads: 50
    shell: 
        """
        mkdir -p {timeStamp}/sim_posterior/
        cd {timeStamp}/sim_posterior/
        {Sc}/python.sif python3 {core_path}/submit_priorgen_gof.py \
            locus_datafile={timeStamp}/locus_datafile locus_write=False global_write=True \
            priorfile={input.posteriors} binpath={core_path} nMultilocus={nmultilocus} 
        split -l {split_size_locus} exec.sh sub_
        chmod a+x sub_*
        list=(sub_*)
        for ll in ${{list[@]}}; do echo ./$ll >> tmp_exec.sh; done
        {Sc}/scrm_py.sif parallel -a tmp_exec.sh -j 4 --halt now,fail=1
        rm exec.sh tmp_exec.sh sub* 
        sed -i '2,${{/dataset/d}}' ABC*.txt
        """

rule gof_posterior:
    input:
        obs= "{timeStamp}/ABCstat_global.txt",
        posteriors = "{timeStamp}/posterior.txt",
        sim="{timeStamp}/sim_posterior/ABCstat_global.txt"
    output:
        "{timeStamp}/gof_posterior.txt"
    shell:
        """
        {Sc}/R.sif Rscript {core_path}/gof_estimate.R obs_dir={timeStamp} \
                sim_dir={timeStamp}/sim_posterior mode=single \
                output=gof_posterior.txt nb_replicate={nPosterior_locus} type=post
        """
########################### Generate locus data ####################

rule simulation_locus: 
    input:
        post = "{timeStamp}/posterior.txt",
        locus_datafile = "{timeStamp}/locus_datafile"
    output:
        "{timeStamp}/sim_locus/ABCstat_locus.txt",
        "{timeStamp}/sim_locus/priorfile_locus.txt"
    threads: 50
    shell:
        """
        mkdir -p {timeStamp}/sim_locus
        cd {timeStamp}/sim_locus
        {Sc}/python.sif python3 {core_path}/submit_priorgen_locus.py  \
            locus_datafile={input.locus_datafile} \
            config_yaml={config_yaml} binpath={core_path} nMultilocus=1000 priorfile={input.post} locus_write=True global_write=False 
        split exec.sh -l{split_size_locus} sub_exec
        list_sub=(sub_exec*)
        for sub in ${{list_sub[@]}}; do echo sh $sub >> tmp_exec.sh; done
        {Sc}/scrm_py.sif parallel -a tmp_exec.sh -j 50
        rm tmp_exec.sh sub_exec* exec.sh
        sed -i '2,${{/dataset/d}}' ABCstat_locus.txt
        """

######################### Detect barrier locus ########################
rule barrier_detection:
    input:
        "{timeStamp}/sim_locus/ABCstat_locus.txt",
        "{timeStamp}/sim_locus/priorfile_locus.txt"
    output:
        '{timeStamp}/Pbarrier.txt',
        '{timeStamp}/report_barrier_detection.txt'
    threads: 8
    shell:
        """
        {Sc}/R.sif Rscript {core_path}/barrier_detection.R  ncores=8\
                ntree={ntree} obs_dir={timeStamp}/ \
                sim_dir={wildcards.timeStamp}/sim_locus/ mode=normal \
                posterior={timeStamp}/posterior.txt
        """
