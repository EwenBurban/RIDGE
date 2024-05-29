import os
# links to the codes
binpath = config['binpath']
core_path = binpath + '/core'
# light or normal mode
lightMode = config['lightMode']

# general property
if lightMode==False:
    nmultilocus = 250 # number of multilocus simulations per iteration (500)
    nPosterior_locus = 1000
    split_size=int(nmultilocus/4)
    split_size_locus=int(nmultilocus/50)
    nCPU_R = 8 # number of CPUs for the model comp for the model forest R functions (8)
    ntree = 1000 # number of tree for the random forest (RF) model comparison (1000)
    nIterations_model_comp = 40# number of subdirectories for the simulations used in the RF model comparison
else:
    nmultilocus = 250 # number of multilocus simulations per iteration (500)
    nPosterior_locus = 1000
    split_size=int(nmultilocus/4)
    split_size_locus=int(nmultilocus/50)
    nCPU_R = 8 # number of CPUs for the model comp for the model forest R functions (8)
    ntree = 1000 # number of tree for the random forest (RF) model comparison (1000)
    nIterations_model_comp = 10 # number of subdirectories for the simulations used in the RF model comparison
ITERATIONS_MODEL_COMP = range(nIterations_model_comp)
MODELS_COMP = ['SC_1M_1N', 'SC_1M_2N', 'SC_2M_1N', 'SC_2M_2N', 'AM_1M_1N', 'AM_1M_2N', 'AM_2M_1N', 'AM_2M_2N', 'IM_1M_1N', 'IM_1M_2N', 'IM_2M_1N', 'IM_2M_2N', 'SI_1N', 'SI_2N']


# informations from the config.yaml file
nameA = config['nameA']
nameB = config['nameB']
timeStamp = config['work_dir']
Nref = config['Nref'] # Nref is the mid point of the prior
window_size = config['window_size']
config_yaml = timeStamp + '/'  +config['config_yaml']
popfile= timeStamp + '/'  +  config['popfile']
contig_file=timeStamp + '/' + config['contig_data']
rec_rate_map=timeStamp + '/' + config['rec_rate_map']
vcf_file=timeStamp + '/' + config['vcf_file']
ploidy=config['ploidy']
mu=config['mu']
homo_rec=config['homo_rec']
homo_rec_rate=config['homo_rec_rate']
hetero_theta=config['hetero_theta']
mode=config['mode']
nLoci=config['nLoci']
mask_file=timeStamp + '/' +config['mask_file']
############# singularity parametrisaiton #########
container_path =binpath + '/container' 
Sc='singularity exec --bind {0},{1} {2}'.format(binpath,timeStamp,container_path)
############### wildcards management ############## 
# Prevent Snakemake to use regex behavior on wildcards
wildcard_constraints:
    timeStamp=timeStamp,
    binpath=binpath,
    model='|'.join(MODELS_COMP)
############### sanity check  ######################
if not os.path.isfile(mask_file):
    mask_file='None'
    print('mask file not found, RIDGE continue without considering the mask file')

if not os.path.isfile(contig_file):
    print('contig file not found, RIDGE stop !')
    exit()

if not os.path.isfile(popfile):
    print('popfile not found, RIDGE stop !')
    exit()

if not os.path.isfile(vcf_file):
    print('vcf_file not found, RIDGE stop !')
    exit()

############## End of Pipeline & Targets ############
if mode=='scan' : 
    expected_output=['ABCstat_locus.txt','prior_bound_suggestion.txt']
elif mode=='prior' :
    expected_output=['ABCstat_global.txt','locus_datafile']
elif mode=='test' : 
    expected_output=['ABCstat_global.txt','ABCstat_locus.txt','gof_prior.txt','QC_plot/QC_prior_density.pdf','QC_plot/QC_prior_acp.pdf']
    nIterations_model_comp = 1
    ITERATIONS_MODEL_COMP = range(nIterations_model_comp)
    nmultilocus = 10 
    nPosterior_locus = 100
elif mode=='all_no_vis' : 
    expected_output=['ABCstat_global.txt','ABCstat_locus.txt','gof_prior.txt',
            'gof_posterior.txt','posterior.txt','model_weight.txt',
            'Pbarrier.txt','barrier_proportion_and_ratio.txt']
else :
    expected_output=['ABCstat_global.txt','ABCstat_locus.txt','gof_prior.txt',
            'gof_posterior.txt','posterior.txt','model_weight.txt','QC_plot/QC_prior_density.pdf','QC_plot/QC_prior_acp.pdf',
            'Pbarrier.txt','barrier_proportion_and_ratio.txt','QC_plot/QC_posterior_density.pdf','QC_plot/QC_posterior_acp.pdf','visual_posterior.pdf','visual_model_weight.pdf']

rule targets: # edit at the end 
    input:
        outputs=expand('{timeStamp}/{expected_output}',timeStamp=timeStamp,expected_output=expected_output)
    shell:
            """
            echo 'done'
            """
#### generate subsample #######

rule generate_segmentation:
    input:
        contig_file = contig_file
    output:
        bed_full_locus = '{timeStamp}/genome_segmentation.txt'
    shell:
        """
            {Sc}/R.sif Rscript {core_path}/generate_genome_segmentation.R contig_file={contig_file}\
                    window_size={window_size} nLoci={nLoci}\
                    output={output.bed_full_locus}
        """

rule generate_locus_datafile:
    input:
        contig_file = contig_file,
        popfile = popfile,
        abc_data='{timeStamp}/ABCstat_locus.txt'
    output:
        datafile = '{timeStamp}/locus_datafile',
        bed_global='{timeStamp}/bed_global_dataset.txt'
    shell:
        """
            {Sc}/R.sif Rscript {core_path}/generate_locus_datafile.R \
                    rho_map={rec_rate_map} homo_rec={homo_rec} homo_rec_rate={homo_rec_rate} contig_data={contig_file}\
                    Nref={Nref} window_size={window_size} nameA={nameA} nameB={nameB} nLoci={nLoci}\
                    popfile={popfile} output={output.datafile} ploidy={ploidy} output_bed={output.bed_global} \
                    abc_data={input.abc_data} mu={mu} hetero_theta={hetero_theta}
        """

rule generate_abc_stat_global: 
    input:
        vcf_file = vcf_file,
        popfile = popfile,
        bed_global = '{timeStamp}/bed_global_dataset.txt'
    output:
        '{timeStamp}/ABCstat_global.txt'
    shell:
        """
            {Sc}/scrm_py.sif python3 {core_path}/vcf2abc.py  data={input.vcf_file} bed_file={input.bed_global}\
                    popfile={popfile} nameA={nameA} nameB={nameB} window_size={window_size}\
                    locus_write="False" global_write="True" output_dir={timeStamp} ploidy={ploidy} mask_file={mask_file}
        """

rule generate_abc_stat_locus: 
    input:
        vcf_file = vcf_file,
        popfile = popfile,
        bed_locus = '{timeStamp}/genome_segmentation.txt'
    output:
        '{timeStamp}/ABCstat_locus.txt'
    shell:
        """
            {Sc}/scrm_py.sif python3 {core_path}/vcf2abc.py  data={input.vcf_file} bed_file={input.bed_locus}\
                    popfile={popfile} nameA={nameA} nameB={nameB} window_size={window_size}\
                    locus_write="True" global_write="False" output_dir={timeStamp} ploidy={ploidy} mask_file={mask_file}
        """

rule generate_suggestion:
    input:
        abc_locus='{timeStamp}/ABCstat_locus.txt'
    output:
        '{timeStamp}/prior_bound_suggestion.txt'
    shell:
        """
            {Sc}/R.sif Rscript {core_path}/make_prior_bound_suggestion.R data={input.abc_locus}\
                    output={output} mu={mu} window_size={window_size} nameA={nameA} nameB={nameB}\
                    popfile={popfile} ploidy={ploidy} hetero_theta={hetero_theta}
        """

############################ generating global data ############

rule simulationsModelComp:
    params:
        nmultilocus={nmultilocus}
    input:
        locus_datafile = "{timeStamp}/locus_datafile"
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
        check_token = "{timeStamp}/check_modelComp_done.log"
    output:
        "{timeStamp}/posterior.txt",
        "{timeStamp}/model_weight.txt",
        "{timeStamp}/sim_posterior/ABCstat_global.txt"
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
            config_yaml={config_yaml} binpath={core_path} nMultilocus=1000 priorfile={input.post} locus_write=True global_write=False Nref={Nref}
        split exec.sh -l{split_size_locus} sub_exec
        list_sub=(sub_exec*)
        for sub in ${{list_sub[@]}}; do echo sh $sub >> tmp_exec.sh; done
        {Sc}/scrm_py.sif parallel -a tmp_exec.sh
        rm tmp_exec.sh sub_exec* exec.sh
        sed -i '2,${{/dataset/d}}' ABCstat_locus.txt
        """

######################### Detect barrier locus ########################
rule barrier_detection_current:
    input:
        "{timeStamp}/sim_locus/ABCstat_locus.txt",
        "{timeStamp}/ABCstat_locus.txt",
        "{timeStamp}/sim_locus/priorfile_locus.txt"
    output:
        '{timeStamp}/Pbarrier.txt',
        '{timeStamp}/barrier_proportion_and_ratio.txt'
    threads: 8
    shell:
        """
        {Sc}/R.sif Rscript {core_path}/barrier_detection.R  ncores=8\
                ntree={ntree} obs_dir={timeStamp}/ \
                sim_dir={wildcards.timeStamp}/sim_locus/ mode=normal \
                posterior={timeStamp}/posterior.txt
        """


####################### QC_plot ############################

rule QC_prior_sim:
    input:
        '{timeStamp}/ABCstat_global.txt',
        '{timeStamp}/gof_prior.txt'
    output:
        '{timeStamp}/QC_plot/QC_prior_density.pdf',
        '{timeStamp}/QC_plot/QC_prior_acp.pdf'
    shell:
        """
        {Sc}/R_visual.sif Rscript {core_path}/QC_prior_sim.R\
            dir={timeStamp} output_dir={timeStamp}/QC_plot
        """


rule QC_posterior:
    input:
        '{timeStamp}/ABCstat_global.txt',
        '{timeStamp}/posterior.txt'
    output:
        '{timeStamp}/QC_plot/QC_posterior_density.pdf',
        '{timeStamp}/QC_plot/QC_posterior_acp.pdf'
    shell:
        """
        {Sc}/R_visual.sif Rscript {core_path}/QC_posterior.R\
            dir={timeStamp} output_dir={timeStamp}/QC_plot
        """

rule plot_posterior:
    input:
        '{timeStamp}/posterior.txt'
    output:
        '{timeStamp}/visual_posterior.pdf'
    shell:
        """
        {Sc}/R_visual.sif Rscript {core_path}/plot_posterior.R\
            posterior={input} prior_dir={timeStamp}/modelComp output={output} Nref={Nref}
        """

rule plot_model_weight:
    input:
        '{timeStamp}/model_weight.txt'
    output:
        '{timeStamp}/visual_model_weight.pdf'
    shell:
        """
        {Sc}/R_visual.sif Rscript {core_path}/plot_model_weight.R\
            model_weight={input} prior_dir={timeStamp}/modelComp output={output} 
         """
###### checkpoints  ###############################
rule check_modelComp:
    input:
        sim =expand("{timeStamp}/modelComp/{model}_{i}/ABCstat_global.txt",timeStamp=timeStamp,model=MODELS_COMP,i=ITERATIONS_MODEL_COMP)
    output:
        check_token = temp(expand("{timeStamp}/check_modelComp_done.log",timeStamp=timeStamp))
    shell:
        """
        {Sc}/R.sif Rscript {core_path}/check_modelComp.R dir={timeStamp}/modelComp
        touch {output.check_token}
        """


