
binpath = config['binpath']
visual_path = binpath + '/visualisation_module'
lightMode = config['lightMode']
timeStamp = config['timeStamp']
if lightMode == False:
    last_est_model_dir = '{}/estimation_param_5'.format(timeStamp)
else:
    last_est_model_dir = '{}/estimation_param_2'.format(timeStamp)
############# singularity parametrisaiton #########
container_path =binpath + '/' + 'container' 
Sc='singularity exec --bind {0},{1} {2}'.format(binpath,timeStamp,container_path)


rule targets:
    input:
        expand('{timeStamp}/visualisation/raw_boxplot.pdf',timeStamp=timeStamp),
        expand('{timeStamp}/visualisation/first_round_plot.pdf',timeStamp=timeStamp),
        expand('{timeStamp}/visualisation/last_round_plot.pdf',timeStamp=timeStamp)
    shell:
        """
        echo visualisation done
        """

rule raw_boxplot:
    input:
        locus_data = '{timeStamp}/ABCstat_locus.txt',
        obs_data = '{timeStamp}/ABCstat_global.txt'
    output:
        '{timeStamp}/visualisation/raw_boxplot.pdf'
    shell:
        """
        {Sc}/R_visual.sif Rscript {visual_path}/raw_boxplot.R locus_data={input.locus_data}\
                obs_data={input.obs_data} output={output}
        """

rule first_round:
    input:
        obs_data = '{timeStamp}/ABCstat_global.txt',
        first_round_dir = directory('{timeStamp}/modelComp'),
        est_model_dir = directory('{timeStamp}/estimation_param_1')
    output:
        '{timeStamp}/visualisation/first_round_plot.pdf'
    shell:
        """
        {Sc}/R_visual.sif Rscript {visual_path}/first_round_plots.R obs_data={input.obs_data}\
            first_round_dir={input.first_round_dir} est_model_dir={input.est_model_dir} output={output}
        """

rule last_round:
    input:
        obs_data = '{timeStamp}/ABCstat_global.txt',
        first_round_dir = directory('{timeStamp}/modelComp'),
        est_model_dir = directory(last_est_model_dir),
        last_round_dir = directory(last_est_model_dir),
        model_weight = '{timeStamp}/models_weight.txt',
        average_posterior_file = '{timeStamp}/locus_posteriors_mw.txt'
    output:
        '{timeStamp}/visualisation/last_round_plot.pdf'
    shell:
        """
        {Sc}/R_visual.sif Rscript {visual_path}/last_round_plots.R obs_data={input.obs_data}\
            first_round_dir={input.first_round_dir} est_model_dir={input.est_model_dir}\
            last_round_dir={input.last_round_dir} average_posterior_file={input.average_posterior_file}\
            output={output} model_weight={input.model_weight}
        """
