# links to the codes
binpath = config['binpath']
input_path = binpath + '/input_module'
# informations from the config.yaml file
nameA = config['nameA']
nameB = config['nameB']
timeStamp = config['timeStamp']
Nref = (0 + config['N_max'])/2.0 # Nref is the mid point of the prior
mu = config['mu']
window_size = config['window_size']
popfile= timeStamp + '/'  +  config['popfile']
contig_file=  timeStamp + '/'  + config['contig_file']
rec_rate_map =  timeStamp + '/'  + config['rec_rate_map']
Nthread_shape_it=8
############# vcfile gather ######################
vcf_list_files = os.listdir(timeStamp + '/rawdata')
pattern= re.compile('.vcf')
vcf_list_files = [ x for x in vcf_list_files if pattern.search(x) ]
############# singularity parametrisaiton #########
container_path =binpath + '/' + 'container' 
Sc='singularity exec --bind {0},{1} {2}'.format(binpath,timeStamp,container_path)
############### wildcards management ############## 
# Prevent Snakemake to use regex behavior on wildcards
wildcard_constraints:
    timeStamp=timeStamp,
    binpath=binpath

############## End of Pipeline & Targets ############
rule targets: # edit at the end 
    input:
        bpfile = expand("{timeStamp}/locus_datafile", timeStamp=timeStamp),
        ABCstat_global = expand("{timeStamp}/ABCstat_global.txt",timeStamp=timeStamp),
        ABCstat_locus = expand("{timeStamp}/ABCstat_locus.txt",timeStamp=timeStamp),
    shell:
            """
            echo 'done'
            """

rule format_vcf:
    input:
        '{timeStamp}/rawdata/{vcf}'
    output:
        '{timeStamp}/formated/{vcf}'
    shell:
        """
        {Sc}/bcftools.sif bcftools annotate -x 'FORMAT' {input} | cat > {output}        
        """

rule filtering:
    input:
        '{timeStamp}/formated/{vcf}'
    output:
        '{timeStamp}/filtered/{vcf}'
    shell:
        """
        {Sc}/bcftools.sif bcftools view --max-alleles 2 --exclude-types indels {input} | cat > {output}
        """

rule phasing:
    input:
        '{timeStamp}/filtered/{vcf}'
    output:
        '{timeStamp}/phased/{vcf}'
    threads: Nthread_shape_it
    shell:
        """
        {Sc}/shapeit.sif shapeit --input-vcf {input} -O {output} -T {Nthread_shape_it} 
        {Sc}/shapeit.sif shapeit -convert --input-haps {output} --output-vcf {output}
        """

rule aggregate:
    input: expand('{timeStamp}/phased/{vcf}',timeStamp=timeStamp,vcf=vcf_list_files)
    output: 
        expand('{timeStamp}/all_phased.vcf',timeStamp=timeStamp)
    shell:
        """
        a=({timeStamp}/phased/*)
        grep -h '#' $a > {output}
        grep -hv '#' {timeStamp}/phased/*.vcf >> {output}
        """


rule haplotyping:
    input:
        expand('{timeStamp}/all_phased.vcf',timeStamp=timeStamp)
    output:
        expand('{timeStamp}/haplotyped.vcf',timeStamp=timeStamp)
    resources:
        mem_mb=2*len(vcf_list_files)*2000
    shell:
        """
        {Sc}/python.sif python3 {input_path}/make_haplo_vcf.py file={input} output={output} popfile={popfile}
        """

rule generate_bed_file:
    input:
        contig_file = contig_file
    params:
        nloci_per_chr_global=100,
        nloci_per_chr_locus=1000,
        nloci_per_chr_full_locus=-1 # here -1 mean there is no sampling, all locus are analysed. 
    output:
        bed_global = '{timeStamp}/bed_global_dataset.txt',
        bed_locus = '{timeStamp}/bed_locus_dataset.txt',
        bed_full_locus = '{timeStamp}/bed_full_locus_dataset.txt'
    shell:
        """
            {Sc}/R.sif Rscript {input_path}/generate_bed_sample.R contig_file={contig_file}\
                    window_size={window_size} nLoci_per_chr={params.nloci_per_chr_global}\
                    output={output.bed_global}
            {Sc}/R.sif Rscript {input_path}/generate_bed_sample.R contig_file={contig_file}\
                    window_size={window_size} nLoci_per_chr={params.nloci_per_chr_locus}\
                    output={output.bed_locus}
            {Sc}/R.sif Rscript {input_path}/generate_bed_sample.R contig_file={contig_file}\
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
            {Sc}/R.sif Rscript {input_path}/get_rec_rate.R contig_file={contig_file}\
                    window_size={window_size} rho_map={rec_rate_map} output={output}
        """

rule generate_locus_datafile:
    input:
        bed_rec_rate = '{timeStamp}/bed_rec_rate.txt',
        contig_file = contig_file,
        popfile = popfile,
        bed_global = '{timeStamp}/bed_global_dataset.txt',
        bed_locus = '{timeStamp}/bed_locus_dataset.txt',
    output:
        datafile = '{timeStamp}/locus_datafile',
        datafile_locussp = '{timeStamp}/locus_datafile_locussp'
    shell:
        """
            {Sc}/R.sif Rscript {input_path}/generate_locus_datafile.R \
                    rho_map={input.bed_rec_rate} bedfile={input.bed_global} mu={mu}\
                    Nref={Nref} window_size={window_size} nameA={nameA} nameB={nameB}\
                    popfile={popfile} output={output.datafile}
            {Sc}/R.sif Rscript {input_path}/generate_locus_datafile.R \
                    rho_map={input.bed_rec_rate} bedfile={input.bed_locus} mu={mu}\
                    Nref={Nref} window_size={window_size} nameA={nameA} nameB={nameB}\
                    popfile={popfile} output={output.datafile_locussp}
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
            {Sc}/scrm_py.sif python3 {input_path}/vcf2abc.py  data={input.vcf_file} bed_file={input.bed_global}\
                    popfile={popfile} nameA={nameA} nameB={nameB} window_size={window_size}\
                    locus_write="False" global_write="True" output_dir={timeStamp}
            {Sc}/scrm_py.sif python3 {input_path}/vcf2abc.py  data={input.vcf_file} bed_file={input.bed_locus}\
                    popfile={popfile} nameA={nameA} nameB={nameB} window_size={window_size}\
                    locus_write="True" global_write="False" output_dir={timeStamp}
        """
#
