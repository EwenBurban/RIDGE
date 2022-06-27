## import module and info from config file
import os
import sys
import re
import pandas as pd
binpath=config['binpath'] 
gwscan_path = binpath + '/gwscan'

nameA = config['nameA']
nameB = config['nameB']
timeStamp = config['timeStamp']
window_size = config['window_size']

config_yaml = timeStamp + '/'  +config['config_yaml']
popfile= timeStamp + '/'  +  config['popfile']
contig_file=  timeStamp + '/'  + config['contig_file']
rec_rate_map =  timeStamp + '/'  + config['rec_rate_map']
gff_file=timeStamp + '/' + config['gff_file']
Nthread_shape_it=8
tag=config['tag']
if tag == 'NULL':
    tag =""
## Singularity config ##
container_path =binpath + '/container' 
Sc='singularity exec --bind {0},{1} {2}'.format(binpath,timeStamp,container_path)

### get contig names ###
contig_data = pd.read_csv(contig_file,sep='\t')
contig_name = contig_data['contig_name'].to_list()
## pipeline ##
wildcard_constraints:
    timeStamp=timeStamp,
    tag=tag,
    nameA=nameA,
    nameB=nameB,
    contig_name = '|'.join(contig_name),
    pop_name='|'.join([nameA,nameB])


rule targets:
    input:
        haplotyped = expand('{timeStamp}/haplotyped.vcf',timeStamp=timeStamp),
        abc = expand('{timeStamp}/gwscan/ABCstat_locus{tag}.txt',timeStamp=timeStamp,tag=tag),
        gwscan = expand('{timeStamp}/gwscan/gwscan{tag}.txt',timeStamp=timeStamp,tag=tag)
    shell:
        """
        echo 'done'
        """


rule calculate_globalstat:
    input:
        vcf = expand('{timeStamp}/haplotyped.vcf',timeStamp=timeStamp),
        popfile = popfile,
        contig_file = contig_file
    output:
        '{timeStamp}/gwscan/ABCstat_locus{tag}.txt'
    shell:
        """
         {Sc}/scrm_py.sif python3 {gwscan_path}/vcf2abc.py data={input.vcf} \
            popfile={input.popfile} contig_file={input.contig_file} nameA={nameA}\
            nameB={nameB} window_size={window_size} output={output}
        """
########### detect sweeps with RAiSD program ######
rule detect_sweep:
    input:
        vcf= expand('{timeStamp}/haplotyped.vcf',timeStamp=timeStamp),
        popfile = popfile
    output:
        expand('{{timeStamp}}/gwscan/raisd_files/RAiSD_Report.{pop_name}.{contig_name}',pop_name=[nameA,nameB],contig_name=contig_name) 
    shell:
        """
        cd {timeStamp}/gwscan/raisd_files
        IFS=',' read -r -a pops <<< $(head -n1 {input.popfile})
        sed '1d' {input.popfile} | cut -f 1 -d ',' > tmp_${{pops[0]}}
        {Sc}/raisd.sif /home/RAiSD -n ${{pops[0]}} -I {input.vcf} -S tmp_${{pops[0]}} -y 1 -s -t -f
        sed '1d' {input.popfile} | cut -f 2 -d ',' > tmp_${{pops[1]}}
        {Sc}/raisd.sif /home/RAiSD -n ${{pops[1]}} -I {input.vcf} -S tmp_${{pops[1]}} -y 1 -s -t -f
        """

rule aggregate_sweep:
    input:
        raisd_files = expand('{{timeStamp}}/gwscan/raisd_files/RAiSD_Report.{{pop_name}}.{contig_name}',contig_name=contig_name),
        contig_file = contig_file
    output:
        expand('{{timeStamp}}/gwscan/{{pop_name}}_raisd{tag}.txt',tag=tag)
    shell:
        """
        {Sc}/R.sif Rscript {gwscan_path}/aggr2window_raisd.R pop={wildcards.pop_name} input_dir={wildcards.timeStamp}/gwscan/raisd_files\
                output={output} window_size={window_size} contig_file={input.contig_file}
        """
########## detect locus contributing to differentiate populations with
########## pcadapt package on R
rule pcadapt_detect:
    input:
        vcf = expand('{timeStamp}/haplotyped.vcf',timeStamp=timeStamp),
        contig_file = contig_file,
        popfile = popfile
    output:
        expand('{timeStamp}/gwscan/pcadapt_detect{tag}.txt',timeStamp=timeStamp,tag=tag)
    shell:
        """
        {Sc}/pcadapt.sif Rscript {gwscan_path}/windowed_pcadapt.R vcf={input.vcf} window_size={window_size} output={output}\
                contig_file={input.contig_file} popfile={input.popfile} nameA={nameA} nameB={nameB}
        """

########## calculating coding rate based on exon position
rule get_coding_rate:
    input:
        gff_file = gff_file,
        contig_file = contig_file
    output:
        expand('{timeStamp}/gwscan/coding_rate{tag}.txt',timeStamp=timeStamp,tag=tag)
    shell:
        """
        {Sc}/R.sif Rscript {gwscan_path}/get_coding_rate.R gff_file={input.gff_file}\
                contig_file={input.contig_file} output={output} window_size={window_size}
        """
########## get mean rec rate per basis from a recombination map

rule get_rec_rate:
    input:
        rec_rate_map = rec_rate_map,
        contig_file =contig_file
    output:
        expand('{timeStamp}/gwscan/rec_rate{tag}.txt',timeStamp=timeStamp,tag=tag)
    shell:
        """
        {Sc}/R.sif Rscript {gwscan_path}/get_rec_rate.R rho_map={input.rec_rate_map} contig_file={input.contig_file}\
                output={output} window_size={window_size}
        """
##### merge all the data together
rule merge_all:
    input:
        abc = expand('{timeStamp}/gwscan/ABCstat_locus{tag}.txt',timeStamp=timeStamp,tag=tag),
        raisd_popA = expand('{timeStamp}/gwscan/{nameA}_raisd{tag}.txt',nameA=nameA,timeStamp=timeStamp,tag=tag),
        raisd_popB = expand('{timeStamp}/gwscan/{nameB}_raisd{tag}.txt',nameB=nameB,timeStamp=timeStamp,tag=tag),
        pcadapt = expand('{timeStamp}/gwscan/pcadapt_detect{tag}.txt',timeStamp=timeStamp,tag=tag),
        rec_rate=expand('{timeStamp}/gwscan/rec_rate{tag}.txt',timeStamp=timeStamp,tag=tag),
        coding_rate=expand('{timeStamp}/gwscan/coding_rate{tag}.txt',timeStamp=timeStamp,tag=tag)
    output: 
        expand('{timeStamp}/gwscan/gwscan{tag}.txt',timeStamp=timeStamp,tag=tag)
    shell:
        """
        {Sc}/R.sif Rscript {gwscan_path}/data_merging.R abc={input.abc} rec_rate={input.rec_rate}\
                coding_rate={input.coding_rate} sweep_A={input.raisd_popA} sweep_B={input.raisd_popB}\
                pcadapt={input.pcadapt} output={output}  
        """
