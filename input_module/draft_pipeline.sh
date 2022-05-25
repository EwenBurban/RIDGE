#!/usr/bin/env sh

. /local/env/envsingularity-3.8.5.sh
wdir="/scratch/eburban/DomIsol/african_rice_test2"
contig_data=${wdir}"/contig_data.txt"
vcf_file=${wdir}"/african_rice_haplotyped.vcf"
binpath="/home/genouest/cnrs_umr6553/eburban/RIDGE/input_module"
container_path="/home/genouest/cnrs_umr6553/eburban/RIDGE/container"
Sc="singularity exec --bind /scratch/eburban ${container_path}"
window_size=10000
Nref=85000
mu=3e-8
popfile=${wdir}"/popfile.csv"
nameA="wild"
nameB="dom"

#${Sc}/R.sif Rscript ${binpath}/generate_bed_sample.R contig_file=${contig_data} window_size=${window_size} nLoci_per_chr=100 output="${wdir}/bed_global_dataset.txt" &&\
#${Sc}/R.sif Rscript ${binpath}/generate_bed_sample.R contig_file=${contig_data} window_size=${window_size} nLoci_per_chr=1000 output="${wdir}/bed_locus_dataset.txt" &&\
${Sc}/R.sif Rscript ${binpath}/get_rec_rate.R contig_file=${contig_data} window_size=${window_size} rho_map=${wdir}/rho_map_african_rice.txt output=${wdir}"/bed_rec_rate.txt" &&\
${Sc}/R.sif Rscript ${binpath}/generate_locus_datafile.R rho_map="${wdir}/bed_rec_rate.txt" bedfile="${wdir}/bed_global_dataset.txt" mu=${mu} Nref=${Nref} window_size=${window_size} nameA=${nameA} nameB=${nameB} popfile=${popfile} output="${wdir}/locus_datafile" &&\
${Sc}/R.sif Rscript ${binpath}/generate_locus_datafile.R rho_map="${wdir}/bed_rec_rate.txt" bedfile="${wdir}/bed_locus_dataset.txt" mu=${mu} Nref=${Nref} window_size=${window_size} nameA=${nameA} nameB=${nameB} popfile=${popfile} output="${wdir}/locus_datafile_locussp" &&\
${Sc}/scrm_py.sif python3 vcf2abc.py  data=${vcf_file} bed_file=${wdir}/bed_global_dataset.txt popfile=${popfile} nameA=${nameA} nameB=${nameB} window_size=${window_size} locus_write="False" global_write="True" output_dir=${wdir}&&\
${Sc}/scrm_py.sif python3 vcf2abc.py  data=${vcf_file} bed_file=${wdir}/bed_locus_dataset.txt popfile=${popfile} nameA=${nameA} nameB=${nameB} window_size=${window_size} locus_write="True" global_write="False" output_dir=${wdir}



