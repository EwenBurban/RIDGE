#!/usr/bin/bash
config=$1
. /local/env/envsnakemake-6.0.5.sh
. /local/env/envsingularity-3.8.0.sh
binpath='/home/genouest/cnrs_umr6553/eburban/RIDGE/faker_module'
snakemake --snakefile ${binpath}/Faker_pipe -p -j 35 --configfile $config --cluster "sbatch -p ecobio --mem-per-cpu=5000" --latency-wait 160 --nolock --rerun-incomplete  # --use-singularity --singularity-args '--bind ${binpath},${worpath}' --singularity-prefix ${container_path} --use-conda # --unlock



