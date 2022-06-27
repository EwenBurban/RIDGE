#!/usr/bin/bash
#SBATCH --mail-user=ewen.burban@univ-rennes1.fr
#SBATCH --mail-type=all
#SBATCH --time=20-60:00:00
## launch DILS for 2 populations
## the provided argument is for --configfile, expecting the yaml file
#. /local/env/envconda.sh
#conda activate DILS
. /local/env/envsnakemake-6.0.5.sh
. /local/env/envsingularity-3.8.0.sh
binpath='/home/genouest/cnrs_umr6553/eburban/RIDGE/core'


snakemake --snakefile ${binpath}/RIDGE_pipe -p -j 70 --until "model_filtering" --configfile ${1} --cluster-config ${binpath}/cluster_2pop.json --cluster "sbatch --nodes={cluster.node} --ntasks={cluster.ntasks} --cpus-per-task={cluster.cpusPerTask} --time={cluster.time} --mem-per-cpu={cluster.memPerCpu} -p {cluster.p}"  --latency-wait 60  --nolock --restart-times 3  --rerun-incomplete
snakemake --snakefile ${binpath}/RIDGE_pipe  -p -j 70 --configfile ${1} --cluster-config ${binpath}/cluster_2pop.json --cluster "sbatch --nodes={cluster.node} --ntasks={cluster.ntasks} --cpus-per-task={cluster.cpusPerTask} --time={cluster.time} --mem-per-cpu={cluster.memPerCpu} -p {cluster.p} "  --latency-wait 60  --nolock 


