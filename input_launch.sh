#!/usr/bin/bash
## the provided argument is for --configfile, expecting the yaml file

## load the necessary environement (used for Genouest cluster)
path=$(realpath $0)
binpath="$( cd -- "$(dirname $path)" >/dev/null 2>&1 ; pwd -P )"
. ${binpath}/config/config.sh
## directory of the code
if [[ $mode != cluster ]]; then 
## without slurm
snakemake --snakefile ${binpath}/input_module/pipeline.py -p -j ${ntask_load} --configfile ${1} --config binpath=${binpath} --latency-wait 60 --nolock 
else
## cluster version with slurm | If you use the cluster version, you must edit the file cluster.json to adapt it to your cluster
snakemake --snakefile ${binpath}/input_module/pipeline.py -p -j ${ntask_load} --configfile ${1}  --config binpath=${binpath} --cluster-config ${binpath}/input_module/cluster.json --cluster "sbatch --nodes={cluster.node} --ntasks={cluster.n} --cpus-per-task={cluster.cpusPerTask} --time={cluster.time} --mem-per-cpu={cluster.memPerCpu} -p {cluster.p} "  --latency-wait 60 --nolock  
fi



