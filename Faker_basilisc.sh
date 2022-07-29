#!/usr/bin/bash
## the provided argument is for --configfile, expecting the yaml file

## load the necessary environement (used for Genouest cluster)
if [ -n $SLURM_JOB_ID ] ; then
path=$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}')
scontrol show job $SLURM_JOBID
IFS=' '
read -ra path_arr <<< "$path"
binpath="$( cd -- "$(dirname ${path_arr[0]})" >/dev/null 2>&1 ; pwd -P )"
else
path=$(realpath $0)
binpath="$( cd -- "$(dirname $path)" >/dev/null 2>&1 ; pwd -P )"
fi
. ${binpath}/config/config.sh
cat ${binpath}/logo
## directory of the code
if [[ $mode != cluster ]]; then 
## without slurm
snakemake --snakefile ${binpath}/core/RIDGE_pipeline.py -p -j ${ntask_load} --until "model_filtering" --configfile ${1}  --config binpath=${binpath}  mode=single --latency-wait 60 --nolock
snakemake --snakefile ${binpath}/core/RIDGE_pipeline.py   -p -j ${ntask_load} --configfile ${1}   --config binpath=${binpath}  mode=single --latency-wait 60 --nolock
else
## cluster version with slurm | If you use the cluster version, you must edit the file cluster.json to adapt it to your cluster
snakemake --snakefile ${binpath}/core/RIDGE_pipeline.py -p -j ${ntask_load} --until "model_filtering"  --config binpath=${binpath}  mode=single --configfile ${1} --cluster-config ${binpath}/core/cluster.json --cluster "sbatch --nodes={cluster.node} --ntasks={cluster.ntasks} --cpus-per-task={cluster.cpusPerTask} --time={cluster.time} --mem-per-cpu={cluster.memPerCpu} -p {cluster.p}"  --latency-wait 60  --nolock --restart-times 3  --rerun-incomplete
snakemake --snakefile ${binpath}/core/RIDGE_pipeline.py  -p -j ${ntask_load} --configfile ${1}  --config binpath=${binpath}  mode=single --cluster-config ${binpath}/core/cluster.json --cluster "sbatch --nodes={cluster.node} --ntasks={cluster.ntasks} --cpus-per-task={cluster.cpusPerTask} --time={cluster.time} --mem-per-cpu={cluster.memPerCpu} -p {cluster.p} "  --latency-wait 60  --nolock 
fi



