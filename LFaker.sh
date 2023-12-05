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
snakemake --snakefile ${binpath}/faker_module/Faker_pipe -p -j ${ntask_load} \
	--configfile ${1} --config binpath=${binpath} --cluster-config ${binpath}/config/cluster.json \
	--cluster "sbatch -A ${ID} --nodes={cluster.node} --ntasks={cluster.ntasks} --cpus-per-task={cluster.cpusPerTask} --time={cluster.time} --mem-per-cpu={cluster.memPerCpu} -p {cluster.p}" \
	--latency-wait 60 --nolock -r 
