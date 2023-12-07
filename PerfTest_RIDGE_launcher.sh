#!/usr/bin/bash
## the provided argument is for --configfile, expecting the yaml file
if [[ -n $SLURM_JOB_ID ]] ; then
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
snakemake --snakefile ${binpath}/core/RIDGE_perftest_pipe.py   -p -j ${ntask_load} --configfile ${1}   --config binpath=${binpath}   --latency-wait 60 --nolock
else
snakemake --snakefile ${binpath}/core/RIDGE_perftest_pipe.py -p -j ${ntask_load}\
	--config binpath=${binpath}  --configfile ${1} mode=${2} --cluster-config ${binpath}/config/cluster.json \
	--cluster "sbatch -A ${ID} --nodes={cluster.node} --ntasks={cluster.ntasks} --cpus-per-task={cluster.cpusPerTask} --time={cluster.time} --mem-per-cpu={cluster.memPerCpu} -p {cluster.p}" \
	--latency-wait 60  --nolock -r  
fi

