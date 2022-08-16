#!/usr/bin/bash
## the provided argument is for --configfile, expecting the yaml file

path=$(realpath $0)
binpath="$( cd -- "$(dirname $path)" >/dev/null 2>&1 ; pwd -P )"

. ${binpath}/config/config.sh
cat ${binpath}/logo
## directory of the code
if [[ $mode != cluster ]]; then 
## without slurm
snakemake --snakefile ${binpath}/input_module/pipeline.py -p -j ${ntask_load} --configfile ${1} --config binpath=${binpath} --latency-wait 60 --nolock 
snakemake --snakefile ${binpath}/gwscan/pipeline.py -p -j ${ntask_load} --configfile ${1} --config binpath=${binpath} --latency-wait 60 --nolock
snakemake --snakefile ${binpath}/core/RIDGE_pipeline.py -p -j ${ntask_load} --until "model_filtering" --configfile ${1}  --config binpath=${binpath}  mode=single --latency-wait 60 --nolock
snakemake --snakefile ${binpath}/core/RIDGE_pipeline.py   -p -j ${ntask_load} --configfile ${1}   --config binpath=${binpath}  mode=single --latency-wait 60 --nolock
snakemake --snakefile ${binpath}/visualisation_module/pipeline.py -p -j ${ntask_load} --configfile ${1} --config binpath=${binpath} --latency-wait 60 --nolock 
else
## cluster version with slurm | If you use the cluster version, you must edit the file cluster.json to adapt it to your cluster
snakemake --snakefile ${binpath}/input_module/pipeline.py -p -j ${ntask_load} --configfile ${1}  --config binpath=${binpath} --cluster-config ${binpath}/input_module/cluster.json --cluster "sbatch --nodes={cluster.node} --ntasks={cluster.n} --cpus-per-task={cluster.cpusPerTask} --time={cluster.time} --mem-per-cpu={cluster.memPerCpu} -p {cluster.p} "  --latency-wait 60 --nolock  
snakemake --snakefile ${binpath}/gwscan/pipeline.py -p -j ${ntask_load} --configfile ${1}  --config binpath=${binpath} --cluster-config ${binpath}/gwscan/cluster.json --cluster "sbatch --nodes={cluster.node} --ntasks={cluster.n} --cpus-per-task={cluster.cpusPerTask} --time={cluster.time} --mem-per-cpu={cluster.memPerCpu} -p {cluster.p} "  --latency-wait 60 --nolock  
snakemake --snakefile ${binpath}/core/RIDGE_pipeline.py -p -j ${ntask_load} --until "model_filtering"  --config binpath=${binpath}  mode=single --configfile ${1} --cluster-config ${binpath}/core/cluster.json --cluster "sbatch --nodes={cluster.node} --ntasks={cluster.ntasks} --cpus-per-task={cluster.cpusPerTask} --time={cluster.time} --mem-per-cpu={cluster.memPerCpu} -p {cluster.p}"  --latency-wait 60  --nolock --restart-times 3  --rerun-incomplete
snakemake --snakefile ${binpath}/core/RIDGE_pipeline.py  -p -j ${ntask_load} --configfile ${1}  --config binpath=${binpath}  mode=single --cluster-config ${binpath}/core/cluster.json --cluster "sbatch --nodes={cluster.node} --ntasks={cluster.ntasks} --cpus-per-task={cluster.cpusPerTask} --time={cluster.time} --mem-per-cpu={cluster.memPerCpu} -p {cluster.p} "  --latency-wait 60  --nolock 
snakemake --snakefile ${binpath}/visualisation_module/pipeline.py -p -j ${ntask_load} --configfile ${1}  --config binpath=${binpath} --cluster-config ${binpath}/visualisation_module/cluster.json --cluster "sbatch --nodes={cluster.node} --ntasks={cluster.n} --cpus-per-task={cluster.cpusPerTask} --time={cluster.time} --mem-per-cpu={cluster.memPerCpu} -p {cluster.p} "  --latency-wait 60 --nolock  
fi



