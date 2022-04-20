#!/usr/bin/bash
#SBATCH --mail-user=ewen.burban@univ-rennes1.fr
#SBATCH --mail-type=all
#SBATCH --time=20-60:00:00
## launch DILS for 2 populations
## the provided argument is for --configfile, expecting the yaml file
#. /local/env/envconda.sh
#conda activate DILS
. /local/env/envsnakemake-6.0.5.sh
#. /local/env/envr-4.0.3.sh
. /local/env/envr-3.6.2.sh
. /local/env/envsingularity-3.8.0.sh
binpath='/home/genouest/cnrs_umr6553/eburban/GW_DILS/bin'


snakemake --snakefile ${binpath}/BetaTSnakepipe_2pop_locus -p -j 40 --until "model_filtering" --configfile ${1} --cluster-config ${binpath}/cluster_glados.json --cluster "sbatch --nodes={cluster.node} --ntasks={cluster.ntasks} --cpus-per-task={cluster.cpusPerTask} --time={cluster.time} --mem-per-cpu={cluster.memPerCpu}"  --latency-wait 60  --nolock --restart-times 3  --rerun-incomplete && \
for i in {1..18}; do
    snakemake --snakefile ${binpath}/BetaTSnakepipe_2pop_locus -p -j 40 --until "model_averaging" --configfile ${1} --cluster-config ${binpath}/cluster_glados.json --cluster "sbatch --nodes={cluster.node} --ntasks={cluster.ntasks} --cpus-per-task={cluster.cpusPerTask} --time={cluster.time} --mem-per-cpu={cluster.memPerCpu}"  --latency-wait 60  --nolock --restart-times 3 && \
    snakemake --snakefile ${binpath}/BetaTSnakepipe_2pop_locus -p -j 40 --configfile ${1} --cluster-config ${binpath}/cluster_glados.json --cluster "sbatch --nodes={cluster.node} --ntasks={cluster.ntasks} --cpus-per-task={cluster.cpusPerTask} --time={cluster.time} --mem-per-cpu={cluster.memPerCpu} "  --latency-wait 60  --nolock --restart-times 3   # --dag | dot -Tpdf > dag.pdf # -k # --unlock
done


