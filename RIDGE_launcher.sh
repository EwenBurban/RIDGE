#!/usr/bin/bash
#SBATCH --mail-user=ewen.burban@univ-rennes1.fr
#SBATCH --mail-type=all
#SBATCH --time=12-60:00:00
## launch DILS for 2 populations
## the provided argument is for --configfile, expecting the yaml file
#. /local/env/envconda.sh
#conda activate DILS
. /local/env/envsnakemake-6.0.5.sh
#. /local/env/envr-4.0.3.sh
. /local/env/envr-3.6.2.sh
. /local/env/envsingularity-3.8.0.sh
binpath='/home/genouest/cnrs_umr6553/eburban/GW_DILS/bin'
#snakemake --snakefile ${binpath}/Best_Snake_2pop_locus -p -j 35 --configfile ${1} --cluster-config ${binpath}/cluster_2pop.json --cluster "sbatch --nodes={cluster.node} --ntasks={cluster.n} --cpus-per-task={cluster.cpusPerTask} --time={cluster.time} --mem-per-cpu={cluster.memPerCpu} -p {cluster.p} "  --latency-wait 60   --rerun-incomplete --nolock  -n # --unlock


#snakemake --snakefile ${binpath}/TSnakepipe_2pop_locus -p -j 135 --configfile ${1} --cluster-config ${binpath}/cluster_2pop.json --cluster "sbatch --nodes={cluster.node} --ntasks={cluster.n} --cpus-per-task={cluster.cpusPerTask} --time={cluster.time} --mem-per-cpu={cluster.memPerCpu} -p {cluster.p} "  --latency-wait 60   --rerun-incomplete --nolock -n # -k # --unlock


snakemake --snakefile ${binpath}/BetaTSnakepipe_2pop_locus -p -j 40 --until "model_filtering" --configfile ${1} --cluster-config ${binpath}/cluster_2pop.json --cluster "sbatch --nodes={cluster.node} --ntasks={cluster.ntasks} --cpus-per-task={cluster.cpusPerTask} --time={cluster.time} --mem-per-cpu={cluster.memPerCpu} -p {cluster.p} "  --latency-wait 60  --nolock --restart-times 1  --rerun-incomplete && \
snakemake --snakefile ${binpath}/BetaTSnakepipe_2pop_locus -p -j 40 --until "model_averaging" --configfile ${1} --cluster-config ${binpath}/cluster_2pop.json --cluster "sbatch --nodes={cluster.node} --ntasks={cluster.ntasks} --cpus-per-task={cluster.cpusPerTask} --time={cluster.time} --mem-per-cpu={cluster.memPerCpu} -p {cluster.p} "  --latency-wait 60  --nolock --restart-times 1 --rerun-incomplete 
snakemake --snakefile ${binpath}/BetaTSnakepipe_2pop_locus -p -j 40 --configfile ${1} --cluster-config ${binpath}/cluster_2pop.json --cluster "sbatch --nodes={cluster.node} --ntasks={cluster.ntasks} --cpus-per-task={cluster.cpusPerTask} --time={cluster.time} --mem-per-cpu={cluster.memPerCpu} -p {cluster.p} "  --latency-wait 60  --nolock --restart-times 1   # --dag | dot -Tpdf > dag.pdf # -k # --unlock



