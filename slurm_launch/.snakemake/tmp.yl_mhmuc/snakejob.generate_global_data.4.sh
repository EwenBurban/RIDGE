#!/bin/sh
# properties = {"type": "single", "rule": "generate_global_data", "local": false, "input": ["/scratch/eburban/DomIsol/test_RIDGE/locus_datafile", "/scratch/eburban/DomIsol/test_RIDGE/config.yaml"], "output": ["/scratch/eburban/DomIsol/test_RIDGE/simdata/IM_2M_2N_0_Ts~10000_Pb~0.1_Mb~10/ABCstat_global.txt", "/scratch/eburban/DomIsol/test_RIDGE/simdata/IM_2M_2N_0_Ts~10000_Pb~0.1_Mb~10/ABCstat_locus.txt", "/scratch/eburban/DomIsol/test_RIDGE/simdata/IM_2M_2N_0_Ts~10000_Pb~0.1_Mb~10/priorfile.txt", "/scratch/eburban/DomIsol/test_RIDGE/simdata/IM_2M_2N_0_Ts~10000_Pb~0.1_Mb~10/priorfile_locus.txt"], "wildcards": {"timeStamp": "/scratch/eburban/DomIsol/test_RIDGE", "sim_model": "IM_2M_2N", "n": "0", "Tsplit": "10000", "Pbar": "0.1", "Mbasal": "10"}, "params": {}, "log": [], "threads": 4, "resources": {"mem_mb": 10000}, "jobid": 4, "cluster": {}}
 cd /home/genouest/cnrs_umr6553/eburban/RIDGE/slurm_launch && \
PATH='/local/miniconda3/envs/snakemake-6.0.5/bin':$PATH /local/miniconda3/envs/snakemake-6.0.5/bin/python3.9 \
-m snakemake '/scratch/eburban/DomIsol/test_RIDGE/simdata/IM_2M_2N_0_Ts~10000_Pb~0.1_Mb~10/priorfile.txt' --snakefile /home/genouest/cnrs_umr6553/eburban/RIDGE/faker_module/Faker_pipe \
--force -j --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /home/genouest/cnrs_umr6553/eburban/RIDGE/slurm_launch/.snakemake/tmp.yl_mhmuc /scratch/eburban/DomIsol/test_RIDGE/locus_datafile /scratch/eburban/DomIsol/test_RIDGE/config.yaml --latency-wait 160 \
 --attempt 1 --force-use-threads --scheduler ilp \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /home/genouest/cnrs_umr6553/eburban/RIDGE/config/test_RIDGE.yaml -p --allowed-rules generate_global_data --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /home/genouest/cnrs_umr6553/eburban/RIDGE/slurm_launch/.snakemake/tmp.yl_mhmuc/4.jobfinished || (touch /home/genouest/cnrs_umr6553/eburban/RIDGE/slurm_launch/.snakemake/tmp.yl_mhmuc/4.jobfailed; exit 1)

