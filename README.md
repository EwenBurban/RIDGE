# Table of content
---

	1- RIDGE
	2- Installation
	3- Usage
	4- Troubleshooting
	4- License
# 1. RIDGE
---

RIDGE is a tool that identifies genomic regions resistant to gene flow between two populations. It utilizes an ABC approach to infer the demographic history of both populations and subsequently detect barriers. The software, RIDGE, necessitates substantial computational resources and heavily relies on parallelization to reduce its computation duration. Thus, it is strongly recommended to install RIDGE on a cluster, preferably in a Linux environment. To facilitate program management and ensure reproducibility, RIDGE employs *Singularity* container technology and *Snakemake* to handle the workflow. 
For each locus (also known as a fragment of the genome), RIDGE outputs a Bayes factor indicating the likelihood of it acting as a gene flow barrier.
# 2. Installation
---

## 2.1. Get the code

Download the code with the following command 

```

git clone https://github.com/EwenBurban/RIDGE.git

```

## 2.2. Install containers

First, ensure Singularity is installed in your cluster. After that, proceed with: Next, create a remote token on Singularity's website (https://cloud.sylabs.io/tokens) and register it in your cluster by using the command `singularity remote login`.  

```
bash cluster_configure.sh
```

It will take sometime to build up all container

## 2.3. Configuration

Before launching anything, ensure that the config/ folder contains at least one file named `config.sh`. If your cluster uses SLURM as the job manager, you must also fill the `cluster_config.yaml` file. Simply refer to the `example/config` folder and adjust it to match your cluster's configuration. 
For more details, consult the `user_manual.pdf`.

## 2.4 Check installation

TO-DO
# 3. Usage
---
Please refer to the `user_manual.pdf` for information on this section. Note that the file to launch RIDGE is `RIDGE.sh`, and `core/RIDGE_V2_pipeline.py` manages all of the code.

# 4. Troubleshooting
---

If you encounter any issues, report them on GitHub.

# 5. License
---

You can freely use all code in this project under the conditions of the GNU GPL Version 3 or later.

