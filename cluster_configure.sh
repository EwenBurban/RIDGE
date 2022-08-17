#!/bin/bash
singularity build --remote container/bcftools.sif container/bcftools.def 
singularity build --remote container/python.sif container/python.def 
singularity build --remote container/shapeit.sif container/shapeit.def 
singularity build --remote container/R.sif container/R.def 
singularity build --remote container/scrm_py.sif container/scrm_py.def
singularity build --remote container/raisd.sif container/raisd.def
singularity build --remote container/pcadapt.sif container/pcadapt.def
singularity build --remote container/R_visual.sif container/R_visual.def
singularity build --remote container/R_visual.sif docker://biocontainers/plink1.9:v1.90b6.6-181012-1-deb_cv1
