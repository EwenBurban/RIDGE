#!/bin/bash
sudo singularity build container/bcftools.sif container/bcftools.def 
sudo singularity build container/python.sif container/python.def 
sudo singularity build container/shapeit.sif container/shapeit.def 
sudo singularity build container/R.sif container/R.def 
sudo singularity build container/scrm_py.sif container/scrm_py.def
sudo singularity build container/raisd.sif container/raisd.def
sudo singularity build container/pcadapt.sif container/pcadapt.def
sudo singularity build container/R_visual.sif container/R_visual.def
sudo singularity build container/plink.sif docker://biocontainers/plink1.9:v1.90b6.6-181012-1-deb_cv1
