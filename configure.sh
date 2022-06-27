#!/bin/bash
sudo singularity build container/bcftools.sif container/bcftools.def 
sudo singularity build container/python.sif container/python.def 
sudo singularity build container/shapeit.sif container/shapeit.def 
sudo singularity build container/R.sif container/R.def 
sudo singularity build container/scrm_py.sif container/scrm_py.def
sudo singularity build container/raisd.sif container/raisd.def
sudo singularity build container/pcadapt.sif container/pcadapt.def
