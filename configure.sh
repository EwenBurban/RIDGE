#!/bin/bash
sudo singularity build container/python.sif container/python.def 
sudo singularity build container/R.sif container/R.def 
sudo singularity build container/scrm_py.sif container/scrm_py.def
sudo singularity build container/R_visual.sif container/R_visual.def
