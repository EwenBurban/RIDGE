#!/bin/bash
singularity build --remote container/python.sif container/python.def 
singularity build --remote container/R.sif container/R.def 
singularity build --remote container/scrm_py.sif container/scrm_py.def
singularity build --remote container/R_visual.sif container/R_visual.def
