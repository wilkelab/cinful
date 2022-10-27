#!/bin/bash

<< note
Script is for manual install of the cinful conda environment.
###### This method is not recommended!! #####
See https://github.com/wilkelab/cinful for recommended install.
note


mamba create -n cinful 
mamba env update -n cinful --file cinful_conda.yml
