#!/bin/sh

conda create --name cinful python=3.8.13 pip
conda env update -n cinful --file cinful_conda.yml
