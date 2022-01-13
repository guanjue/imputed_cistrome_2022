#!/bin/bash

#SBATCH --mem=200G
#SBATCH --cpus-per-task=1

source ~/.bashrc

cd $1

echo ok 
