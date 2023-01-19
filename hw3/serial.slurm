#!/bin/bash

# This script will request one GPU device and 1 CPU core

#SBATCH --job-name=serial
#SBATCH --output=serial.out
#SBATCH --error=serial.err
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cluster=gpu
#SBATCH --partition=a100_multi
#SBATCH --gres=gpu:1

module purge
module load gcc/5.4.0
module load cuda

./serial -n 2000 -o serial.txt

