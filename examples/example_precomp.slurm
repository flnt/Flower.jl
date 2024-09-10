#!/bin/bash
#SBATCH --job-name=seq
#SBATCH --output=%x.o%j
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --mem=80G                       #for nx=512
#SBATCH --partition=cpu_short           # (see available partitions)

# To clean and load modules defined at the compile and link phases
module purge
module restore pdi4 #or load...

# echo of commands
set -x

# To compute in the submission directory
cd ${SLURM_SUBMIT_DIR}


precompile = true
#precompile = false

if [ "$precompile" = true ] ; then
    echo 'Precompiling Flower'
    julia +1.10.5 --project=../flower/Flower.jl --threads=1 --trace-compile=dc_precompile.jl ../flower/Flower.jl/examples/main_current_folder.jl half_circle_pressure_only.yml
else
    echo 'Using precompiled Flower'
    julia +1.10.5 --project=../flower/Flower.jl --threads=1 --sysimage=../flower/Flower.jl/sys_img.so --trace-compile=stderr ../flower/Flower.jl/examples/main_current_folder.jl half_circle.yml
fi
