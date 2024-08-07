#!/bin/bash
#SBATCH --job-name=seq
#SBATCH --output=%x.o%j
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --partition=cpu_short           # (see available partitions)

# To clean and load modules defined at the compile and link phases
module purge
module load ...

# echo of commands
set -x

# To compute in the submission directory
cd ${SLURM_SUBMIT_DIR}

# execution
julia +1.10.0 --project=../Flower.jl --threads=1 ../Flower.jl/examples/main_current_folder.jl ../Flower.jl/examples/validation.yml

# Post-processing
# python3 -c "import plot_flower; plot_flower.plot_all_h5()" yaml_file_name


# srun --nodes=1 --time=00:30:00 -p cpu_short --pty /bin/bash
# source /gpfs/workdir/shared/pdi-deisa/pdi_tuto_env
# module load pdiplugin-pycall/1.6.0/gcc-11.2.0

# Maybe necesary:
# julia --project -e 'using Pkg; Pkg.add("MPIPreferences")'
# julia --project -e 'using MPIPreferences; MPIPreferences.use_system_binary()'


# julia +1.10.4 --project=$HOME/flower/Flower.jl --threads=1 $HOME/flower/Flower.jl/examples/main_current_folder.jl $HOME/flower/Flower.jl/examples/validation.yml


