# Flower




## On HPC

## Installation
```bash
curl -fsSL https://install.julialang.org | sh
```
For example, in $HOME/flower
Clone Flower.jl
Clone DDM.jl

what is added in project in $HOME is not used in $WORKDIR, so in $WORKDIR do:
julia +1.10.4 –project=$HOME/flower/Flower.jl 
import Pkg ; Pkg.add(package)

### Run
In interactive session:
```bash
julia +1.10.4 --project=$HOME/flower/Flower.jl --threads=1 $HOME/flower/Flower.jl/examples/main_current_folder.jl $HOME/flower/Flower.jl/examples/validation.yml
```

The YAML file contains the main parameters for the simulation, IO and post-processing.

## Post-processing 
In the folder where the results are stored:
```bash
python3 -c "import plot_flower; plot_flower.plot_all_fig()" half_circle.yml
```


<!-- # Testing -->

