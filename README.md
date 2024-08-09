# Flower




## On HPC center

### Installation
```bash
curl -fsSL https://install.julialang.org | sh
```
For example, in $HOME/flower
* Clone Flower.jl
*  Clone DDM.jl

what is added in project in $HOME is not used in $WORKDIR, so in $WORKDIR do:
```bash
julia +1.10.4 –project=$HOME/flower/Flower.jl 
import Pkg ; Pkg.add(package)
```

### Run
In interactive session:
```bash
julia +1.10.4 --project=$HOME/flower/Flower.jl --threads=1 $HOME/flower/Flower.jl/examples/main_current_folder.jl $HOME/flower/Flower.jl/examples/validation.yml
```

The YAML file contains the main parameters for the simulation, IO and post-processing.

## Post-processing 
In the folder where the results are stored:
To plot all figures specified in the yaml file:
```bash
python3 -c "import plot_flower; plot_flower.plot_all_fig()" yaml_file
```
To produce all films:
```bash
python3 -c "import plot_flower; plot_flower.plot_all_films()" yaml_file
```
To plot the radius of the interface
```bash
python3 -c "import plot_flower; plot_flower.plot_radius_from_h5()" yaml_file
```



<!-- # Testing -->

