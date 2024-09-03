# Flower

## Precompiling Flower
cf https://vsert.com/posts/precompile/
Use runtests_precomp.jl instead of runtests.jl otherwise the precompilation is much longer
```bash
julia +1.10.4 --project=../Flower.jl --threads=1 --trace-compile=dc_precompile.jl ../Flower.jl/test/runtests_precomp.jl
julia --project=../Flower.jl --trace-compile=dc_precompile.jl ../Flower.jl/src/PackageCompiler_setup_image.jl

time julia --project=../Flower.jl --sysimage=sys_img.so -e 'using Flower; dostuff()'

time julia +1.10.4 --project=../Flower.jl --threads=1 --sysimage=sys_img.so ../Flower.jl/examples/main_current_folder.jl similar_to_Khalighi.yml
```
Example:
```bash
julia --project=../Flower.jl --trace-compile=dc_precompile.jl ../Flower.jl/src/PackageCompiler_setup_image.jl
✔ [02m:08s] PackageCompiler: compiling incremental system image

time julia +1.10.4 --project=../Flower.jl --threads=1 --sysimage=sys_img.so ../Flower.jl/examples/main_current_folder.jl similar_to_Khalighi.yml
```
General example

```bash
time julia --project=../Flower.jl --threads=1 --trace-compile=dc_precompile.jl ../Flower.jl/examples/example.jl

cd ../Flower.jl

julia +1.10.4 --project=. src/PackageCompiler_setup_image.jl

time julia +1.10.4 --project=../Flower.jl --threads=1 --sysimage=../Flower.jl/sys_img.so --trace-compile=stderr ../Flower.jl/examples/example.jl
```

Example
```bash
time julia +1.10.4 --project=../Flower.jl --threads=1 --trace-compile=dc_precompile.jl ../Flower.jl/examples/main_current_folder.jl similar_to_Khalighi.yml

cd ../Flower.jl

julia +1.10.4 --project=. src/PackageCompiler_setup_image.jl

time julia +1.10.4 --project=../Flower.jl --threads=1 --sysimage=../Flower.jl/sys_img.so --trace-compile=stderr ../Flower.jl/examples/main_current_folder.jl similar_to_Khalighi.yml
```
With --sysimage (typical execution times on a laptop): 
```bash
real	0m21,392s
user	0m25,334s
sys	0m18,700s
```

<!-- real	0m14,879s
user	0m19,732s
sys	0m16,808s -->

Without
```bash
real	3m18,956s
user	3m22,081s
sys	0m17,694s
```


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
python3 -c "import plot_flower; plot_flower.plot_all_films_func()" yaml_file
```
To plot the radius of the interface
```bash
python3 -c "import plot_flower; plot_flower.plot_radius_from_h5()" yaml_file
```



# Testing
In test : runtests.jl
