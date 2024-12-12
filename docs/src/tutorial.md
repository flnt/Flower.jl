# Tutorial

## Generate documentation

```bash
julia +1.10.5 --project=../Flower.jl docs/make.jl
```

## Improving type declarations

Profiling with BenchmarkTools:ctrl shif p in VScode, @bprofile and using ProfileSVG

```bash
module testprofile
ARGS = String["levelset_Butler_two_LS.yml"]
include("../Flower.jl/examples/profile_compile.jl")
end
```

Profile Flower.jl with /examples/profile.jl: run the function "run_forward!" once to precompile so that it does not appear in flamegraph and then rerun with @profview
TODO: get more details on the time "run_forward!" takes to compile 

```bash
module testprofile
ARGS = String["levelset_Butler.yml"]
include("../Flower.jl/examples/profile.jl")
end
```

```bash
using ProfileView
using Cthulhu
include("file.jl")
```
Then click on flamegraph, descend_clicked() and look for non-defined types
cf. [tutorial on YouTube](https://www.youtube.com/watch?v=pvduxLowpPY&t=1s) 
and [tutorial on GitHub](https://github.com/timholy/ProfileView.jl) 

## Precompiling issues
If fails to install a library like Cairo_jll 

Failed to precompile [...]
ERROR: LoadError: InitError: could not load library [...]lib/libgobject-2.0.so"
[...]/lib/libgobject-2.0.so: undefined symbol: g_dir_unref


do
```bash
LD_LIBRARY_PATH="" julia
```

If you want to load external libraries and julia throws an error when LD_LIBRARY_PATH is set (not to ""):

do
```bash
LD_LIBRARY_PATH="\$HOME/.julia/artifacts/20c009b8faa6b86ae9ecc8002f2772cd4724774d/lib/:/usr/lib/x86_64-linux-gnu" julia +1.10.5 --project=../Flower.jl --threads=1 --sysimage=../Flower.jl/sys_img.so ../Flower.jl/examples/main_current_folder.jl levelset_Butler.yml
```
with :/usr/lib/x86_64-linux-gnu the place where the external libraries are stored

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

### Ruche doc

[Ruche doc](https://mesocentre.pages.centralesupelec.fr/mesocenter_training/main.pdf)

[Julia on Ruche](https://gitlab-research.centralesupelec.fr/mesocentre-public/ruche_examples/-/tree/master/hello_julia
)

### Installation
```bash
curl -fsSL https://install.julialang.org | sh
```
For example, in in \$WORKDIR/flower (or \$HOME/flower if not possible but slower)
* Clone Flower.jl
*  Clone DDM.jl

<!-- what is added in project in \$HOME is not used in \$WORKDIR, so in \$WORKDIR do:
```bash
julia +1.10.4 –project=\$HOME/flower/Flower.jl 
import Pkg ; Pkg.add(package)
``` -->

### Run
In interactive session:
```bash
julia +1.10.4 --project=\$HOME/flower/Flower.jl --threads=1 \$HOME/flower/Flower.jl/examples/main_current_folder.jl \$HOME/flower/Flower.jl/examples/validation.yml
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


# Passing file in REPL mode

```bash
julia> module testyaml
       ARGS = String["half_circle_imposed_radius.yml"]
       include("../Flower.jl/examples/main_current_folder.jl")
       end
```



# Testing

You can run tests with:

```bash
julia +1.10.5 --project=../Flower.jl --threads=1 --trace-compile=dc_precompile.jl ../Flower.jl/test/runtests.jl
```


or in test : runtests.jl is used when doing:
```bash
]
test
```
Currently, the test half_circle_imposed_radius.yml has an error when launched with test mode instead of the classic command.

# Debugging/logging

## Things to check when the code gets stuck:

JULIA_DEBUG="all" julia... to display @debug prints

if stuck between Debug: After Numerical and after PDI init:

Check that the variables declared in PDI in the YAML are the same in Flower.

## Quick debugging with PDI
use PYCALL to do quick tests (print for example) on variables that are shared in PDI without modifying and recompiling the code

<!-- % \section{Installing Julia (attempt)}
% \url{https://github.com/JuliaLang/juliaup}
% juliaup add 1.9
% version 1.9
% in the same folder: git clone Flower and git clone DDM
% julia +1.9 --project=Flower.jl
% ] add DDM.jl
%
%
% 1.9.0
% using Flower
% instantiate
% resolve
% ctrl d
% julia +1.9.0 --project=Flower.jl
% using Flower
% update
% add StatsFuns@0.9.3 -->

## Debugging

AllocCheck.jl to check if the code allocates

```julia
check_allocs()
```



jet.jl

```julia
./julia ~/src/JETCLI/jet.jl stdlib/REPL
```


semgrep  

```julia
Base.infer_effects
```

```julia
code_typed(sum,(Vector{Float64,}))
```

