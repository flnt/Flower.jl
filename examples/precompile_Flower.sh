cp dc_precompile.jl ../Flower.jl/
cd ../Flower.jl/
LD_LIBRARY_PATH="" julia +1.10.5 --project=. src/PackageCompiler_setup_image.jl

# cp dc_precompile.jl $WORKDIR/flower/Flower.jl/
# cd $WORKDIR/flower/Flower.jl/
