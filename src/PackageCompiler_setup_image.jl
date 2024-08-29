# cf https://vsert.com/posts/precompile/

# julia --project=. --trace-compile=dc_precompile.jl src/precompile.jl

# This will generate ‘dc_precompile.jl’, a Julia file with all the statements to precompile. Create ‘src/setupimage.jl’

using Flower
using Pkg
using PackageCompiler
Pkg.activate(".") #TODO OK to comment this ? 
# otherwise :
# Activating new project at `~/flower/result_similar_to_Khalighi`
# ERROR: LoadError: could not find project at "/local/home/pr277828/flower/result_similar_to_Khalighi"
# )
create_sysimage(sysimage_path="sys_img.so", 
include_transitive_dependencies=false, 
cpu_target="generic", precompile_statements_file="dc_precompile.jl")

# and execute

# julia --project=. src/setupimage.jl

# This will create a ‘system image’, a binary with all the needed code pre-compiled, and using a generic architecture we’re sure this runs everywhere.

# Caution If you do not set cpu_target, it will default to your architecture, if you then run on an older generation CPU, you will crash.

# Now we can use this to speedup execution

# time julia --project=. --sysimage=sys_img.so -e 'using MyPackage; dostuff()'
