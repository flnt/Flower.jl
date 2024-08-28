module PrecompileTools_startup

using Flower
using PrecompileTools

@compile_workload begin
    # inside here, put a "toy example" of everything you want to be fast
    include("validation_test.jl")
end

end