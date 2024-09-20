module testprofile
using Profile
Profile.Allocs.clear()
ARGS = String["levelset_Butler.yml"]
@time Profile.Allocs.@profile sample_rate=1 include("../Flower.jl/examples/main_current_folder.jl")
end

#using PProf
#