function save_field(path::String, num::Numerical, fwd::Forward; step=0)
    if step==0
        JLD.save(path, "u", fwd.uL, "v", fwd.vL, "p", fwd.pL)
    else
        JLD.save(path, "u", fwd.usave[step,:,:], "v", fwd.vsave[step,:,:], "p", fwd.psave[step,:,:])
    end
end

load_field(path::String) = load(path)