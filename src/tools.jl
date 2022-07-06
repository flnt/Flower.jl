function save_field(path::String, num::Numerical, fwd::Forward)
    JLD.save(path, "u", fwd.uL, "v", fwd.vL, "gx", fwd.Gxm1L, "gy", fwd.Gym1L, "p", fwd.pL)
end

load_field(path::String) = load(path)