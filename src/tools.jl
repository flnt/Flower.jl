function save_field(path::String, num::Numerical, gp::Mesh, ph::Phase, fwdPh::ForwardPhase, fwd::Forward)
    JLD.save(path, "num", num, "gp", gp, "ph", ph, "fwdPh", fwdPh, "fwd", fwd)
end

function save_field(path::String, num::Numerical, gp::Mesh, ph::Phase)
    JLD.save(path, "num", num, "gp", gp, "ph", ph)
end

function load_phase!(data, ph::Phase)
    ph.u .= data["u"]
    ph.v .= data["v"]
    ph.p .= data["p"]
    ph.Gxm1 .= data["gx"]
    ph.Gym1 .= data["gy"]
end

function stretching(n::Int, dn0::Float64, dn1::Float64, ds::Float64, ws=12, we=12, maxs=0.04)
    ns = ceil(Int, ds÷dn0)
    ne = ns + log(dn1 / dn0) / log(1 + maxs)

    s = [maxs * 0.25 * (1 + erf(6 * (i - ns) / (ws))) * (1 - erf(6 * (i - ne) / we)) for i = 0:n-1]

    f_ = zeros(size(s))
    f_[1] = dn0
    for k = 2:length(f_)
        f_[k] = f_[k - 1] * (1 + s[k])
    end
    
    f = zeros(size(s))
    f[1] = 0.0
    for k = 2:length(f)
        f[k] = f[k - 1] + f_[k]
    end

    return f
end

function force_coefficients!(num, grid, grid_u, grid_v, op, fwd, ph; A=1., p0=0., step=size(fwd.psave,1), saveCoeffs = true)
    @unpack Re = num
    @unpack nx, ny, ind, LS = grid
    @unpack E11, E12_x, E12_y, E22 = op
    @unpack Cd, Cl = fwd

    D_p = 0.
    D_ν = 0.
    L_p = 0.
    L_ν = 0.

    p = ph.p[:,:]
    u = ph.u[:,:]
    v = ph.v[:,:]

    strain_rate!(dir, E11, E12_x, E12_y, E22, grid_u.LS[1].geoL.dcap, grid_v.LS[1].geoL.dcap,
                 ny, ind.all_indices, ind.inside)

    τ11 = reshape(2 ./ Re .* (E11 * vec(u)), (ny, nx))
    τ12 = reshape(2 ./ Re .* (E12_x * vec(u) .+ E12_y * vec(v)), (ny, nx))
    τ22 = reshape(2 ./ Re .* (E22 * vec(v)), (ny, nx))

    @inbounds for II in ind.inside
        # pressure forces
        D_p += -(p[II] - p0) * (LS[1].geoL.dcap[II,3] - LS[1].geoL.dcap[II,1])
        L_p += -(p[II] - p0) * (LS[1].geoL.dcap[II,4] - LS[1].geoL.dcap[II,2])

        # friction forces (diagonal terms)
        D_ν += τ11[II] * (grid_u.LS[1].geoL.dcap[δx⁺(II),6] - grid_u.LS[1].geoL.dcap[II,6])
        L_ν += τ22[II] * (grid_v.LS[1].geoL.dcap[δy⁺(II),7] - grid_v.LS[1].geoL.dcap[II,7])
    end
    @inbounds for II in ind.all_indices[1:end-1,2:end]
        # friction forces (off-diagonal terms)
        D_ν += τ12[II] * (grid_u.LS[1].geoL.dcap[δy⁺(II),7] - grid_u.LS[1].geoL.dcap[II,7])
        L_ν += τ12[II] * (grid_v.LS[1].geoL.dcap[δy⁺(II),6] - grid_v.LS[1].geoL.dcap[δx⁻(δy⁺(II)),6])
    end

    D = D_p + D_ν
    L = L_p + L_ν

    if saveCoeffs
        Cd[step] = 2D / A
        Cl[step] = 2L / A
    end

    return 2D / A, 2L / A, D, L
end

function vorticity(grid, ph)
    ω = similar(ph.p)
    ω .= (ph.u[:,1:end-1] .- ph.u[:,2:end]) ./ grid.dy
    ω .+= (ph.v[2:end,:] .- ph.v[1:end-1,:]) ./ grid.dx
    
    return ω
end