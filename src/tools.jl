function save_field(path::String, num::Numerical, fwd::Forward)
    JLD.save(path, "u", fwd.uL, "v", fwd.vL, "gx", fwd.Gxm1L, "gy", fwd.Gym1L, "p", fwd.pL)
end

load_field(path::String) = load(path)

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

function force_coefficients!(num, grid, grid_u, grid_v, op, fwd; A=1., p0=0., step=size(fwd.psave,1))
    @unpack Δ, Re = num
    @unpack nx, ny, ind, geoL = grid
    @unpack E11, E12_x, E12_y, E22 = op
    @unpack psave, Uxsave, Uysave, Cd, Cl = fwd

    D_p = 0.
    D_ν = 0.
    L_p = 0.
    L_ν = 0.

    p = psave[step,:,:]
    u = Uxsave[step,:,:]
    v = Uysave[step,:,:]

    strain_rate!(dir, E11, E12_x, E12_y, E22, grid_u.geoL.cap, grid_v.geoL.cap,
                 ny, Δ, ind.all_indices, ind.inside)

    τ11 = reshape(2 ./ Re .* E11 * vec(u), (ny, nx))
    τ12 = reshape(2 ./ Re .* (E12_x * vec(u) .+ E12_y * vec(v)), (ny, nx))
    τ22 = reshape(2 ./ Re .* E22 * vec(v), (ny, nx))

    @inbounds for II in ind.inside
        # pressure forces
        D_p += -(p[II] - p0) * (geoL.cap[II,3] - geoL.cap[II,1]) * Δ
        L_p += -(p[II] - p0) * (geoL.cap[II,4] - geoL.cap[II,2]) * Δ

        # friction forces (diagonal terms)
        D_ν += τ11[II] * (grid_u.geoL.cap[δx⁺(II),6] - grid_u.geoL.cap[II,6]) * Δ
        L_ν += τ22[II] * (grid_v.geoL.cap[δy⁺(II),7] - grid_v.geoL.cap[II,7]) * Δ
    end
    @inbounds for II in ind.all_indices[1:end-1,2:end]
        # friction forces (off-diagonal terms)
        D_ν += τ12[II] * (grid_u.geoL.cap[δy⁺(II),7] - grid_u.geoL.cap[II,7]) * Δ
        L_ν += τ12[II] * (grid_v.geoL.cap[δy⁺(II),6] - grid_v.geoL.cap[δx⁻(δy⁺(II)),6]) * Δ
    end

    D = D_p + D_ν
    L = L_p + L_ν

    Cd[step] = 2D / A
    Cl[step] = 2L / A

    return Cd[step], Cl[step]
end