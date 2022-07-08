function save_field(path::String, num::Numerical, fwd::Forward)
    JLD.save(path, "u", fwd.uL, "v", fwd.vL, "gx", fwd.Gxm1L, "gy", fwd.Gym1L, "p", fwd.pL)
end

load_field(path::String) = load(path)

function force_coefficients(fwd, tmp, num, idx; A=1., p0=0., step=size(fwd.psave,1))
    D_p = 0.
    D_ν = 0.
    L_p = 0.
    L_ν = 0.

    p = fwd.psave[step,:,:]
    u = fwd.Uxsave[step,:,:]
    v = fwd.Uysave[step,:,:]

    strain_rate!(dir, tmp.E11, tmp.E12_x, tmp.E12_y, tmp.E22, tmp.LIQu, tmp.LIQv,
                 num.n, num.Δ, idx.all_indices, idx.inside)

    τ11 = reshape(2 ./ num.Re .* tmp.E11 * vec(u), (num.n, num.n))
    τ12 = reshape(2 ./ num.Re .* (tmp.E12_x * vec(u) .+ tmp.E12_y * vec(v)), (num.n, num.n))
    τ22 = reshape(2 ./ num.Re .* tmp.E22 * vec(v), (num.n, num.n))

    @inbounds for II in idx.inside
        # pressure forces
        D_p += -(p[II] - p0) * (tmp.LIQ[II,3] - tmp.LIQ[II,1]) * num.Δ
        L_p += -(p[II] - p0) * (tmp.LIQ[II,4] - tmp.LIQ[II,2]) * num.Δ

        #viscous forces (diagonal term)
        D_ν += τ11[II] * (tmp.LIQu[δx⁺(II),6] - tmp.LIQu[II,6]) * num.Δ
        L_ν += τ22[II] * (tmp.LIQv[δy⁺(II),7] - tmp.LIQv[II,7]) * num.Δ
    end
    @inbounds for II in idx.all_indices[1:end-1,2:end]
        #viscous forces (off-diagonal term)
        D_ν += τ12[II] * (tmp.LIQu[δy⁺(II),7] - tmp.LIQu[II,7]) * num.Δ
        L_ν += τ12[II] * (tmp.LIQv[δy⁺(II),6] - tmp.LIQv[δx⁻(δy⁺(II)),6]) * num.Δ
    end

    # println("$D_p $D_ν")

    D = D_p + D_ν
    L = L_p + L_ν

    Cd = 2D / A
    Cl = 2L / A

    return Cd, Cl
end