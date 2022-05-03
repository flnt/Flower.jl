function set_stokes!(bcpx, bcpy, BC_p, Lp, CUTp, CAP, n, Δ,
                inside, empty, b_left, b_bottom, b_right, b_top,
                Hu, centroidu, mid_pointu, Du, bcu, BC_u, Lu, CUTu, CAPu, 
                inside_u, empty_u, b_left_u, b_bottom_u, b_right_u, b_top_u,
                Hv, centroidv, mid_pointv, Dv, bcv, BC_v, Lv, CUTv, CAPv, 
                inside_v, empty_v, b_left_v, b_bottom_v, b_right_v, b_top_v,
                Dxu, Dyv, CUTDx, CUTDy, all_indices,
                Gxp, Gyp)
    bcpx .= 0.
    bcpy .= 0.
    bcpx, bcpy = set_bc_bnds(neu, bcpx, bcpy, BC_p)

    laplacian!(neu, Lp, CUTp, bcpx, bcpy, CAP, n, Δ, BC_p, inside, empty,
                b_left[1], b_bottom[1], b_right[1], b_top[1])

    Hu .= 0.
    for II in vcat(b_left_u[1], b_bottom_u[1], b_right_u[1], b_top_u[1])
        Hu[II] = distance(mid_pointu[II], centroidu[II]) * Δ
    end

    Du .= 0. # set to V
    bcu .= Du
    bcux, bcuy = set_bc_bnds(dir, bcu, Hu, BC_u)

    laplacian!(dir, Lu, CUTu, bcux, bcuy, CAPu, n, Δ, BC_u, inside_u, empty_u,
                b_left_u[1], b_bottom_u[1], b_right_u[1], b_top_u[1])

    Hv .= 0.
    for II in vcat(b_left_v[1], b_bottom_v[1], b_right_v[1], b_top_v[1])
        Hv[II] = distance(mid_pointv[II], centroidv[II]) * Δ
    end

    Dv .= 0. # set to V
    bcv .= Dv
    bcvx, bcvy = set_bc_bnds(dir, bcv, Hv, BC_v)

    laplacian!(dir, Lv, CUTv, bcvx, bcvy, CAPv, n+1, Δ, BC_v, inside_v, empty_v,
                b_left_v[1], b_bottom_v[1], b_right_v[1], b_top_v[1])

    divergence!(dir, Dxu, Dyv, CUTDx, CUTDy, bcux, bcvy, CAP, CAPu, CAPv, n, Δ, all_indices)
    gradient!(neu, Gxp, Gyp, Dxu, Dyv, CAP, CAPu, CAPv, n, Δ, BC_p, b_left_u[1], b_bottom_v[1], b_right_u[1], b_top_v[1], b_left[1], b_bottom[1], b_right[1], b_top[1])
    divergence_boundaries(dir, Dxu, Dyv, bcux, bcvy, CAP, CAPu, CAPv, n, Δ, BC_u, BC_v, b_left[1], b_bottom[1], b_right[1], b_top[1])
    
    return nothing
end

function pressure_projection!(p, u, v, Ap, Au, Av, Lp, Lu, Lv,
                            CUTp, CUTu, CUTv,
                            kspp, kspu, kspv, ns, ns_vec,
                            Gxp, Gyp, Dxu, Dyv, CUTDx, CUTDy,
                            Mp, iMp, Mu, Mv, iMGx, iMGy, iMDx, iMDy,
                            τ, iRe, n, b_right)
    Ap .= -τ .* Lp
    Au .= Mu .- 0.5.*iRe.*τ.*Lu
    Av .= Mv .- 0.5.*iRe.*τ.*Lv

    PETSc.destroy(ns)
    ns = update_ksp_solver!(kspp, Ap, true, ns_vec)
    update_ksp_solver!(kspu, Au)
    update_ksp_solver!(kspv, Av)

    Lu = iRe .* (Lu * vec(u) .+ CUTu)
    Lv = iRe .* (Lv * vec(v) .+ CUTv)

    Gx = Mu * iMGx * Gxp * vec(p)
    Gy = Mv * iMGy * Gyp * vec(p)

    Bδucorr = τ .* (Lu .- Gx)
    Bδvcorr = τ .* (Lv .- Gy)

    δucorr = reshape(kspu \ Bδucorr, (n, n+1))
    δvcorr = reshape(kspv \ Bδvcorr, (n+1, n))
    ucorr = δucorr .+ u
    vcorr = δvcorr .+ v

    Duv = Mp * (iMDx * (Dxu * vec(ucorr) .+ CUTDx) .+ iMDy * (Dyv * vec(vcorr) .+ CUTDy))
    Bϕ = -Duv .+ τ .* CUTp
    sum_Bϕ = sum(Bϕ)
    @inbounds @threads for II in b_right[1]
        pII = lexicographic(II, n)
        @inbounds Bϕ[pII] -= sum_Bϕ/n
    end
    vecseq = PETSc.VecSeq(Bϕ)
    PETSc.MatNullSpaceRemove!(ns, vecseq)

    ϕ = reshape(kspp \ vecseq, (n,n))
    PETSc.destroy(vecseq)

    p .+= ϕ .- iRe.*τ.*0.5.*reshape(iMp * (Lp * vec(ϕ) .+ CUTp), (n,n))
    u .= ucorr .- τ.*reshape(iMGx * Gxp * vec(ϕ), (n,n+1))
    v .= vcorr .- τ.*reshape(iMGy * Gyp * vec(ϕ), (n+1,n))

    return nothing
end