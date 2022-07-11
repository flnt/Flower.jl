function no_slip_condition!(V, Vu, Vv, αu, αv, B, BT, n, inside, b_left, b_bottom, b_right, b_top)
    interpolate_scalar!(V, Vu, Vv, B, BT, n, inside, b_left, b_bottom, b_right, b_top)

    normalx = cos.(αu)
    normaly = sin.(αv)

    Vu .*= normalx
    Vv .*= normaly

    replace!(Vu, NaN=>0.0)
    replace!(Vv, NaN=>0.0)

    return nothing
end

function set_stokes!(bcpx, bcpy, BC_p, Lp, Lpm1, CUTp, CUTpm1, CAP, n, Δ, ns_advection,
                inside, empty, b_left, b_bottom, b_right, b_top,
                Hu, centroidu, mid_pointu, Du, Vu, bcu, BC_u, Lu, Lum1, CUTu, CUTum1, CAPu, 
                inside_u, empty_u, b_left_u, b_bottom_u, b_right_u, b_top_u,
                Hv, centroidv, mid_pointv, Dv, Vv, bcv, BC_v, Lv, Lvm1, CUTv, CUTvm1, CAPv, 
                inside_v, empty_v, b_left_v, b_bottom_v, b_right_v, b_top_v,
                Dxu, Dyv, CUTDx, CUTDy, all_indices,
                Gxp, Gyp, Gxpm1, Gypm1,
                CUTCu, CUTCv, Cu, Cv, u, v,
                current_i, ns_vec, MIXED_u, MIXED_v
    )
    Gxpm1 .= Gxp
    Gypm1 .= Gyp
    Lpm1 .= Lp
    Lum1 .= Lu
    Lvm1 .= Lv
    CUTpm1 .= CUTp
    CUTum1 .= CUTu
    CUTvm1 .= CUTv

    bcpx .= 0.
    bcpy .= 0.
    bcpx, bcpy = set_bc_bnds(neu, bcpx, bcpy, BC_p)

    laplacian!(neu, Lp, CUTp, bcpx, bcpy, CAP, n, Δ, BC_p, inside, empty,
                ns_vec, b_left[1], b_bottom[1], b_right[1], b_top[1])

    Hu .= 0.
    for II in vcat(b_left_u[1], b_bottom_u[1], b_right_u[1], b_top_u[1])
        Hu[II] = distance(mid_pointu[II], centroidu[II]) * Δ
    end

    Du .= Vu
    bcu .= Du
    bcux, bcuy = set_bc_bnds(dir, bcu, Hu, BC_u)

    laplacian!(dir, Lu, CUTu, bcux, bcuy, CAPu, n, Δ, BC_u, inside_u, empty_u,
                MIXED_u, b_left_u[1], b_bottom_u[1], b_right_u[1], b_top_u[1])

    Hv .= 0.
    for II in vcat(b_left_v[1], b_bottom_v[1], b_right_v[1], b_top_v[1])
        Hv[II] = distance(mid_pointv[II], centroidv[II]) * Δ
    end

    Dv .= Vv
    bcv .= Dv
    bcvx, bcvy = set_bc_bnds(dir, bcv, Hv, BC_v)

    laplacian!(dir, Lv, CUTv, bcvx, bcvy, CAPv, n+1, Δ, BC_v, inside_v, empty_v,
                MIXED_v, b_left_v[1], b_bottom_v[1], b_right_v[1], b_top_v[1])

    divergence!(dir, Dxu, Dyv, CUTDx, CUTDy, bcux, bcvy, CAP, CAPu, CAPv, n, Δ, all_indices)
    gradient!(neu, Gxp, Gyp, Dxu, Dyv, CAP, CAPu, CAPv, n, Δ, BC_p, b_left_u[1], b_bottom_v[1], b_right_u[1], b_top_v[1], b_left[1], b_bottom[1], b_right[1], b_top[1])
    divergence_boundaries(dir, Dxu, Dyv, bcux, bcvy, CAP, CAPu, CAPv, n, Δ, BC_u, BC_v, b_left[1], b_bottom[1], b_right[1], b_top[1])

    if ns_advection
        bcuCu1_x, bcuCu1_y, bcuCu2_x, bcuCu2_y, bcvCu_x, bcvCu_y = set_bc_bnds(dir, gfcx, Du, Dv, Hu, Hv, u, v, BC_u, BC_v)
        bcvCv1_x, bcvCv1_y, bcvCv2_x, bcvCv2_y, bcuCv_x, bcuCv_y = set_bc_bnds(dir, gfcy, Dv, Du, Hv, Hu, v, u, BC_v, BC_u)

        vector_convection!(dir, gfcx, Cu, CUTCu, u, v, bcuCu1_x, bcuCu1_y, bcuCu2_x, bcuCu2_y, bcvCu_x, bcvCu_y, CAP, n, Δ, BC_u, inside_u, b_left_u[1], b_bottom_u[1], b_right_u[1], b_top_u[1])
        vector_convection!(dir, gfcy, Cv, CUTCv, u, v, bcuCv_x, bcuCv_y, bcvCv1_x, bcvCv1_y, bcvCv2_x, bcvCv2_y, CAP, n, Δ, BC_v, inside_v, b_left_v[1], b_bottom_v[1], b_right_v[1], b_top_v[1])
    end
    
    return nothing
end

function pressure_projection!(p, ϕ, u, v, ns_advection,
                            Gxm1, Gym1, Ap, Au, Av, Lp, Lu, Lv, Lpm1, Lum1, Lvm1,
                            CUTp, CUTu, CUTv, CUTpm1, CUTum1, CUTvm1,
                            Cu, Cv, CUTCu, CUTCv, Cum1, Cvm1,
                            kspp, kspu, kspv, ns, ns_vec,
                            Gxp, Gyp, Gxpm1, Gypm1, Dxu, Dyv, CUTDx, CUTDy,
                            Mp, iMp, Mu, Mv, iMGx, iMGy, iMDx, iMDy,
                            iMu, iMv, iMpm1, iMum1, iMvm1, iMGxm1, iMGym1,
                            τ, iRe, Δ, n,
                            Vu, Vv, projection, projectionu, projectionv,
                            MIXED, MIXED_u, MIXED_v, FULL, EMPTY, EMPTY_u, EMPTY_v,
                            FRESH, FRESH_u, FRESH_v, WAS_MIXED_u, WAS_MIXED_v)
    mat_op!(Ap, Lp, x->-x)
    Au .= Mu .- 0.5.*iRe.*τ.*Lu
    Av .= Mv .- 0.5.*iRe.*τ.*Lv

    PETSc.destroy(ns)
    ns = update_ksp_solver!(kspp, Ap, true, ns_vec)

    Δm1um1 = Mu * (iMum1 * (Lum1 * vec(u) .+ CUTum1))
    Δm1vm1 = Mv * (iMvm1 * (Lvm1 * vec(v) .+ CUTvm1))

    Gxm1 .= Mu * iMGxm1 * Gxpm1 * vec(p)
    Gym1 .= Mv * iMGym1 * Gypm1 * vec(p)

    # Gxϕ = iMGxm1 * (Gxpm1 * vec(ϕ))
    # Gyϕ = iMGym1 * (Gypm1 * vec(ϕ))
    # Gxm1 .= Mu * (Gxm1 .+ Gxϕ .- iRe.*τ.*0.5 .* (iMum1 * (Lum1 * Gxϕ)))
    # Gym1 .= Mv * (Gym1 .+ Gyϕ .- iRe.*τ.*0.5 .* (iMvm1 * (Lvm1 * Gyϕ)))

    if ns_advection
        Convu = 1.5 .* (Cu * vec(u) .+ CUTCu) .- 0.5 .* Cum1
        Convv = 1.5 .* (Cv * vec(v) .+ CUTCv) .- 0.5 .* Cvm1
    else
        Convu = zeros(n*(n+1))
        Convv = zeros(n*(n+1))
    end

    Cum1 .= Cu * vec(u) .+ CUTCu
    Cvm1 .= Cv * vec(v) .+ CUTCv

    init_fresh_cells!(p, projection, FRESH)
    init_fresh_cells!(u, Vu, projectionu, FRESH_u)
    init_fresh_cells!(v, Vv, projectionv, FRESH_v)
    # init_fresh_cells!(Gxm1, projectionu, FRESH_u, n)
    # init_fresh_cells!(Gym1, projectionv, FRESH_v, n+1)
    kill_dead_cells!(p, Lp, EMPTY, MIXED, n)
    kill_dead_cells!(u, Lu, EMPTY_u, MIXED_u, n)
    kill_dead_cells!(v, Lv, EMPTY_v, MIXED_v, n+1)
    kill_dead_cells!(Gxm1, Lu, EMPTY_u, MIXED_u, n)
    kill_dead_cells!(Gym1, Lv, EMPTY_v, MIXED_v, n+1)

    Δum1 = Lu * vec(u) .+ CUTu
    Δvm1 = Lv * vec(v) .+ CUTv

    Δu = 0.5 .* (Δum1 .+ Δm1um1)
    Δv = 0.5 .* (Δvm1 .+ Δm1vm1)

    # Δu = Δum1
    # Δv = Δvm1

    Bδucorr = τ .* (iRe .* Δu .- Gxm1 .- Convu)
    Bδvcorr = τ .* (iRe .* Δv .- Gym1 .- Convv)

    @inbounds @threads for II in EMPTY_u
        pII = lexicographic(II, n)
        Bδucorr[pII] = 0.
    end
    @inbounds @threads for II in EMPTY_v
        pII = lexicographic(II, n+1)
        Bδvcorr[pII] = 0.
    end

    update_ksp_solver!(kspu, Au)
    update_ksp_solver!(kspv, Av)

    δucorr = reshape(kspu \ Bδucorr, (n, n+1))
    δvcorr = reshape(kspv \ Bδvcorr, (n+1, n))
    ucorr = δucorr .+ u
    vcorr = δvcorr .+ v

    Duv = Mp * (iMDx * (Dxu * vec(ucorr) .+ CUTDx) .+ iMDy * (Dyv * vec(vcorr) .+ CUTDy))
    Bϕ = -Duv./τ .+ CUTp

    sum_Bϕ = sum(Bϕ)
    sumMp = sum(Mp[diagind(Mp)])

    non_empty = vcat(FULL, MIXED)
    n_non_empty = length(non_empty)

    @inbounds @threads for II in non_empty
        pII = lexicographic(II, n)
        # @inbounds Bϕ[pII] -= sum_Bϕ/n_non_empty
        @inbounds Bϕ[pII] -= sum_Bϕ * Mp[pII,pII] / sumMp
    end
    vecseq = PETSc.VecSeq(Bϕ)
    PETSc.MatNullSpaceRemove!(ns, vecseq)

    ϕ .= reshape(kspp \ vecseq, (n,n))
    PETSc.destroy(vecseq)

    Δϕ = reshape(iMp * (Lp * vec(ϕ) .+ CUTp), (n,n))
    p .+= ϕ .- iRe.*τ.*0.5.*Δϕ
    Gx = reshape(iMGx * (Gxp * vec(ϕ)), (n,n+1))
    Gy = reshape(iMGy * (Gyp * vec(ϕ)), (n+1,n))

    u .= ucorr .- τ .* Gx
    v .= vcorr .- τ .* Gy

    return nothing
end