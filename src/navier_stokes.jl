function no_slip_condition!(grid, grid_u, grid_v)
    interpolate_scalar!(grid, grid_u, grid_v, grid.V, grid_u.V, grid_v.V)

    normalx = cos.(grid_u.α)
    normaly = sin.(grid_v.α)

    grid_u.V .*= normalx
    grid_v.V .*= normaly

    replace!(grid_u.V, NaN=>0.0)
    replace!(grid_v.V, NaN=>0.0)

    return nothing
end

function set_stokes!(grid, geo, grid_u, geo_u, grid_v, geo_v, op, ph,
                    bcpx, bcpy, BC_p, Lpm1, CUTpm1, Gxpm1, Gypm1, empty,
                    Hu, bcu, BC_u, Lum1, CUTum1, empty_u,
                    Hv, bcv, BC_v, Lvm1, CUTvm1, empty_v,
                    ns_vec, MIXED, MIXED_u, MIXED_v, advection, ns_advection
    )
    @unpack Lp, CUTp, Lu, CUTu, Lv, CUTv, Dxu, CUTDx, Dyv, CUTDy, Gxp, Gyp, Cu, CUTCu, Cv, CUTCv = op
    @unpack u, v, Du, Dv = ph

    Gxpm1 .= Gxp
    Gypm1 .= Gyp
    Lpm1 .= Lp
    Lum1 .= Lu
    Lvm1 .= Lv
    CUTpm1 .= CUTp
    CUTum1 .= CUTu
    CUTvm1 .= CUTv

    empty_laplacian(grid, Lpm1, empty, MIXED)
    empty_laplacian(grid_u, Lum1, empty_u, MIXED_u)
    empty_laplacian(grid_v, Lvm1, empty_v, MIXED_v)

    bcpx .= 0.
    bcpy .= 0.
    bcpx, bcpy = set_bc_bnds(neu, bcpx, bcpy, BC_p, grid.dx, grid.dy)

    if advection
        laplacian!(neu, Lp, CUTp, bcpx, bcpy, geo.dcap, grid.dx, grid.dy, grid.ny, BC_p, grid.ind.inside, empty, MIXED,
                ns_vec, grid.ind.b_left[1], grid.ind.b_bottom[1], grid.ind.b_right[1], grid.ind.b_top[1])
        # laplacian!(neu, Lp, CUTp, bcpx, bcpy, geo.dcap, geo_u.dcap, geo_v.dcap, grid.ny, BC_p, grid.ind.inside, empty, MIXED,
        #         ns_vec, grid.ind.b_left[1], grid.ind.b_bottom[1], grid.ind.b_right[1], grid.ind.b_top[1])
    end

    Hu .= 0.
    for II in vcat(grid_u.ind.b_left[1], grid_u.ind.b_bottom[1], grid_u.ind.b_right[1], grid_u.ind.b_top[1])
        Hu[II] = distance(grid_u.mid_point[II], geo_u.centroid[II], grid_u.dx[II], grid_u.dy[II])
    end

    Du .= grid_u.V
    bcu .= Du
    bcux, bcuy = set_bc_bnds(dir, bcu, Hu, BC_u)

    if advection
        laplacian!(dir, Lu, CUTu, bcux, bcuy, geo_u.dcap, grid_u.ny, BC_u, grid_u.ind.inside, empty_u,
                MIXED_u, grid_u.ind.b_left[1], grid_u.ind.b_bottom[1], grid_u.ind.b_right[1], grid_u.ind.b_top[1])
    end

    Hv .= 0.
    for II in vcat(grid_v.ind.b_left[1], grid_v.ind.b_bottom[1], grid_v.ind.b_right[1], grid_v.ind.b_top[1])
        Hv[II] = distance(grid_v.mid_point[II], geo_v.centroid[II], grid_v.dx[II], grid_v.dy[II])
    end

    Dv .= grid_v.V
    bcv .= Dv
    bcvx, bcvy = set_bc_bnds(dir, bcv, Hv, BC_v)

    if advection
        laplacian!(dir, Lv, CUTv, bcvx, bcvy, geo_v.dcap, grid_v.ny, BC_v, grid_v.ind.inside, empty_v,
                MIXED_v, grid_v.ind.b_left[1], grid_v.ind.b_bottom[1], grid_v.ind.b_right[1], grid_v.ind.b_top[1])

        divergence!(dir, Dxu, Dyv, CUTDx, CUTDy, bcux, bcvy, geo.dcap, geo_u.dcap, geo_v.dcap, grid.ny, 
                grid.ind.all_indices)
        gradient!(neu, Gxp, Gyp, Dxu, Dyv, geo.dcap, geo_u.dcap, geo_v.dcap, grid.ny, BC_p,
                grid_u.ind.b_left[1], grid_v.ind.b_bottom[1], grid_u.ind.b_right[1], grid_v.ind.b_top[1],
                grid.ind.b_left[1], grid.ind.b_bottom[1], grid.ind.b_right[1], grid.ind.b_top[1])
        divergence_boundaries(dir, Dxu, Dyv, bcux, bcvy, geo.dcap, geo_u.dcap, geo_v.dcap, grid.ny, BC_u, BC_v,
                grid.ind.b_left[1], grid.ind.b_bottom[1], grid.ind.b_right[1], grid.ind.b_top[1])
    end

    if ns_advection
        bcuCu1_x, bcuCu1_y, bcuCu2_x, bcuCu2_y, bcvCu_x, bcvCu_y = set_bc_bnds(dir, GridFCx, Du, Dv, Hu, Hv, u, v, BC_u, BC_v)
        bcvCv1_x, bcvCv1_y, bcvCv2_x, bcvCv2_y, bcuCv_x, bcuCv_y = set_bc_bnds(dir, GridFCy, Dv, Du, Hv, Hu, v, u, BC_v, BC_u)

        vector_convection!(dir, GridFCx, Cu, CUTCu, u, v, bcuCu1_x, bcuCu1_y, bcuCu2_x, bcuCu2_y, bcvCu_x, bcvCu_y,
                geo.dcap, grid.ny, BC_u, grid_u.ind.inside,
                grid_u.ind.b_left[1], grid_u.ind.b_bottom[1], grid_u.ind.b_right[1], grid_u.ind.b_top[1])
        vector_convection!(dir, GridFCy, Cv, CUTCv, u, v, bcuCv_x, bcuCv_y, bcvCv1_x, bcvCv1_y, bcvCv2_x, bcvCv2_y,
                geo.dcap, grid.ny, BC_v, grid_v.ind.inside,
                grid_v.ind.b_left[1], grid_v.ind.b_bottom[1], grid_v.ind.b_right[1], grid_v.ind.b_top[1])
    end
    
    return nothing
end

function pressure_projection!(num, grid, geo, grid_u, geo_u, grid_v, geo_v, op, ph,
                            Lum1, Lvm1, CUTum1, CUTvm1, Cum1, Cvm1,
                            kspp, kspu, kspv, ns, ns_vec, Gxpm1, Gypm1,
                            Mp, iMp, Mu, Mv, iMGx, iMGy, iMDx, iMDy,
                            iMum1, iMvm1, iMGxm1, iMGym1,
                            MIXED, MIXED_u, MIXED_v, FULL, EMPTY, EMPTY_u, EMPTY_v,
                            FRESH, FRESH_u, FRESH_v, ns_advection
    )
    @unpack Re, τ = num
    @unpack Lp, CUTp, Lu, CUTu, Lv, CUTv, Dxu, CUTDx, Dyv, CUTDy, Ap, Au, Av, Gxp, Gyp, Cu, CUTCu, Cv, CUTCv = op
    @unpack p, ϕ, Gxm1, Gym1, u, v, Du, Dv = ph

    iRe = 1/Re

    mat_op!(Ap, Lp, x->-x)
    Au .= Mu - 0.5*iRe*τ*Lu
    Av .= Mv - 0.5*iRe*τ*Lv

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
        Convu = zeros(grid_u.nx*grid_u.ny)
        Convv = zeros(grid_v.nx*grid_v.ny)
    end

    Cum1 .= Cu * vec(u) .+ CUTCu
    Cvm1 .= Cv * vec(v) .+ CUTCv

    init_fresh_cells!(grid, p, geo.projection, FRESH)
    init_fresh_cells!(grid_u, u, grid_u.V, geo_u.projection, FRESH_u)
    init_fresh_cells!(grid_v, v, grid_v.V, geo_v.projection, FRESH_v)
    # init_fresh_cells!(Gxm1, projectionu, FRESH_u, grid_u.ny)
    # init_fresh_cells!(Gym1, projectionv, FRESH_v, grid_v.ny)
    kill_dead_cells!(p, Lp, EMPTY, MIXED, grid.ny)
    kill_dead_cells!(u, Lu, EMPTY_u, MIXED_u, grid_u.ny)
    kill_dead_cells!(v, Lv, EMPTY_v, MIXED_v, grid_v.ny)
    kill_dead_cells!(Gxm1, Lu, EMPTY_u, MIXED_u, grid_u.ny)
    kill_dead_cells!(Gym1, Lv, EMPTY_v, MIXED_v, grid_v.ny)

    Δum1 = Lu * vec(u) .+ CUTu
    Δvm1 = Lv * vec(v) .+ CUTv

    Δu = 0.5 .* (Δum1 .+ Δm1um1)
    Δv = 0.5 .* (Δvm1 .+ Δm1vm1)

    # Δu = Δum1
    # Δv = Δvm1

    Bδucorr = τ .* (iRe .* Δu .- Gxm1 .- Convu)
    Bδvcorr = τ .* (iRe .* Δv .- Gym1 .- Convv)

    @inbounds @threads for II in EMPTY_u
        pII = lexicographic(II, grid_u.ny)
        Bδucorr[pII] = 0.
    end
    @inbounds @threads for II in EMPTY_v
        pII = lexicographic(II, grid_v.ny)
        Bδvcorr[pII] = 0.
    end

    update_ksp_solver!(kspu, Au)
    update_ksp_solver!(kspv, Av)

    δucorr = reshape(kspu \ Bδucorr, (grid_u.ny, grid_u.nx))
    δvcorr = reshape(kspv \ Bδvcorr, (grid_v.ny, grid_v.nx))
    ucorr = δucorr .+ u
    vcorr = δvcorr .+ v

    Duv = Mp * (iMDx * (Dxu * vec(ucorr) .+ CUTDx) .+ iMDy * (Dyv * vec(vcorr) .+ CUTDy))
    Bϕ = -Duv./τ .+ CUTp

    sum_Bϕ = sum(Bϕ)
    sumMp = sum(Mp[diagind(Mp)])

    non_empty = vcat(FULL, MIXED)
    n_non_empty = length(non_empty)

    @inbounds @threads for II in non_empty
        pII = lexicographic(II, grid.ny)
        # @inbounds Bϕ[pII] -= sum_Bϕ/n_non_empty
        @inbounds Bϕ[pII] -= sum_Bϕ * Mp[pII,pII] / sumMp
    end
    vecseq = PETSc.VecSeq(Bϕ)
    PETSc.MatNullSpaceRemove!(ns, vecseq)

    ϕ .= reshape(kspp \ vecseq, (grid.ny,grid.nx))
    PETSc.destroy(vecseq)

    Δϕ = reshape(iMp * (Lp * vec(ϕ) .+ CUTp), (grid.ny,grid.nx))
    p .+= ϕ #.- iRe*0.5.*reshape(Duv, (grid.ny,grid.nx))
    Gx = reshape(iMGx * (Gxp * vec(ϕ)), (grid_u.ny,grid_u.nx))
    Gy = reshape(iMGy * (Gyp * vec(ϕ)), (grid_v.ny,grid_v.nx))

    u .= ucorr .- τ .* Gx
    v .= vcorr .- τ .* Gy

    return nothing
end