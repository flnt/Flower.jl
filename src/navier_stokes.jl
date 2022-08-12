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

function extrapolate_boundary!(cap, grid_u, dcapu, u, projection, inside)
    @inbounds @threads for II in inside
        if cap[II,1] < 1e-8 && dcapu[II,5] > 1e-10
            u_1, u_2 = interpolated_temperature(grid_u, projection[II].angle, projection[II].point1, projection[II].point2, u, II)
            if π/4 <= projection[II].angle <= 3π/4 || -π/4 >= projection[II].angle >= -3π/4
                u[II] = y_extrapolation(u_1, u_2, projection[II].point1, projection[II].point2, projection[II].mid_point)
            else
                u[II] = x_extrapolation(u_1, u_2, projection[II].point1, projection[II].point2, projection[II].mid_point)
            end
        end
    end
end

function set_stokes!(grid, geo, grid_u, geo_u, grid_v, geo_v, op, ph,
                    Hx, Hy, bcgpx, bcgpy, bcϕx, bcϕy, bcgϕx, bcgϕy, BC_p,
                    Lpm1, CUTpm1, Gxpm1, Gypm1, empty,
                    Hu, bcu, BC_u, Lum1, CUTum1, empty_u,
                    Hv, bcv, BC_v, Lvm1, CUTvm1, empty_v,
                    ns_vec, MIXED, MIXED_u, MIXED_v,
                    iMu, iMv, FRESH_u, FRESH_v,
                    Re, advection, ns_advection
    )
    @unpack Lp, CUTp, Lu, CUTu, Lv, CUTv, Dxu, CUTDx, Dyv, CUTDy, Gxp, CUTGxp, Gyp, CUTGyp, Gxϕ, CUTGxϕ, Gyϕ, CUTGyϕ, Cu, CUTCu, Cv, CUTCv, utp, vtp = op
    @unpack u, v, Du, Dv, tmp = ph

    iRe = 1 / Re

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

    Hx .= 0.
    Hy .= 0.
    for II in vcat(grid.ind.all_indices)
        Hx[II] = distance(grid.mid_point[II], geo.centroid[II], grid.dx[II], grid.dy[II])
        Hy[II] = distance(grid.mid_point[II], geo.centroid[II], grid.dx[II], grid.dy[II])
    end

    bcϕx .= 0.
    bcϕy .= 0.
    bcϕx, bcϕy = set_bc_bnds(neu, bcϕx, bcϕy, BC_p, grid.dx, grid.dy)

    bcgϕx .= 0.
    bcgϕy .= 0.
    bcgϕx, bcgϕy = set_bc_bnds(neu, bcgϕx, bcgϕy, BC_p, grid.dx, grid.dy, Hx, Hy)

    Hu .= 0.
    for II in vcat(grid_u.ind.b_left[1], grid_u.ind.b_bottom[1], grid_u.ind.b_right[1], grid_u.ind.b_top[1])
        Hu[II] = distance(grid_u.mid_point[II], geo_u.centroid[II], grid_u.dx[II], grid_u.dy[II])
    end

    Du .= grid_u.V
    bcu .= Du
    bcux, bcuy = set_bc_bnds(dir, bcu, Hu, BC_u)

    Hv .= 0.
    for II in vcat(grid_v.ind.b_left[1], grid_v.ind.b_bottom[1], grid_v.ind.b_right[1], grid_v.ind.b_top[1])
        Hv[II] = distance(grid_v.mid_point[II], geo_v.centroid[II], grid_v.dx[II], grid_v.dy[II])
    end

    Dv .= grid_v.V
    bcv .= Dv
    bcvx, bcvy = set_bc_bnds(dir, bcv, Hv, BC_v)

    if advection
        laplacian!(neu, Lp, CUTp, bcϕx, bcϕy, bcgϕx, bcgϕy, geo.dcap, grid.dx, grid.dy, grid.ny, BC_p, grid.ind.inside, empty, MIXED,
                ns_vec, grid.ind.b_left[1], grid.ind.b_bottom[1], grid.ind.b_right[1], grid.ind.b_top[1])

        laplacian!(dir, Lu, CUTu, bcux, bcuy, geo_u.dcap, grid_u.ny, BC_u, grid_u.ind.inside, empty_u,
                MIXED_u, grid_u.ind.b_left[1], grid_u.ind.b_bottom[1], grid_u.ind.b_right[1], grid_u.ind.b_top[1])

        laplacian!(dir, Lv, CUTv, bcvx, bcvy, geo_v.dcap, grid_v.ny, BC_v, grid_v.ind.inside, empty_v,
                MIXED_v, grid_v.ind.b_left[1], grid_v.ind.b_bottom[1], grid_v.ind.b_right[1], grid_v.ind.b_top[1])
    end

    uv_to_p!(utp, vtp, geo.dcap, grid.dx, grid.dy, grid.ny, grid.ind.all_indices)

    u_tmp = copy(u)
    v_tmp = copy(v)

    init_fresh_cells!(grid_u, u_tmp, geo_u.projection, FRESH_u)
    init_fresh_cells!(grid_v, v_tmp, geo_v.projection, FRESH_v)
    kill_dead_cells!(u_tmp, Lu, empty_u, MIXED_u, grid_u.ny)
    kill_dead_cells!(v_tmp, Lv, empty_v, MIXED_v, grid_v.ny)

    Δu = iMu * (Lu * vec(u_tmp) .+ CUTu)
    Δv = iMv * (Lv * vec(v_tmp) .+ CUTv)

    bcgpx .= Hx .* reshape(iRe .* (utp * Δu .+ vtp * Δv), (grid.ny, grid.nx))
    bcgpy .= bcgpx
    # bcgpx .= 0.
    # bcgpy .= 0.
    bcgpx, bcgpy = set_bc_bnds(neu, bcgpx, bcgpy, BC_p, grid.dx, grid.dy, Hx, Hy)

    if advection
        divergence!(dir, Dxu, Dyv, CUTDx, CUTDy, bcux, bcvy, geo.dcap, geo_u.dcap, geo_v.dcap, grid.ny, 
                grid.ind.all_indices)

        gradient!(neu, Gxp, Gyp, CUTGxp, CUTGyp, bcgpx, bcgpy, Dxu, Dyv, geo.cap, geo.dcap, geo_u.dcap, geo_v.dcap,
                grid_u.dx, grid_v.dy, grid.ny, BC_p, grid.ind.all_indices,
                grid_u.ind.b_left[1], grid_v.ind.b_bottom[1], grid_u.ind.b_right[1], grid_v.ind.b_top[1],
                grid.ind.b_left[1], grid.ind.b_bottom[1], grid.ind.b_right[1], grid.ind.b_top[1])

        gradient!(neu, Gxϕ, Gyϕ, CUTGxϕ, CUTGyϕ, bcgϕx, bcgϕy, Dxu, Dyv, geo.cap, geo.dcap, geo_u.dcap, geo_v.dcap,
                grid_u.dx, grid_v.dy, grid.ny, BC_p, grid.ind.all_indices,
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
                            FRESH, FRESH_u, FRESH_v, nullspace, ns_advection
    )
    @unpack Re, τ = num
    @unpack Lp, CUTp, Lu, CUTu, Lv, CUTv, Dxu, CUTDx, Dyv, CUTDy, Ap, Au, Av, Gxp, CUTGxp, Gyp, CUTGyp, Gxϕ, CUTGxϕ, Gyϕ, CUTGyϕ, Cu, CUTCu, Cv, CUTCv = op
    @unpack p, ϕ, Gxm1, Gym1, u, v, ucorr, vcorr, Du, Dv, tmp = ph

    iRe = 1/Re

    mat_op!(Ap, Lp, x->-x)
    Au .= Mu - iRe*τ*Lu
    Av .= Mv - iRe*τ*Lv

    if nullspace
        PETSc.destroy(ns)
    end
    ns = update_ksp_solver!(kspp, Ap, nullspace, ns_vec)

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
    init_fresh_cells!(grid_u, u, geo_u.projection, FRESH_u)
    init_fresh_cells!(grid_v, v, geo_v.projection, FRESH_v)
    kill_dead_cells!(p, Lp, EMPTY, MIXED, grid.ny)
    kill_dead_cells!(u, Lu, EMPTY_u, MIXED_u, grid_u.ny)
    kill_dead_cells!(v, Lv, EMPTY_v, MIXED_v, grid_v.ny)

    Gxm1 .= Mu * iMGx * (Gxp * vec(p) .+ CUTGxp)
    Gym1 .= Mv * iMGy * (Gyp * vec(p) .+ CUTGyp)

    Δu = Lu * vec(u) .+ CUTu
    Δv = Lv * vec(v) .+ CUTv

    Bδucorr = τ .* (iRe .* Δu .- Gxm1 .- Convu)
    Bδvcorr = τ .* (iRe .* Δv .- Gym1 .- Convv)
    # Bδucorr = τ .* (iRe .* Δu .- Convu)
    # Bδvcorr = τ .* (iRe .* Δv .- Convv)

    @inbounds @threads for II in EMPTY_u
        pII = lexicographic(II, grid_u.ny)
        Bδucorr[pII] = 0.
    end
    @inbounds @threads for II in EMPTY_v
        pII = lexicographic(II, grid_v.ny)
        Bδvcorr[pII] = 0.
    end

    # update_ksp_solver!(kspu, Au)
    # update_ksp_solver!(kspv, Av)
    # δucorr = reshape(kspu \ Bδucorr, (grid_u.ny, grid_u.nx))
    # δvcorr = reshape(kspv \ Bδvcorr, (grid_v.ny, grid_v.nx))
    δucorr = reshape(cg(Au, Bδucorr), (grid_u.ny, grid_u.nx))
    δvcorr = reshape(cg(Av, Bδvcorr), (grid_v.ny, grid_v.nx))
    ucorr .= δucorr .+ u
    vcorr .= δvcorr .+ v

    Duv = iMDx * (Dxu * vec(ucorr) .+ CUTDx) .+ iMDy * (Dyv * vec(vcorr) .+ CUTDy)
    Bϕ = - 1. ./ τ .* Mp * Duv .+ CUTp

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
    if nullspace
        PETSc.MatNullSpaceRemove!(ns, vecseq)
    end

    ϕ .= reshape(kspp \ vecseq, (grid.ny,grid.nx))
    PETSc.destroy(vecseq)

    Δϕ = reshape(iMp * (Lp * vec(ϕ) .+ CUTp), (grid.ny,grid.nx))
    p .+= ϕ #.- iRe.*reshape(Duv, (grid.ny,grid.nx))#τ.*Δϕ
    Gx = reshape(iMGx * (Gxϕ * vec(ϕ) .+ CUTGxϕ), (grid_u.ny,grid_u.nx))
    Gy = reshape(iMGy * (Gyϕ * vec(ϕ) .+ CUTGyϕ), (grid_v.ny,grid_v.nx))

    u .= ucorr .- τ .* Gx
    v .= vcorr .- τ .* Gy

    return nothing
end