function no_slip_condition!(grid, grid_u, grid_v)
    interpolate_scalar!(grid, grid_u, grid_v, grid.V, grid_u.V, grid_v.V)

    @inbounds grid_u.α[(abs.(grid_u.α) .- π) .< 1e-8] .= π
    @inbounds grid_v.α[(abs.(grid_v.α) .- π) .< 1e-8] .= π

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

function periodic_velocity!(grid_u, grid_v, Lu, Lv, Bu, Bv, periodic_x, periodic_y)
    if periodic_x
        @inbounds for II in grid_u.ind.b_left[1]
            pII = lexicographic(II, grid_u.ny)

            Lu[pII,pII] = -4.0
            
            Lu[pII+grid_u.ny,end-grid_v.ny+pII] = Lu[pII+grid_u.ny,pII]
            Lu[end-grid_v.ny+pII,pII+grid_u.ny] = Lu[end-grid_v.ny+pII,pII]

            Lu[pII,pII+grid_u.ny] = 0.0
            Lu[pII,end-grid_v.ny+pII] = 0.0
            
            Lu[pII+grid_u.ny,pII] = 0.0
            Lu[end-grid_v.ny+pII,pII] = 0.0


            if II[1] != 1
                Lu[pII-1,pII] = 0.0
                Lu[pII,pII-1] = 0.0
            end
            if II[1] != grid_u.nx
                Lu[pII+1,pII] = 0.0
                Lu[pII,pII+1] = 0.0
            end

            Bu[pII] = 0.0
        end
    end

    if periodic_y
        @inbounds for II in grid_v.ind.b_bottom[1]
            pII = lexicographic(II, grid_v.ny)

            Lv[pII,pII] = -4.0
            
            Lv[pII+1,pII+grid_v.ny-1] = Lv[pII+1,pII]
            Lv[pII+grid_v.ny-1,pII+1] = Lv[pII+grid_v.ny-1,pII]

            Lv[pII,pII+1] = 0.0
            Lv[pII,pII+grid_v.ny-1] = 0.0
            
            Lv[pII+1,pII] = 0.0
            Lv[pII+grid_v.ny-1,pII] = 0.0


            if II[2] != 1
                Lv[pII-grid_v.ny,pII] = 0.0
                Lv[pII,pII-grid_v.ny] = 0.0
            end
            if II[2] != grid_v.nx
                Lv[pII+grid_v.ny,pII] = 0.0
                Lv[pII,pII+grid_v.ny] = 0.0
            end

            Bv[pII] = 0.0
        end
    end

    return nothing
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
                            iMum1, iMvm1, iMGxm1, iMGym1, Mum1, Mvm1,
                            MIXED, MIXED_u, MIXED_v, FULL, EMPTY, EMPTY_u, EMPTY_v,
                            FRESH, FRESH_u, FRESH_v, nullspace, ns_advection, Ra,
                            periodic_x, periodic_y
    )
    @unpack Re, τ = num
    @unpack Lp, CUTp, Lu, CUTu, Lv, CUTv, Dxu, CUTDx, Dyv, CUTDy, Ap, Au, Av, Gxp, CUTGxp, Gyp, CUTGyp, Gxϕ, CUTGxϕ, Gyϕ, CUTGyϕ, Cu, CUTCu, Cv, CUTCv = op
    @unpack p, ϕ, Gxm1, Gym1, u, v, ucorr, vcorr, Du, Dv, tmp = ph

    iRe = 1/Re

    mat_op!(Ap, Lp, x->-x)

    if nullspace
        PETSc.destroy(ns)
    end
    PETSc.destroy(kspp)
    kspp, ns = init_ksp_solver(Ap, grid.nx, nullspace, ns_vec)
    # ns = update_ksp_solver!(kspp, Ap, nullspace, ns_vec)

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
    # Hardfix to ensure periodic borders are initialized
    if periodic_y
        @inbounds @threads for II in FRESH
            if II[1] == 1
                p[II] = p[δy⁺(II)]
            elseif II[1] == grid.ny
                p[II] = p[δy⁻(II)]
            end
        end
    end
    if periodic_x
        @inbounds @threads for II in FRESH
            if II[2] == 1
                p[II] = p[δx⁺(II)]
            elseif II[2] == grid.nx
                p[II] = p[δx⁻(II)]
            end
        end
    end
    init_fresh_cells!(grid_u, u, geo_u.projection, FRESH_u)
    init_fresh_cells!(grid_v, v, geo_v.projection, FRESH_v)
    kill_dead_cells!(p, Lp, EMPTY, MIXED, grid.ny)
    kill_dead_cells!(u, Lu, EMPTY_u, MIXED_u, grid_u.ny)
    kill_dead_cells!(v, Lv, EMPTY_v, MIXED_v, grid_v.ny)

    Gxm1 .= Mu * iMGx * (Gxp * vec(p) .+ CUTGxp)
    Gym1 .= Mv * iMGy * (Gyp * vec(p) .+ CUTGyp)
    # Gxm1 .= Mum1 * iMGxm1 * Gxpm1 * vec(p)
    # Gym1 .= Mvm1 * iMGym1 * Gypm1 * vec(p)

    Δu = Lu * vec(u) .+ CUTu
    Δv = Lv * vec(v) .+ CUTv

    Grav = Ra*vcat(0*mean(ph.T)*ones(grid_u.ny), vec(ph.T))

    @inbounds @threads for II in findall(geo_u.emptied)
        pII = lexicographic(II, grid_u.ny)
        Grav[pII] = 0.
    end

    Bucorr = Mum1 * vec(u) .+ τ .* (iRe.*CUTu .-Gxm1 .- Convu .+ Grav)
    Bvcorr = Mvm1 * vec(v) .+ τ .* (iRe.*CUTv .-Gym1 .- Convv)
    
    @inbounds @threads for II in EMPTY_u
        pII = lexicographic(II, grid_u.ny)
        Bucorr[pII] = 0.
    end
    @inbounds @threads for II in EMPTY_v
        pII = lexicographic(II, grid_v.ny)
        Bvcorr[pII] = 0.
    end

    Au .= Mu - iRe*τ*Lu
    Av .= Mv - iRe*τ*Lv

    ucorr .= reshape(cg(Au, Bucorr), (grid_u.ny, grid_u.nx))
    vcorr .= reshape(cg(Av, Bvcorr), (grid_v.ny, grid_v.nx))

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