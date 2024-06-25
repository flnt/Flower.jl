@inline anisotropy(ϵ, m, θ, θ₀) = ϵ*(1 + 0.4*((8/3)*sin(0.5*m*(θ - θ₀))^4 - 1))

@inline function apply_curvature(num, grid, κ, bc, all_indices)
    @unpack ϵ_κ, ϵ_V = num
    @unpack V = grid

    @inbounds @threads for II in all_indices
        @inbounds bc[II] = bc[II] - ϵ_κ*κ[II] - ϵ_V*V[II]
    end
    return nothing
end

@inline function apply_anisotropy(num, grid, κ, bc, MIXED, sol_projection)
    @unpack ϵ_κ, ϵ_V, m, θ₀ = num
    @unpack V = grid

    @inbounds @threads for II in MIXED
        ϵ_c = anisotropy(ϵ_κ, m, sol_projection[II].angle, θ₀)
        ϵ_v = anisotropy(ϵ_V, m, sol_projection[II].angle, θ₀)
        @inbounds bc[II] = bc[II] - ϵ_c*κ[II] - ϵ_v*V[II]
    end
    return nothing
end

function set_heat_borders!(grid, a0, a1, b, BC_T, per_x, per_y)
    @unpack ny, ind = grid
    
    if per_y
        ind_x = 2:nx-1
        ind_y = 1:ny
    else
        ind_x = 1:nx
        ind_y = 2:ny-1
    end

    @inbounds a0[ind_y,1] .= BC_T.left.val
    if is_dirichlet(BC_T.left)
        @inbounds a1[ind_y,1] .= -1.
        @inbounds b[ind_y,1] .= 0.
    elseif is_neumann(BC_T.left)
        @inbounds a1[ind_y,1] .= 0.
        @inbounds b[ind_y,1] .= 1.
    elseif is_robin(BC_T.left)
        @inbounds a1[ind_y,1] .= -1.
        @inbounds b[ind_y,1] .= 1.
    elseif is_periodic(BC_T.left)
        nothing
    else
        @error ("Not implemented yet")
    end
    @inbounds a0[1,ind_x] .= BC_T.bottom.val
    if is_dirichlet(BC_T.bottom)
        @inbounds a1[1,ind_x] .= -1.
        @inbounds b[1,ind_x] .= 0.
    elseif is_neumann(BC_T.bottom)
        @inbounds a1[1,ind_x] .= 0.
        @inbounds b[1,ind_x] .= 1.
    elseif is_robin(BC_T.bottom)
        @inbounds a1[1,ind_x] .= -1.
        @inbounds b[1,ind_x] .= 1.
    elseif is_periodic(BC_T.bottom)
        nothing
    else
        @error ("Not implemented yet")
    end
    @inbounds a0[ind_y,end] .= BC_T.right.val
    if is_dirichlet(BC_T.right)
        @inbounds a1[ind_y,end] .= -1.
        @inbounds b[ind_y,end] .= 0.
    elseif is_neumann(BC_T.right)
        @inbounds a1[ind_y,end] .= 0.
        @inbounds b[ind_y,end] .= 1.
    elseif is_robin(BC_T.right)
        @inbounds a1[ind_y,end] .= -1.
        @inbounds b[ind_y,end] .= 1.
    elseif is_periodic(BC_T.right)
        nothing
    else
        @error ("Not implemented yet")
    end
    @inbounds a0[end,ind_x] .= BC_T.top.val
    if is_dirichlet(BC_T.top)
        @inbounds a1[end,ind_x] .= -1.
        @inbounds b[end,ind_x] .= 0.
    elseif is_neumann(BC_T.top)
        @inbounds a1[end,ind_x] .= 0.
        @inbounds b[end,ind_x] .= 1.
    elseif is_robin(BC_T.top)
        @inbounds a1[end,ind_x] .= -1.
        @inbounds b[end,ind_x] .= 1.
    elseif is_periodic(BC_T.top)
        nothing
    else
        @error ("Not implemented yet")
    end

    return nothing
end

function set_heat!(bc_type, num, grid, op, geo, ph, θd, BC_T, MIXED, projection,
    A, B,
    op_conv, grid_u, geo_u, grid_v, geo_v,
    periodic_x, periodic_y, convection, ls_advection, BC_int)
    @unpack τ, aniso = num
    @unpack nx, ny, dx, dy, ind  = grid
    @unpack all_indices, inside, b_left, b_bottom, b_right, b_top = ind
    @unpack Bx, By, BxT, ByT, Hx, Hy, HxT, HyT, M, iMx, iMy, χ = op
    @unpack CT, CUTCT = op_conv
    @unpack u, v, uD, vD = ph

    ni = nx * ny
    nb = 2 * nx + 2 * ny
    nt = 2 * ni + nb

    if is_dirichlet(bc_type)
        __a0 = bc_type.val
        __a1 = -1.0
        __b = 0.0
    elseif is_neumann(bc_type)
        __a0 = bc_type.val
        __a1 = 0.0
        __b = 1.0
    elseif is_robin(bc_type)
        __a0 = bc_type.val
        __a1 = -1.0
        __b = 1.0
    elseif is_stefan(bc_type)
        __a0 = θd
        __a1 = -1.0
        __b = 0.0
    elseif is_wall(bc_type)
        __a0 = bc_type.Tval
        __a1 = -1.0
        __b = 0.0
    else
        __a0 = bc_type.val
        __a1 = -1.0
        __b = 0.0
    end

    # Flags with BCs
    a0 = ones(grid) .* __a0
    if aniso
        apply_anisotropy(num, grid, grid.LS[1].κ, a0, MIXED, projection)
    else
        apply_curvature(num, grid, grid.LS[1].κ, a0, all_indices)
    end
    _a1 = ones(grid) .* __a1
    a1 = Diagonal(vec(_a1))
    _b = ones(grid) .* __b
    b = Diagonal(vec(_b))

    a0_b = zeros(nb)
    _a1_b = zeros(nb)
    _b_b = zeros(nb)
    set_borders!(grid, grid.LS[1].cl, grid.LS[1].u, a0_b, _a1_b, _b_b, BC_T, num.n_ext_cl)
    a1_b = Diagonal(vec(_a1_b))
    b_b = Diagonal(vec(_b_b))

    if convection
        HT = zeros(grid)
        @inbounds @threads for II in vcat(b_left[1], b_bottom[1], b_right[1], b_top[1])
            HT[II] = distance(grid.LS[1].mid_point[II], geo.centroid[II], dx[II], dy[II])
        end    
        bcTx, bcTy = set_bc_bnds(dir, a0, HT, BC_T)
    
        # bnds_u = [grid_u.ind.b_left[1], grid_u.ind.b_bottom[1], grid_u.ind.b_right[1], grid_u.ind.b_top[1]]
        # bnds_v = [grid_v.ind.b_left[1], grid_v.ind.b_bottom[1], grid_v.ind.b_right[1], grid_v.ind.b_top[1]]
        # Δu = [grid_u.dx[1,1], grid_u.dy[1,1], grid_u.dx[end,end], grid_u.dy[end,end]] .* 0.5
        # Δv = [grid_v.dx[1,1], grid_v.dy[1,1], grid_v.dx[end,end], grid_v.dy[end,end]] .* 0.5
        
        # Hu = zeros(grid_u)
        # for i in eachindex(bnds_u)
        #     for II in bnds_u[i]
        #         Hu[II] = Δu[i]
        #     end
        # end

        # Hv = zeros(grid_v)
        # for i in eachindex(bnds_v)
        #     for II in bnds_v[i]
        #         Hv[II] = Δv[i]
        #     end
        # end
    
        bcU = zeros(grid_u)
        bcU .= reshape(vec2(uD,grid_u), grid_u)
        bcU[1,:] .= vecb_B(uD, grid_u)
        bcU[end,:] .= vecb_T(uD, grid_u)
        bcU[:,1] .= vecb_L(uD, grid_u)
        bcU[:,end] .= vecb_R(uD, grid_u)
        
        bcV = zeros(grid_v)
        bcV .= reshape(vec2(vD,grid_v), grid_v)
        bcV[:,1] .= vecb_L(vD, grid_v)
        bcV[:,end] .= vecb_R(vD, grid_v)
        bcV[1,:] .= vecb_B(vD, grid_v)
        bcV[end,:] .= vecb_T(vD, grid_v)

        scalar_convection!(dir, CT, CUTCT, u, v, bcTx, bcTy, bcU, bcV, geo.dcap, ny, 
            BC_T, inside, b_left[1], b_bottom[1], b_right[1], b_top[1]
        )
    end

    if ls_advection
        update_all_ls_data(num, grid, grid_u, grid_v, BC_int, periodic_x, periodic_y, false)

        # Mass matrices
        M.diag .= vec(geo.dcap[:,:,5])
        Mx = zeros(ny,nx+1)
        for II in ind.all_indices
            Mx[II] = geo.dcap[II,8]
        end
        for II in ind.b_right[1]
            Mx[δx⁺(II)] = geo.dcap[II,10]
        end
        My = zeros(ny+1,nx)
        for II in ind.all_indices
            My[II] = geo.dcap[II,9]
        end
        for II in ind.b_top[1]
            My[δy⁺(II)] = geo.dcap[II,11]
        end
        iMx.diag .= 1. ./ (vec(Mx) .+ eps(0.01))
        iMy.diag .= 1. ./ (vec(My) .+ eps(0.01))

        # Discrete gradient and divergence operators
        divergence_B!(BxT, ByT, geo.dcap, ny, ind.all_indices)
        mat_assign!(Bx, sparse(-BxT'))
        mat_assign!(By, sparse(-ByT'))

        # Matrices for interior BCs
        for iLS in 1:num.nLS
            bc_matrix!(grid, Hx[iLS], Hy[iLS], geo.dcap, geo.dcap, ny, ind.all_indices)

            mat_assign_T!(HxT[iLS], sparse(Hx[iLS]'))
            mat_assign_T!(HyT[iLS], sparse(Hy[iLS]'))

            periodic_bcs!(grid, Bx, By, Hx[iLS], Hy[iLS], periodic_x, periodic_y)

            χx = (geo.dcap[:,:,3] .- geo.dcap[:,:,1]) .^ 2
            χy = (geo.dcap[:,:,4] .- geo.dcap[:,:,2]) .^ 2
            χ[iLS].diag .= sqrt.(vec(χx .+ χy))
        end
        mat_assign!(BxT, sparse(-Bx'))
        mat_assign!(ByT, sparse(-By'))

        # Matrices for borders BCs
        set_boundary_indicator!(grid, geo, geo, op)
        mass_matrix_borders!(ind, op.iMx_b, op.iMy_b, op.iMx_bd, op.iMy_bd, geo.dcap, ny)
        bc_matrix_borders!(grid, ind, op.Hx_b, op.Hy_b, geo.dcap)
        mat_assign_T!(op.HxT_b, sparse(op.Hx_b'))
        mat_assign_T!(op.HyT_b, sparse(op.Hy_b'))
        periodic_bcs_borders!(grid, op.Hx_b, op.Hy_b, periodic_x, periodic_y)
    end

    LT = BxT * iMx * Bx .+ ByT * iMy * By
    LD = BxT * iMx * Hx[1] .+ ByT * iMy * Hy[1]
    LD_b = BxT * op.iMx_b * op.Hx_b .+ ByT * op.iMy_b * op.Hy_b

    # Implicit part of heat equation
    A[1:ni,1:ni] = pad_crank_nicolson(M .- 0.5 .* τ .* LT, grid, τ)
    A[1:ni,ni+1:2*ni] = - 0.5 .* τ .* LD
    A[1:ni,end-nb+1:end] = - 0.5 .* τ .* LD_b

    # Interior BC
    A[ni+1:2*ni,1:ni] = b * (HxT[1] * iMx * Bx .+ HyT[1] * iMy * By)
    A[ni+1:2*ni,ni+1:2*ni] = pad(b * (HxT[1] * iMx * Hx[1] .+ HyT[1] * iMy * Hy[1]) .- χ[1] * a1)

    # Border BCs
    A[end-nb+1:end,1:ni] = b_b * (op.HxT_b * op.iMx_b' * Bx .+ op.HyT_b * op.iMy_b' * By)
    A[end-nb+1:end,ni+1:2*ni] = b_b * (op.HxT_b * op.iMx_b' * Hx[1] .+ op.HyT_b * op.iMy_b' * op.Hy[1])
    A[end-nb+1:end,end-nb+1:end] = pad(b_b * (op.HxT_b * op.iMx_bd * op.Hx_b .+ op.HyT_b * op.iMy_bd * op.Hy_b) .- op.χ_b * a1_b, 4.0)

    # Explicit part of heat equation
    B[1:ni,1:ni] = M .+ 0.5 .* τ .* LT .- τ .* CT
    B[1:ni,ni+1:2*ni] = 0.5 .* τ .* LD
    B[1:ni,end-nb+1:end] = 0.5 .* τ .* LD_b

    rhs = fnzeros(grid, num)
    if convection
        vec1(rhs,grid) .-= τ .* CUTCT
    end
    vec2(rhs,grid) .+= χ[1] * vec(a0)
    vecb(rhs,grid) .+= op.χ_b * vec(a0_b)

    return rhs
end
