function set_borders!(grid, a0, a1, b, BC, per_x, per_y)
    @unpack nx, ny, ind = grid

    if per_y
        ind_x = 2:nx-1
        ind_y = 1:ny
    else
        ind_x = 1:nx
        ind_y = 2:ny-1
    end
    
    idx = CartesianIndices((ind_y,1))
    @inbounds a0[idx] .= (BC.left.val .* ones(grid))[idx]
    if is_dirichlet(BC.left)
        @inbounds a1[idx] .= -1.0
        @inbounds b[idx] .= 0.0
    elseif is_neumann(BC.left)
        @inbounds a1[idx] .= 0.0
        @inbounds b[idx] .= 1.0
    elseif is_robin(BC.left)
        @inbounds a1[idx] .= -1.0
        @inbounds b[idx] .= 1.0
    elseif is_periodic(BC.left)
        nothing
    elseif is_navier(BC.left)
        @inbounds a1[idx] .= -1.0
        @inbounds b[idx] .= BC.left.λ
    elseif is_navier_cl(BC.left)
        @inbounds a1[intersect(idx, ind.cl)] .= -1.0
        @inbounds a1[symdiff(idx, ind.cl)] .= -1.0
        @inbounds b[intersect(idx, ind.cl)] .= BC.left.λ
        @inbounds b[symdiff(idx, ind.cl)] .= 0.0
    else
        @error ("Not implemented yet")
    end

    idx = CartesianIndices((1,ind_x))
    @inbounds a0[idx] .= (BC.bottom.val .* ones(grid))[idx]
    if is_dirichlet(BC.bottom)
        @inbounds a1[idx] .= -1.0
        @inbounds b[idx] .= 0.0
    elseif is_neumann(BC.bottom)
        @inbounds a1[idx] .= 0.0
        @inbounds b[idx] .= 1.0
    elseif is_robin(BC.bottom)
        @inbounds a1[idx] .= -1.0
        @inbounds b[idx] .= 1.0
    elseif is_periodic(BC.bottom)
        nothing
    elseif is_navier(BC.bottom)
        @inbounds a1[idx] .= -1.0
        @inbounds b[idx] .= BC.bottom.λ
    elseif is_navier_cl(BC.bottom)
        @inbounds a1[intersect(idx, ind.cl)] .= -1.0
        @inbounds a1[symdiff(idx, ind.cl)] .= -1.0
        @inbounds b[intersect(idx, ind.cl)] .= BC.bottom.λ
        @inbounds b[symdiff(idx, ind.cl)] .= 0.0
    else
        @error ("Not implemented yet")
    end

    idx = CartesianIndices((ind_y,nx:nx))
    @inbounds a0[idx] .= (BC.right.val .* ones(grid))[idx]
    if is_dirichlet(BC.right)
        @inbounds a1[idx] .= -1.0
        @inbounds b[idx] .= 0.0
    elseif is_neumann(BC.right)
        @inbounds a1[idx] .= 0.0
        @inbounds b[idx] .= 1.0
    elseif is_robin(BC.right)
        @inbounds a1[idx] .= -1.0
        @inbounds b[idx] .= 1.0
    elseif is_periodic(BC.right)
        nothing
    elseif is_navier(BC.right)
        @inbounds a1[idx] .= -1.0
        @inbounds b[idx] .= BC.right.λ
    elseif is_navier_cl(BC.right)
        @inbounds a1[intersect(idx, ind.cl)] .= -1.0
        @inbounds a1[symdiff(idx, ind.cl)] .= -1.0
        @inbounds b[intersect(idx, ind.cl)] .= BC.right.λ
        @inbounds b[symdiff(idx, ind.cl)] .= 0.0
    else
        @error ("Not implemented yet")
    end

    idx = CartesianIndices((ny:ny,ind_x))
    @inbounds a0[idx] .= (BC.top.val .* ones(grid))[idx]
    if is_dirichlet(BC.top)
        @inbounds a1[idx] .= -1.0
        @inbounds b[idx] .= 0.0
    elseif is_neumann(BC.top)
        @inbounds a1[idx] .= 0.0
        @inbounds b[idx] .= 1.0
    elseif is_robin(BC.top)
        @inbounds a1[idx] .= -1.0
        @inbounds b[idx] .= 1.0
    elseif is_periodic(BC.top)
        nothing
    elseif is_navier(BC.top)
        @inbounds a1[idx] .= -1.0
        @inbounds b[idx] .= BC.top.λ
    elseif is_navier_cl(BC.top)
        @inbounds a1[intersect(idx, ind.cl)] .= -1.0
        @inbounds a1[symdiff(idx, ind.cl)] .= -1.0
        @inbounds b[intersect(idx, ind.cl)] .= BC.top.λ
        @inbounds b[symdiff(idx, ind.cl)] .= 0.0
    else
        @error ("Not implemented yet")
    end

    return nothing
end

function set_borders!(grid, a0, a1, b, BC, n_ext)
    @unpack nx, ny, x, y, dx, dy, u, ind = grid

    idx = 1:ny
    @inbounds a0[idx] .= BC.left.val .* ones(ny)
    if is_dirichlet(BC.left)
        @inbounds a1[idx] .= -1.0
        @inbounds b[idx] .= 0.0
    elseif is_neumann(BC.left)
        @inbounds a1[idx] .= 0.0
        @inbounds b[idx] .= 1.0
    elseif is_robin(BC.left)
        @inbounds a1[idx] .= -1.0
        @inbounds b[idx] .= 1.0
    elseif is_periodic(BC.left)
        nothing
    elseif is_navier(BC.left)
        @inbounds a1[idx] .= -1.0
        @inbounds b[idx] .= BC.left.λ
    elseif is_navier_cl(BC.left)
        idx_cl = intersect(CartesianIndices((idx,1)), ind.cl)
        idx_no = symdiff(CartesianIndices((idx,1)), ind.cl)
        ϵb = zeros(grid)
        ϵb[idx_cl] .= n_ext .* dx[idx_cl]

        @inbounds @threads for II in idx_cl
            @inbounds a1[II[1]] = -1.0
            @inbounds b[II[1]] = BC.left.λ * bell_function2(u[II], ϵb[II])
        end
        @inbounds @threads for II in idx_no
            @inbounds a1[II[1]] = -1.0
            @inbounds b[II[1]] = 0.0
        end
    elseif is_gnbc(BC.left)
        # idx_cl = intersect(CartesianIndices((idx,1)), ind.cl)
        # idx_no = symdiff(CartesianIndices((idx,1)), ind.cl)

        # bell = bell_function(grid, BC.ϵ)
        # @inbounds a0[idx] .+= bell[idx] .* BC.σ ./ BC.μ .* (cos(θd) .- cos(BC.θe))
        # @inbounds @threads for II in idx_cl
        #     @inbounds a1[II[1]] = -1.0
        #     @inbounds b[II[1]] = BC.left.λ
        # end
        # @inbounds @threads for II in idx_no
        #     @inbounds a1[II[1]] = -1.0
        #     @inbounds b[II[1]] = 0.0
        # end
        nothing
    else
        @error ("Not implemented yet")
    end

    _idx = ny+1:ny+nx
    idx = 1:nx
    @inbounds a0[_idx] .= BC.bottom.val .* ones(nx)
    if is_dirichlet(BC.bottom)
        @inbounds a1[_idx] .= -1.0
        @inbounds b[_idx] .= 0.0
    elseif is_neumann(BC.bottom)
        @inbounds a1[_idx] .= 0.0
        @inbounds b[_idx] .= 1.0
    elseif is_robin(BC.bottom)
        @inbounds a1[_idx] .= -1.0
        @inbounds b[_idx] .= 1.0
    elseif is_periodic(BC.bottom)
        nothing
    elseif is_navier(BC.bottom)
        @inbounds a1[_idx] .= -1.0
        @inbounds b[_idx] .= BC.bottom.λ
    elseif is_navier_cl(BC.bottom)
        idx_cl = intersect(CartesianIndices((1,idx)), ind.cl)
        idx_no = symdiff(CartesianIndices((1,idx)), ind.cl)
        ϵb = zeros(grid)
        ϵb[idx_cl] .= n_ext .* dy[idx_cl]

        @inbounds @threads for II in idx_cl
            @inbounds a1[II[2]+ny] = -1.0
            @inbounds b[II[2]+ny] = BC.bottom.λ * bell_function2(u[II], ϵb[II])
        end
        @inbounds @threads for II in idx_no
            @inbounds a1[II[2]+ny] = -1.0
            @inbounds b[II[2]+ny] = 0.0
        end
    elseif is_gnbc(BC.bottom)
        # @inbounds a0[idx] .+= bell[idx] .* BC.σ ./ BC.μ .* (cos(θd) .- cos(BC.θe))
        # @inbounds a1[intersect(idx, ind.cl)] .= -1.0
        # @inbounds a1[symdiff(idx, ind.cl)] .= -1.0
        # @inbounds b[intersect(idx, ind.cl)] .= BC.bottom.λ
        # @inbounds b[symdiff(idx, ind.cl)] .= 0.0
        nothing
    else
        @error ("Not implemented yet")
    end

    _idx = ny+nx+1:2*ny+nx
    idx = 1:ny
    @inbounds a0[_idx] .= BC.right.val .* ones(ny)
    if is_dirichlet(BC.right)
        @inbounds a1[_idx] .= -1.0
        @inbounds b[_idx] .= 0.0
    elseif is_neumann(BC.right)
        @inbounds a1[_idx] .= 0.0
        @inbounds b[_idx] .= 1.0
    elseif is_robin(BC.right)
        @inbounds a1[_idx] .= -1.0
        @inbounds b[_idx] .= 1.0
    elseif is_periodic(BC.right)
        nothing
    elseif is_navier(BC.right)
        @inbounds a1[_idx] .= -1.0
        @inbounds b[_idx] .= BC.right.λ
    elseif is_navier_cl(BC.right)
        idx_cl = intersect(CartesianIndices((idx,nx)), ind.cl)
        idx_no = symdiff(CartesianIndices((idx,nx)), ind.cl)
        ϵb = zeros(grid)
        ϵb[idx_cl] .= n_ext .* dx[idx_cl]

        @inbounds @threads for II in idx_cl
            @inbounds a1[II[1]+ny+nx] = -1.0
            @inbounds b[II[1]+ny+nx] = BC.right.λ * bell_function2(u[II], ϵb[II])
        end
        @inbounds @threads for II in idx_no
            @inbounds a1[II[1]+ny+nx] = -1.0
            @inbounds b[II[1]+ny+nx] = 0.0
        end
    elseif is_gnbc(BC.right)
        # @inbounds a0[idx] .+= bell[idx] .* BC.σ ./ BC.μ .* (cos(θd) .- cos(BC.θe))
        # @inbounds a1[intersect(idx, ind.cl)] .= -1.0
        # @inbounds a1[symdiff(idx, ind.cl)] .= -1.0
        # @inbounds b[intersect(idx, ind.cl)] .= BC.right.λ
        # @inbounds b[symdiff(idx, ind.cl)] .= 0.0
        nothing
    else
        @error ("Not implemented yet")
    end

    _idx = 2*ny+nx+1:2*ny+2*nx
    idx = 1:nx
    @inbounds a0[_idx] .= BC.top.val .* ones(nx)
    if is_dirichlet(BC.top)
        @inbounds a1[_idx] .= -1.0
        @inbounds b[_idx] .= 0.0
    elseif is_neumann(BC.top)
        @inbounds a1[_idx] .= 0.0
        @inbounds b[_idx] .= 1.0
    elseif is_robin(BC.top)
        @inbounds a1[_idx] .= -1.0
        @inbounds b[_idx] .= 1.0
    elseif is_periodic(BC.top)
        nothing
    elseif is_navier(BC.top)
        @inbounds a1[_idx] .= -1.0
        @inbounds b[_idx] .= BC.top.λ
    elseif is_navier_cl(BC.top)
        idx_cl = intersect(CartesianIndices((1,idx)), ind.cl)
        idx_no = symdiff(CartesianIndices((1,idx)), ind.cl)
        ϵb = zeros(grid)
        ϵb[idx_cl] .= n_ext .* dy[idx_cl]

        @inbounds @threads for II in idx_cl
            @inbounds a1[II[2]+2*ny+nx] = -1.0
            @inbounds b[II[2]+2*ny+nx] = BC.top.λ * bell_function2(u[II], ϵb[II])
        end
        @inbounds @threads for II in idx_no
            @inbounds a1[II[2]+2*ny+nx] = -1.0
            @inbounds b[II[2]+2*ny+nx] = 0.0
        end
    elseif is_gnbc(BC.top)
        # @inbounds a0[idx] .+= bell[idx] .* BC.σ ./ BC.μ .* (cos(θd) .- cos(BC.θe))
        # @inbounds a1[intersect(idx, ind.cl)] .= -1.0
        # @inbounds a1[symdiff(idx, ind.cl)] .= -1.0
        # @inbounds b[intersect(idx, ind.cl)] .= BC.top.λ
        # @inbounds b[symdiff(idx, ind.cl)] .= 0.0
        nothing
    else
        @error ("Not implemented yet")
    end

    return nothing
end

function update_dirichlet_field!(grid, bv, v, BC)
    @unpack nx, ny = grid
    tmp = zeros(size(v))

    val = (BC.left.val .* ones(ny))[ind_y]
    if is_dirichlet(BC.left)
        @inbounds tmp[ind_y,1] .= val
    elseif is_neumann(BC.left)
        @inbounds tmp[ind_y,1] .= v[ind_y,1] .- val .* 0.5 .* grid.dx[ind_y,1]
    end

    val = (BC.bottom.val .* ones(nx))[ind_x]
    if is_dirichlet(BC.bottom)
        @inbounds tmp[1,ind_x] .= val
    elseif is_neumann(BC.bottom)
        @inbounds tmp[1,ind_x] .= v[1,ind_x] .- val .* 0.5 .* grid.dy[1,ind_x]
    end

    val = (BC.right.val .* ones(ny))[ind_y]
    if is_dirichlet(BC.right)
        @inbounds tmp[ind_y,end] .= val
    elseif is_neumann(BC.right)
        @inbounds tmp[ind_y,end] .= val .* 0.5 .* grid.dx[ind_y,end] .+ v[ind_y,end]
    end

    val = (BC.top.val .* ones(nx))[ind_x]
    if is_dirichlet(BC.top)
        @inbounds tmp[end,ind_x] .= val
    elseif is_neumann(BC.top)
        @inbounds tmp[end,ind_x] .= val .* 0.5 .* grid.dy[end,ind_x]  .+ v[end,ind_x]
    end

    bv[grid.ny*grid.nx+1:end] .= vec(tmp)

    return nothing
end

function set_cutcell_matrices!(grid, geo, geo_p, opC, periodic_x, periodic_y)
    @unpack nx, ny, ind = grid
    @unpack AxT, AyT, Bx, By, BxT, ByT, Hx, Hy, HxT, HyT, M, iMx, iMy, χ = opC

    χx = (geo.dcap[:,:,3] .- geo.dcap[:,:,1]) .^ 2
    χy = (geo.dcap[:,:,4] .- geo.dcap[:,:,2]) .^ 2
    χ.diag .= sqrt.(vec(χx .+ χy)) 

    M.diag .= vec(geo.dcap[:,:,5])

    Mx = zeros(ny,nx+1)
    for II in ind.all_indices
        Mx[II] = geo.dcap[II,8]
    end
    for II in ind.b_right[1]
        Mx[δx⁺(II)] = geo.dcap[II,10]
    end
    iMx.diag .= 1. ./ (vec(Mx) .+ eps(0.01))

    My = zeros(ny+1,nx)
    for II in ind.all_indices
        My[II] = geo.dcap[II,9]
    end
    for II in ind.b_top[1]
        My[δy⁺(II)] = geo.dcap[II,11]
    end
    iMy.diag .= 1. ./ (vec(My) .+ eps(0.01))

    # Discrete gradient and divergence operators
    divergence_A!(grid, AxT, AyT, geo.dcap, ny, ind.all_indices, periodic_x, periodic_y)
    divergence_B!(BxT, ByT, geo.dcap, ny, ind.all_indices)

    mat_assign!(Bx, sparse(-BxT'))
    mat_assign!(By, sparse(-ByT'))

    # Matrices for BCs
    bc_matrix!(grid, Hx, Hy, geo.dcap, geo_p.dcap, ny, ind.all_indices)

    mat_assign_T!(HxT, sparse(Hx'))
    mat_assign_T!(HyT, sparse(Hy'))

    periodic_bcs!(grid, Bx, By, Hx, Hy, periodic_x, periodic_y)

    mat_assign!(BxT, sparse(-Bx'))
    mat_assign!(ByT, sparse(-By'))

    return nothing
end

function set_other_cutcell_matrices!(
    grid, geo, geo_u, geo_v,
    opC_p, opC_u, opC_v,
    periodic_x, periodic_y
    )
    @unpack nx, ny, ind = grid
    @unpack Bx, By, Gx, Gy = opC_p

    bc_matrix!(grid, opC_u.Gx, opC_v.Gy, geo.dcap, geo_u.dcap, geo_v.dcap, ny, ind.all_indices)

    mat_assign_T!(Gx, sparse(opC_u.Gx'))
    mat_assign_T!(Gy, sparse(opC_v.Gy'))

    periodic_bcs!(grid, Bx, By, opC_u.Gx, opC_v.Gy, periodic_x, periodic_y)
    periodic_bcs_R!(grid, opC_u.Rx, opC_v.Ry, periodic_x, periodic_y)

    return nothing
end

function set_boundary_indicator!(grid::Mesh{GridCC,T,N}, geo, geo_p, opC) where {T,N}
    @unpack nx, ny, ind = grid
    @inbounds @threads for i in 1:ny
        II = ind.b_left[1][i]
        opC.χ_b[i, i] = geo.dcap[II,1]
    end
    @inbounds @threads for i in 1:nx
        II = ind.b_bottom[1][i]
        opC.χ_b[i+ny, i+ny] = geo.dcap[II,2]
    end
    @inbounds @threads for i in 1:ny
        II = ind.b_right[1][i]
        opC.χ_b[i+ny+nx, i+ny+nx] = geo.dcap[II,3]
    end
    @inbounds @threads for i in 1:nx
        II = ind.b_top[1][i]
        opC.χ_b[i+2*ny+nx, i+2*ny+nx] = geo.dcap[II,4]
    end

    return nothing
end

function set_boundary_indicator!(grid::Mesh{GridFCx,T,N}, geo, geo_p, opC) where {T,N}
    @unpack nx, ny, ind = grid
    @inbounds @threads for i in 1:ny
        II = ind.b_left[1][i]
        opC.χ_b[i, i] = geo_p.dcap[II,1]
    end
    @inbounds @threads for i in 1:nx
        II = ind.b_bottom[1][i]
        opC.χ_b[i+ny, i+ny] = geo.dcap[II,2]
    end
    @inbounds @threads for i in 1:ny
        II = ind.b_right[1][i]
        opC.χ_b[i+ny+nx, i+ny+nx] = geo_p.dcap[δx⁻(II),3]
    end
    @inbounds @threads for i in 1:nx
        II = ind.b_top[1][i]
        opC.χ_b[i+2*ny+nx, i+2*ny+nx] = geo.dcap[II,4]
    end

    return nothing
end

function set_boundary_indicator!(grid::Mesh{GridFCy,T,N}, geo, geo_p, opC) where {T,N}
    @unpack nx, ny, ind = grid
    @inbounds @threads for i in 1:ny
        II = ind.b_left[1][i]
        opC.χ_b[i, i] = geo.dcap[II,1]
    end
    @inbounds @threads for i in 1:nx
        II = ind.b_bottom[1][i]
        opC.χ_b[i+ny, i+ny] = geo_p.dcap[II,2]
    end
    @inbounds @threads for i in 1:ny
        II = ind.b_right[1][i]
        opC.χ_b[i+ny+nx, i+ny+nx] = geo.dcap[II,3]
    end
    @inbounds @threads for i in 1:nx
        II = ind.b_top[1][i]
        opC.χ_b[i+2*ny+nx, i+2*ny+nx] = geo_p.dcap[δy⁻(II),4]
    end

    return nothing
end

function set_border_matrices!(
    grid, geo, grid_u, geo_u, grid_v, geo_v,
    opC_p, opC_u, opC_v,
    periodic_x, periodic_y
    )

    set_boundary_indicator!(grid, geo, geo, opC_p)
    set_boundary_indicator!(grid_u, geo_u, geo, opC_u)
    set_boundary_indicator!(grid_v, geo_v, geo, opC_v)

    mass_matrix_borders!(grid.ind, opC_p.iMx_b, opC_p.iMy_b, opC_p.iMx_bd, opC_p.iMy_bd, geo.dcap, grid.ny)
    mass_matrix_borders!(grid_u.ind, opC_u.iMx_b, opC_u.iMy_b, opC_u.iMx_bd, opC_u.iMy_bd, geo_u.dcap, grid_u.ny)
    mass_matrix_borders!(grid_v.ind, opC_v.iMx_b, opC_v.iMy_b, opC_v.iMx_bd, opC_v.iMy_bd, geo_v.dcap, grid_v.ny)

    bc_matrix_borders!(grid, grid.ind, opC_p.Hx_b, opC_p.Hy_b, geo.dcap)
    mat_assign_T!(opC_p.HxT_b, sparse(opC_p.Hx_b'))
    mat_assign_T!(opC_p.HyT_b, sparse(opC_p.Hy_b'))

    bc_matrix_borders!(grid_u, grid.ind, grid_u.ind, opC_u.Hx_b, opC_u.Hy_b, geo.dcap, geo_u.dcap)
    mat_assign_T!(opC_u.HxT_b, sparse(opC_u.Hx_b'))
    mat_assign_T!(opC_u.HyT_b, sparse(opC_u.Hy_b'))

    bc_matrix_borders!(grid_v, grid.ind, grid_v.ind, opC_v.Hx_b, opC_v.Hy_b, geo.dcap, geo_v.dcap)
    mat_assign_T!(opC_v.HxT_b, sparse(opC_v.Hx_b'))
    mat_assign_T!(opC_v.HyT_b, sparse(opC_v.Hy_b'))

    bc_matrix_borders!(grid, opC_u.Gx_b, opC_v.Gy_b, opC_p.Gx_b, opC_p.Gy_b, geo.dcap)

    periodic_bcs_borders!(grid, opC_p.Hx_b, opC_p.Hy_b, periodic_x, periodic_y)
    periodic_bcs_borders!(grid_u, opC_u.Hx_b, opC_u.Hy_b, periodic_x, periodic_y)
    periodic_bcs_borders!(grid_v, opC_v.Hx_b, opC_v.Hy_b, periodic_x, periodic_y)

    return nothing
end

function laplacian(opC)
    @unpack Bx, By, BxT, ByT, iMx, iMy, tmp_x, tmp_y = opC

    mul!(tmp_x, iMx, Bx)
    L = BxT * tmp_x
    mul!(tmp_y, iMy, By)
    L = L .+ ByT * tmp_y

    return L
end

function laplacian_bc(opC)
    @unpack Bx, By, BxT, ByT, Hx, Hy, iMx, iMy, tmp_x, tmp_y, Hx_b, Hy_b, iMx_b, iMy_b = opC

    mul!(tmp_x, iMx, Hx)
    bc_L = BxT * tmp_x
    mul!(tmp_y, iMy, Hy)
    bc_L = bc_L .+ ByT * tmp_y

    bc_L_b = (BxT * iMx_b * Hx_b .+ ByT * iMy_b * Hy_b)

    return bc_L, bc_L_b
end

function set_matrices!(
    grid, geo, grid_u, geo_u, grid_v, geo_v,
    opC_p, opC_u, opC_v,
    periodic_x, periodic_y
    )
    @unpack ny, ind = grid

    set_other_cutcell_matrices!(
        grid, geo, geo_u, geo_v,
        opC_p, opC_u, opC_v,
        periodic_x, periodic_y
    )

    set_cutcell_matrices!(grid, geo, geo, opC_p, periodic_x, periodic_y)
    set_cutcell_matrices!(grid_u, geo_u, geo, opC_u, periodic_x, periodic_y)
    set_cutcell_matrices!(grid_v, geo_v, geo, opC_v, periodic_x, periodic_y)

    Lp = laplacian(opC_p)
    Lu = laplacian(opC_u)
    Lv = laplacian(opC_v)

    set_border_matrices!(
        grid, geo, grid_u, geo_u, grid_v, geo_v,
        opC_p, opC_u, opC_v,
        periodic_x, periodic_y
    )

    bc_Lp, bc_Lp_b = laplacian_bc(opC_p)
    bc_Lu, bc_Lu_b = laplacian_bc(opC_u)
    bc_Lv, bc_Lv_b = laplacian_bc(opC_v)

    return Lp, bc_Lp, bc_Lp_b, Lu, bc_Lu, bc_Lu_b, Lv, bc_Lv, bc_Lv_b
end

function strain_rate(opC_u, opC_v)
    GxT = opC_u.Gx'
    GyT = opC_v.Gy'

    data = Matrix{SparseMatrixCSC{Float64, Int64}}(undef, 2, 2)
    data[1,1] = 2 .* GxT * (opC_u.HxT * opC_u.iMx * opC_u.Bx .+ opC_u.HyT * opC_u.iMy * opC_u.By)
    data[1,2] = 2 .* GxT * (opC_u.HxT * opC_u.iMx * opC_u.Hx .+ opC_u.HyT * opC_u.iMy * opC_u.Hy)
    data[2,1] = 2 .* GyT * (opC_v.HyT * opC_v.iMy * opC_v.By .+ opC_v.HxT * opC_v.iMx * opC_v.Bx)
    data[2,2] = 2 .* GyT * (opC_v.HyT * opC_v.iMy * opC_v.Hy .+ opC_v.HxT * opC_v.iMx * opC_v.Hx)

    return data
end

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

function set_convection!(
    grid, geo, grid_u, geo_u, grid_v, geo_v,
    u, v, op, ph, BC_u, BC_v
    )
    @unpack Cu, CUTCu, Cv, CUTCv = op
    @unpack uD, vD = ph

    Du_x = zeros(grid_u)
    Du_y = zeros(grid_u)
    Du_x .= reshape(veci(uD,grid_u,1), grid_u)
    Du_y .= reshape(veci(uD,grid_u,1), grid_u)
    Du_x[grid_u.ind.MIXED] .= reshape(veci(uD,grid_u,2), grid_u)[grid_u.ind.MIXED]
    Du_y[grid_u.ind.MIXED] .= reshape(veci(uD,grid_u,2), grid_u)[grid_u.ind.MIXED]
    Du_x[:,1] .= vec3_L(uD,grid_u)
    Du_x[:,2] .= u[:,2]
    Du_y[1,:] .= vec3_B(uD,grid_u)
    Du_y[2,:] .= u[2,:]
    Du_x[:,end] .= vec3_R(uD,grid_u)
    Du_x[:,end-1] .= u[:,end-1]
    Du_y[end,:] .= vec3_T(uD,grid_u)
    Du_y[end-1,:] .= u[end-1,:]

    Dv_x = zeros(grid_v)
    Dv_y = zeros(grid_v)
    Dv_x .= reshape(veci(vD,grid_v,1), grid_v)
    Dv_y .= reshape(veci(vD,grid_v,1), grid_v)
    Dv_x[grid_v.ind.MIXED] .= reshape(veci(vD,grid_v,2), grid_v)[grid_v.ind.MIXED]
    Dv_y[grid_v.ind.MIXED] .= reshape(veci(vD,grid_v,2), grid_v)[grid_v.ind.MIXED]
    Dv_x[:,1] .= vec3_L(vD,grid_v)
    Dv_x[:,2] .= v[:,2]
    Dv_y[1,:] .= vec3_B(vD,grid_v)
    Dv_y[2,:] .= v[2,:]
    Dv_x[:,end] .= vec3_R(vD,grid_v)
    Dv_x[:,end-1] .= v[:,end-1]
    Dv_y[end,:] .= vec3_T(vD,grid_v)
    Dv_y[end-1,:] .= v[end-1,:]

    bnds_u = [grid_u.ind.b_left[1], grid_u.ind.b_bottom[1], grid_u.ind.b_right[1], grid_u.ind.b_top[1]]
    bnds_v = [grid_v.ind.b_left[1], grid_v.ind.b_bottom[1], grid_v.ind.b_right[1], grid_v.ind.b_top[1]]
    Δu = [grid_u.dx[1,1], grid_u.dy[1,1], grid_u.dx[end,end], grid_u.dy[end,end]] .* 0.5
    Δv = [grid_v.dx[1,1], grid_v.dy[1,1], grid_v.dx[end,end], grid_v.dy[end,end]] .* 0.5

    Hu = zeros(grid_u)
    for i in eachindex(bnds_u)
        for II in bnds_u[i]
            Hu[II] = Δu[i]
        end
    end

    Hv = zeros(grid_v)
    for i in eachindex(bnds_v)
        for II in bnds_v[i]
            Hv[II] = Δv[i]
        end
    end

    set_bc_bnds(dir, GridFCx, Du_x, Du_y, Dv_x, Dv_y, Hu, Hv, u, v, BC_u, BC_v)
    set_bc_bnds(dir, GridFCy, Dv_x, Dv_y, Du_x, Du_y, Hv, Hu, v, u, BC_v, BC_u)

    # Du_x .= 0.0
    # Du_y .= 0.0
    # Dv_x .= 0.0
    # Dv_y .= 0.0

    vector_convection!(dir, GridFCx, Cu, CUTCu, u, v, Du_x, Du_y, Dv_x, Dv_y,
            geo.dcap, grid.nx, grid.ny, BC_u, grid_u.ind.inside,
            grid_u.ind.b_left[1], grid_u.ind.b_bottom[1], grid_u.ind.b_right[1], grid_u.ind.b_top[1])
    vector_convection!(dir, GridFCy, Cv, CUTCv, u, v, Du_x, Du_y, Dv_x, Dv_y,
            geo.dcap, grid.nx, grid.ny, BC_v, grid_v.ind.inside,
            grid_v.ind.b_left[1], grid_v.ind.b_bottom[1], grid_v.ind.b_right[1], grid_v.ind.b_top[1])
    
    return nothing
end

function CN_set_momentum(
    bc_type, num, grid, opC,
    A, B,
    L, bc_L, bc_L_b,
    Lm1, bc_Lm1, bc_Lm1_b, Mm1, BC,
    ls_advection
    )
    @unpack τ = num
    @unpack Bx, By, Hx, Hy, HxT, HyT, χ, M, iMx, iMy, Hx_b, Hy_b, HxT_b, HyT_b, iMx_b, iMy_b, iMx_bd, iMy_bd, χ_b = opC

    if is_dirichlet(bc_type)
        vel = copy(grid.V)
        __a1 = -1.0
        __b = 0.0
    elseif is_neumann(bc_type)
        vel = 0.0
        __a1 = 0.0
        __b = 1.0
    elseif is_robin(bc_type)
        vel = 0.0
        __a1 = -1.0
        __b = 1.0
    end

    # Number of points in outer boundaries (corners duplicated)
    nb = 2 * grid.nx + 2 * grid.ny

    a0 = ones(grid) .* vel
    _a1 = ones(grid) .* __a1
    a1 = Diagonal(vec(_a1))
    _b = ones(grid) .* __b
    b = Diagonal(vec(_b))

    a0_b = zeros(nb)
    _a1_b = zeros(nb)
    _b_b = zeros(nb)
    set_borders!(grid, a0_b, _a1_b, _b_b, BC, num.n_ext_cl)
    a1_b = Diagonal(vec(_a1_b))
    b_b = Diagonal(vec(_b_b))

    τ_2 = 0.5 * τ

    if ls_advection
        data_A = Matrix{SparseMatrixCSC{Float64, Int64}}(undef, 3, 3)
        # Implicit part of viscous term
        data_A[1,1] = pad_crank_nicolson(M .- τ_2 .* L, grid, τ)
        # Contribution to implicit part of viscous term from inner boundaries
        data_A[1,2] = - τ_2 .* bc_L
        # Contribution to implicit part of viscous term from outer boundaries
        data_A[1,3] = - τ_2 .* bc_L_b
        # Boundary conditions for inner boundaries
        data_A[2,1] = b * (HxT * iMx * Bx .+ HyT * iMy * By)
        data_A[2,2] = pad(b * (HxT * iMx * Hx .+ HyT * iMy * Hy) .- χ * a1)
        data_A[2,3] = spdiagm(grid.nx*grid.ny, nb, zeros(nb))
        # Boundary conditions for outer boundaries
        data_A[3,1] = b_b * (HxT_b * iMx_b' * Bx .+ HyT_b * iMy_b' * By)
        data_A[3,2] = spdiagm(nb, grid.nx*grid.ny, zeros(nb))
        data_A[3,3] = pad(b_b * (HxT_b * iMx_bd * Hx_b .+ HyT_b * iMy_bd * Hy_b) .- χ_b * a1_b)
        A  = [data_A[1,1] data_A[1,2] data_A[1,3];
            data_A[2,1] data_A[2,2] data_A[2,3];
            data_A[3,1] data_A[3,2] data_A[3,3]]

        data_B = Matrix{SparseMatrixCSC{Float64, Int64}}(undef, 3, 3)
        # Explicit part of viscous term
        data_B[1,1] = Mm1 .+ τ_2 .* Lm1
        # Contribution to implicit part of viscous term from inner boundaries
        data_B[1,2] = τ_2 .* bc_Lm1
        # Contribution to implicit part of viscous term from outer boundaries
        data_B[1,3] = τ_2 .* bc_Lm1_b
        data_B[2,1] = spdiagm(0 => fzeros(grid))
        data_B[2,2] = spdiagm(0 => fzeros(grid))
        data_B[2,3] = spdiagm(grid.nx*grid.ny, nb, zeros(nb))
        data_B[3,1] = spdiagm(nb, grid.nx*grid.ny, zeros(nb))
        data_B[3,2] = spdiagm(nb, grid.nx*grid.ny, zeros(nb))
        data_B[3,3] = spdiagm(nb, nb, zeros(nb))
        B  = [data_B[1,1] data_B[1,2] data_B[1,3];
            data_B[2,1] data_B[2,2] data_B[2,3];
            data_B[3,1] data_B[3,2] data_B[3,3]]
    end

    rhs = f3zeros(grid)
    veci(rhs,grid,2) .= χ * vec(a0)
    vec3(rhs,grid) .= χ_b * vec(a0_b)

    return A, B, rhs
end

function FE_set_momentum(
    bc_type, num, grid, opC,
    A, B,
    L, bc_L, bc_L_b, Mm1, BC,
    ls_advection
    )
    @unpack τ = num
    @unpack Bx, By, Hx, Hy, HxT, HyT, χ, M, iMx, iMy, Hx_b, Hy_b, HxT_b, HyT_b, iMx_b, iMy_b, iMx_bd, iMy_bd, χ_b = opC

    if is_dirichlet(bc_type)
        vel = copy(grid.V)
        __a1 = -1.0
        __b = 0.0
    elseif is_neumann(bc_type)
        vel = 0.0
        __a1 = 0.0
        __b = 1.0
    elseif is_robin(bc_type)
        vel = 0.0
        __a1 = -1.0
        __b = 1.0
    end

    # Number of points in outer boundaries (corners duplicated)
    nb = 2 * grid.nx + 2 * grid.ny

    a0 = ones(grid) .* vel
    _a1 = ones(grid) .* __a1
    a1 = Diagonal(vec(_a1))
    _b = ones(grid) .* __b
    b = Diagonal(vec(_b))

    a0_b = zeros(nb)
    _a1_b = zeros(nb)
    _b_b = zeros(nb)
    set_borders!(grid, a0_b, _a1_b, _b_b, BC, num.n_ext_cl)
    a1_b = Diagonal(vec(_a1_b))
    b_b = Diagonal(vec(_b_b))

    if ls_advection
        data_A = Matrix{SparseMatrixCSC{Float64, Int64}}(undef, 3, 3)
        # Implicit part of viscous term
        data_A[1,1] = pad_crank_nicolson(M .- τ .* L, grid, τ)
        # Contribution to implicit part of viscous term from inner boundaries
        data_A[1,2] = - τ .* bc_L
        # Contribution to implicit part of viscous term from outer boundaries
        data_A[1,3] = - τ .* bc_L_b
        # Boundary conditions for inner boundaries
        data_A[2,1] = b * (HxT * iMx * Bx .+ HyT * iMy * By)
        data_A[2,2] = pad(b * (HxT * iMx * Hx .+ HyT * iMy * Hy) .- χ * a1)
        data_A[2,3] = spdiagm(grid.nx*grid.ny, nb, zeros(nb))
        # Boundary conditions for outer boundaries
        data_A[3,1] = b_b * (HxT_b * iMx_b' * Bx .+ HyT_b * iMy_b' * By)
        data_A[3,2] = spdiagm(nb, grid.nx*grid.ny, zeros(nb))
        data_A[3,3] = pad(b_b * (HxT_b * iMx_bd * Hx_b .+ HyT_b * iMy_bd * Hy_b) .- χ_b * a1_b)
        A  = [data_A[1,1] data_A[1,2] data_A[1,3];
            data_A[2,1] data_A[2,2] data_A[2,3];
            data_A[3,1] data_A[3,2] data_A[3,3]]

        data_B = Matrix{SparseMatrixCSC{Float64, Int64}}(undef, 3, 3)
        data_B[1,1] = Mm1
        data_B[1,2] = spdiagm(0 => fzeros(grid))
        data_B[1,3] = spdiagm(grid.nx*grid.ny, nb, zeros(nb))
        data_B[2,1] = spdiagm(0 => fzeros(grid))
        data_B[2,2] = spdiagm(0 => fzeros(grid))
        data_B[2,3] = spdiagm(grid.nx*grid.ny, nb, zeros(nb))
        data_B[3,1] = spdiagm(nb, grid.nx*grid.ny, zeros(nb))
        data_B[3,2] = spdiagm(nb, grid.nx*grid.ny, zeros(nb))
        data_B[3,3] = spdiagm(nb, nb, zeros(nb))
        B  = [data_B[1,1] data_B[1,2] data_B[1,3];
            data_B[2,1] data_B[2,2] data_B[2,3];
            data_B[3,1] data_B[3,2] data_B[3,3]]
    end

    rhs = f3zeros(grid)
    veci(rhs,grid,2) .= χ * vec(a0)
    vec3(rhs,grid) .= χ_b * vec(a0_b)
    
    return A, B, rhs
end

function set_poisson(
    bc_type, num, grid, a0, opC, opC_u, opC_v,
    A, L, bc_L, bc_L_b, BC,
    ls_advection)
    @unpack Bx, By, Hx, Hy, HxT, HyT, χ, M, iMx, iMy, Hx_b, Hy_b, HxT_b, HyT_b, iMx_b, iMy_b, iMx_bd, iMy_bd, χ_b = opC

    if is_dirichlet(bc_type)
        __a1 = -1.0
        __a2 = 0.0
        __b = 0.0
    elseif is_neumann(bc_type)
        __a1 = 0.0
        __a2 = 0.0
        __b = 1.0
    elseif is_robin(bc_type)
        __a1 = -1.0
        __a2 = 0.0
        __b = 1.0
    elseif is_fs(bc_type)
        __a1 = 0.0
        __a2 = 1.0
        __b = 0.0
    end

    # Number of points in outer boundaries (corners duplicated)
    nb = 2 * grid.nx + 2 * grid.ny

    _a1 = ones(grid) .* __a1
    a1 = Diagonal(vec(_a1))
    _a2 = ones(grid) .* __a2
    a2 = Diagonal(vec(_a2))
    _b = ones(grid) .* __b
    b = Diagonal(vec(_b))

    a0_b = zeros(nb)
    _a1_b = zeros(nb)
    _b_b = zeros(nb)
    set_borders!(grid, a0_b, _a1_b, _b_b, BC, num.n_ext_cl)
    a1_b = Diagonal(vec(_a1_b))
    b_b = Diagonal(vec(_b_b))

    if ls_advection
        GxT = opC_u.Gx'
        GyT = opC_v.Gy'
        
        data_A = Matrix{SparseMatrixCSC{Float64, Int64}}(undef, 3, 3)
        # Poisson equation
        data_A[1,1] = pad(L, -4.0)
        data_A[1,2] = bc_L
        data_A[1,3] = bc_L_b
        # Boundary conditions for inner boundaries
        data_A[2,1] = -b * (HxT * iMx * Bx .+ HyT * iMy * By)
        data_A[2,2] = -pad(
            b * (HxT * iMx * Hx .+ HyT * iMy * Hy) .- χ * a1 .+
            a2 * (GxT * opC_u.Gx .+ GyT * opC_v.Gy), 4.0
        )
        data_A[2,3] = spdiagm(grid.nx*grid.ny, nb, zeros(nb))
        # Boundary conditions for outer boundaries
        data_A[3,1] = -b_b * (HxT_b * iMx_b' * Bx .+ HyT_b * iMy_b' * By)
        data_A[3,2] = -b_b * (HxT_b * iMx_b' * Hx .+ HyT_b * iMy_b' * Hy)
        data_A[3,3] = -pad(b_b * (HxT_b * iMx_bd * Hx_b .+ HyT_b * iMy_bd * Hy_b) .- χ_b * a1_b, 4.0)
        A  = [data_A[1,1] data_A[1,2] data_A[1,3];
            data_A[2,1] data_A[2,2] data_A[2,3];
            data_A[3,1] data_A[3,2] data_A[3,3]]
    end

    rhs = f3zeros(grid)
    veci(rhs,grid,2) .= -χ * vec(a0)
    vec3(rhs,grid) .= -χ_b * vec(a0_b)
    
    return A, rhs
end

function set_CN!(
    bc_type_u, bc_type_p, num, grid, geo, grid_u, geo_u, grid_v, geo_v,
    opC_p, opC_u, opC_v, BC_p, BC_u, BC_v,
    Au, Bu, Av, Bv, Aϕ,
    Lpm1, bc_Lpm1, bc_Lpm1_b, Lum1, bc_Lum1, bc_Lum1_b, Lvm1, bc_Lvm1, bc_Lvm1_b, 
    Mum1, Mvm1, iRe, op_conv, ph,
    periodic_x, periodic_y, advection, ls_advection
    )

    if advection
        set_convection!(grid, geo, grid_u, geo_u, grid_v, geo_v, ph.u, ph.v, op_conv, ph, BC_u, BC_v)
    end

    if ls_advection
        update_ls_data(num, grid, grid_u, grid_v, grid.u, grid.κ, periodic_x, periodic_y, false)

        laps = set_matrices!(
            grid, geo, grid_u, geo_u, grid_v, geo_v,
            opC_p, opC_u, opC_v,
            periodic_x, periodic_y
        )
    else
        laps = Lpm1, bc_Lpm1, bc_Lpm1_b, Lum1, bc_Lum1, bc_Lum1_b, Lvm1, bc_Lvm1, bc_Lvm1_b
    end
    Lp, bc_Lp, bc_Lp_b, Lu, bc_Lu, bc_Lu_b, Lv, bc_Lv, bc_Lv_b = laps

    Au, Bu, rhs_u = CN_set_momentum(
        bc_type_u, num, grid_u, opC_u,
        Au, Bu,
        iRe.*Lu, iRe.*bc_Lu, iRe.*bc_Lu_b,
        iRe.*Lum1, iRe.*bc_Lum1, iRe.*bc_Lum1_b, Mum1, BC_u,
        ls_advection
    )
    Av, Bv, rhs_v = CN_set_momentum(
        bc_type_u, num, grid_v, opC_v,
        Av, Bv,
        iRe.*Lv, iRe.*bc_Lv, iRe.*bc_Lv_b,
        iRe.*Lvm1, iRe.*bc_Lvm1, iRe.*bc_Lvm1_b, Mvm1, BC_v,
        ls_advection
    )
    a0_p = zeros(grid)
    Aϕ, rhs_ϕ = set_poisson(
        bc_type_p, num, grid, a0_p, opC_p, opC_u, opC_v,
        Aϕ, Lp, bc_Lp, bc_Lp_b, BC_p,
        ls_advection
    )

    return Au, Bu, rhs_u, Av, Bv, rhs_v, Aϕ, rhs_ϕ, Lp, bc_Lp, bc_Lp_b, Lu, bc_Lu, bc_Lu_b, Lv, bc_Lv, bc_Lv_b
end

function set_FE!(
    bc_type_u, bc_type_p, num, grid, geo, grid_u, geo_u, grid_v, geo_v,
    opC_p, opC_u, opC_v, BC_p, BC_u, BC_v,
    Au, Bu, Av, Bv, Aϕ,
    Lpm1, bc_Lpm1, bc_Lpm1_b, Lum1, bc_Lum1, bc_Lum1_b, Lvm1, bc_Lvm1, bc_Lvm1_b,
    Mum1, Mvm1, iRe, op_conv, ph,
    periodic_x, periodic_y, advection, ls_advection
    )

    if advection
        set_convection!(grid, geo, grid_u, geo_u, grid_v, geo_v, ph.u, ph.v, op_conv, ph, BC_u, BC_v)
    end

    if ls_advection
        update_ls_data(num, grid, grid_u, grid_v, grid.u, grid.κ, periodic_x, periodic_y, false)

        laps = set_matrices!(
            grid, geo, grid_u, geo_u, grid_v, geo_v,
            opC_p, opC_u, opC_v,
            periodic_x, periodic_y
        )
    else
        laps = Lpm1, bc_Lpm1, bc_Lpm1_b, Lum1, bc_Lum1, bc_Lum1_b, Lvm1, bc_Lvm1, bc_Lvm1_b
    end
    Lp, bc_Lp, bc_Lp_b, Lu, bc_Lu, bc_Lu_b, Lv, bc_Lv, bc_Lv_b = laps

    Au, Bu, rhs_u = FE_set_momentum(
        bc_type_u, num, grid_u, opC_u,
        Au, Bu,
        iRe.*Lu, iRe.*bc_Lu, iRe.*bc_Lu_b, Mum1, BC_u,
        ls_advection
    )
    Av, Bv, rhs_v = FE_set_momentum(
        bc_type_u, num, grid_v, opC_v,
        Av, Bv,
        iRe.*Lv, iRe.*bc_Lv, iRe.*bc_Lv_b, Mvm1, BC_v,
        ls_advection
    )
    a0_p = zeros(grid)
    Aϕ, rhs_ϕ = set_poisson(
        bc_type_p, num, grid, a0_p, opC_p, opC_u, opC_v,
        Aϕ, Lp, bc_Lp, bc_Lp_b, BC_p,
        ls_advection
    )

    return Au, Bu, rhs_u, Av, Bv, rhs_v, Aϕ, rhs_ϕ, Lp, bc_Lp, bc_Lp_b, Lu, bc_Lu, bc_Lu_b, Lv, bc_Lv, bc_Lv_b
end

function pressure_projection!(
    time_scheme, bc_type_u, bc_type_p,
    num, grid, geo, grid_u, geo_u, grid_v, geo_v, ph,
    BC_u, BC_v, BC_p,
    opC_p, opC_u, opC_v, op_conv,
    Au, Bu, Av, Bv, Aϕ,
    Lpm1, bc_Lpm1, bc_Lpm1_b, Lum1, bc_Lum1, bc_Lum1_b, Lvm1, bc_Lvm1, bc_Lvm1_b,
    Cum1, Cvm1, Mum1, Mvm1,
    periodic_x, periodic_y, advection, ls_advection, current_i
    )
    @unpack Re, τ, σ, g, β = num
    @unpack p, pD, ϕ, ϕD, u, v, ucorrD, vcorrD, uD, vD, ucorr, vcorr = ph
    @unpack Cu, Cv, CUTCu, CUTCv = op_conv

    iRe = 1.0 / Re
    iτ = 1.0 / τ

    if is_FE(time_scheme)
        Au, Bu, rhs_u, Av, Bv, rhs_v, Aϕ, rhs_ϕ, Lp, bc_Lp, bc_Lp_b, Lu, bc_Lu, bc_Lu_b, Lv, bc_Lv, bc_Lv_b = set_FE!(
            bc_type_u, bc_type_p, num, grid, geo, grid_u, geo_u, grid_v, geo_v,
            opC_p, opC_u, opC_v, BC_p, BC_u, BC_v,
            Au, Bu, Av, Bv, Aϕ,
            Lpm1, bc_Lpm1, bc_Lpm1_b, Lum1, bc_Lum1, bc_Lum1_b, Lvm1, bc_Lvm1, bc_Lvm1_b,
            Mum1, Mvm1, iRe, op_conv, ph,
            periodic_x, periodic_y, advection, ls_advection
        )
    elseif is_CN(time_scheme)
        Au, Bu, rhs_u, Av, Bv, rhs_v, Aϕ, rhs_ϕ, Lp, bc_Lp, bc_Lp_b, Lu, bc_Lu, bc_Lu_b, Lv, bc_Lv, bc_Lv_b = set_CN!(
            bc_type_u, bc_type_p, num, grid, geo, grid_u, geo_u, grid_v, geo_v,
            opC_p, opC_u, opC_v, BC_p, BC_u, BC_v,
            Au, Bu, Av, Bv, Aϕ,
            Lpm1, bc_Lpm1, bc_Lpm1_b, Lum1, bc_Lum1, bc_Lum1_b, Lvm1, bc_Lvm1, bc_Lvm1_b,
            Mum1, Mvm1, iRe, op_conv, ph,
            periodic_x, periodic_y, advection, ls_advection
        )
    end

    grav_x = g .* sin(β) .* opC_u.M * fones(grid_u)
    grav_y = g .* cos(β) .* opC_v.M * fones(grid_v)

    Convu = fzeros(grid_u)
    Convv = fzeros(grid_v)
    Cui = Cu * vec(u) .+ CUTCu
    Cvi = Cv * vec(v) .+ CUTCv
    if advection
        if current_i == 1
            Convu .+= Cui
            Convv .+= Cvi
        else
            Convu .+= 1.5 .* Cui .- 0.5 .* Cum1
            Convv .+= 1.5 .* Cvi .- 0.5 .* Cvm1
        end
    end

    if is_dirichlet(bc_type_u)
        veci(uD,grid_u,1) .= vec(u)
        # update_dirichlet_field!(grid_u, uD, u, BC_u)
        veci(rhs_u,grid_u,1) .+= -τ .* (opC_u.AxT * opC_u.Rx * veci(pD,grid,1) .+ opC_u.Gx * veci(pD,grid,2) .+ opC_u.Gx_b * vec3(pD,grid))
    end
    mul!(rhs_u, Bu, uD, 1.0, 1.0)
    veci(rhs_u,grid_u,1) .+= τ .* grav_x
    veci(rhs_u,grid_u,1) .-= τ .* Convu
    kill_dead_cells!(veci(rhs_u,grid_u,1), grid_u, geo_u)
    kill_dead_cells!(veci(rhs_u,grid_u,2), grid_u, geo_u)
    # blocks = DDM.decompose(Au, grid_u.domdec, grid_u.domdec)
    # bicgstabl!(ucorrD, Au, rhs_u, Pl=ras(blocks,grid_u.pou), log=true)
    # @time bicgstabl!(ucorrD, Au, rhs_u, log=true)
    try
        # @time bicgstabl!(ucorrD, Au, rhs_u, Pl=Diagonal(Au), log=true)
        @time ucorrD .= Au \ rhs_u
    catch e
        ucorrD .= Inf
        println(e)
    end
    # @mytime _, ch = bicgstabl!(ucorrD, Au, rhs_u, Pl=ras(blocks,grid_u.pou), log=true)
    # println(ch)
    ucorr .= reshape(veci(ucorrD,grid_u,1), grid_u)

    if is_dirichlet(bc_type_u)
        veci(vD,grid_v,1) .= vec(v)
        # update_dirichlet_field!(grid_v, vD, v, BC_v)
        veci(rhs_v,grid_v,1) .+= -τ .* (opC_v.AyT * opC_v.Ry * veci(pD,grid,1) .+ opC_v.Gy * veci(pD,grid,2) .+ opC_v.Gy_b * vec3(pD,grid))
    end
    mul!(rhs_v, Bv, vD, 1.0, 1.0)
    veci(rhs_v,grid_v,1) .+= - τ .* grav_y
    veci(rhs_v,grid_v,1) .-= τ .* Convv
    kill_dead_cells!(veci(rhs_v,grid_v,1), grid_v, geo_v)
    kill_dead_cells!(veci(rhs_v,grid_v,2), grid_v, geo_v)
    # blocks = DDM.decompose(Av, grid_v.domdec, grid_v.domdec)
    # bicgstabl!(vcorrD, Av, rhs_v, Pl=ras(blocks,grid_v.pou), log=true)
    # bicgstabl!(vcorrD, Av, rhs_v, log=true)
    try
        # @time bicgstabl!(vcorrD, Av, rhs_v, Pl=Diagonal(Av), log=true)
        @time vcorrD .= Av \ rhs_v
    catch e
        vcorrD .= Inf
        println(e)
    end
    # @mytime _, ch = bicgstabl!(vcorrD, Av, rhs_v, Pl=ras(blocks,grid_v.pou), log=true)
    # println(ch)
    vcorr .= reshape(veci(vcorrD,grid_v,1), grid_v)

    Duv = opC_p.AxT * veci(ucorrD,grid_u,1) .+ opC_p.Gx * veci(ucorrD,grid_u,2) .+ opC_p.Gx_b * vec3(ucorrD,grid_u) .+
          opC_p.AyT * veci(vcorrD,grid_v,1) .+ opC_p.Gy * veci(vcorrD,grid_v,2) .+ opC_p.Gy_b * vec3(vcorrD,grid_v)
    veci(rhs_ϕ,grid,1) .= iτ .* Duv

    if is_fs(bc_type_p)
        Smat = strain_rate(opC_u, opC_v)
        S = Smat[1,1] * veci(ucorrD,grid_u,1) .+ Smat[1,2] * veci(ucorrD,grid_u,2) .+
            Smat[2,1] * veci(vcorrD,grid_v,1) .+ Smat[2,2] * veci(vcorrD,grid_v,2)

        GxT = opC_u.Gx'
        GyT = opC_v.Gy'
        veci(rhs_ϕ,grid,2) .= -iRe .* S .+ σ .* (GxT * opC_u.Gx .+ GyT * opC_v.Gy) * vec(grid.κ)
    end
    # Remove nullspace by adding small quantity to main diagonal
    @inbounds @threads for i in 1:Aϕ.m
        @inbounds Aϕ[i,i] += 1e-10
    end
    # blocks = DDM.decompose(Aϕ, grid.domdec, grid.domdec)
    # bicgstabl!(ϕD, Aϕ, rhs_ϕ, Pl = ras(blocks,grid.pou), log = true)
    # @time bicgstabl!(ϕD, Aϕ, rhs_ϕ, Pl = Diagonal(Aϕ), log = true)
    @time ϕD .= Aϕ \ rhs_ϕ
    # @mytime _, ch = bicgstabl!(ϕD, Aϕ, rhs_ϕ, Pl = ras(blocks,grid.pou), log = true)
    # println(ch)
    ϕ .= reshape(veci(ϕD,grid,1), grid)

    iMu = Diagonal(1 ./ (opC_u.M.diag .+ eps(0.01)))
    iMv = Diagonal(1 ./ (opC_v.M.diag .+ eps(0.01)))

    ∇ϕ_x = opC_u.AxT * opC_u.Rx * vec(ϕ) .+ opC_u.Gx * veci(ϕD,grid,2) .+ opC_u.Gx_b * vec3(ϕD,grid)
    ∇ϕ_y = opC_v.AyT * opC_v.Ry * vec(ϕ) .+ opC_v.Gy * veci(ϕD,grid,2) .+ opC_v.Gy_b * vec3(ϕD,grid)

    iM = Diagonal(1. ./ (vec(geo.dcap[:,:,5]) .+ eps(0.01)))
    if is_fs(bc_type_p)
        veci(pD,grid,1) .= vec(ϕ) #.- iRe .* reshape(iM * Duv, grid))
        veci(pD,grid,2) .= veci(ϕD,grid,2)
        vec3(pD,grid) .= vec3(ϕD,grid)
        p .= reshape(veci(pD,grid,1), grid)
    else
        veci(pD,grid,1) .= vec(p) .+ vec(ϕ) #.- iRe .* iM * Duv
        veci(pD,grid,2) .+= veci(ϕD,grid,2)
        vec3(pD,grid) .+= vec3(ϕD,grid)
        p .= reshape(veci(pD,grid,1), grid)
    end
    
    u .= ucorr .- τ .* reshape(iMu * ∇ϕ_x, grid_u)
    v .= vcorr .- τ .* reshape(iMv * ∇ϕ_y, grid_v)

    veci(uD,grid_u,1) .= vec(u)
    veci(uD,grid_u,2) .= veci(ucorrD,grid_u,2)
    vec3(uD,grid_u) .= vec3(ucorrD,grid_u)
    veci(vD,grid_v,1) .= vec(v)
    veci(vD,grid_v,2) .= veci(vcorrD,grid_v,2)
    vec3(vD,grid_v) .= vec3(vcorrD,grid_v)

    return Aϕ, Lp, bc_Lp, bc_Lp_b, Au, Bu, Lu, bc_Lu, bc_Lu_b, Av, Bv, Lv, bc_Lv, bc_Lv_b, opC_p.M, opC_u.M, opC_v.M, Cui, Cvi
end

"""
    linear_advection!(
        num, grid, geo, grid_u, geo_u, grid_v, geo_v, ph,
        BC_u, BC_v, op_conv
    )

Solves the linear advection equation for a vector field (u, v) with slip BCs.

Convective term is solved implicitly using the midopint rule coupled with a Newton
algorithm.
"""
function linear_advection!(
    num, grid, geo, grid_u, geo_u, grid_v, geo_v, ph,
    BC_u, BC_v, op_conv
    )
    @unpack τ = num
    @unpack u, v, uD, vD = ph
    @unpack Cu, Cv, CUTCu, CUTCv = op_conv

    Convu = fzeros(grid_u)
    rhs_u = f3zeros(grid_u)
    Convv = fzeros(grid_v)
    rhs_v = f3zeros(grid_v)

    res = vcat(fones(grid_u), fones(grid_v))
    dres = vcat(fzeros(grid_u), fzeros(grid_v))

    vec1(uD, grid_u) .= vec(u)
    kill_dead_cells!(vec1(uD, grid_u), grid_u, geo_u)
    kill_dead_cells!(vec2(uD, grid_u), grid_u, geo_u)
    u .= reshape(vec1(uD, grid_u), grid_u)
    u_guess = copy(u)

    vec1(vD, grid_v) .= vec(v)
    kill_dead_cells!(vec1(vD, grid_v), grid_v, geo_v)
    kill_dead_cells!(vec2(vD, grid_v), grid_v, geo_v)
    v .= reshape(vec1(vD, grid_v), grid_v)
    v_guess = copy(v)

    # Newton algorithm to solve the implicit midpoint method
    i = 0
    while sum(res) > 1e-8
        ϵ = 1e-10
        
        res .= residual(u_guess, v_guess, num, grid, geo, grid_u, geo_u, grid_v, geo_v, u, v, op_conv, ph, BC_u, BC_v)
        dres .= dresidual(u_guess, v_guess, res, ϵ, num, grid, geo, grid_u, geo_u, grid_v, geo_v, u, v, op_conv, ph, BC_u, BC_v)

        u_guess .= u_guess .- reshape((res ./ dres)[1:grid_u.nx*grid_u.ny], grid_u)
        v_guess .= v_guess .- reshape((res ./ dres)[grid_u.nx*grid_u.ny+1:end], grid_v)

        i += 1
    end

    u_midp = 0.5 .* (u .+ u_guess)
    v_midp = 0.5 .* (v .+ v_guess)
    set_convection!(grid, geo, grid_u, geo_u, grid_v, geo_v, u_midp, v_midp, op_conv, ph, BC_u, BC_v)

    Convu .= Cu * vec(u_midp) .+ CUTCu
    vec1(rhs_u, grid_u) .-= τ .* Convu
    kill_dead_cells!(vec1(rhs_u, grid_u), grid_u, geo_u)
    kill_dead_cells!(vec2(rhs_u, grid_u), grid_u, geo_u)

    Convv .= Cv * vec(v_midp) .+ CUTCv
    vec1(rhs_v, grid_v) .-= τ .* Convv
    kill_dead_cells!(vec1(rhs_v, grid_v), grid_v, geo_v)
    kill_dead_cells!(vec2(rhs_v, grid_v), grid_v, geo_v)
    
    uD .= uD .+ rhs_u
    u .= reshape(vec1(uD, grid_u), grid_u)
    vD .= vD .+ rhs_v
    v .= reshape(vec1(vD, grid_v), grid_v)

    return nothing
end

function residual(u_guess, v_guess, num, grid, geo, grid_u, geo_u, grid_v, geo_v, u, v, op_conv, ph, BC_u, BC_v)
    @unpack τ = num
    @unpack uD, vD = ph
    @unpack Cu, Cv, CUTCu, CUTCv = op_conv

    res = vcat(fones(grid_u), fones(grid_v))

    Convu = fzeros(grid_u)
    rhs_u = f3zeros(grid_u)
    _u = copy(u)
    _uD = copy(uD)

    Convv = fzeros(grid_v)
    rhs_v = f3zeros(grid_v)
    _v = copy(v)
    _vD = copy(vD)

    u_midp = 0.5 .* (u .+ u_guess)
    v_midp = 0.5 .* (v .+ v_guess)

    set_convection!(grid, geo, grid_u, geo_u, grid_v, geo_v, u_midp, v_midp, op_conv, ph, BC_u, BC_v)
    Convu .= Cu * vec(u_midp) .+ CUTCu
    Convv .= Cv * vec(v_midp) .+ CUTCv
    vec1(rhs_u, grid_u) .-= τ .* Convu
    vec1(rhs_v, grid_v) .-= τ .* Convv
    
    kill_dead_cells!(vec1(rhs_u, grid_u), grid_u, geo_u)
    kill_dead_cells!(vec2(rhs_u, grid_u), grid_u, geo_u)
    kill_dead_cells!(vec1(rhs_v, grid_v), grid_v, geo_v)
    kill_dead_cells!(vec2(rhs_v, grid_v), grid_v, geo_v)

    _uD .= uD .+ rhs_u
    _u .= reshape(vec1(_uD, grid_u), grid_u)
    _vD .= vD .+ rhs_v
    _v .= reshape(vec1(_vD, grid_v), grid_v)

    res .= vcat(vec(abs.(_u .- u_guess)), vec(abs.(_v .- v_guess)))

    return res
end

function dresidual(u_guess, v_guess, res0, eps, num, grid, geo, grid_u, geo_u, grid_v, geo_v, u, v, op_conv, ph, BC_u, BC_v)
    res1 = residual(u_guess .+ eps, v_guess .+ eps, num, grid, geo, grid_u, geo_u, grid_v, geo_v, u, v, op_conv, ph, BC_u, BC_v)
    dres = (res1 .- res0) ./ eps

    return dres
end