function set_borders!(grid, a0, a1, b, BC)
    @unpack ny, ind = grid
    
    @inbounds a0[2:end-1,1] .= BC.left.val
    if is_dirichlet(BC.left.t)
        @inbounds a1[2:end-1,1] .= -1.
        @inbounds b[2:end-1,1] .= 0.
    elseif is_neumann(BC.left.t)
        @inbounds a1[2:end-1,1] .= 0.
        @inbounds b[2:end-1,1] .= 1.
    elseif is_robin(BC.left.t)
        @inbounds a1[2:end-1,1] .= -1.
        @inbounds b[2:end-1,1] .= 1.
    elseif is_periodic(BC.left.t)
        # @inbounds a1[:,1] .= 0.
        # @inbounds b[:,1] .= 0.
    else
        @error ("Not implemented yet")
    end
    @inbounds a0[1,:] .= BC.bottom.val
    if is_dirichlet(BC.bottom.t)
        @inbounds a1[1,:] .= -1.
        @inbounds b[1,:] .= 0.
    elseif is_neumann(BC.bottom.t)
        @inbounds a1[1,:] .= 0.
        @inbounds b[1,:] .= 1.
    elseif is_robin(BC.bottom.t)
        @inbounds a1[1,:] .= -1.
        @inbounds b[1,:] .= 1.
    elseif is_periodic(BC.bottom.t)
        # @inbounds a1[1,2:end-1] .= 0.
        # @inbounds b[1,2:end-1] .= 0.
    else
        @error ("Not implemented yet")
    end
    @inbounds a0[2:end-1,end] .= BC.right.val
    if is_dirichlet(BC.right.t)
        @inbounds a1[2:end-1,end] .= -1.
        @inbounds b[2:end-1,end] .= 0.
    elseif is_neumann(BC.right.t)
        @inbounds a1[2:end-1,end] .= 0.
        @inbounds b[2:end-1,end] .= 1.
    elseif is_robin(BC.right.t)
        @inbounds a1[2:end-1,end] .= -1.
        @inbounds b[2:end-1,end] .= 1.
    elseif is_periodic(BC.right.t)
        # @inbounds a1[:,end] .= 0.
        # @inbounds b[:,end] .= 0.
    else
        @error ("Not implemented yet")
    end
    @inbounds a0[end,:] .= BC.top.val
    if is_dirichlet(BC.top.t)
        @inbounds a1[end,:] .= -1.
        @inbounds b[end,:] .= 0.
    elseif is_neumann(BC.top.t)
        @inbounds a1[end,:] .= 0.
        @inbounds b[end,:] .= 1.
    elseif is_robin(BC.top.t)
        @inbounds a1[end,:] .= -1.
        @inbounds b[end,:] .= 1.
    elseif is_periodic(BC.top.t)
        # @inbounds a1[end,2:end-1] .= 0.
        # @inbounds b[end,2:end-1] .= 0.
    else
        @error ("Not implemented yet")
    end

    return nothing
end

function set_borders!(grid, a0, a1, b0, b1, BC)
    @unpack ny, ind = grid
    
    @inbounds a0[:,1] .= BC.left.val
    if is_dirichlet(BC.left.t)
        @inbounds a1[2:end-1,1] .= -1.
        @inbounds b0[2:end-1,1] .= 0.
        @inbounds b1[2:end-1,1] .= 0.
    elseif is_neumann(BC.left.t)
        @inbounds a1[2:end-1,1] .= 0.
        @inbounds b0[2:end-1,1] .= 0.
        @inbounds b1[2:end-1,1] .= 1.
    elseif is_robin(BC.left.t)
        @inbounds a1[2:end-1,1] .= -1.
        @inbounds b0[2:end-1,1] .= 0.
        @inbounds b1[2:end-1,1] .= 1.
    elseif is_periodic(BC.left.t)
        # @inbounds a1[:,1] .= 0.
        # @inbounds b0[:,1] .= 0.
        # @inbounds b1[:,1] .= 0.
    else
        @error ("Not implemented yet")
    end
    @inbounds a0[1,2:end-1] .= BC.bottom.val
    if is_dirichlet(BC.bottom.t)
        @inbounds a1[1,:] .= -1.
        @inbounds b0[1,:] .= 0.
        @inbounds b1[1,:] .= 0.
    elseif is_neumann(BC.bottom.t)
        @inbounds a1[1,:] .= 0.
        @inbounds b0[1,:] .= 0.
        @inbounds b1[1,:] .= 1.
    elseif is_robin(BC.bottom.t)
        @inbounds a1[1,:] .= -1.
        @inbounds b0[1,:] .= 0.
        @inbounds b1[1,:] .= 1.
    elseif is_periodic(BC.bottom.t)
        # @inbounds a1[1,2:end-1] .= 0.
        # @inbounds b0[1,2:end-1] .= 0.
        # @inbounds b1[1,2:end-1] .= 0.
    else
        @error ("Not implemented yet")
    end
    @inbounds a0[:,end] .= BC.right.val
    if is_dirichlet(BC.right.t)
        @inbounds a1[2:end-1,end] .= -1.
        @inbounds b0[2:end-1,end] .= 0.
        @inbounds b1[2:end-1,end] .= 0.
    elseif is_neumann(BC.right.t)
        @inbounds a1[2:end-1,end] .= 0.
        @inbounds b0[2:end-1,end] .= 0.
        @inbounds b1[2:end-1,end] .= 1.
    elseif is_robin(BC.right.t)
        @inbounds a1[2:end-1,end] .= -1.
        @inbounds b0[2:end-1,end] .= 0.
        @inbounds b1[2:end-1,end] .= 1.
    elseif is_periodic(BC.right.t)
        # @inbounds a1[:,end] .= 0.
        # @inbounds b0[:,end] .= 0.
        # @inbounds b1[:,end] .= 0.
    else
        @error ("Not implemented yet")
    end
    @inbounds a0[end,2:end-1] .= BC.top.val
    if is_dirichlet(BC.top.t)
        @inbounds a1[end,:] .= -1.
        @inbounds b0[end,:] .= 0.
        @inbounds b1[end,:] .= 0.
    elseif is_neumann(BC.top.t)
        @inbounds a1[end,:] .= 0.
        @inbounds b0[end,:] .= 0.
        @inbounds b1[end,:] .= 1.
    elseif is_robin(BC.top.t)
        @inbounds a1[end,:] .= -1.
        @inbounds b0[end,:] .= 0.
        @inbounds b1[end,:] .= 1.
    elseif is_periodic(BC.top.t)
        # @inbounds a1[end,2:end-1] .= 0.
        # @inbounds b0[end,2:end-1] .= 0.
        # @inbounds b1[end,2:end-1] .= 0.
    else
        @error ("Not implemented yet")
    end

    return nothing
end

function update_dirichlet_field!(grid, bv, v, BC)
    tmp = zeros(size(v))

    if is_dirichlet(BC.left.t)
        @inbounds tmp[2:end-1,1] .= BC.left.val
    elseif is_neumann(BC.left.t)
        @inbounds tmp[2:end-1,1] .= v[2:end-1,1] .- BC.left.val .* 0.5 .* grid.dx[2:end-1,1]
    end

    if is_dirichlet(BC.bottom.t)
        @inbounds tmp[1,:] .= BC.bottom.val
    elseif is_neumann(BC.bottom.t)
        @inbounds tmp[1,:] .= v[1,:] .- BC.bottom.val .* 0.5 .* grid.dy[1,:]
    end

    if is_dirichlet(BC.right.t)
        @inbounds tmp[2:end-1,end] .= BC.right.val
    elseif is_neumann(BC.right.t)
        @inbounds tmp[2:end-1,end] .= BC.right.val .* 0.5 .* grid.dx[2:end-1,end] .+ v[2:end-1,end]
    end

    if is_dirichlet(BC.top.t)
        @inbounds tmp[end,:] .= BC.top.val
    elseif is_neumann(BC.top.t)
        @inbounds tmp[end,:] .= BC.top.val .* 0.5 .* grid.dy[end,:]  .+ v[end,:]
    end

    bv[grid.ny*grid.nx+1:end] .= vec(tmp)

    return nothing
end

function set_cutcell_matrices!(grid, geo, opC, periodic_x, periodic_y)
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
    divergence_A!(AxT, AyT, geo.dcap, ny, ind.all_indices)
    divergence_B!(BxT, ByT, geo.dcap, ny, ind.all_indices)

    mat_assign!(Bx, sparse(-BxT'))
    mat_assign!(By, sparse(-ByT'))

    # Matrices for BCs
    bc_matrix!(Hx, Hy, geo.dcap, ny, ind.all_indices)

    mat_assign_T!(HxT, sparse(Hx'))
    mat_assign_T!(HyT, sparse(Hy'))

    periodic_bcs!(grid, Bx, By, Hx, Hy, periodic_x, periodic_y)

    return nothing
end

function set_other_cutcell_matrices(grid, geo, geo_u, geo_v, opC_p, opC_u, opC_v, periodic_x, periodic_y)
    @unpack nx, ny, ind = grid
    @unpack Bx, By, Gx, Gy = opC_p

    bc_matrix!(opC_u.Gx, opC_v.Gy, geo.dcap, geo_u.dcap, geo_v.dcap, ny, ind.all_indices)

    mat_assign_T!(Gx, sparse(opC_u.Gx'))
    mat_assign_T!(Gy, sparse(opC_v.Gy'))

    periodic_bcs!(grid, Bx, By, opC_u.Gx, opC_v.Gy, periodic_x, periodic_y)
    periodic_bcs_R!(grid, opC_u.Rx, opC_v.Ry, periodic_x, periodic_y)
end

function laplacian(opC)
    @unpack Bx, By, BxT, ByT, Hx, Hy, HxT, HyT, iMx, iMy, tmp_x, tmp_y = opC

    mul!(tmp_x, iMx, Bx)
    L = BxT * tmp_x
    mul!(tmp_y, iMy, By)
    L = L .+ ByT * tmp_y
    # mul!(L, ByT, tmp, 1.0, 1.0) # BxT * iMx * Bx .+ ByT * iMy * By

    mul!(tmp_x, iMx, Hx)
    bc_L = BxT * tmp_x
    mul!(tmp_y, iMy, Hy)
    bc_L = bc_L .+ ByT * tmp_y
    # mul!(bc_L, ByT, tmp, 1.0, 1.0) # BxT * iMx * Hx .+ ByT * iMy * Hy

    return L, bc_L
end

function laplacian_fs(opC)
    @unpack AxT, AyT, Bx, By, Hx, Hy, iMx, iMy, tmp_x, tmp_y = opC

    mul!(tmp_x, iMx, Bx)
    L = AxT * tmp_x
    mul!(tmp_y, iMy, By)
    L = L .+ AyT * tmp_y
    # mul!(L, ByT, tmp, 1.0, 1.0) # BxT * iMx * Bx .+ ByT * iMy * By

    mul!(tmp_x, iMx, Hx)
    bc_L = AxT * tmp_x
    mul!(tmp_y, iMy, Hy)
    bc_L = bc_L .+ AyT * tmp_y
    # mul!(bc_L, ByT, tmp, 1.0, 1.0) # BxT * iMx * Hx .+ ByT * iMy * Hy

    return L, bc_L
end

function set_laplacians!(grid, geo, grid_u, geo_u, grid_v, geo_v,
                         opC_p, opC_u, opC_v,
                         periodic_x, periodic_y, fs=false
    )
    @unpack ny, ind = grid

    set_cutcell_matrices!(grid, geo, opC_p, periodic_x, periodic_y)
    set_cutcell_matrices!(grid_u, geo_u, opC_u, periodic_x, periodic_y)
    set_cutcell_matrices!(grid_v, geo_v, opC_v, periodic_x, periodic_y)

    set_other_cutcell_matrices(grid, geo, geo_u, geo_v, opC_p, opC_u, opC_v, periodic_x, periodic_y)

    Lp, bc_Lp = laplacian(opC_p)
    Lu, bc_Lu = laplacian(opC_u)
    Lv, bc_Lv = laplacian(opC_v)

    if fs
        Lp_fs, bc_Lp_fs = laplacian_fs(opC_p)
        Lu_fs, bc_Lu_fs = laplacian_fs(opC_u)
        Lv_fs, bc_Lv_fs = laplacian_fs(opC_v)

        return Lp, bc_Lp, Lu, bc_Lu, Lv, bc_Lv, Lp_fs, bc_Lp_fs, Lu_fs, bc_Lu_fs, Lv_fs, bc_Lv_fs
    else
        return Lp, bc_Lp, Lu, bc_Lu, Lv, bc_Lv
    end
end

function set_crank_nicolson_block(bc_type, num, grid, opC, L, bc_L, Lm1, bc_Lm1, Mm1, BC)
    @unpack τ = num
    @unpack nx, ny = grid
    @unpack Bx, By, Hx, Hy, HxT, HyT, χ, M, iMx, iMy = opC

    if bc_type == dir
        __a1 = -1.
        __b = 0.
    elseif bc_type == neu
        __a1 = 0.
        __b = 1.
    elseif bc_type == rob
        __a1 = -1.
        __b = 1.
    end

    a0 = ones(ny, nx) .* grid.V
    _a1 = ones(ny, nx) .* __a1
    _b = ones(ny, nx) .* __b
    set_borders!(grid, a0, _a1, _b, BC)
    a1 = Diagonal(vec(_a1))
    b = Diagonal(vec(_b))

    τ_2 = 0.5 * τ

    data_A = Matrix{SparseMatrixCSC{Float64, Int64}}(undef, 2, 2)
    data_A[1,1] = pad_crank_nicolson(M .- τ_2 .* L, grid, τ)
    data_A[1,2] = - τ_2 .* bc_L
    data_A[2,1] = b * (HxT * iMx * Bx .+ HyT * iMy * By)
    data_A[2,2] = pad(b * (HxT * iMx * Hx .+ HyT * iMy * Hy) .- χ * a1)
    A  = [data_A[1,1] data_A[1,2];
          data_A[2,1] data_A[2,2]]

    data_B = Matrix{SparseMatrixCSC{Float64, Int64}}(undef, 2, 2)
    data_B[1,1] = Mm1 .+ τ_2 .* Lm1
    data_B[1,2] = τ_2 .* bc_Lm1
    data_B[2,1] = spdiagm(0 => zeros(nx*ny))
    data_B[2,2] = spdiagm(0 => zeros(nx*ny))
    B  = [data_B[1,1] data_B[1,2];
          data_B[2,1] data_B[2,2]]

    rhs = zeros(2*nx*ny)
    rhs[nx*ny+1:end] .= χ * vec(a0)
    
    return A, B, rhs
end

function set_poisson_block(bc_type, grid, a0, opC, L, bc_L, BC)
    @unpack nx, ny = grid
    @unpack Bx, By, Hx, Hy, HxT, HyT, χ, M, iMx, iMy = opC

    if bc_type == dir
        __a1 = -1.
        __b = 0.
    elseif bc_type == neu
        __a1 = 0.
        __b = 1.
    elseif bc_type == rob
        __a1 = -1.
        __b = 1.
    end

    _a1 = ones(ny, nx) .* __a1
    _b = ones(ny, nx) .* __b
    set_borders!(grid, a0, _a1, _b, BC)
    a1 = Diagonal(vec(_a1))
    b = Diagonal(vec(_b))

    data_A = Matrix{SparseMatrixCSC{Float64, Int64}}(undef, 2, 2)
    data_A[1,1] = pad(L, -4.0)
    data_A[1,2] = bc_L
    data_A[2,1] = b * (HxT * iMx * Bx .+ HyT * iMy * By)
    data_A[2,2] = pad(b * (HxT * iMx * Hx .+ HyT * iMy * Hy) .- χ * a1)
    A  = [data_A[1,1] data_A[1,2];
          data_A[2,1] data_A[2,2]]

    rhs = zeros(2*nx*ny)
    rhs[nx*ny+1:end] .= χ * vec(a0)
    
    return A, rhs
end

function projection_no_slip!(num, grid, geo, grid_u, geo_u, grid_v, geo_v, ph,
                            BC_u, BC_v, BC_p,
                            opC_p, opC_u, opC_v,
                            Lum1, bc_Lum1, Lvm1, bc_Lvm1, Mum1, Mvm1,
                            FULL, MIXED, periodic_x, periodic_y
    )
    @unpack Re, τ = num
    @unpack p, pD, ϕ, ϕD, Gxm1, Gym1, u, uD, v, vD, ucorr, vcorr = ph

    iRe = 1.0 / Re
    iτ = 1.0 / τ

    Lp, bc_Lp, Lu, bc_Lu, Lv, bc_Lv = set_laplacians!(grid, geo, grid_u, geo_u, grid_v, geo_v,
                                                      opC_p, opC_u, opC_v,
                                                      periodic_x, periodic_y)

    A_u, B_u, rhs_u = set_crank_nicolson_block(dir, num, grid_u, opC_u, iRe.*Lu, iRe.*bc_Lu, iRe.*Lum1, iRe.*bc_Lum1, Mum1, BC_u)
    A_v, B_v, rhs_v = set_crank_nicolson_block(dir, num, grid_v, opC_v, iRe.*Lv, iRe.*bc_Lv, iRe.*Lvm1, iRe.*bc_Lvm1, Mvm1, BC_v)
    a0_p = zeros(grid.ny, grid.nx)
    A_p, rhs_p = set_poisson_block(neu, grid, a0_p, opC_p, Lp, bc_Lp, BC_p)

    uD[1:grid_u.ny*grid_u.nx] .= vec(u)
    update_dirichlet_field!(grid_u, uD, u, BC_u)
    mul!(rhs_u, B_u, uD, 1.0, 1.0)
    rhs_u[1:grid_u.ny*grid_u.nx] .+= -τ .* (opC_p.Bx * pD[1:grid.ny*grid.nx] .+ opC_p.Hx * pD[grid.ny*grid.nx+1:end])
    blocks = DDM.decompose(A_u, grid_u.domdec, grid_u.domdec)
    @mytime _, ch = bicgstabl!(uD, A_u, rhs_u, Pl=ras(blocks,grid_u.pou), log=true)
    println(ch)
    ucorr .= reshape(uD[1:grid_u.ny*grid_u.nx], (grid_u.ny, grid_u.nx))

    vD[1:grid_v.ny*grid_v.nx] .= vec(v)
    update_dirichlet_field!(grid_v, vD, v, BC_v)
    mul!(rhs_v, B_v, vD, 1.0, 1.0)
    rhs_v[1:grid_v.ny*grid_v.nx] .+= -τ .* (opC_p.By * pD[1:grid.ny*grid.nx] .+ opC_p.Hy * pD[grid.ny*grid.nx+1:end])
    blocks = DDM.decompose(A_v, grid_v.domdec, grid_v.domdec)
    @mytime _, ch = bicgstabl!(vD, A_v, rhs_v, Pl=ras(blocks,grid_v.pou), log=true)
    println(ch)
    vcorr .= reshape(vD[1:grid_v.ny*grid_v.nx], (grid_v.ny, grid_v.nx))

    Duv = opC_p.AxT * uD[1:grid_u.ny*grid_u.nx] .+ opC_p.HxT * uD[grid_u.ny*grid_u.nx+1:end] .+
          opC_p.AyT * vD[1:grid_v.ny*grid_v.nx] .+ opC_p.HyT * vD[grid_v.ny*grid_v.nx+1:end]
    rhs_p[1:grid.ny*grid.nx] .=  iτ .* Duv

    sum_rhs = sum(rhs_p[1:grid.ny*grid.nx])
    sum_Mp = sum(geo.dcap[:,:,5])
    non_empty = vcat(FULL, MIXED)
    @inbounds @threads for II in non_empty
        pII = lexicographic(II, grid.ny)
        @inbounds rhs_p[1:grid.ny*grid.nx][pII] -= sum_rhs * geo.dcap[II,5] / sum_Mp
    end

    ϕD[1:grid.ny*grid.nx] .= 0.
    ϕD[grid.ny*grid.nx+1:end] .= 0.
    blocks = DDM.decompose(A_p, grid.domdec, grid.domdec)
    @mytime _, ch = bicgstabl!(ϕD, A_p, rhs_p, Pl=ras(blocks,grid.pou), log=true)
    println(ch)
    ϕ .= reshape(ϕD[1:grid.ny*grid.nx], (grid.ny, grid.nx))

    iM = Diagonal(1. ./ (vec(geo.dcap[:,:,5]) .+ eps(0.01)))
    pD[1:grid.ny*grid.nx] .= vec(p)
    pD[1:grid.ny*grid.nx] .+= vec(ϕ) .- 0.5 .* iRe .* iM * Duv
    pD[grid.ny*grid.nx+1:end] .+= ϕD[grid.ny*grid.nx+1:end]
    p .= reshape(pD[1:grid.ny*grid.nx], (grid.ny, grid.nx))

    u .= ucorr .- τ .* reshape(opC_p.iMx * (opC_p.Bx * vec(ϕ) .+ opC_p.Hx * ϕD[grid.ny*grid.nx+1:end]), (grid_u.ny, grid_u.nx))
    v .= vcorr .- τ .* reshape(opC_p.iMy * (opC_p.By * vec(ϕ) .+ opC_p.Hy * ϕD[grid.ny*grid.nx+1:end]), (grid_v.ny, grid_v.nx))

    return Lu, bc_Lu, Lv, bc_Lv, opC_u.M, opC_v.M
end

function set_crank_nicolson_block(bc_type, num, 
    grid, opC_p, Lp, bc_Lp, Lp_fs, bc_Lp_fs, BC_p,
    grid_u, opC_u, Lu, bc_Lu, Lu_fs, bc_Lu_fs, BC_u,
    grid_v, opC_v, Lv, bc_Lv, Lv_fs, bc_Lv_fs, BC_v)
    @unpack Re, τ = num

    iRe = 1 / Re

    if bc_type == dir
        __a1 = -1.
        __b0 = 0.
        __b1 = 0.
    elseif bc_type == neu
        __a1 = 0.
        # __b0 = 1.
        # __b1 = 0.
        __b0 = 0.
        __b1 = 1.
    elseif bc_type == rob
        __a1 = -1.
        __b0 = 1.
        __b1 = 0.
    end

    a0_u = zeros(grid_u.ny, grid_u.nx)
    _a1_u = ones(grid_u.ny, grid_u.nx) .* __a1
    _b0_u = ones(grid_u.ny, grid_u.nx) .* __b0
    _b1_u = ones(grid_u.ny, grid_u.nx) .* __b1
    a0_v = zeros(grid_v.ny, grid_v.nx)
    _a1_v = ones(grid_v.ny, grid_v.nx) .* __a1
    _b0_v = ones(grid_v.ny, grid_v.nx) .* __b0
    _b1_v = ones(grid_v.ny, grid_v.nx) .* __b1
    a0_p = zeros(grid.ny, grid.nx)
    _a1_p = ones(grid.ny, grid.nx) .* (-1.)
    _b0_p = ones(grid.ny, grid.nx) .* 0.
    _b1_p = ones(grid.ny, grid.nx) .* 0.
    set_borders!(grid_u, a0_u, _a1_u, _b0_u, _b1_u, BC_u)
    set_borders!(grid_v, a0_v, _a1_v, _b0_v, _b1_v, BC_v)
    set_borders!(grid, a0_p, _a1_p, _b0_p, _b1_p, BC_p)
    a1_u = Diagonal(vec(_a1_u))
    b0_u = Diagonal(vec(_b0_u))
    b1_u = Diagonal(vec(_b1_u))
    a1_v = Diagonal(vec(_a1_v))
    b0_v = Diagonal(vec(_b0_v))
    b1_v = Diagonal(vec(_b1_v))
    a1_p = Diagonal(vec(_a1_p))
    b0_p = Diagonal(vec(_b0_p))
    b1_p = Diagonal(vec(_b1_p))

    τ_2 = 0.5 * τ

    if grid_u.nx*grid_u.ny < grid_v.nx*grid_v.ny
        size_zeros = grid_u.nx*grid_u.ny
    else
        size_zeros = grid_v.nx*grid_v.ny
    end

    data_A = Matrix{SparseMatrixCSC{Float64, Int64}}(undef, 6, 6)
    data_A[1,1] = pad_crank_nicolson(opC_u.M .- τ .* Lu, grid_u, τ)
    data_A[1,2] = - τ .* bc_Lu

    data_A[2,1] = #=-b0_u * opC_u.χ .+=#
                  (b0_u .+ b1_u) * (opC_u.HxT * opC_u.iMx * opC_u.Bx .+
                          opC_u.HyT * opC_u.iMy * opC_u.By)
    data_A[2,2] = pad(#=b0_u * opC_u.χ .+=#
                      (b0_u .+ b1_u) * (opC_u.HxT * opC_u.iMx * opC_u.Hx .+
                              opC_u.HyT * opC_u.iMy * opC_u.Hy) .-
                      opC_u.χ * a1_u)

    data_A[3,3] = pad_crank_nicolson(opC_v.M .- τ .* Lv, grid_v, τ)
    data_A[3,4] = - τ .* bc_Lv

    data_A[4,3] = #=-b0_v * opC_v.χ .+=#
                  (b0_v .+ b1_v) * (opC_v.HyT * opC_v.iMy * opC_v.By .+
                          opC_v.HxT * opC_v.iMx * opC_v.Bx)
    data_A[4,4] = pad(#=b0_v * opC_v.χ .+=#
                      (b0_v .+ b1_v) * (opC_v.HyT * opC_v.iMy * opC_v.Hy .+
                              opC_v.HxT * opC_v.iMx * opC_v.Hx) .-
                      opC_v.χ * a1_v)

    data_A[5,5] = pad(Lp, -4.0)
    data_A[5,6] = bc_Lp

    # data_A[6,5] = b0_p * opC_p.χ * Lp .+
    #               b1_p * (opC_p.HxT * opC_p.iMx * opC_p.Bx .+
    #                       opC_p.HyT * opC_p.iMy * opC_p.By)
    # data_A[6,6] = pad(b0_p * opC_p.χ * bc_Lp .+
    #                   b1_p * (opC_p.HxT * opC_p.iMx * opC_p.Hx .+
    #                           opC_p.HyT * opC_p.iMy * opC_p.Hy) .-
    #                   opC_p.χ * a1_p)
    GxT = opC_u.Gx'
    GyT = opC_v.Gy'
    data_A[6,5] = b1_p * (opC_p.HxT * opC_p.iMx * opC_p.Bx .+
                          opC_p.HyT * opC_p.iMy * opC_p.By)
    data_A[6,6] = pad(b0_p * (GxT * opC_u.Gx .+ GyT * opC_v.Gy) .+
                      b1_p * (opC_p.HxT * opC_p.iMx * opC_p.Hx .+
                              opC_p.HyT * opC_p.iMy * opC_p.Hy) .-
                      opC_p.χ * a1_p)

    data_rhs = Vector{Vector{Float64}}(undef, 6)
    data_rhs[1] = zeros(grid_u.nx*grid_u.ny)
    data_rhs[2] = opC_u.χ * vec(a0_u)
    data_rhs[3] = zeros(grid_v.nx*grid_v.ny)
    data_rhs[4] = opC_v.χ * vec(a0_v)
    data_rhs[5] = zeros(grid.nx*grid.ny)
    data_rhs[6] = opC_p.χ * vec(a0_p)
    
    return data_A, data_rhs
end

function projection_fs!(num, grid, geo, grid_u, geo_u, grid_v, geo_v, ph,
                        BC_u, BC_v, BC_p,
                        opC_p, opC_u, opC_v,
                        Lum1, bc_Lum1, Lvm1, bc_Lvm1, Mum1, Mvm1,
                        FRESH, FRESH_u, FRESH_v,
                        FULL, MIXED, periodic_x, periodic_y, current_i
    )
    @unpack Re, τ, σ, g, β = num
    @unpack p, pD, ϕ, ϕD, Gxm1, Gym1, u, v, uD, vD, uvD, uvϕD, ucorr, vcorr = ph

    iRe = 1.0 / Re
    iτ = 1.0 / τ

    opC_u_m1 = deepcopy(opC_u)
    opC_v_m1 = deepcopy(opC_v)

    Lp, bc_Lp, Lu, bc_Lu, Lv, bc_Lv, Lp_fs, bc_Lp_fs, Lu_fs, bc_Lu_fs, Lv_fs, bc_Lv_fs = set_laplacians!(grid, geo, grid_u, geo_u, grid_v, geo_v,
        opC_p, opC_u, opC_v,
        periodic_x, periodic_y, true)

    A_uvϕ, rhs_uvϕ = set_crank_nicolson_block(neu, num,
        grid, opC_p, Lp, bc_Lp, Lp_fs, bc_Lp_fs, BC_p,
        grid_u, opC_u, iRe.*Lu, iRe.*bc_Lu, iRe.*Lu_fs, iRe.*bc_Lu_fs, BC_u,
        grid_v, opC_v, iRe.*Lv, iRe.*bc_Lv, iRe.*Lv_fs, iRe.*bc_Lv_fs, BC_v)

    a0_u = zeros(grid_u.ny, grid_u.nx)
    _a1_u = zeros(grid_u.ny, grid_u.nx)
    _b0_u = ones(grid_u.ny, grid_u.nx)
    _b1_u = zeros(grid_u.ny, grid_u.nx)
    a0_v = zeros(grid_v.ny, grid_v.nx)
    _a1_v = zeros(grid_v.ny, grid_v.nx)
    _b0_v = ones(grid_v.ny, grid_v.nx)
    _b1_v = zeros(grid_v.ny, grid_v.nx)
    a0_p = zeros(grid.ny, grid.nx)
    _a1_p = zeros(grid.ny, grid.nx)
    _b0_p = ones(grid.ny, grid.nx)
    _b1_p = zeros(grid.ny, grid.nx)
    set_borders!(grid_u, a0_u, _a1_u, _b0_u, _b1_u, BC_u)
    set_borders!(grid_v, a0_v, _a1_v, _b0_v, _b1_v, BC_v)
    set_borders!(grid, a0_p, _a1_p, _b0_p, _b1_p, BC_p)
    b0_u = Diagonal(vec(_b0_u))
    b0_v = Diagonal(vec(_b0_v))
    b0_p = Diagonal(vec(_b0_p))

    no_fs_u = I - b0_u
    no_fs_v = I - b0_v

    rhs_uvϕ[1] .+= Mum1 * vec(u)
    rhs_uvϕ[3] .+= Mvm1 * vec(v)

    grav_x = g .* sin(β) .* Mum1 * ones(grid_u.nx * grid_u.ny)
    grav_y = g .* cos(β) .* Mvm1 * ones(grid_v.nx * grid_v.ny)

    rhs_uvϕ[1] .+= τ .* grav_x
    rhs_uvϕ[3] .+= - τ .* grav_y

    kill_dead_cells!(rhs_uvϕ[1], grid_u, geo_u)
    kill_dead_cells!(rhs_uvϕ[2], grid_u, geo_u)
    kill_dead_cells!(rhs_uvϕ[3], grid_v, geo_v)
    kill_dead_cells!(rhs_uvϕ[4], grid_v, geo_v)

    Au = [A_uvϕ[1,1] A_uvϕ[1,2];
          A_uvϕ[2,1] A_uvϕ[2,2]]
    rhs_u = vcat(rhs_uvϕ[1], rhs_uvϕ[2])
    blocks = DDM.decompose(Au, grid_u.domdec, grid_u.domdec)
    @mytime _, ch = bicgstabl!(uD, Au, rhs_u, Pl=ras(blocks,grid_u.pou), log=true)
    println(ch)

    Av = [A_uvϕ[3,3] A_uvϕ[3,4];
          A_uvϕ[4,3] A_uvϕ[4,4]]
    rhs_v = vcat(rhs_uvϕ[3], rhs_uvϕ[4])
    blocks = DDM.decompose(Av, grid_v.domdec, grid_v.domdec)
    @mytime _, ch = bicgstabl!(vD, Av, rhs_v, Pl=ras(blocks,grid_v.pou), log=true)
    println(ch)

    rhs_ϕ = vcat(rhs_uvϕ[5], rhs_uvϕ[6])

    Duv = opC_p.AxT * uD[1:grid_u.ny*grid_u.nx] .+ opC_p.Gx * uD[grid_u.ny*grid_u.nx+1:end] .+
          opC_p.AyT * vD[1:grid_v.ny*grid_v.nx] .+ opC_p.Gy * vD[grid_v.ny*grid_v.nx+1:end]
    rhs_ϕ[1:grid.ny*grid.nx] .= iτ .* Duv

    GxT = opC_u.Gx'
    GyT = opC_v.Gy'
    S = GxT * ((2 .* opC_u.HxT * opC_u.iMx * opC_u.Bx .+ opC_u.HyT * opC_u.iMy * opC_u.By) * uD[1:grid_u.ny*grid_u.nx] .+
               (2 .* opC_u.HxT * opC_u.iMx * opC_u.Hx .+ opC_u.HyT * opC_u.iMy * opC_u.Hy) * uD[grid_u.ny*grid_u.nx+1:end]) .+
        GyT * ((2 .* opC_v.HyT * opC_v.iMy * opC_v.By .+ opC_v.HxT * opC_v.iMx * opC_v.Bx) * vD[1:grid_v.ny*grid_v.nx] .+
               (2 .* opC_v.HyT * opC_v.iMy * opC_v.Hy .+ opC_v.HxT * opC_v.iMx * opC_v.Hx) * vD[grid_v.ny*grid_v.nx+1:end])

    rhs_ϕ[grid.ny*grid.nx+1:end] .= b0_p * (iRe .* S .- σ .* (GxT * opC_u.Gx .+ GyT * opC_v.Gy) * vec(grid.κ))
    Aϕ = [A_uvϕ[5,5] A_uvϕ[5,6];
          A_uvϕ[6,5] A_uvϕ[6,6]]
    blocks = DDM.decompose(Aϕ, grid.domdec, grid.domdec)
    @mytime _, ch = bicgstabl!(ϕD, Aϕ, rhs_ϕ, Pl=ras(blocks,grid.pou), log=true)
    println(ch)

    ucorr .= reshape(uD[1:grid_u.ny*grid_u.nx], (grid_u.ny, grid_u.nx))
    vcorr .= reshape(vD[1:grid_v.ny*grid_v.nx], (grid_v.ny, grid_v.nx))
    ϕ .= reshape(ϕD[1:grid.ny*grid.nx], (grid.ny, grid.nx))

    iMu = Diagonal(1 ./ (opC_u.M.diag .+ eps(0.01)))
    iMv = Diagonal(1 ./ (opC_v.M.diag .+ eps(0.01)))

    ∇ϕ_x = opC_u.AxT * opC_u.Rx * ϕD[1:grid.ny*grid.nx] .+ opC_u.Gx * ϕD[grid.ny*grid.nx+1:end]
    ∇ϕ_y = opC_v.AyT * opC_v.Ry * ϕD[1:grid.ny*grid.nx] .+ opC_v.Gy * ϕD[grid.ny*grid.nx+1:end]

    p .= ϕ
    Gxm1 .= ∇ϕ_x
    Gym1 .= ∇ϕ_y

    kill_dead_cells!(Gxm1, grid_u, geo_u)
    kill_dead_cells!(Gym1, grid_v, geo_v)

    u .= ucorr .- τ .* reshape(iMu * ∇ϕ_x, (grid_u.ny, grid_u.nx))
    v .= vcorr .- τ .* reshape(iMv * ∇ϕ_y, (grid_v.ny, grid_v.nx))

    return Lu, bc_Lu, Lv, bc_Lv, opC_p.M, opC_u.M, opC_v.M
end