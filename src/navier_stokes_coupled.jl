function set_borders!(grid, num, a0, a1, b, BC)
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
    elseif is_navier(BC.left.t)
        @inbounds a1[2:end-1,1] .= -1.
        @inbounds b[2:end-1,1] .= num.λCA
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
    elseif is_navier(BC.bottom.t)
        @inbounds a1[1,:] .= -1.
        @inbounds b[1,:] .= num.λCA
        @inbounds a0[1,:] .= BC.bottom.val .+ compute_young_stress(grid, num, grid.ind.b_bottom[1])
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
    elseif is_navier(BC.right.t)
        @inbounds a1[2:end-1,end] .= -1.
        @inbounds b[2:end-1,end] .= num.λCA
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
    elseif is_navier(BC.top.t)
        @inbounds a1[end,:] .= -1.
        @inbounds b[end,:] .= num.λCA
    else
        @error ("Not implemented yet")
    end

    return nothing
end

function set_borders!(grid, num, a0, a1, b0, b1, BC)
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
    elseif is_navier(BC.left.t)
        @inbounds a1[2:end-1,1] .= -1.
        @inbounds b0[2:end-1,1] .= 0.
        @inbounds b1[2:end-1,1] .= num.λCA
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
    elseif is_navier(BC.bottom.t)
        @inbounds a1[1,:] .= -1.
        @inbounds b0[1,:] .= 0.
        @inbounds b1[1,:] .= num.λCA
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
    elseif is_navier(BC.right.t)
        @inbounds a1[2:end-1,end] .= -1.
        @inbounds b0[2:end-1,end] .= 0.
        @inbounds b1[2:end-1,end] .= num.λCA
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
    elseif is_navier(BC.top.t)
        @inbounds a1[end,:] .= -1.
        @inbounds b0[end,:] .= 0.
        @inbounds b1[end,:] .= num.λCA
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

function strain_rate(opC_u, opC_v)
    GxT = opC_u.Gx'
    GyT = opC_v.Gy'

    data = Matrix{SparseMatrixCSC{Float64, Int64}}(undef, 2, 2)
    data[1,1] = GxT * (2 .* opC_u.HxT * opC_u.iMx * opC_u.Bx .+ opC_u.HyT * opC_u.iMy * opC_u.By)
    data[1,2] = GxT * (2 .* opC_u.HxT * opC_u.iMx * opC_u.Hx .+ opC_u.HyT * opC_u.iMy * opC_u.Hy)
    data[2,1] = GyT * (2 .* opC_v.HyT * opC_v.iMy * opC_v.By .+ opC_v.HxT * opC_v.iMx * opC_v.Bx)
    data[2,2] = GyT * (2 .* opC_v.HyT * opC_v.iMy * opC_v.Hy .+ opC_v.HxT * opC_v.iMx * opC_v.Hx)

    return data
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
    set_borders!(grid, num, a0, _a1, _b, BC)
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

function set_poisson_block(bc_type, grid, num, a0, opC, L, bc_L, BC)
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
    set_borders!(grid, num, a0, _a1, _b, BC)
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
    @unpack p, pD, ϕ, ϕD, Gxm1, Gym1, u, uD, v, vD, ucorr, vcorr, ucorrD, vcorrD = ph

    iRe = 1.0 / Re
    iτ = 1.0 / τ

    Lp, bc_Lp, Lu, bc_Lu, Lv, bc_Lv = set_laplacians!(grid, geo, grid_u, geo_u, grid_v, geo_v,
                                                      opC_p, opC_u, opC_v,
                                                      periodic_x, periodic_y)

    A_u, B_u, rhs_u = set_crank_nicolson_block(dir, num, grid_u, opC_u, iRe.*Lu, iRe.*bc_Lu, iRe.*Lum1, iRe.*bc_Lum1, Mum1, BC_u)
    A_v, B_v, rhs_v = set_crank_nicolson_block(dir, num, grid_v, opC_v, iRe.*Lv, iRe.*bc_Lv, iRe.*Lvm1, iRe.*bc_Lvm1, Mvm1, BC_v)
    a0_p = zeros(grid.ny, grid.nx)
    A_p, rhs_p = set_poisson_block(neu, grid, num, a0_p, opC_p, Lp, bc_Lp, BC_p)

    update_dirichlet_field!(grid_u, uD, u, BC_u)
    mul!(rhs_u, B_u, uD, 1.0, 1.0)
    veci(rhs_u,grid_u,1) .+= -τ .* (opC_p.Bx * veci(pD,grid,1) .+ opC_p.Hx * veci(pD,grid,2))
    blocks = DDM.decompose(A_u, grid_u.domdec, grid_u.domdec)
     _, ch = bicgstabl!(ucorrD, A_u, rhs_u, Pl=ras(blocks,grid_u.pou), log=false,verbose=false)
    #println(ch)
    ucorr .= reshape(veci(ucorrD,grid_u,1), (grid_u.ny, grid_u.nx))

    veci(vD,grid_v,1) .= vec(v)
    update_dirichlet_field!(grid_v, vD, v, BC_v)
    mul!(rhs_v, B_v, vD, 1.0, 1.0)
    veci(rhs_v,grid_v,1) .+= -τ .* (opC_p.By * veci(pD,grid,1) .+ opC_p.Hy * veci(pD,grid,2))
    blocks = DDM.decompose(A_v, grid_v.domdec, grid_v.domdec)
     _, ch = bicgstabl!(vcorrD, A_v, rhs_v, Pl=ras(blocks,grid_v.pou), log=false,verbose=false)
    #println(ch)
    vcorr .= reshape(veci(vcorrD,grid_v,1), (grid_v.ny, grid_v.nx))

    Duv = opC_p.AxT * veci(ucorrD,grid_u,1) .+ opC_p.HxT * veci(ucorrD,grid_u,2) .+
          opC_p.AyT * veci(vcorrD,grid_v,1) .+ opC_p.HyT * veci(vcorrD,grid_v,2)
    veci(rhs_p,grid,1) .=  iτ .* Duv

    sum_rhs = sum(veci(rhs_p,grid,1))
    sum_Mp = sum(geo.dcap[:,:,5])
    non_empty = vcat(FULL, MIXED)
    @inbounds @threads for II in non_empty
        pII = lexicographic(II, grid.ny)
        @inbounds veci(rhs_p,grid,1)[pII] -= sum_rhs * geo.dcap[II,5] / sum_Mp
    end

    veci(ϕD,grid,1) .= 0.
    veci(ϕD,grid,2) .= 0.
    blocks = DDM.decompose(A_p, grid.domdec, grid.domdec)
     _, ch = bicgstabl!(ϕD, A_p, rhs_p, Pl=ras(blocks,grid.pou), log=false,verbose=false)
    #println(ch)
    ϕ .= reshape(veci(ϕD,grid,1), (grid.ny, grid.nx))

    iM = Diagonal(1. ./ (vec(geo.dcap[:,:,5]) .+ eps(0.01)))
    veci(pD,grid,1) .= vec(p)
    veci(pD,grid,1) .+= vec(ϕ) .- 0.5 .* iRe .* iM * Duv
    veci(pD,grid,2) .+= veci(ϕD,grid,2)
    p .= reshape(veci(pD,grid,1), (grid.ny, grid.nx))

    u .= ucorr .- τ .* reshape(opC_p.iMx * (opC_p.Bx * vec(ϕ) .+ opC_p.Hx * veci(ϕD,grid,2)), (grid_u.ny, grid_u.nx))
    v .= vcorr .- τ .* reshape(opC_p.iMy * (opC_p.By * vec(ϕ) .+ opC_p.Hy * veci(ϕD,grid,2)), (grid_v.ny, grid_v.nx))

    veci(uD,grid_u,1) .= vec(u)
    veci(uD,grid_u,2) .= veci(ucorrD,grid_u,2)
    veci(vD,grid_v,1) .= vec(v)
    veci(vD,grid_v,2) .= veci(vcorrD,grid_v,2)

    return Lu, bc_Lu, Lv, bc_Lv, opC_u.M, opC_v.M
end

function set_crank_nicolson_block(bc_type, num, 
    grid, opC_p, Lp, bc_Lp, Lp_fs, bc_Lp_fs, BC_p,
    grid_u, opC_u, Lu, bc_Lu, Mum1, BC_u,
    grid_v, opC_v, Lv, bc_Lv, Mvm1, BC_v)
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

    a0_u = zeros(grid_u.ny, grid_u.nx) # U + stress curvature here
    _a1_u = ones(grid_u.ny, grid_u.nx) .* __a1
    _b0_u = ones(grid_u.ny, grid_u.nx) .* __b0
    _b1_u = ones(grid_u.ny, grid_u.nx) .* __b1
    a0_v = zeros(grid_v.ny, grid_v.nx)  # + stress
    _a1_v = ones(grid_v.ny, grid_v.nx) .* __a1
    _b0_v = ones(grid_v.ny, grid_v.nx) .* __b0
    _b1_v = ones(grid_v.ny, grid_v.nx) .* __b1
    a0_p = zeros(grid.ny, grid.nx)
    _a1_p = ones(grid.ny, grid.nx) .* (-1.)
    _b0_p = ones(grid.ny, grid.nx) .* 0.
    _b1_p = ones(grid.ny, grid.nx) .* 0.
    set_borders!(grid_u, num, a0_u, _a1_u, _b0_u, _b1_u, BC_u)
    set_borders!(grid_v, num, a0_v, _a1_v, _b0_v, _b1_v, BC_v)
    set_borders!(grid, num, a0_p, _a1_p, _b0_p, _b1_p, BC_p)
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

    data_A = Matrix{SparseMatrixCSC{Float64, Int64}}(undef, 2, 2)
    data_A[1,1] = pad_crank_nicolson(opC_u.M .- τ_2 .* Lu, grid_u, τ)
    data_A[1,2] = - τ_2 .* bc_Lu
    data_A[2,1] = #=-b0_u * opC_u.χ .+=#
                  (b0_u .+ b1_u) * (opC_u.HxT * opC_u.iMx * opC_u.Bx .+
                          opC_u.HyT * opC_u.iMy * opC_u.By)
    data_A[2,2] = pad(#=b0_u * opC_u.χ .+=#
                      (b0_u .+ b1_u) * (opC_u.HxT * opC_u.iMx * opC_u.Hx .+
                              opC_u.HyT * opC_u.iMy * opC_u.Hy) .-
                      opC_u.χ * a1_u)
    Au = [data_A[1,1] data_A[1,2];
          data_A[2,1] data_A[2,2]]

    data_B = Matrix{SparseMatrixCSC{Float64, Int64}}(undef, 2, 2)
    data_B[1,1] = Mum1 #.+ τ_2 .* Lum1

    data_B[1,2] = spdiagm(0 => zeros(grid_u.nx*grid_u.ny))
    data_B[2,1] = spdiagm(0 => zeros(grid_u.nx*grid_u.ny))
    data_B[2,2] = spdiagm(0 => zeros(grid_u.nx*grid_u.ny))
    Bu = [data_B[1,1] data_B[1,2];
          data_B[2,1] data_B[2,2]]

    data_A[1,1] = pad_crank_nicolson(opC_v.M .- τ_2 .* Lv, grid_v, τ)
    data_A[1,2] = - τ_2 .* bc_Lv
    data_A[2,1] = #=-b0_v * opC_v.χ .+=#
                  (b0_v .+ b1_v) * (opC_v.HyT * opC_v.iMy * opC_v.By .+
                          opC_v.HxT * opC_v.iMx * opC_v.Bx)
    data_A[2,2] = pad(#=b0_v * opC_v.χ .+=#
                      (b0_v .+ b1_v) * (opC_v.HyT * opC_v.iMy * opC_v.Hy .+
                              opC_v.HxT * opC_v.iMx * opC_v.Hx) .-
                      opC_v.χ * a1_v)
    Av = [data_A[1,1] data_A[1,2];
          data_A[2,1] data_A[2,2]]

    data_B[1,1] = Mvm1 #.+ τ_2 .* Lvm1

    data_B[1,2] = spdiagm(0 => zeros(grid_v.nx*grid_v.ny))
    data_B[2,1] = spdiagm(0 => zeros(grid_v.nx*grid_v.ny))
    data_B[2,2] = spdiagm(0 => zeros(grid_v.nx*grid_v.ny))
    Bv = [data_B[1,1] data_B[1,2];
          data_B[2,1] data_B[2,2]]
    
    data_A[1,1] = pad(Lp, -4.0)
    data_A[1,2] = bc_Lp
    # data_A[6,5] = b0_p * opC_p.χ * Lp .+
    #               b1_p * (opC_p.HxT * opC_p.iMx * opC_p.Bx .+
    #                       opC_p.HyT * opC_p.iMy * opC_p.By)
    # data_A[6,6] = pad(b0_p * opC_p.χ * bc_Lp .+
    #                   b1_p * (opC_p.HxT * opC_p.iMx * opC_p.Hx .+
    #                           opC_p.HyT * opC_p.iMy * opC_p.Hy) .-
    #                   opC_p.χ * a1_p)
    GxT = opC_u.Gx'
    GyT = opC_v.Gy'
    data_A[2,1] = b1_p * (opC_p.HxT * opC_p.iMx * opC_p.Bx .+
                          opC_p.HyT * opC_p.iMy * opC_p.By)
    data_A[2,2] = pad(b0_p * (GxT * opC_u.Gx .+ GyT * opC_v.Gy) .+
                      b1_p * (opC_p.HxT * opC_p.iMx * opC_p.Hx .+
                              opC_p.HyT * opC_p.iMy * opC_p.Hy) .-
                      opC_p.χ * a1_p)

    Ap = [data_A[1,1] data_A[1,2];
          data_A[2,1] data_A[2,2]]
    
    rhs_u = zeros(2*grid_u.nx*grid_u.ny)
    veci(rhs_u,grid_u,2) .= opC_u.χ * vec(a0_u)
    rhs_v = zeros(2*grid_v.nx*grid_v.ny)
    veci(rhs_v,grid_v,2) .= opC_v.χ * vec(a0_v)
    rhs_p = zeros(2*grid.nx*grid.ny)
    veci(rhs_p,grid,2) .= opC_p.χ * vec(a0_p)
    
    return Au, Bu, rhs_u, Av, Bv, rhs_v, Ap, rhs_p
end

function set_navier_stokes(num, grid, geo, grid_u, geo_u, grid_v, geo_v,
                           opC_p, opC_u, opC_v, BC_p, BC_u, BC_v,
                           Mum1, Mvm1, iRe,
                           periodic_x, periodic_y)
    Lp, bc_Lp, Lu, bc_Lu, Lv, bc_Lv, Lp_fs, bc_Lp_fs, Lu_fs, bc_Lu_fs, Lv_fs, bc_Lv_fs = set_laplacians!(grid, geo, grid_u, geo_u, grid_v, geo_v,
        opC_p, opC_u, opC_v,
        periodic_x, periodic_y, true)

    Au, Bu, rhs_u, Av, Bv, rhs_v, Aϕ, rhs_ϕ = set_crank_nicolson_block(neu, num,
        grid, opC_p, Lp, bc_Lp, Lp_fs, bc_Lp_fs, BC_p,
        grid_u, opC_u, iRe.*Lu, iRe.*bc_Lu, Mum1, BC_u,
        grid_v, opC_v, iRe.*Lv, iRe.*bc_Lv, Mvm1, BC_v)

    return Au, Bu, rhs_u, Av, Bv, rhs_v, Aϕ, rhs_ϕ
end

function projection_fs!(num, grid, geo, grid_u, geo_u, grid_v, geo_v, ph,
                        BC_u, BC_v, BC_p,
                        opC_p, opC_u, opC_v,
                        Lum1, bc_Lum1, Lvm1, bc_Lvm1, Mum1, Mvm1,
                        FULL, MIXED, periodic_x, periodic_y, current_i
    )
    @unpack Re, τ, σ, g, β = num
    @unpack p, pD, ϕ, ϕD, Gxm1, Gym1, u, v, ucorrD, vcorrD, uD, vD, ucorr, vcorr = ph

    iRe = 1.0 / Re
    iτ = 1.0 / τ

    Au, Bu, rhs_u, Av, Bv, rhs_v, Aϕ, rhs_ϕ = set_navier_stokes(num, grid, geo, grid_u, geo_u, grid_v, geo_v,
                                                                opC_p, opC_u, opC_v, BC_p, BC_u, BC_v,
                                                                Mum1, Mvm1, iRe,
                                                                periodic_x, periodic_y)

    a0_u = zeros(grid_u)
    _a1_u = zeros(grid_u)
    _b0_u = ones(grid_u)
    _b1_u = zeros(grid_u)
    a0_v = zeros(grid_v)
    _a1_v = zeros(grid_v)
    _b0_v = ones(grid_v)
    _b1_v = zeros(grid_v)
    a0_p = zeros(grid)
    _a1_p = zeros(grid)
    _b0_p = ones(grid)
    _b1_p = zeros(grid)
    set_borders!(grid_u, a0_u, _a1_u, _b0_u, _b1_u, BC_u)
    set_borders!(grid_v, a0_v, _a1_v, _b0_v, _b1_v, BC_v)
    set_borders!(grid, a0_p, _a1_p, _b0_p, _b1_p, BC_p)
    b0_u = Diagonal(vec(_b0_u))
    b0_v = Diagonal(vec(_b0_v))
    b0_p = Diagonal(vec(_b0_p))

    mul!(rhs_u, Bu, uD, 1.0, 1.0)
    mul!(rhs_v, Bv, vD, 1.0, 1.0)

    grav_x = g .* sin(β) .* opC_u.M * fones(grid_u)
    grav_y = g .* cos(β) .* opC_v.M * fones(grid_v)

    veci(rhs_u,grid_u,1) .+= τ .* grav_x
    veci(rhs_v,grid_v,1) .+= - τ .* grav_y

    kill_dead_cells!(veci(rhs_u,grid_u,1), grid_u, geo_u)
    kill_dead_cells!(veci(rhs_u,grid_u,2), grid_u, geo_u)
    kill_dead_cells!(veci(rhs_v,grid_v,1), grid_v, geo_v)
    kill_dead_cells!(veci(rhs_v,grid_v,2), grid_v, geo_v)

    blocks = DDM.decompose(Au, grid_u.domdec, grid_u.domdec)
     _, ch = bicgstabl!(ucorrD, Au, rhs_u, Pl=ras(blocks,grid_u.pou), log=false,verbose=false)
    #println(ch)

    blocks = DDM.decompose(Av, grid_v.domdec, grid_v.domdec)
     _, ch = bicgstabl!(vcorrD, Av, rhs_v, Pl=ras(blocks,grid_v.pou), log=false,verbose=false)
    #println(ch)

    Duv = opC_p.AxT * veci(ucorrD,grid_u,1) .+ opC_p.Gx * veci(ucorrD,grid_u,2) .+
          opC_p.AyT * veci(vcorrD,grid_v,1) .+ opC_p.Gy * veci(vcorrD,grid_v,2)
    veci(rhs_ϕ,grid,1) .= iτ .* Duv

    Smat = strain_rate(opC_u, opC_v)
    S = Smat[1,1] * veci(ucorrD,grid_u,1) .+ Smat[1,2] * veci(ucorrD,grid_u,2) .+
        Smat[2,1] * veci(vcorrD,grid_v,1) .+ Smat[2,2] * veci(vcorrD,grid_v,2)

    GxT = opC_u.Gx'
    GyT = opC_v.Gy'
    veci(rhs_ϕ,grid,2) .= b0_p * (iRe .* S .- σ .* (GxT * opC_u.Gx .+ GyT * opC_v.Gy) * vec(grid.κ))
    blocks = DDM.decompose(Aϕ, grid.domdec, grid.domdec)
     _, ch = bicgstabl!(ϕD, Aϕ, rhs_ϕ, Pl=ras(blocks,grid.pou), log=true)
    #println(ch)

    ucorr .= reshape(veci(ucorrD,grid_u,1), (grid_u.ny, grid_u.nx))
    vcorr .= reshape(veci(vcorrD,grid_v,1), (grid_v.ny, grid_v.nx))
    ϕ .= reshape(veci(ϕD,grid,1), (grid.ny, grid.nx))

    iMu = Diagonal(1 ./ (opC_u.M.diag .+ eps(0.01)))
    iMv = Diagonal(1 ./ (opC_v.M.diag .+ eps(0.01)))

    ∇ϕ_x = opC_u.AxT * opC_u.Rx * veci(ϕD,grid,1) .+ opC_u.Gx * veci(ϕD,grid,2)
    ∇ϕ_y = opC_v.AyT * opC_v.Ry * veci(ϕD,grid,1) .+ opC_v.Gy * veci(ϕD,grid,2)

    iM = Diagonal(1. ./ (vec(geo.dcap[:,:,5]) .+ eps(0.01)))
    p .= ϕ
    Gxm1 .= ∇ϕ_x
    Gym1 .= ∇ϕ_y

    kill_dead_cells!(Gxm1, grid_u, geo_u)
    kill_dead_cells!(Gym1, grid_v, geo_v)

    u .= ucorr .- τ .* reshape(iMu * ∇ϕ_x, (grid_u.ny, grid_u.nx))
    v .= vcorr .- τ .* reshape(iMv * ∇ϕ_y, (grid_v.ny, grid_v.nx))

    veci(uD,grid_u,1) .= vec(u)
    veci(uD,grid_u,2) .= veci(ucorrD,grid_u,2)
    veci(vD,grid_v,1) .= vec(v)
    veci(vD,grid_v,2) .= veci(vcorrD,grid_v,2)

    return opC_p.M, opC_u.M, opC_v.M
end
