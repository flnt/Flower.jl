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

    bv[2] .= vec(tmp)

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
    A = blockarray(data_A)

    data_B = Matrix{SparseMatrixCSC{Float64, Int64}}(undef, 2, 2)
    data_B[1,1] = Mm1 .+ τ_2 .* Lm1
    data_B[1,2] = τ_2 .* bc_Lm1
    data_B[2,1] = spdiagm(0 => zeros(nx*ny))
    data_B[2,2] = spdiagm(0 => zeros(nx*ny))
    B = blockarray(data_B)

    data_rhs = Vector{Vector{Float64}}(undef, 2)
    data_rhs[1] = zeros(nx*ny)
    data_rhs[2] = χ * vec(a0)
    rhs = blockarray(data_rhs)
    
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
    A = blockarray(data_A)

    data_rhs = Vector{Vector{Float64}}(undef, 2)
    data_rhs[1] = zeros(nx*ny)
    data_rhs[2] = χ * vec(a0)
    rhs = blockarray(data_rhs)
    
    return A, rhs
end

function projection_no_slip!(num, grid, geo, grid_u, geo_u, grid_v, geo_v, ph,
                            BC_u, BC_v, BC_p,
                            opC_p, opC_u, opC_v,
                            ws_p, history_p, ws_u, history_u, ws_v, history_v,
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

    uD.data[1] .= vec(u)
    update_dirichlet_field!(grid_u, uD.data, u, BC_u)
    mul!(rhs_u, B_u, uD, 1.0, 1.0)
    rhs_u.data[1] .+= -τ .* (opC_p.Bx * pD.data[1] .+ opC_p.Hx * pD.data[2])
    @time solved, tired, broken, it = YAK.bicgstab!(uD, A_u, rhs_u, ws_u; Pl=I, Pr=Diagonal(A_u), rtol=1e-12, atol=1e-9, history = history_u)
    ucorr .= reshape(uD.data[1], (grid_u.ny, grid_u.nx))

    println("solved: $solved | tired: $tired | broken: $broken")
    monitor("None", history_u, it)

    vD.data[1] .= vec(v)
    update_dirichlet_field!(grid_v, vD.data, v, BC_v)
    mul!(rhs_v, B_v, vD, 1.0, 1.0)
    rhs_v.data[1] .+= -τ .* (opC_p.By * pD.data[1] .+ opC_p.Hy * pD.data[2])
    @time solved, tired, broken, it = YAK.bicgstab!(vD, A_v, rhs_v, ws_v; Pl=I, Pr=Diagonal(A_v), rtol=1e-12, atol=1e-9, history = history_v)
    vcorr .= reshape(vD.data[1], (grid_v.ny, grid_v.nx))

    println("solved: $solved | tired: $tired | broken: $broken")
    monitor("None", history_v, it)

    Duv =  opC_p.AxT * uD.data[1] .+ opC_p.HxT * uD.data[2] .+ opC_p.AyT * vD.data[1] .+ opC_p.HyT * vD.data[2]
    rhs_p.data[1] .=  iτ .* Duv

    sum_rhs = sum(rhs_p.data[1])
    sum_Mp = sum(geo.dcap[:,:,5])
    non_empty = vcat(FULL, MIXED)
    @inbounds @threads for II in non_empty
        pII = lexicographic(II, grid.ny)
        @inbounds rhs_p.data[1][pII] -= sum_rhs * geo.dcap[II,5] / sum_Mp
    end

    data = Matrix{Union{Diagonal{Float64,Vector{Float64}},SuiteSparse.UMFPACK.UmfpackLU{Float64, Int64}}}(undef, 2, 2)
    data[1,1] = lu(pad(Lp, -4.0))
    data[2,1] = Diagonal(zeros(grid.nx * grid.ny))
    data[2,2] = Diagonal(ones(grid.nx * grid.ny))
    LTdata = LowerTriangular(data)
    Pr = blockarray(LTdata)

    ϕD.data[1] .= 0.
    ϕD.data[2] .= 0.
    @time solved, tired, broken, it = YAK.bicgstab!(ϕD, A_p, rhs_p, ws_p; Pl=I, Pr=Pr, rtol=1e-12, atol=1e-9, history = history_p)
    # @time solved, tired, broken, it = YAK.bicgstab!(ϕD, A_p, rhs_p, ws_p; Pl=I, Pr=Diagonal(A_p), rtol=1e-12, atol=1e-9, history = history_p)
    # @time solved, tired, broken, it = YAK.bicgstab!(ϕD, A_p, rhs_p, ws_p; Pl=I, Pr=I, rtol=1e-12, atol=1e-8, history = history_p)
    ϕ .= reshape(ϕD.data[1], (grid.ny, grid.nx))

    println("solved: $solved | tired: $tired | broken: $broken")
    monitor("None", history_p, it)

    iM = Diagonal(1. ./ (vec(geo.dcap[:,:,5]) .+ eps(0.01)))
    pD.data[1] .= vec(p)
    pD.data[1] .+= ϕD.data[1] .- 0.5 .* iRe .* iM * Duv
    pD.data[2] .+= ϕD.data[2]
    p .= reshape(pD.data[1], (grid.ny, grid.nx))

    u .= ucorr .- τ .* reshape(opC_p.iMx * (opC_p.Bx * vec(ϕ) .+ opC_p.Hx * ϕD.data[2]), (grid_u.ny, grid_u.nx))
    v .= vcorr .- τ .* reshape(opC_p.iMy * (opC_p.By * vec(ϕ) .+ opC_p.Hy * ϕD.data[2]), (grid_v.ny, grid_v.nx))

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
    # _a1_p = ones(grid.ny, grid.nx) .* __a1
    # _b0_p = ones(grid.ny, grid.nx) .* __b0
    # _b1_p = ones(grid.ny, grid.nx) .* __b1
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
    data_A[1,3] = spdiagm(grid_u.nx*grid_u.ny, grid_v.nx*grid_v.ny, 0 => zeros(size_zeros))
    data_A[1,4] = spdiagm(grid_u.nx*grid_u.ny, grid_v.nx*grid_v.ny, 0 => zeros(size_zeros))
    data_A[1,5] = spdiagm(grid_u.nx*grid_u.ny, grid.nx*grid.ny, 0 => zeros(grid.nx*grid.ny))
    data_A[1,6] = spdiagm(grid_u.nx*grid_u.ny, grid.nx*grid.ny, 0 => zeros(grid.nx*grid.ny))

    data_A[2,1] = #=-b0_u * opC_u.χ .+=#
                  (b0_u .+ b1_u) * (opC_u.HxT * opC_u.iMx * opC_u.Bx .+
                          opC_u.HyT * opC_u.iMy * opC_u.By)
    data_A[2,2] = pad(#=b0_u * opC_u.χ .+=#
                      (b0_u .+ b1_u) * (opC_u.HxT * opC_u.iMx * opC_u.Hx .+
                              opC_u.HyT * opC_u.iMy * opC_u.Hy) .-
                      opC_u.χ * a1_u)
    data_A[2,3] = spdiagm(grid_u.nx*grid_u.ny, grid_v.nx*grid_v.ny, 0 => zeros(size_zeros))
    data_A[2,4] = spdiagm(grid_u.nx*grid_u.ny, grid_v.nx*grid_v.ny, 0 => zeros(size_zeros))
    data_A[2,5] = spdiagm(grid_u.nx*grid_u.ny, grid.nx*grid.ny, 0 => zeros(grid.nx*grid.ny))
    data_A[2,6] = spdiagm(grid_u.nx*grid_u.ny, grid.nx*grid.ny, 0 => zeros(grid.nx*grid.ny))

    data_A[3,1] = spdiagm(grid_v.nx*grid_v.ny, grid_u.nx*grid_u.ny, 0 => zeros(size_zeros))
    data_A[3,2] = spdiagm(grid_v.nx*grid_v.ny, grid_u.nx*grid_u.ny, 0 => zeros(size_zeros))
    data_A[3,3] = pad_crank_nicolson(opC_v.M .- τ .* Lv, grid_v, τ)
    data_A[3,4] = - τ .* bc_Lv
    data_A[3,5] = spdiagm(grid_v.nx*grid_v.ny, grid.nx*grid.ny, 0 => zeros(grid.nx*grid.ny))
    data_A[3,6] = spdiagm(grid_v.nx*grid_v.ny, grid.nx*grid.ny, 0 => zeros(grid.nx*grid.ny))

    data_A[4,1] = spdiagm(grid_v.nx*grid_v.ny, grid_u.nx*grid_u.ny, 0 => zeros(size_zeros))
    data_A[4,2] = spdiagm(grid_v.nx*grid_v.ny, grid_u.nx*grid_u.ny, 0 => zeros(size_zeros))
    data_A[4,3] = #=-b0_v * opC_v.χ .+=#
                  (b0_v .+ b1_v) * (opC_v.HyT * opC_v.iMy * opC_v.By .+
                          opC_v.HxT * opC_v.iMx * opC_v.Bx)
    data_A[4,4] = pad(#=b0_v * opC_v.χ .+=#
                      (b0_v .+ b1_v) * (opC_v.HyT * opC_v.iMy * opC_v.Hy .+
                              opC_v.HxT * opC_v.iMx * opC_v.Hx) .-
                      opC_v.χ * a1_v)
    data_A[4,5] = spdiagm(grid_v.nx*grid_v.ny, grid.nx*grid.ny, 0 => zeros(grid.nx*grid.ny))
    data_A[4,6] = spdiagm(grid_v.nx*grid_v.ny, grid.nx*grid.ny, 0 => zeros(grid.nx*grid.ny))

    data_A[5,1] = spdiagm(grid.nx*grid.ny, grid_u.nx*grid_u.ny, 0 => zeros(grid.nx*grid.ny))
    data_A[5,2] = spdiagm(grid.nx*grid.ny, grid_u.nx*grid_u.ny, 0 => zeros(grid.nx*grid.ny))
    data_A[5,3] = spdiagm(grid.nx*grid.ny, grid_v.nx*grid_v.ny, 0 => zeros(grid.nx*grid.ny))
    data_A[5,4] = spdiagm(grid.nx*grid.ny, grid_v.nx*grid_v.ny, 0 => zeros(grid.nx*grid.ny))
    data_A[5,5] = pad(Lp, -4.0)
    data_A[5,6] = bc_Lp

    data_A[6,1] = spdiagm(grid.nx*grid.ny, grid_u.nx*grid_u.ny, 0 => zeros(grid.nx*grid.ny))
    data_A[6,2] = spdiagm(grid.nx*grid.ny, grid_u.nx*grid_u.ny, 0 => zeros(grid.nx*grid.ny))
    data_A[6,3] = spdiagm(grid.nx*grid.ny, grid_v.nx*grid_v.ny, 0 => zeros(grid.nx*grid.ny))
    data_A[6,4] = spdiagm(grid.nx*grid.ny, grid_v.nx*grid_v.ny, 0 => zeros(grid.nx*grid.ny))
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

    A = blockarray(data_A)

    data_rhs = Vector{Vector{Float64}}(undef, 6)
    data_rhs[1] = zeros(grid_u.nx*grid_u.ny)
    data_rhs[2] = opC_u.χ * vec(a0_u)
    data_rhs[3] = zeros(grid_v.nx*grid_v.ny)
    data_rhs[4] = opC_v.χ * vec(a0_v)
    data_rhs[5] = zeros(grid.nx*grid.ny)
    data_rhs[6] = opC_p.χ * vec(a0_p)

    rhs = blockarray(data_rhs)
    
    return A, rhs
end

function projection_fs!(num, grid, geo, grid_u, geo_u, grid_v, geo_v, ph,
                        BC_u, BC_v, BC_p,
                        opC_p, opC_u, opC_v,
                        ws_uvϕ, history_uvϕ,
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

    # init_fresh_cells!(grid, p, geo.projection, FRESH, periodic_x, periodic_y)
    # init_fresh_cells!(grid, pD.data[1], geo.projection, FRESH, periodic_x, periodic_y)
    # init_fresh_cells!(grid, pD.data[2], pD.data[1], geo.projection, FRESH, periodic_x, periodic_y)
    # # Hardfix to ensure periodic borders are initialized
    # if periodic_y
    #     @inbounds @threads for II in FRESH
    #         pII = lexicographic(II, grid.ny)
    #         if II[1] == 1
    #             p[II] = p[δy⁺(II)]
    #             pD.data[1][pII] = pD.data[1][pII+1]
    #             pD.data[2][pII] = pD.data[2][pII+1]
    #         elseif II[1] == grid.ny
    #             p[II] = p[δy⁻(II)]
    #             pD.data[1][pII] = pD.data[1][pII-1]
    #             pD.data[2][pII] = pD.data[2][pII-1]
    #         end
    #     end
    # end
    # if periodic_x
    #     @inbounds @threads for II in FRESH
    #         pII = lexicographic(II, grid.ny)
    #         if II[2] == 1
    #             p[II] = p[δx⁺(II)]
    #             pD.data[1][pII] = pD.data[1][pII+grid.ny]
    #             pD.data[2][pII] = pD.data[2][pII+grid.ny]
    #         elseif II[2] == grid.nx
    #             p[II] = p[δx⁻(II)]
    #             pD.data[1][pII] = pD.data[1][pII-grid.ny]
    #             pD.data[2][pII] = pD.data[2][pII-grid.ny]
    #         end
    #     end
    # end

    # kill_dead_cells!(p, grid, geo)
    # kill_dead_cells!(pD.data[1], grid, geo)
    # kill_dead_cells!(pD.data[2], grid, geo)

    # init_fresh_cells!(grid_u, u, geo_u.projection, FRESH_u, periodic_x, periodic_y)
    # init_fresh_cells!(grid_v, v, geo_v.projection, FRESH_v, periodic_x, periodic_y)
    # kill_dead_cells!(u, grid_u, geo_u)
    # kill_dead_cells!(v, grid_v, geo_v)

    # rhs_uvϕ.data[1] .+= opC_u.M * vec(u)
    # rhs_uvϕ.data[3] .+= opC_v.M * vec(v)

    rhs_uvϕ.data[1] .+= Mum1 * vec(u)
    rhs_uvϕ.data[3] .+= Mvm1 * vec(v)

    # rhs_uvϕ.data[1] .+= -τ .* Gxm1
    # rhs_uvϕ.data[1] .+= τ .* σ .* b0_u * opC_u_m1.Gx * vec(grid.κ)
    # rhs_uvϕ.data[3] .+= -τ .* Gym1
    # rhs_uvϕ.data[3] .+= τ .* σ .* b0_v * opC_v_m1.Gy * vec(grid.κ)

    grav_x = g .* sin(β) .* Mum1 * ones(grid_u.nx * grid_u.ny)
    grav_y = g .* cos(β) .* Mvm1 * ones(grid_v.nx * grid_v.ny)

    rhs_uvϕ.data[1] .+= τ .* grav_x
    rhs_uvϕ.data[3] .+= - τ .* grav_y

    kill_dead_cells!(rhs_uvϕ.data[1], grid_u, geo_u)
    kill_dead_cells!(rhs_uvϕ.data[2], grid_u, geo_u)
    kill_dead_cells!(rhs_uvϕ.data[3], grid_v, geo_v)
    kill_dead_cells!(rhs_uvϕ.data[4], grid_v, geo_v)

    # data = Matrix{Union{Diagonal{Float64,Vector{Float64}},SuiteSparse.UMFPACK.UmfpackLU{Float64, Int64}}}(undef, 6, 6)
    # for i = 1:6, j = 1:6
    #     zu = zeros(grid_u.nx * grid_u.ny)
    #     if i != j
    #         if i < 3
    #             data[i,j] = Diagonal()
    #         elseif 
    #     end
    # end
    # data[1,1] = lu(pad(Lp, -4.0))
    # data[2,1] = Diagonal(zeros(grid.nx * grid.ny))
    # data[2,2] = Diagonal(ones(grid.nx * grid.ny))
    # LTdata = LowerTriangular(data_p)
    # Pr = blockarray(LTdata)

    # @time solved, tired, broken, it = YAK.bicgstab!(uvϕD, A_uvϕ, rhs_uvϕ, ws_uvϕ; Pl=I, Pr=Diagonal(A_uvϕ), rtol=1e-12, atol=1e-9, history = history_uvϕ, itmax=2000)
    # @time solved, tired, broken, it = YAK.bicgstab!(uvϕD, A_uvϕ, rhs_uvϕ, ws_uvϕ; Pl=I, Pr=Pr, rtol=1e-12, atol=1e-9, history = history_uvϕ, itmax=2000)

    # println("solved: $solved | tired: $tired | broken: $broken")
    # monitor("Diagonal", history_uvϕ, it)

    # Au = blockarray(A_uvϕ.data[1:2,1:2])
    # rhs_u = blockarray(rhs_uvϕ.data[1:2])
    # ws_u = bicgstabws(uD, Au, rhs_u)
    # history_u = eltype(ws_u)[]
    # @time solved, tired, broken, it = YAK.bicgstab!(uD, Au, rhs_u, ws_u; Pl=I, Pr=Diagonal(Au), rtol=1e-12, atol=1e-9, history = history_u, itmax=2000)
    # println("solved: $solved | tired: $tired | broken: $broken")
    # monitor("Diagonal", history_u, it)

    # Av = blockarray(A_uvϕ.data[3:4,3:4])
    # rhs_v = blockarray(rhs_uvϕ.data[3:4])
    # ws_v = bicgstabws(vD, Av, rhs_v)
    # history_v = eltype(ws_v)[]
    # @time solved, tired, broken, it = YAK.bicgstab!(vD, Av, rhs_v, ws_v; Pl=I, Pr=Diagonal(Av), rtol=1e-12, atol=1e-9, history = history_v, itmax=2000)
    # println("solved: $solved | tired: $tired | broken: $broken")
    # monitor("Diagonal", history_v, it)

    Auv = blockarray(A_uvϕ.data[1:4,1:4])
    rhs_uv = blockarray(rhs_uvϕ.data[1:4])
    ws_uv = bicgstabws(uvD, Auv, rhs_uv)
    history_uv = eltype(ws_uv)[]
    @time solved, tired, broken, it = YAK.bicgstab!(uvD, Auv, rhs_uv, ws_uv; Pl=I, Pr=Diagonal(Auv), rtol=1e-12, atol=1e-9, history = history_uv, itmax=2000)
    println("solved: $solved | tired: $tired | broken: $broken")
    monitor("Diagonal", history_uv, it)

    Aϕ = blockarray(A_uvϕ.data[5:6,5:6])
    rhs_ϕ = blockarray(rhs_uvϕ.data[5:6])
    ws_ϕ = bicgstabws(ϕD, Aϕ, rhs_ϕ)
    history_ϕ = eltype(ws_ϕ)[]

    Duv = opC_p.AxT * uvD.data[1] .+ opC_p.Gx * uvD.data[2] .+ opC_p.AyT * uvD.data[3] .+ opC_p.Gy * uvD.data[4]
    rhs_ϕ.data[1] .= iτ .* Duv

    GxT = opC_u.Gx'
    GyT = opC_v.Gy'
    S = GxT * ((2 .* opC_u.HxT * opC_u.iMx * opC_u.Bx .+ opC_u.HyT * opC_u.iMy * opC_u.By) * uD.data[1] .+
               (2 .* opC_u.HxT * opC_u.iMx * opC_u.Hx .+ opC_u.HyT * opC_u.iMy * opC_u.Hy) * uD.data[2]) .+
        GyT * ((2 .* opC_v.HyT * opC_v.iMy * opC_v.By .+ opC_v.HxT * opC_v.iMx * opC_v.Bx) * vD.data[1] .+
               (2 .* opC_v.HyT * opC_v.iMy * opC_v.Hy .+ opC_v.HxT * opC_v.iMx * opC_v.Hx) * vD.data[2])

    rhs_ϕ.data[2] .= b0_p * (iRe .* S .- σ .* (GxT * opC_u.Gx .+ GyT * opC_v.Gy) * vec(grid.κ))
    
    data = Matrix{Union{Diagonal{Float64,Vector{Float64}},SuiteSparse.UMFPACK.UmfpackLU{Float64, Int64}}}(undef, 2, 2)
    data[1,1] = lu(pad(Lp, -4.0))
    data[2,1] = Diagonal(zeros(grid.nx * grid.ny))
    data[2,2] = Diagonal(ones(grid.nx * grid.ny))
    LTdata = LowerTriangular(data)
    Pr = blockarray(LTdata)

    # ϕD.data[1] .= 0.
    # ϕD.data[2] .= 0.
    
    kill_dead_cells!(rhs_ϕ.data[1], grid, geo)
    kill_dead_cells!(rhs_ϕ.data[2], grid, geo)

    # @time solved, tired, broken, it = YAK.bicgstab!(ϕD, Aϕ, rhs_ϕ, ws_ϕ; Pl=I, Pr=Diagonal(Aϕ), rtol=1e-12, atol=1e-9, history = history_ϕ, itmax=2000)
    @time solved, tired, broken, it = YAK.bicgstab!(ϕD, Aϕ, rhs_ϕ, ws_ϕ; Pl=I, Pr=Pr, rtol=1e-12, atol=1e-9, history = history_ϕ, itmax=2000)
    println("solved: $solved | tired: $tired | broken: $broken")
    monitor("Diagonal", history_ϕ, it)

    # _A = [Aϕ.data[1,1] Aϕ.data[1,2]; Aϕ.data[2,1] Aϕ.data[2,2]]
    # _ϕD = _A \ vcat(rhs_ϕ.data[1], rhs_ϕ.data[1])
    # ϕD.data[1] .= _ϕD[1:grid.ny*grid.nx]
    # ϕD.data[2] .= _ϕD[grid.ny*grid.nx+1:end]

    ucorr .= reshape(uvD.data[1], (grid_u.ny, grid_u.nx))
    vcorr .= reshape(uvD.data[3], (grid_v.ny, grid_v.nx))
    ϕ .= reshape(ϕD.data[1], (grid.ny, grid.nx))

    # kill_dead_cells!(p, grid, geo)
    # kill_dead_cells!(pD.data[1], grid, geo)
    # kill_dead_cells!(pD.data[2], grid, geo)

    # pD.data[1] .= vec(p)
    # pD.data[1] .+= uvϕD.data[5]
    # pD.data[2] .+= uvϕD.data[6]
    # p .= reshape(pD.data[1], (grid.ny, grid.nx))

    iMu = Diagonal(1 ./ (opC_u.M.diag .+ eps(0.01)))
    iMv = Diagonal(1 ./ (opC_v.M.diag .+ eps(0.01)))

    ∇ϕ_x = opC_u.AxT * opC_u.Rx * ϕD.data[1] .+ opC_u.Gx * ϕD.data[2]
    ∇ϕ_y = opC_v.AyT * opC_v.Ry * ϕD.data[1] .+ opC_v.Gy * ϕD.data[2]

    # Gxm1 .= Gxm1 + ∇ϕ_x
    # Gym1 .= Gym1 + ∇ϕ_y

    p .= ϕ
    Gxm1 .= ∇ϕ_x
    Gym1 .= ∇ϕ_y

    kill_dead_cells!(Gxm1, grid_u, geo_u)
    kill_dead_cells!(Gym1, grid_v, geo_v)

    u .= ucorr .- τ .* reshape(iMu * ∇ϕ_x, (grid_u.ny, grid_u.nx))
    v .= vcorr .- τ .* reshape(iMv * ∇ϕ_y, (grid_v.ny, grid_v.nx))

    return Aϕ, Lu, bc_Lu, Lv, bc_Lv, opC_p.M, opC_u.M, opC_v.M, reshape(Duv, (grid.ny, grid.nx))
end