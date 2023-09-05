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

function set_borders!(grid, geo, a0, a1, b, BC, per_x, per_y)
    @unpack nx, ny, x, y, u, ind = grid
    
    if per_y
        ind_x = 2:nx-1
        ind_y = 1:ny
    else
        ind_x = 1:nx
        ind_y = 2:ny-1
    end

    χx = abs.(geo.dcap[:,:,3] .- geo.dcap[:,:,1])
    χy = abs.(geo.dcap[:,:,4] .- geo.dcap[:,:,2])

    idx = CartesianIndices((ind_y,1))
    @inbounds a0[idx] .= χx[idx] .* (BC.left.val .* ones(grid))[idx]
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
        idx_cl = intersect(idx, ind.cl)
        ϵb = 5.0 .* sqrt.((x[δx⁺.(idx_cl)] .- x[idx_cl]).^2 .+ (y[δx⁺.(idx_cl)] .- y[idx_cl]).^2)

        @inbounds a1[idx_cl] .= -1.0
        @inbounds a1[symdiff(idx, ind.cl)] .= -1.0
        @inbounds b[idx_cl] .= BC.left.λ .* bell_function2.(u[idx_cl], ϵb)
        @inbounds b[symdiff(idx, ind.cl)] .= 0.0
    elseif is_gnbc(BC.left)
        bell = bell_function(grid, BC.ϵ)
        @inbounds a0[idx] .+= bell[idx] .* BC.σ ./ BC.μ .* (cos(θd) .- cos(BC.θe))
        @inbounds a1[intersect(idx, ind.cl)] .= -1.0
        @inbounds a1[symdiff(idx, ind.cl)] .= -1.0
        @inbounds b[intersect(idx, ind.cl)] .= BC.left.λ
        @inbounds b[symdiff(idx, ind.cl)] .= 0.0
    else
        @error ("Not implemented yet")
    end

    idx = CartesianIndices((1,ind_x))
    @inbounds a0[idx] .= χy[idx] .* (BC.bottom.val .* ones(grid))[idx]
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
        idx_cl = intersect(idx, ind.cl)
        ϵb = 5.0 .* sqrt.((x[δy⁺.(idx_cl)] .- x[idx_cl]).^2 .+ (y[δy⁺.(idx_cl)] .- y[idx_cl]).^2)

        @inbounds a1[idx_cl] .= -1.0
        @inbounds a1[symdiff(idx, ind.cl)] .= -1.0
        @inbounds b[idx_cl] .= BC.bottom.λ .* bell_function2.(u[idx_cl], ϵb)
        @inbounds b[symdiff(idx, ind.cl)] .= 0.0
    elseif is_gnbc(BC.bottom)
        @inbounds a0[idx] .+= bell[idx] .* BC.σ ./ BC.μ .* (cos(θd) .- cos(BC.θe))
        @inbounds a1[intersect(idx, ind.cl)] .= -1.0
        @inbounds a1[symdiff(idx, ind.cl)] .= -1.0
        @inbounds b[intersect(idx, ind.cl)] .= BC.bottom.λ
        @inbounds b[symdiff(idx, ind.cl)] .= 0.0
    else
        @error ("Not implemented yet")
    end

    idx = CartesianIndices((ind_y,nx:nx))
    @inbounds a0[idx] .= χx[idx] .* (BC.right.val .* ones(grid))[idx]
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
        idx_cl = intersect(idx, ind.cl)
        ϵb = 5.0 .* sqrt.((x[δx⁻.(idx_cl)] .- x[idx_cl]).^2 .+ (y[δx⁻.(idx_cl)] .- y[idx_cl]).^2)

        @inbounds a1[idx_cl] .= -1.0
        @inbounds a1[symdiff(idx, ind.cl)] .= -1.0
        @inbounds b[idx_cl] .= BC.right.λ .* bell_function2.(u[idx_cl], ϵb)
        @inbounds b[symdiff(idx, ind.cl)] .= 0.0
    elseif is_gnbc(BC.right)
        @inbounds a0[idx] .+= bell[idx] .* BC.σ ./ BC.μ .* (cos(θd) .- cos(BC.θe))
        @inbounds a1[intersect(idx, ind.cl)] .= -1.0
        @inbounds a1[symdiff(idx, ind.cl)] .= -1.0
        @inbounds b[intersect(idx, ind.cl)] .= BC.right.λ
        @inbounds b[symdiff(idx, ind.cl)] .= 0.0
    else
        @error ("Not implemented yet")
    end

    idx = CartesianIndices((ny:ny,ind_x))
    @inbounds a0[idx] .= χy[idx] .* (BC.top.val .* ones(grid))[idx]
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
        idx_cl = intersect(idx, ind.cl)
        ϵb = 5.0 .* sqrt.((x[δy⁻.(idx_cl)] .- x[idx_cl]).^2 .+ (y[δy⁻.(idx_cl)] .- y[idx_cl]).^2)

        @inbounds a1[idx_cl] .= -1.0
        @inbounds a1[symdiff(idx, ind.cl)] .= -1.0
        @inbounds b[idx_cl] .= BC.top.λ .* bell_function2.(u[idx_cl], ϵb)
        @inbounds b[symdiff(idx, ind.cl)] .= 0.0
    elseif is_gnbc(BC.top)
        @inbounds a0[idx] .+= bell[idx] .* BC.σ ./ BC.μ .* (cos(θd) .- cos(BC.θe))
        @inbounds a1[intersect(idx, ind.cl)] .= -1.0
        @inbounds a1[symdiff(idx, ind.cl)] .= -1.0
        @inbounds b[intersect(idx, ind.cl)] .= BC.top.λ
        @inbounds b[symdiff(idx, ind.cl)] .= 0.0
    else
        @error ("Not implemented yet")
    end

    return nothing
end

function update_dirichlet_field!(grid, bv, v, BC, per_x, per_y)
    @unpack nx, ny = grid
    tmp = zeros(size(v))

    if per_y
        ind_x = 2:nx-1
        ind_y = 1:ny
    else
        ind_x = 1:nx
        ind_y = 2:ny-1
    end

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

    mul!(tmp_x, iMx, Hx)
    bc_L = BxT * tmp_x
    mul!(tmp_y, iMy, Hy)
    bc_L = bc_L .+ ByT * tmp_y

    return L, bc_L
end

function set_laplacians!(grid, geo, grid_u, geo_u, grid_v, geo_v,
                         opC_p, opC_u, opC_v,
                         periodic_x, periodic_y
    )
    @unpack ny, ind = grid

    set_cutcell_matrices!(grid, geo, opC_p, periodic_x, periodic_y)
    set_cutcell_matrices!(grid_u, geo_u, opC_u, periodic_x, periodic_y)
    set_cutcell_matrices!(grid_v, geo_v, opC_v, periodic_x, periodic_y)

    set_other_cutcell_matrices(grid, geo, geo_u, geo_v, opC_p, opC_u, opC_v, periodic_x, periodic_y)

    Lp, bc_Lp = laplacian(opC_p)
    Lu, bc_Lu = laplacian(opC_u)
    Lv, bc_Lv = laplacian(opC_v)

    return Lp, bc_Lp, Lu, bc_Lu, Lv, bc_Lv
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

function set_convection!(grid, geo, grid_u, geo_u, grid_v, geo_v,
                         op, ph, BC_u, BC_v
    )
    @unpack Cu, CUTCu, Cv, CUTCv = op
    @unpack u, v, uD, vD = ph

    Du = reshape(veci(uD,grid_u,2), (grid_u.ny, grid_u.nx))
    Dv = reshape(veci(vD,grid_v,2), (grid_v.ny, grid_v.nx))

    Hu = zeros(grid_u)
    for II in vcat(grid_u.ind.b_left[1], grid_u.ind.b_bottom[1], grid_u.ind.b_right[1], grid_u.ind.b_top[1])
        Hu[II] = distance(grid_u.mid_point[II], geo_u.centroid[II], grid_u.dx[II], grid_u.dy[II])
    end

    Hv = zeros(grid_v)
    for II in vcat(grid_v.ind.b_left[1], grid_v.ind.b_bottom[1], grid_v.ind.b_right[1], grid_v.ind.b_top[1])
        Hv[II] = distance(grid_v.mid_point[II], geo_v.centroid[II], grid_v.dx[II], grid_v.dy[II])
    end

    bcuCu_x, bcuCu_y, bcvCu_x, bcvCu_y = set_bc_bnds(dir, GridFCx, Du, Dv, Hu, Hv, u, v, BC_u, BC_v)
    bcvCv_x, bcvCv_y, bcuCv_x, bcuCv_y = set_bc_bnds(dir, GridFCy, Dv, Du, Hv, Hu, v, u, BC_v, BC_u)

    vector_convection!(dir, GridFCx, Cu, CUTCu, u, v, bcuCu_x, bcuCu_y, bcvCu_x, bcvCu_y,
            geo.dcap, grid.ny, BC_u, grid_u.ind.inside,
            grid_u.ind.b_left[1], grid_u.ind.b_bottom[1], grid_u.ind.b_right[1], grid_u.ind.b_top[1])
    vector_convection!(dir, GridFCy, Cv, CUTCv, u, v, bcuCv_x, bcuCv_y, bcvCv_x, bcvCv_y,
            geo.dcap, grid.ny, BC_v, grid_v.ind.inside,
            grid_v.ind.b_left[1], grid_v.ind.b_bottom[1], grid_v.ind.b_right[1], grid_v.ind.b_top[1])
    
    return nothing
end

function CN_set_momentum(
    bc_type, num, grid, geo, opC,
    L, bc_L, Lm1, bc_Lm1, Mm1, BC,
    per_x, per_y
    )
    @unpack τ = num
    @unpack Bx, By, Hx, Hy, HxT, HyT, χ, M, iMx, iMy = opC

    if bc_type == dir
        vel = copy(grid.V)
        __a1 = -1.0
        __b = 0.0
    elseif bc_type == neu
        vel = 0.0
        __a1 = 0.0
        __b = 0.0
    elseif bc_type == rob
        vel = 0.0
        __a1 = -1.0
        __b = 1.0
    end

    a0 = ones(grid) .* vel
    _a1 = ones(grid) .* __a1
    _b = ones(grid) .* __b
    set_borders!(grid, geo, a0, _a1, _b, BC, per_x, per_y)
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
    data_B[2,1] = spdiagm(0 => fzeros(grid))
    data_B[2,2] = spdiagm(0 => fzeros(grid))
    B  = [data_B[1,1] data_B[1,2];
          data_B[2,1] data_B[2,2]]

    rhs = f2zeros(grid)
    veci(rhs,grid,2) .= vec(a0)
    
    return A, B, rhs
end

function FE_set_momentum(
    bc_type, num, grid, geo, opC,
    L, bc_L, Mm1, BC,
    per_x, per_y
    )
    @unpack τ = num
    @unpack Bx, By, Hx, Hy, HxT, HyT, χ, M, iMx, iMy = opC

    if bc_type == dir
        vel = copy(grid.V)
        __a1 = -1.0
        __b = 0.0
    elseif bc_type == neu
        vel = 0.0
        __a1 = 0.0
        __b = 0.0
    elseif bc_type == rob
        vel = 0.0
        __a1 = -1.0
        __b = 1.0
    end

    a0 = ones(grid) .* vel
    _a1 = ones(grid) .* __a1
    _b = ones(grid) .* __b
    set_borders!(grid, geo, a0, _a1, _b, BC, per_x, per_y)
    a1 = Diagonal(vec(_a1))
    b = Diagonal(vec(_b))

    data_A = Matrix{SparseMatrixCSC{Float64, Int64}}(undef, 2, 2)
    data_A[1,1] = pad_crank_nicolson(M .- τ .* L, grid, τ)
    data_A[1,2] = - τ .* bc_L
    data_A[2,1] = b * (HxT * iMx * Bx .+ HyT * iMy * By)
    data_A[2,2] = pad(b * (HxT * iMx * Hx .+ HyT * iMy * Hy) .- χ * a1)
    A  = [data_A[1,1] data_A[1,2];
          data_A[2,1] data_A[2,2]]

    data_B = Matrix{SparseMatrixCSC{Float64, Int64}}(undef, 2, 2)
    data_B[1,1] = Mm1
    data_B[1,2] = spdiagm(0 => fzeros(grid))
    data_B[2,1] = spdiagm(0 => fzeros(grid))
    data_B[2,2] = spdiagm(0 => fzeros(grid))
    B  = [data_B[1,1] data_B[1,2];
          data_B[2,1] data_B[2,2]]

    rhs = f2zeros(grid)
    veci(rhs,grid,2) .= vec(a0)
    
    return A, B, rhs
end

function set_poisson(bc_type, grid, geo, a0, opC, L, bc_L, BC, per_x, per_y)
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

    _a1 = ones(grid) .* __a1
    _b = ones(grid) .* __b
    set_borders!(grid, geo, a0, _a1, _b, BC, per_x, per_y)
    a1 = Diagonal(vec(_a1))
    b = Diagonal(vec(_b))

    data_A = Matrix{SparseMatrixCSC{Float64, Int64}}(undef, 2, 2)
    data_A[1,1] = pad(L, -4.0)
    data_A[1,2] = bc_L
    data_A[2,1] = b * (HxT * iMx * Bx .+ HyT * iMy * By)
    data_A[2,2] = pad(b * (HxT * iMx * Hx .+ HyT * iMy * Hy) .- χ * a1)
    A  = [data_A[1,1] data_A[1,2];
          data_A[2,1] data_A[2,2]]

    rhs = f2zeros(grid)
    veci(rhs,grid,2) .= vec(a0)
    
    return A, rhs
end

function set_CN!(
    bc_type_u, bc_type_p, num, grid, geo, grid_u, geo_u, grid_v, geo_v,
    opC_p, opC_u, opC_v, BC_p, BC_u, BC_v,
    Lum1, bc_Lum1, Lvm1, bc_Lvm1, Mum1, Mvm1, iRe,
    op_conv, ph,
    periodic_x, periodic_y, advection
    )

    laps = set_laplacians!(
        grid, geo, grid_u, geo_u, grid_v, geo_v,
        opC_p, opC_u, opC_v,
        periodic_x, periodic_y
    )
    Lp, bc_Lp, Lu, bc_Lu, Lv, bc_Lv = laps
    
    if advection
        set_convection!(grid, geo, grid_u, geo_u, grid_v, geo_v, op_conv, ph, BC_u, BC_v)
    end

    Au, Bu, rhs_u = CN_set_momentum(
        bc_type_u, num, grid_u, geo_u, opC_u,
        iRe.*Lu, iRe.*bc_Lu, iRe.*Lum1, iRe.*bc_Lum1, Mum1, BC_u,
        periodic_x, periodic_y
    )
    Av, Bv, rhs_v = CN_set_momentum(
        bc_type_u, num, grid_v, geo_v, opC_v,
        iRe.*Lv, iRe.*bc_Lv, iRe.*Lvm1, iRe.*bc_Lvm1, Mvm1, BC_v,
        periodic_x, periodic_y
    )
    a0_p = zeros(grid)
    Aϕ, rhs_ϕ = set_poisson(bc_type_p, grid, geo, a0_p, opC_p, Lp, bc_Lp, BC_p, periodic_x, periodic_y)

    return Au, Bu, rhs_u, Av, Bv, rhs_v, Aϕ, rhs_ϕ, Lu, bc_Lu, Lv, bc_Lv
end

function set_FE!(bc_type_u, bc_type_p, num, grid, geo, grid_u, geo_u, grid_v, geo_v,
                           opC_p, opC_u, opC_v, BC_p, BC_u, BC_v,
                           Mum1, Mvm1, iRe,
                           op_conv, ph,
                           periodic_x, periodic_y, advection)

    laps = set_laplacians!(
        grid, geo, grid_u, geo_u, grid_v, geo_v,
        opC_p, opC_u, opC_v,
        periodic_x, periodic_y
    )
    Lp, bc_Lp, Lu, bc_Lu, Lv, bc_Lv = laps

    if advection
        set_convection!(grid, geo, grid_u, geo_u, grid_v, geo_v, op_conv, ph, BC_u, BC_v)
    end

    Au, Bu, rhs_u = FE_set_momentum(
        bc_type_u, num, grid_u, geo_u, opC_u,
        iRe.*Lu, iRe.*bc_Lu, Mum1, BC_u,
        periodic_x, periodic_y
    )
    Av, Bv, rhs_v = FE_set_momentum(
        bc_type_u, num, grid_v, geo_v, opC_v,
        iRe.*Lv, iRe.*bc_Lv, Mvm1, BC_v,
        periodic_x, periodic_y
    )
    a0_p = zeros(grid)
    Aϕ, rhs_ϕ = set_poisson(bc_type_p, grid, geo, a0_p, opC_p, Lp, bc_Lp, BC_p, periodic_x, periodic_y)

    return Au, Bu, rhs_u, Av, Bv, rhs_v, Aϕ, rhs_ϕ, Lu, bc_Lu, Lv, bc_Lv
end

function pressure_projection!(
    time_scheme, bc_type_u, bc_type_p,
    num, grid, geo, grid_u, geo_u, grid_v, geo_v, ph,
    BC_u, BC_v, BC_p,
    opC_p, opC_u, opC_v, op_conv,
    Lum1, bc_Lum1, Lvm1, bc_Lvm1,
    Cum1, Cvm1, Mum1, Mvm1,
    periodic_x, periodic_y, advection
    )
    @unpack Re, τ, σ, g, β = num
    @unpack p, pD, ϕ, ϕD, u, v, ucorrD, vcorrD, uD, vD, ucorr, vcorr = ph
    @unpack Cu, Cv, CUTCu, CUTCv = op_conv

    iRe = 1.0 / Re
    iτ = 1.0 / τ

    if is_FE(time_scheme)
        Au, Bu, rhs_u, Av, Bv, rhs_v, Aϕ, rhs_ϕ, Lu, bc_Lu, Lv, bc_Lv = set_FE!(
            bc_type_u, bc_type_p, num, grid, geo, grid_u, geo_u, grid_v, geo_v,
            opC_p, opC_u, opC_v, BC_p, BC_u, BC_v,
            Mum1, Mvm1, iRe,
            op_conv, ph,
            periodic_x, periodic_y, advection
        )
    elseif is_CN(time_scheme)
        Au, Bu, rhs_u, Av, Bv, rhs_v, Aϕ, rhs_ϕ, Lu, bc_Lu, Lv, bc_Lv = set_CN!(
            bc_type_u, bc_type_p, num, grid, geo, grid_u, geo_u, grid_v, geo_v,
            opC_p, opC_u, opC_v, BC_p, BC_u, BC_v,
            Lum1, bc_Lum1, Lvm1, bc_Lvm1, Mum1, Mvm1, iRe,
            op_conv, ph,
            periodic_x, periodic_y, advection
        )
    end

    a0_p = zeros(grid)
    _a1_p = zeros(grid)
    _b_p = ones(grid)
    set_borders!(grid, geo, a0_p, _a1_p, _b_p, BC_p, periodic_x, periodic_y)
    b_p = Diagonal(vec(_b_p))

    grav_x = g .* sin(β) .* opC_u.M * fones(grid_u)
    grav_y = g .* cos(β) .* opC_v.M * fones(grid_v)

    Convu = fzeros(grid_u)
    Convv = fzeros(grid_v)
    Cui = Cu * vec(u) .+ CUTCu
    Cvi = Cv * vec(v) .+ CUTCv
    if advection
        Convu .+= 1.5 .* Cui .- 0.5 .* Cum1
        Convv .+= 1.5 .* Cvi .- 0.5 .* Cvm1
    end

    if is_dirichlet(bc_type_u)
        veci(uD,grid_u,1) .= vec(u)
        update_dirichlet_field!(grid_u, uD, u, BC_u, periodic_x, periodic_y)
        veci(rhs_u,grid_u,1) .+= -τ .* (opC_u.AxT * opC_u.Rx * veci(pD,grid,1) .+ opC_u.Gx * veci(pD,grid,2))
    end
    mul!(rhs_u, Bu, uD, 1.0, 1.0)
    veci(rhs_u,grid_u,1) .+= τ .* grav_x
    veci(rhs_u,grid_u,1) .-= τ .* Convu
    kill_dead_cells!(veci(rhs_u,grid_u,1), grid_u, geo_u)
    kill_dead_cells!(veci(rhs_u,grid_u,2), grid_u, geo_u)
    blocks = DDM.decompose(Au, grid_u.domdec, grid_u.domdec)
    bicgstabl!(ucorrD, Au, rhs_u, Pl=ras(blocks,grid_u.pou), log=true)
    # @mytime _, ch = bicgstabl!(ucorrD, Au, rhs_u, Pl=ras(blocks,grid_u.pou), log=true)
    # println(ch)
    ucorr .= reshape(veci(ucorrD,grid_u,1), grid_u)

    if is_dirichlet(bc_type_u)
        veci(vD,grid_v,1) .= vec(v)
        update_dirichlet_field!(grid_v, vD, v, BC_v, periodic_x, periodic_y)
        veci(rhs_v,grid_v,1) .+= -τ .* (opC_v.AyT * opC_v.Ry * veci(pD,grid,1) .+ opC_v.Gy * veci(pD,grid,2))
    end
    mul!(rhs_v, Bv, vD, 1.0, 1.0)
    veci(rhs_v,grid_v,1) .+= - τ .* grav_y
    veci(rhs_v,grid_v,1) .-= τ .* Convv
    kill_dead_cells!(veci(rhs_v,grid_v,1), grid_v, geo_v)
    kill_dead_cells!(veci(rhs_v,grid_v,2), grid_v, geo_v)
    blocks = DDM.decompose(Av, grid_v.domdec, grid_v.domdec)
    bicgstabl!(vcorrD, Av, rhs_v, Pl=ras(blocks,grid_v.pou), log=true)
    # @mytime _, ch = bicgstabl!(vcorrD, Av, rhs_v, Pl=ras(blocks,grid_v.pou), log=true)
    # println(ch)
    vcorr .= reshape(veci(vcorrD,grid_v,1), grid_v)

    Duv = opC_p.AxT * veci(ucorrD,grid_u,1) .+ opC_p.Gx * veci(ucorrD,grid_u,2) .+
          opC_p.AyT * veci(vcorrD,grid_v,1) .+ opC_p.Gy * veci(vcorrD,grid_v,2)
    veci(rhs_ϕ,grid,1) .= iτ .* Duv

    if is_dirichlet(bc_type_p)
        Smat = strain_rate(opC_u, opC_v)
        S = Smat[1,1] * veci(ucorrD,grid_u,1) .+ Smat[1,2] * veci(ucorrD,grid_u,2) .+
            Smat[2,1] * veci(vcorrD,grid_v,1) .+ Smat[2,2] * veci(vcorrD,grid_v,2)

        GxT = opC_u.Gx'
        GyT = opC_v.Gy'
        veci(rhs_ϕ,grid,2) .= b_p * (iRe .* S .- σ .* (GxT * opC_u.Gx .+ GyT * opC_v.Gy) * vec(grid.κ))
    end
    blocks = DDM.decompose(Aϕ, grid.domdec, grid.domdec)
    bicgstabl!(ϕD, Aϕ, rhs_ϕ, Pl=ras(blocks,grid.pou), log=true)
    # @mytime _, ch = bicgstabl!(ϕD, Aϕ, rhs_ϕ, Pl=ras(blocks,grid.pou), log=true)
    # println(ch)
    ϕ .= reshape(veci(ϕD,grid,1), grid)

    iMu = Diagonal(1 ./ (opC_u.M.diag .+ eps(0.01)))
    iMv = Diagonal(1 ./ (opC_v.M.diag .+ eps(0.01)))

    ∇ϕ_x = opC_u.AxT * opC_u.Rx * vec(ϕ) .+ opC_u.Gx * veci(ϕD,grid,2)
    ∇ϕ_y = opC_v.AyT * opC_v.Ry * vec(ϕ) .+ opC_v.Gy * veci(ϕD,grid,2)

    iM = Diagonal(1. ./ (vec(geo.dcap[:,:,5]) .+ eps(0.01)))
    if is_dirichlet(bc_type_p)
        veci(pD,grid,1) .= vec(ϕ) #.- iRe .* reshape(iM * Duv, grid))
        veci(pD,grid,2) .+= veci(ϕD,grid,2)
        p .= reshape(veci(pD,grid,1), grid)
    else
        veci(pD,grid,1) .= vec(p) .+ vec(ϕ) #.- iRe .* iM * Duv
        veci(pD,grid,2) .+= veci(ϕD,grid,2)
        p .= reshape(veci(pD,grid,1), grid)
    end
    
    u .= ucorr .- τ .* reshape(iMu * ∇ϕ_x, grid_u)
    v .= vcorr .- τ .* reshape(iMv * ∇ϕ_y, grid_v)

    veci(uD,grid_u,1) .= vec(u)
    veci(uD,grid_u,2) .= veci(ucorrD,grid_u,2)
    veci(vD,grid_v,1) .= vec(v)
    veci(vD,grid_v,2) .= veci(vcorrD,grid_v,2)

    return Lu, bc_Lu, Lv, bc_Lv, opC_p.M, opC_u.M, opC_v.M, Cui, Cvi
end