@inline function apply_curvature(num, grid, bc, all_indices)
    @unpack ϵ_κ, ϵ_V = num
    @unpack κ, V = grid

    @inbounds @threads for II in all_indices
        @inbounds bc[II] = bc[II] - ϵ_κ*κ[II] - ϵ_V*V[II]
    end
    return nothing
end

@inline function apply_anisotropy(num, grid, bc, MIXED, sol_projection)
    @unpack ϵ_κ, ϵ_V, m, θ₀ = num
    @unpack κ, V = grid

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
    op_conv, grid_u, geo_u, grid_v, geo_v,
    periodic_x, periodic_y, convection)
    @unpack τ, aniso = num
    @unpack nx, ny, dx, dy, ind, mid_point = grid
    @unpack all_indices, inside, b_left, b_bottom, b_right, b_top = ind
    @unpack Bx, By, BxT, ByT, Hx, Hy, HxT, HyT, iMx, iMy, χ = op
    @unpack CT, CUTCT = op_conv
    @unpack u, v, uD, vD = ph

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

    # Flags with BCs
    a0 = ones(grid) .* θd
    if aniso
        apply_anisotropy(num, grid, a0, MIXED, projection)
    else
        apply_curvature(num, grid, a0, all_indices)
    end
    _a1 = ones(ny, nx) .* __a1
    _b = ones(ny, nx) .* __b
    set_borders!(grid, a0, _a1, _b, BC_T, periodic_x, periodic_y)
    a1 = Diagonal(vec(_a1))
    b = Diagonal(vec(_b))

    if convection
        HT = zeros(grid)
        @inbounds @threads for II in vcat(b_left[1], b_bottom[1], b_right[1], b_top[1])
            HT[II] = distance(mid_point[II], geo.centroid[II], dx[II], dy[II])
        end    
        bcTx, bcTy = set_bc_bnds(dir, a0, HT, BC_T)
    
        Hu = zeros(grid_u)
        for II in vcat(grid_u.ind.b_left[1], grid_u.ind.b_bottom[1], grid_u.ind.b_right[1], grid_u.ind.b_top[1])
            Hu[II] = distance(grid_u.mid_point[II], geo_u.centroid[II], grid_u.dx[II], grid_u.dy[II])
        end
    
        Hv = zeros(grid_v) 
        for II in vcat(grid_v.ind.b_left[1], grid_v.ind.b_bottom[1], grid_v.ind.b_right[1], grid_v.ind.b_top[1])
            Hv[II] = distance(grid_v.mid_point[II], geo_v.centroid[II], grid_v.dx[II], grid_v.dy[II])
        end
    
        bcU = reshape(veci(uD,grid_u,2), (grid_u.ny, grid_u.nx))
        bcV = reshape(veci(vD,grid_v,2), (grid_v.ny, grid_v.nx))
        scalar_convection!(dir, CT, CUTCT, u, v, bcTx, bcTy, bcU, bcV, geo.dcap, ny, BC_T, inside, b_left[1], b_bottom[1], b_right[1], b_top[1])
    end

    χx = (geo.dcap[:,:,3] .- geo.dcap[:,:,1]) .^ 2
    χy = (geo.dcap[:,:,4] .- geo.dcap[:,:,2]) .^ 2
    χ .= Diagonal(sqrt.(vec(χx .+ χy))) 

    # Mass matrices
    M = Diagonal(vec(geo.dcap[:,:,5]))
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

    # Matrices for BCs
    bc_matrix!(Hx, Hy, geo.dcap, ny, ind.all_indices)

    mat_assign_T!(HxT, sparse(Hx'))
    mat_assign_T!(HyT, sparse(Hy'))

    periodic_bcs!(grid, Bx, By, Hx, Hy, periodic_x, periodic_y)
    
    LT = BxT * iMx * Bx .+ ByT * iMy * By
    LD = BxT * iMx * Hx .+ ByT * iMy * Hy

    dataA = Matrix{SparseMatrixCSC{Float64, Int64}}(undef, 2, 2)
    dataA[1,1] = pad_crank_nicolson(M .- 0.5 .* τ .* LT, grid, τ)
    dataA[1,2] = - 0.5 .* τ .* LD
    dataA[2,1] = b * (HxT * iMx * Bx .+ HyT * iMy * By)
    dataA[2,2] = pad(b * (HxT * iMx * Hx .+ HyT * iMy * Hy) .- χ * a1)
    A  = [dataA[1,1] dataA[1,2];
          dataA[2,1] dataA[2,2]]

    dataB = Matrix{SparseMatrixCSC{Float64, Int64}}(undef, 2, 2)
    dataB[1,1] = M .+ 0.5 .* τ .* LT
    if convection
       dataB[1,1] .-= τ .* CT
    end
    dataB[1,2] = 0.5 .* τ .* LD
    dataB[2,1] = spdiagm(0 => zeros(nx*ny))
    dataB[2,2] = spdiagm(0 => zeros(nx*ny))
    B  = [dataB[1,1] dataB[1,2];
          dataB[2,1] dataB[2,2]]

    rhs = zeros(2*nx*ny)
    if convection
        veci(rhs,grid,1) .-= τ .* CUTCT
    end
    veci(rhs,grid,2) .+= χ * vec(a0)

    return A, B, rhs
end
