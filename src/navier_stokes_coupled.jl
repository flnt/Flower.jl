"""
Prepare boundary conditions

!!! danger "Sign of BC"
    * BC.left.val .* ones(ny)
    * BC.bottom.val .* ones(nx)
    * BC.right.val .* ones(ny)
    * BC.bottom.val .* ones(nx)
    In the current implementation, the sign needs to be checked.
    Dirichlet: a1 = -1
    Neumann: b =1
    Robin: a1 = -1 b = 1
    a0 : BC value 
"""
function set_borders!(grid, cl, u, a0, a1, b, BC, n_ext)
    @unpack nx, ny, x, y, dx, dy, ind = grid

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
        idx_cl = intersect(CartesianIndices((idx,1)), cl)
        idx_no = setdiff(CartesianIndices((idx,1)), cl)
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
        # idx_cl = intersect(CartesianIndices((idx,1)), cl)
        # idx_no = setdiff(CartesianIndices((idx,1)), cl)

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
        idx_cl = intersect(CartesianIndices((1,idx)), cl)
        idx_no = setdiff(CartesianIndices((1,idx)), cl)
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
        # @inbounds a1[intersect(idx, cl)] .= -1.0
        # @inbounds a1[setdiff(idx, cl)] .= -1.0
        # @inbounds b[intersect(idx, cl)] .= BC.bottom.λ
        # @inbounds b[setdiff(idx, cl)] .= 0.0
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
        idx_cl = intersect(CartesianIndices((idx,nx)), cl)
        idx_no = setdiff(CartesianIndices((idx,nx)), cl)
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
        # @inbounds a1[intersect(idx, cl)] .= -1.0
        # @inbounds a1[setdiff(idx, cl)] .= -1.0
        # @inbounds b[intersect(idx, cl)] .= BC.right.λ
        # @inbounds b[setdiff(idx, cl)] .= 0.0
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
        idx_cl = intersect(CartesianIndices((1,idx)), cl)
        idx_no = setdiff(CartesianIndices((1,idx)), cl)
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
        # @inbounds a1[intersect(idx, cl)] .= -1.0
        # @inbounds a1[setdiff(idx, cl)] .= -1.0
        # @inbounds b[intersect(idx, cl)] .= BC.top.λ
        # @inbounds b[setdiff(idx, cl)] .= 0.0
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


"""
set iMx, Ax, Ay,... indicator χ..., Bx, By...

The mass matrix M is set by:
```julia
    M.diag .= vec(geo[end].dcap[:,:,5])
```
Mx is of size (ny,nx+1)
Mx: geo[end].dcap[II,8]
At right: Mx[δx⁺(II)] = geo[end].dcap[II,10]

My is of size (ny+1,nx)
My: geo[end].dcap[II,9]
At top: My[δy⁺(II)] = geo[end].dcap[II,11]

divergence_A!(grid, AxT, AyT, geo[end].dcap, ny, ind.all_indices, periodic_x, periodic_y)
divergence_B!(BxT, ByT, geo[end].dcap, ny, ind.all_indices)

mat_assign!(Bx, sparse(-BxT'))
mat_assign!(By, sparse(-ByT'))


bc_matrix!(grid, Hx[iLS], Hy[iLS], geo[iLS].dcap, geo_p[iLS].dcap, ny, ind.all_indices)

mat_assign_T!

periodic_bcs!(grid, Bx, By, Hx[iLS], Hy[iLS], periodic_x, periodic_y)

"""
function set_cutcell_matrices!(num, grid, geo, geo_p, opC, periodic_x, periodic_y)
    @unpack nx, ny, ind = grid
    @unpack AxT, AyT, Bx, By, BxT, ByT, Hx, Hy, HxT, HyT, M, iMx, iMy, χ = opC

    M.diag .= vec(geo[end].dcap[:,:,5])
    Mx = zeros(ny,nx+1)
    for II in ind.all_indices
        Mx[II] = geo[end].dcap[II,8]
        pII = lexicographic(II, grid.ny)
        iMx.diag[pII] = inv_weight_eps(num,Mx[II])
    end
    for II in ind.b_right[1]
        Mx[δx⁺(II)] = geo[end].dcap[II,10]
        pII = lexicographic(δx⁺(II), grid.ny)
        iMx.diag[pII] = inv_weight_eps(num,Mx[δx⁺(II)])
    end


    My = zeros(ny+1,nx)
    for II in ind.all_indices
        My[II] = geo[end].dcap[II,9]
        pII = lexicographic(II, grid.ny + 1)
        iMy.diag[pII] = inv_weight_eps(num,My[II])
    end
    for II in ind.b_top[1]
        My[δy⁺(II)] = geo[end].dcap[II,11]
        pII = lexicographic(δy⁺(II), grid.ny + 1)
        iMy.diag[pII] = inv_weight_eps(num,My[δy⁺(II)])
    end   

    # Discrete gradient and divergence operators
    divergence_A!(grid, AxT, AyT, geo[end].dcap, ny, ind.all_indices, periodic_x, periodic_y)
    divergence_B!(BxT, ByT, geo[end].dcap, ny, ind.all_indices)

    mat_assign!(Bx, sparse(-BxT'))
    mat_assign!(By, sparse(-ByT'))

    # Matrices for BCs
    for iLS in 1:num.nLS
        bc_matrix!(grid, Hx[iLS], Hy[iLS], geo[iLS].dcap, geo_p[iLS].dcap, ny, ind.all_indices)

        mat_assign_T!(HxT[iLS], sparse(Hx[iLS]'))
        mat_assign_T!(HyT[iLS], sparse(Hy[iLS]'))

        periodic_bcs!(grid, Bx, By, Hx[iLS], Hy[iLS], periodic_x, periodic_y)

        χx = (geo[iLS].dcap[:,:,3] .- geo[iLS].dcap[:,:,1]) .^ 2
        χy = (geo[iLS].dcap[:,:,4] .- geo[iLS].dcap[:,:,2]) .^ 2
        χ[iLS].diag .= sqrt.(vec(χx .+ χy))
    end

    mat_assign!(BxT, sparse(-Bx'))
    mat_assign!(ByT, sparse(-By'))

    return nothing
end


"""
set_other_cutcell_matrices!

set Bx, By, Gx, Gy, Rx, Ry...
"""
function set_other_cutcell_matrices!(
    num, grid, geo, geo_u, geo_v,
    opC_p, opC_u, opC_v,
    periodic_x, periodic_y
    )
    @unpack nx, ny, ind = grid
    @unpack Bx, By, Gx, Gy = opC_p

    for iLS in 1:num.nLS
        bc_matrix!(grid, opC_u.Gx[iLS], opC_v.Gy[iLS], geo[iLS].dcap, geo_u[iLS].dcap, geo_v[iLS].dcap, ny, ind.all_indices)

        mat_assign_T!(Gx[iLS], sparse(opC_u.Gx[iLS]'))
        mat_assign_T!(Gy[iLS], sparse(opC_v.Gy[iLS]'))

        periodic_bcs!(grid, Bx, By, opC_u.Gx[iLS], opC_v.Gy[iLS], periodic_x, periodic_y)
    end
    periodic_bcs_R!(grid, opC_u.Rx, opC_v.Ry, periodic_x, periodic_y)

    return nothing
end

"""
set_boundary_indicator!(grid::Mesh{GridCC,T,N}, geo, geo_p, opC) where {T,N}

set indicator*length for Robin BC for scalar grid (p grid)
"""
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


"""
set_boundary_indicator!(grid::Mesh{GridFCx,T,N}, geo, geo_p, opC) where {T,N}

set indicator*length for Robin BC for u grid
"""
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


"""

set indicator*length for Robin BC for v grid
on the left: geo.dcap[II,1]
"""
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


"""
set boundary indicators, heights for the wall
"""
function set_border_matrices!(num,
    grid, geo, grid_u, geo_u, grid_v, geo_v,
    opC_p, opC_u, opC_v,
    periodic_x, periodic_y
    )

    set_boundary_indicator!(grid, geo, geo, opC_p)
    set_boundary_indicator!(grid_u, geo_u, geo, opC_u)
    set_boundary_indicator!(grid_v, geo_v, geo, opC_v)

    mass_matrix_borders!(num,grid.ind, opC_p.iMx_b, opC_p.iMy_b, opC_p.iMx_bd, opC_p.iMy_bd, geo.dcap, grid.ny)
    mass_matrix_borders!(num,grid_u.ind, opC_u.iMx_b, opC_u.iMy_b, opC_u.iMx_bd, opC_u.iMy_bd, geo_u.dcap, grid_u.ny)
    mass_matrix_borders!(num,grid_v.ind, opC_v.iMx_b, opC_v.iMy_b, opC_v.iMx_bd, opC_v.iMy_bd, geo_v.dcap, grid_v.ny)

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


"""
```julia
L = BxT * iMx * Bx .+ ByT * iMy * By
```
"""
function laplacian(opC::Operators{Float64, Int64})
    @unpack Bx, By, BxT, ByT, iMx, iMy, tmp_x, tmp_y = opC

    mul!(tmp_x, iMx, Bx)
    L = BxT * tmp_x
    mul!(tmp_y, iMy, By)
    L = L .+ ByT * tmp_y

    return L
end


"""
```julia
    bc_L[iLS]= BxT * iMx * Hx[iLS] .+ ByT * iMy * Hy[iLS]
    bc_L_b = (BxT * iMx_b * Hx_b .+ ByT * iMy_b * Hy_b)
```
"""
function laplacian_bc(opC::Operators{Float64, Int64}, nLS::Int64)
    @unpack BxT, ByT, Hx, Hy, iMx, iMy, Hx_b, Hy_b, iMx_b, iMy_b = opC

    bc_L = []
    for iLS in 1:nLS
        push!(bc_L, BxT * iMx * Hx[iLS] .+ ByT * iMy * Hy[iLS])
    end

    bc_L_b = (BxT * iMx_b * Hx_b .+ ByT * iMy_b * Hy_b)

    return bc_L, bc_L_b
end


"""    
set Laplacian and BC for p, u, v

### Inputs

- `num::Numerical{Float64, Int64}`: A structure containing numerical parameters for the simulation.
- `grid::Mesh{Flower.GridCC, Float64, Int64}`: The mesh for the pressure (cut-cell grid).
- `geo::Array{Flower.GeometricInfo{Float64}, 1}`: Geometric information for the pressure grid.
- `grid_u::Mesh{Flower.GridFCx, Float64, Int64}`: The mesh for the x-velocity component (face-centered grid in x-direction).
- `geo_u::Array{Flower.GeometricInfo{Float64}, 1}`: Geometric information for the x-velocity grid.
- `grid_v::Mesh{Flower.GridFCy, Float64, Int64}`: The mesh for the y-velocity component (face-centered grid in y-direction).
- `geo_v::Array{Flower.GeometricInfo{Float64}, 1}`: Geometric information for the y-velocity grid.
- `opC_p::Operators{Float64, Int64}`: Operators for the pressure grid.
- `opC_u::Operators{Float64, Int64}`: Operators for the x-velocity grid.
- `opC_v::Operators{Float64, Int64}`: Operators for the y-velocity grid.
- `periodic_x::Bool`: Flag indicating whether the domain is periodic in the x-direction.
- `periodic_y::Bool`: Flag indicating whether the domain is periodic in the y-direction.

"""
function set_matrices!(
    num::Numerical{Float64, Int64},
    grid::Mesh{Flower.GridCC, Float64, Int64},
    geo::Array{Flower.GeometricInfo{Float64}, 1}, 
    grid_u::Mesh{Flower.GridFCx, Float64, Int64}, 
    geo_u::Array{Flower.GeometricInfo{Float64}, 1}, 
    grid_v::Mesh{Flower.GridFCy, Float64, Int64}, 
    geo_v::Array{Flower.GeometricInfo{Float64}, 1},
    opC_p::Operators{Float64, Int64}, 
    opC_u::Operators{Float64, Int64}, 
    opC_v::Operators{Float64, Int64},
    periodic_x::Bool, 
    periodic_y::Bool
    )
    @unpack ny, ind = grid

    set_other_cutcell_matrices!(
        num, grid, geo, geo_u, geo_v,
        opC_p, opC_u, opC_v,
        periodic_x, periodic_y
    )

    set_cutcell_matrices!(num, grid, geo, geo, opC_p, periodic_x, periodic_y)
    set_cutcell_matrices!(num, grid_u, geo_u, geo, opC_u, periodic_x, periodic_y)
    set_cutcell_matrices!(num, grid_v, geo_v, geo, opC_v, periodic_x, periodic_y)

    Lp = laplacian(opC_p)
    Lu = laplacian(opC_u)
    Lv = laplacian(opC_v)

    set_border_matrices!(
        num,grid, geo[end], grid_u, geo_u[end], grid_v, geo_v[end],
        opC_p, opC_u, opC_v,
        periodic_x, periodic_y
    )

    bc_Lp, bc_Lp_b = laplacian_bc(opC_p, num.nLS)
    bc_Lu, bc_Lu_b = laplacian_bc(opC_u, num.nLS)
    bc_Lv, bc_Lv_b = laplacian_bc(opC_v, num.nLS)

    return Lp, bc_Lp, bc_Lp_b, Lu, bc_Lu, bc_Lu_b, Lv, bc_Lv, bc_Lv_b
end

function strain_rate(iLS, opC_u, opC_v, opC_p)
    data = Matrix{SparseMatrixCSC{Float64, Int64}}(undef, 2, 2)

    data[1,1] = opC_p.HxT[iLS] * (opC_u.HxT[iLS] * opC_u.iMx * opC_u.Bx .+ opC_u.HyT[iLS] * opC_u.iMy * opC_u.By)
    data[1,2] = opC_p.HxT[iLS] * (opC_u.HxT[iLS] * opC_u.iMx * opC_u.Hx[iLS] .+ opC_u.HyT[iLS] * opC_u.iMy * opC_u.Hy[iLS])
    data[2,1] = opC_p.HyT[iLS] * (opC_v.HyT[iLS] * opC_v.iMy * opC_v.By .+ opC_v.HxT[iLS] * opC_v.iMx * opC_v.Bx)
    data[2,2] = opC_p.HyT[iLS] * (opC_v.HyT[iLS] * opC_v.iMy * opC_v.Hy[iLS] .+ opC_v.HxT[iLS] * opC_v.iMx * opC_v.Hx[iLS])

    return data
end

"""
no_slip_condition!
    
"""
function no_slip_condition!(num, grid, grid_u, LS_u, grid_v, LS_v, periodic_x, periodic_y)
    #interpolate velocity from scalar grid to u and v grids
    interpolate_scalar!(grid, grid_u, grid_v, grid.V, grid_u.V, grid_v.V)

    normalx = cos.(LS_u.α)
    normaly = sin.(LS_v.α)

    grid_u.V .*= normalx
    grid_v.V .*= normaly

    if any(isnan, grid_u.V) || any(isnan, grid_v.V)

        @error("NaN no_slip_condition!")
        replace!(grid_u.V, NaN=>0.0)
        replace!(grid_v.V, NaN=>0.0)
    end

    i_u_ext, l_u_ext, b_u_ext, r_u_ext, t_u_ext = indices_extension(grid_u, LS_u, grid_u.ind.inside, periodic_x, periodic_y)
    i_v_ext, l_v_ext, b_v_ext, r_v_ext, t_v_ext = indices_extension(grid_v, LS_v, grid_v.ind.inside, periodic_x, periodic_y)

    field_extension!(grid_u, LS_u.u, grid_u.V, i_u_ext, l_u_ext, b_u_ext, r_u_ext, t_u_ext, num.NB, periodic_x, periodic_y)
    field_extension!(grid_v, LS_v.u, grid_v.V, i_v_ext, l_v_ext, b_v_ext, r_v_ext, t_v_ext, num.NB, periodic_x, periodic_y)

    return nothing
end


"""
uses vector_convection

set BC
"""
function set_convection!(
    num, grid, geo, grid_u, LS_u, grid_v, LS_v,
    u, v, op, ph, BC_u, BC_v,opC_p, opC_u, opC_v
    )
    @unpack Cu, CUTCu, Cv, CUTCv = op
    @unpack uD, vD = ph

    Du_x = zeros(grid_u)
    Du_y = zeros(grid_u)

    for iLS in 1:num.nLS
        Du_x[LS_u[iLS].MIXED] .= reshape(veci(uD,grid_u,iLS+1), grid_u)[LS_u[iLS].MIXED] #TODO PmIII
        Du_y[LS_u[iLS].MIXED] .= reshape(veci(uD,grid_u,iLS+1), grid_u)[LS_u[iLS].MIXED]
    end

    Du_x[:,1] .= vecb_L(uD,grid_u) 
    Du_y[1,:] .= vecb_B(uD,grid_u)
    Du_x[:,end] .= vecb_R(uD,grid_u)
    Du_y[end,:] .= vecb_T(uD,grid_u)

    Dv_x = zeros(grid_v)
    Dv_y = zeros(grid_v)

    for iLS in 1:num.nLS
        Dv_x[LS_v[iLS].MIXED] .= reshape(veci(vD,grid_v,iLS+1), grid_v)[LS_v[iLS].MIXED] #TODO PmIII
        Dv_y[LS_v[iLS].MIXED] .= reshape(veci(vD,grid_v,iLS+1), grid_v)[LS_v[iLS].MIXED]
    end
   
    Dv_x[:,1] .= vecb_L(vD,grid_v)
    Dv_y[1,:] .= vecb_B(vD,grid_v)
    Dv_x[:,end] .= vecb_R(vD,grid_v)
    Dv_y[end,:] .= vecb_T(vD,grid_v)

    if num.prediction == "PmIII" #cf Brown 2001
        # and we first consider boundary conditions
        # ```math
        # \begin{align}
        # \hat{\mathbf{n}} \cdot \mathbf{u}^*|_{\partial \Omega} &= \hat{\mathbf{n}} \cdot \mathbf{u}_b^{n+1} \\
        # \hat{\mathbf{t}} \cdot \mathbf{u}^*|_{\partial \Omega} &= \hat{\mathbf{t}} \cdot (\mathbf{u}_b^{n+1} + \Delta t \nabla_h \phi^n)|_{\partial \Omega}.
        # \end{align}
        # ```
        
        #region Compute gradient
        grad_x = zeros(grid_u)
        grad_y = zeros(grid_v)
        compute_grad_T_x_T_y_array_u_v_capacities!(num, grid, grid_u, grid_v, opC_u, opC_v, grad_x, grad_y, ph.pD)
        #endregion Compute gradient

        # printstyled(color=:red, @sprintf "\n grad min max x %.2e %.2e y %.2e %.2e\n" minimum(grd_x) maximum(grd_x) minimum(grd_y) maximum(grd_y))
    
        # printstyled(color=:red, @sprintf "\n set_convection B %.2e T %.2e L %.2e R %.2e\n" maximum(abs.(vecb_B(∇ϕ_x,grid_u))) maximum(abs.(vecb_T(∇ϕ_x,grid_u))) maximum(abs.(vecb_L(∇ϕ_y,grid_v))) maximum(abs.(vecb_R(∇ϕ_y,grid_v))))
        # printstyled(color=:red, @sprintf "\n set_convection B %.2e T %.2e L %.2e R %.2e\n" maximum(grd_x[end,:]) maximum(grd_x[1,:]) maximum(grd_y[:,1]) maximum(grd_y[:,end]))

        # print("\n dt ", num.τ)
        dt = num.τ

        Du_y[1,:] .+= dt* grad_x[1,:] #vecb_B(uD,grid_u) + 
        Du_y[end,:] .+= dt* grad_x[end,:] #vecb_T(uD,grid_u) + 
        Dv_x[:,1] .+= dt* grad_y[:,1] #vecb_L(vD,grid_v) +
        Dv_x[:,end] .+= dt* grad_y[:,end] #vecb_R(vD,grid_v) + 
        
        # ∇ϕ_x .= 0.0
        # ∇ϕ_y .= 0.0

    end



    # Du_x .= reshape(vec1(uD,grid_u), grid_u)
    # Du_y .= reshape(vec1(uD,grid_u), grid_u)

    # Du_x .= reshape(vec2(uD,grid_u), grid_u)
    # Du_y .= reshape(vec2(uD,grid_u), grid_u)

    # Du_x[:,2] .= u[:,2]
    # Du_y[2,:] .= u[2,:]

    # Du_x[:,end-1] .= u[:,end-1]
    # Du_y[end-1,:] .= u[end-1,:]

    # Dv_x .= reshape(vec1(vD,grid_v), grid_v)
    # Dv_y .= reshape(vec1(vD,grid_v), grid_v)

    # Dv_x .= reshape(vec2(vD,grid_v), grid_v)
    # Dv_y .= reshape(vec2(vD,grid_v), grid_v)

    # Dv_x[:,2] .= v[:,2]
    # Dv_y[2,:] .= v[2,:]
    # Dv_x[:,end-1] .= v[:,end-1]
    # Dv_y[end-1,:] .= v[end-1,:]


    # bnds_u = [grid_u.ind.b_left[1], grid_u.ind.b_bottom[1], grid_u.ind.b_right[1], grid_u.ind.b_top[1]]
    # bnds_v = [grid_v.ind.b_left[1], grid_v.ind.b_bottom[1], grid_v.ind.b_right[1], grid_v.ind.b_top[1]]
    # Δu = [grid_u.dx[1,1] * 0.25, grid_u.dy[1,1] * 0.5, grid_u.dx[end,end] * 0.25, grid_u.dy[end,end] * 0.5]
    # Δv = [grid_v.dx[1,1] * 0.5, grid_v.dy[1,1] * 0.25, grid_v.dx[end,end] * 0.5, grid_v.dy[end,end] * 0.25]

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

    # set_bc_bnds(dir, GridFCx, Du_x, Du_y, Dv_x, Dv_y, Hu, Hv, u, v, BC_u, BC_v)
    # set_bc_bnds(dir, GridFCy, Dv_x, Dv_y, Du_x, Du_y, Hv, Hu, v, u, BC_v, BC_u)

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

"""
    FE_set_momentum_coupled(
    bc_type, num, grid, opC,
    A, B,
    L, bc_L, bc_L_b, Mm1, BC,
    ls_advection
    )

Set `u` and `v` system matrices for Forward-Euler scheme in the diffusive term in 
presence of a Navier slip BC.

```julia
    opu.HxT[iLS] * opu.iMx_b * opu.Hx_b
```

!!! todo "change bc_type to BCu and BCv"

"""
function FE_set_momentum_coupled(
    bc_type, num, gp, gu, gv,
    opp, opu, opv,
    A, B,
    Lu, bc_Lu, bc_Lu_b, Mum1, BCu,
    Lv, bc_Lv, bc_Lv_b, Mvm1, BCv,
    ls_advection
    )
    @unpack τ, Re, nLS, nNavier, visc_coeff = num

    # iRe = 1.0 / Re
    iRe = visc_coeff

    nip = gp.nx * gp.ny
    nbp = 2 * gp.nx + 2 * gp.ny

    niu = gu.nx * gu.ny
    nbu = 2 * gu.nx + 2 * gu.ny
    ntu = (nLS - nNavier + 1) * niu + nbu

    niv = gv.nx * gv.ny
    nbv = 2 * gv.nx + 2 * gv.ny
    ntv = (nLS - nNavier + 1) * niv + nbv

    rhs = zeros(ntu + ntv + nNavier * nip)

    a0_bu = zeros(nbu)
    _a1_bu = zeros(nbu)
    _b_bu = zeros(nbu)
    for iLS in 1:num.nLS
        set_borders!(gu, gu.LS[iLS].cl, gu.LS[iLS].u, a0_bu, _a1_bu, _b_bu, BCu, num.n_ext_cl)
    end
    a1_bu = Diagonal(vec(_a1_bu))
    b_bu = Diagonal(vec(_b_bu))

    a0_bv = zeros(nbv)
    _a1_bv = zeros(nbv)
    _b_bv = zeros(nbv)
    for iLS in 1:num.nLS
        set_borders!(gv, gv.LS[iLS].cl, gv.LS[iLS].u, a0_bv, _a1_bv, _b_bv, BCv, num.n_ext_cl)
    end
    a1_bv = Diagonal(vec(_a1_bv))
    b_bv = Diagonal(vec(_b_bv))

    if ls_advection
        A.nzval .= 0.0
        # Implicit part of viscous term
        A[1:niu,1:niu] = pad_crank_nicolson(opu.M .- τ .* Lu, gu, τ)
        # Contribution to implicit part of viscous term from outer boundaries
        A[1:niu,ntu-nbu+1:ntu] = - τ .* bc_Lu_b
        # Boundary conditions for outer boundaries
        A[ntu-nbu+1:ntu,1:niu] = b_bu * (opu.HxT_b * opu.iMx_b' * opu.Bx .+ opu.HyT_b * opu.iMy_b' * opu.By)
        A[ntu-nbu+1:ntu,ntu-nbu+1:ntu] = pad(b_bu * (
            opu.HxT_b * opu.iMx_bd * opu.Hx_b .+ 
            opu.HyT_b * opu.iMy_bd * opu.Hy_b
        ) .- opu.χ_b * a1_bu)

        # Implicit part of viscous term
        A[ntu+1:ntu+niv,ntu+1:ntu+niv] = pad_crank_nicolson(opv.M .- τ .* Lv, gv, τ)
        # Contribution to implicit part of viscous term from outer boundaries
        A[ntu+1:ntu+niv,ntu+ntv-nbv+1:ntu+ntv] = - τ .* bc_Lv_b
        # Boundary conditions for outer boundaries
        A[ntu+ntv-nbv+1:ntu+ntv,ntu+1:ntu+niv] = b_bv * (opv.HxT_b * opv.iMx_b' * opv.Bx .+ opv.HyT_b * opv.iMy_b' * opv.By)
        A[ntu+ntv-nbv+1:ntu+ntv,ntu+ntv-nbv+1:ntu+ntv] = pad(b_bv * (
            opv.HxT_b * opv.iMx_bd * opv.Hx_b .+ 
            opv.HyT_b * opv.iMy_bd * opv.Hy_b
        ) .- opv.χ_b * a1_bv)

        B[1:niu,1:niu] = Mum1
        B[ntu+1:ntu+niv,ntu+1:ntu+niv] = Mvm1
    end

    nNav1 = 0
    _iLS = 1
    for iLS in 1:num.nLS
        #TODO can be improved for readability / compilation:
        # "_a1u = ones(gu) .* __a1u
        # a1u = Diagonal(vec(_a1u))
        # _bu = ones(gu) .* __bu
        # bu = Diagonal(vec(_bu)) " is the same


        if is_dirichlet(bc_type[iLS])
            velu = copy(gu.V)
            a0u = ones(gu) .* velu
            __a1u = -1.0
            __bu = 0.0
            _a1u = ones(gu) .* __a1u
            a1u = Diagonal(vec(_a1u))
            _bu = ones(gu) .* __bu
            bu = Diagonal(vec(_bu))

            velv = copy(gv.V)
            a0v = ones(gv) .* velv
            __a1v = -1.0
            __bv = 0.0
            _a1v = ones(gv) .* __a1v
            a1v = Diagonal(vec(_a1v))
            _bv = ones(gv) .* __bv
            bv = Diagonal(vec(_bv))
        elseif is_neumann(bc_type[iLS])
            velu = 0.0
            a0u = ones(gu) .* velu
            __a1u = 0.0
            __bu = 1.0
            _a1u = ones(gu) .* __a1u
            a1u = Diagonal(vec(_a1u))
            _bu = ones(gu) .* __bu
            bu = Diagonal(vec(_bu))

            velv = 0.0
            a0v = ones(gv) .* velv
            __a1v = 0.0
            __bv = 1.0
            _a1v = ones(gv) .* __a1v
            a1v = Diagonal(vec(_a1v))
            _bv = ones(gv) .* __bv
            bv = Diagonal(vec(_bv))
        elseif is_robin(bc_type[iLS])
            velu = 0.0
            a0u = ones(gu) .* velu
            __a1u = -1.0
            __bu = 1.0
            _a1u = ones(gu) .* __a1u
            a1u = Diagonal(vec(_a1u))
            _bu = ones(gu) .* __bu
            bu = Diagonal(vec(_bu))

            velv = 0.0
            a0v = ones(gv) .* velv
            __a1v = -1.0
            __bv = 1.0
            _a1v = ones(gv) .* __a1v
            a1v = Diagonal(vec(_a1v))
            _bv = ones(gv) .* __bv
            bv = Diagonal(vec(_bv))
        elseif is_fs(bc_type[iLS])
            velu = 0.0
            a0u = ones(gu) .* velu
            __a1u = 0.0
            __bu = 1.0
            _a1u = ones(gu) .* __a1u
            a1u = Diagonal(vec(_a1u))
            _bu = ones(gu) .* __bu
            bu = Diagonal(vec(_bu))

            velv = 0.0
            a0v = ones(gv) .* velv
            __a1v = 0.0
            __bv = 1.0
            _a1v = ones(gv) .* __a1v
            a1v = Diagonal(vec(_a1v))
            _bv = ones(gv) .* __bv
            bv = Diagonal(vec(_bv))
        elseif is_wall_no_slip(bc_type[iLS])
            velu = bc_type[iLS].val
            a0u = ones(gu) .* velu
            __a1u = -1.0
            __bu = 0.0
            _a1u = ones(gu) .* __a1u
            a1u = Diagonal(vec(_a1u))
            _bu = ones(gu) .* __bu
            bu = Diagonal(vec(_bu))

            velv = bc_type[iLS].val
            a0v = ones(gv) .* velv
            __a1v = -1.0
            __bv = 0.0
            _a1v = ones(gv) .* __a1v
            a1v = Diagonal(vec(_a1v))
            _bv = ones(gv) .* __bv
            bv = Diagonal(vec(_bv))
        elseif is_navier_cl(bc_type[iLS])
            velp = bc_type[iLS].val
            a0p = ones(gp) .* velp

            # Assumes FreeSurface() is the first interface!!! (LS[1])
            ϵb = num.n_ext_cl .* sqrt.(gp.dx.^2 .+ gp.dy.^2)
            bp = Diagonal(vec(bc_type[iLS].λ * bell_function2.(gp.LS[1].u, ϵb)))
        elseif is_navier(bc_type[iLS])
            velp = bc_type[iLS].val
            a0p = ones(gp) .* velp

            λ = bc_type[iLS].λ .* ones(gp)
            bp = Diagonal(vec(λ))
        else
            velu = bc_type[iLS].val
            a0u = ones(gu) .* velu
            __a1u = -1.0
            __bu = 0.0
            _a1u = ones(gu) .* __a1u
            a1u = Diagonal(vec(_a1u))
            _bu = ones(gu) .* __bu
            bu = Diagonal(vec(_bu))

            velv = bc_type[iLS].val
            a0v = ones(gv) .* velv
            __a1v = -1.0
            __bv = 0.0
            _a1v = ones(gv) .* __a1v
            a1v = Diagonal(vec(_a1v))
            _bv = ones(gv) .* __bv
            bv = Diagonal(vec(_bv))
        end

        sbu = _iLS*niu+1:(_iLS+1)*niu
        sbv = ntu+_iLS*niv+1:ntu+(_iLS+1)*niv

        if ls_advection
            if !is_navier_cl(bc_type[iLS]) && !is_navier(bc_type[iLS])
                # Contribution to implicit part of viscous term from inner boundaries
                A[1:niu,sbu] = - τ .* bc_Lu[iLS]
                A[ntu+1:ntu+niv,sbv] = - τ .* bc_Lv[iLS]
                # Boundary conditions for inner boundaries
                A[sbu,1:niu] = bu * (opu.HxT[iLS] * opu.iMx * opu.Bx .+ opu.HyT[iLS] * opu.iMy * opu.By)
                A[sbv,ntu+1:ntu+niv] = bv * (opv.HxT[iLS] * opv.iMx * opv.Bx .+ opv.HyT[iLS] * opv.iMy * opv.By)
                # Contribution to Neumann BC from other boundaries
                nNav2 = 0
                for i in 1:num.nLS
                    if i != iLS && (!is_navier_cl(bc_type[i]) && !is_navier(bc_type[i]))
                        A[sbu,i*niu+1:(i+1)*niu] = bu * (
                            opu.HxT[iLS] * opu.iMx * opu.Hx[i] .+
                            opu.HyT[iLS] * opu.iMy * opu.Hy[i]
                        )
                        A[sbv,ntu+i*niv+1:ntu+(i+1)*niv] = bv * (
                            opv.HxT[iLS] * opv.iMx * opv.Hx[i] .+
                            opv.HyT[iLS] * opv.iMy * opv.Hy[i]
                        )
                    elseif i != iLS
                        sinα = Diagonal(vec(sin.(gp.LS[i].α)))
                        cosα = Diagonal(vec(cos.(gp.LS[i].α)))

                        if any(isnan, sinα.diag) || any(isnan, cosα.diag)
                            @error("NaN FE_set_momentum_coupled")
                            replace!(sinα.diag, NaN=>0.0)
                            replace!(cosα.diag, NaN=>0.0)
                        end

                        cap1 = hcat(
                            zeros(gp.ny),
                            gp.LS[i].geoL.cap[:,1:end-1,3] .- gp.LS[i].geoL.cap[:,1:end-1,1],
                            zeros(gp.ny)
                        )
                        cap2 = vcat(
                            zeros(1,gp.nx),
                            gp.LS[i].geoL.cap[1:end-1,:,4] .- gp.LS[i].geoL.cap[1:end-1,:,2],
                            zeros(1,gp.nx)
                        )
                        cap3 = hcat(
                            zeros(gp.ny),
                            gp.LS[i].geoL.cap[:,2:end,3] .- gp.LS[i].geoL.cap[:,2:end,1],
                            zeros(gp.ny)
                        )
                        cap4 = vcat(
                            zeros(1,gp.nx),
                            gp.LS[i].geoL.cap[2:end,:,4] .- gp.LS[i].geoL.cap[2:end,:,2],
                            zeros(1,gp.nx)
                        )
                        capx = cap1 + cap3
                        capy = cap2 + cap4

                        avgx = copy(opp.Bx)
                        avgx.nzval .= 0.0
                        @inbounds @threads for II in gu.ind.all_indices[:,2:end-1]
                            pII = lexicographic(II, gu.ny)
                            if capx[II] > num.epsilon_dist
                                avgx[pII,pII-gu.ny] = cap1[II] * inv_weight_eps(num,capx[II])
                                avgx[pII,pII] = cap3[II] * inv_weight_eps(num,capx[II])
                            else
                                avgx[pII,pII-gu.ny] = 0.5
                                avgx[pII,pII] = 0.5
                            end
                        end
                        @inbounds @threads for II in gu.ind.all_indices[:,1]
                            pII = lexicographic(II, gu.ny)
                            avgx[pII,pII] = 0.5
                        end
                        @inbounds @threads for II in gu.ind.all_indices[:,end]
                            pII = lexicographic(II, gu.ny)
                            avgx[pII,pII-gu.ny] = 0.5
                        end
                        avgy = copy(opp.By)
                        avgy.nzval .= 0.0
                        @inbounds @threads for II in gv.ind.all_indices[2:end-1,:]
                            pII = lexicographic(II, gp.ny)
                            pJJ = lexicographic(II, gv.ny)
                            if capy[II] > num.epsilon_dist
                                avgy[pJJ,pII-1] = cap2[II] * inv_weight_eps(num,capy[II])
                                avgy[pJJ,pII] = cap4[II] * inv_weight_eps(num,capy[II])
                            else
                                avgy[pJJ,pII-1] = 0.5
                                avgy[pJJ,pII] = 0.5
                            end
                        end
                        @inbounds @threads for II in gv.ind.all_indices[1,:]
                            pII = lexicographic(II, gp.ny)
                            pJJ = lexicographic(II, gv.ny)
                            avgy[pJJ,pII] = 0.5
                        end
                        @inbounds @threads for II in gv.ind.all_indices[end,:]
                            pII = lexicographic(II, gp.ny)
                            pJJ = lexicographic(II, gv.ny)
                            avgy[pJJ,pII-1] = 0.5
                        end
                        A[sbu,ntu+ntv+1+nNav2*nip:ntu+ntv+(nNav2+1)*nip] = bu * (
                            opu.HxT[iLS] * opu.iMx * opu.Hx[i] .+
                            opu.HyT[iLS] * opu.iMy * opu.Hy[i]
                        ) * avgx * sinα
                        A[sbv,ntu+ntv+1+nNav2*nip:ntu+ntv+(nNav2+1)*nip] = bv * (
                            opv.HxT[iLS] * opv.iMx * opv.Hx[i] .+
                            opv.HyT[iLS] * opv.iMy * opv.Hy[i]
                        ) * avgy * (-cosα)

                        nNav2 += 1
                    end
                end
                A[sbu,sbu] = pad(bu * (
                    opu.HxT[iLS] * opu.iMx * opu.Hx[iLS] .+
                    opu.HyT[iLS] * opu.iMy * opu.Hy[iLS]
                ) .- opu.χ[iLS] * a1u)
                A[sbv,sbv] = pad(bv * (
                    opv.HxT[iLS] * opv.iMx * opv.Hx[iLS] .+
                    opv.HyT[iLS] * opv.iMy * opv.Hy[iLS]
                ) .- opv.χ[iLS] * a1v)

                A[sbu,ntu-nbu+1:ntu] = bu * (
                    opu.HxT[iLS] * opu.iMx_b * opu.Hx_b .+ opu.HyT[iLS] * opu.iMy_b * opu.Hy_b
                )
                A[sbv,ntu+ntv-nbv+1:ntu+ntv] = bv * (
                    opv.HxT[iLS] * opv.iMx_b * opv.Hx_b .+ opv.HyT[iLS] * opv.iMy_b * opv.Hy_b
                )
                # Boundary conditions for outer boundaries
                A[ntu-nbu+1:ntu,sbu] = b_bu * (
                    opu.HxT_b * opu.iMx_b' * opu.Hx[iLS] .+ opu.HyT_b * opu.iMy_b' * opu.Hy[iLS]
                )
                A[ntu+ntv-nbv+1:ntu+ntv,sbv] = b_bv * (
                    opv.HxT_b * opv.iMx_b' * opv.Hx[iLS] .+ opv.HyT_b * opv.iMy_b' * opv.Hy[iLS]
                )

                _iLS += 1
            else # Tangential component of velocity if Navier BC #if !is_navier_cl(bc_type[iLS]) && !is_navier(bc_type[iLS])
                sinα_p = Diagonal(vec(sin.(gp.LS[iLS].α)))
                # replace!(sinα_p.diag, NaN=>0.0)
                cosα_p = Diagonal(vec(cos.(gp.LS[iLS].α)))
                # replace!(cosα_p.diag, NaN=>0.0)

                sinα_u = Diagonal(vec(sin.(gu.LS[iLS].α)))
                # replace!(sinα_u.diag, NaN=>0.0)
                cosα_v = Diagonal(vec(cos.(gv.LS[iLS].α)))
                # replace!(cosα_v.diag, NaN=>0.0)

                if any(isnan, sinα_p.diag) || any(isnan, cosα_p.diag) || any(isnan, sinα_u.diag) || any(isnan, cosα_v.diag)
                    @error("NaN FE_set_momentum_coupled")
                    replace!(sinα.diag, NaN=>0.0)
                    replace!(cosα.diag, NaN=>0.0)
                    replace!(sinα_u.diag, NaN=>0.0)
                    replace!(cosα_v.diag, NaN=>0.0)
                end


                # Contribution to implicit part of viscous term from inner boundaries
                cap1 = hcat(
                    zeros(gp.ny),
                    gp.LS[iLS].geoL.cap[:,1:end-1,3] .- gp.LS[iLS].geoL.cap[:,1:end-1,1],
                    zeros(gp.ny)
                )
                cap2 = vcat(
                    zeros(1,gp.nx),
                    gp.LS[iLS].geoL.cap[1:end-1,:,4] .- gp.LS[iLS].geoL.cap[1:end-1,:,2],
                    zeros(1,gp.nx)
                )
                cap3 = hcat(
                    zeros(gp.ny),
                    gp.LS[iLS].geoL.cap[:,2:end,3] .- gp.LS[iLS].geoL.cap[:,2:end,1],
                    zeros(gp.ny)
                )
                cap4 = vcat(
                    zeros(1,gp.nx),
                    gp.LS[iLS].geoL.cap[2:end,:,4] .- gp.LS[iLS].geoL.cap[2:end,:,2],
                    zeros(1,gp.nx)
                )
                capx = cap1 + cap3
                capy = cap2 + cap4

                #Compute the averaging coefficients for 
                #if empty: 0.5 otherwise infty
                avgx = copy(opp.Bx)
                avgx.nzval .= 0.0
                @inbounds @threads for II in gu.ind.all_indices[:,2:end-1]
                    pII = lexicographic(II, gu.ny)
                    if capx[II] > num.epsilon_dist
                        avgx[pII,pII-gu.ny] = cap1[II] * inv_weight_eps(num,capx[II])
                        avgx[pII,pII] = cap3[II] * inv_weight_eps(num,capx[II])
                    else
                        if gp.LS[iLS].geoL.cap[δx⁻(II),5] > num.epsilon_vol && gp.LS[iLS].geoL.cap[II,5] > num.epsilon_vol
                            if !(δx⁻(II) in gp.LS[1].MIXED)
                                avgx[pII,pII] = 1.0
                            elseif !(II in gp.LS[1].MIXED)
                                avgx[pII,pII-gu.ny] = 1.0 #averaging coefficients
                            else
                                avgx[pII,pII-gu.ny] = 0.5
                                avgx[pII,pII] = 0.5
                            end
                        elseif gp.LS[iLS].geoL.cap[δx⁻(II),5] > num.epsilon_vol
                            avgx[pII,pII-gu.ny] = 1.0
                        else
                            avgx[pII,pII] = 1.0
                        end
                    end
                end
                @inbounds @threads for II in gu.ind.all_indices[:,1]
                    pII = lexicographic(II, gu.ny)
                    avgx[pII,pII] = 0.5
                end
                @inbounds @threads for II in gu.ind.all_indices[:,end]
                    pII = lexicographic(II, gu.ny)
                    avgx[pII,pII-gu.ny] = 0.5
                end
                avgy = copy(opp.By)
                avgy.nzval .= 0.0
                @inbounds @threads for II in gv.ind.all_indices[2:end-1,:]
                    pII = lexicographic(II, gp.ny)
                    pJJ = lexicographic(II, gv.ny)
                    if capy[II] > num.epsilon_dist
                        avgy[pJJ,pII-1] = cap2[II] * inv_weight_eps(num,capy[II])
                        avgy[pJJ,pII] = cap4[II] * inv_weight_eps(num,capy[II])
                    else
                        if gp.LS[iLS].geoL.cap[δy⁻(II),5] > num.epsilon_vol && gp.LS[iLS].geoL.cap[II,5] > num.epsilon_vol
                            avgy[pJJ,pII-1] = 0.5
                            avgy[pJJ,pII] = 0.5
                            if !(δy⁻(II) in gp.LS[1].MIXED)
                                avgy[pJJ,pII] = 1.0
                            elseif !(II in gp.LS[1].MIXED)
                                avgy[pJJ,pII-1] = 1.0
                            else
                                avgy[pJJ,pII-1] = 0.5
                                avgy[pJJ,pII] = 0.5
                            end
                        elseif gp.LS[iLS].geoL.cap[δy⁻(II),5] > num.epsilon_vol
                            avgy[pJJ,pII-1] = 1.0
                        else
                            avgy[pJJ,pII] = 1.0
                        end
                    end
                end
                @inbounds @threads for II in gv.ind.all_indices[1,:]
                    pII = lexicographic(II, gp.ny)
                    pJJ = lexicographic(II, gv.ny)
                    avgy[pJJ,pII] = 0.5
                end
                @inbounds @threads for II in gv.ind.all_indices[end,:]
                    pII = lexicographic(II, gp.ny)
                    pJJ = lexicographic(II, gv.ny)
                    avgy[pJJ,pII-1] = 0.5
                end
                A[1:niu,ntu+ntv+1+nNav1*nip:ntu+ntv+(nNav1+1)*nip] = - iRe * τ .* (
                    opu.BxT * opu.iMx * opu.Hx[iLS] .+
                    opu.ByT * opu.iMy * opu.Hy[iLS]
                ) * sinα_u * avgx
                A[ntu+1:ntu+niv,ntu+ntv+1+nNav1*nip:ntu+ntv+(nNav1+1)*nip] = - iRe * τ .* (
                    opv.BxT * opv.iMx * opv.Hx[iLS] .+
                    opv.ByT * opv.iMy * opv.Hy[iLS]
                ) * (-cosα_v) * avgy

                # Boundary conditions for inner boundaries
                cap1 = gu.LS[iLS].geoL.cap[:,1:end-1,5]
                cap2 = gv.LS[iLS].geoL.cap[1:end-1,:,5]
                cap3 = gu.LS[iLS].geoL.cap[:,2:end,5]
                cap4 = gv.LS[iLS].geoL.cap[2:end,:,5]
                capx = cap1 + cap3
                capy = cap2 + cap4

                avgu = copy(opp.BxT)
                avgu.nzval .= 0.0
                avgv = copy(opp.ByT)
                avgv.nzval .= 0.0
                @inbounds @threads for II in gp.ind.all_indices
                    pII = lexicographic(II, gp.ny)
                    pJJ = lexicographic(II, gv.ny)

                    avgu[pII,pII] = cap1[II] * inv_weight_eps(num,capx[II])
                    avgu[pII,pII+gu.ny] = cap3[II] * inv_weight_eps(num,capx[II])

                    avgv[pII,pJJ] = cap2[II] * inv_weight_eps(num,capy[II])
                    avgv[pII,pJJ+1] = cap4[II] * inv_weight_eps(num,capy[II])
                end
                A[ntu+ntv+1+nNav1*nip:ntu+ntv+(nNav1+1)*nip,1:niu] = bp * (
                    opp.HxT[iLS] * opp.iMx * opp.Bx .+
                    opp.HyT[iLS] * opp.iMy * opp.By
                ) * sinα_p * avgu
                A[ntu+ntv+1+nNav1*nip:ntu+ntv+(nNav1+1)*nip,ntu+1:ntu+niv] = bp * (
                    opp.HxT[iLS] * opp.iMx * opp.Bx .+
                    opp.HyT[iLS] * opp.iMy * opp.By
                ) * (-cosα_p) * avgv

                for i in 1:num.nLS
                    if i != iLS && (!is_navier_cl(bc_type[i]) && !is_navier(bc_type[i]))
                        cap1 = gu.LS[i].geoL.cap[:,1:end-1,3] .- gu.LS[i].geoL.cap[:,1:end-1,1]
                        cap2 = gv.LS[i].geoL.cap[1:end-1,:,4] .- gv.LS[i].geoL.cap[1:end-1,:,2]
                        cap3 = gu.LS[i].geoL.cap[:,2:end,3] .- gu.LS[i].geoL.cap[:,2:end,1]
                        cap4 = gv.LS[i].geoL.cap[2:end,:,4] .- gv.LS[i].geoL.cap[2:end,:,2]
                        capx = cap1 + cap3
                        capy = cap2 + cap4

                        avgu = copy(opp.BxT)
                        avgu.nzval .= 0.0
                        avgv = copy(opp.ByT)
                        avgv.nzval .= 0.0
                        @inbounds @threads for II in gp.ind.all_indices
                            pII = lexicographic(II, gp.ny)
                            pJJ = lexicographic(II, gv.ny)

                            avgu[pII,pII] = cap1[II] * inv_weight_eps(num,capx[II])
                            avgu[pII,pII+gu.ny] = cap3[II] * inv_weight_eps(num,capx[II])

                            avgv[pII,pJJ] = cap2[II] * inv_weight_eps(num,capy[II])
                            avgv[pII,pJJ+1] = cap4[II] * inv_weight_eps(num,capy[II])
                        end
                        A[ntu+ntv+1+nNav1*nip:ntu+ntv+(nNav1+1)*nip,i*niu+1:(i+1)*niu] = bp * (
                            opp.HxT[iLS] * opp.iMx * opp.Hx[i] .+
                            opp.HyT[iLS] * opp.iMy * opp.Hy[i]
                        ) * sinα_p * avgu
                        A[ntu+ntv+1+nNav1*nip:ntu+ntv+(nNav1+1)*nip,ntu+i*niv+1:ntu+(i+1)*niv] = bp * (
                            opp.HxT[iLS] * opp.iMx * opp.Hx[i] .+
                            opp.HyT[iLS] * opp.iMy * opp.Hy[i]
                        ) * (-cosα_p) * avgv
                    end
                end
                A[ntu+ntv+1+nNav1*nip:ntu+ntv+(nNav1+1)*nip,ntu+ntv+1+nNav1*nip:ntu+ntv+(nNav1+1)*nip] = pad(bp * (
                    opp.HxT[iLS] * opp.iMx * opp.Hx[iLS] .+
                    opp.HyT[iLS] * opp.iMy * opp.Hy[iLS]
                ) .+ opp.χ[iLS])

                sinα_b = Diagonal(zeros(nbp))
                sinα_b.diag[1:gp.ny] .= sin.(gp.LS[iLS].α[:,1])
                sinα_b.diag[gp.ny+1:gp.ny+gp.nx] .= sin.(gp.LS[iLS].α[1,:])
                sinα_b.diag[gp.ny+gp.nx+1:2gp.ny+gp.nx] .= sin.(gp.LS[iLS].α[:,end])
                sinα_b.diag[2gp.ny+gp.nx+1:end] .= sin.(gp.LS[iLS].α[end,:])
                cosα_b = Diagonal(zeros(nbp))
                cosα_b.diag[1:gp.ny] .= cos.(gp.LS[iLS].α[:,1])
                cosα_b.diag[gp.ny+1:gp.ny+gp.nx] .= cos.(gp.LS[iLS].α[1,:])
                cosα_b.diag[gp.ny+gp.nx+1:2gp.ny+gp.nx] .= cos.(gp.LS[iLS].α[:,end])
                cosα_b.diag[2gp.ny+gp.nx+1:end] .= cos.(gp.LS[iLS].α[end,:])

                avgu = spdiagm(nbp, nbu, 0 => zeros(nbp), 1 => zeros(nbp-1))
                for ii in 1:gp.ny
                    avgu[ii,ii] = 1.0
                    avgu[gp.ny+gp.nx+ii,gu.ny+gu.nx+ii] = 1.0
                end
                for ii in 1:gp.nx
                    avgu[ii+gp.ny,ii+gu.ny] = 0.5
                    avgu[ii+gp.ny,ii+gu.ny+1] = 0.5
                    avgu[ii+2gp.ny+gp.nx,ii+2gu.ny+gu.nx] = 0.5
                    avgu[ii+2gp.ny+gp.nx,ii+2gu.ny+gu.nx+1] = 0.5
                end
                avgv = spdiagm(nbp, nbv, 0 => zeros(nbp), 1 => zeros(nbp-1))
                for ii in 1:gp.ny
                    avgv[ii,ii] = 0.5
                    avgv[ii,ii+1] = 0.5
                    avgv[gp.ny+gp.nx+ii,gv.ny+gv.nx+ii] = 0.5
                    avgv[gp.ny+gp.nx+ii,gv.ny+gv.nx+ii+1] = 0.5
                end
                for ii in 1:gp.nx
                    avgv[ii+gp.ny,ii+gv.ny] = 1.0
                    avgv[ii+2gp.ny+gp.nx,ii+2gv.ny+gv.nx] = 1.0
                end
                A[ntu+ntv+1+nNav1*nip:ntu+ntv+(nNav1+1)*nip,ntu-nbu+1:ntu] = bp * (
                    opp.HxT[iLS] * opp.iMx_b * opp.Hx_b .+ 
                    opp.HyT[iLS] * opp.iMy_b * opp.Hy_b
                ) * sinα_b * avgu
                A[ntu+ntv+1+nNav1*nip:ntu+ntv+(nNav1+1)*nip,ntu+ntv-nbv+1:ntu+ntv] = bp * (
                    opp.HxT[iLS] * opp.iMx_b * opp.Hx_b .+ 
                    opp.HyT[iLS] * opp.iMy_b * opp.Hy_b
                ) * (-cosα_b) * avgv
                
                # Boundary conditions for outer boundaries
                avgu = spdiagm(nbu, nbp, 0 => zeros(nbp), 1 => zeros(nbp-1))
                for ii in 1:gu.ny
                    avgu[ii,ii] = 1.0
                    avgu[gu.ny+gu.nx+ii,gp.ny+gp.nx+ii] = 1.0
                end
                avgu[gu.ny+1,gp.ny+1]  = 1.0
                avgu[gu.ny+gu.nx,gp.ny+gp.nx]  = 1.0
                avgu[2gu.ny+gu.nx+1,2gp.ny+gp.nx+1]  = 1.0
                avgu[end,end]  = 1.0
                for ii in 2:(gu.nx-1)
                    avgu[ii+gu.ny,ii+gp.ny-1] = 0.5
                    avgu[ii+gu.ny,ii+gp.ny] = 0.5
                    avgu[ii+2gu.ny+gu.nx,ii+2gp.ny+gp.nx-1] = 0.5
                    avgu[ii+2gu.ny+gu.nx,ii+2gp.ny+gp.nx] = 0.5
                end
                avgv = spdiagm(nbv, nbp, 0 => zeros(nbp), 1 => zeros(nbp-1))
                avgv[1,1]  = 1.0
                avgv[gv.ny,gp.ny]  = 1.0
                avgv[gv.ny+gv.nx+1,gp.ny+gp.nx+1]  = 1.0
                avgv[2gv.ny+gv.nx,2gp.ny+gp.nx]  = 1.0
                for ii in 2:(gv.ny-1)
                    avgv[ii,ii-1] = 0.5
                    avgv[ii,ii] = 0.5
                    avgv[gv.ny+gv.nx+ii,gp.ny+gp.nx+ii] = 0.5
                    avgv[gv.ny+gv.nx+ii,gp.ny+gp.nx+ii] = 0.5
                end
                for ii in 1:gv.nx
                    avgv[ii+gv.ny,ii+gp.ny] = 1.0
                    avgv[ii+2gv.ny+gv.nx,ii+2gp.ny+gp.nx] = 1.0
                end
                A[ntu-nbu+1:ntu,ntu+ntv+1+nNav1*nip:ntu+ntv+(nNav1+1)*nip] = b_bu * (
                    opu.HxT_b * opu.iMx_b' * opu.Hx[iLS] .+ opu.HyT_b * opu.iMy_b' * opu.Hy[iLS]
                ) * avgx * sinα_p
                A[ntu+ntv-nbv+1:ntu+ntv,ntu+ntv+1+nNav1*nip:ntu+ntv+(nNav1+1)*nip] = b_bv * (
                    opv.HxT_b * opv.iMx_b' * opv.Hx[iLS] .+ opv.HyT_b * opv.iMy_b' * opv.Hy[iLS]
                ) * avgy * (-cosα_p)
            end # if !is_navier_cl(bc_type[iLS]) && !is_navier(bc_type[iLS])
        end

        if !is_navier_cl(bc_type[iLS]) && !is_navier(bc_type[iLS])
            @inbounds rhs[sbu] .= opu.χ[iLS] * vec(a0u)
            @inbounds rhs[sbv] .= opv.χ[iLS] * vec(a0v)
        else
            @inbounds rhs[ntu+ntv+1+nNav1*nip:ntu+ntv+(nNav1+1)*nip] .= opp.χ[iLS] * vec(a0p)
            nNav1 += 1
        end
    end

    @inbounds rhs[ntu-nbu+1:ntu] .= opu.χ_b * vec(a0_bu)
    @inbounds rhs[ntu+ntv-nbv+1:ntu+ntv] .= opv.χ_b * vec(a0_bv)
    
    return rhs
end

"""
    FE_set_momentum_coupled(
    bc_type, num, grid, opC,
    A, B,
    L, bc_L, bc_L_b, Mm1, BC,
    ls_advection
    )

Set `u` and `v` system matrices for Forward-Euler scheme in the diffusive term in 
presence of a Navier slip BC.
"""
#    rhs = zeros(ntu + ntv + nNavier * nip)

function FE_set_momentum_coupled2(
    bc_type, num, gp, gu, gv,
    opp, opu, opv,
    A, B,
    rhs,
    Lu, bc_Lu, bc_Lu_b, Mum1, BCu,
    Lv, bc_Lv, bc_Lv_b, Mvm1, BCv,
    ls_advection
    )
    @unpack τ, Re, nLS, nNavier, visc_coeff = num

    # iRe = 1.0 / Re
    iRe = visc_coeff

    nip = gp.nx * gp.ny
    nbp = 2 * gp.nx + 2 * gp.ny

    niu = gu.nx * gu.ny
    nbu = 2 * gu.nx + 2 * gu.ny
    ntu = (nLS - nNavier + 1) * niu + nbu

    niv = gv.nx * gv.ny
    nbv = 2 * gv.nx + 2 * gv.ny
    ntv = (nLS - nNavier + 1) * niv + nbv

    #Reset to zero
    rhs = 0.0

    a0_bu = zeros(nbu)
    _a1_bu = zeros(nbu)
    _b_bu = zeros(nbu)
    for iLS in 1:num.nLS
        set_borders!(gu, gu.LS[iLS].cl, gu.LS[iLS].u, a0_bu, _a1_bu, _b_bu, BCu, num.n_ext_cl)
    end
    a1_bu = Diagonal(vec(_a1_bu))
    b_bu = Diagonal(vec(_b_bu))

    a0_bv = zeros(nbv)
    _a1_bv = zeros(nbv)
    _b_bv = zeros(nbv)
    for iLS in 1:num.nLS
        set_borders!(gv, gv.LS[iLS].cl, gv.LS[iLS].u, a0_bv, _a1_bv, _b_bv, BCv, num.n_ext_cl)
    end
    a1_bv = Diagonal(vec(_a1_bv))
    b_bv = Diagonal(vec(_b_bv))

    if ls_advection
        A.nzval .= 0.0
        # Implicit part of viscous term
        A[1:niu,1:niu] = pad_crank_nicolson(opu.M .- τ .* Lu, gu, τ)
        # Contribution to implicit part of viscous term from outer boundaries
        A[1:niu,ntu-nbu+1:ntu] = - τ .* bc_Lu_b
        # Boundary conditions for outer boundaries
        A[ntu-nbu+1:ntu,1:niu] = b_bu * (opu.HxT_b * opu.iMx_b' * opu.Bx .+ opu.HyT_b * opu.iMy_b' * opu.By)
        A[ntu-nbu+1:ntu,ntu-nbu+1:ntu] = pad(b_bu * (
            opu.HxT_b * opu.iMx_bd * opu.Hx_b .+ 
            opu.HyT_b * opu.iMy_bd * opu.Hy_b
        ) .- opu.χ_b * a1_bu)

        # Implicit part of viscous term
        A[ntu+1:ntu+niv,ntu+1:ntu+niv] = pad_crank_nicolson(opv.M .- τ .* Lv, gv, τ)
        # Contribution to implicit part of viscous term from outer boundaries
        A[ntu+1:ntu+niv,ntu+ntv-nbv+1:ntu+ntv] = - τ .* bc_Lv_b
        # Boundary conditions for outer boundaries
        A[ntu+ntv-nbv+1:ntu+ntv,ntu+1:ntu+niv] = b_bv * (opv.HxT_b * opv.iMx_b' * opv.Bx .+ opv.HyT_b * opv.iMy_b' * opv.By)
        A[ntu+ntv-nbv+1:ntu+ntv,ntu+ntv-nbv+1:ntu+ntv] = pad(b_bv * (
            opv.HxT_b * opv.iMx_bd * opv.Hx_b .+ 
            opv.HyT_b * opv.iMy_bd * opv.Hy_b
        ) .- opv.χ_b * a1_bv)

        B[1:niu,1:niu] = Mum1
        B[ntu+1:ntu+niv,ntu+1:ntu+niv] = Mvm1
    end

    nNav1 = 0
    _iLS = 1
    for iLS in 1:num.nLS
        #TODO can be improved for readability / compilation:
        # "_a1u = ones(gu) .* __a1u
        # a1u = Diagonal(vec(_a1u))
        # _bu = ones(gu) .* __bu
        # bu = Diagonal(vec(_bu)) " is the same


        _a1u = ones(gu) .* __a1u
        a1u = Diagonal(vec(_a1u))
        _bu = ones(gu) .* __bu
        bu = Diagonal(vec(_bu))

        # allocates ones...
        # better tot do mapping function if cte... else if matrix ...


        if is_dirichlet(bc_type[iLS])
            velu = copy(gu.V) #TODO really reed to copy ? allocates a whole new 
            a0u = ones(gu) .* velu
            __a1u = -1.0
            __bu = 0.0
            

            velv = copy(gv.V)
            a0v = ones(gv) .* velv
            __a1v = -1.0
            __bv = 0.0
            _a1v = ones(gv) .* __a1v
            a1v = Diagonal(vec(_a1v))
            _bv = ones(gv) .* __bv
            bv = Diagonal(vec(_bv))
        elseif is_neumann(bc_type[iLS])
            velu = 0.0
            a0u = ones(gu) .* velu
            __a1u = 0.0
            __bu = 1.0
            _a1u = ones(gu) .* __a1u
            a1u = Diagonal(vec(_a1u))
            _bu = ones(gu) .* __bu
            bu = Diagonal(vec(_bu))

            velv = 0.0
            a0v = ones(gv) .* velv
            __a1v = 0.0
            __bv = 1.0
            _a1v = ones(gv) .* __a1v
            a1v = Diagonal(vec(_a1v))
            _bv = ones(gv) .* __bv
            bv = Diagonal(vec(_bv))
        elseif is_robin(bc_type[iLS])
            velu = 0.0
            a0u = ones(gu) .* velu
            __a1u = -1.0
            __bu = 1.0
            _a1u = ones(gu) .* __a1u
            a1u = Diagonal(vec(_a1u))
            _bu = ones(gu) .* __bu
            bu = Diagonal(vec(_bu))

            velv = 0.0
            a0v = ones(gv) .* velv
            __a1v = -1.0
            __bv = 1.0
            _a1v = ones(gv) .* __a1v
            a1v = Diagonal(vec(_a1v))
            _bv = ones(gv) .* __bv
            bv = Diagonal(vec(_bv))
        elseif is_fs(bc_type[iLS])
            velu = 0.0
            a0u = ones(gu) .* velu
            __a1u = 0.0
            __bu = 1.0
            _a1u = ones(gu) .* __a1u
            a1u = Diagonal(vec(_a1u))
            _bu = ones(gu) .* __bu
            bu = Diagonal(vec(_bu))

            velv = 0.0
            a0v = ones(gv) .* velv
            __a1v = 0.0
            __bv = 1.0
            _a1v = ones(gv) .* __a1v
            a1v = Diagonal(vec(_a1v))
            _bv = ones(gv) .* __bv
            bv = Diagonal(vec(_bv))
        elseif is_wall_no_slip(bc_type[iLS])
            velu = bc_type[iLS].val
            a0u = ones(gu) .* velu
            __a1u = -1.0
            __bu = 0.0
            _a1u = ones(gu) .* __a1u
            a1u = Diagonal(vec(_a1u))
            _bu = ones(gu) .* __bu
            bu = Diagonal(vec(_bu))

            velv = bc_type[iLS].val
            a0v = ones(gv) .* velv
            __a1v = -1.0
            __bv = 0.0
            _a1v = ones(gv) .* __a1v
            a1v = Diagonal(vec(_a1v))
            _bv = ones(gv) .* __bv
            bv = Diagonal(vec(_bv))
        elseif is_navier_cl(bc_type[iLS])
            velp = bc_type[iLS].val
            a0p = ones(gp) .* velp

            # Assumes FreeSurface() is the first interface!!! (LS[1])
            ϵb = num.n_ext_cl .* sqrt.(gp.dx.^2 .+ gp.dy.^2)
            bp = Diagonal(vec(bc_type[iLS].λ * bell_function2.(gp.LS[1].u, ϵb)))
        elseif is_navier(bc_type[iLS])
            velp = bc_type[iLS].val
            a0p = ones(gp) .* velp

            λ = bc_type[iLS].λ .* ones(gp)
            bp = Diagonal(vec(λ))
        else
            velu = bc_type[iLS].val
            a0u = ones(gu) .* velu
            __a1u = -1.0
            __bu = 0.0
            _a1u = ones(gu) .* __a1u
            a1u = Diagonal(vec(_a1u))
            _bu = ones(gu) .* __bu
            bu = Diagonal(vec(_bu))

            velv = bc_type[iLS].val
            a0v = ones(gv) .* velv
            __a1v = -1.0
            __bv = 0.0
            _a1v = ones(gv) .* __a1v
            a1v = Diagonal(vec(_a1v))
            _bv = ones(gv) .* __bv
            bv = Diagonal(vec(_bv))
        end

        sbu = _iLS*niu+1:(_iLS+1)*niu
        sbv = ntu+_iLS*niv+1:ntu+(_iLS+1)*niv

        if ls_advection
            if !is_navier_cl(bc_type[iLS]) && !is_navier(bc_type[iLS])
                # Contribution to implicit part of viscous term from inner boundaries
                A[1:niu,sbu] = - τ .* bc_Lu[iLS]
                A[ntu+1:ntu+niv,sbv] = - τ .* bc_Lv[iLS]
                # Boundary conditions for inner boundaries
                A[sbu,1:niu] = bu * (opu.HxT[iLS] * opu.iMx * opu.Bx .+ opu.HyT[iLS] * opu.iMy * opu.By)
                A[sbv,ntu+1:ntu+niv] = bv * (opv.HxT[iLS] * opv.iMx * opv.Bx .+ opv.HyT[iLS] * opv.iMy * opv.By)
                # Contribution to Neumann BC from other boundaries
                nNav2 = 0
                for i in 1:num.nLS
                    if i != iLS && (!is_navier_cl(bc_type[i]) && !is_navier(bc_type[i]))
                        A[sbu,i*niu+1:(i+1)*niu] = bu * (
                            opu.HxT[iLS] * opu.iMx * opu.Hx[i] .+
                            opu.HyT[iLS] * opu.iMy * opu.Hy[i]
                        )
                        A[sbv,ntu+i*niv+1:ntu+(i+1)*niv] = bv * (
                            opv.HxT[iLS] * opv.iMx * opv.Hx[i] .+
                            opv.HyT[iLS] * opv.iMy * opv.Hy[i]
                        )
                    elseif i != iLS
                        sinα = Diagonal(vec(sin.(gp.LS[i].α)))
                        # replace!(sinα.diag, NaN=>0.0)
                        cosα = Diagonal(vec(cos.(gp.LS[i].α)))
                        # replace!(cosα.diag, NaN=>0.0)

                        if any(isnan, sinα.diag) || any(isnan, cosα.diag)
                            @error("NaN FE_set_momentum_coupled")
                            replace!(sinα.diag, NaN=>0.0)
                            replace!(cosα.diag, NaN=>0.0)
                        end

                        #TODO hcat can be improved by new 
                        cap1 = hcat(
                            zeros(gp.ny),
                            gp.LS[i].geoL.cap[:,1:end-1,3] .- gp.LS[i].geoL.cap[:,1:end-1,1],
                            zeros(gp.ny)
                        )
                        cap2 = vcat(
                            zeros(1,gp.nx),
                            gp.LS[i].geoL.cap[1:end-1,:,4] .- gp.LS[i].geoL.cap[1:end-1,:,2],
                            zeros(1,gp.nx)
                        )
                        cap3 = hcat(
                            zeros(gp.ny),
                            gp.LS[i].geoL.cap[:,2:end,3] .- gp.LS[i].geoL.cap[:,2:end,1],
                            zeros(gp.ny)
                        )
                        cap4 = vcat(
                            zeros(1,gp.nx),
                            gp.LS[i].geoL.cap[2:end,:,4] .- gp.LS[i].geoL.cap[2:end,:,2],
                            zeros(1,gp.nx)
                        )
                        capx = cap1 + cap3
                        capy = cap2 + cap4

                        avgx = copy(opp.Bx)
                        avgx.nzval .= 0.0
                        @inbounds @threads for II in gu.ind.all_indices[:,2:end-1]
                            pII = lexicographic(II, gu.ny)
                            if capx[II] > num.epsilon_dist
                                avgx[pII,pII-gu.ny] = cap1[II] * inv_weight_eps(num,capx[II])
                                avgx[pII,pII] = cap3[II] * inv_weight_eps(num,capx[II])
                            else
                                avgx[pII,pII-gu.ny] = 0.5
                                avgx[pII,pII] = 0.5
                            end
                        end
                        @inbounds @threads for II in gu.ind.all_indices[:,1]
                            pII = lexicographic(II, gu.ny)
                            avgx[pII,pII] = 0.5
                        end
                        @inbounds @threads for II in gu.ind.all_indices[:,end]
                            pII = lexicographic(II, gu.ny)
                            avgx[pII,pII-gu.ny] = 0.5
                        end
                        avgy = copy(opp.By)
                        avgy.nzval .= 0.0
                        @inbounds @threads for II in gv.ind.all_indices[2:end-1,:]
                            pII = lexicographic(II, gp.ny)
                            pJJ = lexicographic(II, gv.ny)
                            if capy[II] > num.epsilon_dist
                                avgy[pJJ,pII-1] = cap2[II] * inv_weight_eps(num,capy[II])
                                avgy[pJJ,pII] = cap4[II] * inv_weight_eps(num,capy[II])
                            else
                                avgy[pJJ,pII-1] = 0.5
                                avgy[pJJ,pII] = 0.5
                            end
                        end
                        @inbounds @threads for II in gv.ind.all_indices[1,:]
                            pII = lexicographic(II, gp.ny)
                            pJJ = lexicographic(II, gv.ny)
                            avgy[pJJ,pII] = 0.5
                        end
                        @inbounds @threads for II in gv.ind.all_indices[end,:]
                            pII = lexicographic(II, gp.ny)
                            pJJ = lexicographic(II, gv.ny)
                            avgy[pJJ,pII-1] = 0.5
                        end
                        A[sbu,ntu+ntv+1+nNav2*nip:ntu+ntv+(nNav2+1)*nip] = bu * (
                            opu.HxT[iLS] * opu.iMx * opu.Hx[i] .+
                            opu.HyT[iLS] * opu.iMy * opu.Hy[i]
                        ) * avgx * sinα
                        A[sbv,ntu+ntv+1+nNav2*nip:ntu+ntv+(nNav2+1)*nip] = bv * (
                            opv.HxT[iLS] * opv.iMx * opv.Hx[i] .+
                            opv.HyT[iLS] * opv.iMy * opv.Hy[i]
                        ) * avgy * (-cosα)

                        nNav2 += 1
                    end
                end
                A[sbu,sbu] = pad(bu * (
                    opu.HxT[iLS] * opu.iMx * opu.Hx[iLS] .+
                    opu.HyT[iLS] * opu.iMy * opu.Hy[iLS]
                ) .- opu.χ[iLS] * a1u)
                A[sbv,sbv] = pad(bv * (
                    opv.HxT[iLS] * opv.iMx * opv.Hx[iLS] .+
                    opv.HyT[iLS] * opv.iMy * opv.Hy[iLS]
                ) .- opv.χ[iLS] * a1v)

                A[sbu,ntu-nbu+1:ntu] = bu * (
                    opu.HxT[iLS] * opu.iMx_b * opu.Hx_b .+ opu.HyT[iLS] * opu.iMy_b * opu.Hy_b
                )
                A[sbv,ntu+ntv-nbv+1:ntu+ntv] = bv * (
                    opv.HxT[iLS] * opv.iMx_b * opv.Hx_b .+ opv.HyT[iLS] * opv.iMy_b * opv.Hy_b
                )
                # Boundary conditions for outer boundaries
                A[ntu-nbu+1:ntu,sbu] = b_bu * (
                    opu.HxT_b * opu.iMx_b' * opu.Hx[iLS] .+ opu.HyT_b * opu.iMy_b' * opu.Hy[iLS]
                )
                A[ntu+ntv-nbv+1:ntu+ntv,sbv] = b_bv * (
                    opv.HxT_b * opv.iMx_b' * opv.Hx[iLS] .+ opv.HyT_b * opv.iMy_b' * opv.Hy[iLS]
                )

                _iLS += 1
            else # Tangential component of velocity if Navier BC
                sinα_p = Diagonal(vec(sin.(gp.LS[iLS].α)))
                # replace!(sinα_p.diag, NaN=>0.0)
                cosα_p = Diagonal(vec(cos.(gp.LS[iLS].α)))
                # replace!(cosα_p.diag, NaN=>0.0)

                sinα_u = Diagonal(vec(sin.(gu.LS[iLS].α)))
                # replace!(sinα_u.diag, NaN=>0.0)
                cosα_v = Diagonal(vec(cos.(gv.LS[iLS].α)))
                # replace!(cosα_v.diag, NaN=>0.0)

                if any(isnan, sinα_p.diag) || any(isnan, cosα_p.diag) || any(isnan, sinα_u.diag) || any(isnan, cosα_v.diag)
                    @error("NaN FE_set_momentum_coupled")
                    replace!(sinα.diag, NaN=>0.0)
                    replace!(cosα.diag, NaN=>0.0)
                    replace!(sinα_u.diag, NaN=>0.0)
                    replace!(cosα_v.diag, NaN=>0.0)
                end

                # Contribution to implicit part of viscous term from inner boundaries
                cap1 = hcat(
                    zeros(gp.ny),
                    gp.LS[iLS].geoL.cap[:,1:end-1,3] .- gp.LS[iLS].geoL.cap[:,1:end-1,1],
                    zeros(gp.ny)
                )
                cap2 = vcat(
                    zeros(1,gp.nx),
                    gp.LS[iLS].geoL.cap[1:end-1,:,4] .- gp.LS[iLS].geoL.cap[1:end-1,:,2],
                    zeros(1,gp.nx)
                )
                cap3 = hcat(
                    zeros(gp.ny),
                    gp.LS[iLS].geoL.cap[:,2:end,3] .- gp.LS[iLS].geoL.cap[:,2:end,1],
                    zeros(gp.ny)
                )
                cap4 = vcat(
                    zeros(1,gp.nx),
                    gp.LS[iLS].geoL.cap[2:end,:,4] .- gp.LS[iLS].geoL.cap[2:end,:,2],
                    zeros(1,gp.nx)
                )
                capx = cap1 + cap3
                capy = cap2 + cap4

                avgx = copy(opp.Bx)
                avgx.nzval .= 0.0
                @inbounds @threads for II in gu.ind.all_indices[:,2:end-1]
                    pII = lexicographic(II, gu.ny)
                    if capx[II] > num.epsilon_dist
                        avgx[pII,pII-gu.ny] = cap1[II] * inv_weight_eps(num,capx[II])
                        avgx[pII,pII] = cap3[II] * inv_weight_eps(num,capx[II])
                    else
                        if gp.LS[iLS].geoL.cap[δx⁻(II),5] > num.epsilon_vol && gp.LS[iLS].geoL.cap[II,5] > num.epsilon_vol
                            if !(δx⁻(II) in gp.LS[1].MIXED)
                                avgx[pII,pII] = 1.0
                            elseif !(II in gp.LS[1].MIXED)
                                avgx[pII,pII-gu.ny] = 1.0
                            else
                                avgx[pII,pII-gu.ny] = 0.5
                                avgx[pII,pII] = 0.5
                            end
                        elseif gp.LS[iLS].geoL.cap[δx⁻(II),5] > num.epsilon_vol
                            avgx[pII,pII-gu.ny] = 1.0
                        else
                            avgx[pII,pII] = 1.0
                        end
                    end
                end
                @inbounds @threads for II in gu.ind.all_indices[:,1]
                    pII = lexicographic(II, gu.ny)
                    avgx[pII,pII] = 0.5
                end
                @inbounds @threads for II in gu.ind.all_indices[:,end]
                    pII = lexicographic(II, gu.ny)
                    avgx[pII,pII-gu.ny] = 0.5
                end
                avgy = copy(opp.By)
                avgy.nzval .= 0.0
                @inbounds @threads for II in gv.ind.all_indices[2:end-1,:]
                    pII = lexicographic(II, gp.ny)
                    pJJ = lexicographic(II, gv.ny)
                    if capy[II] > num.epsilon_dist
                        avgy[pJJ,pII-1] = cap2[II] * inv_weight_eps(num,capy[II])
                        avgy[pJJ,pII] = cap4[II] * inv_weight_eps(num,capy[II])
                    else
                        if gp.LS[iLS].geoL.cap[δy⁻(II),5] > num.epsilon_vol && gp.LS[iLS].geoL.cap[II,5] > num.epsilon_vol
                            avgy[pJJ,pII-1] = 0.5
                            avgy[pJJ,pII] = 0.5
                            if !(δy⁻(II) in gp.LS[1].MIXED)
                                avgy[pJJ,pII] = 1.0
                            elseif !(II in gp.LS[1].MIXED)
                                avgy[pJJ,pII-1] = 1.0
                            else
                                avgy[pJJ,pII-1] = 0.5
                                avgy[pJJ,pII] = 0.5
                            end
                        elseif gp.LS[iLS].geoL.cap[δy⁻(II),5] > num.epsilon_vol
                            avgy[pJJ,pII-1] = 1.0
                        else
                            avgy[pJJ,pII] = 1.0
                        end
                    end
                end
                @inbounds @threads for II in gv.ind.all_indices[1,:]
                    pII = lexicographic(II, gp.ny)
                    pJJ = lexicographic(II, gv.ny)
                    avgy[pJJ,pII] = 0.5
                end
                @inbounds @threads for II in gv.ind.all_indices[end,:]
                    pII = lexicographic(II, gp.ny)
                    pJJ = lexicographic(II, gv.ny)
                    avgy[pJJ,pII-1] = 0.5
                end
                A[1:niu,ntu+ntv+1+nNav1*nip:ntu+ntv+(nNav1+1)*nip] = - iRe * τ .* (
                    opu.BxT * opu.iMx * opu.Hx[iLS] .+
                    opu.ByT * opu.iMy * opu.Hy[iLS]
                ) * sinα_u * avgx
                A[ntu+1:ntu+niv,ntu+ntv+1+nNav1*nip:ntu+ntv+(nNav1+1)*nip] = - iRe * τ .* (
                    opv.BxT * opv.iMx * opv.Hx[iLS] .+
                    opv.ByT * opv.iMy * opv.Hy[iLS]
                ) * (-cosα_v) * avgy

                # Boundary conditions for inner boundaries
                cap1 = gu.LS[iLS].geoL.cap[:,1:end-1,5]
                cap2 = gv.LS[iLS].geoL.cap[1:end-1,:,5]
                cap3 = gu.LS[iLS].geoL.cap[:,2:end,5]
                cap4 = gv.LS[iLS].geoL.cap[2:end,:,5]
                capx = cap1 + cap3
                capy = cap2 + cap4

                avgu = copy(opp.BxT)
                avgu.nzval .= 0.0
                avgv = copy(opp.ByT)
                avgv.nzval .= 0.0
                @inbounds @threads for II in gp.ind.all_indices
                    pII = lexicographic(II, gp.ny)
                    pJJ = lexicographic(II, gv.ny)

                    avgu[pII,pII] = cap1[II] * inv_weight_eps(num,capx[II])
                    avgu[pII,pII+gu.ny] = cap3[II] * inv_weight_eps(num,capx[II])

                    avgv[pII,pJJ] = cap2[II] * inv_weight_eps(num,capy[II])
                    avgv[pII,pJJ+1] = cap4[II] * inv_weight_eps(num,capy[II])
                end
                A[ntu+ntv+1+nNav1*nip:ntu+ntv+(nNav1+1)*nip,1:niu] = bp * (
                    opp.HxT[iLS] * opp.iMx * opp.Bx .+
                    opp.HyT[iLS] * opp.iMy * opp.By
                ) * sinα_p * avgu
                A[ntu+ntv+1+nNav1*nip:ntu+ntv+(nNav1+1)*nip,ntu+1:ntu+niv] = bp * (
                    opp.HxT[iLS] * opp.iMx * opp.Bx .+
                    opp.HyT[iLS] * opp.iMy * opp.By
                ) * (-cosα_p) * avgv

                for i in 1:num.nLS
                    if i != iLS && (!is_navier_cl(bc_type[i]) && !is_navier(bc_type[i]))
                        cap1 = gu.LS[i].geoL.cap[:,1:end-1,3] .- gu.LS[i].geoL.cap[:,1:end-1,1]
                        cap2 = gv.LS[i].geoL.cap[1:end-1,:,4] .- gv.LS[i].geoL.cap[1:end-1,:,2]
                        cap3 = gu.LS[i].geoL.cap[:,2:end,3] .- gu.LS[i].geoL.cap[:,2:end,1]
                        cap4 = gv.LS[i].geoL.cap[2:end,:,4] .- gv.LS[i].geoL.cap[2:end,:,2]
                        capx = cap1 + cap3
                        capy = cap2 + cap4

                        avgu = copy(opp.BxT)
                        avgu.nzval .= 0.0
                        avgv = copy(opp.ByT)
                        avgv.nzval .= 0.0
                        @inbounds @threads for II in gp.ind.all_indices
                            pII = lexicographic(II, gp.ny)
                            pJJ = lexicographic(II, gv.ny)

                            avgu[pII,pII] = cap1[II] * inv_weight_eps(num,capx[II])
                            avgu[pII,pII+gu.ny] = cap3[II] * inv_weight_eps(num,capx[II])

                            avgv[pII,pJJ] = cap2[II] * inv_weight_eps(num,capy[II])
                            avgv[pII,pJJ+1] = cap4[II] * inv_weight_eps(num,capy[II])
                        end
                        A[ntu+ntv+1+nNav1*nip:ntu+ntv+(nNav1+1)*nip,i*niu+1:(i+1)*niu] = bp * (
                            opp.HxT[iLS] * opp.iMx * opp.Hx[i] .+
                            opp.HyT[iLS] * opp.iMy * opp.Hy[i]
                        ) * sinα_p * avgu
                        A[ntu+ntv+1+nNav1*nip:ntu+ntv+(nNav1+1)*nip,ntu+i*niv+1:ntu+(i+1)*niv] = bp * (
                            opp.HxT[iLS] * opp.iMx * opp.Hx[i] .+
                            opp.HyT[iLS] * opp.iMy * opp.Hy[i]
                        ) * (-cosα_p) * avgv
                    end
                end
                A[ntu+ntv+1+nNav1*nip:ntu+ntv+(nNav1+1)*nip,ntu+ntv+1+nNav1*nip:ntu+ntv+(nNav1+1)*nip] = pad(bp * (
                    opp.HxT[iLS] * opp.iMx * opp.Hx[iLS] .+
                    opp.HyT[iLS] * opp.iMy * opp.Hy[iLS]
                ) .+ opp.χ[iLS])

                sinα_b = Diagonal(zeros(nbp))
                sinα_b.diag[1:gp.ny] .= sin.(gp.LS[iLS].α[:,1])
                sinα_b.diag[gp.ny+1:gp.ny+gp.nx] .= sin.(gp.LS[iLS].α[1,:])
                sinα_b.diag[gp.ny+gp.nx+1:2gp.ny+gp.nx] .= sin.(gp.LS[iLS].α[:,end])
                sinα_b.diag[2gp.ny+gp.nx+1:end] .= sin.(gp.LS[iLS].α[end,:])
                cosα_b = Diagonal(zeros(nbp))
                cosα_b.diag[1:gp.ny] .= cos.(gp.LS[iLS].α[:,1])
                cosα_b.diag[gp.ny+1:gp.ny+gp.nx] .= cos.(gp.LS[iLS].α[1,:])
                cosα_b.diag[gp.ny+gp.nx+1:2gp.ny+gp.nx] .= cos.(gp.LS[iLS].α[:,end])
                cosα_b.diag[2gp.ny+gp.nx+1:end] .= cos.(gp.LS[iLS].α[end,:])

                avgu = spdiagm(nbp, nbu, 0 => zeros(nbp), 1 => zeros(nbp-1))
                for ii in 1:gp.ny
                    avgu[ii,ii] = 1.0
                    avgu[gp.ny+gp.nx+ii,gu.ny+gu.nx+ii] = 1.0
                end
                for ii in 1:gp.nx
                    avgu[ii+gp.ny,ii+gu.ny] = 0.5
                    avgu[ii+gp.ny,ii+gu.ny+1] = 0.5
                    avgu[ii+2gp.ny+gp.nx,ii+2gu.ny+gu.nx] = 0.5
                    avgu[ii+2gp.ny+gp.nx,ii+2gu.ny+gu.nx+1] = 0.5
                end
                avgv = spdiagm(nbp, nbv, 0 => zeros(nbp), 1 => zeros(nbp-1))
                for ii in 1:gp.ny
                    avgv[ii,ii] = 0.5
                    avgv[ii,ii+1] = 0.5
                    avgv[gp.ny+gp.nx+ii,gv.ny+gv.nx+ii] = 0.5
                    avgv[gp.ny+gp.nx+ii,gv.ny+gv.nx+ii+1] = 0.5
                end
                for ii in 1:gp.nx
                    avgv[ii+gp.ny,ii+gv.ny] = 1.0
                    avgv[ii+2gp.ny+gp.nx,ii+2gv.ny+gv.nx] = 1.0
                end
                A[ntu+ntv+1+nNav1*nip:ntu+ntv+(nNav1+1)*nip,ntu-nbu+1:ntu] = bp * (
                    opp.HxT[iLS] * opp.iMx_b * opp.Hx_b .+ 
                    opp.HyT[iLS] * opp.iMy_b * opp.Hy_b
                ) * sinα_b * avgu
                A[ntu+ntv+1+nNav1*nip:ntu+ntv+(nNav1+1)*nip,ntu+ntv-nbv+1:ntu+ntv] = bp * (
                    opp.HxT[iLS] * opp.iMx_b * opp.Hx_b .+ 
                    opp.HyT[iLS] * opp.iMy_b * opp.Hy_b
                ) * (-cosα_b) * avgv
                
                # Boundary conditions for outer boundaries
                avgu = spdiagm(nbu, nbp, 0 => zeros(nbp), 1 => zeros(nbp-1))
                for ii in 1:gu.ny
                    avgu[ii,ii] = 1.0
                    avgu[gu.ny+gu.nx+ii,gp.ny+gp.nx+ii] = 1.0
                end
                avgu[gu.ny+1,gp.ny+1]  = 1.0
                avgu[gu.ny+gu.nx,gp.ny+gp.nx]  = 1.0
                avgu[2gu.ny+gu.nx+1,2gp.ny+gp.nx+1]  = 1.0
                avgu[end,end]  = 1.0
                for ii in 2:(gu.nx-1)
                    avgu[ii+gu.ny,ii+gp.ny-1] = 0.5
                    avgu[ii+gu.ny,ii+gp.ny] = 0.5
                    avgu[ii+2gu.ny+gu.nx,ii+2gp.ny+gp.nx-1] = 0.5
                    avgu[ii+2gu.ny+gu.nx,ii+2gp.ny+gp.nx] = 0.5
                end
                avgv = spdiagm(nbv, nbp, 0 => zeros(nbp), 1 => zeros(nbp-1))
                avgv[1,1]  = 1.0
                avgv[gv.ny,gp.ny]  = 1.0
                avgv[gv.ny+gv.nx+1,gp.ny+gp.nx+1]  = 1.0
                avgv[2gv.ny+gv.nx,2gp.ny+gp.nx]  = 1.0
                for ii in 2:(gv.ny-1)
                    avgv[ii,ii-1] = 0.5
                    avgv[ii,ii] = 0.5
                    avgv[gv.ny+gv.nx+ii,gp.ny+gp.nx+ii] = 0.5
                    avgv[gv.ny+gv.nx+ii,gp.ny+gp.nx+ii] = 0.5
                end
                for ii in 1:gv.nx
                    avgv[ii+gv.ny,ii+gp.ny] = 1.0
                    avgv[ii+2gv.ny+gv.nx,ii+2gp.ny+gp.nx] = 1.0
                end
                A[ntu-nbu+1:ntu,ntu+ntv+1+nNav1*nip:ntu+ntv+(nNav1+1)*nip] = b_bu * (
                    opu.HxT_b * opu.iMx_b' * opu.Hx[iLS] .+ opu.HyT_b * opu.iMy_b' * opu.Hy[iLS]
                ) * avgx * sinα_p
                A[ntu+ntv-nbv+1:ntu+ntv,ntu+ntv+1+nNav1*nip:ntu+ntv+(nNav1+1)*nip] = b_bv * (
                    opv.HxT_b * opv.iMx_b' * opv.Hx[iLS] .+ opv.HyT_b * opv.iMy_b' * opv.Hy[iLS]
                ) * avgy * (-cosα_p)
            end
        end

        if !is_navier_cl(bc_type[iLS]) && !is_navier(bc_type[iLS])
            @inbounds rhs[sbu] .= opu.χ[iLS] * vec(a0u)
            @inbounds rhs[sbv] .= opv.χ[iLS] * vec(a0v)
        else
            @inbounds rhs[ntu+ntv+1+nNav1*nip:ntu+ntv+(nNav1+1)*nip] .= opp.χ[iLS] * vec(a0p)
            nNav1 += 1
        end
    end

    @inbounds rhs[ntu-nbu+1:ntu] .= opu.χ_b * vec(a0_bu)
    @inbounds rhs[ntu+ntv-nbv+1:ntu+ntv] .= opv.χ_b * vec(a0_bv)
    
    return rhs
end


"""
    CN_set_momentum(
    bc_type, num, grid, opC,
    A, B,
    L, bc_L, bc_L_b,
    Lm1, bc_Lm1, bc_Lm1_b, Mm1, BC,
    ls_advection
    )

Set `u` or `v` system matrices for Crank-Nicholson scheme in the diffusive term.
"""
function CN_set_momentum(
    bc_type, num, grid, opC,
    A, B,
    L, bc_L, bc_L_b,
    Lm1, bc_Lm1, bc_Lm1_b, Mm1, BC,
    ls_advection
    )
    @unpack τ = num
    @unpack Bx, By, Hx, Hy, HxT, HyT, χ, M, iMx, iMy, Hx_b, Hy_b, HxT_b, HyT_b, iMx_b, iMy_b, iMx_bd, iMy_bd, χ_b = opC

    τ2 = 0.5 * τ

    ni = grid.nx * grid.ny
    nb = 2 * grid.nx + 2 * grid.ny

    rhs = fnzeros(grid, num)

    a0_b = zeros(nb)
    _a1_b = zeros(nb)
    _b_b = zeros(nb)
    for iLS in 1:num.nLS
        set_borders!(grid, grid.LS[iLS].cl, grid.LS[iLS].u, a0_b, _a1_b, _b_b, BC, num.n_ext_cl)
    end
    a1_b = Diagonal(vec(_a1_b))
    b_b = Diagonal(vec(_b_b))

    if ls_advection
        # Implicit part of viscous term
        A[1:ni,1:ni] = pad_crank_nicolson(M .- τ2 .* L, grid, τ)
        # Contribution to implicit part of viscous term from outer boundaries
        A[1:ni,end-nb+1:end] = - τ2 .* bc_L_b
        # Boundary conditions for outer boundaries
        A[end-nb+1:end,1:ni] = b_b * (HxT_b * iMx_b' * Bx .+ HyT_b * iMy_b' * By)
        A[end-nb+1:end,end-nb+1:end] = pad(b_b * (HxT_b * iMx_bd * Hx_b .+ HyT_b * iMy_bd * Hy_b) .- χ_b * a1_b)

        # Explicit part of viscous term
        B[1:ni,1:ni] = Mm1 .+ τ2 .* Lm1
        # Contribution to implicit part of viscous term from outer boundaries
        B[1:ni,end-nb+1:end] = τ2 .* bc_Lm1_b
    end

    for iLS in 1:num.nLS
        if is_dirichlet(bc_type[iLS])
            vel = copy(grid.V)
            __a1 = -1.0
            __b = 0.0
        elseif is_neumann(bc_type[iLS])
            vel = 0.0
            __a1 = 0.0
            __b = 1.0
        elseif is_robin(bc_type[iLS])
            vel = 0.0
            __a1 = -1.0
            __b = 1.0
        elseif is_fs(bc_type[iLS])
            vel = 0.0
            __a1 = 0.0
            __b = 1.0
        elseif is_wall_no_slip(bc_type[iLS])
            vel = bc_type[iLS].val
            __a1 = -1.0
            __b = 0.0
                else
            vel = bc_type[iLS].val
            __a1 = -1.0
            __b = 0.0
        end

        a0 = ones(grid) .* vel
        _a1 = ones(grid) .* __a1
        a1 = Diagonal(vec(_a1))
        _b = ones(grid) .* __b
        b = Diagonal(vec(_b))

        sb = iLS*ni+1:(iLS+1)*ni

        if ls_advection
            # Contribution to implicit part of viscous term from inner boundaries
            A[1:ni,sb] = - τ2 .* bc_L[iLS]
            # Boundary conditions for inner boundaries
            A[sb,1:ni] = b * (HxT[iLS] * iMx * Bx .+ HyT[iLS] * iMy * By)
            for i in 1:num.nLS
                if i != iLS
                    A[sb,i*ni+1:(i+1)*ni] = b * (HxT[iLS] * iMx * Hx[i] .+ HyT[iLS] * iMy * Hy[i])
                end
            end
            A[sb,sb] = pad(b * (HxT[iLS] * iMx * Hx[iLS] .+ HyT[iLS] * iMy * Hy[iLS]) .- χ[iLS] * a1)
            A[sb,end-nb+1:end] = b * (HxT[iLS] * iMx_b * Hx_b .+ HyT[iLS] * iMy_b * Hy_b)
            # Boundary conditions for outer boundaries
            A[end-nb+1:end,sb] = b_b * (HxT_b * iMx_b' * Hx[iLS] .+ HyT_b * iMy_b' * Hy[iLS])

            # Contribution to implicit part of viscous term from inner boundaries
            B[1:ni,sb] = τ2 .* bc_Lm1[iLS]
        end

        veci(rhs,grid,iLS+1) .= χ[iLS] * vec(a0)
    end

    vecb(rhs,grid) .= χ_b * vec(a0_b)

    return rhs
end


"""

Set `u` or `v` system matrices for Forward Euler scheme in the diffusive term.

Forward Euler
# Arguments
- bc_type: BC for interface, num, grid
- opC
- A: out
- B: out,  B[1:ni,1:ni] = Mm1
- L: Laplacian
- bc_L: Laplacian BC
- bc_L_b: Laplacian BC?
- Mm1: ?
- BC: BC for wall
- ls_advection
    FE_set_momentum(
    bc_type, num, grid, opC,
    A, B,
    L, bc_L, bc_L_b, Mm1, BC,
    ls_advection
    )

"""
function FE_set_momentum(
    num, grid, opC,
    A, B,
    L, bc_L, bc_L_b, Mm1, BC,
    ls_advection
    )
    @unpack τ = num
    @unpack Bx, By, Hx, Hy, HxT, HyT, χ, M, iMx, iMy, Hx_b, Hy_b, HxT_b, HyT_b, iMx_b, iMy_b, iMx_bd, iMy_bd, χ_b = opC

    ni = grid.nx * grid.ny
    nb = 2 * grid.nx + 2 * grid.ny

    rhs = fnzeros(grid, num)

    a0_b = zeros(nb)
    _a1_b = zeros(nb)
    _b_b = zeros(nb)
    for iLS in 1:num.nLS
        set_borders!(grid, grid.LS[iLS].cl, grid.LS[iLS].u, a0_b, _a1_b, _b_b, BC, num.n_ext_cl)
    end
    a1_b = Diagonal(vec(_a1_b))
    b_b = Diagonal(vec(_b_b))

    if ls_advection
        # Implicit part of viscous term
        A[1:ni,1:ni] = pad_crank_nicolson(M .- τ .* L, grid, τ)
        # Contribution to implicit part of viscous term from outer boundaries
        A[1:ni,end-nb+1:end] = - τ .* bc_L_b
        # Boundary conditions for outer boundaries
        A[end-nb+1:end,1:ni] = b_b * (HxT_b * iMx_b' * Bx .+ HyT_b * iMy_b' * By)
        A[end-nb+1:end,end-nb+1:end] = pad(b_b * (HxT_b * iMx_bd * Hx_b .+ HyT_b * iMy_bd * Hy_b) .- χ_b * a1_b)

        B[1:ni,1:ni] = Mm1
    end

    vecb(rhs,grid) .= χ_b * vec(a0_b)

    for iLS in 1:num.nLS
        if is_dirichlet(BC.LS[iLS])
            vel = copy(grid.V)
            __a1 = -1.0
            __b = 0.0
        elseif is_neumann(BC.LS[iLS])
            vel = 0.0
            __a1 = 0.0
            __b = 1.0
        elseif is_robin(BC.LS[iLS])
            vel = 0.0
            __a1 = -1.0
            __b = 1.0
        elseif is_fs(BC.LS[iLS])
            vel = 0.0
            __a1 = 0.0
            __b = 1.0
        elseif is_wall_no_slip(BC.LS[iLS])
            vel = BC.LS[iLS].val
            __a1 = -1.0
            __b = 0.0
        else
            vel = BC.LS[iLS].val
            __a1 = -1.0
            __b = 0.0
        end

        a0 = ones(grid) .* vel
        _a1 = ones(grid) .* __a1
        a1 = Diagonal(vec(_a1))
        _b = ones(grid) .* __b
        b = Diagonal(vec(_b))

        sb = iLS*ni+1:(iLS+1)*ni

        if ls_advection
            # Contribution to implicit part of viscous term from inner boundaries
            A[1:ni,sb] = - τ .* bc_L[iLS]
            # Boundary conditions for inner boundaries
            A[sb,1:ni] = b * (HxT[iLS] * iMx * Bx .+ HyT[iLS] * iMy * By)
            # Contribution to Neumann BC from other boundaries
            for i in 1:num.nLS
                if i != iLS
                    A[sb,i*ni+1:(i+1)*ni] = b * (HxT[iLS] * iMx * Hx[i] .+ HyT[iLS] * iMy * Hy[i])
                end
            end
            A[sb,sb] = pad(b * (HxT[iLS] * iMx * Hx[iLS] .+ HyT[iLS] * iMy * Hy[iLS]) .- χ[iLS] * a1)
            A[sb,end-nb+1:end] = b * (HxT[iLS] * iMx_b * Hx_b .+ HyT[iLS] * iMy_b * Hy_b)
            # Boundary conditions for outer boundaries
            A[end-nb+1:end,sb] = b_b * (HxT_b * iMx_b' * Hx[iLS] .+ HyT_b * iMy_b' * Hy[iLS])
        end

        veci(rhs,grid,iLS+1) .= χ[iLS] * vec(a0)

        if num.io_pdi>0
            try
                # printstyled(color=:magenta, @sprintf "\n PDI v1       FE_set_momentum %.5i %.5i \n" num.current_i num.nLS)
                #in YAML file: save only if iscal ==1 for example
                PDI_status = @ccall "libpdi".PDI_multi_expose("Navier_Stokes_set_momentum"::Cstring,
                "current_nx"::Cstring, grid.nx::Ref{Clonglong}, PDI_OUT::Cint,
                "current_ny"::Cstring, grid.ny::Ref{Clonglong}, PDI_OUT::Cint,
                "current_rhs_1D"::Cstring, rhs::Ptr{Cdouble}, PDI_OUT::Cint,
                "a0_2D"::Cstring, a0::Ptr{Cdouble}, PDI_OUT::Cint,
                C_NULL::Ptr{Cvoid})::Cint
    
            catch error
                printstyled(color=:red, @sprintf "\n PDI error \n")
                print(error)
                # print("\n PDI_status ",PDI_status)
                printstyled(color=:red, @sprintf "\n PDI error \n")
            end
        end #if io_pdi

    end
    
    #a0 not defined outside of iLS loop so cannot expose to PDI a0 after end "#if io_pdi"

    return rhs
end


function FE_set_momentum_old(
    bc_type, num, grid, opC,
    A, B,
    L, bc_L, bc_L_b, Mm1, BC,
    ls_advection
    )
    @unpack τ = num
    @unpack Bx, By, Hx, Hy, HxT, HyT, χ, M, iMx, iMy, Hx_b, Hy_b, HxT_b, HyT_b, iMx_b, iMy_b, iMx_bd, iMy_bd, χ_b = opC

    ni = grid.nx * grid.ny
    nb = 2 * grid.nx + 2 * grid.ny

    rhs = fnzeros(grid, num)

    a0_b = zeros(nb)
    _a1_b = zeros(nb)
    _b_b = zeros(nb)
    for iLS in 1:num.nLS
        set_borders!(grid, grid.LS[iLS].cl, grid.LS[iLS].u, a0_b, _a1_b, _b_b, BC, num.n_ext_cl)
    end
    a1_b = Diagonal(vec(_a1_b))
    b_b = Diagonal(vec(_b_b))

    if ls_advection
        # Implicit part of viscous term
        A[1:ni,1:ni] = pad_crank_nicolson(M .- τ .* L, grid, τ)
        # Contribution to implicit part of viscous term from outer boundaries
        A[1:ni,end-nb+1:end] = - τ .* bc_L_b
        # Boundary conditions for outer boundaries
        A[end-nb+1:end,1:ni] = b_b * (HxT_b * iMx_b' * Bx .+ HyT_b * iMy_b' * By)
        A[end-nb+1:end,end-nb+1:end] = pad(b_b * (HxT_b * iMx_bd * Hx_b .+ HyT_b * iMy_bd * Hy_b) .- χ_b * a1_b)

        B[1:ni,1:ni] = Mm1
    end

    for iLS in 1:num.nLS
        if is_dirichlet(bc_type[iLS])
            vel = copy(grid.V)
            __a1 = -1.0
            __b = 0.0
        elseif is_neumann(bc_type[iLS])
            vel = 0.0
            __a1 = 0.0
            __b = 1.0
        elseif is_robin(bc_type[iLS])
            vel = 0.0
            __a1 = -1.0
            __b = 1.0
        elseif is_fs(bc_type[iLS])
            vel = 0.0
            __a1 = 0.0
            __b = 1.0
        elseif is_wall_no_slip(bc_type[iLS])
            vel = bc_type[iLS].val
            __a1 = -1.0
            __b = 0.0
                else
            vel = bc_type[iLS].val
            __a1 = -1.0
            __b = 0.0
        end

        a0 = ones(grid) .* vel
        _a1 = ones(grid) .* __a1
        a1 = Diagonal(vec(_a1))
        _b = ones(grid) .* __b
        b = Diagonal(vec(_b))

        sb = iLS*ni+1:(iLS+1)*ni

        if ls_advection
            # Contribution to implicit part of viscous term from inner boundaries
            A[1:ni,sb] = - τ .* bc_L[iLS]
            # Boundary conditions for inner boundaries
            A[sb,1:ni] = b * (HxT[iLS] * iMx * Bx .+ HyT[iLS] * iMy * By)
            # Contribution to Neumann BC from other boundaries
            for i in 1:num.nLS
                if i != iLS
                    A[sb,i*ni+1:(i+1)*ni] = b * (HxT[iLS] * iMx * Hx[i] .+ HyT[iLS] * iMy * Hy[i])
                end
            end
            A[sb,sb] = pad(b * (HxT[iLS] * iMx * Hx[iLS] .+ HyT[iLS] * iMy * Hy[iLS]) .- χ[iLS] * a1)
            A[sb,end-nb+1:end] = b * (HxT[iLS] * iMx_b * Hx_b .+ HyT[iLS] * iMy_b * Hy_b)
            # Boundary conditions for outer boundaries
            A[end-nb+1:end,sb] = b_b * (HxT_b * iMx_b' * Hx[iLS] .+ HyT_b * iMy_b' * Hy[iLS])
        end

        veci(rhs,grid,iLS+1) .= χ[iLS] * vec(a0)
    end

    vecb(rhs,grid) .= χ_b * vec(a0_b)
    
    return rhs
end

 
@doc raw"""
# Arguments
- bc_type: BC for interface, num, grid, 
- a0, 
- opC, 
- opC_u, 
- pC_v,
- A, system matrix
- L, 
- bc_L, 
- bc_L_b, 
- BC: BC for wall (aka border)
- ls_advection



cf. [`(Rodriguez et al. 2024)`](https://link.springer.com/article/10.1007/s00707-024-04133-4) for a 1D expression in the i-th cell, x component 


```math
\begin{aligned}
-\mathcal{B}_{x,i} [&\mathcal{W}_{x,i+1} (\mathcal{B}_{x,i+1} p^\omega_{i+1} - \mathcal{B}_{x,i} p^\omega _i ) - \mathcal{W}_{x,i} (\mathcal{B}_{x,i} p^\omega_i - \mathcal{B}_{x,i-1} p^\omega_{i-1} )] \\
-\mathcal{B}_{x,i} \{&\mathcal{W}_{x,i+1} [(\mathcal{B}_{x,i+1}  - \mathcal{A}_{x,i+1} )p^\gamma_{i+1} - (\mathcal{A}_{x,i+1} - \mathcal{B}_{x,i} ) p^\gamma_i ] \\
&-\mathcal{W}_{x,i} [(\mathcal{B}_{x,i} - A_{x,i} ) p^\gamma_i + (A_{x,i} - \mathcal{B}_{x,i-1}) p^\gamma_{i-1} ] \} \\
= V_i f^\omega_i&
\end{aligned}
```

"""
function set_poisson(
    bc_type, num, grid, a0, opC, opC_u, opC_v,
    A, L, bc_L, bc_L_b, BC,
    ls_advection)
    @unpack Bx, By, Hx, Hy, HxT, HyT, χ, M, iMx, iMy, Hx_b, Hy_b, HxT_b, HyT_b, iMx_b, iMy_b, iMx_bd, iMy_bd, χ_b = opC

    ni = grid.nx * grid.ny
    nb = 2 * grid.nx + 2 * grid.ny

    rhs = fnzeros(grid, num)

    a0_b = zeros(nb)
    _a1_b = zeros(nb)
    _b_b = zeros(nb)
    for iLS in 1:num.nLS
        set_borders!(grid, grid.LS[iLS].cl, grid.LS[iLS].u, a0_b, _a1_b, _b_b, BC, num.n_ext_cl)
    end
    a1_b = Diagonal(vec(_a1_b))
    b_b = Diagonal(vec(_b_b))

    if ls_advection
        # Poisson equation
        A[1:ni,1:ni] = pad(L, -4.0)
        A[1:ni,end-nb+1:end] = bc_L_b

        # Boundary conditions for outer boundaries
        A[end-nb+1:end,1:ni] = -b_b * (HxT_b * iMx_b' * Bx .+ HyT_b * iMy_b' * By)
        A[end-nb+1:end,end-nb+1:end] = -pad(b_b * (HxT_b * iMx_bd * Hx_b .+ HyT_b * iMy_bd * Hy_b) .- χ_b * a1_b, 4.0)
    end

    for iLS in 1:num.nLS
        if ls_advection
            if is_dirichlet(bc_type[iLS])
                __a1 = -1.0
                __a2 = 0.0
                __b = 0.0
            elseif is_neumann(bc_type[iLS])
                __a1 = 0.0
                __a2 = 0.0
                __b = 1.0
            elseif is_robin(bc_type[iLS])
                __a1 = -1.0
                __a2 = 0.0
                __b = 1.0
            elseif is_fs(bc_type[iLS])
                __a1 = 0.0
                __a2 = 1.0
                __b = 0.0
            elseif is_wall_no_slip(bc_type[iLS])
                __a1 = 0.0
                __a2 = 0.0
                __b = 1.0
            elseif is_navier(bc_type[iLS])
                __a1 = 0.0
                __a2 = 0.0
                __b = 1.0
            elseif is_navier_cl(bc_type[iLS])
                __a1 = 0.0
                __a2 = 0.0
                __b = 1.0
            else
                __a1 = 0.0
                __a2 = 0.0
                __b = 1.0
            end
    
            _a1 = ones(grid) .* __a1
            a1 = Diagonal(vec(_a1))
            _a2 = ones(grid) .* __a2
            a2 = Diagonal(vec(_a2))
            _b = ones(grid) .* __b
            b = Diagonal(vec(_b))

            fs_mat = HxT[iLS] * Hx[iLS] .+ HyT[iLS] * Hy[iLS]

            sb = iLS*ni+1:(iLS+1)*ni
            
            # Poisson equation
            A[1:ni,sb] = bc_L[iLS]
            # Boundary conditions for inner boundaries
            A[sb,1:ni] = -b * (HxT[iLS] * iMx * Bx .+ HyT[iLS] * iMy * By)
            # Contribution to Neumann BC from other boundaries
            for i in 1:num.nLS
                if i != iLS
                    A[sb,i*ni+1:(i+1)*ni] = -b * (HxT[iLS] * iMx * Hx[i] .+ HyT[iLS] * iMy * Hy[i])
                end
            end
            A[sb,sb] = -pad(
                b * (HxT[iLS] * iMx * Hx[iLS] .+ HyT[iLS] * iMy * Hy[iLS]) .- χ[iLS] * a1 .+
                a2 * Diagonal(diag(fs_mat)), 4.0
            )
            A[sb,end-nb+1:end] = b * (HxT[iLS] * iMx_b * Hx_b .+ HyT[iLS] * iMy_b * Hy_b)
            # Boundary conditions for outer boundaries
            A[end-nb+1:end,sb] = -b_b * (HxT_b * iMx_b' * Hx[iLS] .+ HyT_b * iMy_b' * Hy[iLS])
        end

        veci(rhs,grid,iLS+1) .= -χ[iLS] * vec(a0[iLS])
    end

    vecb(rhs,grid) .= -χ_b * vec(a0_b)
    
    return rhs
end


"""

"""
function set_Crank_Nicolson!(
    bc_int, num, grid, geo, grid_u, geo_u, grid_v, geo_v,
    opC_p, opC_u, opC_v, BC_p, BC_u, BC_v,
    Au, Bu, Av, Bv, Aϕ,
    Lpm1, bc_Lpm1, bc_Lpm1_b, Lum1, bc_Lum1, bc_Lum1_b, Lvm1, bc_Lvm1, bc_Lvm1_b, 
    Mum1, Mvm1, iRe, op_conv, ph,
    periodic_x, periodic_y, advection, ls_advection
    )

    # @unpack rho1=num
    # irho1=1.0/rho1

    if advection
        set_convection!(num, grid, geo[end], grid_u, grid_u.LS, grid_v, grid_v.LS, ph.u, ph.v, op_conv, ph, BC_u, BC_v, opC_p, opC_u, opC_v)
    end

    if ls_advection
        update_all_ls_data(num, grid, grid_u, grid_v, bc_int, periodic_x, periodic_y, false)

        laps = set_matrices!(
            num, grid, geo, grid_u, geo_u, grid_v, geo_v,
            opC_p, opC_u, opC_v,
            periodic_x, periodic_y
        )
    else
        laps = Lpm1, bc_Lpm1, bc_Lpm1_b, Lum1, bc_Lum1, bc_Lum1_b, Lvm1, bc_Lvm1, bc_Lvm1_b
    end
    Lp, bc_Lp, bc_Lp_b, Lu, bc_Lu, bc_Lu_b, Lv, bc_Lv, bc_Lv_b = laps

    rhs_u = CN_set_momentum(
        bc_int, num, grid_u, opC_u,
        Au, Bu,
        iRe.*Lu, iRe.*bc_Lu, iRe.*bc_Lu_b,
        iRe.*Lum1, iRe.*bc_Lum1, iRe.*bc_Lum1_b, Mum1, BC_u,
        ls_advection
    )
    rhs_v = CN_set_momentum(
        bc_int, num, grid_v, opC_v,
        Av, Bv,
        iRe.*Lv, iRe.*bc_Lv, iRe.*bc_Lv_b,
        iRe.*Lvm1, iRe.*bc_Lvm1, iRe.*bc_Lvm1_b, Mvm1, BC_v,
        ls_advection
    )
    a0_p = []
    for i in 1:num.nLS
        push!(a0_p, zeros(grid))
    end
    # rhs_ϕ = set_poisson(
    #     bc_int, num, grid, a0_p, opC_p, opC_u, opC_v,
    #     Aϕ, Lp, bc_Lp, bc_Lp_b, BC_p,
    #     ls_advection
    # )

    # vecb(rhs,grid) .= +χ_b * vec(a0_b) #was - in set_poisson , -a0, a1 = -1,...

    # In solve_poisson, the equation a \frac{\partial p}{\partial n} + bp = g , 
    # if inhomogeneous Neumann: -1 at bottom and left
    # +1 sign at top and right
    
    rhs_ϕ = solve_poisson(
        bc_int, num, grid, a0_p, opC_p, opC_u, opC_v,
        Aϕ, Lp, bc_Lp, bc_Lp_b, BC_p,
        ls_advection
    )

    return rhs_u, rhs_v, rhs_ϕ, Lp, bc_Lp, bc_Lp_b, Lu, bc_Lu, bc_Lu_b, Lv, bc_Lv, bc_Lv_b
end

"""
    set_Forward_Euler!(
        bc_int, num, grid, geo, grid_u, geo_u, grid_v, geo_v,
        opC_p, opC_u, opC_v, BC_p, BC_u, BC_v,
        Au, Bu, Av, Bv, Aϕ, Auv, Buv,
        Lpm1, bc_Lpm1, bc_Lpm1_b, Lum1, bc_Lum1, bc_Lum1_b, Lvm1, bc_Lvm1, bc_Lvm1_b,
        Mum1, Mvm1, iRe, op_conv, ph,
        periodic_x, periodic_y, advection, ls_advection, navier
    )

Sets up the matrices and right-hand side (RHS) for Forward Euler (FE) 
 for the Navier-Stokes equations, optionally including advection and coupling terms.

### Arguments

- `bc_int`: Boundary conditions for the interface.
- `num`: Numerical parameters structure.
- `grid`: Grid structure.
- `geo`: Geometry structure.
- `grid_u`, `geo_u`: Grid for the x-component of velocity (u).
- `grid_v`, `geo_v`: Grid for the y-component (v).
- `opC_p`, `opC_u`, `opC_v`: Operator structures for pressure, u-component, and v-component.
- `BC_p`, `BC_u`, `BC_v`: Boundary conditions for pressure, u-component, and v-component.
- `Au`, `Bu`, `Av`, `Bv`: Matrices for the u and v components.
- `Aϕ`: Matrix for the pressure.
- `Auv`, `Buv`: Matrices for the coupled system.
- `Lpm1`, `bc_Lpm1`, `bc_Lpm1_b`: Laplacian matrix and boundary conditions for pressure.
- `Lum1`, `bc_Lum1`, `bc_Lum1_b`: Laplacian matrix and boundary conditions for u-component.
- `Lvm1`, `bc_Lvm1`, `bc_Lvm1_b`: Laplacian matrix and boundary conditions for v-component.
- `Mum1`, `Mvm1`: Mass matrices for the u and v components.
- `iRe`: Inverse of the Reynolds number.
- `op_conv`: Operator for convection.
- `ph`: Solution structure.
- `periodic_x`, `periodic_y`: Flags indicating periodic boundary conditions in x and y directions.
- `advection`: Flag indicating whether advection terms should be included.
- `ls_advection`: Flag indicating whether advection is activated.
- `navier`: Flag indicating whether the system is Navier-Stokes (coupled) or Stokes (decoupled).

### Returns

- `rhs_u`: Right-hand side vector for the u-component.
- `rhs_v`: Right-hand side vector for the v-component.
- `rhs_ϕ`: Right-hand side vector for the pressure.
- `rhs_uv`: Right-hand side vector for the coupled system (if `navier` is true).
- `Lp`, `bc_Lp`, `bc_Lp_b`: Laplacian matrix and boundary conditions for pressure.
- `Lu`, `bc_Lu`, `bc_Lu_b`: Laplacian matrix and boundary conditions for the u-component.
- `Lv`, `bc_Lv`, `bc_Lv_b`: Laplacian matrix and boundary conditions for the v-component.

### Description

1. **Advection Setup**: If advection is enabled, the convection terms are set up using the `set_convection!` function.
2. **Laplacian Matrices**: If advection is enabled, the Laplacian matrices are updated using the `set_matrices!` function. Otherwise, the Laplacian matrices are used.
3. **Right-Hand Side Vectors**:
   - For the Stokes system (`navier` is false), the right-hand side vectors for the u and v components are computed using the `FE_set_momentum` function.
   - For the Navier-Stokes system (`navier` is true), the right-hand side vector for the coupled system is computed using the `FE_set_momentum_coupled` function.
4. **Pressure Poisson Equation**: The right-hand side vector for the pressure Poisson equation is computed using the `set_poisson` function.

"""
function set_Forward_Euler!(
    bc_int, num, grid, geo, grid_u, geo_u, grid_v, geo_v,
    opC_p, opC_u, opC_v, BC_p, BC_u, BC_v,
    Au, Bu, Av, Bv, Aϕ, Auv, Buv,
    Lpm1, bc_Lpm1, bc_Lpm1_b, Lum1, bc_Lum1, bc_Lum1_b, Lvm1, bc_Lvm1, bc_Lvm1_b,
    Mum1, Mvm1, iRe, op_conv, ph,
    periodic_x, periodic_y, advection, ls_advection, navier
    )

    if advection
        set_convection!(num, grid, geo[end], grid_u, grid_u.LS, grid_v, grid_v.LS, ph.u, ph.v, op_conv, ph, BC_u, BC_v,opC_p, opC_u, opC_v)
    end

    if ls_advection
        update_all_ls_data(num, grid, grid_u, grid_v, bc_int, periodic_x, periodic_y, false)

        laps = set_matrices!(
            num, grid, geo, grid_u, geo_u, grid_v, geo_v,
            opC_p, opC_u, opC_v,
            periodic_x, periodic_y
        )
    else
        laps = Lpm1, bc_Lpm1, bc_Lpm1_b, Lum1, bc_Lum1, bc_Lum1_b, Lvm1, bc_Lvm1, bc_Lvm1_b
    end
    Lp, bc_Lp, bc_Lp_b, Lu, bc_Lu, bc_Lu_b, Lv, bc_Lv, bc_Lv_b = laps

    if !navier
        rhs_u = FE_set_momentum(
            num, grid_u, opC_u,
            Au, Bu,
            iRe.*Lu, iRe.*bc_Lu, iRe.*bc_Lu_b, Mum1, BC_u,
            ls_advection
        )
        rhs_v = FE_set_momentum(
            num, grid_v, opC_v,
            Av, Bv,
            iRe.*Lv, iRe.*bc_Lv, iRe.*bc_Lv_b, Mvm1, BC_v,
            ls_advection
        )
        rhs_uv = nothing
    else
        rhs_u = nothing
        rhs_v = nothing
        rhs_uv = FE_set_momentum_coupled(
            bc_int, num, grid, grid_u, grid_v,
            opC_p, opC_u, opC_v,
            Auv, Buv,
            iRe.*Lu, iRe.*bc_Lu, iRe.*bc_Lu_b, Mum1, BC_u,
            iRe.*Lv, iRe.*bc_Lv, iRe.*bc_Lv_b, Mvm1, BC_v,
            ls_advection
        )
    end
    a0_p = []
    for i in 1:num.nLS
        push!(a0_p, zeros(grid))
    end
    # rhs_ϕ = set_poisson(
    #     bc_int, num, grid, a0_p, opC_p, opC_u, opC_v,
    #     Aϕ, Lp, bc_Lp, bc_Lp_b, BC_p,
    #     ls_advection
    # )

    # vecb(rhs,grid) .= +χ_b * vec(a0_b) #was - in set_poisson , -a0, a1 = -1,...

    # In solve_poisson, the equation a \frac{\partial p}{\partial n} + bp = g , 
    # if inhomogeneous Neumann: -1 at bottom and left
    # +1 sign at top and right
    
    rhs_ϕ = solve_poisson(
        bc_int, num, grid, a0_p, opC_p, opC_u, opC_v,
        Aϕ, Lp, bc_Lp, bc_Lp_b, BC_p,
        ls_advection
    )

    return rhs_u, rhs_v, rhs_ϕ, rhs_uv, Lp, bc_Lp, bc_Lp_b, Lu, bc_Lu, bc_Lu_b, Lv, bc_Lv, bc_Lv_b
end


"""
solves Navier-Stokes equations with a pressure projection method. 


#### Variables and Data Structures
- `vec1(ucorrD, grid_u)`: Velocity correction for the horizontal grid.
- `vec1(vcorrD, grid_v)`: Velocity correction for the vertical grid.
- `vec1(rhs_ϕ, grid)`: Right-hand side of the Poisson equation.
- `vec1(pD, grid)`: Pressure correction.
- `vec1(uD, grid_u)`: Updated horizontal velocity.
- `vec1(vD, grid_v)`: Updated vertical velocity.
- `vec1(ϕD, grid)`: Pressure correction potential.
- `ϕ`: Pressure correction potential.
- `u`: Updated horizontal velocity.
- `v`: Updated vertical velocity.
- `p`: Pressure.
- `opC_p`, `opC_u`, `opC_v`: Operator matrices for pressure, horizontal velocity, and vertical velocity, respectively.
- `geo`, `geo_u`, `geo_v`: Geometric data for the grid.
- `bc_int`: Internal boundary conditions.
- `nLS`: Number of levelsets.
- `ntu`, `ntv`, `niu`, `niv`: Grid dimensions.
- `nbv`: Number of boundary cells in the vertical direction.
- `nNav`: Counter for Navier boundary conditions.
- `iLS`: index of levelset
- `iRe`: Reynolds number.
- `ρ1`, `ρ2`: Densities.
- `σ`: Surface tension coefficient.
- `mass_flux`: Mass flux.
- `pres_free_suface`: Free surface pressure.
- `diff_inv_rho`: Difference in inverse densities.
- `jump_mass_flux`: Flag for mass flux jump.
- `τ`: Time step.
- `Aϕ`: Matrix for the Poisson equation.
- `num`: Numerical parameters.
- `epsilon_mode`, `epsilon_vol`: parameters for epsilon handling.
- `strain_rate`: Function to compute strain rate.
- `Diagonal`: Function to create a diagonal matrix.
- `inv_weight_eps2`: Function to compute inverse weights.
- `iMu`, `iMv`: Inverse weight matrices for horizontal and vertical velocities, respectively.
- `∇ϕ_x`, `∇ϕ_y`: Gradients of the pressure correction potential.
- `iM`: Inverse weight matrix for the pressure correction.

#### Functions and Operations
1. **Initialization and Updates**
   - `vecb(vcorrD, grid_v) .= uvD[ntu+ntv-nbv+1:ntu+ntv]`: Updates the vertical velocity correction.
   - `kill_dead_cells!(vec1(vcorrD,grid_v), grid_v, geo_v[end])`: Removes dead cells from the vertical velocity correction grid.
   - `vcorr .= reshape(vec1(vcorrD,grid_v), grid_v)`: Reshapes the vertical velocity correction.

2. **Navier and Non-Navier Boundary Conditions**
   - Loop through linear solvers (`iLS`) to apply boundary conditions:
     - If not Navier or Navier-CL, update and apply boundary conditions for horizontal and vertical velocities.
     - If Navier or Navier-CL, update the Navier matrix.

3. **Divergence Calculation**
   - Calculate the divergence of the velocity corrections (`Duv`).
   - Add contributions from internal boundary conditions.

4. **Poisson Equation**
   - Set the right-hand side of the Poisson equation (`rhs_ϕ`).
   - Handle free surface conditions and Marangoni effects if `jump_mass_flux` is true.
   - Remove nullspace from the matrix `Aϕ`.
   - Apply boundary conditions and solve the Poisson equation using `Aϕ / rhs_ϕ`.

5. **Pressure Correction**
   - Update the pressure correction potential (`ϕ`).
   - Compute the gradients of the pressure correction potential (`∇ϕ_x`, `∇ϕ_y`).

6. **Velocity Correction**
   - Update the horizontal and vertical velocities (`u`, `v`) using the pressure correction gradients.
   - Apply boundary conditions and remove dead cells.

7. **Return Values**
   - Return various matrices and boundary conditions for further use.

#### Notes
- The code includes handling for free surfaces and Navier boundary conditions.
- The Poisson equation is solved using a linear solver (`Aϕ / rhs_ϕ`).
- The pressure correction is applied to update the velocities.
- Dead cells are removed from the grids to maintain numerical stability.

---

This documentation provides a high-level overview of the code's functionality and the key operations performed. For detailed implementation of specific functions or operations, refer to the corresponding sections of the code.
"""
function pressure_projection!(
    time_scheme, bc_int,
    num, grid, geo, grid_u, geo_u, grid_v, geo_v, ph,
    BC_u, BC_v, BC_p,
    opC_p, opC_u, opC_v, op_conv,
    Au, Bu, Av, Bv, Aϕ, Auv, Buv,
    Lpm1, bc_Lpm1, bc_Lpm1_b, Lum1, bc_Lum1, bc_Lum1_b, Lvm1, bc_Lvm1, bc_Lvm1_b,
    Cum1, Cvm1, Mum1, Mvm1,
    periodic_x, periodic_y, advection, ls_advection, current_i, Ra, navier, pres_free_suface,jump_mass_flux,mass_flux
    )
    @unpack Re, τ, σ, g, β, nLS, nNavier = num
    @unpack p, pD, ϕ, ϕD, u, v, ucorrD, vcorrD, uD, vD, ucorr, vcorr, uT = ph
    @unpack Cu, Cv, CUTCu, CUTCv = op_conv
    @unpack rho1,rho2,visc_coeff = num

    iRe = visc_coeff
    iτ = 1.0 / τ
    irho1 = 1.0/rho1
    mu1_over_rho1 = num.mu1 / num.rho1 

    II = CartesianIndex(div(grid_v.ny,2),1)
    pII = lexicographic(II,grid_v.ny)
    print("\nLv[pII,:] ",Lvm1[pII,:])
    print("\bc_Lvm1_b[pII,:] ",bc_Lvm1_b[pII,:])



    PDI_status = @ccall "libpdi".PDI_multi_expose("print_before_prediction"::Cstring,
    "u_1D"::Cstring, uD::Ptr{Cdouble}, PDI_OUT::Cint,
    "v_1D"::Cstring, vD::Ptr{Cdouble}, PDI_OUT::Cint,
    "p_1D"::Cstring, ph.pD::Ptr{Cdouble}, PDI_OUT::Cint,
    C_NULL::Ptr{Cvoid})::Cint

    #region prediction 
   
    #region add gradient of pressure to prediction
    # Compute gradient of pressure localized on u and v grids , times volume
    # $ \nabla p^{n-1/2} $
    if num.prediction == "PmI" || num.prediction == "PmII" || num.prediction == "PmIIimposedpressure" || num.prediction == "PmIIimposedpressureBCincrement"
        #cf Brown 2001

        ∇ϕ_x = opC_u.AxT * opC_u.Rx * vec1(pD,grid) .+ opC_u.Gx_b * vecb(pD,grid)
        ∇ϕ_y = opC_v.AyT * opC_v.Ry * vec1(pD,grid) .+ opC_v.Gy_b * vecb(pD,grid)
        for iLS in 1:nLS
            ∇ϕ_x .+= opC_u.Gx[iLS] * veci(pD,grid,iLS+1)
            ∇ϕ_y .+= opC_v.Gy[iLS] * veci(pD,grid,iLS+1)
        end

        ph.Gxm1 .= 0.0 #TODO
        ph.Gym1 .= 0.0
        
        # gradient \times volume of cell (cells at border: volume = dx*dy/2)
        ph.Gxm1 .= copy(∇ϕ_x)
        ph.Gym1 .= copy(∇ϕ_y)

        #TODO add gradient check
    
        # ph.Gxm1 .= ∇ϕ_x
        # ph.Gym1 .= ∇ϕ_y


        PDI_status = @ccall "libpdi".PDI_multi_expose("print_pressure_gradient_in_prediction"::Cstring,
        "grad_x_1D"::Cstring, ph.Gxm1::Ptr{Cdouble}, PDI_OUT::Cint,
        "grad_y_1D"::Cstring, ph.Gym1::Ptr{Cdouble}, PDI_OUT::Cint,
        "p_1D"::Cstring, ph.pD::Ptr{Cdouble}, PDI_OUT::Cint,
        C_NULL::Ptr{Cvoid})::Cint

        ∇ϕ_x .= 0.0
        ∇ϕ_y .= 0.0

        grad_x = zeros(grid_u)
        grad_y = zeros(grid_v)
        compute_grad_T_x_T_y_array_u_v_capacities!(num, grid, grid_u, grid_v, opC_u, opC_v, grad_x, grad_y, ph.pD)
    
        #TODO divergence level
      
        PDI_status = @ccall "libpdi".PDI_multi_expose("check_pressure_velocity_end"::Cstring,
        # "grad_x"::Cstring,grad_x::Ptr{Cdouble}, PDI_OUT::Cint,
        # "grad_y"::Cstring, grad_y::Ptr{Cdouble}, PDI_OUT::Cint,
        "grad_u"::Cstring,grad_x::Ptr{Cdouble}, PDI_OUT::Cint,
        "grad_v"::Cstring, grad_y::Ptr{Cdouble}, PDI_OUT::Cint,
        "u_1D"::Cstring, ucorrD::Ptr{Cdouble}, PDI_OUT::Cint,
        "v_1D"::Cstring, vcorrD::Ptr{Cdouble}, PDI_OUT::Cint,
        "p_1D"::Cstring, ph.pD::Ptr{Cdouble}, PDI_OUT::Cint,
        C_NULL::Ptr{Cvoid})::Cint
    
        
    end
    #endregion add gradient to prediction

    nip = grid.nx * grid.ny

    niu = grid_u.nx * grid_u.ny
    nbu = 2 * grid_u.nx + 2 * grid_u.ny
    ntu = (nLS - nNavier + 1) * niu + nbu

    niv = grid_v.nx * grid_v.ny
    nbv = 2 * grid_v.nx + 2 * grid_v.ny
    ntv = (nLS - nNavier + 1) * niv + nbv

    if num.prediction == "PmIIimposedpressure" || num.prediction == "PmIIimposedpressureBCincrement" 
        BC_Poisson = Boundaries() #Neumann everywhere
    else
        BC_Poisson = copy(BC_p) 
    end

    if is_Forward_Euler(time_scheme)
        rhs_u, rhs_v, rhs_ϕ, rhs_uv, Lp, bc_Lp, bc_Lp_b, Lu, bc_Lu, bc_Lu_b, Lv, bc_Lv, bc_Lv_b = set_Forward_Euler!(
            bc_int, num, grid, geo, grid_u, geo_u, grid_v, geo_v,
            opC_p, opC_u, opC_v, BC_Poisson, BC_u, BC_v,
            Au, Bu, Av, Bv, Aϕ, Auv, Buv,
            Lpm1, bc_Lpm1, bc_Lpm1_b, Lum1, bc_Lum1, bc_Lum1_b, Lvm1, bc_Lvm1, bc_Lvm1_b,
            Mum1, Mvm1, mu1_over_rho1, op_conv, ph,
            periodic_x, periodic_y, advection, ls_advection, navier
        )
    elseif is_Crank_Nicolson(time_scheme)
        rhs_u, rhs_v, rhs_ϕ, Lp, bc_Lp, bc_Lp_b, Lu, bc_Lu, bc_Lu_b, Lv, bc_Lv, bc_Lv_b = set_Crank_Nicolson!(
            bc_int, num, grid, geo, grid_u, geo_u, grid_v, geo_v,
            opC_p, opC_u, opC_v, BC_Poisson, BC_u, BC_v,
            Au, Bu, Av, Bv, Aϕ,
            Lpm1, bc_Lpm1, bc_Lpm1_b, Lum1, bc_Lum1, bc_Lum1_b, Lvm1, bc_Lvm1, bc_Lvm1_b,
            Mum1, Mvm1, mu1_over_rho1, op_conv, ph,
            periodic_x, periodic_y, advection, ls_advection
        )
    end

    ra_x = Ra .* sin(β) .* opC_u.M * vec(hcat(zeros(grid_u.ny), ph.T))
    ra_y = Ra .* cos(β) .* opC_v.M * vec(vcat(zeros(1,grid_v.nx), ph.T))

    grav_x = g .* sin(β) .* opC_u.M * fones(grid_u)
    grav_y = g .* cos(β) .* opC_v.M * fones(grid_v)

    Convu = fzeros(grid_u)
    Convv = fzeros(grid_v)
    Cui = Cu * vec(u) .+ CUTCu
    Cvi = Cv * vec(v) .+ CUTCv

    if advection
        # scheme
        if current_i == 1
            Convu .+= Cui
            Convv .+= Cvi
        else
            Convu .+= 1.5 .* Cui .- 0.5 .* Cum1 #Cui returned at the end of function to Cum1
            Convv .+= 1.5 .* Cvi .- 0.5 .* Cvm1
        end
    end

    
    # TODO PDI_multi_expose() #Cu u CUTCu

    # u and v are coupled if a Navier slip BC is employed inside, otherwise they are uncoupled
    if !navier
        # if is_wall_no_slip(bc_int)
        #     vec1(uD,grid_u) .= vec(u)
        #     # update_dirichlet_field!(grid_u, uD, u, BC_u)
        #     vec1(rhs_u,grid_u) .+= -τ .* (opC_u.AxT * opC_u.Rx * vec1(pD,grid) .+ opC_u.Gx_b * vecb(pD,grid))
        #     for iLS in 1:nLS
        #         vec1(rhs_u,grid_u) .+= -τ .* (opC_u.Gx[iLS] * veci(pD,grid,iLS+1))
        #     end
        # end
        mul!(rhs_u, Bu, uD, 1.0, 1.0)
        vec1(rhs_u,grid_u) .+= τ .* grav_x
        vec1(rhs_u,grid_u) .-= τ .* Convu
        vec1(rhs_u,grid_u) .+= τ .* ra_x
        vec1(rhs_u,grid_u) .-= τ .* irho1 .* ph.Gxm1 
        
        kill_dead_cells!(vec1(rhs_u,grid_u), grid_u, geo_u[end])
        for iLS in 1:nLS
            kill_dead_cells!(veci(rhs_u,grid_u,iLS+1), grid_u, geo_u[end])
        end
        # @time bicgstabl!(ucorrD, Au, rhs_u, log=true)
        try
            # @time bicgstabl!(ucorrD, Au, rhs_u, Pl=Diagonal(Au), log=true)
            @time ucorrD .= Au \ rhs_u
        catch e
            ucorrD .= Inf
            println(e)
        end

        kill_dead_cells!(vec1(ucorrD,grid_u), grid_u, geo_u[end])
        for iLS in 1:nLS
            kill_dead_cells!(veci(ucorrD,grid_u,iLS+1), grid_u, geo_u[end])
        end
        ucorr .= reshape(vec1(ucorrD,grid_u), grid_u)

        # if is_wall_no_slip(bc_int)
        #     vec1(vD,grid_v) .= vec(v)
        #     # update_dirichlet_field!(grid_v, vD, v, BC_v)
        #     vec1(rhs_v,grid_v) .+= -τ .* (opC_v.AyT * opC_v.Ry * vec1(pD,grid) .+opC_v.Gy_b * vecb(pD,grid))
        #     for iLS in 1:nLS
        #         vec1(rhs_v,grid_v) .+= -τ .* (opC_v.Gy[iLS] * veci(pD,grid,iLS+1))
        #     end
        # end

        PDI_status = @ccall "libpdi".PDI_multi_expose("rhs_v"::Cstring,
        "v_1D"::Cstring, rhs_v::Ptr{Cdouble}, PDI_OUT::Cint,
        C_NULL::Ptr{Cvoid})::Cint
        
        # add M * vD , i.e. dx*dy*v:
        mul!(rhs_v, Bv, vD, 1.0, 1.0) 

        PDI_status = @ccall "libpdi".PDI_multi_expose("rhs_v"::Cstring,
        "v_1D"::Cstring, rhs_v::Ptr{Cdouble}, PDI_OUT::Cint,
        C_NULL::Ptr{Cvoid})::Cint

        vec1(rhs_v,grid_v) .+= - τ .* grav_y #TODO - minus sign here + sign there

        PDI_status = @ccall "libpdi".PDI_multi_expose("rhs_v"::Cstring,
        "v_1D"::Cstring, rhs_v::Ptr{Cdouble}, PDI_OUT::Cint,
        C_NULL::Ptr{Cvoid})::Cint

        vec1(rhs_v,grid_v) .-= τ .* Convv
        
        PDI_status = @ccall "libpdi".PDI_multi_expose("rhs_v"::Cstring,
        "v_1D"::Cstring, rhs_v::Ptr{Cdouble}, PDI_OUT::Cint,
        C_NULL::Ptr{Cvoid})::Cint

        vec1(rhs_v,grid_v) .+= τ .* ra_y

        PDI_status = @ccall "libpdi".PDI_multi_expose("rhs_v"::Cstring,
        "v_1D"::Cstring, rhs_v::Ptr{Cdouble}, PDI_OUT::Cint,
        C_NULL::Ptr{Cvoid})::Cint

        vec1(rhs_v,grid_v) .-= τ .* irho1 .* ph.Gym1

        PDI_status = @ccall "libpdi".PDI_multi_expose("rhs_v"::Cstring,
        "v_1D"::Cstring, rhs_v::Ptr{Cdouble}, PDI_OUT::Cint,
        C_NULL::Ptr{Cvoid})::Cint

        kill_dead_cells!(vec1(rhs_v,grid_v), grid_v, geo_v[end])
        for iLS in 1:nLS
            kill_dead_cells!(veci(rhs_v,grid_v,iLS+1), grid_v, geo_v[end])
        end
        # bicgstabl!(vcorrD, Av, rhs_v, log=true)
        

        try
            # @time bicgstabl!(vcorrD, Av, rhs_v, Pl=Diagonal(Av), log=true)
            @time vcorrD .= Av \ rhs_v
        catch e
            vcorrD .= Inf
            println(e)
        end


        PDI_status = @ccall "libpdi".PDI_multi_expose("print_velocity_prediction"::Cstring,
        "u_1D"::Cstring, ucorrD::Ptr{Cdouble}, PDI_OUT::Cint,
        "v_1D"::Cstring, vcorrD::Ptr{Cdouble}, PDI_OUT::Cint,
        "p_1D"::Cstring, ph.pD::Ptr{Cdouble}, PDI_OUT::Cint,
        C_NULL::Ptr{Cvoid})::Cint

        # II = CartesianIndex(div(grid_v.ny,2),1)
        # pII = lexicographic(II,grid_v.ny)
        # print("\n A coeff ",Av[pII,:])

        # PDI_status = @ccall "libpdi".PDI_multi_expose("print_velocity_prediction"::Cstring,
        # "u_1D"::Cstring, ucorrD::Ptr{Cdouble}, PDI_OUT::Cint,
        # "v_1D"::Cstring, vcorrD::Ptr{Cdouble}, PDI_OUT::Cint,
        # "p_1D"::Cstring, ph.pD::Ptr{Cdouble}, PDI_OUT::Cint,
        # C_NULL::Ptr{Cvoid})::Cint

        # print("\n proj v ",reshape(vec1(vcorrD,grid_v),grid_v)[div(grid_v.ny,2),:]," size ",size(reshape(vec1(vcorrD,grid_v),grid_v)[div(grid_v.ny,2),:]))
        # print("\n proj v vecb_L",vecb_L(vcorrD,grid_v)," size ",size(vecb_L(vcorrD,grid_v))," size ",size(vecb_B(vcorrD,grid_v)))


        kill_dead_cells!(vec1(vcorrD,grid_v), grid_v, geo_v[end])
        for iLS in 1:nLS
            kill_dead_cells!(veci(vcorrD,grid_v,iLS+1), grid_v, geo_v[end])
        end
        vcorr .= reshape(vec1(vcorrD,grid_v), grid_v)
    else #navier
        uvm1 = zeros(ntu + ntv + nNavier * nip)
        uvm1[1:niu] .= vec1(uD,grid_u)
        uvm1[ntu+1:ntu+niv] .= vec1(vD,grid_v)
        uvm1[ntu-nbu+1:ntu] .= vecb(uD,grid_u)
        uvm1[ntu+ntv-nbv+1:ntu+ntv] .= vecb(vD,grid_v)
        _iLS = 1
        for iLS in 1:num.nLS
            if !is_navier(bc_int[iLS]) && !is_navier_cl(bc_int[iLS])
                uvm1[_iLS*niu+1:(_iLS+1)*niu] .= veci(uD,grid_u,iLS+1)
                uvm1[ntu+_iLS*niv+1:ntu+(_iLS+1)*niv] .= veci(vD,grid_v,iLS+1)
                _iLS += 1
            end
        end

        rhs_uv .+=  Buv * uvm1

        rhs_uv[1:niu] .+= τ .* grav_x
        rhs_uv[1:niu] .-= τ .* Convu
        rhs_uv[1:niu] .+= τ .* ra_x
        rhs_uv[1:niu] .-= τ .* irho1 .* ph.Gxm1 

        rhs_uv[ntu+1:ntu+niv] .+= τ .* grav_y
        rhs_uv[ntu+1:ntu+niv] .-= τ .* Convv
        rhs_uv[ntu+1:ntu+niv] .+= τ .* ra_y
        rhs_uv[ntu+1:ntu+niv] .-= τ .* irho1 .* ph.Gym1 

        @views kill_dead_cells!(rhs_uv[1:niu], grid_u, geo_u[end])
        @views kill_dead_cells!(rhs_uv[ntu+1:ntu+niv], grid_v, geo_v[end])
        _iLS = 1
        for iLS in 1:nLS
            sbu = _iLS*niu+1:(_iLS+1)*niu
            sbv = ntu+_iLS*niv+1:ntu+(_iLS+1)*niv
            if !is_navier(bc_int[iLS]) && !is_navier_cl(bc_int[iLS])
                @views kill_dead_cells!(rhs_uv[sbu], grid_u, geo_u[end])
                @views kill_dead_cells!(rhs_uv[sbv], grid_v, geo_v[end])
                _iLS += 1
            end
        end

        uvD = ones(ntu + ntv + nNavier * nip)
        try
            @time uvD .= Auv \ rhs_uv
        catch e
            uvD .= Inf
            println(e)
        end

        vec1(ucorrD, grid_u) .= uvD[1:niu]
        vecb(ucorrD, grid_u) .= uvD[ntu-nbu+1:ntu]
        kill_dead_cells!(vec1(ucorrD,grid_u), grid_u, geo_u[end])
        ucorr .= reshape(vec1(ucorrD,grid_u), grid_u)

        vec1(vcorrD, grid_v) .= uvD[ntu+1:ntu+niv]
        vecb(vcorrD, grid_v) .= uvD[ntu+ntv-nbv+1:ntu+ntv]
        kill_dead_cells!(vec1(vcorrD,grid_v), grid_v, geo_v[end])
        vcorr .= reshape(vec1(vcorrD,grid_v), grid_v)

        nNav = 0
        _iLS = 1
        for iLS in 1:nLS
            if !is_navier(bc_int[iLS]) && !is_navier_cl(bc_int[iLS])
                veci(ucorrD,grid_u,iLS+1) .= uvD[_iLS*niu+1:(_iLS+1)*niu]
                kill_dead_cells!(veci(ucorrD,grid_u,iLS+1), grid_u, geo_u[end])

                veci(vcorrD,grid_v,iLS+1) .= uvD[ntu+_iLS*niv+1:ntu+(_iLS+1)*niv]
                kill_dead_cells!(veci(vcorrD,grid_v,iLS+1), grid_v, geo_v[end])
                _iLS += 1
            else
                @inbounds uT[nNav+1,:] .= vec(uvD[ntu+ntv+1+nNav*nip:ntu+ntv+(nNav+1)*nip])
                nNav += 1
            end
        end
    end

    PDI_status = @ccall "libpdi".PDI_multi_expose("print_velocity_prediction"::Cstring,
    "u_1D"::Cstring, ucorrD::Ptr{Cdouble}, PDI_OUT::Cint,
    "v_1D"::Cstring, vcorrD::Ptr{Cdouble}, PDI_OUT::Cint,
    "p_1D"::Cstring, ph.pD::Ptr{Cdouble}, PDI_OUT::Cint,
    C_NULL::Ptr{Cvoid})::Cint

    II = CartesianIndex(div(grid_v.ny,2),1)
    pII = lexicographic(II,grid_v.ny)
    # print("\n A coeff ",Av[pII,:])
    # print("\nM ",opC_p.iMy.diag[pII])

    print("\nAv[pII,:]./mu1_over_rho1 ",Av[pII,:]./mu1_over_rho1)

    # print("\nmu1_over_rho1 ",mu1_over_rho1)


    # # Test analytical vel
    # ucorrD = copy(uD)
    # vcorrD = copy(vD)
    # print("\n test analytical vel ")

    # PDI_status = @ccall "libpdi".PDI_multi_expose("print_velocity_prediction"::Cstring,
    # "u_1D"::Cstring, ucorrD::Ptr{Cdouble}, PDI_OUT::Cint,
    # "v_1D"::Cstring, vcorrD::Ptr{Cdouble}, PDI_OUT::Cint,
    # "p_1D"::Cstring, ph.pD::Ptr{Cdouble}, PDI_OUT::Cint,
    # C_NULL::Ptr{Cvoid})::Cint


    #endregion prediction 

    #region correction 

    # Compute divergence of velocity
    Duv = opC_p.AxT * vec1(ucorrD,grid_u) .+ opC_p.Gx_b * vecb(ucorrD,grid_u) .+
          opC_p.AyT * vec1(vcorrD,grid_v) .+ opC_p.Gy_b * vecb(vcorrD,grid_v)
    for iLS in 1:nLS
        if !is_navier(bc_int[iLS]) && !is_navier_cl(bc_int[iLS])
            Duv .+= opC_p.Gx[iLS] * veci(ucorrD,grid_u,iLS+1) .+ 
                    opC_p.Gy[iLS] * veci(vcorrD,grid_v,iLS+1)
        end
    end


    #TODO replace

    # Poisson equation: source term
    # divergence of velocity / dt
    vec1(rhs_ϕ,grid) .= iτ .* Duv

    #region needs to be corrected/documented for the signs
    
    # pres_free_suface = 0.0
    #TODO Marangoni
    #TODO phase change
    diff_inv_rho = 1.0/rho1 - 1.0/rho2
    # jump_mass_flux = 0.0 #TODO

    if jump_mass_flux
        for iLS in 1:nLS
            if is_fs(bc_int[iLS])
                Smat = strain_rate(iLS, opC_u, opC_v, opC_p)
                S = Smat[1,1] * vec1(ucorrD,grid_u) .+ Smat[1,2] * veci(ucorrD,grid_u,iLS+1) .+
                    Smat[2,1] * vec1(vcorrD,grid_v) .+ Smat[2,2] * veci(vcorrD,grid_v,iLS+1)
    
                fs_mat = opC_p.HxT[iLS] * opC_p.Hx[iLS] .+ opC_p.HyT[iLS] * opC_p.Hy[iLS]
                veci(rhs_ϕ,grid,iLS+1) .= -2.0 .* mu1_over_rho1 .* S .+ Diagonal(diag(fs_mat)) * ( σ .* vec(grid.LS[iLS].κ) .- pres_free_suface .- diff_inv_rho * mass_flux ^ 2)
            end
        end
    else
        for iLS in 1:nLS
            if is_fs(bc_int[iLS])
                Smat = strain_rate(iLS, opC_u, opC_v, opC_p)
                S = Smat[1,1] * vec1(ucorrD,grid_u) .+ Smat[1,2] * veci(ucorrD,grid_u,iLS+1) .+
                    Smat[2,1] * vec1(vcorrD,grid_v) .+ Smat[2,2] * veci(vcorrD,grid_v,iLS+1)

                fs_mat = opC_p.HxT[iLS] * opC_p.Hx[iLS] .+ opC_p.HyT[iLS] * opC_p.Hy[iLS]
                veci(rhs_ϕ,grid,iLS+1) .= -2.0 .* mu1_over_rho1 .* S .+ Diagonal(diag(fs_mat)) * ( σ .* vec(grid.LS[iLS].κ) .- pres_free_suface )
            end
        end
    end
    # Remove nullspace by adding small quantity to main diagonal
    if num.null_space == 0
        @inbounds @threads for i in 1:Aϕ.m
            @inbounds Aϕ[i,i] += 1e-10
        end
    end
    kill_dead_cells!(vec1(rhs_ϕ,grid), grid, geo[end])
    for iLS in 1:nLS
        kill_dead_cells!(veci(rhs_ϕ,grid,iLS+1), grid, geo[end])
    end
    # @time bicgstabl!(ϕD, Aϕ, rhs_ϕ, Pl = Diagonal(Aϕ), log = true)

    # rhs_ϕ .*= rho1 #TODO #TODO not BC

    # vec1(rhs_ϕ,grid) .*= rho1 

    # Aϕ .*= irho1

    # vecb(rhs_ϕ,grid) .*= irho1

    #endregion needs to be corrected/documented for the signs

    # Solve Poisson equation

    # \phi^{n+1}: ϕD
    @time ϕD .= Aϕ \ rhs_ϕ
    kill_dead_cells!(vec1(ϕD,grid), grid, geo[end])
    for iLS in 1:nLS
        kill_dead_cells!(veci(ϕD,grid,iLS+1), grid, geo[end])
    end
    ϕ .= reshape(vec1(ϕD,grid), grid)

    iMu = Diagonal(inv_weight_eps2.(num.epsilon_mode,num.epsilon_vol,opC_u.M.diag))
    iMv = Diagonal(inv_weight_eps2.(num.epsilon_mode,num.epsilon_vol,opC_v.M.diag))
    # Gradient of pressure, eq. 17 in 
    #"A Conservative Cartesian Cut-Cell Method for Mixed Boundary Conditions and the Incompressible Navier-Stokes Equations on Staggered Meshes"
    ∇ϕ_x = opC_u.AxT * opC_u.Rx * vec(ϕ) .+ opC_u.Gx_b * vecb(ϕD,grid)
    ∇ϕ_y = opC_v.AyT * opC_v.Ry * vec(ϕ) .+ opC_v.Gy_b * vecb(ϕD,grid)
    for iLS in 1:nLS
        ∇ϕ_x .+= opC_u.Gx[iLS] * veci(ϕD,grid,iLS+1)
        ∇ϕ_y .+= opC_v.Gy[iLS] * veci(ϕD,grid,iLS+1)
    end

    # ∇ϕ_x = irho1 .* opC_u.AxT * opC_u.Rx * vec(ϕ) .+ opC_u.Gx_b * vecb(ϕD,grid)
    # ∇ϕ_y = irho1 .* opC_v.AyT * opC_v.Ry * vec(ϕ) .+ opC_v.Gy_b * vecb(ϕD,grid)
    # for iLS in 1:nLS
    #     ∇ϕ_x .+= irho1 .* opC_u.Gx[iLS] * veci(ϕD,grid,iLS+1)
    #     ∇ϕ_y .+= irho1 .* opC_v.Gy[iLS] * veci(ϕD,grid,iLS+1)
    # end

    # if num.prediction == 1 already done
    #     ph.Gxm1 .+= ∇ϕ_x
    #     ph.Gym1 .+= ∇ϕ_y
    # end


    # iM = Diagonal(1. ./ (vec(geo[end].dcap[:,:,5]) .+ eps(0.01)))

    # iM = Diagonal(inv_weight_eps.(num,geo[end].dcap[:,:,5]))

    iM = Diagonal(inv_weight_eps2.(num.epsilon_mode,num.epsilon_vol,vec(geo[end].dcap[:,:,5])))

    # iM = Diagonal(1. ./ (vec(geo[end].dcap[:,:,5]) ))

    # if is_fs(bc_int)
    # p^{n-1/2} is the pressure of the previous timestep (for pressure, which is lagging by dt/2 wrt other wariables)
    if num.prediction == "PmI"
        # \nabla_h p^{n+1/2} = \nabla_h p^{n-1/2} + \nabla_h \phi^{n+1}.
        # "not consistent with a second-order discretization of the Navier–Stokes equations since, 
        # due to Eq. (72), the normal component of the pressure gradient will remain constant in time at the boundary"
        # Brown 2001
        vec1(pD,grid) .+= vec(ϕ) #no τ  since div u not rho1

    elseif num.prediction == "PmII" || num.prediction == "PmIIimposedpressure" || num.prediction == "PmIIimposedpressureBCincrement"
        # \nabla_h p^{n+1/2} = \nabla_h p^{n-1/2} + \nabla_h \phi^{n+1} - 
        #nu dt/2
        # \Delta t \nabla_h^2 \phi^{n+1} = \nabla_h \cdot \mathbf{u}^{*} \quad \text{in } \Omega
        

        print("\n max p vec1 ",maximum(vec1(pD,grid)))
        print("\n max p",minimum(ϕ),maximum(ϕ))
        print("\n max p",minimum(num.mu_cin1./2 .* reshape(iM * Duv,grid))," ",maximum(num.mu_cin1./2 .* reshape(iM * Duv,grid)))

        # todo zero neumann bc !!!!
        # vec1(pD,grid) .+= vec(ϕ .- num.mu_cin1./2 .* reshape(iM * Duv,grid))
        # or,
        # for better readability
        vec1(pD,grid) .= vec1(pD,grid) .+ vec(ϕ .- num.mu_cin1./2 .* reshape(iM * Duv,grid)) 
        # vec1(pD,grid) .+= vec(ϕ .- num.mu_cin1./2 .* reshape(iM * Duv,grid)) 

        print("\n max p vec1 ",maximum(vec1(pD,grid)))
        print("\n max p",maximum(ϕ))
        print("\n max p",minimum(num.mu_cin1./2 .* reshape(iM * Duv,grid))," ",maximum(num.mu_cin1./2 .* reshape(iM * Duv,grid)))



    elseif num.prediction == "PmIII"
        # \Delta t \nabla_h^2 \phi^{n+1} = \nabla_h \cdot \mathbf{u}^{*} \quad \text{in } \Omega
        vec1(pD,grid) .= vec(ϕ .- num.mu_cin1./2 .* reshape(iM * Duv,grid)) #no contribution from p^{n-1/2}
    
    elseif num.prediction == "Flower" #occursin("Flower",num.prediction)
        # \nabla_h p^{n+1/2} = \nabla_h \phi^{n+1} # TODO: does not correspond to any formula in Brown 2001 ?
        vec1(pD,grid) .= vec(ϕ) #.- mu1_over_rho1 .* reshape(iM * Duv, grid))
    
    else
        @error("wrong prediction method, does not exist")
    end

    # PDI_status = @ccall "libpdi".PDI_multi_expose("print_pressure_projection"::Cstring,
    # "u_1D"::Cstring, ucorrD::Ptr{Cdouble}, PDI_OUT::Cint,
    # "v_1D"::Cstring, vcorrD::Ptr{Cdouble}, PDI_OUT::Cint,
    # "p_1D"::Cstring, ph.pD::Ptr{Cdouble}, PDI_OUT::Cint,
    # C_NULL::Ptr{Cvoid})::Cint

    #region interfacial pressure is overwritten
    # TODO check and document: for pressure-imposed Poiseuille, not sure
    #  
    if num.prediction == "PmIIimposedpressure"

    elseif num.prediction == "PmIIimposedpressureBCincrement"
        #TODO reapply BC 
        #TODO
        #init p for Poiseuille : grad_x = 0
        #increment : grad_x=0
        # so grad_x still zero ? even with div ?
        tmp_vec_p = zeros(grid) 
        # tmp_vec_p .= 0.0
        get_height!(grid.LS[1],grid.ind,grid.dx,grid.dy,grid.LS[end].geoS,tmp_vec_p) #here tmp_vec_p solid

        init_fields_multiple_levelsets!(num,ph.pD,ph.p,tmp_vec_p,BC_p,grid,num.pres_intfc,"pL")

    else #update pressure at boundaries
        for iLS in 1:nLS
            veci(pD,grid,iLS+1) .= veci(ϕD,grid,iLS+1)
        end
        vecb(pD,grid) .= vecb(ϕD,grid) 
    end


    #endregion interfacial pressure is overwritten

    #TODO reapply Neumann BC for pressure boundary values ?


    p .= reshape(vec1(pD,grid), grid)

    #TODO
    # compute_grad_p!(num,grid, grid_u, grid_v, pD, opC_p, opC_u, opC_v)


    # else
    #     vec1(pD,grid) .= vec(p) .+ vec(ϕ) #.- mu1_over_rho1 .* iM * Duv
    #     vec2(pD,grid) .+= vec2(ϕD,grid)
    #     vecb(pD,grid) .+= vecb(ϕD,grid)
    #     p .= reshape(vec1(pD,grid), grid)
    # end

    # vec1(∇ϕ_x,grid) .*= irho1 
    # vec1(∇ϕ_y,grid) .*= irho1

    
    # u .= ucorr .- τ .* reshape(iMu * ∇ϕ_x, grid_u)
    # v .= vcorr .- τ .* reshape(iMv * ∇ϕ_y, grid_v)

    u .= ucorr .- τ .* irho1 .* reshape(iMu * ∇ϕ_x, grid_u)
    v .= vcorr .- τ .* irho1 .* reshape(iMv * ∇ϕ_y, grid_v)

    kill_dead_cells!(u, grid_u, geo_u[end])
    kill_dead_cells!(v, grid_v, geo_v[end])

    vec1(uD,grid_u) .= vec(u)
    vecb(uD,grid_u) .= vecb(ucorrD,grid_u)
    vec1(vD,grid_v) .= vec(v)
    vecb(vD,grid_v) .= vecb(vcorrD,grid_v)
    for iLS in 1:nLS
        if !is_navier(bc_int[iLS]) && !is_navier_cl(bc_int[iLS])
            veci(uD,grid_u,iLS+1) .= veci(ucorrD,grid_u,iLS+1)
            veci(vD,grid_v,iLS+1) .= veci(vcorrD,grid_v,iLS+1)
        end
        # if is_fs(bc_int[iLS])
        #     @inbounds for II in grid_u.ind.all_indices
        #         pII = lexicographic(II, grid_u.ny)
        #         if abs(veci(ucorrD,grid_u,iLS+1)[pII]) > 1e-12
        #             veci(ucorrD,grid_u,iLS+1)[pII] -= (τ .* iMu * ∇ϕ_x)[pII]
        #         end
        #     end
        #     @inbounds for II in grid_v.ind.all_indices
        #         pII = lexicographic(II, grid_v.ny)
        #         if abs(veci(vcorrD,grid_v,iLS+1)[pII]) > 1e-12
        #             veci(vcorrD,grid_v,iLS+1)[pII] -= (τ .* iMv * ∇ϕ_y)[pII]
        #         end
        #     end
        # end
    end


    grad_x = zeros(grid_u)
    grad_y = zeros(grid_v)
    compute_grad_T_x_T_y_array_u_v_capacities!(num, grid, grid_u, grid_v, opC_u, opC_v, grad_x, grad_y, ph.pD)

    #TODO divergence level
  
    PDI_status = @ccall "libpdi".PDI_multi_expose("check_pressure_velocity_end"::Cstring,
    # "grad_x"::Cstring,grad_x::Ptr{Cdouble}, PDI_OUT::Cint,
    # "grad_y"::Cstring, grad_y::Ptr{Cdouble}, PDI_OUT::Cint,
    "grad_u"::Cstring,grad_x::Ptr{Cdouble}, PDI_OUT::Cint,
    "grad_v"::Cstring, grad_y::Ptr{Cdouble}, PDI_OUT::Cint,
    "u_1D"::Cstring, ucorrD::Ptr{Cdouble}, PDI_OUT::Cint,
    "v_1D"::Cstring, vcorrD::Ptr{Cdouble}, PDI_OUT::Cint,
    "p_1D"::Cstring, ph.pD::Ptr{Cdouble}, PDI_OUT::Cint,
    C_NULL::Ptr{Cvoid})::Cint

    # PDI_status = @ccall "libpdi".PDI_multi_expose("print_velocity_correction"::Cstring,
    # "u_1D"::Cstring, ucorrD::Ptr{Cdouble}, PDI_OUT::Cint,
    # "v_1D"::Cstring, vcorrD::Ptr{Cdouble}, PDI_OUT::Cint,
    # "p_1D"::Cstring, ph.pD::Ptr{Cdouble}, PDI_OUT::Cint,
    # C_NULL::Ptr{Cvoid})::Cint

    #endregion correction 


    return Lp, bc_Lp, bc_Lp_b, Lu, bc_Lu, bc_Lu_b, Lv, bc_Lv, bc_Lv_b, opC_p.M, opC_u.M, opC_v.M, Cui, Cvi
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
    # i = 0
    # while sum(res) > 1e-8
    # for i in 1:1
    #     ϵ = 1e-10
        
    #     res .= residual(u_guess, v_guess, num, grid, geo, grid_u, geo_u, grid_v, geo_v, u, v, op_conv, ph, BC_u, BC_v)
    #     dres .= dresidual(u_guess, v_guess, res, ϵ, num, grid, geo, grid_u, geo_u, grid_v, geo_v, u, v, op_conv, ph, BC_u, BC_v)

    #     u_guess .= u_guess .- reshape((res ./ dres)[1:grid_u.nx*grid_u.ny], grid_u)
    #     v_guess .= v_guess .- reshape((res ./ dres)[grid_u.nx*grid_u.ny+1:end], grid_v)

    #     i += 1
    #     println("sum: $(sum(res)) | max: $(maximum(res))")
    # end

    # u_midp = 0.5 .* (u .+ u_guess)
    # v_midp = 0.5 .* (v .+ v_guess)
    u_midp = copy(u)
    v_midp = copy(v)
    set_convection!(num, grid, geo, grid_u, grid_u.LS, grid_v, grid_v.LS, u_midp, v_midp, op_conv, ph, BC_u, BC_v,opC_p, opC_u, opC_v)

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

    set_convection!(num, grid, geo, grid_u, grid_u.LS, grid_v, grid_v.LS, u_midp, v_midp, op_conv, ph, BC_u, BC_v,opC_p, opC_u, opC_v)
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