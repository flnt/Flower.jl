function fill_empty_rows!(num, grid, geo, O)
    @unpack τ = num
    @unpack nx, ny, dx, dy, ind = grid

    @inbounds @threads for II in ind.all_indices
        pII = lexicographic(II, ny)
        if geo.cap[II,5] < 1e-12
            O[pII,pII] = dx[II] * dy[II] + 0.5 * τ * 4.0
        end
        if (2*dx[II]+2*dy[II] - sum(geo.dcap[II,1:4])) < 1e-12
            O[pII+nx*ny,pII+nx*ny] = 1.
        end
    end

    return nothing
end


"""
pad(A, a=1.0) adds a to the diagonal of A where it is null
"""
function pad(A, a=1.0)
    d = collect(diag(A))
    for i in eachindex(d)
        d[i] = ifelse(iszero(d[i]), a*one(d[i]), zero(d[i]))
    end
    A + Diagonal(d)
end


"""
pad_crank_nicolson

adds dx[II] * dy[II] + 0.5 * τ * 4.0 to the diagonal where it is null, for ind.all_indices
"""
function pad_crank_nicolson(A, grid, τ)
    @unpack ny, dx, dy, ind = grid

    d = collect(diag(A))
    for II in ind.all_indices
        pII = lexicographic(II, ny)
        pad = dx[II] * dy[II] + 0.5 * τ * 4.0
        d[pII] = ifelse(iszero(d[pII]), pad, zero(d[pII]))
    end
    A + Diagonal(d)
end

# @inline function get_capacities(cap, II)
#     @inbounds ret = (cap[II,1], cap[II,2], cap[II,3], cap[II,4], cap[II,6], cap[II,7],
#            cap[II,8] + eps(0.01), cap[II,9] + eps(0.01), cap[II,10] + eps(0.01), cap[II,11] + eps(0.01))
#     return ret
# end

# @inline function get_capacities_x(cap, II)
#     @inbounds ret = cap[II,1], cap[II,3], cap[II,6]
#     return ret
# end

# @inline function get_capacities_y(cap, II)
#     @inbounds ret = cap[II,2], cap[II,4], cap[II,7]
#     return ret
# end

"""
```julia
    return cap[II,1], cap[II,2], cap[II,3], cap[II,4], cap[II,6], cap[II,7]
```
"""
@inline function surface_capacities(cap, II)
    @inbounds return cap[II,1], cap[II,2], cap[II,3], cap[II,4], cap[II,6], cap[II,7]
end


@doc raw"""
    computes the capacity for the divergence, in the article, this corresponds to:

```math
-G ^ \top = \left [ \begin{array}{>{\centering\arraybackslash$} p{1.2cm} <{$} >{\centering\arraybackslash$} p{1.2cm} <{$}}
B _ x D _ x ^ + & B _ y D _ y ^ +
\end{array} \right ]
```

fills along x -dcap[II,6] +dcap[II,6] , along y -dcap[II,7] +dcap[II,7]
"""
function divergence_B!(Ox, Oy, dcap, n, all_indices)
    @inbounds @threads for II in all_indices
        pII_1 = lexicographic(II, n)
        pII_2 = lexicographic(δx⁺(II), n)
        
        B1 = dcap[II,6]
        B2 = dcap[II,7]

        pJJ_1 = lexicographic(II, n+1)
        pJJ_2 = lexicographic(δy⁺(II), n+1)
        
        @inbounds Ox[pII_1,pII_1] = -B1
        @inbounds Ox[pII_1,pII_2] = B1

        @inbounds Oy[pII_1,pJJ_1] = -B2
        @inbounds Oy[pII_1,pJJ_2] = B2
    end

    return nothing
end


"""
fills with  dcap 1, 2, 3, 4 i.e. A1 = dcap[II,1]
"""
function divergence_A!(grid, Ox, Oy, dcap, n, all_indices, per_x, per_y)
    @inbounds @threads for II in all_indices
        pII_1 = lexicographic(II, n)
        pII_2 = lexicographic(δx⁺(II), n)
        
        A1 = dcap[II,1]
        A2 = dcap[II,2]
        A3 = dcap[II,3]
        A4 = dcap[II,4]

        pJJ_1 = lexicographic(II, n+1)
        pJJ_2 = lexicographic(δy⁺(II), n+1)
        
        @inbounds Ox[pII_1,pII_1] = -A1
        @inbounds Ox[pII_1,pII_2] = A3

        @inbounds Oy[pII_1,pJJ_1] = -A2
        @inbounds Oy[pII_1,pJJ_2] = A4
    end
    if !per_x
        @inbounds @threads for II in all_indices[:,1]
            pII_1 = lexicographic(II, n)
            
            @inbounds Ox[pII_1,pII_1] = 0.0
        end
        @inbounds @threads for II in all_indices[:,end]
            pII_1 = lexicographic(II, n)
            pII_2 = lexicographic(δx⁺(II), n)
            
            @inbounds Ox[pII_1,pII_2] = 0.0
        end
    elseif isFCx(grid)
        @inbounds @threads for II in all_indices[:,1]
            pII_1 = lexicographic(II, n)
            pII_2 = lexicographic(δx⁺(II), n)
            pJJ_1 = lexicographic(II + CartesianIndex(0,grid.nx-1), n)
            pJJ_2 = lexicographic(II + CartesianIndex(0,grid.nx), n)
            
            @inbounds Ox[pII_1,pII_1] = Ox[pJJ_1,pJJ_1]
            @inbounds Ox[pJJ_1,pJJ_2] = Ox[pII_1,pII_2]
        end
    end
    if !per_y
        @inbounds @threads for II in all_indices[1,:]
            pII_1 = lexicographic(II, n)
            pJJ_1 = lexicographic(II, n+1)
            
            @inbounds Oy[pII_1,pJJ_1] = 0.0
        end
        @inbounds @threads for II in all_indices[end,:]
            pII_1 = lexicographic(II, n)
            pJJ_2 = lexicographic(δy⁺(II), n+1)
            
            @inbounds Oy[pII_1,pJJ_2] = 0.0
        end
    elseif isFCy(grid)
        @inbounds @threads for II in all_indices[1,:]
            pII_1 = lexicographic(II, n)
            pII_2 = lexicographic(δy⁺(II), n)
            pJJ_1 = lexicographic(II + CartesianIndex(grid.ny-1,0), n)
            pJJ_2 = lexicographic(II + CartesianIndex(grid.ny,0), n)
            
            @inbounds Oy[pII_1,pII_1] = Oy[pJJ_1,pJJ_1]
            @inbounds Oy[pJJ_1,pJJ_2] = Oy[pII_1,pII_2]
        end
    end

    return nothing
end

function bc_matrix!(grid::Mesh{GridCC,T,N}, Hx, Hy, dcap, dcap_p, n, all_indices) where {T,N}
    @inbounds @threads for II in @view all_indices[2:end,2:end]
        pII = lexicographic(II, n)
        pJJ = lexicographic(II, n+1)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, II)
        
        @inbounds Hx[pII, pII-n] = -(dcap[δx⁻(II),3] - dcap[δx⁻(II),6])
        @inbounds Hx[pII, pII] = -(B1 - A1)
        @inbounds Hy[pJJ, pII-1] = -(dcap[δy⁻(II),4] - dcap[δy⁻(II),7])
        @inbounds Hy[pJJ, pII] = -(B2 - A2)
    end

    @inbounds for II in @view all_indices[:,1]
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, II)
        
        @inbounds Hx[pII, pII] = -(B1 - A1)
    end
    @inbounds for II in @view all_indices[:,end]
        pII = lexicographic(δx⁺(II), n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, II)
        
        @inbounds Hx[pII, pII-n] = -(A3 - B1)
    end
    @inbounds for II in @views vcat(all_indices[1,2:end], all_indices[end,2:end])
        pII = lexicographic(II, n)
        pJJ = lexicographic(II, n+1)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, II)
        
        @inbounds Hx[pII, pII-n] = -(dcap[δx⁻(II),3] - dcap[δx⁻(II),6])
        @inbounds Hx[pII, pII] = -(B1 - A1)
    end

    @inbounds @threads for II in @view all_indices[1,:]
        pII = lexicographic(II, n)
        pJJ = lexicographic(II, n+1)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, II)
        
        @inbounds Hy[pJJ, pII] = -(B2 - A2)
    end
    @inbounds @threads for II in @view all_indices[end,:]
        pII = lexicographic(II, n)
        pJJ = lexicographic(δy⁺(II), n+1)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, II)
        
        @inbounds Hy[pJJ, pII] = -(A4 - B2)
    end
    @inbounds @threads for II in @views vcat(all_indices[2:end,1], all_indices[2:end,end])
        pII = lexicographic(II, n)
        pJJ = lexicographic(II, n+1)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, II)
        
        @inbounds Hy[pJJ, pII-1] = -(dcap[δy⁻(II),4] - dcap[δy⁻(II),7])
        @inbounds Hy[pJJ, pII] = -(B2 - A2)
    end

    return nothing
end

function bc_matrix!(grid::Mesh{GridFCx,T,N}, Hx, Hy, dcap, dcap_p, n, all_indices) where {T,N}
    @inbounds @threads for II in @view all_indices[2:end,2:end]
        pII = lexicographic(II, n)
        pJJ = lexicographic(II, n+1)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, II)
        
        @inbounds Hx[pII, pII-n] = -(dcap[δx⁻(II),3] - dcap[δx⁻(II),6])
        @inbounds Hx[pII, pII] = -(B1 - A1)
        @inbounds Hy[pJJ, pII-1] = -(dcap[δy⁻(II),4] - dcap[δy⁻(II),7])
        @inbounds Hy[pJJ, pII] = -(B2 - A2)
    end
    @inbounds for II in @view all_indices[:,1]
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, II)
        
        @inbounds Hx[pII, pII] = -(B1 - dcap_p[II,1])
    end
    @inbounds for II in @view all_indices[:,end]
        pII = lexicographic(δx⁺(II), n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, II)
        
        @inbounds Hx[pII, pII-n] = -(dcap_p[δx⁻(II),3] - B1)
    end
    @inbounds for II in @views vcat(all_indices[1,2:end], all_indices[end,2:end])
        pII = lexicographic(II, n)
        pJJ = lexicographic(II, n+1)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, II)
        
        @inbounds Hx[pII, pII-n] = -(dcap[δx⁻(II),3] - dcap[δx⁻(II),6])
        @inbounds Hx[pII, pII] = -(B1 - A1)
    end

    @inbounds @threads for II in @view all_indices[1,:]
        pII = lexicographic(II, n)
        pJJ = lexicographic(II, n+1)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, II)
        
        @inbounds Hy[pJJ, pII] = -(B2 - A2)
    end
    @inbounds @threads for II in @view all_indices[end,:]
        pII = lexicographic(II, n)
        pJJ = lexicographic(δy⁺(II), n+1)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, II)
        
        @inbounds Hy[pJJ, pII] = -(A4 - B2)
    end
    @inbounds @threads for II in @views vcat(all_indices[2:end,1], all_indices[2:end,end])
        pII = lexicographic(II, n)
        pJJ = lexicographic(II, n+1)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, II)
        
        @inbounds Hy[pJJ, pII-1] = -(dcap[δy⁻(II),4] - dcap[δy⁻(II),7])
        @inbounds Hy[pJJ, pII] = -(B2 - A2)
    end

    return nothing
end

function bc_matrix!(grid::Mesh{GridFCy,T,N}, Hx, Hy, dcap, dcap_p, n, all_indices) where {T,N}
    @inbounds @threads for II in @view all_indices[2:end,2:end]
        pII = lexicographic(II, n)
        pJJ = lexicographic(II, n+1)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, II)
        
        @inbounds Hx[pII, pII-n] = -(dcap[δx⁻(II),3] - dcap[δx⁻(II),6])
        @inbounds Hx[pII, pII] = -(B1 - A1)
        @inbounds Hy[pJJ, pII-1] = -(dcap[δy⁻(II),4] - dcap[δy⁻(II),7])
        @inbounds Hy[pJJ, pII] = -(B2 - A2)
    end

    @inbounds for II in @view all_indices[:,1]
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, II)
        
        @inbounds Hx[pII, pII] = -(B1 - A1)
    end
    @inbounds for II in @view all_indices[:,end]
        pII = lexicographic(δx⁺(II), n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, II)
        
        @inbounds Hx[pII, pII-n] = -(A3 - B1)
    end
    @inbounds for II in @views vcat(all_indices[1,2:end], all_indices[end,2:end])
        pII = lexicographic(II, n)
        pJJ = lexicographic(II, n+1)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, II)
        
        @inbounds Hx[pII, pII-n] = -(dcap[δx⁻(II),3] - dcap[δx⁻(II),6])
        @inbounds Hx[pII, pII] = -(B1 - A1)
    end

    @inbounds @threads for II in @view all_indices[1,:]
        pII = lexicographic(II, n)
        pJJ = lexicographic(II, n+1)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, II)
        
        @inbounds Hy[pJJ, pII] = -(B2 - dcap_p[II,2])
    end
    @inbounds @threads for II in @view all_indices[end,:]
        pII = lexicographic(II, n)
        pJJ = lexicographic(δy⁺(II), n+1)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, II)
        
        @inbounds Hy[pJJ, pII] = -(dcap_p[δy⁻(II),4] - B2)
    end
    @inbounds @threads for II in @views vcat(all_indices[2:end,1], all_indices[2:end,end])
        pII = lexicographic(II, n)
        pJJ = lexicographic(II, n+1)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, II)
        
        @inbounds Hy[pJJ, pII-1] = -(dcap[δy⁻(II),4] - dcap[δy⁻(II),7])
        @inbounds Hy[pJJ, pII] = -(B2 - A2)
    end

    return nothing
end

function bc_matrix!(grid, Hx, Hy, dcap, dcap_u, dcap_v, n, all_indices)
    @inbounds @threads for II in @view all_indices[2:end,2:end]
        pII = lexicographic(II, n)
        pJJ = lexicographic(II, n+1)
        @inbounds A1, A2 = dcap[II,1], dcap[II,2]
        @inbounds Ax1, Ax3 = dcap_u[II,1], dcap_u[II,3]
        @inbounds Ay2, Ay4 = dcap_v[II,2], dcap_v[II,4]
        
        @inbounds Hx[pII, pII-n] = -(A1 - Ax1)
        @inbounds Hx[pII, pII] = -(Ax3 - A1)
        @inbounds Hy[pJJ, pII-1] = -(A2 - Ay2)
        @inbounds Hy[pJJ, pII] = -(Ay4 - A2)
    end

    @inbounds for II in @view all_indices[:,1]
        pII = lexicographic(II, n)
        @inbounds A1 = dcap[II,1]
        @inbounds Ax3 = dcap_u[II,3]
        
        @inbounds Hx[pII, pII] = -(Ax3 - A1)
    end
    @inbounds for II in @view all_indices[:,end]
        pII = lexicographic(δx⁺(II), n)
        @inbounds A3 = dcap[II,3]
        @inbounds Ax3 = dcap_u[II,3]
        
        @inbounds Hx[pII, pII-n] = -(A3 - Ax3)
    end
    @inbounds for II in @views vcat(all_indices[1,2:end], all_indices[end,2:end])
        pII = lexicographic(II, n)
        @inbounds A1 = dcap[II,1]
        @inbounds Ax1, Ax3 = dcap_u[II,1], dcap_u[II,3]
        
        @inbounds Hx[pII, pII-n] = -(A1 - Ax1)
        @inbounds Hx[pII, pII] = -(Ax3 - A1)
    end

    @inbounds @threads for II in @view all_indices[1,:]
        pII = lexicographic(II, n)
        pJJ = lexicographic(II, n+1)
        @inbounds A2 = dcap[II,2]
        @inbounds Ay4 = dcap_v[II,4]
        
        @inbounds Hy[pJJ, pII] = -(Ay4 - A2)
    end
    @inbounds @threads for II in @view all_indices[end,:]
        pII = lexicographic(II, n)
        pJJ = lexicographic(δy⁺(II), n+1)
        @inbounds A4 = dcap[II,4]
        @inbounds Ay4 = dcap_v[II,4]
        
        @inbounds Hy[pJJ, pII] = -(A4 - Ay4)
    end
    @inbounds @threads for II in @views vcat(all_indices[2:end,1], all_indices[2:end,end])
        pII = lexicographic(II, n)
        pJJ = lexicographic(II, n+1)
        @inbounds A2 = dcap[II,2]
        @inbounds Ay2, Ay4 = dcap_v[II,2], dcap_v[II,4]
        
        @inbounds Hy[pJJ, pII-1] = -(A2 - Ay2)
        @inbounds Hy[pJJ, pII] = -(Ay4 - A2)
    end

    return nothing
end

function periodic_bcs!(grid, Gx, Gy, Hx, Hy, periodic_x, periodic_y)
    @unpack nx, ny, ind = grid

    if periodic_x
        @inbounds @threads for i in eachindex(grid.ind.periodic_x[1])
            II = grid.ind.periodic_x[1][i]
            JJ = grid.ind.periodic_x[2][i]
            pII = lexicographic(II, ny)
            pJJ = lexicographic(JJ, ny)

            Gx[pII,pJJ] = Gx[pJJ+ny,pJJ]
            Gx[pJJ+ny,pII] = Gx[pII,pII]

            Hx[pII,pJJ] = Hx[pJJ+ny,pJJ]
            Hx[pJJ+ny,pII] = Hx[pII,pII]
        end
    end
    if periodic_y
        @inbounds @threads for i in eachindex(grid.ind.periodic_y[1])
            II = grid.ind.periodic_y[1][i]
            JJ = grid.ind.periodic_y[2][i]
            pII = lexicographic(II, ny)
            pJJ = lexicographic(JJ, ny)
            pII_y = lexicographic(II, ny+1)
            pJJ_y = lexicographic(δy⁺(JJ), ny+1)

            Gy[pII_y,pJJ] = Gy[pJJ_y,pJJ]
            Gy[pJJ_y,pII] = Gy[pII_y,pII]

            Hy[pII_y,pJJ] = Hy[pJJ_y,pJJ]
            Hy[pJJ_y,pII] = Hy[pII_y,pII]
        end
    end

    return nothing
end

function periodic_bcs_R!(grid, Rx, Ry, periodic_x, periodic_y)
    @unpack nx, ny, ind = grid

    if periodic_x
        @inbounds @threads for i in eachindex(grid.ind.periodic_x[1])
            II = grid.ind.periodic_x[1][i]
            JJ = grid.ind.periodic_x[2][i]
            pII = lexicographic(II, ny)
            pJJ = lexicographic(JJ, ny)

            Rx[pII,pJJ] = 1.0
            Rx[pJJ+2*ny,pII] = 1.0
        end
    end
    if periodic_y
        for i in eachindex(grid.ind.periodic_y[1])
            II = grid.ind.periodic_y[1][i]
            JJ = grid.ind.periodic_y[2][i]
            pII = lexicographic(II, ny)
            pII_y = lexicographic(δy⁺(δy⁺(δx⁻(JJ))), ny+2)+1
            pJJ = lexicographic(JJ, ny)
            pJJ_y = lexicographic(δy⁺(δy⁺(JJ)), ny+2)

            Ry[pII_y, pJJ] = 1.0
            Ry[pJJ_y, pII] = 1.0
        end
    end

    return nothing
end


"""
Sets ```iMx_b, iMx_bd, iMy_b, iMy_bd``` with an epsilon parameter and capacities 8, 9, 10, 11 at cell II i.e. dcap[II,8]
"""
function mass_matrix_borders!(num,ind, iMx_b, iMy_b, iMx_bd, iMy_bd, dcap, n)
    @unpack b_left, b_bottom, b_right, b_top = ind

    idx = 1
    @inbounds for II in b_left[1]
        pII = lexicographic(II, n)
        # printstyled(color=:red, @sprintf "\n mass: %.5i mass %.4e eps %.4e mass + eps %.4e inv %.4e inv(m+eps) %.4e\n" idx dcap[II,8] eps(0.01) dcap[II,8]+eps(0.01) 1/dcap[II,8] 1/(dcap[II,8]+eps(0.01)))
        @inbounds iMx_b[pII, idx] = inv_weight_eps(num,dcap[II,8]) 
        @inbounds iMx_bd[idx, idx] = inv_weight_eps(num,dcap[II,8])
        idx += 1
    end
    @inbounds for II in b_bottom[1]
        pII = lexicographic(II, n+1)
        @inbounds iMy_b[pII, idx] = inv_weight_eps(num,dcap[II,9])
        @inbounds iMy_bd[idx, idx] = inv_weight_eps(num,dcap[II,9])
        idx += 1
    end
    @inbounds for II in b_right[1]
        pII = lexicographic(δx⁺(II), n)
        @inbounds iMx_b[pII, idx] = inv_weight_eps(num,dcap[II,10])
        @inbounds iMx_bd[idx, idx] = inv_weight_eps(num,dcap[II,10])
        idx += 1
    end
    @inbounds for II in b_top[1]
        pII = lexicographic(δy⁺(II), n+1)
        @inbounds iMy_b[pII, idx] = inv_weight_eps(num,dcap[II,11])
        @inbounds iMy_bd[idx, idx] = inv_weight_eps(num,dcap[II,11])
        idx += 1
    end

    return nothing
end


"""
Compute Hx and Hy, left: -A1 bottom: -A2 right: A3 top: A4
"""
function bc_matrix_borders!(grid::Mesh{GridCC,T,N}, ind, Hx, Hy, dcap) where {T,N}
    @unpack nx, ny = grid
    @unpack b_left, b_bottom, b_right, b_top = ind

    @inbounds @threads for idx in 1:ny
        # II = CartesianIndex(idx,1)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, b_left[1][idx])
        @inbounds Hx[idx, idx] = -A1
    end
    @inbounds @threads for idx in 1:nx
        # II = CartesianIndex(1,idx)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, b_bottom[1][idx])
        @inbounds Hy[idx+ny, idx+ny] = -A2
    end
    @inbounds @threads for idx in 1:ny
        # II = CartesianIndex(idx,nx)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, b_right[1][idx])
        @inbounds Hx[idx+ny+nx, idx+ny+nx] = A3
    end
    @inbounds @threads for idx in 1:nx
        # II = CartesianIndex(ny,idx)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, b_top[1][idx])
        @inbounds Hy[idx+2*ny+nx, idx+2*ny+nx] = A4
    end

    return nothing
end


function bc_matrix_borders!(grid::Mesh{GridFCx,T,N}, ind, ind_u, Hx, Hy, dcap, dcap_u) where {T,N}
    @unpack nx, ny = grid
    @unpack b_left, b_right = ind
    @unpack b_bottom, b_top = ind_u

    @inbounds @threads for idx in 1:ny
        # II = CartesianIndex(idx,1)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, b_left[1][idx])
        @inbounds Hx[idx, idx] = -A1
    end
    @inbounds @threads for idx in 1:nx
        # II = CartesianIndex(1,idx)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap_u, b_bottom[1][idx])
        @inbounds Hy[idx+ny, idx+ny] = -A2
    end
    @inbounds @threads for idx in 1:ny
        # II = CartesianIndex(idx,nx)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, b_right[1][idx])
        @inbounds Hx[idx+ny+nx, idx+ny+nx] = A3
    end
    @inbounds @threads for idx in 1:nx
        # II = CartesianIndex(ny,idx)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap_u, b_top[1][idx])
        @inbounds Hy[idx+2*ny+nx, idx+2*ny+nx] = A4
    end

    return nothing
end


function bc_matrix_borders!(grid::Mesh{GridFCy,T,N}, ind, ind_v, Hx, Hy, dcap, dcap_v) where {T,N}
    @unpack nx, ny = grid
    @unpack b_bottom, b_top = ind
    @unpack b_left, b_right = ind_v

    @inbounds @threads for idx in 1:ny
        # II = CartesianIndex(idx,1)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap_v, b_left[1][idx])
        @inbounds Hx[idx, idx] = -A1
    end
    @inbounds @threads for idx in 1:nx
        # II = CartesianIndex(1,idx)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, b_bottom[1][idx])
        @inbounds Hy[idx+ny, idx+ny] = -A2
    end
    @inbounds @threads for idx in 1:ny
        # II = CartesianIndex(idx,nx)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap_v, b_right[1][idx])
        @inbounds Hx[idx+ny+nx, idx+ny+nx] = A3
    end
    @inbounds @threads for idx in 1:nx
        # II = CartesianIndex(ny,idx)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, b_top[1][idx])
        @inbounds Hy[idx+2*ny+nx, idx+2*ny+nx] = A4
    end

    return nothing
end


"""
    bc_matrix_borders!(grid, Hx_u, Hy_v, Hx_p, Hy_p, dcap)
    
"""
function bc_matrix_borders!(grid, Hx_u, Hy_v, Hx_p, Hy_p, dcap)
    @unpack nx, ny, ind = grid
    @unpack b_left, b_bottom, b_right, b_top = ind

    @inbounds @threads for i in 1:ny
        II = CartesianIndex(i,1)
        pII = lexicographic(II, ny)
        @inbounds A1 = dcap[II,1]

        @inbounds Hx_u[pII, i] = -A1
        @inbounds Hx_p[pII, i] = -A1
    end
    @inbounds @threads for i in 1:nx
        II = CartesianIndex(1,i)
        pII = lexicographic(II, ny+1)
        pJJ = lexicographic(II, ny)
        @inbounds A2 = dcap[II,2]
        
        @inbounds Hy_v[pII, i+ny] = -A2
        @inbounds Hy_p[pJJ, i+ny+1] = -A2
    end
    @inbounds @threads for i in 1:ny
        II = CartesianIndex(i,nx)
        pII = lexicographic(δx⁺(II), ny)
        pJJ = lexicographic(II, ny)
        @inbounds A3 = dcap[II,3]
        
        @inbounds Hx_u[pII, i+ny+nx] = A3
        @inbounds Hx_p[pJJ, i+ny+nx+1] = A3
    end
    @inbounds @threads for i in 1:nx
        II = CartesianIndex(ny,i)
        pII = lexicographic(δy⁺(II), ny+1)
        pJJ = lexicographic(II, ny)
        @inbounds A4 = dcap[II,4]
        
        @inbounds Hy_v[pII, i+2*ny+nx] = A4
        @inbounds Hy_p[pJJ, i+2*(ny+1)+nx] = A4
    end

    return nothing
end


"""
Modifies Hx and Hy
"""
function periodic_bcs_borders!(grid, Hx, Hy, periodic_x, periodic_y)
    @unpack nx, ny = grid

    if periodic_x
        @inbounds @threads for idx in 1:ny
            @inbounds Hx[idx,idx+nx+ny] = Hx[idx+nx+ny,idx]
        end
    end
    if periodic_y
        @inbounds @threads for idx in ny+1:ny+nx
            @inbounds Hy[idx,idx+nx+ny] = Hy[idx+nx+ny,idx]
        end
    end

    return nothing
end

# function harmonic_average(W4, W3)
#     # Harmonic average of volume capacities
#     if W4 < 1e-8 || W3 < 1e-8
#         Ŵ =  W4 + W3
#     else
#         Ŵ = 2 * W4 * W3 / (W4 + W3)
#     end

#     return Ŵ
# end

# # Didn't implement the inhomogeneous BC term for the moment
# # since we only need this operator to compute forces inside
# # the domain and we have the no-slip condition at the wall
# function strain_rate!(::Dirichlet, O11, O12_x, O12_y, O22, cap_x, cap_y, n, all_indices, inside)
#     @inbounds @threads for II in inside
#         pII = lexicographic(II, n)

#         JJ = δx⁺(II)
#         pJJ = lexicographic(JJ, n)
#         A1_1, A2_1, A3_1, A4_1, B1_1, B2_1, W1_1, W2_1, W3_1, W4_1 = get_capacities(cap_x, II)
#         A1_2, A2_2, A3_2, A4_2, B1_2, B2_2, W1_2, W2_2, W3_2, W4_2 = get_capacities(cap_x, JJ)

#         @inbounds O11[pII,pII] = -B1_1 / W3_1
#         @inbounds O11[pII,pJJ] = B1_2 / W3_1

#         JJ = δy⁺(II)
#         _pII = lexicographic(II, n+1)
#         pJJ = lexicographic(JJ, n+1)
#         A1_1, A2_1, A3_1, A4_1, B1_1, B2_1, W1_1, W2_1, W3_1, W4_1 = get_capacities(cap_y, II)
#         A1_2, A2_2, A3_2, A4_2, B1_2, B2_2, W1_2, W2_2, W3_2, W4_2 = get_capacities(cap_y, JJ)

#         @inbounds O22[pII,_pII] = -B2_1 / W4_1
#         @inbounds O22[pII,pJJ] = B2_2 / W4_1
#     end
#     @inbounds @threads for II in all_indices[1:end-1,2:end]
#         pII = lexicographic(II, n)

#         JJ_1 = II
#         JJ_2 = δy⁺(II)
#         pJJ_1 = lexicographic(JJ_1, n)
#         pJJ_2 = lexicographic(JJ_2, n)
#         A1_1_x, A2_1_x, A3_1_x, A4_1_x, B1_1_x, B2_1_x, W1_1_x, W2_1_x, W3_1_x, W4_1_x = get_capacities(cap_x, JJ_1)
#         A1_2_x, A2_2_x, A3_2_x, A4_2_x, B1_2_x, B2_2_x, W1_2_x, W2_2_x, W3_2_x, W4_2_x = get_capacities(cap_x, JJ_2)

#         @inbounds O12_x[pII,pJJ_1] = -B2_1_x / 2W4_1_x
#         @inbounds O12_x[pII,pJJ_2] = B2_2_x / 2W4_1_x

#         KK_1 = δy⁺(δx⁻(II))
#         KK_2 = δy⁺(II)
#         pKK_1 = lexicographic(KK_1, n+1)
#         pKK_2 = lexicographic(KK_2, n+1)
#         A1_1_y, A2_1_y, A3_1_y, A4_1_y, B1_1_y, B2_1_y, W1_1_y, W2_1_y, W3_1_y, W4_1_y = get_capacities(cap_y, KK_1)
#         A1_2_y, A2_2_y, A3_2_y, A4_2_y, B1_2_y, B2_2_y, W1_2_y, W2_2_y, W3_2_y, W4_2_y = get_capacities(cap_y, KK_2)

#         Ŵ = harmonic_average(W4_1_x, W3_1_y)

#         @inbounds O12_y[pII,pKK_1] = -B1_1_y / 2Ŵ
#         @inbounds O12_y[pII,pKK_2] = B1_2_y / 2Ŵ
#     end

#     return nothing
# end

# function set_bc_bnds(::Dirichlet, Du, Dv, Hu, Hv, u, v, BC_u, BC_v)
#     Dx = copy(Du)
#     Dy = copy(Dv)

#     if is_neumann(BC_u.left.t)
#         @inbounds Dx[:,1] .= u[:,1] .+ Hu[:,1] .* BC_u.left.val
#     elseif is_dirichlet(BC_u.left.t)
#         @inbounds Dx[:,1] .= BC_u.left.val
#     elseif is_periodic(BC_u.left.t)
#         @inbounds Dx[:,1] .= u[:,end]
#     end
#     if is_neumann(BC_v.bottom.t)
#         @inbounds Dy[1,:] .= v[1,:] .+ Hv[1,:] .* BC_v.bottom.val 
#     elseif is_dirichlet(BC_v.bottom.t)
#         @inbounds Dy[1,:] .= BC_v.bottom.val
#     elseif is_periodic(BC_v.bottom.t)
#         @inbounds Dy[1,:] .= v[end,:]
#     end
#     if is_neumann(BC_u.right.t)
#         @inbounds Dx[:,end] .= u[:,end] .+ Hu[:,end] .* BC_u.right.val 
#     elseif is_dirichlet(BC_u.right.t)
#         @inbounds Dx[:,end] .= BC_u.right.val
#     elseif is_periodic(BC_u.right.t)
#         @inbounds Dx[:,end] .= u[:,1]
#     end
#     if is_neumann(BC_v.top.t)
#         @inbounds Dy[end,:] .= v[end,:] .+ Hv[end,:] .* BC_v.top.val 
#     elseif is_dirichlet(BC_v.top.t)
#         @inbounds Dy[end,:] .= BC_v.top.val
#     elseif is_periodic(BC_v.top.t)
#         @inbounds Dy[end,:] .= v[1,:]
#     end

#     return Dx, Dy
# end

# @inline function set_sca_conv_bnd!(::Dirichlet, ::Dirichlet, O, fun, A1, A2, B1, D, n, b_indices, b_periodic)
#     return nothing
# end

# @inline function set_sca_conv_bnd!(::Dirichlet, ::Neumann, O, fun, A1, A2, B1, D, n, b_indices, b_periodic)
#     @inbounds @threads for II in b_indices
#         pII = lexicographic(II, n)
#         @inbounds O[pII,pII] += -0.5 * ((A2[II] - B1[II]) * D[fun(II)] + (B1[II] - A1[II]) * D[II])
#     end
#     return nothing
# end

# @inline function set_sca_conv_bnd!(::Dirichlet, ::Periodic, O, fun, A1, A2, B1, D, n, b_indices, b_periodic)
#     @inbounds for (II, JJ) in zip(b_indices, b_periodic)
#         pII = lexicographic(II, n)
#         pJJ = lexicographic(JJ, n)
#         @inbounds O[pII,pJJ] += -0.5 * ((A2[II] - B1[II]) * D[fun(II)] + (B1[II] - A1[II]) * D[II])
#     end
#     return nothing
# end

# function scalar_convection!(::Dirichlet, O, B, u, v, Dx, Dy, Du, Dv, cap, n, BC, inside, b_left, b_bottom, b_right, b_top)
#     B .= 0.0
#     @inbounds @threads for II in inside
#         pII = lexicographic(II, n)
#         A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II)
#         u1, v2, u3, v4 = u[II], v[II], u[δx⁺(II)], v[δy⁺(II)]

#         @inbounds O[pII,pII] = 0.5 * (A3 * u3 - A1 * u1 + A4 * v4 - A2 * v2)
#         @inbounds O[pII,pII+n] = 0.5 * A3 * u3
#         @inbounds O[pII,pII-n] = -0.5 * A1 * u1
#         @inbounds O[pII,pII+1] = 0.5 * A4 * v4
#         @inbounds O[pII,pII-1] = -0.5 * A2 * v2

#         @inbounds O[pII,pII] += -0.5 * ((A3 - B1) * Du[δx⁺(II)] + (B1 - A1) * Du[II])
#         @inbounds O[pII,pII] += -0.5 * ((A4 - B2) * Dv[δy⁺(II)] + (B2 - A2) * Dv[II])

#         @inbounds B[pII] += -0.5 * Dx[II] * ((A3 - B1) * Du[δx⁺(II)] + (B1 - A1) * Du[II])
#         @inbounds B[pII] += -0.5 * Dy[II] * ((A4 - B2) * Dv[δy⁺(II)] + (B2 - A2) * Dv[II])
#     end

#     @inbounds @threads for II in vcat(b_left, b_bottom[2:end-1], b_right, b_top[2:end-1])
#         pII = lexicographic(II, n)
#         A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II)
#         u1, v2, u3, v4 = u[II], v[II], u[δx⁺(II)], v[δy⁺(II)]

#         @inbounds O[pII,pII] = 0.5 * (A3 * u3 - A1 * u1 + A4 * v4 - A2 * v2)

#         @inbounds O[pII,pII] += -0.5 * ((A3 - B1) * Du[δx⁺(II)] + (B1 - A1) * Du[II])
#         @inbounds O[pII,pII] += -0.5 * ((A4 - B2) * Dv[δy⁺(II)] + (B2 - A2) * Dv[II])

#         @inbounds B[pII] += -0.5 * Dx[II] * ((A3 - B1) * Du[δx⁺(II)] + (B1 - A1) * Du[II])
#         @inbounds B[pII] += -0.5 * Dy[II] * ((A4 - B2) * Dv[δy⁺(II)] + (B2 - A2) * Dv[II])
#     end
#     @inbounds @threads for II in vcat(b_left, b_bottom[2:end-1], b_top[2:end-1])
#         pII = lexicographic(II, n)
#         A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II)

#         @inbounds O[pII,pII+n] = 0.5 * A3 * u[δx⁺(II)]
#     end
#     @inbounds @threads for II in vcat(b_bottom[2:end-1], b_right, b_top[2:end-1])
#         pII = lexicographic(II, n)
#         A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II)

#         @inbounds O[pII,pII-n] = -0.5 * A1 * u[II]
#     end
#     @inbounds @threads for II in vcat(b_left[2:end-1], b_bottom, b_right[2:end-1])
#         pII = lexicographic(II, n)
#         A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II)

#         @inbounds O[pII,pII+1] = 0.5 * A4 * v[δy⁺(II)]
#     end
#     @inbounds @threads for II in vcat(b_left[2:end-1], b_right[2:end-1], b_top)
#         pII = lexicographic(II, n)
#         A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II)

#         @inbounds O[pII,pII-1] = -0.5 * A2 * v[II]
#     end

#     if is_periodic(BC.left.t) && is_periodic(BC.right.t)
#         @inbounds for (II,JJ) in zip(b_right, b_left)
#             pII = lexicographic(II, n)
#             pJJ = lexicographic(JJ, n)
#             A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II)
    
#             @inbounds O[pII,pJJ] = 0.5 * A3 * u[δx⁺(II)]
#         end
#         @inbounds for (II,JJ) in zip(b_left, b_right)
#             pII = lexicographic(II, n)
#             pJJ = lexicographic(JJ, n)
#             A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II)
    
#             @inbounds O[pII,pJJ] = -0.5 * A1 * u[II]
#         end
#     end
#     if is_periodic(BC.bottom.t) && is_periodic(BC.top.t)
#         @inbounds for (II,JJ) in zip(b_top, b_bottom)
#             pII = lexicographic(II, n)
#             pJJ = lexicographic(JJ, n)
#             A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II)
    
#             @inbounds O[pII,pJJ] = 0.5 * A4 * v[δy⁺(II)]
#         end
#         @inbounds for (II,JJ) in zip(b_bottom, b_top)
#             pII = lexicographic(II, n)
#             pJJ = lexicographic(JJ, n)
#             A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II)
    
#             @inbounds O[pII,pJJ] = -0.5 * A2 * v[II]
#         end
#     end

#     @inbounds _A1 = @view cap[:,:,1]
#     @inbounds _A2 = @view cap[:,:,2]
#     @inbounds _A3 = @view cap[:,:,3]
#     @inbounds _A4 = @view cap[:,:,4]
#     @inbounds _B1 = @view cap[:,:,6]
#     @inbounds _B2 = @view cap[:,:,7]

#     set_sca_conv_bnd!(dir, BC.left.t, O, δx⁺, _A1, _A3, _B1, Du, n, b_left, b_right)
#     set_sca_conv_bnd!(dir, BC.bottom.t, O, δy⁺, _A2, _A4, _B2, Dv, n, b_bottom, b_top)
#     set_sca_conv_bnd!(dir, BC.right.t, O, δx⁺, _A1, _A3, _B1, Du, n, b_right, b_left)
#     set_sca_conv_bnd!(dir, BC.top.t, O, δy⁺, _A2, _A4, _B2, Dv, n, b_top, b_bottom)

#     return nothing
# end

# function set_bc_bnds(::Dirichlet, ::Union{Type{GridFCx},Type{GridFCy}}, Du, Dv, Hu, Hv, u, v, BC_u, BC_v)
#     Du1_x = copy(Du)
#     Du1_y = copy(Du)
#     Du2_x = copy(Du)
#     Du2_y = copy(Du)
#     Dv_x = copy(Dv)
#     Dv_y = copy(Dv)

#     if is_neumann(BC_u.left.t)
#         @inbounds Du1_x[:,1] .= u[:,1] .+ Hu[:,1] .* BC_u.left.val
#         @inbounds Du1_x[:,2] .= u[:,2]
#         @inbounds Du2_x[:,1] .= Hu[:,1] .* BC_u.left.val
#         @inbounds Du2_x[:,2] .= u[:,2]
#     elseif is_dirichlet(BC_u.left.t)
#         @inbounds Du1_x[:,1] .= BC_u.left.val
#         @inbounds Du1_x[:,2] .= BC_u.left.val
#         @inbounds Du2_x[:,1] .= BC_u.left.val
#         @inbounds Du2_x[:,2] .= BC_u.left.val
#     elseif is_periodic(BC_u.left.t)
#         @inbounds Du1_x[:,1] .= u[:,end]
#         @inbounds Du1_x[:,2] .= u[:,2]
#         @inbounds Du2_x[:,1] .= 0.0
#         @inbounds Du2_x[:,2] .= u[:,2]
#     end
#     if is_neumann(BC_u.bottom.t)
#         @inbounds Du1_y[1,:] .= u[1,:] .+ Hu[1,:] .* BC_u.bottom.val
#         @inbounds Du1_y[2,:] .= u[2,:]
#         @inbounds Du2_y[1,:] .= Hu[1,:] .* BC_u.bottom.val
#         @inbounds Du2_y[2,:] .= u[2,:]
#     elseif is_dirichlet(BC_u.bottom.t)
#         @inbounds Du1_y[1,:] .= BC_u.bottom.val
#         @inbounds Du1_y[2,:] .= BC_u.bottom.val
#         @inbounds Du2_y[1,:] .= BC_u.bottom.val
#         @inbounds Du2_y[2,:] .= BC_u.bottom.val
#     elseif is_periodic(BC_u.bottom.t)
#         @inbounds Du1_y[1,:] .= u[end,:]
#         @inbounds Du1_y[2,:] .= u[2,:]
#         @inbounds Du2_y[1,:] .= 0.0
#         @inbounds Du2_y[2,:] .= u[2,:]
#     end
#     if is_neumann(BC_u.right.t)
#         @inbounds Du1_x[:,end] .= u[:,end] .+ Hu[:,end] .* BC_u.right.val 
#         @inbounds Du1_x[:,end-1] .= u[:,end-1]
#         @inbounds Du2_x[:,end] .= Hu[:,end] .* BC_u.right.val
#         @inbounds Du2_x[:,end-1] .= u[:,end-1]
#     elseif is_dirichlet(BC_u.right.t)
#         @inbounds Du1_x[:,end] .= BC_u.right.val
#         @inbounds Du1_x[:,end-1] .= BC_u.right.val
#         @inbounds Du2_x[:,end] .= BC_u.right.val
#         @inbounds Du2_x[:,end-1] .= BC_u.right.val
#     elseif is_periodic(BC_u.right.t)
#         @inbounds Du1_x[:,end] .= u[:,1]
#         @inbounds Du1_x[:,end-1] .= u[:,end-1]
#         @inbounds Du2_x[:,end] .= 0.0
#         @inbounds Du2_x[:,end-1] .= u[:,end-1]
#     end
#     if is_neumann(BC_u.top.t)
#         @inbounds Du1_y[end,:] .= u[end,:] .+ Hu[end,:] .* BC_u.top.val
#         @inbounds Du1_y[end-1,:] .= u[end-1,:]
#         @inbounds Du2_y[end,:] .= Hu[end,:] .* BC_u.top.val
#         @inbounds Du2_y[end-1,:] .= u[end-1,:]
#     elseif is_dirichlet(BC_u.top.t)
#         @inbounds Du1_y[end,:] .= BC_u.top.val
#         @inbounds Du1_y[end-1,:] .= BC_u.top.val
#         @inbounds Du2_y[end,:] .= BC_u.top.val
#         @inbounds Du2_y[end-1,:] .= BC_u.top.val
#     elseif is_periodic(BC_u.top.t)
#         @inbounds Du1_y[end,:] .= u[1,:]
#         @inbounds Du1_y[end-1,:] .= u[end-1,:]
#         @inbounds Du2_y[end,:] .= 0.0
#         @inbounds Du2_y[end-1,:] .= u[end-1,:]
#     end

#     if is_neumann(BC_v.left.t)
#         @inbounds Dv_x[:,1] .= v[:,1] .+ Hv[:,1] .* BC_v.left.val 
#         @inbounds Dv_x[:,2] .= v[:,2]
#     elseif is_dirichlet(BC_v.left.t)
#         @inbounds Dv_x[:,1] .= BC_v.left.val
#         @inbounds Dv_x[:,2] .= BC_v.left.val
#     elseif is_periodic(BC_v.left.t)
#         @inbounds Dv_x[:,1] .= v[:,end]
#         @inbounds Dv_x[:,2] .= v[:,2]
#     end
#     if is_neumann(BC_v.bottom.t)
#         @inbounds Dv_y[1,:] .= v[1,:] .+ Hv[1,:] .* BC_v.bottom.val 
#         @inbounds Dv_y[2,:] .= v[2,:]
#     elseif is_dirichlet(BC_v.bottom.t)
#         @inbounds Dv_y[1,:] .= BC_v.bottom.val
#         @inbounds Dv_y[2,:] .= BC_v.bottom.val
#     elseif is_periodic(BC_v.bottom.t)
#         @inbounds Dv_y[1,:] .= v[end,:]
#         @inbounds Dv_y[2,:] .= v[2,:]
#     end
#     if is_neumann(BC_v.right.t)
#         @inbounds Dv_x[:,end] .= v[:,end] .+ Hv[:,end] .* BC_v.right.val 
#         @inbounds Dv_x[:,end-1] .= v[:,end-1]
#     elseif is_dirichlet(BC_v.right.t)
#         @inbounds Dv_x[:,end] .= BC_v.right.val
#         @inbounds Dv_x[:,end-1] .= BC_v.right.val
#     elseif is_periodic(BC_v.right.t)
#         @inbounds Dv_x[:,end] .= v[:,1]
#         @inbounds Dv_x[:,end-1] .= v[:,end-1]
#     end
#     if is_neumann(BC_v.top.t)
#         @inbounds Dv_y[end,:] .= v[end,:] .+ Hv[end,:] .* BC_v.top.val 
#         @inbounds Dv_y[end-1,:] .= v[end-1,:]
#     elseif is_dirichlet(BC_v.top.t)
#         @inbounds Dv_y[end,:] .= BC_v.top.val
#         @inbounds Dv_y[end-1,:] .= BC_v.top.val
#     elseif is_periodic(BC_v.top.t)
#         @inbounds Dv_y[end,:] .= v[1,:]
#         @inbounds Dv_y[end-1,:] .= v[end-1,:]
#     end

#     return Du1_x, Du1_y, Du2_x, Du2_y, Dv_x, Dv_y
# end

# @inline function set_vec_conv_bnd!(::Dirichlet, ::Dirichlet, O, fun_cap, fun1, fun2, fun3, fun4, A1, A2, A3, A4, B1, B2, Du, Dv, n, b_indices, b_periodic)
#     return nothing
# end

# @inline function set_vec_conv_bnd!(::Dirichlet, ::Neumann, O, fun_cap, fun1, fun2, fun3, fun4, A3, B1, A1, A4, B2, A2, Du, Dv, n, b_indices, b_periodic)
#     @inbounds @threads for II in b_indices
#         pII = lexicographic(II, n)
#         @inbounds O[pII,pII] += -0.25 * (A3[fun_cap(II)] - B1[fun_cap(II)]) * Du[fun1(II)]
#         @inbounds O[pII,pII] += -0.25 * (B1[fun_cap(II)] - A1[fun_cap(II)]) * Du[fun2(II)]
#         @inbounds O[pII,pII] += -0.25 * (A4[fun_cap(II)] - B2[fun_cap(II)]) * Dv[fun3(II)]
#         @inbounds O[pII,pII] += -0.25 * (B2[fun_cap(II)] - A2[fun_cap(II)]) * Dv[fun4(II)]
#     end
#     return nothing
# end

# @inline function set_vec_conv_bnd!(::Dirichlet, ::Periodic, O, fun_cap, fun1, fun2, fun3, fun4, A1, A2, A3, A4, B1, B2, Du, Dv, n, b_indices, b_periodic)
#     @inbounds for (II,JJ) in zip(b_indices, b_periodic)
#         pII = lexicographic(II, n)
#         pJJ = lexicographic(JJ, n)
#         @inbounds O[pII,pJJ] = -0.25 * (A1[fun_cap(II)] - B1[fun_cap(II)]) * D[fun1(II)]
#         @inbounds O[pII,pJJ] = -0.25 * (B1[fun_cap(II)] - A2[fun_cap(II)]) * D[fun2(II)]
#     end
#     return nothing
# end

# function fill_inside_conv!(::Type{GridFCx}, O, B, u, v, Du1_x, Du1_y, Dv_y, cap, n, II)
#     pII = lexicographic(II, n)
#     A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(II))
#     A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)

#     Auim1, Aui, Auip1 = A1_1 * u[δx⁻(II)], A1_2 * u[II], A3_2 * u[δx⁺(II)]
#     Avim1jm1, Avim1jp1 = A2_1 * v[δx⁻(II)], A4_1 * v[δy⁺(δx⁻(II))]
#     Avip1jm1, Avip1jp1 = A2_2 * v[II], A4_2 * v[δy⁺(II)]

#     Au1 = 0.5 * (Auim1 + Aui)
#     Au2 = 0.5 * (Avim1jm1 + Avip1jm1)
#     Au3 = 0.5 * (Aui + Auip1)
#     Au4 = 0.5 * (Avim1jp1 + Avip1jp1)

#     @inbounds O[pII,pII] = 0.5 * (Au3 - Au1 + Au4 - Au2)
#     @inbounds O[pII,pII+n] = 0.5 * Au3
#     @inbounds O[pII,pII-n] = -0.5 * Au1
#     @inbounds O[pII,pII+1] = 0.5 * Au4
#     @inbounds O[pII,pII-1] = -0.5 * Au2

#     @inbounds O[pII,pII] += -0.25 * (A3_2 - B1_2) * Du1_x[δx⁺(II)]
#     @inbounds O[pII,pII] += -0.25 * (B1_2 - B1_1) * Du1_x[II]
#     @inbounds O[pII,pII] += -0.25 * (B1_1 - A1_1) * Du1_x[δx⁻(II)]

#     @inbounds O[pII,pII] += -0.25 * (A4_1 - B2_1) * Dv_y[δx⁻(δy⁺(II))]
#     @inbounds O[pII,pII] += -0.25 * (B2_1 - A2_1) * Dv_y[δx⁻(II)]
#     @inbounds O[pII,pII] += -0.25 * (A4_2 - B2_2) * Dv_y[δy⁺(II)]
#     @inbounds O[pII,pII] += -0.25 * (B2_2 - A2_2) * Dv_y[II]

#     # @inbounds B[pII] += -0.25 * Du2_x[II] * (A3_2 - B1_2) * Du1_x[δx⁺(II)]
#     # @inbounds B[pII] += -0.25 * Du2_x[II] * (B1_2 - B1_1) * Du1_x[II]
#     # @inbounds B[pII] += -0.25 * Du2_x[II] * (B1_1 - A1_1) * Du1_x[δx⁻(II)]

#     # @inbounds B[pII] += -0.25 * Du2_y[II] * (A4_1 - B2_1) * Dv_y[δx⁻(δy⁺(II))]
#     # @inbounds B[pII] += -0.25 * Du2_y[II] * (B2_1 - A2_1) * Dv_y[δx⁻(II)]
#     # @inbounds B[pII] += -0.25 * Du2_y[II] * (A4_2 - B2_2) * Dv_y[δy⁺(II)]
#     # @inbounds B[pII] += -0.25 * Du2_y[II] * (B2_2 - A2_2) * Dv_y[II]

#     @inbounds B[pII] += -0.25 * Du1_x[II] * (A3_2 - B1_2) * Du1_x[δx⁺(II)]
#     @inbounds B[pII] += -0.25 * Du1_x[II] * (B1_2 - B1_1) * Du1_x[II]
#     @inbounds B[pII] += -0.25 * Du1_x[II] * (B1_1 - A1_1) * Du1_x[δx⁻(II)]

#     @inbounds B[pII] += -0.25 * Du1_y[II] * (A4_1 - B2_1) * Dv_y[δx⁻(δy⁺(II))]
#     @inbounds B[pII] += -0.25 * Du1_y[II] * (B2_1 - A2_1) * Dv_y[δx⁻(II)]
#     @inbounds B[pII] += -0.25 * Du1_y[II] * (A4_2 - B2_2) * Dv_y[δy⁺(II)]
#     @inbounds B[pII] += -0.25 * Du1_y[II] * (B2_2 - A2_2) * Dv_y[II]

#     return nothing
# end

# function vector_convection!(::Dirichlet, ::Type{GridFCx}, O, B, u, v, Du1_x, Du1_y, Du2_x, Du2_y, Dv_x, Dv_y, cap, n, BC, inside, b_left, b_bottom, b_right, b_top)
#     B .= 0.0
#     @inbounds @threads for II in inside
#         fill_inside_conv!(GridFCx, O, B, u, v, Du1_x, Du1_y, Dv_y, cap, n, II)
#     end

#     @inbounds @threads for II in vcat(b_left, b_bottom[2:end-1], b_right, b_top[2:end-1])
#         pII = lexicographic(II, n)
#         @inbounds O[pII,pII] = 0.0
#     end
#     bnds = (b_left, b_bottom[2:end-1], b_top[2:end-1])
#     bc = ((Du1_x, Du1_x, Dv_x), (Du1_y, Du1_y, Dv_y), (Du1_y, Du1_y, Dv_y))
#     for (bnd, (Du1, Du2, Dv)) in zip(bnds, bc)
#         @inbounds for II in bnd
#             pII = lexicographic(II, n)
#             A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
            
#             Aui, Auip1 = A1_2 * u[II], A3_2 * u[δx⁺(II)]

#             Au3 = 0.5 * (Aui + Auip1)

#             @inbounds O[pII,pII] += 0.5 * Au3
#             @inbounds O[pII,pII+n] = 0.5 * Au3

#             @inbounds O[pII,pII] += -0.25 * (A3_2 - B1_2) * Du1[δx⁺(II)]
#             @inbounds O[pII,pII] += -0.25 * (B1_2 - A1_2) * Du1[II]

#             @inbounds O[pII,pII] += -0.25 * (A4_2 - B2_2) * Dv[δy⁺(II)]
#             @inbounds O[pII,pII] += -0.25 * (B2_2 - A2_2) * Dv[II]

#             @inbounds B[pII] += -0.25 * Du2[II] * (A3_2 - B1_2) * Du1[δx⁺(II)]
#             @inbounds B[pII] += -0.25 * Du2[II] * (B1_2 - A1_2) * Du1[II]

#             @inbounds B[pII] += -0.25 * Du2[II] * (A4_2 - B2_2) * Dv[δy⁺(II)]
#             @inbounds B[pII] += -0.25 * Du2[II] * (B2_2 - A2_2) * Dv[II]
#         end
#     end
#     bnds = (b_bottom[2:end-1], b_right, b_top[2:end-1])
#     bc = ((Du1_y, Du1_y, Dv_y), (Du1_x, Du1_x, Dv_x), (Du1_y, Du1_y, Dv_y))
#     for (bnd, (Du1, Du2, Dv)) in zip(bnds, bc)
#         @inbounds for II in bnd
#             pII = lexicographic(II, n)
#             A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(II))
            
#             Auim1, Aui = A1_1 * u[δx⁻(II)], A3_1 * u[II]

#             Au1 = 0.5 * (Auim1 + Aui)

#             @inbounds O[pII,pII] += -0.5 * Au1
#             @inbounds O[pII,pII-n] = -0.5 * Au1

#             @inbounds O[pII,pII] += -0.25 * (A3_1 - B1_1) * Du1[II]
#             @inbounds O[pII,pII] += -0.25 * (B1_1 - A1_1) * Du1[δx⁻(II)]

#             @inbounds O[pII,pII] += -0.25 * (A4_1 - B2_1) * Dv[δx⁻(δy⁺(II))]
#             @inbounds O[pII,pII] += -0.25 * (B2_1 - A2_1) * Dv[δx⁻(II)]

#             @inbounds B[pII] += -0.25 * Du2[II] * (A3_1 - B1_1) * Du1[II]
#             @inbounds B[pII] += -0.25 * Du2[II] * (B1_1 - A1_1) * Du1[δx⁻(II)]

#             @inbounds B[pII] += -0.25 * Du2[II] * (A4_1 - B2_1) * Dv[δx⁻(δy⁺(II))]
#             @inbounds B[pII] += -0.25 * Du2[II] * (B2_1 - A2_1) * Dv[δx⁻(II)]
#         end
#     end
#     @inbounds for II in b_bottom[2:end-1]
#         pII = lexicographic(II, n)
#         A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(II))
#         A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
        
#         Avim1jp1 = A4_1 * v[δy⁺(δx⁻(II))]
#         Avip1jp1 = A4_2 * v[δy⁺(II)]

#         Au4 = 0.5 * (Avim1jp1 + Avip1jp1)

#         @inbounds O[pII,pII] += 0.5 * Au4
#         @inbounds O[pII,pII+1] = 0.5 * Au4
#     end
#     @inbounds for II in b_top[2:end-1]
#         pII = lexicographic(II, n)
#         A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(II))
#         A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
        
#         Avim1jm1 = A2_1 * v[δx⁻(II)]
#         Avip1jm1 = A2_2 * v[II]

#         Au2 = 0.5 * (Avim1jm1 + Avip1jm1)

#         @inbounds O[pII,pII] += -0.5 * Au2
#         @inbounds O[pII,pII-1] = -0.5 * Au2
#     end
#     @inbounds for II in b_left[2:end]
#         pII = lexicographic(II, n)
#         A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
        
#         Avip1jm1 = A2_2 * v[II]

#         Au2 = 0.5 * Avip1jm1

#         if is_periodic(BC.left.t)
#             JJ = II + CartesianIndex(0, n)
#             A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, JJ)
#             Avim1jm1 = A2_1 * v[δx⁻(JJ)]
#             Au2 += 0.5 * Avim1jm1
#         end

#         @inbounds O[pII,pII] += -0.5 * Au2
#         @inbounds O[pII,pII-1] = -0.5 * Au2
#     end
#     @inbounds for II in b_left[1:end-1]
#         pII = lexicographic(II, n)
#         A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
        
#         Avip1jp1 = A4_2 * v[δy⁺(II)]

#         Au4 = 0.5 * Avip1jp1

#         if is_periodic(BC.left.t)
#             JJ = II + CartesianIndex(0, n)
#             A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, JJ)
#             Avim1jp1 = A4_1 * v[δy⁺(δx⁻(JJ))]
#             Au4 += 0.5 * Avim1jp1
#         end

#         @inbounds O[pII,pII] += 0.5 * Au4
#         @inbounds O[pII,pII+1] = 0.5 * Au4
#     end
#     @inbounds for II in b_right[2:end]
#         pII = lexicographic(II, n)
#         A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(II))
        
#         Avim1jm1 = A2_1 * v[δx⁻(II)]

#         Au2 = 0.5 * Avim1jm1

#         if is_periodic(BC.right.t)
#             JJ = II + CartesianIndex(0, -n)
#             A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ)
#             Avip1jm1 = A2_2 * v[JJ]
#             Au2 += 0.5 * Avip1jm1
#         end

#         @inbounds O[pII,pII] += -0.5 * Au2
#         @inbounds O[pII,pII-1] = -0.5 * Au2
#     end
#     @inbounds for II in b_right[1:end-1]
#         pII = lexicographic(II, n)
#         A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(II))
        
#         Avim1jp1 = A4_1 * v[δy⁺(δx⁻(II))]

#         Au4 = 0.5 * Avim1jp1

#         if is_periodic(BC.right.t)
#             JJ = II + CartesianIndex(0, -n)
#             A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ)
#             Avip1jp1 = A4_2 * v[δy⁺(JJ)]
#             Au4 += 0.5 * Avip1jp1
#         end

#         @inbounds O[pII,pII] += 0.5 * Au4
#         @inbounds O[pII,pII+1] = 0.5 * Au4
#     end

#     if is_periodic(BC.left.t) && is_periodic(BC.right.t)
#         (Du1, Du2, Dv) = (Du1_x, Du1_x, Dv_x)
#         @inbounds for (II,JJ) in zip(b_left, b_right)
#             pII = lexicographic(II, n)
#             pJJ = lexicographic(JJ, n)
#             A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(JJ))
            
#             Auim1, Aui = A1_1 * u[JJ], A3_1 * u[II]

#             Au1 = 0.5 * (Auim1 + Aui)

#             @inbounds O[pII,pII] += -0.5 * Au1
#             @inbounds O[pII,pJJ] = -0.5 * Au1

#             @inbounds O[pII,pII] += -0.25 * (A3_1 - B1_1) * Du1[II]
#             @inbounds O[pII,pII] += -0.25 * (B1_1 - A1_1) * Du1[JJ]

#             @inbounds O[pII,pII] += -0.25 * (A4_1 - B2_1) * Dv[δy⁺(JJ)]
#             @inbounds O[pII,pII] += -0.25 * (B2_1 - A2_1) * Dv[JJ]

#             @inbounds B[pII] += -0.25 * Du2[II] * (A3_1 - B1_1) * Du1[II]
#             @inbounds B[pII] += -0.25 * Du2[II] * (B1_1 - A1_1) * Du1[JJ]

#             @inbounds B[pII] += -0.25 * Du2[II] * (A4_1 - B2_1) * Dv[δy⁺(JJ)]
#             @inbounds B[pII] += -0.25 * Du2[II] * (B2_1 - A2_1) * Dv[JJ]
#         end
#         (Du1, Du2, Dv) = (Du1_x, Du1_x, Dv_x)
#         @inbounds for (II, JJ) in zip(b_right, b_left)
#             pII = lexicographic(II, n)
#             pJJ = lexicographic(JJ, n)
#             A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ)
            
#             Aui, Auip1 = A1_2 * u[II], A3_2 * u[JJ]

#             Au3 = 0.5 * (Aui + Auip1)

#             @inbounds O[pII,pII] += 0.5 * Au3
#             @inbounds O[pII,pJJ] = 0.5 * Au3

#             @inbounds O[pII,pII] += -0.25 * (A3_2 - B1_2) * Du1[JJ]
#             @inbounds O[pII,pII] += -0.25 * (B1_2 - A1_2) * Du1[II]

#             @inbounds O[pII,pII] += -0.25 * (A4_2 - B2_2) * Dv[δy⁺(II)]
#             @inbounds O[pII,pII] += -0.25 * (B2_2 - A2_2) * Dv[II]

#             @inbounds B[pII] += -0.25 * Du2[II] * (A3_2 - B1_2) * Du1[JJ]
#             @inbounds B[pII] += -0.25 * Du2[II] * (B1_2 - A1_2) * Du1[II]

#             @inbounds B[pII] += -0.25 * Du2[II] * (A4_2 - B2_2) * Dv[δy⁺(II)]
#             @inbounds B[pII] += -0.25 * Du2[II] * (B2_2 - A2_2) * Dv[II]
#         end
#     end
#     if is_periodic(BC.bottom.t) && is_periodic(BC.top.t)
#         @inbounds for (II,JJ) in zip(b_bottom[2:end-1], b_top[2:end-1])
#             pII = lexicographic(II, n)
#             pJJ = lexicographic(JJ, n)
#             A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(II))
#             A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
            
#             Avim1jm1 = A2_1 * v[δx⁻(II)]
#             Avip1jm1 = A2_2 * v[II]
    
#             Au2 = 0.5 * (Avim1jm1 + Avip1jm1)
    
#             @inbounds O[pII,pII] += -0.5 * Au2
#             @inbounds O[pII,pJJ] = -0.5 * Au2
#         end
#         @inbounds for (II,JJ) in zip(b_top[2:end-1], b_bottom[2:end-1])
#             pII = lexicographic(II, n)
#             pJJ = lexicographic(JJ, n)
#             A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(δx⁻(II)))
#             A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, δy⁻(II))
            
#             Avim1jp1 = A4_1 * v[δx⁻(JJ)]
#             Avip1jp1 = A4_2 * v[JJ]

#             Au4 = 0.5 * (Avim1jp1 + Avip1jp1)
    
#             @inbounds O[pII,pII] += 0.5 * Au4
#             @inbounds O[pII,pJJ] = 0.5 * Au4
#         end

#         ii = b_left[1]
#         pii = lexicographic(ii, n)
#         A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, ii)
        
#         Avip1jm1 = A2_2 * v[ii]

#         Au2 = 0.5 * Avip1jm1

#         if is_periodic(BC.left.t)
#             JJ = ii + CartesianIndex(0, n)
#             A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, JJ)
#             Avim1jm1 = A2_1 * v[δx⁻(JJ)]
#             Au2 += 0.5 * Avim1jm1
#         end

#         JJ = ii + CartesianIndex(n-1, 0)
#         pJJ = lexicographic(JJ, n)
#         @inbounds O[pii,pii] += -0.5 * Au2
#         @inbounds O[pii,pJJ] = -0.5 * Au2
        
#         ii = b_left[end]
#         pii = lexicographic(ii, n)
#         A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, ii)
        
#         Avip1jp1 = A4_2 * v[δy⁺(ii)]

#         Au4 = 0.5 * Avip1jp1

#         if is_periodic(BC.left.t)
#             JJ = ii + CartesianIndex(0, n)
#             A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, JJ)
#             Avim1jp1 = A4_1 * v[δy⁺(δx⁻(JJ))]
#             Au4 += 0.5 * Avim1jp1
#         end

#         JJ = ii + CartesianIndex(-n+1, 0)
#         pJJ = lexicographic(JJ, n)
#         @inbounds O[pii,pii] += 0.5 * Au4
#         @inbounds O[pii,pJJ] = 0.5 * Au4

#         ii = b_right[1]
#         pii = lexicographic(ii, n)
#         A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(ii))
        
#         Avim1jm1 = A2_1 * v[δx⁻(ii)]

#         Au2 = 0.5 * Avim1jm1

#         if is_periodic(BC.right.t)
#             JJ = ii + CartesianIndex(0, -n)
#             A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ)
#             Avip1jm1 = A2_2 * v[JJ]
#             Au2 += 0.5 * Avip1jm1
#         end

#         JJ = ii + CartesianIndex(n-1, 0)
#         pJJ = lexicographic(JJ, n)
#         @inbounds O[pii,pii] += -0.5 * Au2
#         @inbounds O[pii,pJJ] = -0.5 * Au2
        
#         ii = b_right[end]
#         pii = lexicographic(ii, n)
#         A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(ii))
        
#         Avim1jp1 = A4_1 * v[δy⁺(δx⁻(ii))]

#         Au4 = 0.5 * Avim1jp1

#         if is_periodic(BC.right.t)
#             JJ = ii + CartesianIndex(0, -n)
#             A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ)
#             Avip1jp1 = A4_2 * v[δy⁺(JJ)]
#             Au4 += 0.5 * Avip1jp1
#         end
        
#         JJ = ii + CartesianIndex(-n+1, 0)
#         pJJ = lexicographic(JJ, n)
#         @inbounds O[pii,pii] += 0.5 * Au4
#         @inbounds O[pii,pJJ] = 0.5 * Au4
#     end

#     return nothing
# end

# function fill_inside_conv!(::Type{GridFCy}, O, B, u, v, Du_x, Dv1_x, Dv1_y, cap, n, II)
#     pII = lexicographic(II, n+1)
#     A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(II))
#     A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
    
#     Avim1, Avi, Avip1 = A2_1 * v[δy⁻(II)], A2_2 * v[II], A4_2 * v[δy⁺(II)]
#     Auim1jm1, Auip1jm1 = A1_1 * u[δy⁻(II)], A3_1 * u[δx⁺(δy⁻(II))]
#     Auim1jp1, Auip1jp1 = A1_2 * u[II], A3_2 * u[δx⁺(II)]

#     Au1 = 0.5 * (Auim1jm1 + Auim1jp1)
#     Au2 = 0.5 * (Avim1 + Avi)
#     Au3 = 0.5 * (Auip1jm1 + Auip1jp1)
#     Au4 = 0.5 * (Avi + Avip1)

#     @inbounds O[pII,pII] = 0.5 * (Au3 - Au1 + Au4 - Au2)
#     @inbounds O[pII,pII+n+1] = 0.5 * Au3
#     @inbounds O[pII,pII-n-1] = -0.5 * Au1
#     @inbounds O[pII,pII+1] = 0.5 * Au4
#     @inbounds O[pII,pII-1] = -0.5 * Au2

#     @inbounds O[pII,pII] += -0.25 * (A4_2 - B2_2) * Dv1_y[δy⁺(II)]
#     @inbounds O[pII,pII] += -0.25 * (B2_2 - B2_1) * Dv1_y[II]
#     @inbounds O[pII,pII] += -0.25 * (B2_1 - A2_1) * Dv1_y[δy⁻(II)]

#     @inbounds O[pII,pII] += -0.25 * (A3_1 - B1_1) * Du_x[δy⁻(δx⁺(II))]
#     @inbounds O[pII,pII] += -0.25 * (B1_1 - A1_1) * Du_x[δy⁻(II)]
#     @inbounds O[pII,pII] += -0.25 * (A3_2 - B1_2) * Du_x[δx⁺(II)]
#     @inbounds O[pII,pII] += -0.25 * (B1_2 - A1_2) * Du_x[II]

#     # @inbounds B[pII] += -0.25 * Dv2_y[II] * (A4_2 - B2_2) * Dv1_y[δy⁺(II)]
#     # @inbounds B[pII] += -0.25 * Dv2_y[II] * (B2_2 - B2_1) * Dv1_y[II]
#     # @inbounds B[pII] += -0.25 * Dv2_y[II] * (B2_1 - A2_1) * Dv1_y[δy⁻(II)]

#     # @inbounds B[pII] += -0.25 * Dv2_x[II] * (A3_1 - B1_1) * Du_x[δy⁻(δx⁺(II))]
#     # @inbounds B[pII] += -0.25 * Dv2_x[II] * (B1_1 - A1_1) * Du_x[δy⁻(II)]
#     # @inbounds B[pII] += -0.25 * Dv2_x[II] * (A3_2 - B1_2) * Du_x[δx⁺(II)]
#     # @inbounds B[pII] += -0.25 * Dv2_x[II] * (B1_2 - A1_2) * Du_x[II]

#     @inbounds B[pII] += -0.25 * Dv1_y[II] * (A4_2 - B2_2) * Dv1_y[δy⁺(II)]
#     @inbounds B[pII] += -0.25 * Dv1_y[II] * (B2_2 - B2_1) * Dv1_y[II]
#     @inbounds B[pII] += -0.25 * Dv1_y[II] * (B2_1 - A2_1) * Dv1_y[δy⁻(II)]

#     @inbounds B[pII] += -0.25 * Dv1_x[II] * (A3_1 - B1_1) * Du_x[δy⁻(δx⁺(II))]
#     @inbounds B[pII] += -0.25 * Dv1_x[II] * (B1_1 - A1_1) * Du_x[δy⁻(II)]
#     @inbounds B[pII] += -0.25 * Dv1_x[II] * (A3_2 - B1_2) * Du_x[δx⁺(II)]
#     @inbounds B[pII] += -0.25 * Dv1_x[II] * (B1_2 - A1_2) * Du_x[II]
# end

# function vector_convection!(::Dirichlet, ::Type{GridFCy}, O, B, u, v, Du_x, Du_y, Dv1_x, Dv1_y, Dv2_x, Dv2_y, cap, n, BC, inside, b_left, b_bottom, b_right, b_top)
#     B .= 0.0
#     @inbounds @threads for II in inside
#         fill_inside_conv!(GridFCy, O, B, u, v, Du_x, Dv1_x, Dv1_y, cap, n, II)
#     end

#     @inbounds @threads for II in vcat(b_left, b_bottom[2:end-1], b_right, b_top[2:end-1])
#         pII = lexicographic(II, n+1)
#         @inbounds O[pII,pII] = 0.0
#     end
#     bnds = (b_left[2:end-1], b_bottom, b_right[2:end-1])
#     bc = ((Du_x, Dv1_x, Dv1_x), (Du_y, Dv1_y, Dv1_y), (Du_x, Dv1_x, Dv1_x))
#     for (bnd, (Du, Dv1, Dv2)) in zip(bnds, bc)
#         @inbounds for II in bnd
#             pII = lexicographic(II, n+1)
#             A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
            
#             Avi, Avip1 = A2_2 * v[II], A4_2 * v[δy⁺(II)]

#             Au4 = 0.5 * (Avi + Avip1)

#             @inbounds O[pII,pII] += 0.5 * Au4
#             @inbounds O[pII,pII+1] = 0.5 * Au4

#             @inbounds O[pII,pII] += -0.25 * (A4_2 - B2_2) * Dv1[δy⁺(II)]
#             @inbounds O[pII,pII] += -0.25 * (B2_2 - A2_2) * Dv1[II]

#             @inbounds O[pII,pII] += -0.25 * (A3_2 - B1_2) * Du[δx⁺(II)]
#             @inbounds O[pII,pII] += -0.25 * (B1_2 - A1_2) * Du[II]

#             @inbounds B[pII] += -0.25 * Dv2[II] * (A4_2 - B2_2) * Dv1[δy⁺(II)]
#             @inbounds B[pII] += -0.25 * Dv2[II] * (B2_2 - A2_2) * Dv1[II]

#             @inbounds B[pII] += -0.25 * Dv2[II] * (A3_2 - B1_2) * Du[δx⁺(II)]
#             @inbounds B[pII] += -0.25 * Dv2[II] * (B1_2 - A1_2) * Du[II]
#         end
#     end
#     bnds = (b_left[2:end-1], b_right[2:end-1], b_top)
#     bc = ((Du_x, Dv1_x, Dv1_x), (Du_x, Dv1_x, Dv1_x), (Du_y, Dv1_y, Dv1_y))
#     for (bnd, (Du, Dv1, Dv2)) in zip(bnds, bc)
#         @inbounds for II in bnd
#             pII = lexicographic(II, n+1)
#             A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(II))
            
#             Avim1, Avi = A2_1 * v[δy⁻(II)], A4_1 * v[II]

#             Au2 = 0.5 * (Avim1 + Avi)

#             @inbounds O[pII,pII] += -0.5 * Au2
#             @inbounds O[pII,pII-1] = -0.5 * Au2

#             @inbounds O[pII,pII] += -0.25 * (A4_1 - B2_1) * Dv1[II]
#             @inbounds O[pII,pII] += -0.25 * (B2_1 - A2_1) * Dv1[δy⁻(II)]

#             @inbounds O[pII,pII] += -0.25 * (A3_1 - B1_1) * Du[δy⁻(δx⁺(II))]
#             @inbounds O[pII,pII] += -0.25 * (B1_1 - A1_1) * Du[δy⁻(II)]

#             @inbounds B[pII] += -0.25 * Dv2[II] * (A4_1 - B2_1) * Dv1[II]
#             @inbounds B[pII] += -0.25 * Dv2[II] * (B2_1 - A2_1) * Dv1[δy⁻(II)]

#             @inbounds B[pII] += -0.25 * Dv2[II] * (A3_1 - B1_1) * Du[δy⁻(δx⁺(II))]
#             @inbounds B[pII] += -0.25 * Dv2[II] * (B1_1 - A1_1) * Du[δy⁻(II)]
#         end
#     end
#     @inbounds for II in b_left[2:end-1]
#         pII = lexicographic(II, n+1)
#         A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(II))
#         A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
        
#         Auip1jm1 = A3_1 * u[δx⁺(δy⁻(II))]
#         Auip1jp1 = A3_2 * u[δx⁺(II)]

#         Au3 = 0.5 * (Auip1jm1 + Auip1jp1)

#         @inbounds O[pII,pII] += 0.5 * Au3
#         @inbounds O[pII,pII+n+1] = 0.5 * Au3
#     end
#     @inbounds for II in b_right[2:end-1]
#         pII = lexicographic(II, n+1)
#         A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(II))
#         A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
        
#         Auim1jm1 = A1_1 * u[δy⁻(II)]
#         Auim1jp1 = A1_2 * u[II]

#         Au1 = 0.5 * (Auim1jm1 + Auim1jp1)

#         @inbounds O[pII,pII] += -0.5 * Au1
#         @inbounds O[pII,pII-n-1] = -0.5 * Au1
#     end
#     @inbounds for II in b_bottom[2:end]
#         pII = lexicographic(II, n+1)
#         A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
        
#         Auim1jp1 = A1_2 * u[II]

#         Au1 = 0.5 * Auim1jp1

#         if is_periodic(BC.bottom.t)
#             JJ = II + CartesianIndex(n, 0)
#             A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, JJ)
#             Auim1jm1 = A1_1 * u[JJ]
#             Au1 += 0.5 * Auim1jm1
#         end

#         @inbounds O[pII,pII] += -0.5 * Au1
#         @inbounds O[pII,pII-n-1] = -0.5 * Au1
#     end
#     @inbounds for II in b_bottom[1:end-1]
#         pII = lexicographic(II, n+1)
#         A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
        
#         Auip1jp1 = A3_2 * u[δx⁺(II)]

#         Au3 = 0.5 * Auip1jp1

#         if is_periodic(BC.bottom.t)
#             JJ = II + CartesianIndex(n, 0)
#             A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, JJ)
#             Auip1jm1 = A3_1 * u[δx⁺(JJ)]
#             Au3 += 0.5 * Auip1jm1
#         end

#         @inbounds O[pII,pII] += 0.5 * Au3
#         @inbounds O[pII,pII+n+1] = 0.5 * Au3
#     end
#     @inbounds for II in b_top[2:end]
#         pII = lexicographic(II, n+1)
#         A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(II))
        
#         Auim1jm1 = A1_1 * u[δy⁻(II)]

#         Au1 = 0.5 * Auim1jm1

#         if is_periodic(BC.top.t)
#             JJ = II + CartesianIndex(-n, 0)
#             A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ)
#             Auim1jp1 = A1_2 * u[JJ]
#             Au1 += 0.5 * Auim1jp1
#         end

#         @inbounds O[pII,pII] += -0.5 * Au1
#         @inbounds O[pII,pII-n-1] = -0.5 * Au1
#     end
#     @inbounds for II in b_top[1:end-1]
#         pII = lexicographic(II, n+1)
#         A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(II))
        
#         Auip1jm1 = A3_1 * u[δx⁺(δy⁻(II))]

#         Au3 = 0.5 * Auip1jm1

#         if is_periodic(BC.top.t)
#             JJ = II + CartesianIndex(n, 0)
#             A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ)
#             Auip1jp1 = A3_2 * u[δx⁺(JJ)]
#             Au3 += 0.5 * Auip1jp1
#         end

#         @inbounds O[pII,pII] += 0.5 * Au3
#         @inbounds O[pII,pII+n+1] = 0.5 * Au3
#     end

#     if is_periodic(BC.bottom.t) && is_periodic(BC.top.t)
#         (Du, Dv1, Dv2) = (Du_y, Dv1_y, Dv1_y)
#         @inbounds for (II, JJ) in zip(b_bottom, b_top)
#             pII = lexicographic(II, n+1)
#             pJJ = lexicographic(JJ, n+1)
#             A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(JJ))
            
#             Avim1, Avi = A2_1 * v[JJ], A4_1 * v[II]

#             Au2 = 0.5 * (Avim1 + Avi)

#             @inbounds O[pII,pII] += -0.5 * Au2
#             @inbounds O[pII,pJJ] = -0.5 * Au2

#             @inbounds O[pII,pII] += -0.25 * (A4_1 - B2_1) * Dv1[II]
#             @inbounds O[pII,pII] += -0.25 * (B2_1 - A2_1) * Dv1[JJ]

#             @inbounds O[pII,pII] += -0.25 * (A3_1 - B1_1) * Du[δx⁺(JJ)]
#             @inbounds O[pII,pII] += -0.25 * (B1_1 - A1_1) * Du[JJ]

#             @inbounds B[pII] += -0.25 * Dv2[II] * (A4_1 - B2_1) * Dv1[II]
#             @inbounds B[pII] += -0.25 * Dv2[II] * (B2_1 - A2_1) * Dv1[JJ]

#             @inbounds B[pII] += -0.25 * Dv2[II] * (A3_1 - B1_1) * Du[δx⁺(JJ)]
#             @inbounds B[pII] += -0.25 * Dv2[II] * (B1_1 - A1_1) * Du[JJ]
#         end
#         (Du, Dv1, Dv2) = (Du_y, Dv1_y, Dv1_y)
#         @inbounds for (II, JJ) in zip(b_top, b_bottom)
#             pII = lexicographic(II, n+1)
#             pJJ = lexicographic(JJ, n+1)
#             A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ)
            
#             Avi, Avip1 = A2_2 * v[II], A4_2 * v[JJ]

#             Au4 = 0.5 * (Avi + Avip1)

#             @inbounds O[pII,pII] += 0.5 * Au4
#             @inbounds O[pII,pJJ] = 0.5 * Au4

#             @inbounds O[pII,pII] += -0.25 * (A4_2 - B2_2) * Dv1[JJ]
#             @inbounds O[pII,pII] += -0.25 * (B2_2 - A2_2) * Dv1[II]

#             @inbounds O[pII,pII] += -0.25 * (A3_2 - B1_2) * Du[δx⁺(II)]
#             @inbounds O[pII,pII] += -0.25 * (B1_2 - A1_2) * Du[II]

#             @inbounds B[pII] += -0.25 * Dv2[II] * (A4_2 - B2_2) * Dv1[JJ]
#             @inbounds B[pII] += -0.25 * Dv2[II] * (B2_2 - A2_2) * Dv1[II]

#             @inbounds B[pII] += -0.25 * Dv2[II] * (A3_2 - B1_2) * Du[δx⁺(II)]
#             @inbounds B[pII] += -0.25 * Dv2[II] * (B1_2 - A1_2) * Du[II]
#         end
#     end
#     if is_periodic(BC.left.t) && is_periodic(BC.right.t)
#         @inbounds for (II,JJ) in zip(b_left[2:end-1], b_right[2:end-1])
#             pII = lexicographic(II, n+1)
#             pJJ = lexicographic(JJ, n+1)
#             A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(II))
#             A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
            
#             Auim1jm1 = A1_1 * u[δy⁻(II)]
#             Auim1jp1 = A1_2 * u[II]
    
#             Au1 = 0.5 * (Auim1jm1 + Auim1jp1)
    
#             @inbounds O[pII,pII] += -0.5 * Au1
#             @inbounds O[pII,pJJ] = -0.5 * Au1
#         end
#         @inbounds for (II,JJ) in zip(b_right[2:end-1], b_left[2:end-1])
#             pII = lexicographic(II, n+1)
#             pJJ = lexicographic(JJ, n+1)
#             A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(II))
#             A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
            
#             Auip1jm1 = A3_1 * u[δx⁺(δy⁻(II))]
#             Auip1jp1 = A3_2 * u[δx⁺(II)]
    
#             Au3 = 0.5 * (Auip1jm1 + Auip1jp1)
    
#             @inbounds O[pII,pII] += 0.5 * Au3
#             @inbounds O[pII,pJJ] = 0.5 * Au3
#         end

#         ii = b_bottom[1]
#         pii = lexicographic(ii, n+1)
#         A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, ii)
        
#         Auim1jp1 = A1_2 * u[ii]

#         Au1 = 0.5 * Auim1jp1

#         if is_periodic(BC.bottom.t)
#             JJ = ii + CartesianIndex(n, 0)
#             A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, JJ)
#             Auim1jm1 = A1_1 * u[JJ]
#             Au1 += 0.5 * Auim1jm1
#         end

#         JJ = ii + CartesianIndex(0, n-1)
#         pJJ = lexicographic(JJ, n+1)
#         @inbounds O[pii,pii] += -0.5 * Au1
#         @inbounds O[pii,pJJ] = -0.5 * Au1
        
#         ii = b_bottom[end]
#         pii = lexicographic(ii, n+1)
#         A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, ii)
        
#         Auip1jp1 = A3_2 * u[δx⁺(ii)]

#         Au3 = 0.5 * Auip1jp1

#         if is_periodic(BC.bottom.t)
#             JJ = II + CartesianIndex(n, 0)
#             A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, JJ)
#             Auip1jm1 = A3_1 * u[δx⁺(JJ)]
#             Au3 += 0.5 * Auip1jm1
#         end

#         JJ = ii + CartesianIndex(0, -n+1)
#         pJJ = lexicographic(JJ, n+1)
#         @inbounds O[pii,pii] += 0.5 * Au3
#         @inbounds O[pii,pJJ] = 0.5 * Au3
        
#         ii = b_top[1]
#         pii = lexicographic(ii, n+1)
#         A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(ii))
        
#         Auim1jm1 = A1_1 * u[δy⁻(ii)]

#         Au1 = 0.5 * Auim1jm1

#         if is_periodic(BC.top.t)
#             JJ = II + CartesianIndex(-n, 0)
#             A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ)
#             Auim1jp1 = A1_2 * u[JJ]
#             Au1 += 0.5 * Auim1jp1
#         end

#         JJ = ii + CartesianIndex(0, n-1)
#         pJJ = lexicographic(JJ, n+1)
#         @inbounds O[pii,pii] += -0.5 * Au1
#         @inbounds O[pii,pJJ] = -0.5 * Au1
        
#         ii = b_top[end]
#         pii = lexicographic(ii, n+1)
#         A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(ii))
        
#         Auip1jm1 = A3_1 * u[δx⁺(δy⁻(ii))]

#         Au3 = 0.5 * Auip1jm1

#         if is_periodic(BC.top.t)
#             JJ = II + CartesianIndex(n, 0)
#             A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ)
#             Auip1jp1 = A3_2 * u[δx⁺(JJ)]
#             Au3 += 0.5 * Auip1jp1
#         end

#         JJ = ii + CartesianIndex(0, -n+1)
#         pJJ = lexicographic(JJ, n+1)
#         @inbounds O[pii,pii] += 0.5 * Au3
#         @inbounds O[pii,pJJ] = 0.5 * Au3
#     end

#     return nothing
# end
