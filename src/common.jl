"""
    lexicographic(II, n)

Returns the index of type `Int` corresponding to the 1D view of a `n × n` array at the CartesianIndex `II`.

"""
@inline lexicographic(II, n) = muladd(n, II[2]-1, II[1])
@inline within_bounds(i, grid) = 1 <= i <= grid.nx*grid.ny
@inline within_bounds(i::CartesianIndex, grid) = 1 <= i[1] <= grid.ny && 1 <= i[2] <= grid.nx

const newaxis = [CartesianIndex()]

@inline δy⁺(II) = CartesianIndex(II[1]+1, II[2])
@inline δy⁻(II) = CartesianIndex(II[1]-1, II[2])
@inline δx⁺(II) = CartesianIndex(II[1], II[2]+1)
@inline δx⁻(II) = CartesianIndex(II[1], II[2]-1)

@inline δy⁺(II, ny, per) = per && II[1]==ny ? CartesianIndex(1, II[2]) : δy⁺(II)
@inline δy⁻(II, ny, per) = per && II[1]==1 ? CartesianIndex(ny, II[2]) : δy⁻(II)
@inline δx⁺(II, nx, per) = per && II[2]==nx ? CartesianIndex(II[1], 1) : δx⁺(II)
@inline δx⁻(II, nx, per) = per && II[2]==1 ? CartesianIndex(II[1], nx) : δx⁻(II)

# Temporary function to get a certain field from a vector with 
# multiple fields. To be removed when working with decomposed 
# vectors directly
veci(a, g::G, p::Int) where {G<:Grid} = @view a[g.ny*g.nx*(p-1)+1:g.ny*g.nx*p]

function veci(a, g::Vector{G}, p::Int) where {G<:Grid}
    c0 = 1
    c1 = g[1].ny * g[1].nx
    for i = 2:p
        ii = (p - 1) ÷ 2 + 1
        c0 += g[ii].ny * g[ii].nx
        c1 += g[ii].ny * g[ii].nx
    end

    return @view a[c0:c1]
end

@inline opposite(α) = ifelse(α >= 0. ,α - π, α + π)

@inline distance(A::Point, B::Point, dx=1.0, dy=1.0) =
    sqrt(((A.x-B.x)*dx)^2 + ((A.y-B.y)*dy)^2)

@inline midpoint(A::Point, B::Point) = Point((A.x + B.x)/2 - 0.5, (A.y + B.y)/2 - 0.5)

@inline indomain(A::Point{T}, x::Vector{T}, y::Vector{T}) where {T} =
    ifelse(A.x >= x[1] && A.x <= x[end] &&
           A.y >= y[1] && A.y <= y[end], true, false)

@inline <(A::Point{T}, L::T) where {T} = ifelse(abs(A.x) <= L/2 && abs(A.y) <= L/2, true, false)
@inline isnan(A::Point{T}) where {T} = ifelse(isnan(A.x) || isnan(A.y), true, false)
@inline +(A::Point{T}, B::Point{T}) where {T} = Point(A.x + B.x, A.y + B.y)
@inline -(A::Point{T}, B::Point{T}) where {T} = Point(A.x - B.x, A.y - B.y)
@inline *(A::Point{T}, m::N) where {T, N} = Point(m*A.x, m*A.y)
@inline *(m::N, A::Point{T}) where {T, N} = Point(m*A.x, m*A.y)
@inline *(A::Point{T}, B::Point{T}) where {T} = Point(A.x*B.x, A.y*B.y)
@inline abs(A::Point) = Point(abs(A.x), abs(A.y))

@inline mysign(a,b) = a/sqrt(a^2 + b^2)
@inline mysign(a) = ifelse(a >= 0, 1, -1)

@inline ⁺(a) = max(a,0)
@inline ⁻(a) = min(a,0)

@inline c∇x(u, II) = @inbounds u[δx⁺(II)] - u[δx⁻(II)]
@inline c∇y(u, II) = @inbounds u[δy⁺(II)] - u[δy⁻(II)]
@inline c∇x(u, II, nx, per) = @inbounds u[δx⁺(II, nx, per)] - u[δx⁻(II, nx, per)]
@inline c∇y(u, II, ny, per) = @inbounds u[δy⁺(II, ny, per)] - u[δy⁻(II, ny, per)]
@inline c∇x(u, II, h) = @inbounds (u[δx⁺(II)] - u[δx⁻(II)])/2h
@inline c∇y(u, II, h) = @inbounds (u[δy⁺(II)] - u[δy⁻(II)])/2h
@inline c∇x(u, II, h, nx, per) = @inbounds (u[δx⁺(II, nx, per)] - u[δx⁻(II, nx, per)])/h
@inline c∇y(u, II, h, ny, per) = @inbounds (u[δy⁺(II, ny, per)] - u[δy⁻(II, ny, per)])/h

@inline ∇x⁺(u, II) = @inbounds u[δx⁺(II)] - u[II]
@inline ∇y⁺(u, II) = @inbounds u[δy⁺(II)] - u[II]
@inline ∇x⁻(u, II) = @inbounds u[δx⁻(II)] - u[II]
@inline ∇y⁻(u, II) = @inbounds u[δy⁻(II)] - u[II]

@inline ∇x⁺(u, II, nx, per) = @inbounds u[δx⁺(II, nx, per)] - u[II]
@inline ∇y⁺(u, II, ny, per) = @inbounds u[δy⁺(II, ny, per)] - u[II]
@inline ∇x⁻(u, II, nx, per) = @inbounds u[δx⁻(II, nx, per)] - u[II]
@inline ∇y⁻(u, II, ny, per) = @inbounds u[δy⁻(II, ny, per)] - u[II]

@inline normal(u, II) = @SVector [mysign(c∇x(u, II), c∇y(u, II)), mysign(c∇y(u, II), c∇x(u, II))]

@inline minmod(a, b) = ifelse(a*b <= 0, 0.0, myabs(a,b))
@inline myabs(a, b) = ifelse(abs(a) < abs(b), a, b)

@inline Dxx(u, II, dx) = @inbounds (2.0 * (u[δx⁺(II)] - u[II]) / (dx[δx⁺(II)] + dx[II]) -
                                    2.0 * (u[II] - u[δx⁻(II)]) / (dx[δx⁻(II)] + dx[II])) / dx[II]
@inline Dyy(u, II, dy) = @inbounds (2.0 * (u[δy⁺(II)] - u[II]) / (dy[δy⁺(II)] + dy[II]) -
                                    2.0 * (u[II] - u[δy⁻(II)]) / (dy[δy⁻(II)] + dy[II])) / dy[II]

@inline Dxx(u, II, dx, nx, per) = @inbounds (2.0 * (u[δx⁺(II, nx, per)] - u[II]) / (dx[δx⁺(II, nx, per)] + dx[II]) -
                                    2.0 * (u[II] - u[δx⁻(II, nx, per)]) / (dx[δx⁻(II, nx, per)] + dx[II])) / dx[II]
@inline Dyy(u, II, dy, ny, per) = @inbounds (2.0 * (u[δy⁺(II, ny, per)] - u[II]) / (dy[δy⁺(II, ny, per)] + dy[II]) -
                                    2.0 * (u[II] - u[δy⁻(II, ny, per)]) / (dy[δy⁻(II, ny, per)] + dy[II])) / dy[II]

@inline Dxx(u, II) = @inbounds u[δx⁺(II)] - 2u[II] + u[δx⁻(II)]
@inline Dyy(u, II) = @inbounds u[δy⁺(II)] - 2u[II] + u[δy⁻(II)]

@inline Dxx(u, II, nx, per) = @inbounds u[δx⁺(II, nx, per)] - 2u[II] + u[δx⁻(II, nx, per)]
@inline Dyy(u, II, ny, per) = @inbounds u[δy⁺(II, ny, per)] - 2u[II] + u[δy⁻(II, ny, per)]

@inline Dxy(u, II, h) = @inbounds (u[δx⁻(δy⁻(II))] + u[δx⁺(δy⁺(II))] - u[δx⁻(δy⁺(II))] - u[δx⁺(δy⁻(II))])/4h^2

function smekerka_curvature(u, II, h)
    gx = c∇x(u, II, h)
    gy = c∇y(u, II, h)
    gxx = Dxx(u, II, h)
    gyy = Dyy(u, II, h)
    gxy = Dxy(u, II, h)
    k = (gxx + gyy)/(gx^2 + gy^2)^0.5
    return k
end

@inline in_bounds(a, n) = ifelse(a > 2 && a < n-1, true, false)
@inline in_bounds(a, n, per) = per ? ifelse(a > 1 && a < n, true, false) : in_bounds(a, n)

@inline function static_stencil(a, II::CartesianIndex)
   return @inbounds SA_F64[a[II.I[1]-1, II.I[2]-1] a[II.I[1]-1, II.I[2]] a[II.I[1]-1, II.I[2]+1];
               a[II.I[1], II.I[2]-1] a[II.I[1], II.I[2]] a[II.I[1], II.I[2]+1];
               a[II.I[1]+1, II.I[2]-1] a[II.I[1]+1, II.I[2]] a[II.I[1]+1, II.I[2]+1]]
end

@inline function static_stencil(a, II::CartesianIndex, nx, ny, per_x, per_y)
    return @inbounds SA_F64[a[δx⁻(δy⁻(II, ny, per_y), nx, per_x)] a[δy⁻(II, ny, per_y)] a[δx⁺(δy⁻(II, ny, per_y), nx, per_x)];
                a[δx⁻(II, nx, per_x)] a[II] a[δx⁺(II, nx, per_x)];
                a[δx⁻(δy⁺(II, ny, per_y), nx, per_x)] a[δy⁺(II, ny, per_y)] a[δx⁺(δy⁺(II, ny, per_y), nx, per_x)]]
 end

@inline function B_BT(II, x, y, f=x->x)
    B = inv((@SMatrix [((x[f(δx⁻(II))]-x[II])/(x[f(δx⁺(II))] - x[f(δx⁻(II))]))^2 (x[f(δx⁻(II))]-x[II])/(x[f(δx⁺(II))] - x[f(δx⁻(II))]) 1.0;
                      ((x[f(II)]-x[II])/(x[f(δx⁺(II))] - x[f(δx⁻(II))]))^2 (x[f(II)]-x[II])/(x[f(δx⁺(II))] - x[f(δx⁻(II))]) 1.0;
                      ((x[f(δx⁺(II))]-x[II])/(x[f(δx⁺(II))] - x[f(δx⁻(II))]))^2 (x[f(δx⁺(II))]-x[II])/(x[f(δx⁺(II))] - x[f(δx⁻(II))]) 1.0]))

    BT = inv((@SMatrix [((y[f(δy⁻(II))]-y[II])/(y[f(δy⁺(II))] - y[f(δy⁻(II))]))^2 ((y[f(II)]-y[II])/(y[f(δy⁺(II))] - y[f(δy⁻(II))]))^2 ((y[f(δy⁺(II))]-y[II])/(y[f(δy⁺(II))] - y[f(δy⁻(II))]))^2;
                        (y[f(δy⁻(II))]-y[II])/(y[f(δy⁺(II))] - y[f(δy⁻(II))]) (y[f(II)]-y[II])/(y[f(δy⁺(II))] - y[f(δy⁻(II))]) (y[f(δy⁺(II))]-y[II])/(y[f(δy⁺(II))] - y[f(δy⁻(II))]);
                        1.0 1.0 1.0]))
    
    return B, BT
end

@inline function B_BT(II::CartesianIndex, grid::G, per_x, per_y) where {G<:Grid}
    @unpack x, y, nx, ny, dx, dy, α = grid

    _dx = dx[II]
    _dy = dy[II]

    dx⁺ = 0.5 * (dx[II] + dx[δx⁺(II, nx, per_x)])
    dx⁻ = 0.5 * (dx[δx⁻(II, nx, per_x)] + dx[II])

    dy⁺ = 0.5 * (dy[II] + dy[δy⁺(II, ny, per_y)])
    dy⁻ = 0.5 * (dy[δy⁻(II, ny, per_y)] + dy[II])

    B = inv((@SMatrix [(dx⁻/_dx)^2 -dx⁻/_dx 1.0;
        0.0 0.0 1.0;
        (dx⁺/_dx)^2 dx⁺/_dx 1.0]))

    BT = inv((@SMatrix [(dy⁻/_dy)^2 0.0 (dy⁺/_dy)^2;
        -dy⁻/_dy 0.0 dy⁺/_dy;
        1.0 1.0 1.0]))

    return B, BT
end

@inline function B_BT(II::CartesianIndex, grid::G) where {G<:Grid}
    B = inv(@SMatrix [0.25 -0.5 1.0;
                      0.0 0.0 1.0;
                      0.25 0.5 1.0])

    BT = inv(@SMatrix [0.25 -0.5 1.0;
                        0.0 0.0 1.0;
                        0.25 0.5 1.0])
    
    return B, BT
end

@inline big_static_stencil(a, II::CartesianIndex) = @inbounds @SMatrix [a[II.I[1]-2, II.I[2]-2] a[II.I[1]-2, II.I[2]-1] a[II.I[1]-2, II.I[2]] a[II.I[1]-2, II.I[2]+1] a[II.I[1]-2, II.I[2]+2];
   a[II.I[1]-1, II.I[2]-2] a[II.I[1]-1, II.I[2]-1] a[II.I[1]-1, II.I[2]] a[II.I[1]-1, II.I[2]+1] a[II.I[1]-1, II.I[2]+2];
   a[II.I[1], II.I[2]-2] a[II.I[1], II.I[2]-1] a[II.I[1], II.I[2]] a[II.I[1], II.I[2]+1] a[II.I[1], II.I[2]+2];
   a[II.I[1]+1, II.I[2]-2] a[II.I[1]+1, II.I[2]-1] a[II.I[1]+1, II.I[2]] a[II.I[1]+1, II.I[2]+1] a[II.I[1]+1, II.I[2]+2];
   a[II.I[1]+2, II.I[2]-2] a[II.I[1]+2, II.I[2]-1] a[II.I[1]+2, II.I[2]] a[II.I[1]+2, II.I[2]+1] a[II.I[1]+2, II.I[2]+2]]

@inline function mean_curvature(a::StaticArray, h)
    f = @SVector [(a[2,3]-a[2,1])/2h, (a[3,2]-a[1,2])/2h]
    s = @SVector [(a[2,3]-2a[2,2]+a[2,1])/h^2, (a[3,2]-2a[2,2]+a[1,2])/h^2]
    m = @SVector [(a[3,3]-a[1,3]-a[3,1]+a[1,1])/4h^2]
    return (s[1]*(f[2])^2 - 2*f[1]*f[2]*m[1] + s[2]*(f[1])^2)/(f[1]^2 + f[2]^2)^1.5
end

@inline function mean_curvature_indices(a, II, h)
    f = @SVector [(a[δx⁺(II)]-a[δx⁻(II)])/2h, (a[δy⁺(II)]-a[δy⁻(II)])/2h]
    s = @SVector [(a[δx⁻(II)]-2a[II]+a[δx⁺(II)])/h^2, (a[δy⁻(II)]-2a[II]+a[δy⁺(II)])/h^2]
    m = @SVector [(a[δx⁻(δy⁻(II))]+a[δx⁺(δy⁺(II))]-a[δx⁻(δy⁺(II))]+a[δx⁺(δy⁻(II))])/4h^2]
    return (s[1]*(f[2])^2 - 2*f[1]*f[2]*m[1] + s[2]*(f[1])^2)/(1e-10 + f[1]^2 + f[2]^2)^1.5
end

function mean_curvature_interpolated(u, II, h, B, BT, mid)
    tmp = zeros(3,3)
    @threads for i in -1:1
        for j in -1:1
            JJ = CartesianIndex(II[1]+i, II[2]+j)
            st = static_stencil(u, JJ)
            tmp[2+i, 2+j] = mean_curvature(st, h)
        end
    end
    itp = B*tmp*BT
    a = biquadratic(itp, mid.x, mid.y)
    return min(max(mean(tmp),-1/h), 1/h)
end


@inline parabola_fit_curvature(itp, mid_point, dx, dy) =
    @inbounds 2*itp[1,3]/((1+(2*itp[1,3]*mid_point.x/dx + itp[2,3])^2)^1.5)/(dx)^2 +
              2*itp[3,1]/((1+(2*itp[3,1]*mid_point.y/dy + itp[3,2])^2)^1.5)/(dy)^2

@inline function linear_fit(a, b, x)
    a = (a - b)/sign(x+eps(0.1))
    return a*x + b
end

@inline @. parabola_fit(x, p) = p[1]*x^2 + p[2]*x + p[3]

function get_NB_width_indices_base(n)
    x = [CartesianIndex(i,0) for i in -n:n]
    y = [CartesianIndex(0,j) for j in -n:n]
    return vcat(x,y)
end

@inline function normf(field, pos, cap, h)
    AVG = 0.
    RMS = 0.
    VOLUME = 0.
    MAX = 0.
    dv = 0.
    @inbounds for II in pos
        v = abs(field[II])
        dv = cap[II]*h^2
        if dv > 0.
            VOLUME += dv
            AVG += dv*v
            RMS += dv*v^2
            if (v > MAX) MAX = v end
        end
    end
    return AVG/VOLUME, sqrt(RMS/VOLUME), MAX
end

@inline function normf(field, pos, h)
    AVG = 0.
    RMS = 0.
    VOLUME = 0.
    MAX = 0.
    @inbounds for II in pos
        v = abs(field[II])
        VOLUME += h^2
        AVG += (h^2)*v
        RMS += (h^2)*v^2
        if (v > MAX) MAX = v end
    end
    return AVG/VOLUME, sqrt(RMS/VOLUME), MAX
end

@inline function norma(field, pos, cap, h)
    AVG = 0.
    RMS = 0.
    VOLUME = 0.
    MAX = 0.
    dv = 0.
    @inbounds for II in pos
        v = abs(field[II])
        dv = cap[II]*h^2
        # dv = h^2
        if dv > 0.
            VOLUME += dv
            AVG += dv*v
            RMS += dv*v^2
            if (v > MAX) MAX = v end
        end
    end
    return AVG/VOLUME, sqrt(RMS), MAX
end

@inline function Richardson_extrapolation(e, r)
    p = log(abs(e[end-2] - e[end-1])/abs(e[end-1] - e[end]))/log(r)
    ext = e[end] + (e[end] - e[end-1])/(2^p - 1)
    return abs.(e .- ext)
end

function fit_order(x, y)
    coeffs = fit(log.(x), log.(y), 1).coeffs
    return exp(coeffs[1]), -coeffs[2]
end

@inline function arc_length2(POS, ind)
    d_tot = 0.
    d_ = 0.
    ind_ = copy(ind)
    for i in length(ind):-1:1
       d = 1e300
       pop!(ind_)
       for j in axes(ind_,1)
           d_ = distance(abs(POS[ind[i]].pos), abs(POS[ind_[j]].pos))
           if d_ < d
               d = d_
           end
       end
       d_tot += d_
    end
    return d_tot
end

function find_2closest_points(POS, ind, II)
    closest_indices = Vector{Int64}(undef, 2)
    closest_cartesian_indices = Vector{CartesianIndex{2}}(undef, 2)
    all_indices = copy(ind)
    for  i = 1:2
        d = 1e300
        for JJ in axes(all_indices,1)
            d_ = distance(POS[II].pos, POS[all_indices[JJ]].pos)
            if 0. < d_ < d
                d = d_
                closest_indices[i] = JJ
            end
        end
        closest_cartesian_indices[i] = all_indices[closest_indices[i]]
        if i == 1 deleteat!(all_indices, closest_indices[1]) end
    end
    return closest_cartesian_indices
end

function monitor(header, history, it)
    println("****** ", header, " ******")
    if it > 0
        println("res[0] = ", first(history))
        println("res[", it, "] = ", history[it+1]) 
    else
        println("x0 is a solution")
    end
end

@inline is_dirichlet(::Dirichlet) = true
@inline is_dirichlet(::BoundaryCondition) = false

@inline is_neumann(::Neumann) = true
@inline is_neumann(::BoundaryCondition) = false

@inline is_robin(::Robin) = true
@inline is_robin(::BoundaryCondition) = false

@inline is_periodic(::Periodic) = true
@inline is_periodic(::BoundaryCondition) = false

@inline dirichlet(target, Δ, λ, val) = muladd(2.0, val, -target)
@inline neumann(target, Δ, λ, val) = muladd(Δ, val, target)
@inline robin(target, Δ, λ, val) = muladd((2*λ-Δ)/(2*λ+Δ), target, (2*Δ*val)/(Δ+2*λ))
@inline periodic(target, Δ, λ, val) = target

function bcs!(field, BC::Boundary{B, N, T, T}, Δ) where {B, N, T}
    @inbounds for KK in axes(field,3)
        for (II,JJ) in zip(BC.ind[1], BC.ind[2])
            field[II,KK] = BC.f(field[JJ,KK], Δ, BC.λ, BC.val)
        end
    end
    return nothing
end

function bcs!(field, BC::Boundary{B, N, T, Vector{T}}, Δ) where {B, N, T}
    @inbounds for KK in axes(field,3)
        for (II,JJ,LL) in zip(BC.ind[1], BC.ind[2], axes(field,1))
            field[II,KK] = BC.f(field[JJ,KK], Δ, BC.λ, BC.val[LL])
        end
    end
    return nothing
end

function get_fresh_cells!(grid, geo, Mm1, indices)
    @inbounds @threads for II in indices
        pII = lexicographic(II, grid.ny)
        if Mm1.diag[pII]/(grid.dx[II]*grid.dy[II]) < 1e-8 && geo.cap[II,5] > 1e-8
            geo.fresh[II] = true
        end
    end
    return nothing
end

function kill_dead_cells!(T::Matrix, grid, geo)
    @unpack ind = grid

    @inbounds @threads for II in ind.all_indices
        if geo.cap[II,5] < 1e-12
            T[II] = 0.
        end
    end
end

function kill_dead_cells!(T::Vector, grid, geo)
    @unpack ny, ind = grid

    @inbounds @threads for II in ind.all_indices
        pII = lexicographic(II, ny)
        if geo.cap[II,5] < 1e-12
            T[pII] = 0.
        end
    end
end

function kill_dead_cells!(S::SubArray{T,N,P,I,L}, grid, geo) where {T,N,P<:Vector{T},I,L}
    @unpack ny, ind = grid

    @inbounds @threads for II in ind.all_indices
        pII = lexicographic(II, ny)
        if geo.cap[II,5] < 1e-12
            S[pII] = 0.
        end
    end
end

function init_borders!(T::Matrix, BC_T, val=0.0)
    if is_dirichlet(BC_T.left.t)
        T[2:end-1,1] .= BC_T.left.val
    elseif is_periodic(BC_T.left.t)
        T[2:end-1,1] .= val
    end
    if is_dirichlet(BC_T.bottom.t)
        T[1,:] .= BC_T.bottom.val
    elseif is_periodic(BC_T.bottom.t)
        T[1,:] .= val
    end
    if is_dirichlet(BC_T.right.t)
        T[2:end-1,end] .= BC_T.right.val
    elseif is_periodic(BC_T.right.t)
        T[2:end-1,end] .= val
    end
    if is_dirichlet(BC_T.top.t)
        T[end,:] .= BC_T.top.val
    elseif is_periodic(BC_T.top.t)
        T[end,:] .= val
    end
end

"""
    export_all()

Exports every function of an included file.

* @args Noting : nothing

## Notes

Do not use this, use `@export` instead!

"""
function export_all()
    for n in names(@__MODULE__; all=true)
        if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
            @eval export $n
        end
    end
end

function mat_assign!(mat1, mat2)
    mat1.nzval .= 0.
    rows = rowvals(mat2)
    m, n = size(mat2)
    @inbounds @threads for j = 1:n
        for i in nzrange(mat2, j)
            @inbounds row = rows[i]
            mat1[row,j] = mat2[row,j]
        end
    end

    return nothing
end

function mat_assign_T!(mat1, mat2)
    mat1.nzval .= 0.
    rows = rowvals(mat1)
    m, n = size(mat1)
    @inbounds @threads for j = 1:n
        for i in nzrange(mat1, j)
            @inbounds row = rows[i]
            mat1[row,j] = mat2[row,j]
        end
    end

    return nothing
end

function mat_op!(mat1, mat2, op)
    mat1.nzval .= 0.
    rows = rowvals(mat2)
    vals = nonzeros(mat2)
    m, n = size(mat2)
    @inbounds @threads for j = 1:n
        for i in nzrange(mat2, j)
            @inbounds row = rows[i]
            @inbounds val = vals[i]
            @inbounds mat1[row,j] = op(val)
        end
    end

    return nothing
end

function mat_T_op!(mat1, mat2, op)
    mat1.nzval .= 0.
    rows = rowvals(mat2)
    vals = nonzeros(mat2)
    m, n = size(mat2)
    @inbounds @threads for j = 1:n
        for i in nzrange(mat2, j)
            @inbounds row = rows[i]
            @inbounds val = vals[i]
            @inbounds mat1[j,row] = op(val)
        end
    end

    return nothing
end

# Operations on sparse arrays that can be parallelized
for op ∈ (:*, :+, :-)
    @eval begin
        function $op(A::AbstractSparseMatrix{Tv,Ti}, B::Tv) where {Tv<:Number,Ti}
            C = SparseMatrixCSC{Tv,Ti}(A.m, A.n, A.colptr, A.rowval, zeros(length(A.nzval)))
            @inbounds @threads for col in 1:size(A, 2)
                for j in nzrange(A, col)
                    C.nzval[j] = $op(A.nzval[j], B)
                end
            end
            C
        end
        function $op(B::Tv, A::AbstractSparseMatrix{Tv,Ti}) where {Tv<:Number,Ti}
            C = SparseMatrixCSC{Tv,Ti}(A.m, A.n, A.colptr, A.rowval, zeros(length(A.nzval)))
            @inbounds @threads for col in 1:size(A, 2)
                for j in nzrange(A, col)
                    C.nzval[j] = $op(B, A.nzval[j])
                end
            end
            C
        end

        # function $op(A::AbstractSparseMatrix{Tv,Ti}, B::AbstractSparseMatrix{Tv,Ti}) where {Tv<:Number,Ti}
        #     C = SparseMatrixCSC{Tv,Ti}(A.m, A.n, A.colptr, A.rowval, zeros(length(A.nzval)))
        #     @inbounds @threads for col in 1:size(A, 2)
        #         for j in nzrange(A, col)
        #             C.nzval[j] = $op(A.nzval[j], B.nzval[j])
        #         end
        #     end
        #     C
        # end
        # function $op(A::AbstractSparseMatrix{Tv,Ti}, B::Diagonal{Tv,Vector{Tv}}) where {Tv<:Number,Ti}
        #     C = SparseMatrixCSC{Tv,Ti}(A.m, A.n, A.colptr, A.rowval, copy(A.nzval))
        #     b = B.diag
        #     nzv = nonzeros(A)
        #     rv = rowvals(A)
        #     @inbounds @threads for col in 1:size(A, 2)
        #         nz = nzrange(A, col)
        #         j = nz[findfirst(x->rv[x]==col, collect(nz))]
        #         C.nzval[j] = $op(A.nzval[j], b[col])
        #     end
        #     C
        # end
    end
end
function (-)(B::Diagonal{Tv,Vector{Tv}}, A::AbstractSparseMatrix{Tv,Ti}) where {Tv<:Number,Ti}
    C = SparseMatrixCSC{Tv,Ti}(A.m, A.n, A.colptr, A.rowval, -A.nzval)
    b = B.diag
    nzv = nonzeros(A)
    rv = rowvals(A)
    @inbounds @threads for col in 1:size(A, 2)
        nz = nzrange(A, col)
        j = nz[findfirst(x->rv[x]==col, collect(nz))]
        C.nzval[j] += b[col]
    end
    C
end

function mytime_print(elapsedtime, gctime=0)
    timestr = Base.Ryu.writefixed(Float64(elapsedtime/1e9), 6)
    str = sprint() do io
        print(io, length(timestr) < 10 ? (" "^(10 - length(timestr))) : "")
        print(io, timestr, " seconds")
        parens = gctime > 0
        parens && print(io, " (")
        if gctime > 0
            print(io, Base.Ryu.writefixed(Float64(gctime/1e9), 6), " gc time")
        end
        parens && print(io, " )")
    end
    println(str)
    nothing
end

# macro that computes time removing the
# garbage collector time
macro mytime(ex)
    quote
        Base.Experimental.@force_compile
        local stats = Base.gc_num()
        local elapsedtime = Base.time_ns()
        local val = $(esc(ex))
        elapsedtime = Base.time_ns() - elapsedtime
        local diff = Base.GC_Diff(Base.gc_num(), stats)
        local t = elapsedtime - diff.total_time
        mytime_print(t, diff.total_time)
        (time=t/1e9, gctime=diff.total_time/1e9)
    end
end
