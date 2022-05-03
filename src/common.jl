"""
    lexicographic(II, n)

Returns the index of type `Int` corresponding to the 1D view of a `n × n` array at the CartesianIndex `II`.

"""
@inline lexicographic(II, n) = muladd(n, II[2]-1, II[1])
@inline within_bounds(i, n) = 1 <= i <= n^2

const newaxis = [CartesianIndex()]

@inline δy⁺(II) = CartesianIndex(II[1]+1, II[2])
@inline δy⁻(II) = CartesianIndex(II[1]-1, II[2])
@inline δx⁺(II) = CartesianIndex(II[1], II[2]+1)
@inline δx⁻(II) = CartesianIndex(II[1], II[2]-1)

@inline opposite(α) = ifelse(α >= 0. ,α - π, α + π)

@inline distance(A::Point, B::Point) = sqrt((A.x-B.x)^2 + (A.y-B.y)^2)

@inline midpoint(A::Point, B::Point) = Point((A.x + B.x)/2 - 0.5, (A.y + B.y)/2 - 0.5)

@inline <(A::Point{T}, L::T) where {T} = ifelse(abs(A.x) <= L/2 && abs(A.y) <= L/2, true, false)
@inline isnan(A::Point{T}) where {T} = ifelse(isnan(A.x) || isnan(A.y), true, false)
@inline +(A::Point{T}, B::Point{T}) where {T} = Point(A.x + B.x, A.y + B.y)
@inline -(A::Point{T}, B::Point{T}) where {T} = Point(A.x - B.x, A.y - B.y)
@inline *(A::Point{T}, m::N) where {T, N} = Point(m*A.x, m*A.y)
@inline *(m::N, A::Point{T}) where {T, N} = Point(m*A.x, m*A.y)
@inline abs(A::Point) = Point(abs(A.x), abs(A.y))

@inline mysign(a,b) = a/sqrt(a^2 + b^2)
@inline mysign(a) = ifelse(a >= 0, 1, -1)

@inline ⁺(a) = max(a,0)
@inline ⁻(a) = min(a,0)

@inline c∇x(u, II) = @inbounds u[δx⁺(II)] - u[δx⁻(II)]
@inline c∇y(u, II) = @inbounds u[δy⁺(II)] - u[δy⁻(II)]
@inline c∇x(u, II, h) = @inbounds (u[δx⁺(II)] - u[δx⁻(II)])/2h
@inline c∇y(u, II, h) = @inbounds (u[δy⁺(II)] - u[δy⁻(II)])/2h

@inline ∇x⁺(u, II) = @inbounds u[δx⁺(II)] - u[II]
@inline ∇y⁺(u, II) = @inbounds u[δy⁺(II)] - u[II]
@inline ∇x⁻(u, II) = @inbounds u[δx⁻(II)] - u[II]
@inline ∇y⁻(u, II) = @inbounds u[δy⁻(II)] - u[II]

@inline normal(u, II) = @SVector [mysign(c∇x(u, II), c∇y(u, II)), mysign(c∇y(u, II), c∇x(u, II))]

@inline minmod(a, b) = ifelse(a*b <= 0, 0.0, myabs(a,b))
@inline myabs(a, b) = ifelse(abs(a) < abs(b), a, b)

@inline Dxx(u, II, h) = @inbounds (u[δx⁺(II)] - 2u[II] + u[δx⁻(II)])/h^2
@inline Dyy(u, II, h) = @inbounds (u[δy⁺(II)] - 2u[II] + u[δy⁻(II)])/h^2

@inline Dxx(u, II) = @inbounds u[δx⁺(II)] - 2u[II] + u[δx⁻(II)]
@inline Dyy(u, II) = @inbounds u[δy⁺(II)] - 2u[II] + u[δy⁻(II)]

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

@inline function static_stencil(a, II::CartesianIndex)
   return @inbounds SA_F64[a[II.I[1]-1, II.I[2]-1] a[II.I[1]-1, II.I[2]] a[II.I[1]-1, II.I[2]+1];
               a[II.I[1], II.I[2]-1] a[II.I[1], II.I[2]] a[II.I[1], II.I[2]+1];
               a[II.I[1]+1, II.I[2]-1] a[II.I[1]+1, II.I[2]] a[II.I[1]+1, II.I[2]+1]]
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


@inline parabola_fit_curvature(itp, mid_point, h) =
    @inbounds (2*itp[1,3]/((1+(2*itp[1,3]*mid_point.y + itp[2,3])^2)^1.5) + 2*itp[3,1]/((1+(2*itp[3,1]*mid_point.x + itp[3,2])^2)^1.5))/h^2


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

@inline function Richardson_extrapolation(e, r)
    p = log(abs(e[end-2] - e[end-1])/abs(e[end-1] - e[end]))/log(r)
    ext = e[end] + (e[end] - e[end-1])/(2^p - 1)
    return abs.(e .- ext)
end

function fit_order(x, y)
    coeffs = fit(log.(x), log.(y), 1).coeffs
    return exp(coeffs[1]), -coeffs[2]
end

@inline function arc_length2(POS, ind, h)
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

@inline is_dirichlet(::Dirichlet) = true
@inline is_dirichlet(::BoundaryCondition) = false

@inline is_neumann(::Neumann) = true
@inline is_neumann(::BoundaryCondition) = false

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
