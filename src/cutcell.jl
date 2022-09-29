@inline SOUTH_face(itp, g, II_0, p=Point(0.0,-0.5), dx=1.0) = dx/2 - p.x + find_zero(x -> biquadratic(itp, x, p.y/(g.y[δy⁺(II_0)] - g.y[δy⁻(II_0)])), (-dx/2+p.x,dx/2+p.x)./(g.x[δx⁺(II_0)] - g.x[δx⁻(II_0)]), FalsePosition(), maxevals = 10, atol = 1e-9) * (g.x[δx⁺(II_0)] - g.x[δx⁻(II_0)])
@inline WEST_face(itp, g, II_0, p=Point(-0.5,0.0), dy=1.0) = dy/2 - p.y + find_zero(y -> biquadratic(itp, p.x/(g.x[δx⁺(II_0)] - g.x[δx⁻(II_0)]), y), (-dy/2+p.y,dy/2+p.y)./(g.y[δy⁺(II_0)] - g.y[δy⁻(II_0)]), FalsePosition(), maxevals = 10, atol = 1e-9) * (g.y[δy⁺(II_0)] - g.y[δy⁻(II_0)])
@inline NORTH_face(itp, g, II_0, p=Point(0.0,0.5), dx=1.0) = dx/2 - p.x + find_zero(x -> biquadratic(itp, x, p.y/(g.y[δy⁺(II_0)] - g.y[δy⁻(II_0)])), (-dx/2+p.x,dx/2+p.x)./(g.x[δx⁺(II_0)] - g.x[δx⁻(II_0)]), FalsePosition(), maxevals = 10, atol = 1e-9) * (g.x[δx⁺(II_0)] - g.x[δx⁻(II_0)])
@inline EAST_face(itp, g, II_0, p=Point(0.5,0.0), dy=1.0) = dy/2 - p.y + find_zero(y -> biquadratic(itp, p.x/(g.x[δx⁺(II_0)] - g.x[δx⁻(II_0)]), y), (-dy/2+p.y,dy/2+p.y)./(g.y[δy⁺(II_0)] - g.y[δy⁻(II_0)]), FalsePosition(), maxevals = 10, atol = 1e-9) * (g.y[δy⁺(II_0)] - g.y[δy⁻(II_0)])

@inline WE(a, cases, II, l_face, s_face) = ismixed(cases[δx⁻(II)]) ? 0.5*(a[II,1] + a[δx⁻(II), 3]) : (is_liquid(cases[δx⁻(II)]) ? l_face : s_face)
@inline WE_border(a, cases, II, l_face, s_face) = a[II,1]
@inline EW(a, cases, II, l_face, s_face) = ismixed(cases[δx⁺(II)]) ? 0.5*(a[II,3] + a[δx⁺(II), 1]) : (is_liquid(cases[δx⁺(II)]) ? l_face : s_face)
@inline EW_border(a, cases, II, l_face, s_face) = a[II,3]
@inline SN(a, cases, II, l_face, s_face) = ismixed(cases[δy⁻(II)]) ? 0.5*(a[II,2] + a[δy⁻(II), 4]) : (is_liquid(cases[δy⁻(II)]) ? l_face : s_face)
@inline SN_border(a, cases, II, l_face, s_face) = a[II,2]
@inline NS(a, cases, II, l_face, s_face) = ismixed(cases[δy⁺(II)]) ? 0.5*(a[II,4] + a[δy⁺(II), 2]) : (is_liquid(cases[δy⁺(II)]) ? l_face : s_face)
@inline NS_border(a, cases, II, l_face, s_face) = a[II,4]

@inline WE(a, II) = 0.5*(a[II,1] + a[δx⁻(II), 3])
@inline EW(a, II) = 0.5*(a[II,3] + a[δx⁺(II), 1])
@inline SN(a, II) = 0.5*(a[II,2] + a[δy⁻(II), 4])
@inline NS(a, II) = 0.5*(a[II,4] + a[δy⁺(II), 2])

@inline ispositive(a) = min(max(a,0),1)
@inline isovalue(v::SVector{4,Float64}) = 4. * v[1] + 2. * v[2] + v[3] + 0.5 * v[4]

@inline vertices_sign(itp, II_0, II, x, y, dx, dy) = 
    1 .- sign.(@SVector [biquadratic(itp, (x[II] - x[II_0] - dx/2)/(x[δx⁺(II_0)] - x[δx⁻(II_0)]),
                            (y[II] - y[II_0] + dy/2)/(y[δy⁺(II_0)] - y[δy⁻(II_0)])),
                        biquadratic(itp, (x[II] - x[II_0] + dx/2)/(x[δx⁺(II_0)] - x[δx⁻(II_0)]),
                            (y[II] - y[II_0] + dy/2)/(y[δy⁺(II_0)] - y[δy⁻(II_0)])),
                        biquadratic(itp, (x[II] - x[II_0] + dx/2)/(x[δx⁺(II_0)] - x[δx⁻(II_0)]),
                            (y[II] - y[II_0] - dy/2)/(y[δy⁺(II_0)] - y[δy⁻(II_0)])),
                        biquadratic(itp, (x[II] - x[II_0] - dx/2)/(x[δx⁺(II_0)] - x[δx⁻(II_0)]),
                            (y[II] - y[II_0] - dy/2)/(y[δy⁺(II_0)] - y[δy⁻(II_0)]))])

@inline face_pos(II_0, II, x, y, dx, dy) =
    (Point(x[II] - x[II_0] - dx/2, y[II] - y[II_0]),
    Point(x[II] - x[II_0], y[II] - y[II_0] - dy/2),
    Point(x[II] - x[II_0] + dx/2, y[II] - y[II_0]),
    Point(x[II] - x[II_0], y[II] - y[II_0] + dy/2))

@inline is_near_interface(a, st) = @inbounds ifelse(a != sign(st[2,1]) || a != sign(st[1,2]) || a != sign(st[2,3]) || a != sign(st[3,2]), true, false)

@inline is_near_interface(u, II::CartesianIndex) = @inbounds ifelse(u[II]*u[δx⁺(II)] < 0 || u[II]*u[δy⁺(II)] < 0 || u[II]*u[δx⁻(II)] < 0 || u[II]*u[δy⁻(II)] < 0, true, false)
@inline is_near_interface(u, II::CartesianIndex, nx, ny, per_x, per_y) = @inbounds ifelse(u[II]*u[δx⁺(II, nx, per_x)] < 0 || u[II]*u[δy⁺(II, ny, per_y)] < 0 || u[II]*u[δx⁻(II, nx, per_x)] < 0 || u[II]*u[δy⁻(II, ny, per_y)] < 0, true, false)

@inline ismixed(a) = a !== 0.0 && a !== 15.0
@inline is_not_mixed(a) = a == 0.0 || a == 15.0
@inline is_liquid(a) = a == 0.0
@inline is_solid(a) = a == 15.0

@inline biquadratic(m, x, y) = @inbounds m[1,1]*(y^2)*(x^2) + m[1,2]*(y^2)*x + m[2,1]*y*x^2 + m[1,3]*y^2 + m[3,1]*x^2 + m[2,2]*y*x + m[2,3]*y + m[3,2]*x + m[3,3]

@inline check_grid_alignement(u, II, κ, ϵ) = @SVector [abs(u[II] + u[δx⁺(II)]) <= ϵ && sign(u[δy⁺(II)]) == sign(u[δy⁻(II)]) && abs(κ) < ϵ,
    abs(u[II] + u[δx⁻(II)]) <= ϵ && sign(u[δy⁺(II)]) == sign(u[δy⁻(II)]) && abs(κ) < ϵ,
    abs(u[II] + u[δy⁺(II)]) <= ϵ && sign(u[δx⁺(II)]) == sign(u[δx⁻(II)]) && abs(κ) < ϵ,
    abs(u[II] + u[δy⁻(II)]) <= ϵ && sign(u[δx⁺(II)]) == sign(u[δx⁻(II)]) && abs(κ) < ϵ]

@inline function is_aligned_with_grid(a)
    @inbounds for i in a
        if i == true
            return true
        end
    end
    return false
end

@inline function get_centroid(nodes)
   x = 0.; y = 0.; area = 0.; k = 0.;
   c = nodes[end]
   d = nodes[end]
   @inbounds for i in 1:length(nodes)
      c = nodes[i]

      k = c.y*d.x - c.x*d.y
      area += k
      x += (c.x + d.x)*k
      y += (c.y + d.y)*k

      d = c
   end
   area *= 3
   return Point((x/=area) - 0.5, (y/=area) - 0.5)
end

function find_radius(grid, MIXED, h)
    @unpack x, y, mid_point = grid
    H = x[1,:]
    radius = Vector{Float64}(undef,0)
    for II in MIXED
        try
            radius_ = distance(Point(0.,0.), Point(H[II.I[2]]+h*mid_point[II].x, H[II.I[1]]+h*mid_point[II].y))
            if !isnan(radius_) push!(radius, radius_) end
        catch
            @goto endofloop
        end
        @label endofloop
    end
    return mean(radius)
end

function find_position_hill(grid, MIXED, h)
    @unpack x, y, mid_point = grid
    H = x[1,:]
    radius = Vector{Float64}(undef,0)
    for II in MIXED
        try
            radius_ = H[II.I[1]]+h*mid_point[II].y
            if !isnan(radius_) push!(radius, radius_) end
        catch
            @goto endofloop
        end
        @label endofloop
    end
    return mean(radius)
end


@inline function line_intersection(l1::Line, l2::Line)
    x1, y1 = l1.p1.x, l1.p1.y
    x2, y2 = l1.p2.x, l1.p2.y
    x3, y3 = l2.p1.x, l2.p1.y
    x4, y4 = l2.p2.x, l2.p2.y

    d = (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4)
    px = ((x1*y2-y1*x2)*(x3-x4) - (x1-x2)*(x3*y4-y3*x4)) / d
    py = ((x1*y2-y1*x2)*(y3-y4) - (y1-y2)*(x3*y4-y3*x4)) / d

    return Point(px,py)
end

@inline function Bcaps(ip, centroid)
    if isnan(ip) || isinf(ip) || abs(ip) >= 0.5
        cap = 1.0
    elseif sign(ip) == sign(centroid)
        cap = 0.5 - abs(ip)
    else
        cap = 0.5 + abs(ip)
    end
    return cap
end

function clip_large_cell!(grid, geo, II, ϵ, nx, ny)
    if geo.cap[II,5] > (1.0-ϵ)
        geo.cap[II,:] .= vcat(ones(7), 0.5.*ones(4))
        if II[2] > 1
            geo.cap[δx⁻(II),3] = 1.0
        end
        if II[1] > 1
            geo.cap[δy⁻(II),4] = 1.0
        end
        if II[2] < nx
            geo.cap[δx⁺(II),1] = 1.0
        end
        if II[1] < ny
            geo.cap[δy⁺(II),2] = 1.0
        end
        geo.centroid[II] = Point(0.0, 0.0)
        grid.mid_point[II] = Point(0.0, 0.0)
    end
    return nothing
end

function clip_small_cell!(grid, geo, II, ϵ, nx, ny)
    if geo.cap[II,5] < ϵ
        geo.cap[II,:] .= 0.
        geo.emptied[II] = true
        if II[2] > 1
            geo.cap[δx⁻(II),3] = 0.
        end
        if II[1] > 1
            geo.cap[δy⁻(II),4] = 0.
        end
        if II[2] < nx
            geo.cap[δx⁺(II),1] = 0.
        end
        if II[1] < ny
            geo.cap[δy⁺(II),2] = 0.
        end
        geo.centroid[II] = Point(0.0, 0.0)
        grid.mid_point[II] = Point(0.0, 0.0)
    end
    return nothing
end

function clip_cells!(grid, ϵ, periodic_x, periodic_y)
    @unpack nx, ny, ind, geoS, geoL = grid

    @inbounds @threads for II in ind.inside
        clip_large_cell!(grid, geoS, II, ϵ, nx, ny)
        clip_large_cell!(grid, geoL, II, ϵ, nx, ny)
    end
    @inbounds @threads for II in ind.inside
        clip_small_cell!(grid, geoS, II, ϵ, nx, ny)
        clip_small_cell!(grid, geoL, II, ϵ, nx, ny)
    end
    @inbounds @threads for II in @view ind.b_left[1][2:end-1]
        if geoS.cap[II,5] > (1.0-ϵ)
            if periodic_x
                geoS.cap[II,1] = 1.0
            end
            geoS.cap[II,2:11] .= vcat(ones(6), 0.5.*ones(4))
            geoS.cap[δy⁻(II),4] = 1.0
            geoS.cap[δx⁺(II),1] = 1.0
            geoS.cap[δy⁺(II),2] = 1.0
        end
    end
    @inbounds @threads for II in @view ind.b_left[1][2:end-1]
        if geoS.cap[II,5] < ϵ
            geoS.cap[II,:] .= 0.
            geoS.emptied[II] = true
            geoS.cap[δy⁻(II),4] = 0.
            geoS.cap[δx⁺(II),1] = 0.
            geoS.cap[δy⁺(II),2] = 0.
        end
    end
    @inbounds @threads for II in @view ind.b_bottom[1][2:end-1]
        if geoS.cap[II,5] > (1.0-ϵ)
            geoS.cap[II,1] = 1.0
            if periodic_y
                geoS.cap[II,2] = 1.0
            end
            geoS.cap[II,3:11] .= vcat(ones(5), 0.5.*ones(4))
            geoS.cap[δx⁻(II),3] = 1.0
            geoS.cap[δx⁺(II),1] = 1.0
            geoS.cap[δy⁺(II),2] = 1.0
        end
    end
    @inbounds @threads for II in @view ind.b_bottom[1][2:end-1]
        if geoS.cap[II,5] < ϵ
            geoS.cap[II,:] .= 0.
            geoS.emptied[II] = true
            geoS.cap[δx⁻(II),3] = 0.
            geoS.cap[δx⁺(II),1] = 0.
            geoS.cap[δy⁺(II),2] = 0.
        end
    end
    @inbounds @threads for II in @view ind.b_right[1][2:end-1]
        if geoS.cap[II,5] > (1.0-ϵ)
            geoS.cap[II,1:2] .= ones(2)
            if periodic_x
                geoS.cap[II,3] = 1.0
            end
            geoS.cap[II,4:11] .= vcat(ones(4), 0.5.*ones(4))
            geoS.cap[δx⁻(II),3] = 1.0
            geoS.cap[δy⁻(II),4] = 1.0
            geoS.cap[δy⁺(II),2] = 1.0
        end
    end
    @inbounds @threads for II in @view ind.b_right[1][2:end-1]
        if geoS.cap[II,5] < ϵ
            geoS.cap[II,:] .= 0.
            geoS.emptied[II] = true
            geoS.cap[δx⁻(II),3] = 0.
            geoS.cap[δy⁻(II),4] = 0.
            geoS.cap[δy⁺(II),2] = 0.
        end
    end
    @inbounds @threads for II in @view ind.b_top[1][2:end-1]
        if geoS.cap[II,5] > (1.0-ϵ)
            geoS.cap[II,1:3] .= ones(3)
            if periodic_y
                geoS.cap[II,4] = 1.0
            end
            geoS.cap[II,5:11] .= vcat(ones(3), 0.5.*ones(4))
            geoS.cap[δx⁻(II),3] = 1.0
            geoS.cap[δy⁻(II),4] = 1.0
            geoS.cap[δx⁺(II),1] = 1.0
        end
    end
    @inbounds @threads for II in @view ind.b_top[1][2:end-1]
        if geoS.cap[II,5] < ϵ
            geoS.cap[II,:] .= 0.
            geoS.emptied[II] = true
            geoS.cap[δx⁻(II),3] = 0.
            geoS.cap[δy⁻(II),4] = 0.
            geoS.cap[δx⁺(II),1] = 0.
        end
    end
    @inbounds @threads for II in @view ind.b_left[1][2:end-1]
        if geoL.cap[II,5] > (1.0-ϵ)
            if periodic_x
                geoL.cap[II,1] = 1.0
            end
            geoL.cap[II,2:11] .= vcat(ones(6), 0.5.*ones(4))
            geoL.cap[δy⁻(II),4] = 1.0
            geoL.cap[δx⁺(II),1] = 1.0
            geoL.cap[δy⁺(II),2] = 1.0
        end
    end
    @inbounds @threads for II in @view ind.b_left[1][2:end-1]
        if geoL.cap[II,5] < ϵ
            geoL.cap[II,:] .= 0.
            geoL.emptied[II] = true
            geoL.cap[δy⁻(II),4] = 0.
            geoL.cap[δx⁺(II),1] = 0.
            geoL.cap[δy⁺(II),2] = 0.
        end
    end
    @inbounds @threads for II in @view ind.b_bottom[1][2:end-1]
        if geoL.cap[II,5] > (1.0-ϵ)
            geoL.cap[II,1] = 1.0
            if periodic_y
                geoL.cap[II,2] = 1.0
            end
            geoL.cap[II,3:11] .= vcat(ones(5), 0.5.*ones(4))
            geoL.cap[δx⁻(II),3] = 1.0
            geoL.cap[δx⁺(II),1] = 1.0
            geoL.cap[δy⁺(II),2] = 1.0
        end
    end
    @inbounds @threads for II in @view ind.b_bottom[1][2:end-1]
        if geoL.cap[II,5] < ϵ
            geoL.cap[II,:] .= 0.
            geoL.emptied[II] = true
            geoL.cap[δx⁻(II),3] = 0.
            geoL.cap[δx⁺(II),1] = 0.
            geoL.cap[δy⁺(II),2] = 0.
        end
    end
    @inbounds @threads for II in @view ind.b_right[1][2:end-1]
        if geoL.cap[II,5] > (1.0-ϵ)
            geoL.cap[II,1:2] .= ones(2)
            if periodic_x
                geoL.cap[II,3] = 1.0
            end
            geoL.cap[II,4:11] .= vcat(ones(4), 0.5.*ones(4))
            geoL.cap[δx⁻(II),3] = 1.0
            geoL.cap[δy⁻(II),4] = 1.0
            geoL.cap[δy⁺(II),2] = 1.0
        end
    end
    @inbounds @threads for II in @view ind.b_right[1][2:end-1]
        if geoL.cap[II,5] < ϵ
            geoL.cap[II,:] .= 0.
            geoL.emptied[II] = true
            geoL.cap[δx⁻(II),3] = 0.
            geoL.cap[δy⁻(II),4] = 0.
            geoL.cap[δy⁺(II),2] = 0.
        end
    end
    @inbounds @threads for II in @view ind.b_top[1][2:end-1]
        if geoL.cap[II,5] > (1.0-ϵ)
            geoL.cap[II,1:3] .= ones(3)
            if periodic_y
                geoL.cap[II,4] = 1.0
            end
            geoL.cap[II,5:11] .= vcat(ones(3), 0.5.*ones(4))
            geoL.cap[δx⁻(II),3] = 1.0
            geoL.cap[δy⁻(II),4] = 1.0
            geoL.cap[δx⁺(II),1] = 1.0
        end
    end
    @inbounds @threads for II in @view ind.b_top[1][2:end-1]
        if geoL.cap[II,5] < ϵ
            geoL.cap[II,:] .= 0.
            geoL.emptied[II] = true
            geoL.cap[δx⁻(II),3] = 0.
            geoL.cap[δy⁻(II),4] = 0.
            geoL.cap[δx⁺(II),1] = 0.
        end
    end
    return nothing
end

function empty_cell_u!(geo, geo_u, II, ny)
    if geo.cap[II,3] < 1e-10
        geo_u.dcap[δx⁺(II),:] .= 0.
        geo_u.emptied[δx⁺(II)] = true
        geo_u.dcap[II,3] = 0.
        geo_u.dcap[δx⁺(δx⁺(II)),1] = 0.
        if II[1] != 1
            geo_u.dcap[δy⁻(δx⁺(II)),4] = 0.
        end
        if II[1] != ny
            geo_u.dcap[δy⁺(δx⁺(II)),2] = 0.
        end
    end
    return nothing
end

function empty_cell_v!(geo, geo_v, II, nx)
    if geo.cap[II,4] < 1e-10
        geo_v.dcap[δy⁺(II),:] .= 0.
        geo_v.emptied[δy⁺(II)] = true
        geo_v.dcap[II,4] = 0.
        geo_v.dcap[δy⁺(δy⁺(II)),2] = 0.
        if II[2] != 1
            geo_v.dcap[δx⁻(δy⁺(II)),3] = 0.
        end
        if II[2] != nx
            geo_v.dcap[δx⁺(δy⁺(II)),1] = 0.
        end
    end
    return nothing
end

function empty_cells_uv!(grid, grid_u, grid_v)
    @unpack nx, ny, ind, geoS, geoL = grid

    # If A_i = 0, V_i = 0
    @inbounds @threads for II in ind.inside
        empty_cell_u!(geoS, grid_u.geoS, II, ny)
        empty_cell_u!(geoL, grid_u.geoL, II, ny)

        empty_cell_v!(geoS, grid_v.geoS, II, nx)
        empty_cell_v!(geoL, grid_v.geoL, II, nx)
    end
    @inbounds @threads for II in ind.b_left[1]
        if geoS.cap[II,1] < 1e-10
            grid_u.geoS.dcap[II,:] .= 0.
            grid_u.geoS.emptied[II] = true
            grid_u.geoS.dcap[δx⁺(II),1] = 0.
            if II[1] != 1
                grid_u.geoS.dcap[δy⁻(II),4] = 0.
            end
            if II[1] != ny
                grid_u.geoS.dcap[δy⁺(II),2] = 0.
            end
        end
        if geoL.cap[II,1] < 1e-10
            grid_u.geoL.dcap[II,:] .= 0.
            grid_u.geoL.emptied[II] = true
            grid_u.geoL.dcap[δx⁺(II),1] = 0.
            if II[1] != 1
                grid_u.geoL.dcap[δy⁻(II),4] = 0.
            end
            if II[1] != ny
                grid_u.geoL.dcap[δy⁺(II),2] = 0.
            end
        end
    end
    @inbounds @threads for II in ind.b_left[1]
        empty_cell_u!(geoS, grid_u.geoS, II, ny)
        empty_cell_u!(geoL, grid_u.geoL, II, ny)
    end
    @inbounds @threads for II in ind.b_right[1]
        if geoS.cap[II,3] < 1e-10
            grid_u.geoS.dcap[δx⁺(II),:] .= 0.
            grid_u.geoS.emptied[II] = true
            grid_u.geoS.dcap[II,3] = 0.
            if II[1] != 1
                grid_u.geoS.dcap[δy⁻(δx⁺(II)),4] = 0.
            end
            if II[1] != ny
                grid_u.geoS.dcap[δy⁺(δx⁺(II)),2] = 0.
            end
        end
        if geoL.cap[II,3] < 1e-10
            grid_u.geoL.dcap[δx⁺(II),:] .= 0.
            grid_u.geoL.emptied[II] = true
            grid_u.geoL.dcap[II,3] = 0.
            if II[1] != 1
                grid_u.geoL.dcap[δy⁻(δx⁺(II)),4] = 0.
            end
            if II[1] != ny
                grid_u.geoL.dcap[δy⁺(δx⁺(II)),2] = 0.
            end
        end
    end
    @inbounds @threads for II in vcat(ind.b_bottom[1][2:end-1], ind.b_top[1][2:end-1])
        if geoS.cap[II,3] < 1e-10
            grid_u.geoS.dcap[δx⁺(II),:] .= 0.
            grid_u.geoS.emptied[δx⁺(II)] = true
            grid_u.geoS.dcap[II,3] = 0.
            grid_u.geoS.dcap[δx⁺(δx⁺(II)),1] = 0.
        end
        if geoL.cap[II,3] < 1e-10
            grid_u.geoL.dcap[δx⁺(II),:] .= 0.
            grid_u.geoL.emptied[δx⁺(II)] = true
            grid_u.geoL.dcap[II,3] = 0.
            grid_u.geoL.dcap[δx⁺(δx⁺(II)),1] = 0.
        end
    end
    @inbounds @threads for II in ind.b_bottom[1]
        if geoS.cap[II,2] < 1e-10
            grid_v.geoS.dcap[II,:] .= 0.
            grid_v.geoS.emptied[II] = true
            grid_v.geoS.dcap[δy⁺(II),2] = 0.
            if II[2] != 1
                grid_v.geoS.dcap[δx⁻(II),3] = 0.
            end
            if II[2] != nx
                grid_v.geoS.dcap[δx⁺(II),1] = 0.
            end
        end
        if geoL.cap[II,2] < 1e-10
            grid_v.geoL.dcap[II,:] .= 0.
            grid_v.geoL.emptied[II] = true
            grid_v.geoL.dcap[δy⁺(II),2] = 0.
            if II[2] != 1
                grid_v.geoL.dcap[δx⁻(II),3] = 0.
            end
            if II[2] != nx
                grid_v.geoL.dcap[δx⁺(II),1] = 0.
            end
        end
    end
    @inbounds @threads for II in ind.b_bottom[1]
        empty_cell_v!(geoS, grid_v.geoS, II, nx)
        empty_cell_v!(geoL, grid_v.geoL, II, nx)
    end
    @inbounds @threads for II in ind.b_top[1]
        if geoS.cap[II,4] < 1e-10
            grid_v.geoS.dcap[δy⁺(II),:] .= 0.
            grid_v.geoS.emptied[II] = true
            grid_v.geoS.dcap[II,4] = 0.
            if II[2] != 1
                grid_v.geoS.dcap[δx⁻(δy⁺(II)),3] = 0.
            end
            if II[2] != nx
                grid_v.geoS.dcap[δx⁺(δy⁺(II)),1] = 0.
            end
        end
        if geoL.cap[II,4] < 1e-10
            grid_v.geoL.dcap[δy⁺(II),:] .= 0.
            grid_v.geoL.emptied[II] = true
            grid_v.geoL.dcap[II,4] = 0.
            if II[2] != 1
                grid_v.geoL.dcap[δx⁻(δy⁺(II)),3] = 0.
            end
            if II[2] != nx
                grid_v.geoL.dcap[δx⁺(δy⁺(II)),1] = 0.
            end
        end
    end
    @inbounds @threads for II in vcat(ind.b_left[1][2:end-1], ind.b_right[1][2:end-1])
        if geoS.cap[II,4] < 1e-10
            grid_v.geoS.dcap[δy⁺(II),:] .= 0.
            grid_v.geoS.emptied[δy⁺(II)] = true
            grid_v.geoS.dcap[II,4] = 0.
            grid_v.geoS.dcap[δy⁺(δy⁺(II)),2] = 0.
        end
        if geoL.cap[II,4] < 1e-10
            grid_v.geoL.dcap[δy⁺(II),:] .= 0.
            grid_v.geoL.emptied[δy⁺(II)] = true
            grid_v.geoL.dcap[II,4] = 0.
            grid_v.geoL.dcap[δy⁺(δy⁺(II)),2] = 0.
        end
    end
end

function dimensionalize!(grid, geo)
    @unpack dx, dy = grid
    @unpack cap, dcap = geo

    @inbounds dcap[:,:,1] .= view(cap,:,:,1) .* dy
    @inbounds dcap[:,:,2] .= view(cap,:,:,2) .* dx
    @inbounds dcap[:,:,3] .= view(cap,:,:,3) .* dy
    @inbounds dcap[:,:,4] .= view(cap,:,:,4) .* dx
    @inbounds dcap[:,:,5] .= view(cap,:,:,5) .* dx .* dy
    @inbounds dcap[:,:,6] .= view(cap,:,:,6) .* dy
    @inbounds dcap[:,:,7] .= view(cap,:,:,7) .* dx
    @inbounds dcap[:,:,8] .= view(cap,:,:,8) .* dx .* dy
    @inbounds dcap[:,:,9] .= view(cap,:,:,9) .* dx .* dy
    @inbounds dcap[:,:,10] .= view(cap,:,:,10) .* dx .* dy
    @inbounds dcap[:,:,11] .= view(cap,:,:,11) .* dx .* dy

    return nothing
end

function postprocess_grids!(grid, grid_u, grid_v, MIXED, periodic_x, periodic_y, advection, ϵ)
    clip_cells!(grid, ϵ, periodic_x, periodic_y)
    
    dimensionalize!(grid, grid.geoS)
    dimensionalize!(grid, grid.geoL)
    dimensionalize!(grid_u, grid_u.geoS)
    dimensionalize!(grid_u, grid_u.geoL)
    dimensionalize!(grid_v, grid_v.geoS)
    dimensionalize!(grid_v, grid_v.geoL)

    empty_cells_uv!(grid, grid_u, grid_v)

    set_cap_bcs!(grid, periodic_x, periodic_y)
    set_cap_bcs!(grid_u, periodic_x, periodic_y)
    set_cap_bcs!(grid_v, periodic_x, periodic_y)

    Wcapacities!(grid.geoS.dcap, periodic_x, periodic_y)
    Wcapacities!(grid.geoL.dcap, periodic_x, periodic_y)
    Wcapacities!(grid_u.geoS.dcap, periodic_x, periodic_y)
    Wcapacities!(grid_u.geoL.dcap, periodic_x, periodic_y)
    Wcapacities!(grid_v.geoS.dcap, periodic_x, periodic_y)
    Wcapacities!(grid_v.geoL.dcap, periodic_x, periodic_y)

    average_face_capacities!(grid.geoS.dcap)
    average_face_capacities!(grid.geoL.dcap)
    average_face_capacities!(grid_u.geoS.dcap)
    average_face_capacities!(grid_u.geoL.dcap)
    average_face_capacities!(grid_v.geoS.dcap)
    average_face_capacities!(grid_v.geoL.dcap)
end

function marching_squares!(num, grid)
    @unpack x, y, nx, ny, dx, dy, ind, u, iso, faces, geoS, geoL, mid_point, α, κ = grid

    empty_capacities = vcat(zeros(7), zeros(4))
    full_capacities = vcat(ones(7), 0.5.*ones(4))
    κ .= zeros(ny,nx)
    @inbounds @threads for II in ind.all_indices
        if II in ind.inside
            II_0 = II
        elseif II in ind.b_left[1][2:end-1]
            II_0 = δx⁺(II)
        elseif II in ind.b_bottom[1][2:end-1]
            II_0 = δy⁺(II)
        elseif II in ind.b_right[1][2:end-1]
            II_0 = δx⁻(II)
        elseif II in ind.b_top[1][2:end-1]
            II_0 = δy⁻(II)
        elseif II == ind.b_left[1][1]
            II_0 = δy⁺(δx⁺(II))
        elseif II == ind.b_left[1][end]
            II_0 = δy⁻(δx⁺(II))
        elseif II == ind.b_right[1][1]
            II_0 = δy⁺(δx⁻(II))
        elseif II == ind.b_right[1][end]
            II_0 = δy⁻(δx⁻(II))
        end

        st = static_stencil(u, II_0)
        posW, posS, posE, posN = face_pos(II_0, II, x, y, dx[II], dy[II])
        
        a = sign(u[II])
        ISO = ifelse(a > 0, 0., 15.)
        if is_near_interface(a, st)
            κ_ = 0.
            if II in ind.inside
                grid_alignement = check_grid_alignement(u, II, κ_, 0.0)
            else
                grid_alignement = false
            end
            if is_aligned_with_grid(grid_alignement)
                if a > 0
                    ISO = 0.
                    @goto notmixed
                else
                    ISO = -1.
                    κ[II] = 0.
                    geoS.cap[II,:] .= full_capacities
                    geoL.cap[II,:] .= empty_capacities
                end
            else
                B, BT = B_BT(II_0, x, y)
                itp = B * st * BT
                vertices = vertices_sign(itp, II_0, II, x, y, dx[II], dy[II])
                ISO = isovalue(vertices)

                if is_not_mixed(ISO) @goto notmixed end
                face_capacities(grid, itp, ISO, II_0, II, posW, posS, posE, posN)

            end
            if ISO == -1
                geoS.projection[II].flag = false
                geoL.projection[II].flag = false
            end
        end
        @label notmixed
        if is_solid(ISO)
            geoS.cap[II,:] .= full_capacities
            geoL.cap[II,:] .= empty_capacities
        elseif is_liquid(ISO)
            geoS.cap[II,:] .= empty_capacities
            geoL.cap[II,:] .= full_capacities
        end
        iso[II] = ISO
    end
    
    return nothing
end

function get_iterface_location!(grid, indices)
    @unpack x, y, iso, geoS, geoL, mid_point, cut_points, α = grid

    @inbounds @threads for II in indices
        f = average_face_capacities(grid, II)
        geoS.cap[II,:], geoL.cap[II,:], α[II], geoS.centroid[II], geoL.centroid[II], mid_point[II], cut_points[II] = capacities(f, iso[II])
        geoS.projection[II], geoL.projection[II] = projection_2points(grid, II)
    end
    return nothing
end

function get_interface_location_borders!(grid::Mesh{GridFCx,T,N}, periodic_x, periodic_y) where {T,N}
    @unpack nx, ny, ind, u, geoS, geoL, mid_point, cut_points = grid
    @unpack b_left, b_bottom, b_right, b_top = ind

    f = SA_F64[0.5, 0.5]
    empty_capacities = SA_F64[vcat(zeros(7), zeros(4))...]
    capacities_6 = capacities(f, 6.0)
    capacities_9 = capacities(f, 9.0)

    if !periodic_x
        @inbounds @threads for II in b_left[1]
            if u[II] >= 0.0
                @inbounds geoS.cap[II,:] .= empty_capacities
                @inbounds _, geoL.cap[II,:], _, _, geoL.centroid[II], mid_point[II], cut_points[II] = capacities_9
            else
                @inbounds geoS.cap[II,:], _, _, geoS.centroid[II], _, mid_point[II], cut_points[II] = capacities_6
                @inbounds geoL.cap[II,:] .= empty_capacities
            end
        end
        @inbounds @threads for II in b_right[1]
            if u[II] >= 0.0
                @inbounds geoS.cap[II,:] .= empty_capacities
                @inbounds _, geoL.cap[II,:], _, _, geoL.centroid[II], mid_point[II], cut_points[II] = capacities_6
            else
                @inbounds geoS.cap[II,:], _, _, geoS.centroid[II], _, mid_point[II], cut_points[II] = capacities_9
                @inbounds geoL.cap[II,:] .= empty_capacities
            end
        end
    end
    return nothing
end

function get_interface_location_borders!(grid::Mesh{GridFCy,T,N}, periodic_x, periodic_y) where {T,N}
    @unpack ind, u, geoS, geoL, mid_point, cut_points = grid
    @unpack b_left, b_bottom, b_right, b_top = ind

    f = SA_F64[0.5, 0.5]
    empty_capacities = SA_F64[vcat(zeros(7), zeros(4))...]
    capacities_3 = capacities(f, 3.0)
    capacities_12 = capacities(f, 12.0)

    if !periodic_y
        @inbounds @threads for II in b_bottom[1]
            if u[II] >= 0.0
                @inbounds geoS.cap[II,:] .= empty_capacities
                @inbounds _, geoL.cap[II,:], _, _, geoL.centroid[II], mid_point[II], cut_points[II] = capacities_3
            else
                @inbounds geoS.cap[II,:], _, _, geoS.centroid[II], _, mid_point[II], cut_points[II] = capacities_12
                @inbounds geoL.cap[II,:] .= empty_capacities
            end
        end
        @inbounds @threads for II in b_top[1]
            if u[II] >= 0.0
                @inbounds geoS.cap[II,:] .= empty_capacities
                @inbounds _, geoL.cap[II,:], _, _, geoL.centroid[II], mid_point[II], cut_points[II] = capacities_12
            else
                @inbounds geoS.cap[II,:], _, _, geoS.centroid[II], _, mid_point[II], cut_points[II] = capacities_3
                @inbounds geoL.cap[II,:] .= empty_capacities
            end
        end
    end
end

function get_curvature(num, grid, inside)
    @unpack Δ = num
    @unpack x, y, dx, dy, ind, u, geoL, κ = grid

    @inbounds for II in inside
        if II in ind.inside
            II_0 = II
        elseif II in ind.b_left[1][2:end-1]
            II_0 = δx⁺(II)
        elseif II in ind.b_bottom[1][2:end-1]
            II_0 = δy⁺(II)
        elseif II in ind.b_right[1][2:end-1]
            II_0 = δx⁻(II)
        elseif II in ind.b_top[1][2:end-1]
            II_0 = δy⁻(II)
        elseif II == ind.b_left[1][1]
            II_0 = δy⁺(δx⁺(II))
        elseif II == ind.b_left[1][end]
            II_0 = δy⁻(δx⁺(II))
        elseif II == ind.b_right[1][1]
            II_0 = δy⁺(δx⁻(II))
        elseif II == ind.b_right[1][end]
            II_0 = δy⁻(δx⁻(II))
        end
        mid_point = geoL.projection[II].mid_point
        Δx2 = x[δx⁺(II_0)] - x[δx⁻(II_0)]
        Δy2 = y[δy⁺(II_0)] - y[δy⁻(II_0)]
        st = static_stencil(u, II_0)
        B, BT = B_BT(II_0, x, y)
        itp = B * st * BT
        # If someone finds why it works with this scaling factor
        # ( * 0.25 / (dx[II]*dy[II]) ), I pay him a beer (Alex)
        κ[II] = parabola_fit_curvature(itp, mid_point, Δx2, Δy2) * 0.25 / (dx[II]*dy[II])
    end
end

function capacities(F_prev, case)
    # add some little volume so that weird cases don't appear
    F = [F_prev[1], F_prev[2]]
    if F[1] <= 1e-8
        F[1] += 1e-8
    elseif  F[1] >= (1.0-1e-8)
        F[1] -= 1e-8
    end
    if F[2] <= 1e-8
        F[2] += 1e-8
    elseif  F[2] >= (1.0-1e-8)
        F[2] -= 1e-8
    end

    cap_sol = @SVector zeros(11)
    cap_liq = @SVector zeros(11)
    α = 0.0
    sol_centroid = Point(NaN, NaN)
    liq_centroid = Point(NaN, NaN)
    mid_point = Point(NaN, NaN)
    cut_points = [Point(NaN, NaN), Point(NaN, NaN)]
    if case == 1.0 #W S
        s1 = (Point(0.,0.), Point(F[2],0.), Point(0.,F[1]), Point(0.,0.))
        l1 = (Point(1.,0.), Point(F[2],0.), Point(0., F[1]), Point(0.,1.), Point(1.,1.), Point(1., 0.))
        sol_centroid = get_centroid(s1)
        liq_centroid = get_centroid(l1)
        α = atan(F[2], F[1])
        mid_point = midpoint(Point(F[2],0.), Point(0.,F[1]))
        cut_points = [Point(F[2],0.) - Point(0.5,0.5), Point(0.,F[1]) - Point(0.5,0.5)]
        bsol, bliq, s_ipx, s_ipy, l_ipx, l_ipy = Bcapacities(cut_points, sol_centroid, liq_centroid)
        svol = 0.5*F[1]*F[2]
        lvol = 1. - 0.5*F[1]*F[2]

        s_w3 = 0.5*(F[2]-(sol_centroid.x+0.5))*(s_ipy.y+0.5)
        s_w1 = svol - s_w3
        s_w4 = 0.5*(F[1]-(sol_centroid.y+0.5))*(s_ipx.x+0.5)
        s_w2 = svol - s_w4

        if isnan(l_ipy.y) || isinf(l_ipy.y) || abs(l_ipy.y) >= 0.5
            l_w3 = 1.0 - (liq_centroid.x+0.5)
            l_w1 = lvol - l_w3
        else
            l_w1 = 0.5*(liq_centroid.x+0.5)*(1.0-(l_ipy.y+0.5)+1.0-F[1])
            l_w3 = lvol - l_w1
        end
        if isnan(l_ipx.x) || isinf(l_ipx.x) || abs(l_ipx.x) >= 0.5
            l_w4 = 1.0 - (liq_centroid.y+0.5)
            l_w2 = lvol - l_w4
        else
            l_w2 = 0.5*(liq_centroid.y+0.5)*(1.0-(l_ipx.x+0.5)+1.0-F[2])
            l_w4 = lvol - l_w2
        end

        cap_sol = SA_F64[F[1], F[2], 0, 0, svol, bsol..., s_w1, s_w2, s_w3, s_w4]
        cap_liq = SA_F64[1.0 - F[1], 1.0 - F[2], 1.0, 1.0, lvol, bliq..., l_w1, l_w2, l_w3, l_w4]
    elseif case == 2.0 #S E
        s2 = (Point(F[1],0.), Point(1.,0.), Point(1.,F[2]), Point(F[1],0.))
        l2 = (Point(0.,0.), Point(F[1],0.), Point(1., F[2]), Point(1.,1.), Point(0.,1.), Point(0., 0.))
        sol_centroid = get_centroid(s2)
        liq_centroid = get_centroid(l2)
        α = atan(1- F[1], -F[2])
        mid_point = midpoint(Point(F[1],0.), Point(1.,F[2]))
        cut_points = [Point(F[1],0.) - Point(0.5,0.5), Point(1.,F[2]) - Point(0.5,0.5)]
        bsol, bliq, s_ipx, s_ipy, l_ipx, l_ipy = Bcapacities(cut_points, sol_centroid, liq_centroid)
        svol = 0.5*(1 - F[1])*F[2]
        lvol = 1. - 0.5*(1 - F[1])*F[2]

        s_w1 = 0.5*(sol_centroid.x+0.5-F[1])*(s_ipy.y+0.5)
        s_w3 = svol - s_w1
        s_w4 = 0.5*(F[2]-(sol_centroid.y+0.5))*(1.0-(s_ipx.x+0.5))
        s_w2 = svol - s_w4

        if isnan(l_ipy.y) || isinf(l_ipy.y) || abs(l_ipy.y) >= 0.5
            l_w1 = liq_centroid.x+0.5
            l_w3 = lvol - l_w1
        else
            l_w3 = 0.5*(1.0-(liq_centroid.x+0.5))*(1.0-(l_ipy.y+0.5)+1.0-F[2])
            l_w1 = lvol - l_w3
        end
        if isnan(l_ipx.x) || isinf(l_ipx.x) || abs(l_ipx.x) >= 0.5
            l_w4 = 1.0 - (liq_centroid.y+0.5)
            l_w2 = lvol - l_w4
        else
            l_w2 = 0.5*(liq_centroid.y+0.5)*(l_ipx.x+0.5+F[1])
            l_w4 = lvol - l_w2
        end

        cap_sol = SA_F64[0.0, 1 - F[1], F[2], 0.0, svol, bsol..., s_w1, s_w2, s_w3, s_w4]
        cap_liq = SA_F64[1.0, F[1], 1.0 - F[2], 1.0, lvol, bliq..., l_w1, l_w2, l_w3, l_w4]
    elseif case == 3.0 #W E
        s3 = (Point(0.,0.), Point(1.,0.), Point(1.,F[2]), Point(0.,F[1]), Point(0.,0.))
        l3 = (Point(1.,1.), Point(0.,1.), Point(0.,F[1]), Point(1.,F[2]), Point(1.,1.))
        sol_centroid = get_centroid(s3)
        liq_centroid = get_centroid(l3)
        α = atan(1. , F[1]-F[2])
        mid_point = midpoint(Point(1.,F[2]), Point(0.,F[1]))
        cut_points = [Point(1.,F[2]) - Point(0.5,0.5), Point(0.,F[1]) - Point(0.5,0.5)]
        bsol, bliq, s_ipx, s_ipy, l_ipx, l_ipy = Bcapacities(cut_points, sol_centroid, liq_centroid)
        svol = 0.5*(F[1]+F[2])
        lvol = 1.0 - 0.5*(F[1]+F[2])

        s_w1 = 0.5*(sol_centroid.x+0.5)*(F[1]+s_ipy.y+0.5)
        s_w3 = svol - s_w1
        l_w1 = 0.5*(liq_centroid.x+0.5)*(1.0-F[1]+1.0-(l_ipy.y+0.5))
        l_w3 = lvol - l_w1

        if isnan(s_ipx.x) || isinf(s_ipx.x) || abs(s_ipx.x) >= 0.5
            s_w2 = sol_centroid.y+0.5
            s_w4 = svol - s_w2
        elseif F[2] >= F[1]
            s_w4 = 0.5*(1.0-(s_ipx.x+0.5))*(F[2]-(sol_centroid.y+0.5))
            s_w2 = svol - s_w4
        else
            s_w4 = 0.5*(s_ipx.x+0.5)*(F[1]-(sol_centroid.y+0.5))
            s_w2 = svol - s_w4
        end
        if isnan(l_ipx.x) || isinf(l_ipx.x) || abs(l_ipx.x) >= 0.5
            l_w4 = 1.0 - (liq_centroid.y+0.5)
            l_w2 = lvol - l_w4
        elseif F[2] >= F[1]
            l_w2 = 0.5*(liq_centroid.y+0.5-F[1])*(l_ipx.x+0.5)
            l_w4 = lvol - l_w2
        else
            l_w2 = 0.5*(liq_centroid.y+0.5-F[2])*(1.0-(l_ipx.x+0.5))
            l_w4 = lvol - l_w2
        end

        cap_sol = SA_F64[F[1], 1.0, F[2], 0.0, svol, bsol..., s_w1, s_w2, s_w3, s_w4]
        cap_liq = SA_F64[1.0 - F[1], 0.0, 1.0 - F[2], 1.0, lvol, bliq..., l_w1, l_w2, l_w3, l_w4]
    elseif case == 4.0 #E N
        s4 = (Point(1.,1.), Point(F[2],1.), Point(1.,F[1]), Point(1.,1.))
        l4 = (Point(0.,0.), Point(1.,0.), Point(1.,F[1]), Point(F[2],1.), Point(0.,1.), Point(0.,0.))
        sol_centroid = get_centroid(s4)
        liq_centroid = get_centroid(l4)
        α = atan(-1+F[2], -1 +F[1])
        mid_point = midpoint(Point(F[2],1.), Point(1.,F[1]))
        cut_points = [Point(F[2],1.) - Point(0.5,0.5), Point(1.,F[1]) - Point(0.5,0.5)]
        bsol, bliq, s_ipx, s_ipy, l_ipx, l_ipy = Bcapacities(cut_points, sol_centroid, liq_centroid)
        svol = 0.5*(1. - F[1])*(1. - F[2])
        lvol = 1. - 0.5*(1. - F[1])*(1. - F[2])

        s_w1 = 0.5*(sol_centroid.x+0.5-F[2])*(1.0-(s_ipy.y+0.5))
        s_w3 = svol - s_w1
        s_w2 = 0.5*(sol_centroid.y+0.5-F[1])*(1.0-(s_ipx.x+0.5))
        s_w4 = svol - s_w2

        if isnan(l_ipy.y) || isinf(l_ipy.y) || abs(l_ipy.y) >= 0.5
            l_w1 = liq_centroid.x+0.5
            l_w3 = lvol - l_w1
        else
            l_w3 = 0.5*(1.0-(liq_centroid.x+0.5))*(l_ipy.y+0.5+F[1])
            l_w1 = lvol - l_w3
        end
        if isnan(l_ipx.x) || isinf(l_ipx.x) || abs(l_ipx.x) >= 0.5
            l_w2 = liq_centroid.y+0.5
            l_w4 = lvol - l_w2
        else
            l_w4 = 0.5*(1.0-(liq_centroid.y+0.5))*(l_ipx.x+0.5+F[2])
            l_w2 = lvol - l_w4
        end

        cap_liq = SA_F64[1.0, 1.0, F[1], F[2], lvol, bliq..., l_w1, l_w2, l_w3, l_w4]
        cap_sol = SA_F64[0.0, 0.0, 1. - F[1], 1. - F[2], svol, bsol..., s_w1, s_w2, s_w3, s_w4]
    elseif case == 6.0 #S N
        s6 = (Point(1.,1.), Point(F[2],1.), Point(F[1],0.), Point(1.,0.), Point(1.,1.))
        l6 = (Point(0.,0.), Point(F[1],0.), Point(F[2],1.), Point(0.,1.), Point(0.,0.))
        sol_centroid = get_centroid(s6)
        liq_centroid = get_centroid(l6)
        α = atan(F[2]-F[1], -1.)
        mid_point = midpoint(Point(F[2],1.), Point(F[1],0.))
        cut_points = [Point(F[2],1.) - Point(0.5,0.5), Point(F[1],0.) - Point(0.5,0.5)]
        bsol, bliq, s_ipx, s_ipy, l_ipx, l_ipy = Bcapacities(cut_points, sol_centroid, liq_centroid)
        svol = 1.0 - 0.5*(F[1]+F[2])
        lvol = 0.5*(F[1]+F[2])

        s_w2 = 0.5*(sol_centroid.y+0.5)*(1.0-F[1]+1.0-(s_ipx.x+0.5))
        s_w4 = svol - s_w2
        l_w2 = 0.5*(liq_centroid.y+0.5)*(F[1]+l_ipx.x+0.5)
        l_w4 = lvol - l_w2

        if isnan(s_ipy.y) || isinf(s_ipy.y) || abs(s_ipy.y) >= 0.5
            s_w3 = 1.0-(sol_centroid.x+0.5)
            s_w1 = svol - s_w3
        elseif F[2] >= F[1]
            s_w1 = 0.5*(s_ipy.y+0.5)*(sol_centroid.x+0.5-F[1])
            s_w3 = svol - s_w1
        else
            s_w1 = 0.5*(1.0-(s_ipy.y+0.5))*(sol_centroid.x+0.5-F[2])
            s_w3 = svol - s_w1
        end
        if isnan(l_ipy.y) || isinf(l_ipy.y) || abs(l_ipy.y) >= 0.5
            l_w1 = liq_centroid.x+0.5
            l_w3 = lvol - l_w1
        elseif F[2] >= F[1]
            l_w3 = 0.5*(F[2]-(liq_centroid.x+0.5))*(1.0-(l_ipy.y+0.5))
            l_w1 = lvol - l_w3
        else
            l_w3 = 0.5*(F[1]-(liq_centroid.x+0.5))*(l_ipy.y+0.5)
            l_w1 = lvol - l_w3
        end
         
        cap_liq = SA_F64[1.0, F[1], 0.0, F[2], lvol, bliq..., l_w1, l_w2, l_w3, l_w4]
        cap_sol = SA_F64[0.0, 1. - F[1], 1.0, 1. - F[2], svol, bsol..., s_w1, s_w2, s_w3, s_w4]
    elseif case == 7.0 #W N
        s7 = (Point(0.,0.), Point(1.,0.), Point(1.,1.), Point(F[2],1.), Point(0.,F[1]), Point(0.,0.))
        l7 = (Point(0.,1.), Point(0.,F[1]), Point(F[2],1.), Point(0.,1.))
        sol_centroid = get_centroid(s7)
        liq_centroid = get_centroid(l7)
        α = atan(F[2], -1. + F[1])
        mid_point = midpoint(Point(F[2],1.), Point(0.,F[1]))
        cut_points = [Point(F[2],1.) - Point(0.5,0.5), Point(0.,F[1]) - Point(0.5,0.5)]
        bsol, bliq, s_ipx, s_ipy, l_ipx, l_ipy = Bcapacities(cut_points, sol_centroid, liq_centroid)
        svol = 1. - 0.5*(1.0-F[1])*F[2]
        lvol = 0.5*(1.0-F[1])*F[2]

        l_w3 = 0.5*(F[2]-(liq_centroid.x+0.5))*(1.0-(l_ipy.y+0.5))
        l_w1 = lvol - l_w3
        l_w2 = 0.5*(liq_centroid.y+0.5-F[1])*(l_ipx.x+0.5)
        l_w4 = lvol - l_w2

        if isnan(s_ipy.y) || isinf(s_ipy.y) || abs(s_ipy.y) >= 0.5
            s_w3 = 1.0-(sol_centroid.x+0.5)
            s_w1 = svol - s_w3
        else
            s_w1 = 0.5*(sol_centroid.x+0.5)*(s_ipy.y+0.5+F[1])
            s_w3 = svol - s_w1
        end
        if isnan(s_ipx.x) || isinf(s_ipx.x) || abs(s_ipx.x) >= 0.5
            s_w2 = sol_centroid.y+0.5
            s_w4 = svol - s_w2
        else
            s_w4 = 0.5*(1.0-(sol_centroid.y+0.5))*(1.0-(s_ipx.x+0.5)+1.0-F[2])
            s_w2 = svol - s_w4
        end

        cap_sol = SA_F64[F[1], 1.0, 1.0, 1. - F[2], svol, bsol..., s_w1, s_w2, s_w3, s_w4]
        cap_liq = SA_F64[1.0 - F[1], 0.0, 0.0, F[2], lvol, bliq..., l_w1, l_w2, l_w3, l_w4]
    elseif case == 8.0 #W N
        s8 = (Point(0.,1.), Point(0.,F[1]), Point(F[2],1.), Point(0.,1.))
        l8 = (Point(0.,0.), Point(1.,0.), Point(1.,1.), Point(F[2],1.), Point(0.,F[1]), Point(0.,0.))
        sol_centroid = get_centroid(s8)
        liq_centroid = get_centroid(l8)
        α = atan(-F[2], 1. - F[1])
        mid_point = midpoint(Point(0.,F[1]), Point(F[2],1.))
        cut_points = [Point(0.,F[1]) - Point(0.5,0.5), Point(F[2],1.) - Point(0.5,0.5)]
        bsol, bliq, s_ipx, s_ipy, l_ipx, l_ipy = Bcapacities(cut_points, sol_centroid, liq_centroid)
        svol = 0.5*F[2]*(1. -F[1])
        lvol = 1. - 0.5*F[2]*(1.0-F[1])

        s_w3 = 0.5*(F[2]-(sol_centroid.x+0.5))*(1.0-(s_ipy.y+0.5))
        s_w1 = svol - s_w3
        s_w2 = 0.5*(sol_centroid.y+0.5-F[1])*(s_ipx.x+0.5)
        s_w4 = svol - s_w2

        if isnan(l_ipy.y) || isinf(l_ipy.y) || abs(l_ipy.y) >= 0.5
            l_w3 = 1.0-(liq_centroid.x+0.5)
            l_w1 = lvol - l_w3
        else
            l_w1 = 0.5*(liq_centroid.x+0.5)*(l_ipy.y+0.5+F[1])
            l_w3 = lvol - l_w1
        end
        if isnan(l_ipx.x) || isinf(l_ipx.x) || abs(l_ipx.x) >= 0.5
            l_w2 = liq_centroid.y+0.5
            l_w4 = lvol - l_w2
        else
            l_w4 = 0.5*(1.0-(liq_centroid.y+0.5))*(1.0-(l_ipx.x+0.5)+1.0-F[2])
            l_w2 = lvol - l_w4
        end

        cap_sol = SA_F64[1. - F[1], 0.0, 0.0, F[2], svol, bsol..., s_w1, s_w2, s_w3, s_w4]
        cap_liq = SA_F64[F[1], 1.0, 1.0, 1. - F[2], lvol, bliq..., l_w1, l_w2, l_w3, l_w4]
    elseif case == 9.0 #S N
        s9 = (Point(0.,0.), Point(F[1],0.), Point(F[2],1.), Point(0.,1.), Point(0.,0.))
        l9 = (Point(1.,1.), Point(F[2],1.), Point(F[1],0.), Point(1.,0.), Point(1.,1.))
        sol_centroid = get_centroid(s9)
        liq_centroid = get_centroid(l9)
        α = atan(F[1]-F[2], 1.)
        mid_point = midpoint(Point(F[1],0.), Point(F[2],1.))
        cut_points = [Point(F[1],0.) - Point(0.5,0.5), Point(F[2],1.) - Point(0.5,0.5)]
        bsol, bliq, s_ipx, s_ipy, l_ipx, l_ipy = Bcapacities(cut_points, sol_centroid, liq_centroid)
        svol = 0.5*(F[1]+F[2])
        lvol = 1. - 0.5*(F[1]+F[2])

        l_w2 = 0.5*(liq_centroid.y+0.5)*(1.0-F[1]+1.0-(l_ipx.x+0.5))
        l_w4 = lvol - l_w2
        s_w2 = 0.5*(sol_centroid.y+0.5)*(F[1]+s_ipx.x+0.5)
        s_w4 = svol - s_w2

        if isnan(l_ipy.y) || isinf(l_ipy.y) || abs(l_ipy.y) >= 0.5
            l_w3 = 1.0-(liq_centroid.x+0.5)
            l_w1 = lvol - l_w3
        elseif F[2] >= F[1]
            l_w1 = 0.5*(l_ipy.y+0.5)*(liq_centroid.x+0.5-F[1])
            l_w3 = lvol - l_w1
        else
            l_w1 = 0.5*(1.0-(l_ipy.y+0.5))*(liq_centroid.x+0.5-F[2])
            l_w3 = lvol - l_w1
        end
        if isnan(s_ipy.y) || isinf(s_ipy.y) || abs(s_ipy.y) >= 0.5
            s_w1 = sol_centroid.x+0.5
            s_w3 = svol - s_w1
        elseif F[2] >= F[1]
            s_w3 = 0.5*(F[2]-(sol_centroid.x+0.5))*(1.0-(s_ipy.y+0.5))
            s_w1 = svol - s_w3
        else
            s_w3 = 0.5*(F[1]-(sol_centroid.x+0.5))*(s_ipy.y+0.5)
            s_w1 = svol - s_w3
        end

        cap_sol = SA_F64[1.0, F[1], 0.0, F[2], svol, bsol..., s_w1, s_w2, s_w3, s_w4]
        cap_liq = SA_F64[0, 1.0 - F[1], 1.0, 1.0 - F[2], lvol, bliq..., l_w1, l_w2, l_w3, l_w4]
    elseif case == 11.0 #E N
        s11 = (Point(0.,0.), Point(1.,0.), Point(1.,F[1]), Point(F[2],1.), Point(0.,1.), Point(0.,0.))
        l11 = (Point(1.,1.), Point(F[2],1.), Point(1.,F[1]), Point(1.,1.))
        sol_centroid = get_centroid(s11)
        liq_centroid = get_centroid(l11)
        α = atan(1. - F[2], 1. - F[1])
        mid_point = midpoint(Point(1.,F[1]), Point(F[2],1.))
        cut_points = [Point(1.,F[1]) - Point(0.5,0.5), Point(F[2],1.) - Point(0.5,0.5)]
        bsol, bliq, s_ipx, s_ipy, l_ipx, l_ipy = Bcapacities(cut_points, sol_centroid, liq_centroid)
        svol = 1. - 0.5*((1.0 - F[1])*(1.0 - F[2]))
        lvol = 0.5*((1.0 - F[1])*(1.0 - F[2]))

        l_w1 = 0.5*(liq_centroid.x+0.5-F[2])*(1.0-(l_ipy.y+0.5))
        l_w3 = lvol - l_w1
        l_w2 = 0.5*(liq_centroid.y+0.5-F[1])*(1.0-(l_ipx.x+0.5))
        l_w4 = lvol - l_w2

        if isnan(s_ipy.y) || isinf(s_ipy.y) || abs(s_ipy.y) >= 0.5
            s_w1 = sol_centroid.x+0.5
            s_w3 = svol - s_w1
        else
            s_w3 = 0.5*(1.0-(sol_centroid.x+0.5))*(s_ipy.y+0.5+F[1])
            s_w1 = svol - s_w3
        end
        if isnan(s_ipx.x) || isinf(s_ipx.x) || abs(s_ipx.x) >= 0.5
            s_w2 = sol_centroid.y+0.5
            s_w4 = svol - s_w2
        else
            s_w4 = 0.5*(1.0-(sol_centroid.y+0.5))*(s_ipx.x+0.5+F[2])
            s_w2 = svol - s_w4
        end

        cap_sol = SA_F64[1.0, 1.0, F[1], F[2], svol, bsol..., s_w1, s_w2, s_w3, s_w4]
        cap_liq = SA_F64[0.0, 0.0, 1. - F[1], 1. - F[2], lvol, bliq..., l_w1, l_w2, l_w3, l_w4]
    elseif case == 12.0 #W E
        s12 = (Point(1.,1.), Point(0.,1.), Point(0.,F[1]), Point(1.,F[2]), Point(1.,1.))
        l12 = (Point(0.,0.), Point(1.,0.), Point(1.,F[2]), Point(0.,F[1]), Point(0.,0.))
        sol_centroid = get_centroid(s12)
        liq_centroid = get_centroid(l12)
        α = atan(-1., (1. - F[1])-(1. - F[2]))
        mid_point = midpoint(Point(0.,F[1]), Point(1.,F[2]))
        cut_points = [Point(0.,F[1]) - Point(0.5,0.5), Point(1.,F[2]) - Point(0.5,0.5)]
        bsol, bliq, s_ipx, s_ipy, l_ipx, l_ipy = Bcapacities(cut_points, sol_centroid, liq_centroid)
        svol = 1. - 0.5*(F[1]+F[2])
        lvol = 0.5*(F[1]+F[2])

        l_w1 = 0.5*(liq_centroid.x+0.5)*(F[1]+l_ipy.y+0.5)
        l_w3 = lvol - l_w1
        s_w1 = 0.5*(sol_centroid.x+0.5)*(1.0-F[1]+1.0-(s_ipy.y+0.5))
        s_w3 = svol - s_w1

        if isnan(l_ipx.x) || isinf(l_ipx.x) || abs(l_ipx.x) >= 0.5
            l_w2 = liq_centroid.y+0.5
            l_w4 = lvol - l_w2
        elseif F[2] >= F[1]
            l_w4 = 0.5*(1.0-(l_ipx.x+0.5))*(F[2]-(liq_centroid.y+0.5))
            l_w2 = lvol - l_w4
        else
            l_w4 = 0.5*(l_ipx.x+0.5)*(F[1]-(liq_centroid.y+0.5))
            l_w2 = lvol - l_w4
        end
        if isnan(s_ipx.x) || isinf(s_ipx.x) || abs(s_ipx.x) >= 0.5
            s_w4 = 1.0 - (sol_centroid.y+0.5)
            s_w2 = svol - s_w4
        elseif F[2] >= F[1]
            s_w2 = 0.5*(sol_centroid.y+0.5-F[1])*(s_ipx.x+0.5)
            s_w4 = svol - s_w2
        else
            s_w2 = 0.5*(sol_centroid.y+0.5-F[2])*(1.0-(s_ipx.x+0.5))
            s_w4 = svol - s_w2
        end

        cap_sol = SA_F64[1. - F[1], 0.0, 1. - F[2], 1.0, svol, bsol..., s_w1, s_w2, s_w3, s_w4]
        cap_liq = SA_F64[F[1], 1.0, F[2], 0.0, lvol, bliq..., l_w1, l_w2, l_w3, l_w4]
    elseif case == 13.0 #S E
        s13 = (Point(0.,0.), Point(F[1],0.), Point(1.,F[2]), Point(1.,1.), Point(0.,1.), Point(0.,0.))
        l13 = (Point(1.,0.), Point(1.,F[2]), Point(F[1],0.), Point(1.,0.))
        sol_centroid = get_centroid(s13)
        liq_centroid = get_centroid(l13)
        α = atan(-1.0+F[1], F[2])
        mid_point = midpoint(Point(F[1],0.), Point(1.,F[2]))
        cut_points = [Point(F[1],0.) - Point(0.5,0.5), Point(1.,F[2]) - Point(0.5,0.5)]
        bsol, bliq, s_ipx, s_ipy, l_ipx, l_ipy = Bcapacities(cut_points, sol_centroid, liq_centroid)
        svol = 1. - 0.5*(F[2]*(1. - F[1]))
        lvol = 0.5*(F[2]*(1. - F[1]))

        l_w1 = 0.5*(liq_centroid.x+0.5-F[1])*(l_ipy.y+0.5)
        l_w3 = lvol - l_w1
        l_w4 = 0.5*(F[2]-(liq_centroid.y+0.5))*(1.0-(l_ipx.x+0.5))
        l_w2 = lvol - l_w4

        if isnan(s_ipy.y) || isinf(s_ipy.y) || abs(s_ipy.y) >= 0.5
            s_w1 = sol_centroid.x+0.5
            s_w3 = svol - s_w1
        else
            s_w3 = 0.5*(1.0-(sol_centroid.x+0.5))*(1.0-(s_ipy.y+0.5)+1.0-F[2])
            s_w1 = svol - s_w3
        end
        if isnan(s_ipx.x) || isinf(s_ipx.x) || abs(s_ipx.x) >= 0.5
            s_w4 = 1.0 - (sol_centroid.y+0.5)
            s_w2 = svol - s_w4
        else
            s_w2 = 0.5*(sol_centroid.y+0.5)*(s_ipx.x+0.5+F[1])
            s_w4 = svol - s_w2
        end

        cap_sol = SA_F64[1.0, F[1], 1. - F[2], 1.0, svol, bsol..., s_w1, s_w2, s_w3, s_w4]
        cap_liq = SA_F64[0.0, 1.0 - F[1], F[2], 0.0, lvol, bliq..., l_w1, l_w2, l_w3, l_w4]
    elseif case == 14.0 #W S
        s14 = (Point(1.,1.), Point(0.,1.), Point(0.,F[1]), Point(F[2],0.), Point(1.,0.), Point(1.,1.))
        l14 = (Point(0.,0.), Point(F[2],0.), Point(0.,F[1]), Point(0.,0.))
        sol_centroid = get_centroid(s14)
        liq_centroid = get_centroid(l14)
        α = atan(-F[2], -F[1])
        mid_point = midpoint(Point(0.,F[1]), Point(F[2],0.))
        cut_points = [Point(0.,F[1]) - Point(0.5,0.5), Point(F[2],0.) - Point(0.5,0.5)]
        bsol, bliq, s_ipx, s_ipy, l_ipx, l_ipy = Bcapacities(cut_points, sol_centroid, liq_centroid)
        svol = 1. - 0.5*F[1]*F[2]
        lvol = 0.5*F[1]*F[2]

        l_w3 = 0.5*(F[2]-(liq_centroid.x+0.5))*(l_ipy.y+0.5)
        l_w1 = lvol - l_w3
        l_w4 = 0.5*(F[1]-(liq_centroid.y+0.5))*(l_ipx.x+0.5)
        l_w2 = lvol - l_w4

        if isnan(s_ipy.y) || isinf(s_ipy.y) || abs(s_ipy.y) >= 0.5
            s_w3 = 1.0 - (sol_centroid.x+0.5)
            s_w1 = svol - s_w3
        else
            s_w1 = 0.5*(sol_centroid.x+0.5)*(1.0-(s_ipy.y+0.5)+1.0-F[1])
            s_w3 = svol - s_w1
        end
        if isnan(s_ipx.x) || isinf(s_ipx.x) || abs(s_ipx.x) >= 0.5
            s_w4 = 1.0 - (sol_centroid.y+0.5)
            s_w2 = svol - s_w4
        else
            s_w2 = 0.5*(sol_centroid.y+0.5)*(1.0-(s_ipx.x+0.5)+1.0-F[2])
            s_w4 = svol - s_w2
        end

        cap_sol = SA_F64[1.0 - F[1], 1.0 - F[2], 1.0, 1.0, svol, bsol..., s_w1, s_w2, s_w3, s_w4]
        cap_liq = SA_F64[F[1], F[2], 0.0, 0.0, lvol, bliq..., l_w1, l_w2, l_w3, l_w4]
    end
    return cap_sol, cap_liq, min(float(π),max(-float(π),α)), sol_centroid, liq_centroid, mid_point, cut_points
end

function set_cap_bcs!(grid::Mesh{GridCC,T,N}, periodic_x, periodic_y) where {T,N}
    @unpack nx, ny, geoS, geoL, mid_point = grid
    # set A and mid_point at the boundaries to 0 if not periodic in that direction
    if periodic_x
        @inbounds @threads for i = 1:ny
            @inbounds geoS.dcap[i,1,1] = geoS.dcap[i,end,3]
            @inbounds geoL.dcap[i,1,1] = geoL.dcap[i,end,3]
        end
    else
        @inbounds @threads for i = 1:ny
            @inbounds geoS.dcap[i,1,1] = 0.0
            @inbounds geoS.dcap[i,end,3] = 0.0
            @inbounds geoL.dcap[i,1,1] = 0.0
            @inbounds geoL.dcap[i,end,3] = 0.0
            @inbounds mid_point[i,1] += Point(-0.5, 0.0)
            @inbounds mid_point[i,end] += Point(0.5, 0.0)
        end
    end
    if periodic_y
        @inbounds @threads for i = 1:nx
            @inbounds geoS.dcap[1,i,2] = geoS.dcap[end,i,4]
            @inbounds geoL.dcap[1,i,2] = geoL.dcap[end,i,4]
        end
    else
        @inbounds @threads for i = 1:nx
            @inbounds geoS.dcap[1,i,2] = 0.0
            @inbounds geoS.dcap[end,i,4] = 0.0
            @inbounds geoL.dcap[1,i,2] = 0.0
            @inbounds geoL.dcap[end,i,4] = 0.0
            @inbounds mid_point[1,i] += Point(0.0, -0.5)
            @inbounds mid_point[end,i] += Point(0.0, 0.5)
        end
    end

    return nothing
end

function set_cap_bcs!(grid::Mesh{GridFCx,T,N}, periodic_x, periodic_y) where {T,N}
    @unpack nx, ny, ind, u, geoS, geoL, mid_point, cut_points = grid
    @unpack b_left, b_bottom, b_right, b_top = ind

    # set A at the boundaries
    if periodic_y
        @inbounds @threads for i = 1:nx
            @inbounds geoS.dcap[1,i,2] = geoS.dcap[end,i,4]
            @inbounds geoL.dcap[1,i,2] = geoL.dcap[end,i,4]
        end
    else
        @inbounds @threads for i = 1:nx
            @inbounds geoS.dcap[1,i,2] = 0.0
            @inbounds geoS.dcap[end,i,4] = 0.0
            @inbounds geoL.dcap[1,i,2] = 0.0
            @inbounds geoL.dcap[end,i,4] = 0.0
            @inbounds mid_point[1,i] += Point(0.0, -0.5)
            @inbounds mid_point[end,i] += Point(0.0, 0.5)
        end
    end

    # set mid_point in the outer boundaries
    if periodic_x
        @inbounds @threads for i = 1:ny
            @inbounds geoS.dcap[i,1,1] = geoS.dcap[i,end,3]
            @inbounds geoL.dcap[i,1,1] = geoL.dcap[i,end,3]
        end
    end

    return nothing
end

function set_cap_bcs!(grid::Mesh{GridFCy,T,N}, periodic_x, periodic_y) where {T,N}
    @unpack nx, ny, ind, u, geoS, geoL, mid_point, cut_points = grid
    @unpack b_left, b_bottom, b_right, b_top = ind

    # set A at the boundaries to 0 if not periodic in that direction
    if periodic_x
        @inbounds @threads for i = 1:ny
            @inbounds geoS.dcap[i,1,1] = geoS.dcap[i,end,3]
            @inbounds geoL.dcap[i,1,1] = geoL.dcap[i,end,3]
        end
    else
        @inbounds @threads for i = 1:ny
            @inbounds geoS.dcap[i,1,1] = 0.0
            @inbounds geoS.dcap[i,end,3] = 0.0
            @inbounds geoL.dcap[i,1,1] = 0.0
            @inbounds geoL.dcap[i,end,3] = 0.0
            @inbounds mid_point[i,1] += Point(-0.5, 0.0)
            @inbounds mid_point[i,end] += Point(0.5, 0.0)
        end
    end

    # set mid_point in the outer boundaries
    if periodic_y
        @inbounds @threads for i = 1:nx
            @inbounds geoS.dcap[1,i,2] = geoS.dcap[end,i,4]
            @inbounds geoL.dcap[1,i,2] = geoL.dcap[end,i,4]
        end
    end

    return nothing
end

function Bcapacities(cut_points, sol_centroid, liq_centroid)
    l_int = Line(cut_points[1], cut_points[2])

    ly = Line(sol_centroid, sol_centroid+Point(0.0,1.0))
    lx = Line(sol_centroid, sol_centroid+Point(1.0,0.0))
    s_ipy = line_intersection(l_int, ly)
    s_ipx = line_intersection(l_int, lx)

    sol = [Bcaps(s_ipy.y, sol_centroid.y), Bcaps(s_ipx.x, sol_centroid.x)]

    ly = Line(liq_centroid, liq_centroid+Point(0.0,1.0))
    lx = Line(liq_centroid, liq_centroid+Point(1.0,0.0))
    l_ipy = line_intersection(l_int, ly)
    l_ipx = line_intersection(l_int, lx)

    liq = [Bcaps(l_ipy.y, liq_centroid.y), Bcaps(l_ipx.x, liq_centroid.x)]

    return sol, liq, s_ipx, s_ipy, l_ipx, l_ipy
end

function Wcapacities!(cap, periodic_x, periodic_y)
    tmp = copy(cap)

    @inbounds cap[:,2:end,8] .= @views tmp[:,2:end,8] .+ tmp[:,1:end-1,10]
    @inbounds cap[:,1:end-1,10] .= @views tmp[:,1:end-1,10] .+ tmp[:,2:end,8]
    @inbounds cap[2:end,:,9] .= @views tmp[2:end,:,9] .+ tmp[1:end-1,:,11]
    @inbounds cap[1:end-1,:,11] .= @views tmp[1:end-1,:,11] .+ tmp[2:end,:,9]

    if periodic_x
        @inbounds cap[:,1,8] .= @views tmp[:,1,8] .+ tmp[:,end,10]
        @inbounds cap[:,end,10] .= @views tmp[:,end,10] .+ tmp[:,1,8]
    end
    if periodic_y
        @inbounds cap[1,:,9] .= @views tmp[1,:,9] .+ tmp[end,:,11]
        @inbounds cap[end,:,11] .= @views tmp[end,:,11] .+ tmp[1,:,9]
    end

    return nothing
end

function face_capacities(grid, itp, case, II_0, II, posW, posS, posE, posN)
    @unpack dx, dy, faces = grid

    if case == 1.0 || case == 14.0
        faces[II, 1] = ispositive(WEST_face(itp, grid, II_0, posW, dy[II])) / dy[II]
        faces[II, 2] = ispositive(SOUTH_face(itp, grid, II_0, posS, dx[II])) / dx[II]
    elseif case == 2.0 || case == 13.0
        faces[II, 2] = ispositive(SOUTH_face(itp, grid, II_0, posS, dx[II])) / dx[II]
        faces[II, 3] = ispositive(EAST_face(itp, grid, II_0, posE, dy[II])) / dy[II]
    elseif case == 3.0 || case == 12.0
        faces[II, 1] = ispositive(WEST_face(itp, grid, II_0, posW, dy[II])) / dy[II]
        faces[II, 3] = ispositive(EAST_face(itp, grid, II_0, posE, dy[II])) / dy[II]
    elseif case == 4.0 || case == 11.0
        faces[II, 3] = ispositive(EAST_face(itp, grid, II_0, posE, dy[II])) / dy[II]
        faces[II, 4] = ispositive(NORTH_face(itp, grid, II_0, posN, dx[II])) / dx[II]
    elseif case == 6.0 || case == 9.0
        faces[II, 2] = ispositive(SOUTH_face(itp, grid, II_0, posS, dx[II])) / dx[II]
        faces[II, 4] = ispositive(NORTH_face(itp, grid, II_0, posN, dx[II])) / dx[II]
    elseif case == 7.0 || case == 8.0
        faces[II, 1] = ispositive(WEST_face(itp, grid, II_0, posW, dy[II])) / dy[II]
        faces[II, 4] = ispositive(NORTH_face(itp, grid, II_0, posN, dx[II])) / dx[II]
    end
end

function average_face_capacities(grid, II)
    @unpack ind, iso, faces = grid
    @unpack b_left, b_bottom, b_right, b_top = ind

    case = iso[II]

    if II in b_left[1]
        WE_fun = WE_border
    else
        WE_fun = WE
    end
    if II in b_bottom[1]
        SN_fun = SN_border
    else
        SN_fun = SN
    end
    if II in b_right[1]
        EW_fun = EW_border
    else
        EW_fun = EW
    end
    if II in b_top[1]
        NS_fun = NS_border
    else
        NS_fun = NS
    end
    if case == 1.0
        f = SA_F64[WE_fun(faces, iso, II, 0.0, 1.0), SN_fun(faces, iso, II, 0.0, 1.0)]
    elseif case == 14.0
        f = SA_F64[WE_fun(faces, iso, II, 1.0, 0.0), SN_fun(faces, iso, II, 1.0, 0.0)]
    elseif case == 2.0
        f = SA_F64[SN_fun(faces, iso, II, 1.0, 0.0), EW_fun(faces, iso, II, 0.0, 1.0)]
    elseif case == 13.0
        f = SA_F64[SN_fun(faces, iso, II, 0.0, 1.0), EW_fun(faces, iso, II, 1.0, 0.0)]
    elseif case == 3.0
        f = SA_F64[WE_fun(faces, iso, II, 0.0, 1.0), EW_fun(faces, iso, II, 0.0, 1.0)]
    elseif case == 12.0
        f = SA_F64[WE_fun(faces, iso, II, 1.0, 0.0), EW_fun(faces, iso, II, 1.0, 0.0)]
    elseif case == 4.0
        f = SA_F64[EW_fun(faces, iso, II, 1.0, 0.0), NS_fun(faces, iso, II, 1.0, 0.0)]
    elseif case == 11.0
        f = SA_F64[EW_fun(faces, iso, II, 0.0, 1.0), NS_fun(faces, iso, II, 0.0, 1.0)]
    elseif case == 6.0
        f = SA_F64[SN_fun(faces, iso, II, 1.0, 0.0), NS_fun(faces, iso, II, 1.0, 0.0)]
    elseif case == 9.0
        f = SA_F64[SN_fun(faces, iso, II, 0.0, 1.0), NS_fun(faces, iso, II, 0.0, 1.0)]
    elseif case == 7.0
        f = SA_F64[WE_fun(faces, iso, II, 0.0, 1.0), NS_fun(faces, iso, II, 1.0, 0.0)]
    elseif case == 8.0
        f = SA_F64[WE_fun(faces, iso, II, 1.0, 0.0), NS_fun(faces, iso, II, 0.0, 1.0)]
    else
        f = @SVector zeros(2)
    end
    return float(f)
end

function average_face_capacities!(a)
    # average face capacities due to mismatch in iso 
    n = size(a,1)
    @inbounds @threads for i = 1:n
        @inbounds tmp1 = @view a[i,1:end-1,3]
        @inbounds tmp2 = @view a[i,2:end,1]
        @inbounds avg_cap = 0.5*(tmp1 .+ tmp2)
        @inbounds a[i,1:end-1,3] .= avg_cap
        @inbounds a[i,2:end,1] .= avg_cap
    end
    n = size(a,2)
    @inbounds @threads for i = 1:n
        @inbounds tmp1 = @view a[1:end-1,i,4]
        @inbounds tmp2 = @view a[2:end,i,2]
        @inbounds avg_cap = 0.5*(tmp1 .+ tmp2)
        @inbounds a[1:end-1,i,4] .= avg_cap
        @inbounds a[2:end,i,2] .= avg_cap
    end
end

function get_cells_indices(iso, all)
    local M = 1
    local S = 1
    local L = 1
    MIXED = Vector{CartesianIndex{2}}(undef, length(all))
    SOLID = Vector{CartesianIndex{2}}(undef, length(all))
    LIQUID = Vector{CartesianIndex{2}}(undef, length(all))
    @simd for II in all
         @inbounds if is_solid(iso[II])
             SOLID[S] = II
             S+=1
         elseif is_liquid(iso[II])
             LIQUID[L] = II
             L+=1
         else
             MIXED[M] = II
             M+=1
         end
    end
    resize!(MIXED, M-1); resize!(SOLID, S-1); resize!(LIQUID, L-1);
    sizehint!(MIXED, length(MIXED)); sizehint!(SOLID, length(SOLID)); sizehint!(LIQUID, length(LIQUID))
    return MIXED, SOLID, LIQUID
end

function get_cells_indices(iso, all, nx, ny, periodic_x, periodic_y)
    local M = 1
    local S = 1
    local L = 1
    MIXED = Vector{CartesianIndex{2}}(undef, length(all))
    SOLID = Vector{CartesianIndex{2}}(undef, length(all))
    LIQUID = Vector{CartesianIndex{2}}(undef, length(all))
    @simd for II in all
        @inbounds if !periodic_x && (II[2] == 1 || II[2] == nx)

        elseif  !periodic_y && (II[1] == 1 || II[1] == ny)
         
        elseif is_solid(iso[II])
            SOLID[S] = II
            S+=1
        elseif is_liquid(iso[II])
            LIQUID[L] = II
            L+=1
        else
            MIXED[M] = II
            M+=1
        end
    end
    resize!(MIXED, M-1); resize!(SOLID, S-1); resize!(LIQUID, L-1);
    sizehint!(MIXED, length(MIXED)); sizehint!(SOLID, length(SOLID)); sizehint!(LIQUID, length(LIQUID))
    return MIXED, SOLID, LIQUID
end

function get_cells_indices(iso, inside, NB_width)
    local M = 1
    local S = 1
    local L = 1
    local N = 1
    MIXED = Vector{CartesianIndex{2}}(undef, length(inside))
    SOLID = Vector{CartesianIndex{2}}(undef, length(inside))
    LIQUID = Vector{CartesianIndex{2}}(undef, length(inside))
    NB_indices = Vector{CartesianIndex{2}}(undef, length(inside))
    for II in inside
         @inbounds if is_solid(iso[II])
             SOLID[S] = II
             S+=1
         elseif is_liquid(iso[II])
             LIQUID[L] = II
             L+=1
         else
             MIXED[M] = II
             M+=1
             NB_indices[N] = II
             N += 1
             for JJ in NB_width
                 NB_indices[N] = II + JJ
                 N += 1
             end
         end
    end
    resize!(MIXED, M-1); resize!(SOLID, S-1); resize!(LIQUID, L-1); resize!(NB_indices, N-1);
    sizehint!(MIXED, length(MIXED)); sizehint!(SOLID, length(SOLID)); sizehint!(LIQUID, length(LIQUID)); sizehint!(NB_indices, length(NB_indices))
    return MIXED, SOLID, LIQUID, unique!(NB_indices)
end

function get_NB_width(MIXED, NB_indices_base)
        NB_indices = Vector{CartesianIndex{2}}(undef, 0)
        for II in MIXED
            for JJ in NB_indices_base
                push!(NB_indices, II + JJ)
            end
        end
        return unique!(NB_indices)
end


@inline function affine(A, B, C)
    if isnan(C.y)
        a = (A.y-B.y)/(A.x-B.x)
        if (abs(a) === Inf) || (abs(a) === NaN) a = 0. end
        b = A.y - A.x*a
        s = muladd(a, C.x, b)
        return Point(C.x, s)
    elseif isnan(C.x)
        a = (A.x-B.x)/(A.y-B.y)
        if (abs(a) === Inf) || (abs(a) === NaN) a = 0. end
        b = A.x - A.y*a
        s = muladd(a, C.y, b)
        return Point(s, C.y)
    end
end

function projection_2points(grid, II)
    @unpack x_nodes, y_nodes, x, y, nx, ny, dx, dy, mid_point, α = grid

    h = Point(dx[II], dy[II])
    ĥ = sqrt(h.x^2 + h.y^2)
    absolute_position = Point(x[II], y[II])
    mid = mid_point[II]*h
    β = opposite(α[II])
    tmpS = Point(mid.x + cos(α[II])*ĥ, mid.y + sin(α[II])*ĥ)
    tmpL = Point(mid.x + cos(β)*ĥ, mid.y + sin(β)*ĥ)
    S2 = Point(NaN, NaN)
    L2 = Point(NaN, NaN)
    Sflag = false
    Lflag = false
    
    if II[2] < 3
        mx = -dx[II]
        m2x = -2*dx[II]
        px = x[δx⁺(II)] - x[II]
        p2x = 2*(x[δx⁺(II)] - x[II])
    elseif II[2] > nx-2
        mx = x[δx⁻(II)] - x[II]
        m2x = 2*(x[δx⁻(II)] - x[II])
        px = dx[II]
        p2x = 2*dx[II]
    else
        mx = x[δx⁻(II)] - x[II]
        m2x = 2*(x[δx⁻(II)] - x[II])
        px = x[δx⁺(II)] - x[II]
        p2x = 2*(x[δx⁺(II)] - x[II])
    end
    if II[1] < 3
        my = -dy[II]
        m2y = -2*dy[II]
        py = y[δy⁺(II)] - y[II]
        p2y = 2*(y[δy⁺(II)] - y[II])
    elseif II[1] > ny-2
        my = y[δy⁻(II)] - y[II]
        m2y = 2*(y[δy⁻(II)] - y[II])
        py = dy[II]
        p2y = 2*dy[II]
    else
        my = y[δy⁻(II)] - y[II]
        m2y = 2*(y[δy⁻(II)] - y[II])
        py = y[δy⁺(II)] - y[II]
        p2y = 2*(y[δy⁺(II)] - y[II])
    end
    
    if -π/4 < α[II] <= π/4
        S1 = affine(mid, tmpS, Point(mx, NaN))
        S2 = affine(mid, tmpS, Point(m2x, NaN))
        L1 = affine(mid, tmpL, Point(px, NaN))
        L2 = affine(mid, tmpL, Point(p2x, NaN))
    elseif π/4 < α[II] <= 3π/4
        S1 = affine(mid, tmpS, Point(NaN, my))
        S2 = affine(mid, tmpS, Point(NaN, m2y))
        L1 = affine(mid, tmpL, Point(NaN, py))
        L2 = affine(mid, tmpL, Point(NaN, p2y))
    elseif α[II] > 3π/4 || α[II] <= -3π/4
        S1 = affine(mid, tmpS, Point(px, NaN))
        S2 = affine(mid, tmpS, Point(p2x, NaN))
        L1 = affine(mid, tmpL, Point(mx, NaN))
        L2 = affine(mid, tmpL, Point(m2x, NaN))
    elseif -3π/4 < α[II] <= -π/4
        S1 = affine(mid, tmpS, Point(NaN, py))
        S2 = affine(mid, tmpS, Point(NaN, p2y))
        L1 = affine(mid, tmpL, Point(NaN, my))
        L2 = affine(mid, tmpL, Point(NaN, m2y))
    end
    if S2 != Point(NaN, NaN) && L2 != Point(NaN, NaN)
        Sflag = indomain(absolute_position + S2, x_nodes, y_nodes)
        Lflag = indomain(absolute_position + L2, x_nodes, y_nodes)
    end
    pos = absolute_position + mid
    return Gradient(Sflag, β, mid, S1, S2, distance(mid, S1), distance(mid, S2), pos), Gradient(Lflag, α[II], mid, L1, L2, distance(mid, L1), distance(mid, L2), pos)
end

function kill_dead_cells!(T::Matrix, L, EMPTY, MIXED, n)
    @inbounds @threads for II in EMPTY
        pII = lexicographic(II, n)
        sumL = abs(sum(abs.(L[:,pII])) - 4.0)
        if sumL <= 1e-8
            T[II] = 0.
        end
    end
    @inbounds @threads for II in MIXED
        pII = lexicographic(II, n)
        sumL = abs(sum(abs.(L[:,pII])) - 4.0)
        if sumL <= 1e-8
            T[II] = 0.
        end
    end
end

function kill_dead_cells!(T::Vector, L, EMPTY, MIXED, n)
    @inbounds @threads for II in EMPTY
        pII = lexicographic(II, n)
        sumL = abs(sum(abs.(L[:,pII])) - 4.0)
        if sumL <= 1e-8
            T[pII] = 0.
        end
    end
    @inbounds @threads for II in MIXED
        pII = lexicographic(II, n)
        sumL = abs(sum(abs.(L[:,pII])) - 4.0)
        if sumL <= 1e-8
            T[pII] = 0.
        end
    end
end

function init_fresh_cells!(grid, T, projection, FRESH, periodic_x, periodic_y)
    @inbounds @threads for II in FRESH
        if projection[II].flag
            T_1, T_2 = interpolated_temperature(grid, projection[II].angle, projection[II].point1, projection[II].point2, T, II, periodic_x, periodic_y)
            if π/4 <= projection[II].angle <= 3π/4 || -π/4 >= projection[II].angle >= -3π/4
                T[II] = y_extrapolation(T_1, T_2, projection[II].point1, projection[II].point2, projection[II].mid_point)
            else
                T[II] = x_extrapolation(T_1, T_2, projection[II].point1, projection[II].point2, projection[II].mid_point)
            end
        end
    end
end

function init_fresh_cells!(grid, T::Vector, projection, FRESH, n, periodic_x, periodic_y)
    @inbounds @threads for II in FRESH
        if projection[II].flag
            pII = lexicographic(II, n)
            T_1, T_2 = interpolated_temperature(grid, projection[II].angle, projection[II].point1, projection[II].point2, T, II, periodic_x, periodic_y)
            if π/4 <= projection[II].angle <= 3π/4 || -π/4 >= projection[II].angle >= -3π/4
                T[pII] = y_extrapolation(T_1, T_2, projection[II].point1, projection[II].point2, projection[II].mid_point)
            else
                T[pII] = x_extrapolation(T_1, T_2, projection[II].point1, projection[II].point2, projection[II].mid_point)
            end
        end
    end
end

function init_fresh_cells!(grid, u::Matrix, V, projection, FRESH, periodic_x, periodic_y)
    @inbounds @threads for II in FRESH
        if projection[II].flag
            u_1, u_2 = interpolated_temperature(grid, projection[II].angle, projection[II].point1, projection[II].point2, V, II, periodic_x, periodic_y)
            if π/4 <= projection[II].angle <= 3π/4 || -π/4 >= projection[II].angle >= -3π/4
                u[II] = y_extrapolation(u_1, u_2, projection[II].point1, projection[II].point2, projection[II].mid_point)
            else
                u[II] = x_extrapolation(u_1, u_2, projection[II].point1, projection[II].point2, projection[II].mid_point)
            end
        end
    end
end

function x_extrapolation(T_1, T_2, p1, p2, pnew)
    a1 = (T_2 - T_1)/(p2.x - p1.x)
    a0 = T_2 - a1*p2.x
    return a1*pnew.x + a0
end

function y_extrapolation(T_1, T_2, p1, p2, pnew)
    a1 = (T_2 - T_1)/(p2.y - p1.y)
    a0 = T_2 - a1*p2.y
    return a1*pnew.y + a0
end
