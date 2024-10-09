@inline SOUTH_face(itp, p=Point(0.0,-0.5), dx=1.0, dx2=2.0, dy2=2.0) = dx/2 - p.x + find_zero(x -> biquadratic(itp, x, p.y/dy2), (-dx/2+p.x,dx/2+p.x)./dx2, FalsePosition(), maxevals = 10, atol = 1e-9) * dx2
@inline WEST_face(itp, p=Point(-0.5,0.0), dy=1.0, dx2=2.0, dy2=2.0) = dy/2 - p.y + find_zero(y -> biquadratic(itp, p.x/dx2, y), (-dy/2+p.y,dy/2+p.y)./dy2, FalsePosition(), maxevals = 10, atol = 1e-9) * dy2
@inline NORTH_face(itp, p=Point(0.0,0.5), dx=1.0, dx2=2.0, dy2=2.0) = dx/2 - p.x + find_zero(x -> biquadratic(itp, x, p.y/dy2), (-dx/2+p.x,dx/2+p.x)./dx2, FalsePosition(), maxevals = 10, atol = 1e-9) * dx2
@inline EAST_face(itp, p=Point(0.5,0.0), dy=1.0, dx2=2.0, dy2=2.0) = dy/2 - p.y + find_zero(y -> biquadratic(itp, p.x/dx2, y), (-dy/2+p.y,dy/2+p.y)./dy2, FalsePosition(), maxevals = 10, atol = 1e-9) * dy2

@inline WE(g, LS, II, l_face, s_face, per) = ismixed(LS.iso[δx⁻(II, g.nx, per)]) ? 0.5*(LS.faces[II,1] + LS.faces[δx⁻(II, g.nx, per), 3]) : (is_liquid(LS.iso[δx⁻(II, g.nx, per)]) ? l_face : s_face)
@inline WE_border(g, LS, II, l_face, s_face, per) = LS.faces[II,1]
@inline EW(g, LS, II, l_face, s_face, per) = ismixed(LS.iso[δx⁺(II, g.nx, per)]) ? 0.5*(LS.faces[II,3] + LS.faces[δx⁺(II, g.nx, per), 1]) : (is_liquid(LS.iso[δx⁺(II, g.nx, per)]) ? l_face : s_face)
@inline EW_border(g, LS, II, l_face, s_face, per) = LS.faces[II,3]
@inline SN(g, LS, II, l_face, s_face, per) = ismixed(LS.iso[δy⁻(II, g.ny, per)]) ? 0.5*(LS.faces[II,2] + LS.faces[δy⁻(II, g.ny, per), 4]) : (is_liquid(LS.iso[δy⁻(II, g.ny, per)]) ? l_face : s_face)
@inline SN_border(g, LS, II, l_face, s_face, per) = LS.faces[II,2]
@inline NS(g, LS, II, l_face, s_face, per) = ismixed(LS.iso[δy⁺(II, g.ny, per)]) ? 0.5*(LS.faces[II,4] + LS.faces[δy⁺(II, g.ny, per), 2]) : (is_liquid(LS.iso[δy⁺(II, g.ny, per)]) ? l_face : s_face)
@inline NS_border(g, LS, II, l_face, s_face, per) = LS.faces[II,4]

@inline WE(a, II) = 0.5*(a[II,1] + a[δx⁻(II), 3])
@inline EW(a, II) = 0.5*(a[II,3] + a[δx⁺(II), 1])
@inline SN(a, II) = 0.5*(a[II,2] + a[δy⁻(II), 4])
@inline NS(a, II) = 0.5*(a[II,4] + a[δy⁺(II), 2])

@inline ispositive(a) = min(max(a,0),1)
@inline isovalue(v::SVector{4,Float64}) = 4. * v[1] + 2. * v[2] + v[3] + 0.5 * v[4]

@inline vertices_sign(itp, g, II_0, II, dx, dy, dx2, dy2) = 
    1 .- sign.(@SVector [biquadratic(itp, (g.x[II] - g.x[II_0] - dx/2)/dx2,
                            (g.y[II] - g.y[II_0] + dy/2)/dy2),
                        biquadratic(itp, (g.x[II] - g.x[II_0] + dx/2)/dx2,
                            (g.y[II] - g.y[II_0] + dy/2)/dy2),
                        biquadratic(itp, (g.x[II] - g.x[II_0] + dx/2)/dx2,
                            (g.y[II] - g.y[II_0] - dy/2)/dy2),
                        biquadratic(itp, (g.x[II] - g.x[II_0] - dx/2)/dx2,
                            (g.y[II] - g.y[II_0] - dy/2)/dy2)])

@inline face_pos(II_0, II, x, y, dx, dy) =
    (Point(x[II] - x[II_0] - dx/2, y[II] - y[II_0]),
    Point(x[II] - x[II_0], y[II] - y[II_0] - dy/2),
    Point(x[II] - x[II_0] + dx/2, y[II] - y[II_0]),
    Point(x[II] - x[II_0], y[II] - y[II_0] + dy/2))

@inline is_near_interface(a, st) = @inbounds ifelse(a != sign(st[2,1]) || a != sign(st[1,2]) || a != sign(st[2,3]) || a != sign(st[3,2]), true, false)
@inline is_near_interface_l(a, st) = @inbounds ifelse(a != sign(st[1,1]) || a != sign(st[3,1]) || a != sign(st[2,2]), true, false)
@inline is_near_interface_b(a, st) = @inbounds ifelse(a != sign(st[1,1]) || a != sign(st[1,3]) || a != sign(st[2,2]), true, false)
@inline is_near_interface_r(a, st) = @inbounds ifelse(a != sign(st[1,3]) || a != sign(st[3,3]) || a != sign(st[2,2]), true, false)
@inline is_near_interface_t(a, st) = @inbounds ifelse(a != sign(st[3,1]) || a != sign(st[3,3]) || a != sign(st[2,2]), true, false)
@inline is_near_interface_bl(a, st) = @inbounds ifelse(a != sign(st[1,2]) || a != sign(st[2,1]), true, false)
@inline is_near_interface_br(a, st) = @inbounds ifelse(a != sign(st[1,2]) || a != sign(st[2,3]), true, false)
@inline is_near_interface_tr(a, st) = @inbounds ifelse(a != sign(st[3,2]) || a != sign(st[2,3]), true, false)
@inline is_near_interface_tl(a, st) = @inbounds ifelse(a != sign(st[3,2]) || a != sign(st[2,1]), true, false)


@inline is_near_interface(u, II::CartesianIndex) = @inbounds ifelse(u[II]*u[δx⁺(II)] < 0 || u[II]*u[δy⁺(II)] < 0 || u[II]*u[δx⁻(II)] < 0 || u[II]*u[δy⁻(II)] < 0, true, false)
@inline is_near_interface(u, II::CartesianIndex, nx, ny, per_x, per_y) = @inbounds ifelse(u[II]*u[δx⁺(II, nx, per_x)] < 0 || u[II]*u[δy⁺(II, ny, per_y)] < 0 || u[II]*u[δx⁻(II, nx, per_x)] < 0 || u[II]*u[δy⁻(II, ny, per_y)] < 0, true, false)
@inline is_near_interface_l(u, II::CartesianIndex, nx, ny, per_x, per_y) = @inbounds ifelse(u[II]*u[δx⁺(II, nx, per_x)] < 0 || u[II]*u[δy⁺(II, ny, per_y)] < 0 || u[II]*u[δy⁻(II, ny, per_y)] < 0, true, false)
@inline is_near_interface_b(u, II::CartesianIndex, nx, ny, per_x, per_y) = @inbounds ifelse(u[II]*u[δx⁺(II, nx, per_x)] < 0 || u[II]*u[δy⁺(II, ny, per_y)] < 0 || u[II]*u[δx⁻(II, nx, per_x)] < 0, true, false)
@inline is_near_interface_r(u, II::CartesianIndex, nx, ny, per_x, per_y) = @inbounds ifelse(u[II]*u[δy⁺(II, ny, per_y)] < 0 || u[II]*u[δx⁻(II, nx, per_x)] < 0 || u[II]*u[δy⁻(II, ny, per_y)] < 0, true, false)
@inline is_near_interface_t(u, II::CartesianIndex, nx, ny, per_x, per_y) = @inbounds ifelse(u[II]*u[δx⁺(II, nx, per_x)] < 0 || u[II]*u[δx⁻(II, nx, per_x)] < 0 || u[II]*u[δy⁻(II, ny, per_y)] < 0, true, false)
@inline is_near_interface_bl(u, II::CartesianIndex, nx, ny, per_x, per_y) = @inbounds ifelse(u[II]*u[δx⁺(II, nx, per_x)] < 0 || u[II]*u[δy⁺(II, ny, per_y)] < 0, true, false)
@inline is_near_interface_br(u, II::CartesianIndex, nx, ny, per_x, per_y) = @inbounds ifelse(u[II]*u[δy⁺(II, ny, per_y)] < 0 || u[II]*u[δx⁻(II, nx, per_x)] < 0, true, false)
@inline is_near_interface_tr(u, II::CartesianIndex, nx, ny, per_x, per_y) = @inbounds ifelse(u[II]*u[δx⁻(II, nx, per_x)] < 0 || u[II]*u[δy⁻(II, ny, per_y)] < 0, true, false)
@inline is_near_interface_tl(u, II::CartesianIndex, nx, ny, per_x, per_y) = @inbounds ifelse(u[II]*u[δx⁺(II, nx, per_x)] < 0 || u[II]*u[δy⁻(II, ny, per_y)] < 0, true, false)

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
   @inbounds for i in eachindex(nodes)
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

function find_radius(grid, LS)
    @unpack x, y = grid
    @unpack mid_point, MIXED = LS

    radius = Vector{Float64}(undef,0)
    for II in MIXED
        try
            radius_ = distance(Point(0.,0.), Point(x[II]+mid_point[II].x, y[II]+mid_point[II].y))
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

function clip_large_cell!(grid, LS, geo, II, ϵ, neighbours)
    @unpack nx, ny = grid

    if geo.cap[II,5] > (1.0-ϵ)
        geo.cap[II,:] .= vcat(ones(7), 0.5.*ones(4))
        if neighbours
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
        end
        geo.centroid[II] = Point(0.0, 0.0)
        LS.mid_point[II] = Point(0.0, 0.0)
    end
    return nothing
end

function clip_small_cell!(grid, LS, geo, II, ϵ, neighbours)
    @unpack nx, ny = grid

    if geo.cap[II,5] < ϵ
        geo.cap[II,:] .= 0.
        geo.emptied[II] = true
        if neighbours
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
        end
        geo.centroid[II] = Point(0.0, 0.0)
        LS.mid_point[II] = Point(0.0, 0.0)
    end
    return nothing
end

function empty_cell!(grid, LS, geo, II, neighbours = true)
    @unpack nx, ny = grid
    
    geo.cap[II,:] .= 0.0
    geo.emptied[II] = true
    if neighbours
        if II[2] > 1
            geo.cap[δx⁻(II),3] = 0.0
        end
        if II[1] > 1
            geo.cap[δy⁻(II),4] = 0.0
        end
        if II[2] < nx
            geo.cap[δx⁺(II),1] = 0.0
        end
        if II[1] < ny
            geo.cap[δy⁺(II),2] = 0.0
        end
    end
    geo.centroid[II] = Point(0.0, 0.0)
    LS.mid_point[II] = Point(0.0, 0.0)
    LS.iso[II] = 15.0
end

function clip_cells!(grid, LS, ϵ, ϵwall, neighbours, BC_int)
    @unpack nx, ny, ind = grid
    @unpack geoS, geoL = LS

    @inbounds @threads for II in ind.inside
        clip_large_cell!(grid, LS, geoS, II, ϵ, neighbours)
        clip_large_cell!(grid, LS, geoL, II, ϵ, neighbours)
    end
    for iLS in eachindex(BC_int)
        if is_wall(BC_int[iLS])
            @inbounds @threads for II in vcat(grid.LS[iLS].MIXED, grid.LS[iLS].SOLID)
                clip_large_cell!(grid, LS, geoS, II, ϵwall, neighbours)
                clip_large_cell!(grid, LS, geoL, II, ϵwall, neighbours)
            end
        end
    end
    @inbounds @threads for II in ind.inside
        clip_small_cell!(grid, LS, geoS, II, ϵ, neighbours)
        clip_small_cell!(grid, LS, geoL, II, ϵ, neighbours)
    end
    for iLS in eachindex(BC_int)
        if is_wall(BC_int[iLS])
            @inbounds @threads for II in vcat(grid.LS[iLS].MIXED, grid.LS[iLS].SOLID)
                clip_small_cell!(grid, LS, geoS, II, ϵwall, neighbours)
                clip_small_cell!(grid, LS, geoL, II, ϵwall, neighbours)
            end
        end
    end
    @inbounds @threads for II in @view ind.b_left[1][2:end-1]
        if geoS.cap[II,5] > (1.0-ϵ)
            geoS.cap[II,1] = 1.0
            geoS.cap[II,2:11] .= vcat(ones(6), 0.5.*ones(4))
            if neighbours
                geoS.cap[δy⁻(II),4] = 1.0
                geoS.cap[δx⁺(II),1] = 1.0
                geoS.cap[δy⁺(II),2] = 1.0
            end
        end
    end
    @inbounds @threads for II in @view ind.b_left[1][2:end-1]
        if geoS.cap[II,5] < ϵ
            geoS.cap[II,:] .= 0.
            geoS.emptied[II] = true
            if neighbours
                geoS.cap[δy⁻(II),4] = 0.0
                geoS.cap[δx⁺(II),1] = 0.0
                geoS.cap[δy⁺(II),2] = 0.0
            end
        end
    end
    @inbounds @threads for II in @view ind.b_bottom[1][2:end-1]
        if geoS.cap[II,5] > (1.0-ϵ)
            geoS.cap[II,1] = 1.0
            geoS.cap[II,2] = 1.0
            geoS.cap[II,3:11] .= vcat(ones(5), 0.5.*ones(4))
            if neighbours
                geoS.cap[δx⁻(II),3] = 1.0
                geoS.cap[δx⁺(II),1] = 1.0
                geoS.cap[δy⁺(II),2] = 1.0
            end
        end
    end
    @inbounds @threads for II in @view ind.b_bottom[1][2:end-1]
        if geoS.cap[II,5] < ϵ
            geoS.cap[II,:] .= 0.
            geoS.emptied[II] = true
            if neighbours
                geoS.cap[δx⁻(II),3] = 0.0
                geoS.cap[δx⁺(II),1] = 0.0
                geoS.cap[δy⁺(II),2] = 0.0
            end
        end
    end
    @inbounds @threads for II in @view ind.b_right[1][2:end-1]
        if geoS.cap[II,5] > (1.0-ϵ)
            geoS.cap[II,1:2] .= ones(2)
            geoS.cap[II,3] = 1.0
            geoS.cap[II,4:11] .= vcat(ones(4), 0.5.*ones(4))
            if neighbours
                geoS.cap[δx⁻(II),3] = 1.0
                geoS.cap[δy⁻(II),4] = 1.0
                geoS.cap[δy⁺(II),2] = 1.0
            end
        end
    end
    @inbounds @threads for II in @view ind.b_right[1][2:end-1]
        if geoS.cap[II,5] < ϵ
            geoS.cap[II,:] .= 0.
            geoS.emptied[II] = true
            if neighbours
                geoS.cap[δx⁻(II),3] = 0.0
                geoS.cap[δy⁻(II),4] = 0.0
                geoS.cap[δy⁺(II),2] = 0.0
            end
        end
    end
    @inbounds @threads for II in @view ind.b_top[1][2:end-1]
        if geoS.cap[II,5] > (1.0-ϵ)
            geoS.cap[II,1:3] .= ones(3)
            geoS.cap[II,4] = 1.0
            geoS.cap[II,5:11] .= vcat(ones(3), 0.5.*ones(4))
            if neighbours
                geoS.cap[δx⁻(II),3] = 1.0
                geoS.cap[δy⁻(II),4] = 1.0
                geoS.cap[δx⁺(II),1] = 1.0
            end
        end
    end
    @inbounds @threads for II in @view ind.b_top[1][2:end-1]
        if geoS.cap[II,5] < ϵ
            geoS.cap[II,:] .= 0.
            geoS.emptied[II] = true
            if neighbours
                geoS.cap[δx⁻(II),3] = 0.0
                geoS.cap[δy⁻(II),4] = 0.0
                geoS.cap[δx⁺(II),1] = 0.0
            end
        end
    end
    @inbounds @threads for II in @view ind.b_left[1][2:end-1]
        if geoL.cap[II,5] > (1.0-ϵ)
            geoL.cap[II,1] = 1.0
            geoL.cap[II,2:11] .= vcat(ones(6), 0.5.*ones(4))
            if neighbours
                geoL.cap[δy⁻(II),4] = 1.0
                geoL.cap[δx⁺(II),1] = 1.0
                geoL.cap[δy⁺(II),2] = 1.0
            end
        end
    end
    @inbounds @threads for II in @view ind.b_left[1][2:end-1]
        if geoL.cap[II,5] < ϵ
            geoL.cap[II,:] .= 0.
            geoL.emptied[II] = true
            if neighbours
                geoL.cap[δy⁻(II),4] = 0.0
                geoL.cap[δx⁺(II),1] = 0.0
                geoL.cap[δy⁺(II),2] = 0.0
            end
        end
    end
    @inbounds @threads for II in @view ind.b_bottom[1][2:end-1]
        if geoL.cap[II,5] > (1.0-ϵ)
            geoL.cap[II,1] = 1.0
            geoL.cap[II,2] = 1.0
            geoL.cap[II,3:11] .= vcat(ones(5), 0.5.*ones(4))
            if neighbours
                geoL.cap[δx⁻(II),3] = 1.0
                geoL.cap[δx⁺(II),1] = 1.0
                geoL.cap[δy⁺(II),2] = 1.0
            end
        end
    end
    @inbounds @threads for II in @view ind.b_bottom[1][2:end-1]
        if geoL.cap[II,5] < ϵ
            geoL.cap[II,:] .= 0.
            geoL.emptied[II] = true
            if neighbours
                geoL.cap[δx⁻(II),3] = 0.0
                geoL.cap[δx⁺(II),1] = 0.0
                geoL.cap[δy⁺(II),2] = 0.0
            end
        end
    end
    @inbounds @threads for II in @view ind.b_right[1][2:end-1]
        if geoL.cap[II,5] > (1.0-ϵ)
            geoL.cap[II,1:2] .= ones(2)
            geoL.cap[II,3] = 1.0
            geoL.cap[II,4:11] .= vcat(ones(4), 0.5.*ones(4))
            if neighbours
                geoL.cap[δx⁻(II),3] = 1.0
                geoL.cap[δy⁻(II),4] = 1.0
                geoL.cap[δy⁺(II),2] = 1.0
            end
        end
    end
    @inbounds @threads for II in @view ind.b_right[1][2:end-1]
        if geoL.cap[II,5] < ϵ
            geoL.cap[II,:] .= 0.
            geoL.emptied[II] = true
            if neighbours
                geoL.cap[δx⁻(II),3] = 0.0
                geoL.cap[δy⁻(II),4] = 0.0
                geoL.cap[δy⁺(II),2] = 0.0
            end
        end
    end
    @inbounds @threads for II in @view ind.b_top[1][2:end-1]
        if geoL.cap[II,5] > (1.0-ϵ)
            geoL.cap[II,1:3] .= ones(3)
            geoL.cap[II,4] = 1.0
            geoL.cap[II,5:11] .= vcat(ones(3), 0.5.*ones(4))
            if neighbours
                geoL.cap[δx⁻(II),3] = 1.0
                geoL.cap[δy⁻(II),4] = 1.0
                geoL.cap[δx⁺(II),1] = 1.0
            end
        end
    end
    @inbounds @threads for II in @view ind.b_top[1][2:end-1]
        if geoL.cap[II,5] < ϵ
            geoL.cap[II,:] .= 0.
            geoL.emptied[II] = true
            if neighbours
                geoL.cap[δx⁻(II),3] = 0.0
                geoL.cap[δy⁻(II),4] = 0.0
                geoL.cap[δx⁺(II),1] = 0.0
            end
        end
    end
    return nothing
end

function clip_A_acc_to_V(grid, grid_u, grid_v, geo, geo_u, geo_v, ϵ, ϵwall, neighbours, BC_int)
    @unpack ind = grid
    @inbounds @threads for II in ind.all_indices
        if geo.cap[II,5] < ϵ
            geo_u.cap[II,3] = 0.0
            geo_v.cap[II,4] = 0.0
            if neighbours
                geo_u.cap[δx⁺(II),1] = 0.0
                geo_v.cap[δy⁺(II),2] = 0.0
            end
        end
    end
    for iLS in eachindex(BC_int)
        if is_wall(BC_int[iLS])
            @inbounds @threads for II in vcat(grid.LS[iLS].MIXED, grid.LS[iLS].SOLID)
                if geo.cap[II,5] < ϵwall
                    geo_u.cap[II,3] = 0.0
                    geo_v.cap[II,4] = 0.0
                    if neighbours
                        geo_u.cap[δx⁺(II),1] = 0.0
                        geo_v.cap[δy⁺(II),2] = 0.0
                    end
                end
            end
        end
    end
    @inbounds @threads for II in grid_u.ind.all_indices
        if geo_u.cap[II,5] < ϵ
            if neighbours
                if II[2] > 1
                    geo.cap[δx⁻(II),3] = 0.
                end
            end
            if II[2] < grid_u.nx
                geo.cap[II,1] = 0.
            end
        end
    end
    for iLS in eachindex(BC_int)
        if is_wall(BC_int[iLS])
            @inbounds @threads for II in vcat(grid_u.LS[iLS].MIXED, grid_u.LS[iLS].SOLID)
                if geo_u.cap[II,5] < ϵwall
                    if neighbours
                        if II[2] > 1
                            geo.cap[δx⁻(II),3] = 0.
                        end
                    end
                    if II[2] < grid_u.nx
                        geo.cap[II,1] = 0.
                    end
                end
            end
        end
    end
    @inbounds @threads for II in grid_v.ind.all_indices
        if geo_v.cap[II,5] < ϵ
            if neighbours
                if II[1] > 1
                    geo.cap[δy⁻(II),4] = 0.
                end
            end
            if II[1] < grid_v.ny
                geo.cap[II,2] = 0.
            end
        end
    end
    for iLS in eachindex(BC_int)
        if is_wall(BC_int[iLS])
            @inbounds @threads for II in vcat(grid_v.LS[iLS].MIXED, grid_v.LS[iLS].SOLID)
                if geo_v.cap[II,5] < ϵwall
                    if neighbours
                        if II[1] > 1
                            geo.cap[δy⁻(II),4] = 0.
                        end
                    end
                    if II[1] < grid_v.ny
                        geo.cap[II,2] = 0.
                    end
                end
            end
        end
    end

    return nothing
end

"""
    clip_middle_cells!(grid, LS)

Clip the middle cells if the cells on top and bottom or
at the left and right are empty simultaneously.
"""
function clip_middle_cells!(grid, LS)

    @inbounds @threads for II in grid.ind.inside
        if LS.geoL.emptied[δx⁻(II)] && LS.geoL.emptied[δx⁺(II)]
            empty_cell!(grid, LS, LS.geoL, II)
        end
        if LS.geoL.emptied[δy⁻(II)] && LS.geoL.emptied[δy⁺(II)]
            empty_cell!(grid, LS, LS.geoL, II)
        end

        if LS.geoS.emptied[δx⁻(II)] && LS.geoS.emptied[δx⁺(II)]
            empty_cell!(grid, LS, LS.geoS, II)
        end
        if LS.geoS.emptied[δy⁻(II)] && LS.geoS.emptied[δy⁺(II)]
            empty_cell!(grid, LS, LS.geoS, II)
        end
    end

    return nothing
end

"""
    dimensionalize!()

    dcap... = cap*dx...dx*dy...
"""
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
"""
    postprocess_grids1!()

    clip_cells!, clip_A_acc_to_V, set_cap_bcs!
"""
function postprocess_grids1!(num, grid, LS, grid_u, LS_u, grid_v, LS_v, periodic_x, periodic_y, neighbours, empty, BC_int)
    clip_cells!(grid, LS, num.ϵ, num.ϵwall, neighbours, BC_int)
    clip_cells!(grid_u, LS_u, num.ϵ, num.ϵwall, neighbours, BC_int)
    clip_cells!(grid_v, LS_v, num.ϵ, num.ϵwall, neighbours, BC_int)

    clip_A_acc_to_V(grid, grid_u, grid_v, LS.geoS, LS_u.geoS, LS_v.geoS, num.ϵ, num.ϵwall, neighbours, BC_int)
    clip_A_acc_to_V(grid, grid_u, grid_v, LS.geoL, LS_u.geoL, LS_v.geoL, num.ϵ, num.ϵwall, neighbours, BC_int)

    clip_middle_cells!(grid, LS)
    clip_middle_cells!(grid_u, LS_u)
    clip_middle_cells!(grid_v, LS_v)

    # Compute angles this way so that they exist next to clipped cells
    non_mixed = vcat(LS.LIQUID, LS.SOLID)
    LS.α[non_mixed] .= atan.(
        LS.geoL.cap[non_mixed,4] .- LS.geoL.cap[non_mixed,2],
        LS.geoL.cap[non_mixed,3] .- LS.geoL.cap[non_mixed,1]
    )
    non_mixed_u = vcat(LS_u.LIQUID, LS_u.SOLID)
    LS_u.α[non_mixed_u] .= atan.(
        LS_u.geoL.cap[non_mixed_u,4] .- LS_u.geoL.cap[non_mixed_u,2],
        LS_u.geoL.cap[non_mixed_u,3] .- LS_u.geoL.cap[non_mixed_u,1]
    )
    non_mixed_v = vcat(LS_v.LIQUID, LS_v.SOLID)
    LS_v.α[non_mixed_v] .= atan.(
        LS_v.geoL.cap[non_mixed_v,4] .- LS_v.geoL.cap[non_mixed_v,2],
        LS_v.geoL.cap[non_mixed_v,3] .- LS_v.geoL.cap[non_mixed_v,1]
    )
    # try
    #     set_cap_bcs!(grid, num, LS, periodic_x, periodic_y, empty)
    #     set_cap_bcs!(grid_u, num, LS_u, periodic_x, periodic_y, empty)
    #     set_cap_bcs!(grid_v, num, LS_v, periodic_x, periodic_y, empty)

    # catch errorLS
    #     printstyled(color=:red, @sprintf "\n LS not updated: set_cap_bcs! p u v \n")
    #     print(errorLS)
    # end
    set_cap_bcs!(grid, num, LS, periodic_x, periodic_y, empty)
    set_cap_bcs!(grid_u, num, LS_u, periodic_x, periodic_y, empty)
    set_cap_bcs!(grid_v, num, LS_v, periodic_x, periodic_y, empty)

    return nothing
end

"""
    postprocess_grids2!()

    dimensionalize!, Wcapacities!, average_face_capacities!
"""
function postprocess_grids2!(grid::Mesh{Flower.GridCC,Float64,Int64},
    LS::Levelset{Float64, Int64},
    grid_u::Mesh{Flower.GridFCx, Float64, Int64},
    LS_u::Levelset{Float64, Int64},
    grid_v::Mesh{Flower.GridFCy, Float64, Int64},
    LS_v::Levelset{Float64, Int64}, 
    periodic_x::Bool, 
    periodic_y::Bool, 
    average::Bool)
    dimensionalize!(grid, LS.geoS)
    dimensionalize!(grid, LS.geoL)
    dimensionalize!(grid_u, LS_u.geoS)
    dimensionalize!(grid_u, LS_u.geoL)
    dimensionalize!(grid_v, LS_v.geoS)
    dimensionalize!(grid_v, LS_v.geoL)

    Wcapacities!(LS.geoS.dcap, periodic_x, periodic_y)
    Wcapacities!(LS.geoL.dcap, periodic_x, periodic_y)
    Wcapacities!(LS_u.geoS.dcap, periodic_x, periodic_y)
    Wcapacities!(LS_u.geoL.dcap, periodic_x, periodic_y)
    Wcapacities!(LS_v.geoS.dcap, periodic_x, periodic_y)
    Wcapacities!(LS_v.geoL.dcap, periodic_x, periodic_y)

    if average
        average_face_capacities!(LS.geoS.dcap, periodic_x, periodic_y)
        average_face_capacities!(LS.geoL.dcap, periodic_x, periodic_y)
        average_face_capacities!(LS_u.geoS.dcap, periodic_x, periodic_y)
        average_face_capacities!(LS_u.geoL.dcap, periodic_x, periodic_y)
        average_face_capacities!(LS_v.geoS.dcap, periodic_x, periodic_y)
        average_face_capacities!(LS_v.geoL.dcap, periodic_x, periodic_y)
    end

    return nothing
end

"""
    ilp2cap()

Compute the length of a line intersecting a polygon.
"""
function ilp2cap(l, p)
    inters = LibGEOS.intersection(l, p)
    obj = GeoInterface.convert(GeometryBasics, inters)
    if typeof(obj) == Point2{Float64}
        len = 0.0
    else
        len = LibGEOS.distance(obj.points[1][1], obj.points[1][2])
    end

    return len
end

"""
    crossing_2levelsets(grid, LS1, LS2)

Update geometric moments when two interfaces cross.
"""
function crossing_2levelsets!(num, grid, LS1, LS2, BC_int)
    @unpack x, y, nx, ny, dx, dy, ind = grid

    f_left = readgeom("LINESTRING(0.0 0.0, 0.0 1.0)")
    f_bottom = readgeom("LINESTRING(0.0 0.0, 1.0 0.0)")
    f_right = readgeom("LINESTRING(1.0 0.0, 1.0 1.0)")
    f_top = readgeom("LINESTRING(0.0 1.0, 1.0 1.0)")
    cell_faces = [f_left, f_bottom, f_right, f_top]

    # Cells that have been clipped and became full should be removed from mixed cells
    ls1_liquid = copy(LS1.LIQUID)
    ls1_mixed = copy(LS1.MIXED)
    idx_filled = findall(LS1.geoL.cap[LS1.MIXED,5] .> 1 .- num.ϵ)
    append!(ls1_liquid, ls1_mixed[idx_filled])
    deleteat!(ls1_mixed, idx_filled)
    idx_emptied = findall(LS1.geoL.emptied[ls1_mixed])
    deleteat!(ls1_mixed, idx_emptied)

    ls2_liquid = copy(LS2.LIQUID)
    ls2_mixed = copy(LS2.MIXED)
    idx_filled = findall(LS2.geoL.cap[LS2.MIXED,5] .> 1 .- num.ϵ)
    idx_emptied = findall(LS2.geoL.cap[LS2.MIXED,5] .> 1 .- num.ϵ)
    append!(ls2_liquid, ls2_mixed[idx_filled])
    deleteat!(ls2_mixed, idx_filled)
    idx_emptied = findall(LS2.geoL.emptied[ls2_mixed])
    deleteat!(ls2_mixed, idx_emptied)

    full_mixed = intersect(ls1_liquid, ls2_mixed)
    mixed_full = intersect(ls2_liquid, ls1_mixed)
    @inbounds @threads for II in full_mixed
        grid.LS[end].geoL.cap[II,:] .= LS2.geoL.cap[II,:]
        grid.LS[end].geoL.centroid[II] = LS2.geoL.centroid[II]
    end
    @inbounds @threads for II in mixed_full
        grid.LS[end].geoL.cap[II,:] .= LS1.geoL.cap[II,:]
        grid.LS[end].geoL.centroid[II] = LS1.geoL.centroid[II]
    end

    mixed_mixed = intersect(ls1_mixed, ls2_mixed)
    Acap = [3,4,1,2]
    @inbounds for II in mixed_mixed
        if is_fs(BC_int[1])
            append!(LS1.cl, [II])
        elseif is_fs(BC_int[2])
            append!(LS2.cl, [II])
        end
        if (abs(LS1.cut_points[II][1].x) == 0.5 && abs(LS1.cut_points[II][1].y) == 0.5)
            LS1.cut_points[II][1] -= Point(1e-12 * sign(LS1.cut_points[II][1].x), 1e-12 * sign(LS1.cut_points[II][1].y))
        end
        if (abs(LS1.cut_points[II][2].x) == 0.5 && abs(LS1.cut_points[II][2].y) == 0.5)
            LS1.cut_points[II][2] -= Point(1e-12 * sign(LS1.cut_points[II][2].x), 1e-12 * sign(LS1.cut_points[II][2].y))
        end
        if (abs(LS2.cut_points[II][1].x) == 0.5 && abs(LS2.cut_points[II][1].y) == 0.5)
            LS2.cut_points[II][1] -= Point(1e-12 * sign(LS2.cut_points[II][1].x), 1e-12 * sign(LS2.cut_points[II][1].y))
        end
        if (abs(LS2.cut_points[II][2].x) == 0.5 && abs(LS2.cut_points[II][2].y) == 0.5)
            LS2.cut_points[II][2] -= Point(1e-12 * sign(LS2.cut_points[II][2].x), 1e-12 * sign(LS2.cut_points[II][2].y))
        end
        l1 = Line(LS1.cut_points[II][1], LS1.cut_points[II][2])
        l2 = Line(LS2.cut_points[II][1], LS2.cut_points[II][2])
        p = line_intersection(l1, l2)

        poly1 = points2polygon(LS1.geoL.vertices[II])
        vol1 = LibGEOS.area(poly1)
        poly2 = points2polygon(LS2.geoL.vertices[II])
        vol2 = LibGEOS.area(poly2)
        poly = LibGEOS.intersection(poly1, poly2)
        vol = LibGEOS.area(poly)

        # If the intersection is not empty
        if vol > num.ϵwall
            # Add point as contact line point
            _cent = LibGEOS.centroid(poly)
            cent = GeoInterface.convert(GeometryBasics, _cent)

            lx = readgeom(
                "LINESTRING(
                    $(cent.data[1]) $(cent.data[2]-2.0),
                    $(cent.data[1]) $(cent.data[2]+2.0)
                )"
            )
            ly = readgeom(
                "LINESTRING(
                    $(cent.data[1]-2.0) $(cent.data[2]),
                    $(cent.data[1]+2.0) $(cent.data[2])
                )"
            )
            Bx = ilp2cap(lx, poly)
            By = ilp2cap(ly, poly)

            grid.LS[end].geoL.cap[II,1] = minimum([LS1.geoL.cap[II,1], LS2.geoL.cap[II,1]])
            grid.LS[end].geoL.cap[II,2] = minimum([LS1.geoL.cap[II,2], LS2.geoL.cap[II,2]])
            grid.LS[end].geoL.cap[II,3] = minimum([LS1.geoL.cap[II,3], LS2.geoL.cap[II,3]])
            grid.LS[end].geoL.cap[II,4] = minimum([LS1.geoL.cap[II,4], LS2.geoL.cap[II,4]])
            grid.LS[end].geoL.cap[II,5] = vol
            grid.LS[end].geoL.cap[II,6] = Bx
            grid.LS[end].geoL.cap[II,7] = By
            grid.LS[end].geoL.cap[II,8:11] .= 0.5 * vol
            grid.LS[end].geoL.centroid[II] = Point(cent.data[1], cent.data[2])

            if within_cell(p)
                min_face = zeros(Int, 4) # Interface with minimum moment at each face
                double_mixed = zeros(Bool, 4)
                # The smallest keeps its value and the largest the intersection value
                for capn in 1:4
                    # Check special case where both interfaces intersect the cell
                    # but the moments don't overlap
                    face1 = LibGEOS.intersection(poly1, cell_faces[capn])
                    face2 = LibGEOS.intersection(poly2, cell_faces[capn])
                    if ((LS1.geoL.cap[II,capn] > 1e-12) && (LS2.geoL.cap[II,capn] > 1e-12) && 
                        (LS1.geoL.cap[II,capn] < (1.0 - 1e-12)) && (LS2.geoL.cap[II,capn] < (1.0 - 1e-12)) &&
                        (!GeoInterface.contains(face1, face2) && !GeoInterface.contains(face2, face1)))
                        double_mixed[capn] = true
                    end

                    set_A_caps!(capn, LS1, LS2, II, poly1, poly2, p, min_face)
                end
    
                # If the face is not intersected by any interface
                for capn in 1:4
                    if min_face[capn] == 0 && min_face[Acap[capn]] != 0 && !double_mixed[Acap[capn]]
                        set_A_caps_full_mixed!(capn, LS1, LS2, II, poly1, poly2, p, min_face[Acap[capn]])
                    elseif min_face[capn] == 0 && double_mixed[Acap[capn]]
                        set_A_caps_double_mixed_full_empty!(capn, Acap[capn], LS1, LS2, II, poly1, p)
                        if capn == 2 || capn == 4
                            lx = readgeom(
                                "LINESTRING(
                                    $(p.x+0.5) $(p.y-2.0),
                                    $(p.x+0.5) $(p.y+2.0)
                                )"
                            )
                            try                                
                                Bx = ilp2cap(lx, poly)
                            catch
                                face_non_empty = Acap[capn]
                                if face_non_empty == 2
                                    id_cut1 = abs(LS1.cut_points[II][1].y + 0.5) < 1e-8 ? 1 : 2
                                    id_cut2 = abs(LS2.cut_points[II][1].y + 0.5) < 1e-8 ? 1 : 2
                                    if (p.x + 0.5) > LS1.cut_points[II][id_cut1].x
                                        if LS1.cut_points[II][id_cut1].x > LS2.cut_points[II][id_cut2].x
                                            Bx = ilp2cap(lx, poly2)
                                        else
                                            Bx = ilp2cap(lx, poly1)
                                        end
                                    else
                                        if LS1.cut_points[II][id_cut1].x < LS2.cut_points[II][id_cut2].x
                                            Bx = ilp2cap(lx, poly2)
                                        else
                                            Bx = ilp2cap(lx, poly1)
                                        end
                                    end
                                else
                                    id_cut1 = abs(LS1.cut_points[II][1].y - 0.5) < 1e-8 ? 1 : 2
                                    id_cut2 = abs(LS2.cut_points[II][1].y - 0.5) < 1e-8 ? 1 : 2
                                    if (p.x + 0.5) > LS1.cut_points[II][id_cut1].x
                                        if LS1.cut_points[II][id_cut1].x > LS2.cut_points[II][id_cut2].x
                                            Bx = ilp2cap(lx, poly2)
                                        else
                                            Bx = ilp2cap(lx, poly1)
                                        end
                                    else
                                        if LS1.cut_points[II][id_cut1].x < LS2.cut_points[II][id_cut2].x
                                            Bx = ilp2cap(lx, poly2)
                                        else
                                            Bx = ilp2cap(lx, poly1)
                                        end
                                    end
                                end
                            end
                            grid.LS[end].geoL.cap[II,6] = Bx
                        else
                            ly = readgeom(
                                "LINESTRING(
                                    $(p.x-2.0) $(p.y+0.5),
                                    $(p.x+2.0) $(p.y+0.5)
                                )"
                            )
                            try
                                By = ilp2cap(ly, poly)
                            catch
                                face_non_empty = Acap[capn]
                                if face_non_empty == 1
                                    id_cut1 = abs(LS1.cut_points[II][1].x + 0.5) < 1e-8 ? 1 : 2
                                    id_cut2 = abs(LS2.cut_points[II][1].x + 0.5) < 1e-8 ? 1 : 2
                                    if (p.y + 0.5) > LS1.cut_points[II][id_cut1].y
                                        if LS1.cut_points[II][id_cut1].y > LS2.cut_points[II][id_cut2].y
                                            By = ilp2cap(ly, poly2)
                                        else
                                            By = ilp2cap(ly, poly1)
                                        end
                                    else
                                        if LS1.cut_points[II][id_cut1].y < LS2.cut_points[II][id_cut2].y
                                            By = ilp2cap(ly, poly2)
                                        else
                                            By = ilp2cap(ly, poly1)
                                        end
                                    end
                                else
                                    id_cut1 = abs(LS1.cut_points[II][1].x - 0.5) < 1e-8 ? 1 : 2
                                    id_cut2 = abs(LS2.cut_points[II][1].x - 0.5) < 1e-8 ? 1 : 2
                                    if (p.y + 0.5) > LS1.cut_points[II][id_cut1].y
                                        if LS1.cut_points[II][id_cut1].y > LS2.cut_points[II][id_cut2].y
                                            By = ilp2cap(ly, poly2)
                                        else
                                            By = ilp2cap(ly, poly1)
                                        end
                                    else
                                        if LS1.cut_points[II][id_cut1].y < LS2.cut_points[II][id_cut2].y
                                            By = ilp2cap(ly, poly2)
                                        else
                                            By = ilp2cap(ly, poly1)
                                        end
                                    end
                                end
                            end
                            grid.LS[end].geoL.cap[II,7] = By
                        end
                    elseif (min_face[capn] == 0 && 
                        (LS1.geoL.cap[II,capn]*LS2.geoL.cap[II,capn] + LS1.geoL.cap[II,Acap[capn]]*LS2.geoL.cap[II,Acap[capn]]) < 1e-12
                        )
                        if capn == 2 || capn == 4
                            ly = readgeom(
                                "LINESTRING(
                                    $(p.x-2.0) $(p.y+0.5),
                                    $(p.x+2.0) $(p.y+0.5)
                                )"
                            )
                            try
                                By = ilp2cap(ly, poly)
                            catch
                                face_non_empty = Acap[capn]
                                if face_non_empty == 1
                                    id_cut1 = abs(LS1.cut_points[II][1].x + 0.5) < 1e-8 ? 1 : 2
                                    id_cut2 = abs(LS2.cut_points[II][1].x + 0.5) < 1e-8 ? 1 : 2
                                    if (p.y + 0.5) > LS1.cut_points[II][id_cut1].y
                                        if LS1.cut_points[II][id_cut1].y > LS2.cut_points[II][id_cut2].y
                                            By = ilp2cap(ly, poly2)
                                        else
                                            By = ilp2cap(ly, poly1)
                                        end
                                    else
                                        if LS1.cut_points[II][id_cut1].y < LS2.cut_points[II][id_cut2].y
                                            By = ilp2cap(ly, poly2)
                                        else
                                            By = ilp2cap(ly, poly1)
                                        end
                                    end
                                else
                                    id_cut1 = abs(LS1.cut_points[II][1].x - 0.5) < 1e-8 ? 1 : 2
                                    id_cut2 = abs(LS2.cut_points[II][1].x - 0.5) < 1e-8 ? 1 : 2
                                    if (p.y + 0.5) > LS1.cut_points[II][id_cut1].y
                                        if LS1.cut_points[II][id_cut1].y > LS2.cut_points[II][id_cut2].y
                                            By = ilp2cap(ly, poly2)
                                        else
                                            By = ilp2cap(ly, poly1)
                                        end
                                    else
                                        if LS1.cut_points[II][id_cut1].y < LS2.cut_points[II][id_cut2].y
                                            By = ilp2cap(ly, poly2)
                                        else
                                            By = ilp2cap(ly, poly1)
                                        end
                                    end
                                end
                            end
                            grid.LS[end].geoL.cap[II,7] = By
                            if LS1.geoL.cap[II,capn] > 1e-12
                                LS1.geoL.cap[II,capn] = By
                                LS1.geoL.cap[II,Acap[capn]] = 0.0
                                LS2.geoL.cap[II,Acap[capn]] = By
                                LS2.geoL.cap[II,capn] = 0.0
                            else
                                LS1.geoL.cap[II,capn] = 0.0
                                LS1.geoL.cap[II,Acap[capn]] = By
                                LS2.geoL.cap[II,Acap[capn]] = 0.0
                                LS2.geoL.cap[II,capn] = By
                            end
                        else
                            lx = readgeom(
                                "LINESTRING(
                                    $(p.x+0.5) $(p.y-2.0),
                                    $(p.x+0.5) $(p.y+2.0)
                                )"
                            )
                            try                                
                                Bx = ilp2cap(lx, poly)
                            catch
                                face_non_empty = Acap[capn]
                                if face_non_empty == 2
                                    id_cut1 = abs(LS1.cut_points[II][1].y + 0.5) < 1e-8 ? 1 : 2
                                    id_cut2 = abs(LS2.cut_points[II][1].y + 0.5) < 1e-8 ? 1 : 2
                                    if (p.x + 0.5) > LS1.cut_points[II][id_cut1].x
                                        if LS1.cut_points[II][id_cut1].x > LS2.cut_points[II][id_cut2].x
                                            Bx = ilp2cap(lx, poly2)
                                        else
                                            Bx = ilp2cap(lx, poly1)
                                        end
                                    else
                                        if LS1.cut_points[II][id_cut1].x < LS2.cut_points[II][id_cut2].x
                                            Bx = ilp2cap(lx, poly2)
                                        else
                                            Bx = ilp2cap(lx, poly1)
                                        end
                                    end
                                else
                                    id_cut1 = abs(LS1.cut_points[II][1].y - 0.5) < 1e-8 ? 1 : 2
                                    id_cut2 = abs(LS2.cut_points[II][1].y - 0.5) < 1e-8 ? 1 : 2
                                    if (p.x + 0.5) > LS1.cut_points[II][id_cut1].x
                                        if LS1.cut_points[II][id_cut1].x > LS2.cut_points[II][id_cut2].x
                                            Bx = ilp2cap(lx, poly2)
                                        else
                                            Bx = ilp2cap(lx, poly1)
                                        end
                                    else
                                        if LS1.cut_points[II][id_cut1].x < LS2.cut_points[II][id_cut2].x
                                            Bx = ilp2cap(lx, poly2)
                                        else
                                            Bx = ilp2cap(lx, poly1)
                                        end
                                    end
                                end
                            end
                            grid.LS[end].geoL.cap[II,6] = Bx
                            if LS1.geoL.cap[II,capn] > 1e-12
                                LS1.geoL.cap[II,capn] = Bx
                                LS1.geoL.cap[II,Acap[capn]] = 0.0
                                LS2.geoL.cap[II,Acap[capn]] = Bx
                                LS2.geoL.cap[II,capn] = 0.0
                            else
                                LS1.geoL.cap[II,capn] = 0.0
                                LS1.geoL.cap[II,Acap[capn]] = Bx
                                LS2.geoL.cap[II,Acap[capn]] = 0.0
                                LS2.geoL.cap[II,capn] = Bx
                            end
                        end
                    elseif min_face[capn] == 0
                        set_A_caps_full_empty!(capn, Acap[capn], LS1, LS2, II, poly1, p)
                        if capn == 2 || capn == 4
                            lx = readgeom(
                                "LINESTRING(
                                    $(p.x+0.5) $(p.y-2.0),
                                    $(p.x+0.5) $(p.y+2.0)
                                )"
                            )
                            Bx = ilp2cap(lx, poly)
                            grid.LS[end].geoL.cap[II,6] = Bx
                        else
                            ly = readgeom(
                                "LINESTRING(
                                    $(p.x-2.0) $(p.y+0.5),
                                    $(p.x+2.0) $(p.y+0.5)
                                )"
                            )
                            By = ilp2cap(ly, poly)
                            grid.LS[end].geoL.cap[II,7] = By
                        end
                    elseif min_face[capn] != 0 && double_mixed[Acap[capn]]
                        set_A_caps_double_mixed_mixed!(capn, Acap[capn], LS1, LS2, II, poly1, poly2, p, min_face[Acap[capn]])
                        if capn == 2 || capn == 4
                            lx = readgeom(
                                "LINESTRING(
                                    $(p.x+0.5) $(p.y-2.0),
                                    $(p.x+0.5) $(p.y+2.0)
                                )"
                            )
                            try
                                Bx = ilp2cap(lx, poly)    
                            catch
                                if min_face[capn] == 1
                                    Bx = ilp2cap(lx, poly2)
                                else
                                    Bx = ilp2cap(lx, poly1)
                                end
                            end
                            
                            grid.LS[end].geoL.cap[II,6] = Bx
                        else
                            ly = readgeom(
                                "LINESTRING(
                                    $(p.x-2.0) $(p.y+0.5),
                                    $(p.x+2.0) $(p.y+0.5)
                                )"
                            )
                            try
                                By = ilp2cap(ly, poly)
                            catch
                                if min_face[capn] == 1
                                    By = ilp2cap(ly, poly2)
                                else
                                    By = ilp2cap(ly, poly1)
                                end
                            end
                            grid.LS[end].geoL.cap[II,7] = By
                        end
                    end
                end

                for (B, capn) in zip(6:7, [Bx, By])
                    set_B_caps!(capn, LS1, LS2, II, B)
                end
            else
                # The smallest keeps its value and the largest is set to 0
                if vol1 > vol2
                    LS1.geoL.cap[II,:] .= 0.0
                else
                    LS2.geoL.cap[II,:] .= 0.0
                end
            end

            if grid.LS[end].geoL.cap[II,6] < 1e-8
                grid.LS[end].geoL.cap[II,6] = 0.5 * (grid.LS[end].geoL.cap[II,1] + grid.LS[end].geoL.cap[II,3])
            end
            if grid.LS[end].geoL.cap[II,7] < 1e-8
                grid.LS[end].geoL.cap[II,7] = 0.5 * (grid.LS[end].geoL.cap[II,2] + grid.LS[end].geoL.cap[II,4])
            end
        else
            empty_cell!(grid, grid.LS[end], grid.LS[end].geoL, II)
            empty_cell!(grid, LS1, LS1.geoL, II, false)
            empty_cell!(grid, LS2, LS2.geoL, II, false)
            if vol > 1e-12
                grid.LS[end].geoL.double_emptied[II] = true
                LS1.geoL.double_emptied[II] = true
                LS2.geoL.double_emptied[II] = true
            end
        end
    end

    empty_LS1 = findall(LS1.geoL.emptied)
    empty_LS2 = findall(LS2.geoL.emptied)
    @inbounds @threads for II in empty_LS1
        empty_cell!(grid, LS2, LS2.geoL, II, false)
        empty_cell!(grid, grid.LS[end], grid.LS[end].geoL, II)
    end
    @inbounds @threads for II in empty_LS2
        empty_cell!(grid, LS1, LS1.geoL, II, false)
        empty_cell!(grid, grid.LS[end], grid.LS[end].geoL, II)
    end
    
    return nothing
end

function set_A_caps!(capn, LS1, LS2, II, poly1, poly2, p, min_face)
    if ((LS1.geoL.cap[II,capn] < (1.0 - 1e-12)) && (LS1.geoL.cap[II,capn] > 1e-12) ||
        (LS2.geoL.cap[II,capn] < (1.0 - 1e-12)) && (LS2.geoL.cap[II,capn] > 1e-12))
        iLS = findmin([LS1.geoL.cap[II,capn], LS2.geoL.cap[II,capn]])[2]
        @inbounds min_face[capn] = iLS

        if capn == 1 || capn == 3
            l = readgeom(
                "LINESTRING($(p.x+0.5) $(p.y-2.0),
                $(p.x+0.5) $(p.y+2.0))"
            )
        else
            l = readgeom(
                "LINESTRING($(p.x-2.0) $(p.y+0.5),
                $(p.x+2.0) $(p.y+0.5))"
            )
        end

        if iLS == 1
            try
                cap = ilp2cap(l, poly2)
                LS2.geoL.cap[II,capn] = cap
            catch
                cap = ilp2cap(l, poly1)
                LS2.geoL.cap[II,capn] = cap
            end
        elseif iLS == 2
            try                
                cap = ilp2cap(l, poly1)
                LS1.geoL.cap[II,capn] = cap
            catch
                cap = ilp2cap(l, poly2)
                LS1.geoL.cap[II,capn] = cap
            end
        end
    end

    return nothing
end

function set_A_caps_full_mixed!(capn, LS1, LS2, II, poly1, poly2, p, min_face)
    if capn == 1 || capn == 3
        l = readgeom(
            "LINESTRING($(p.x+0.5) $(p.y-2.0),
            $(p.x+0.5) $(p.y+2.0))"
        )
    else
        l = readgeom(
            "LINESTRING($(p.x-2.0) $(p.y+0.5),
            $(p.x+2.0) $(p.y+0.5))"
        )
    end

    if min_face == 1
        try
            cap = ilp2cap(l, poly2)
            LS1.geoL.cap[II,capn] = cap    
        catch
            cap = ilp2cap(l, poly1)
            LS1.geoL.cap[II,capn] = cap
        end
    elseif min_face == 2
        try
            cap = ilp2cap(l, poly1)
            LS2.geoL.cap[II,capn] = cap 
        catch
            cap = ilp2cap(l, poly2)
            LS2.geoL.cap[II,capn] = cap
        end
    end

    return nothing
end

function set_A_caps_double_mixed_full_empty!(capn, capm, LS1, LS2, II, poly1, p)
    if capn == 1 || capn == 3
        l = readgeom(
            "LINESTRING($(p.x+0.5) $(p.y-2.0),
            $(p.x+0.5) $(p.y+2.0))"
        )
    else
        l = readgeom(
            "LINESTRING($(p.x-2.0) $(p.y+0.5),
            $(p.x+2.0) $(p.y+0.5))"
        )
    end

    cap = ilp2cap(l, poly1)
    if LS1.geoL.cap[II,capn] > 0.5
        LS1.geoL.cap[II,capn] = 1.0
        LS1.geoL.cap[II,capm] = cap
        LS2.geoL.cap[II,capn] = cap
        LS2.geoL.cap[II,capm] = 0.0
    else
        LS1.geoL.cap[II,capn] = cap
        LS2.geoL.cap[II,capn] = 1.0 - cap
    end

    return nothing
end

function set_A_caps_full_empty!(capn, capm, LS1, LS2, II, poly1, p)
    if capn == 1 || capn == 3
        l = readgeom(
            "LINESTRING($(p.x+0.5) $(p.y-2.0),
            $(p.x+0.5) $(p.y+2.0))"
        )
    else
        l = readgeom(
            "LINESTRING($(p.x-2.0) $(p.y+0.5),
            $(p.x+2.0) $(p.y+0.5))"
        )
    end

    cap = ilp2cap(l, poly1)
    if LS1.geoL.cap[II,capn] > 0.5
        LS1.geoL.cap[II,capn] = 1.0
        LS1.geoL.cap[II,capm] = cap
        LS2.geoL.cap[II,capn] = cap
        LS2.geoL.cap[II,capm] = 0.0
    else
        LS1.geoL.cap[II,capm] = 1.0
        LS1.geoL.cap[II,capn] = cap
        LS2.geoL.cap[II,capm] = cap
        LS2.geoL.cap[II,capn] = 0.0
    end

    return nothing
end

function set_A_caps_double_mixed_mixed!(capn, capm, LS1, LS2, II, poly1, poly2, p, min_face)
    if capn == 1 || capn == 3
        l = readgeom(
            "LINESTRING($(p.x+0.5) $(p.y-2.0),
            $(p.x+0.5) $(p.y+2.0))"
        )
    else
        l = readgeom(
            "LINESTRING($(p.x-2.0) $(p.y+0.5),
            $(p.x+2.0) $(p.y+0.5))"
        )
    end

    if min_face == 1
        cap = ilp2cap(l, poly1)
        LS1.geoL.cap[II,capm] = cap
        LS2.geoL.cap[II,capm] = 0.0
    elseif min_face == 2
        cap = ilp2cap(l, poly2)
        LS2.geoL.cap[II,capm] = cap
        LS1.geoL.cap[II,capm] = 0.0
    end

    return nothing
end

function set_B_caps!(capn, LS1, LS2, II, B)
    if capn == 6
        if B > LS1.geoL.cap[II,1] && B > LS1.geoL.cap[II,3]
            LS1.geoL.cap[II,capn] = maximum([LS1.geoL.cap[II,1], LS1.geoL.cap[II,3]])
        else
            LS1.geoL.cap[II,capn] = B
        end
        if B > LS2.geoL.cap[II,1] && B > LS2.geoL.cap[II,3]
            LS2.geoL.cap[II,capn] = maximum([LS2.geoL.cap[II,1], LS2.geoL.cap[II,3]])
        else
            LS2.geoL.cap[II,capn] = B
        end
    elseif capn == 7
        if B > LS1.geoL.cap[II,2] && B > LS1.geoL.cap[II,4]
            LS1.geoL.cap[II,capn] = maximum([LS1.geoL.cap[II,2], LS1.geoL.cap[II,4]])
        else
            LS1.geoL.cap[II,capn] = B
        end
        if B > LS2.geoL.cap[II,2] && B > LS2.geoL.cap[II,4]
            LS2.geoL.cap[II,capn] = maximum([LS2.geoL.cap[II,2], LS2.geoL.cap[II,4]])
        else
            LS2.geoL.cap[II,capn] = B
        end
    end

    return nothing
end

function _marching_squares!(grid, LS, u, periodic_x, periodic_y, II, II_0, near_interface)
    @unpack x, y, nx, ny, dx, dy, ind = grid
    @unpack iso, faces, geoS, geoL, mid_point, α = LS

    empty_capacities = vcat(zeros(7), zeros(4))
    full_capacities = vcat(ones(7), 0.5.*ones(4))

    st = static_stencil(u, II_0, nx, ny, periodic_x, periodic_y)
    posW, posS, posE, posN = face_pos(II_0, II, x, y, dx[II], dy[II])
    
    a = sign(u[II])
    ISO = ifelse(a > 0, 0., 15.)
    if near_interface(a, st)
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
                geoS.cap[II,:] .= full_capacities
                geoL.cap[II,:] .= empty_capacities
            end
        else
            B, BT = B_BT(II_0, grid, periodic_x, periodic_y)
            itp = B * st * BT
            vertices = vertices_sign(itp, grid, II_0, II, dx[II], dy[II], dx[II_0], dy[II_0])
            ISO = isovalue(vertices)

            if is_not_mixed(ISO) @goto notmixed end
            face_capacities(grid, faces, itp, ISO, II_0, II, posW, posS, posE, posN)

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

function marching_squares!(grid, LS, u, periodic_x, periodic_y)
    @unpack x, y, nx, ny, dx, dy, ind = grid

    indices = vcat(
        vec(ind.inside), ind.b_left[1][2:end-1], ind.b_bottom[1][2:end-1],
        ind.b_right[1][2:end-1], ind.b_top[1][2:end-1]
    )

    if periodic_x && periodic_y
        _indices = indices
    elseif !periodic_x && periodic_y
        _indices = ind.all_indices[:,2:end-1]
    elseif periodic_x && !periodic_y
        _indices = ind.all_indices[2:end-1,:]
    else
        _indices = ind.all_indices[2:end-1,2:end-1]
    end

    if !periodic_x
        @inbounds @threads for II in ind.b_left[1][2:end-1]
            II_0 = δx⁺(II)
            _marching_squares!(grid, LS, u, periodic_x, periodic_y, II, II_0, is_near_interface_l)
        end
        @inbounds @threads for II in ind.b_right[1][2:end-1]
            II_0 = δx⁻(II)
            _marching_squares!(grid, LS, u, periodic_x, periodic_y, II, II_0, is_near_interface_r)
        end
    end
    if !periodic_y
        @inbounds @threads for II in ind.b_bottom[1][2:end-1]
            II_0 = δy⁺(II)
            _marching_squares!(grid, LS, u, periodic_x, periodic_y, II, II_0, is_near_interface_b)
        end
        @inbounds @threads for II in ind.b_top[1][2:end-1]
            II_0 = δy⁻(II)
            _marching_squares!(grid, LS, u, periodic_x, periodic_y, II, II_0, is_near_interface_t)
        end
    end

    @inbounds @threads for II in _indices
        _marching_squares!(grid, LS, u, periodic_x, periodic_y, II, II, is_near_interface)
    end

    II = ind.b_left[1][1]
    II_0 = δy⁺(δx⁺(II))
    _marching_squares!(grid, LS, u, periodic_x, periodic_y, II, II_0, is_near_interface_bl)
    
    II = ind.b_left[1][end]
    II_0 = δy⁻(δx⁺(II))
    _marching_squares!(grid, LS, u, periodic_x, periodic_y, II, II_0, is_near_interface_tl)
    
    II = ind.b_right[1][1]
    II_0 = δy⁺(δx⁻(II))
    _marching_squares!(grid, LS, u, periodic_x, periodic_y, II, II_0, is_near_interface_br)
    
    II = ind.b_right[1][end]
    II_0 = δy⁻(δx⁻(II))
    _marching_squares!(grid, LS, u, periodic_x, periodic_y, II, II_0, is_near_interface_tr)
    
    return nothing
end

function get_interface_location!(grid, LS, periodic_x, periodic_y)
    @unpack x, y = grid
    @unpack iso, geoS, geoL, mid_point, cut_points, α, MIXED = LS

    @inbounds @threads for II in MIXED
        f = average_face_capacities(grid, LS, iso, II, periodic_x, periodic_y)
        geoS.cap[II,:], geoL.cap[II,:], α[II], geoS.centroid[II], geoL.centroid[II], mid_point[II], cut_points[II], geoS.vertices[II], geoL.vertices[II] = capacities(f, iso[II])
        geoS.projection[II], geoL.projection[II] = projection_2points(grid, LS, II)
    end
    return nothing
end

function get_interface_location_borders!(grid::Mesh{GridFCx,T,N}, u, periodic_x) where {T,N}
    @unpack nx, ny, ind, geoS, geoL, mid_point, cut_points = grid
    @unpack b_left, b_bottom, b_right, b_top = ind

    f = SA_F64[0.5, 0.5]
    empty_capacities = SA_F64[vcat(zeros(7), zeros(4))...]
    capacities_6 = capacities(f, 6.0)
    capacities_9 = capacities(f, 9.0)

    if !periodic_x
        @inbounds @threads for II in b_left[1]
            if u[II] >= 0.0
                @inbounds geoS.cap[II,:] .= empty_capacities
                @inbounds _, geoL.cap[II,:], _, _, _, _, _, _, _ = capacities_9
            else
                @inbounds geoS.cap[II,:], _, _, _, _, _, _, _, _ = capacities_6
                @inbounds geoL.cap[II,:] .= empty_capacities
            end
        end
        @inbounds @threads for II in b_right[1]
            if u[II] >= 0.0
                @inbounds geoS.cap[II,:] .= empty_capacities
                @inbounds _, geoL.cap[II,:], _, _, _, _, _, _, _ = capacities_6
            else
                @inbounds geoS.cap[II,:], _, _, _, _, _, _, _, _ = capacities_9
                @inbounds geoL.cap[II,:] .= empty_capacities
            end
        end
    end

    return nothing
end

function get_interface_location_borders!(grid::Mesh{GridFCy,T,N}, u, periodic_y) where {T,N}
    @unpack ind, geoS, geoL, mid_point, cut_points = grid
    @unpack b_left, b_bottom, b_right, b_top = ind

    f = SA_F64[0.5, 0.5]
    empty_capacities = SA_F64[vcat(zeros(7), zeros(4))...]
    capacities_3 = capacities(f, 3.0)
    capacities_12 = capacities(f, 12.0)

    if !periodic_y
        @inbounds @threads for II in b_bottom[1]
            if u[II] >= 0.0
                @inbounds geoS.cap[II,:] .= empty_capacities
                @inbounds _, geoL.cap[II,:], _, _, _, _, _, _, _ = capacities_3
            else
                @inbounds geoS.cap[II,:], _, _, _, _, _, _, _, _ = capacities_12
                @inbounds geoL.cap[II,:] .= empty_capacities
            end
        end
        @inbounds @threads for II in b_top[1]
            if u[II] >= 0.0
                @inbounds geoS.cap[II,:] .= empty_capacities
                @inbounds _, geoL.cap[II,:], _, _, _, _, _, _, _ = capacities_12
            else
                @inbounds geoS.cap[II,:], _, _, _, _, _, _, _, _ = capacities_3
                @inbounds geoL.cap[II,:] .= empty_capacities
            end
        end
    end

    return nothing
end

function get_curvature(num, grid, geoL, u, κ, inside, per_x, per_y)
    @unpack Δ = num
    @unpack x, y, nx, ny, ind = grid

    κ .= zeros(grid)
    @inbounds @threads for II in inside
        if !per_x && !per_y
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
            else
                continue
            end
        elseif per_x && !per_y
            if II in ind.inside || II in ind.b_left[1][2:end-1] || II in ind.b_right[1][2:end-1]
                II_0 = II
            elseif II in ind.b_bottom[1]
                II_0 = δy⁺(II)
            elseif II in ind.b_top[1]
                II_0 = δy⁻(II)
            else 
                continue
            end
        else
            if II in ind.inside || II in ind.b_bottom[1][2:end-1] || II in ind.b_top[1][2:end-1]
                II_0 = II
            elseif II in ind.b_left[1]
                II_0 = δx⁺(II)
            elseif II in ind.b_right[1]
                II_0 = δx⁻(II)
            else
                continue
            end
        end
        mid_point = geoL.projection[II].mid_point
        st = static_stencil(u, II_0, nx, ny, per_x, per_y)
        B, BT = B_BT(II_0, grid, per_x, per_y)
        itp = B * st * BT
        
        dx = grid.dx[II]
        dy = grid.dy[II]

        κ[II] = parabola_fit_curvature(itp, mid_point, dx, dy)
    end
end

function capacities(F_prev, case)
    # add some little volume so that weird cases don't appear
    F = [F_prev[1], F_prev[2]]
    # if F[1] <= 1e-8
    #     F[1] += 1e-8
    # elseif  F[1] >= (1.0-1e-8)
    #     F[1] -= 1e-8
    # end
    # if F[2] <= 1e-8
    #     F[2] += 1e-8
    # elseif  F[2] >= (1.0-1e-8)
    #     F[2] -= 1e-8
    # end

    cap_sol = @SVector zeros(11)
    cap_liq = @SVector zeros(11)
    α = 0.0
    sol_centroid = Point(NaN, NaN)
    liq_centroid = Point(NaN, NaN)
    mid_point = Point(NaN, NaN)
    cut_points = [Point(NaN, NaN), Point(NaN, NaN)]
    vertS = [Point(NaN, NaN)]
    vertL = [Point(NaN, NaN)]
    if case == 1.0 #W S
        vertS = [Point(0.,0.), Point(F[2],0.), Point(0.,F[1]), Point(0.,0.)]
        vertL = [Point(1.,0.), Point(F[2],0.), Point(0., F[1]), Point(0.,1.), Point(1.,1.), Point(1., 0.)]
        sol_centroid = get_centroid(vertS)
        liq_centroid = get_centroid(vertL)
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
        vertS = [Point(F[1],0.), Point(1.,0.), Point(1.,F[2]), Point(F[1],0.)]
        vertL = [Point(0.,0.), Point(F[1],0.), Point(1., F[2]), Point(1.,1.), Point(0.,1.), Point(0., 0.)]
        sol_centroid = get_centroid(vertS)
        liq_centroid = get_centroid(vertL)
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
        vertS = [Point(0.,0.), Point(1.,0.), Point(1.,F[2]), Point(0.,F[1]), Point(0.,0.)]
        vertL = [Point(1.,1.), Point(0.,1.), Point(0.,F[1]), Point(1.,F[2]), Point(1.,1.)]
        sol_centroid = get_centroid(vertS)
        liq_centroid = get_centroid(vertL)
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
        vertS = [Point(1.,1.), Point(F[2],1.), Point(1.,F[1]), Point(1.,1.)]
        vertL = [Point(0.,0.), Point(1.,0.), Point(1.,F[1]), Point(F[2],1.), Point(0.,1.), Point(0.,0.)]
        sol_centroid = get_centroid(vertS)
        liq_centroid = get_centroid(vertL)
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
        vertS = [Point(1.,1.), Point(F[2],1.), Point(F[1],0.), Point(1.,0.), Point(1.,1.)]
        vertL = [Point(0.,0.), Point(F[1],0.), Point(F[2],1.), Point(0.,1.), Point(0.,0.)]
        sol_centroid = get_centroid(vertS)
        liq_centroid = get_centroid(vertL)
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
        vertS = [Point(0.,0.), Point(1.,0.), Point(1.,1.), Point(F[2],1.), Point(0.,F[1]), Point(0.,0.)]
        vertL = [Point(0.,1.), Point(0.,F[1]), Point(F[2],1.), Point(0.,1.)]
        sol_centroid = get_centroid(vertS)
        liq_centroid = get_centroid(vertL)
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
        vertS = [Point(0.,1.), Point(0.,F[1]), Point(F[2],1.), Point(0.,1.)]
        vertL = [Point(0.,0.), Point(1.,0.), Point(1.,1.), Point(F[2],1.), Point(0.,F[1]), Point(0.,0.)]
        sol_centroid = get_centroid(vertS)
        liq_centroid = get_centroid(vertL)
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
        vertS = [Point(0.,0.), Point(F[1],0.), Point(F[2],1.), Point(0.,1.), Point(0.,0.)]
        vertL = [Point(1.,1.), Point(F[2],1.), Point(F[1],0.), Point(1.,0.), Point(1.,1.)]
        sol_centroid = get_centroid(vertS)
        liq_centroid = get_centroid(vertL)
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
        vertS = [Point(0.,0.), Point(1.,0.), Point(1.,F[1]), Point(F[2],1.), Point(0.,1.), Point(0.,0.)]
        vertL = [Point(1.,1.), Point(F[2],1.), Point(1.,F[1]), Point(1.,1.)]
        sol_centroid = get_centroid(vertS)
        liq_centroid = get_centroid(vertL)
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
        vertS = [Point(1.,1.), Point(0.,1.), Point(0.,F[1]), Point(1.,F[2]), Point(1.,1.)]
        vertL = [Point(0.,0.), Point(1.,0.), Point(1.,F[2]), Point(0.,F[1]), Point(0.,0.)]
        sol_centroid = get_centroid(vertS)
        liq_centroid = get_centroid(vertL)
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
        vertS = [Point(0.,0.), Point(F[1],0.), Point(1.,F[2]), Point(1.,1.), Point(0.,1.), Point(0.,0.)]
        vertL = [Point(1.,0.), Point(1.,F[2]), Point(F[1],0.), Point(1.,0.)]
        sol_centroid = get_centroid(vertS)
        liq_centroid = get_centroid(vertL)
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
        vertS = [Point(1.,1.), Point(0.,1.), Point(0.,F[1]), Point(F[2],0.), Point(1.,0.), Point(1.,1.)]
        vertL = [Point(0.,0.), Point(F[2],0.), Point(0.,F[1]), Point(0.,0.)]
        sol_centroid = get_centroid(vertS)
        liq_centroid = get_centroid(vertL)
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
    return cap_sol, cap_liq, min(float(π),max(-float(π),α)), sol_centroid, liq_centroid, mid_point, cut_points, vertS, vertL
end

function set_cap_bcs!(grid::Mesh{GridCC,T,N}, num, LS, periodic_x, periodic_y, empty = true) where {T,N}
    @unpack nx, ny = grid
    @unpack geoS, geoL, mid_point = LS

    # set A and mid_point at the boundaries to 0 if not periodic in that direction
    if periodic_x
        @inbounds @threads for i = 1:ny
            @inbounds geoS.cap[i,1,1] = geoS.cap[i,end,3]
            @inbounds geoL.cap[i,1,1] = geoL.cap[i,end,3]
        end
    elseif empty
        @inbounds @threads for i = 1:ny
            @inbounds geoS.cap[i,1,1] = 0.0
            @inbounds geoS.cap[i,end,3] = 0.0
            @inbounds geoL.cap[i,1,1] = 0.0
            @inbounds geoL.cap[i,end,3] = 0.0
            @inbounds mid_point[i,1] += Point(-0.5, 0.0)
            @inbounds mid_point[i,end] += Point(0.5, 0.0)
        end
    end
    if periodic_y
        @inbounds @threads for i = 1:nx
            @inbounds geoS.cap[1,i,2] = geoS.cap[end,i,4]
            @inbounds geoL.cap[1,i,2] = geoL.cap[end,i,4]
        end
    elseif empty
        @inbounds @threads for i = 1:nx
            @inbounds geoS.cap[1,i,2] = 0.0
            @inbounds geoS.cap[end,i,4] = 0.0
            @inbounds geoL.cap[1,i,2] = 0.0
            @inbounds geoL.cap[end,i,4] = 0.0
            @inbounds mid_point[1,i] += Point(0.0, -0.5)
            @inbounds mid_point[end,i] += Point(0.0, 0.5)
        end
    end

    return nothing
end

function set_cap_diff_S_L!(num,geoS,geoL,tmpS,tmpL,cap_id,II)

    diffS = geoS.cap[II,cap_id] - tmpS
    if diffS > 0.0
        @inbounds geoS.cap[II,5] = geoS.cap[II,5] * diffS * inv_weight_eps(num,geoS.cap[II,cap_id])
    else
        @inbounds geoS.cap[II,5] = geoS.cap[II,5] * 0.5
    end
    diffL = geoL.cap[II,cap_id] - tmpL
    if diffL > 0.0
        @inbounds geoL.cap[II,5] = geoL.cap[II,5] * diffL * inv_weight_eps(num,geoL.cap[II,cap_id]) 
    else
        @inbounds geoL.cap[II,5] = geoL.cap[II,5] * 0.5
    end
end

function set_cap_bcs!(grid::Mesh{GridFCx,T,N}, num, LS, periodic_x, periodic_y, empty = true) where {T,N}
    @unpack nx, ny, ind = grid
    @unpack b_left, b_bottom, b_right, b_top = ind
    @unpack geoS, geoL, mid_point, cut_points, iso = LS

    # set A at the boundaries
    if periodic_y
        @inbounds @threads for i = 1:nx
            @inbounds geoS.cap[1,i,2] = geoS.cap[end,i,4]
            @inbounds geoL.cap[1,i,2] = geoL.cap[end,i,4]
        end
    elseif empty
        @inbounds @threads for i = 1:nx
            @inbounds geoS.cap[1,i,2] = 0.0
            @inbounds geoS.cap[end,i,4] = 0.0
            @inbounds geoL.cap[1,i,2] = 0.0
            @inbounds geoL.cap[end,i,4] = 0.0
            @inbounds mid_point[1,i] += Point(0.0, -0.5)
            @inbounds mid_point[end,i] += Point(0.0, 0.5)
        end
    end

    # set mid_point in the outer boundaries
    if periodic_x
        @inbounds @threads for i = 1:ny
            @inbounds geoS.cap[i,1,1] = geoS.cap[i,end,3]
            @inbounds geoL.cap[i,1,1] = geoL.cap[i,end,3]
        end
    else
        @inbounds @threads for II in b_left[1]
            @inbounds geoS.cap[II,1] = 0.0
            @inbounds geoL.cap[II,1] = 0.0
            if iso[II] == 0
                @inbounds geoL.cap[II,2] = 0.5
                @inbounds geoL.cap[II,4] = 0.5
                tmpS = 0.0
                tmpL = 0.5
                @inbounds geoL.centroid[II] = geoL.centroid[II] + Point(0.25,0.0)
            elseif iso[II] == 15
                @inbounds geoS.cap[II,2] = 0.5
                @inbounds geoS.cap[II,4] = 0.5
                tmpS = 0.5
                tmpL = 0.0
                @inbounds geoS.centroid[II] = geoS.centroid[II] + Point(0.25,0.0)
            end
            if iso[II] == 1 || iso[II] == 9 || iso[II] == 13
                @inbounds geoS.cap[II,2] = max(geoS.cap[II,2] - 0.5, 0.0)
                @inbounds geoL.cap[II,2] = min(geoL.cap[II,2], 0.5)
                tmpS = max(geoS.cap[II,7] - 0.5, 0.0)
                tmpL = min(geoL.cap[II,7], 0.5)
            elseif iso[II] == 2 || iso[II] == 6 || iso[II] == 14
                @inbounds geoS.cap[II,2] = min(geoS.cap[II,2], 0.5)
                @inbounds geoL.cap[II,2] = max(geoL.cap[II,2] - 0.5, 0.0)
                tmpS = min(geoS.cap[II,7], 0.5)
                tmpL = max(geoL.cap[II,7] - 0.5, 0.0)
            elseif iso[II] == 3 || iso[II] == 4 || iso[II] == 7 || iso[II] == 8 || iso[II] == 11 || iso[II] == 12
                @inbounds geoS.cap[II,2] *= 0.5
                @inbounds geoL.cap[II,2] *= 0.5
            end
            if iso[II] == 4 || iso[II] == 6 || iso[II] == 7
                @inbounds geoS.cap[II,4] = min(geoS.cap[II,4], 0.5)
                @inbounds geoL.cap[II,4] = max(geoL.cap[II,4] - 0.5, 0.0)
                tmpS = min(geoS.cap[II,7], 0.5)
                tmpL = max(geoL.cap[II,7] - 0.5, 0.0)
            elseif iso[II] == 8 || iso[II] == 9 || iso[II] == 11
                @inbounds geoS.cap[II,4] = max(geoS.cap[II,4] - 0.5, 0.0)
                @inbounds geoL.cap[II,4] = min(geoL.cap[II,4], 0.5)
                tmpS = max(geoS.cap[II,7] - 0.5, 0.0)
                tmpL = min(geoL.cap[II,7], 0.5)
            elseif iso[II] == 1 || iso[II] == 2 || iso[II] == 3 || iso[II] == 12 || iso[II] == 13 || iso[II] == 14
                @inbounds geoS.cap[II,4] *= 0.5
                @inbounds geoL.cap[II,4] *= 0.5
            end
            if iso[II] == 3
                tmpS = min(geoS.cap[II,7], 0.5)
                tmpL = max(geoL.cap[II,7] - 0.5, 0.0)
            elseif iso[II] == 12
                tmpS = max(geoS.cap[II,7] - 0.5, 0.0)
                tmpL = min(geoL.cap[II,7], 0.5)
            end

            set_cap_diff_S_L!(num,geoS,geoL,tmpS,tmpL,7,II)

            @inbounds geoS.cap[II,7] = tmpS
            @inbounds geoL.cap[II,7] = tmpL
            @inbounds geoS.cap[II,8:11] .= geoS.cap[II,5] * 0.5
            @inbounds geoL.cap[II,8:11] .= geoL.cap[II,5] * 0.5
            @inbounds geoS.centroid[II] = Point(0.5 * geoS.cap[II,7],0.0)
            @inbounds geoL.centroid[II] = Point(0.5 * geoL.cap[II,7],0.0)
        end
        @inbounds @threads for II in b_right[1]
            @inbounds geoS.cap[II,3] = 0.0
            @inbounds geoL.cap[II,3] = 0.0
            if iso[II] == 0
                @inbounds geoL.cap[II,2] = 0.5
                @inbounds geoL.cap[II,4] = 0.5
                tmpS = 0.0
                tmpL = 0.5
                @inbounds geoL.centroid[II] = geoL.centroid[II] - Point(0.25,0.0)
            elseif iso[II] == 15
                @inbounds geoS.cap[II,2] = 0.5
                @inbounds geoS.cap[II,4] = 0.5
                tmpS = 0.5
                tmpL = 0.0
                @inbounds geoS.centroid[II] = geoS.centroid[II] - Point(0.25,0.0)
            end
            if iso[II] == 1 || iso[II] == 9 || iso[II] == 13
                @inbounds geoS.cap[II,2] = min(geoS.cap[II,2], 0.5)
                @inbounds geoL.cap[II,2] = max(geoL.cap[II,2] - 0.5, 0.0)
                tmpS = min(geoS.cap[II,7], 0.5)
                tmpL = max(geoL.cap[II,7] - 0.5, 0.0)
            elseif iso[II] == 2 || iso[II] == 6 || iso[II] == 14
                @inbounds geoS.cap[II,2] = max(geoS.cap[II,2] - 0.5, 0.0)
                @inbounds geoL.cap[II,2] = min(geoL.cap[II,2], 0.5)
                tmpS = max(geoS.cap[II,7] - 0.5, 0.0)
                tmpL = min(geoL.cap[II,7], 0.5)
            elseif iso[II] == 3 || iso[II] == 4 || iso[II] == 7 || iso[II] == 8 || iso[II] == 11 || iso[II] == 12
                @inbounds geoS.cap[II,2] *= 0.5
                @inbounds geoL.cap[II,2] *= 0.5
            end
            if iso[II] == 4 || iso[II] == 6 || iso[II] == 7
                @inbounds geoS.cap[II,4] = max(geoS.cap[II,4] - 0.5, 0.0)
                @inbounds geoL.cap[II,4] = min(geoL.cap[II,4], 0.5)
                tmpS = max(geoS.cap[II,7] - 0.5, 0.0)
                tmpL = min(geoL.cap[II,7], 0.5)
            elseif iso[II] == 8 || iso[II] == 9 || iso[II] == 11
                @inbounds geoS.cap[II,4] = min(geoS.cap[II,4], 0.5)
                @inbounds geoL.cap[II,4] = max(geoL.cap[II,4] - 0.5, 0.0)
                tmpS = min(geoS.cap[II,7], 0.5)
                tmpL = max(geoL.cap[II,7] - 0.5, 0.0)
            elseif iso[II] == 1 || iso[II] == 2 || iso[II] == 3 || iso[II] == 12 || iso[II] == 13 || iso[II] == 14
                @inbounds geoS.cap[II,4] *= 0.5
                @inbounds geoL.cap[II,4] *= 0.5
            end
            if iso[II] == 3
                tmpS = max(geoS.cap[II,7] - 0.5, 0.0)
                tmpL = min(geoL.cap[II,7], 0.5)
            elseif iso[II] == 12
                tmpS = min(geoS.cap[II,7], 0.5)
                tmpL = max(geoL.cap[II,7] - 0.5, 0.0)
            end

            set_cap_diff_S_L!(num,geoS,geoL,tmpS,tmpL,7,II)

            @inbounds geoS.cap[II,7] = tmpS
            @inbounds geoL.cap[II,7] = tmpL
            @inbounds geoS.cap[II,8:11] .= geoS.cap[II,5] * 0.5
            @inbounds geoL.cap[II,8:11] .= geoL.cap[II,5] * 0.5
            @inbounds geoS.centroid[II] = Point(-0.5 * geoS.cap[II,7],0.0)
            @inbounds geoL.centroid[II] = Point(-0.5 * geoL.cap[II,7],0.0)
        end
    end

    return nothing
end

function set_cap_bcs!(grid::Mesh{GridFCy,T,N}, num, LS, periodic_x, periodic_y, empty = true) where {T,N}
    @unpack nx, ny, ind = grid
    @unpack b_left, b_bottom, b_right, b_top = ind
    @unpack geoS, geoL, mid_point, cut_points, iso = LS


    # set A at the boundaries to 0 if not periodic in that direction
    if periodic_x
        @inbounds @threads for i = 1:ny
            @inbounds geoS.cap[i,1,1] = geoS.cap[i,end,3]
            @inbounds geoL.cap[i,1,1] = geoL.cap[i,end,3]
        end
    elseif empty
        @inbounds @threads for i = 1:ny
            @inbounds geoS.cap[i,1,1] = 0.0
            @inbounds geoS.cap[i,end,3] = 0.0
            @inbounds geoL.cap[i,1,1] = 0.0
            @inbounds geoL.cap[i,end,3] = 0.0
            @inbounds mid_point[i,1] += Point(-0.5, 0.0)
            @inbounds mid_point[i,end] += Point(0.5, 0.0)
        end
    end

    # set mid_point in the outer boundaries
    if periodic_y
        @inbounds @threads for i = 1:nx
            @inbounds geoS.cap[1,i,2] = geoS.cap[end,i,4]
            @inbounds geoL.cap[1,i,2] = geoL.cap[end,i,4]
        end
    else
        @inbounds @threads for II in b_bottom[1]
            @inbounds geoS.cap[II,2] = 0.0
            @inbounds geoL.cap[II,2] = 0.0
            if iso[II] == 0
                @inbounds geoL.cap[II,1] = 0.5
                @inbounds geoL.cap[II,3] = 0.5
                tmpS = 0.0
                tmpL = 0.5
                @inbounds geoL.centroid[II] = geoL.centroid[II] + Point(0.0,0.25)
            elseif iso[II] == 15
                @inbounds geoS.cap[II,1] = 0.5
                @inbounds geoS.cap[II,3] = 0.5
                tmpS = 0.5
                tmpL = 0.0
                @inbounds geoS.centroid[II] = geoS.centroid[II] + Point(0.0,0.25)
            end
            if iso[II] == 1 || iso[II] == 3 || iso[II] == 7
                @inbounds geoS.cap[II,1] = max(geoS.cap[II,1] - 0.5, 0.0)
                @inbounds geoL.cap[II,1] = min(geoL.cap[II,1], 0.5)
                tmpS = max(geoS.cap[II,6] - 0.5, 0.0)
                tmpL = min(geoL.cap[II,6], 0.5)
            elseif iso[II] == 8 || iso[II] == 12 || iso[II] == 14
                @inbounds geoS.cap[II,1] = min(geoS.cap[II,1], 0.5)
                @inbounds geoL.cap[II,1] = max(geoL.cap[II,1] - 0.5, 0.0)
                tmpS = min(geoS.cap[II,6], 0.5)
                tmpL = max(geoL.cap[II,6] - 0.5, 0.0)
            elseif iso[II] == 2 || iso[II] == 4 || iso[II] == 6 || iso[II] == 9 || iso[II] == 11 || iso[II] == 13
                @inbounds geoS.cap[II,1] *= 0.5
                @inbounds geoL.cap[II,1] *= 0.5
            end
            if iso[II] == 4 || iso[II] == 12 || iso[II] == 13
                @inbounds geoS.cap[II,3] = min(geoS.cap[II,3], 0.5)
                @inbounds geoL.cap[II,3] = max(geoL.cap[II,3] - 0.5, 0.0)
                tmpS = min(geoS.cap[II,6], 0.5)
                tmpL = max(geoL.cap[II,6] - 0.5, 0.0)
            elseif iso[II] == 2 || iso[II] == 3 || iso[II] == 11
                @inbounds geoS.cap[II,3] = max(geoS.cap[II,3] - 0.5, 0.0)
                @inbounds geoL.cap[II,3] = min(geoL.cap[II,3], 0.5)
                tmpS = max(geoS.cap[II,6] - 0.5, 0.0)
                tmpL = min(geoL.cap[II,6], 0.5)
            elseif iso[II] == 1 || iso[II] == 6 || iso[II] == 7 || iso[II] == 8 || iso[II] == 9 || iso[II] == 14
                @inbounds geoS.cap[II,3] *= 0.5
                @inbounds geoL.cap[II,3] *= 0.5
            end
            if iso[II] == 9
                tmpS = min(geoS.cap[II,6], 0.5)
                tmpL = max(geoL.cap[II,6] - 0.5, 0.0)
            elseif iso[II] == 6
                tmpS = max(geoS.cap[II,6] - 0.5, 0.0)
                tmpL = min(geoL.cap[II,6], 0.5)
            end

            set_cap_diff_S_L!(num,geoS,geoL,tmpS,tmpL,6,II)

            @inbounds geoS.cap[II,6] = tmpS
            @inbounds geoL.cap[II,6] = tmpL
            @inbounds geoS.cap[II,8:11] .= geoS.cap[II,5] * 0.5
            @inbounds geoL.cap[II,8:11] .= geoL.cap[II,5] * 0.5
            @inbounds geoS.centroid[II] = Point(0.0, 0.5 * geoS.cap[II,6])
            @inbounds geoL.centroid[II] = Point(0.0, 0.5 * geoL.cap[II,6])
        end
        @inbounds @threads for II in b_top[1]
            @inbounds geoS.cap[II,4] = 0.0
            @inbounds geoL.cap[II,4] = 0.0
            if iso[II] == 0
                @inbounds geoL.cap[II,1] = 0.5
                @inbounds geoL.cap[II,3] = 0.5
                tmpS = 0.0
                tmpL = 0.5
                @inbounds geoL.centroid[II] = geoL.centroid[II] - Point(0.0,0.25)
            elseif iso[II] == 15
                @inbounds geoS.cap[II,1] = 0.5
                @inbounds geoS.cap[II,3] = 0.5
                tmpS = 0.5
                tmpL = 0.0
                @inbounds geoS.centroid[II] = geoS.centroid[II] - Point(0.0,0.25)
            end
            if iso[II] == 1 || iso[II] == 3 || iso[II] == 7
                @inbounds geoS.cap[II,1] = min(geoS.cap[II,1], 0.5)
                @inbounds geoL.cap[II,1] = max(geoL.cap[II,1] - 0.5, 0.0)
                tmpS = min(geoS.cap[II,6], 0.5)
                tmpL = max(geoL.cap[II,6] - 0.5, 0.0)
            elseif iso[II] == 8 || iso[II] == 12 || iso[II] == 14
                @inbounds geoS.cap[II,1] = max(geoS.cap[II,1] - 0.5, 0.0)
                @inbounds geoL.cap[II,1] = min(geoL.cap[II,1], 0.5)
                tmpS = max(geoS.cap[II,6] - 0.5, 0.0)
                tmpL = min(geoL.cap[II,6], 0.5)
            elseif iso[II] == 2 || iso[II] == 4 || iso[II] == 6 || iso[II] == 9 || iso[II] == 11 || iso[II] == 13
                @inbounds geoS.cap[II,1] *= 0.5
                @inbounds geoL.cap[II,1] *= 0.5
            end
            if iso[II] == 4 || iso[II] == 12 || iso[II] == 13
                @inbounds geoS.cap[II,3] = max(geoS.cap[II,3] - 0.5, 0.0)
                @inbounds geoL.cap[II,3] = min(geoL.cap[II,3], 0.5)
                tmpS = max(geoS.cap[II,6] - 0.5, 0.0)
                tmpL = min(geoL.cap[II,6], 0.5)
            elseif iso[II] == 2 || iso[II] == 3 || iso[II] == 11
                @inbounds geoS.cap[II,3] = min(geoS.cap[II,3], 0.5)
                @inbounds geoL.cap[II,3] = max(geoL.cap[II,3] - 0.5, 0.0)
                tmpS = min(geoS.cap[II,6], 0.5)
                tmpL = max(geoL.cap[II,6] - 0.5, 0.0)
            elseif iso[II] == 1 || iso[II] == 6 || iso[II] == 7 || iso[II] == 8 || iso[II] == 9 || iso[II] == 14
                @inbounds geoS.cap[II,3] *= 0.5
                @inbounds geoL.cap[II,3] *= 0.5
            end
            if iso[II] == 9
                tmpS = max(geoS.cap[II,6] - 0.5, 0.0)
                tmpL = min(geoL.cap[II,6], 0.5)
            elseif iso[II] == 6
                tmpS = min(geoS.cap[II,6], 0.5)
                tmpL = max(geoL.cap[II,6] - 0.5, 0.0)
            end

            set_cap_diff_S_L!(num,geoS,geoL,tmpS,tmpL,6,II)

            @inbounds geoS.cap[II,6] = tmpS
            @inbounds geoL.cap[II,6] = tmpL
            @inbounds geoS.cap[II,8:11] .= geoS.cap[II,5] * 0.5
            @inbounds geoL.cap[II,8:11] .= geoL.cap[II,5] * 0.5
            @inbounds geoS.centroid[II] = Point(0.0,-0.5 * geoS.cap[II,6])
            @inbounds geoL.centroid[II] = Point(0.0,-0.5 * geoL.cap[II,6])
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

# function Wcapacities!(cap, periodic_x, periodic_y)
#     tmp = copy(cap) #TODO large allocation

#     @inbounds cap[:,2:end,8] .= @views tmp[:,2:end,8] .+ tmp[:,1:end-1,10]
#     @inbounds cap[:,1:end-1,10] .= @views tmp[:,1:end-1,10] .+ tmp[:,2:end,8]
#     @inbounds cap[2:end,:,9] .= @views tmp[2:end,:,9] .+ tmp[1:end-1,:,11]
#     @inbounds cap[1:end-1,:,11] .= @views tmp[1:end-1,:,11] .+ tmp[2:end,:,9]

#     if periodic_x
#         @inbounds cap[:,1,8] .= @views tmp[:,1,8] .+ tmp[:,end,10]
#         @inbounds cap[:,end,10] .= @views tmp[:,end,10] .+ tmp[:,1,8]
#     end
#     if periodic_y
#         @inbounds cap[1,:,9] .= @views tmp[1,:,9] .+ tmp[end,:,11]
#         @inbounds cap[end,:,11] .= @views tmp[end,:,11] .+ tmp[1,:,9]
#     end

#     return nothing
# end

function Wcapacities!(cap, periodic_x, periodic_y)
    #avoid tmp = copy(cap) #TODO large allocation
    #indices do not overlap
    @views @inbounds cap[:,2:end,8] .+= cap[:,1:end-1,10]
    @views @inbounds cap[:,1:end-1,10] .= cap[:,2:end,8]
    @views @inbounds cap[2:end,:,9] .+= cap[1:end-1,:,11]
    @views @inbounds cap[1:end-1,:,11] .=  cap[2:end,:,9]

    if periodic_x
        @views @inbounds cap[:,1,8] .+= cap[:,end,10]
        @views @inbounds cap[:,end,10] .= cap[:,1,8]
    end
    if periodic_y
        @views @inbounds cap[1,:,9] .+= cap[end,:,11]
        @views @inbounds cap[end,:,11] .= cap[1,:,9]
    end

    return nothing
end

function face_capacities(grid, faces, itp, case, II_0, II, posW, posS, posE, posN)
    @unpack nx, ny, dx, dy = grid

    if case == 1.0 || case == 14.0
        faces[II, 1] = ispositive(WEST_face(itp, posW, dy[II], dx[II_0], dy[II_0])) / dy[II]
        faces[II, 2] = ispositive(SOUTH_face(itp, posS, dx[II], dx[II_0], dy[II_0])) / dx[II]
    elseif case == 2.0 || case == 13.0
        faces[II, 2] = ispositive(SOUTH_face(itp, posS, dx[II], dx[II_0], dy[II_0])) / dx[II]
        faces[II, 3] = ispositive(EAST_face(itp, posE, dy[II], dx[II_0], dy[II_0])) / dy[II]
    elseif case == 3.0 || case == 12.0
        faces[II, 1] = ispositive(WEST_face(itp, posW, dy[II], dx[II_0], dy[II_0])) / dy[II]
        faces[II, 3] = ispositive(EAST_face(itp, posE, dy[II], dx[II_0], dy[II_0])) / dy[II]
    elseif case == 4.0 || case == 11.0
        faces[II, 3] = ispositive(EAST_face(itp, posE, dy[II], dx[II_0], dy[II_0])) / dy[II]
        faces[II, 4] = ispositive(NORTH_face(itp, posN, dx[II], dx[II_0], dy[II_0])) / dx[II]
    elseif case == 6.0 || case == 9.0
        faces[II, 2] = ispositive(SOUTH_face(itp, posS, dx[II], dx[II_0], dy[II_0])) / dx[II]
        faces[II, 4] = ispositive(NORTH_face(itp, posN, dx[II], dx[II_0], dy[II_0])) / dx[II]
    elseif case == 7.0 || case == 8.0
        faces[II, 1] = ispositive(WEST_face(itp, posW, dy[II], dx[II_0], dy[II_0])) / dy[II]
        faces[II, 4] = ispositive(NORTH_face(itp, posN, dx[II], dx[II_0], dy[II_0])) / dx[II]
    end
end

function average_face_capacities(grid, LS, iso, II, per_x, per_y)
    @unpack nx, ny, ind = grid
    @unpack b_left, b_bottom, b_right, b_top = ind

    case = iso[II]

    if II in b_left[1] && !per_x
        WE_fun = WE_border
    else
        WE_fun = WE
    end
    if II in b_bottom[1] && !per_y
        SN_fun = SN_border
    else
        SN_fun = SN
    end
    if II in b_right[1] && !per_x
        EW_fun = EW_border
    else
        EW_fun = EW
    end
    if II in b_top[1] && !per_y
        NS_fun = NS_border
    else
        NS_fun = NS
    end
    if case == 1.0
        f = SA_F64[WE_fun(grid, LS, II, 0.0, 1.0, per_x), SN_fun(grid, LS, II, 0.0, 1.0, per_y)]
    elseif case == 14.0
        f = SA_F64[WE_fun(grid, LS, II, 1.0, 0.0, per_x), SN_fun(grid, LS, II, 1.0, 0.0, per_y)]
    elseif case == 2.0
        f = SA_F64[SN_fun(grid, LS, II, 1.0, 0.0, per_y), EW_fun(grid, LS, II, 0.0, 1.0, per_x)]
    elseif case == 13.0
        f = SA_F64[SN_fun(grid, LS, II, 0.0, 1.0, per_y), EW_fun(grid, LS, II, 1.0, 0.0, per_x)]
    elseif case == 3.0
        f = SA_F64[WE_fun(grid, LS, II, 0.0, 1.0, per_x), EW_fun(grid, LS, II, 0.0, 1.0, per_x)]
    elseif case == 12.0
        f = SA_F64[WE_fun(grid, LS, II, 1.0, 0.0, per_x), EW_fun(grid, LS, II, 1.0, 0.0, per_x)]
    elseif case == 4.0
        f = SA_F64[EW_fun(grid, LS, II, 1.0, 0.0, per_x), NS_fun(grid, LS, II, 1.0, 0.0, per_y)]
    elseif case == 11.0
        f = SA_F64[EW_fun(grid, LS, II, 0.0, 1.0, per_x), NS_fun(grid, LS, II, 0.0, 1.0, per_y)]
    elseif case == 6.0
        f = SA_F64[SN_fun(grid, LS, II, 1.0, 0.0, per_y), NS_fun(grid, LS, II, 1.0, 0.0, per_y)]
    elseif case == 9.0
        f = SA_F64[SN_fun(grid, LS, II, 0.0, 1.0, per_y), NS_fun(grid, LS, II, 0.0, 1.0, per_y)]
    elseif case == 7.0
        f = SA_F64[WE_fun(grid, LS, II, 0.0, 1.0, per_x), NS_fun(grid, LS, II, 1.0, 0.0, per_y)]
    elseif case == 8.0
        f = SA_F64[WE_fun(grid, LS, II, 1.0, 0.0, per_x), NS_fun(grid, LS, II, 0.0, 1.0, per_y)]
    else
        f = @SVector zeros(2)
    end
    return float(f)
end

function average_face_capacities!(a, per_x, per_y)
    # average face capacities due to mismatch in iso 
    n = size(a,1)
    @inbounds @threads for i = 1:n
        @inbounds tmp1 = @view a[i,1:end-1,3]
        @inbounds tmp2 = @view a[i,2:end,1]
        @inbounds avg_cap = 0.5*(tmp1 .+ tmp2)
        @inbounds a[i,1:end-1,3] .= avg_cap
        @inbounds a[i,2:end,1] .= avg_cap

        if per_x
            @inbounds avg_border = 0.5*(a[i,1,1] + a[i,end,3])
            @inbounds a[i,1,1] = avg_border
            @inbounds a[i,end,3] = avg_border
        end
    end
    n = size(a,2)
    @inbounds @threads for i = 1:n
        @inbounds tmp1 = @view a[1:end-1,i,4]
        @inbounds tmp2 = @view a[2:end,i,2]
        @inbounds avg_cap = 0.5*(tmp1 .+ tmp2)
        @inbounds a[1:end-1,i,4] .= avg_cap
        @inbounds a[2:end,i,2] .= avg_cap

        if per_y
            @inbounds avg_border = 0.5*(a[1,i,2] + a[end,i,4])
            @inbounds a[1,i,2] = avg_border
            @inbounds a[end,i,4] = avg_border
        end
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

function get_NB_width(grid, MIXED, NB_indices_base)
    @unpack nx, ny = grid
    NB_indices = Vector{CartesianIndex{2}}(undef, 0)
    for II in MIXED
        for JJ in NB_indices_base
            KK = II + JJ
            if KK[1] > 0 && KK[2] > 0 && KK[1] < ny+1 && KK[2] < nx+1
                push!(NB_indices, KK)
            end
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

function projection_2points(grid, LS, II)
    @unpack x_nodes, y_nodes, x, y, nx, ny, dx, dy = grid
    @unpack mid_point, α = LS

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

function init_fresh_cells!(grid, T::Vector, projection, FRESH, periodic_x, periodic_y)
    _T = reshape(T, (grid.ny, grid.nx))
    @inbounds @threads for II in FRESH
        if projection[II].flag
            pII = lexicographic(II, grid.ny)
            T_1, T_2 = interpolated_temperature(grid, projection[II].angle, projection[II].point1, projection[II].point2, _T, II, periodic_x, periodic_y)
            if π/4 <= projection[II].angle <= 3π/4 || -π/4 >= projection[II].angle >= -3π/4
                T[pII] = y_extrapolation(T_1, T_2, projection[II].point1, projection[II].point2, projection[II].mid_point)
            else
                T[pII] = x_extrapolation(T_1, T_2, projection[II].point1, projection[II].point2, projection[II].mid_point)
            end
        end
    end
end

function init_fresh_cells!(grid, u::Matrix, V::Matrix, projection, FRESH, periodic_x, periodic_y)
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

function init_fresh_cells!(grid, u::Vector, V, projection, FRESH, periodic_x, periodic_y)
    _V = reshape(V, (grid.ny, grid.nx))
    @inbounds @threads for II in FRESH
        if projection[II].flag
            pII = lexicographic(II, grid.ny)
            u_1, u_2 = interpolated_temperature(grid, projection[II].angle, projection[II].point1, projection[II].point2, _V, II, periodic_x, periodic_y)
            if π/4 <= projection[II].angle <= 3π/4 || -π/4 >= projection[II].angle >= -3π/4
                u[pII] = y_extrapolation(u_1, u_2, projection[II].point1, projection[II].point2, projection[II].mid_point)
            else
                u[pII] = x_extrapolation(u_1, u_2, projection[II].point1, projection[II].point2, projection[II].mid_point)
            end
        end
    end
end

function init_fresh_cells!(grid, u::SubArray{T,N,P,I,L}, V, projection, FRESH, periodic_x, periodic_y) where {T,N,P<:Vector{T},I,L}
    _V = reshape(V, (grid.ny, grid.nx))
    @inbounds @threads for II in FRESH
        if !isassigned(projection, II)
            continue
        end
        if projection[II].flag
            pII = lexicographic(II, grid.ny)
            u_1, u_2 = interpolated_temperature(grid, projection[II].angle, projection[II].point1, projection[II].point2, _V, II, periodic_x, periodic_y)
            if π/4 <= projection[II].angle <= 3π/4 || -π/4 >= projection[II].angle >= -3π/4
                u[pII] = y_extrapolation(u_1, u_2, projection[II].point1, projection[II].point2, projection[II].mid_point)
            else
                u[pII] = x_extrapolation(u_1, u_2, projection[II].point1, projection[II].point2, projection[II].mid_point)
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
