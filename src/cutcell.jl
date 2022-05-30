@inline SOUTH_face(itp) = 0.5 + find_zero(x -> biquadratic(itp, x, -0.5), (-0.5,0.5), FalsePosition(), maxevals = 10, atol = 1e-9)
@inline WEST_face(itp) = 0.5 + find_zero(y -> biquadratic(itp, -0.5, y), (-0.5,0.5), FalsePosition(), maxevals = 10, atol = 1e-9)
@inline NORTH_face(itp) = 0.5 + find_zero(x -> biquadratic(itp, x, 0.5), (-0.5,0.5), FalsePosition(), maxevals = 10, atol = 1e-9)
@inline EAST_face(itp) = 0.5 + find_zero(y -> biquadratic(itp, 0.5, y), (-0.5,0.5), FalsePosition(), maxevals = 10, atol = 1e-9)

@inline WE(a, cases, II, l_face, s_face) = ismixed(cases[δx⁻(II)]) ? 0.5*(a[II,1] + a[δx⁻(II), 3]) : (is_liquid(cases[δx⁻(II)]) ? l_face : s_face)
@inline EW(a, cases, II, l_face, s_face) = ismixed(cases[δx⁺(II)]) ? 0.5*(a[II,3] + a[δx⁺(II), 1]) : (is_liquid(cases[δx⁺(II)]) ? l_face : s_face)
@inline SN(a, cases, II, l_face, s_face) = ismixed(cases[δy⁻(II)]) ? 0.5*(a[II,2] + a[δy⁻(II), 4]) : (is_liquid(cases[δy⁻(II)]) ? l_face : s_face)
@inline NS(a, cases, II, l_face, s_face) = ismixed(cases[δy⁺(II)]) ? 0.5*(a[II,4] + a[δy⁺(II), 2]) : (is_liquid(cases[δy⁺(II)]) ? l_face : s_face)

@inline WE(a, II) = 0.5*(a[II,1] + a[δx⁻(II), 3])
@inline EW(a, II) = 0.5*(a[II,3] + a[δx⁺(II), 1])
@inline SN(a, II) = 0.5*(a[II,2] + a[δy⁻(II), 4])
@inline NS(a, II) = 0.5*(a[II,4] + a[δy⁺(II), 2])

@inline ispositive(a) = min(max(a,0),1)
@inline isovalue(v::SVector{4,Float64}) = 4. * v[1] + 2. * v[2] + v[3] + 0.5 * v[4]

@inline vertices_sign(itp) = 1 .- sign.(@SVector [biquadratic(itp, -0.5, 0.5),
                            biquadratic(itp, 0.5, 0.5),
                            biquadratic(itp, 0.5, -0.5),
                            biquadratic(itp, -0.5, -0.5)])

@inline is_near_interface(a, st) = @inbounds ifelse(a != sign(st[2,1]) || a != sign(st[1,2]) || a != sign(st[2,3]) || a != sign(st[3,2]), true, false)

@inline is_near_interface(u, II::CartesianIndex) = @inbounds ifelse(u[II]*u[δx⁺(II)] < 0 || u[II]*u[δy⁺(II)] < 0 || u[II]*u[δx⁻(II)] < 0 || u[II]*u[δy⁻(II)] < 0, true, false)

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

function find_radius(u, MIXED, ISO, B, BT, H, inside_indices, Δ)
    radius = Vector{Float64}(undef,0)
    for II in MIXED
        st = static_stencil(u, II)
        itp = B*st*BT
        try
            faces = face_capacities(itp, ISO[II])
            a, b, c, d, e, mid_point = capacities(faces, ISO[II])
            radius_ = distance(Point(0.,0.), Point(H[II.I[2]]+Δ*mid_point.x,  H[II.I[1]]+Δ*mid_point.y))
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

function marching_squares!(H, iso, u, TS, TL, κ, SOL, LIQ, sol_projection, liq_projection, Δ, L0, B, BT, inside_indices, ϵ, nx, ny, faces)
    empty_capacities = vcat(zeros(7), zeros(4))
    full_capacities = vcat(ones(7), 0.5.*ones(4))
    κ .= zeros(ny,nx)
    @inbounds @threads for II in inside_indices
        st = static_stencil(u, II)
        a = sign(u[II])
        ISO = ifelse(a > 0, 0., 15.)
        if is_near_interface(a, st)
            κ_ = 0.
            grid_alignement = check_grid_alignement(u, II, κ_, ϵ)
            if is_aligned_with_grid(grid_alignement)
                if a > 0
                    ISO = 0.
                    @goto notmixed
                else
                    ISO = -1.
                    κ[II] = 0.
                    SOL[II,:] .= full_capacities
                    LIQ[II,:] .= empty_capacities
                    TL[II] = 0.0
                    if grid_alignement[1] == true
                        α = 0.
                        mid_point = Point(0.5, 0.)
                    elseif grid_alignement[2] == true
                        α = float(π)
                        mid_point = Point(-0.5, 0.)
                    elseif grid_alignement[3] == true
                        α = π/2
                        mid_point = Point(0., 0.5)
                    elseif grid_alignement[4] == true
                        α = -π/2
                        mid_point = Point(0., -0.5)
                    end
                end
            else
                itp = B*st*BT
                vertices = vertices_sign(itp)
                ISO = isovalue(vertices)

                if is_not_mixed(ISO) @goto notmixed end
                face_capacities(faces, itp, ISO, II)
            end
            if ISO == -1
                sol_projection[II].flag = false
                liq_projection[II].flag = false
            end
        end
        @label notmixed
        if is_solid(ISO)
            SOL[II,:] .= full_capacities
            LIQ[II,:] .= empty_capacities
            # TL[II] = 0.0
        elseif is_liquid(ISO)
            SOL[II,:] .= empty_capacities
            LIQ[II,:] .= full_capacities
            # TS[II] = 0.0
        end
        iso[II] = ISO
    end
    return nothing
end

function get_iterface_location!(::GridCC, x, y, iso, u, TS, TL, κ, SOL, LIQ, sol_projection, liq_projection, sol_centroid, liq_centroid, mid_point, α, Δ, L0, B, BT, idx, inside_indices, ϵ, n, faces, periodic_x, periodic_y)
    @inbounds @threads for II in inside_indices
        f = average_face_capacities(faces, iso[II], iso, II)
        SOL[II,:], LIQ[II,:], α[II], sol_centroid[II], liq_centroid[II], mid_point[II], cut_points = capacities(f, iso[II])
        absolute_position = Point(x[II], y[II])
        sol_projection[II], liq_projection[II] = projection_2points(absolute_position, mid_point[II], α[II], L0, Δ)
    end
    set_cap_bcs!(gcc, SOL, LIQ, sol_centroid, liq_centroid, mid_point, u, idx.b_left[1], idx.b_bottom[1], idx.b_right[1], idx.b_top[1], periodic_x, periodic_y, n)
    Wcapacities!(SOL, periodic_x, periodic_y)
    Wcapacities!(LIQ, periodic_x, periodic_y)
    average_face_capacities!(SOL)
    average_face_capacities!(LIQ)
    return nothing
end

function get_iterface_location!(::GridFCx, x, y, iso, u, TS, TL, κ, SOL, LIQ, sol_projection, liq_projection, sol_centroid, liq_centroid, mid_point, α, Δ, L0, B, BT, idx, inside_indices, ϵ, n, faces, periodic_x, periodic_y)
    @inbounds @threads for II in inside_indices
        f = average_face_capacities(faces, iso[II], iso, II)
        SOL[II,:], LIQ[II,:], α[II], sol_centroid[II], liq_centroid[II], mid_point[II], cut_points = capacities(f, iso[II])
        absolute_position = Point(x[II], y[II])
        sol_projection[II], liq_projection[II] = projection_2points(absolute_position, mid_point[II], α[II], L0, Δ)
    end
    set_cap_bcs!(gfcx, SOL, LIQ, sol_centroid, liq_centroid, mid_point, u, idx.b_left[1], idx.b_bottom[1], idx.b_right[1], idx.b_top[1], periodic_x, periodic_y, n)
    Wcapacities!(SOL, periodic_x, periodic_y)
    Wcapacities!(LIQ, periodic_x, periodic_y)
    average_face_capacities!(SOL)
    average_face_capacities!(LIQ)
    return nothing
end

function get_iterface_location!(::GridFCy, x, y, iso, u, TS, TL, κ, SOL, LIQ, sol_projection, liq_projection, sol_centroid, liq_centroid, mid_point, α, Δ, L0, B, BT, idx, inside_indices, ϵ, n, faces, periodic_x, periodic_y)
    @inbounds @threads for II in inside_indices
        f = average_face_capacities(faces, iso[II], iso, II)
        SOL[II,:], LIQ[II,:], α[II], sol_centroid[II], liq_centroid[II], mid_point[II], cut_points = capacities(f, iso[II])
        absolute_position = Point(x[II], y[II])
        sol_projection[II], liq_projection[II] = projection_2points(absolute_position, mid_point[II], α[II], L0, Δ)
    end
    set_cap_bcs!(gfcy, SOL, LIQ, sol_centroid, liq_centroid, mid_point, u, idx.b_left[1], idx.b_bottom[1], idx.b_right[1], idx.b_top[1], periodic_x, periodic_y, n)
    Wcapacities!(SOL, periodic_x, periodic_y)
    Wcapacities!(LIQ, periodic_x, periodic_y)
    average_face_capacities!(SOL)
    average_face_capacities!(LIQ)
    return nothing
end

function get_curvature(u, POS, κ, B, BT, h, inside_indices)
    @inbounds for II in inside_indices
        mid_point = POS[II].mid_point
        st = static_stencil(u, II)
        itp = B*st*BT
        κ[II] = parabola_fit_curvature(itp, mid_point, h)
    end
end

function capacities(F, case)
    cap_sol = @SVector zeros(11)
    cap_liq = @SVector zeros(11)
    α = 0.0
    sol_centroid = Point(NaN, NaN)
    liq_centroid = Point(NaN, NaN)
    mid_point = Point(NaN, NaN)
    cut_points = (Point(NaN, NaN), Point(NaN, NaN))
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

function set_cap_bcs!(::GridCC, SOL, LIQ, sol_centroid, liq_centroid, mid_point, u, b_left, b_bottom, b_right, b_top, periodic_x, periodic_y, n)
    empty_capacities = SA_F64[vcat(zeros(7), zeros(4))...]
    full_capacities = SA_F64[vcat(ones(7), 0.5.*ones(4))...]

    # fill boundary capacities
    @inbounds @threads for II in b_left
        if u[II] >= 0.0
            @inbounds SOL[II,:] .= empty_capacities
            @inbounds LIQ[II,:] .= full_capacities
        else
            @inbounds SOL[II,:] .= full_capacities
            @inbounds LIQ[II,:] .= empty_capacities
        end
    end
    @inbounds @threads for II in b_bottom
        if u[II] >= 0.0
            @inbounds SOL[II,:] .= empty_capacities
            @inbounds LIQ[II,:] .= full_capacities
        else
            @inbounds SOL[II,:] .= full_capacities
            @inbounds LIQ[II,:] .= empty_capacities
        end
    end
    @inbounds @threads for II in b_right
        if u[II] >= 0.0
            @inbounds SOL[II,:] .= empty_capacities
            @inbounds LIQ[II,:] .= full_capacities
        else
            @inbounds SOL[II,:] .= full_capacities
            @inbounds LIQ[II,:] .= empty_capacities
        end
    end
    @inbounds @threads for II in b_top
        if u[II] >= 0.0
            @inbounds SOL[II,:] .= empty_capacities
            @inbounds LIQ[II,:] .= full_capacities
        else
            @inbounds SOL[II,:] .= full_capacities
            @inbounds LIQ[II,:] .= empty_capacities
        end
    end

    # set A and mid_point at the boundaries to 0 if not periodic in that direction
    if !periodic_x
        @inbounds @threads for i = 1:n
            @inbounds SOL[i,1,1] = 0.0
            @inbounds SOL[i,end,3] = 0.0
            @inbounds LIQ[i,1,1] = 0.0
            @inbounds LIQ[i,end,3] = 0.0
            @inbounds mid_point[i,1] += Point(-0.5, 0.0)
            @inbounds mid_point[i,end] += Point(0.5, 0.0)
        end
    end
    if !periodic_y
        @inbounds @threads for i = 1:n
            @inbounds SOL[1,i,2] = 0.0
            @inbounds SOL[end,i,4] = 0.0
            @inbounds LIQ[1,i,2] = 0.0
            @inbounds LIQ[end,i,4] = 0.0
            @inbounds mid_point[1,i] += Point(0.0, -0.5)
            @inbounds mid_point[end,i] += Point(0.0, 0.5)
        end
    end
end

function set_cap_bcs!(::GridFCx, SOL, LIQ, sol_centroid, liq_centroid, mid_point, u, b_left, b_bottom, b_right, b_top, periodic_x, periodic_y, n)
    f = SA_F64[0.5, 0.5]
    empty_capacities = SA_F64[vcat(zeros(7), zeros(4))...]
    full_capacities = SA_F64[vcat(ones(7), 0.5.*ones(4))...]
    capacities_6 = capacities(f, 6.0)
    capacities_9 = capacities(f, 9.0)

    # fill boundary capacities
    if periodic_x
        @inbounds @threads for II in b_left
            if u[II] >= 0.0
                @inbounds SOL[II,:] .= empty_capacities
                @inbounds LIQ[II,:] .= full_capacities
            else
                @inbounds SOL[II,:] .= full_capacities
                @inbounds LIQ[II,:] .= empty_capacities
            end
        end
        @inbounds @threads for II in b_right
            if u[II] >= 0.0
                @inbounds SOL[II,:] .= empty_capacities
                @inbounds LIQ[II,:] .= full_capacities
            else
                @inbounds SOL[II,:] .= full_capacities
                @inbounds LIQ[II,:] .= empty_capacities
            end
        end
    else
        @inbounds @threads for II in b_left
            if u[II] >= 0.0
                @inbounds SOL[II,:] .= empty_capacities
                @inbounds _, LIQ[II,:], _, _, liq_centroid[II], mid_point[II], _ = capacities_9
            else
                @inbounds SOL[II,:], _, _, sol_centroid[II], _, mid_point[II], _ = capacities_6
                @inbounds LIQ[II,:] .= empty_capacities
            end
        end
        @inbounds @threads for II in b_right
            if u[II] >= 0.0
                @inbounds SOL[II,:] .= empty_capacities
                @inbounds _, LIQ[II,:], _, _, liq_centroid[II], mid_point[II], _ = capacities_6
            else
                @inbounds SOL[II,:], _, _, sol_centroid[II], _, mid_point[II], _ = capacities_9
                @inbounds LIQ[II,:] .= empty_capacities
            end
        end
    end
    @inbounds @threads for II in b_bottom
        if II ∉ b_left && II ∉ b_right && u[II] >= 0.0
            @inbounds SOL[II,:] .= empty_capacities
            @inbounds LIQ[II,:] .= full_capacities
        elseif II ∉ b_left && II ∉ b_right
            @inbounds SOL[II,:] .= full_capacities
            @inbounds LIQ[II,:] .= empty_capacities
        end
    end
    @inbounds @threads for II in b_top
        if II ∉ b_left && II ∉ b_right && u[II] >= 0.0
            @inbounds SOL[II,:] .= empty_capacities
            @inbounds LIQ[II,:] .= full_capacities
        elseif II ∉ b_left && II ∉ b_right
            @inbounds SOL[II,:] .= full_capacities
            @inbounds LIQ[II,:] .= empty_capacities
        end
    end

    # set A at the boundaries
    if !periodic_y
        @inbounds @threads for i = 1:n+1
            @inbounds SOL[1,i,2] = 0.0
            @inbounds SOL[end,i,4] = 0.0
            @inbounds LIQ[1,i,2] = 0.0
            @inbounds LIQ[end,i,4] = 0.0
        end
    end

    # set mid_point in the outer boundaries
    if !periodic_y
        @inbounds @threads for i = 2:n
            @inbounds mid_point[1,i] += Point(0.0, -0.5)
            @inbounds mid_point[end,i] += Point(0.0, 0.5)
        end
    end
end

function set_cap_bcs!(::GridFCy, SOL, LIQ, sol_centroid, liq_centroid, mid_point, u, b_left, b_bottom, b_right, b_top, periodic_x, periodic_y, n)
    f = SA_F64[0.5, 0.5]
    empty_capacities = SA_F64[vcat(zeros(7), zeros(4))...]
    full_capacities = SA_F64[vcat(ones(7), 0.5.*ones(4))...]
    capacities_3 = capacities(f, 3.0)
    capacities_12 = capacities(f, 12.0)

    # fill boundary capacities
    @inbounds @threads for II in b_left
        if II ∉ b_bottom && II ∉ b_top && u[II] >= 0.0
            @inbounds SOL[II,:] .= empty_capacities
            @inbounds LIQ[II,:] .= full_capacities
        elseif II ∉ b_bottom && II ∉ b_top
            @inbounds SOL[II,:] .= full_capacities
            @inbounds LIQ[II,:] .= empty_capacities
        end
    end
    @inbounds @threads for II in b_right
        if II ∉ b_bottom && II ∉ b_top && u[II] >= 0.0
            @inbounds SOL[II,:] .= empty_capacities
            @inbounds LIQ[II,:] .= full_capacities
        elseif II ∉ b_bottom && II ∉ b_top
            @inbounds SOL[II,:] .= full_capacities
            @inbounds LIQ[II,:] .= empty_capacities
        end
    end
    if periodic_y
        @inbounds @threads for II in b_bottom
            if u[II] >= 0.0
                @inbounds SOL[II,:] .= empty_capacities
                @inbounds LIQ[II,:] .= full_capacities
            else
                @inbounds SOL[II,:] .= full_capacities
                @inbounds LIQ[II,:] .= empty_capacities
            end
        end
        @inbounds @threads for II in b_top
            if u[II] >= 0.0
                @inbounds SOL[II,:] .= empty_capacities
                @inbounds LIQ[II,:] .= full_capacities
            else
                @inbounds SOL[II,:], _, _, sol_centroid[II], _, mid_point[II], _ = capacities_3
                @inbounds LIQ[II,:] .= empty_capacities
            end
        end
    else
        @inbounds @threads for II in b_bottom
            if u[II] >= 0.0
                @inbounds SOL[II,:] .= empty_capacities
                @inbounds _, LIQ[II,:], _, _, liq_centroid[II], mid_point[II], _ = capacities_3
            else
                @inbounds SOL[II,:], _, _, sol_centroid[II], _, mid_point[II], _ = capacities_12
                @inbounds LIQ[II,:] .= empty_capacities
            end
        end
        @inbounds @threads for II in b_top
            if u[II] >= 0.0
                @inbounds SOL[II,:] .= empty_capacities
                @inbounds _, LIQ[II,:], _, _, liq_centroid[II], mid_point[II], _ = capacities_12
            else
                @inbounds SOL[II,:], _, _, sol_centroid[II], _, mid_point[II], _ = capacities_3
                @inbounds LIQ[II,:] .= empty_capacities
            end
        end
    end

    # set A at the boundaries to 0 if not periodic in that direction
    if !periodic_x
        @inbounds @threads for i = 1:n+1
            @inbounds SOL[i,1,1] = 0.0
            @inbounds SOL[i,end,3] = 0.0
            @inbounds LIQ[i,1,1] = 0.0
            @inbounds LIQ[i,end,3] = 0.0
        end
    end

    # set mid_point in the outer boundaries
    if !periodic_y
        @inbounds @threads for i = 2:n
            @inbounds mid_point[1,i] += Point(0.0, -0.5)
            @inbounds mid_point[end,i] += Point(0.0, 0.5)
        end
    end
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

    @inbounds cap[:,2:end,8] .= tmp[:,2:end,8] .+ tmp[:,1:end-1,10]
    @inbounds cap[:,1:end-1,10] .= tmp[:,1:end-1,10] .+ tmp[:,2:end,8]
    @inbounds cap[2:end,:,9] .= tmp[2:end,:,9] .+ tmp[1:end-1,:,11]
    @inbounds cap[1:end-1,:,11] .= tmp[1:end-1,:,11] .+ tmp[2:end,:,9]

    if periodic_x
        @inbounds cap[:,1,8] .= tmp[:,1,8] .+ tmp[:,end,10]
        @inbounds cap[:,end,10] .= tmp[:,end,10] .+ tmp[:,1,8]
    end
    if periodic_y
        @inbounds cap[1,:,9] .= tmp[1,:,9] .+ tmp[end,:,11]
        @inbounds cap[end,:,11] .= tmp[end,:,11] .+ tmp[1,:,9]
    end

    return nothing
end

function face_capacities(itp, case)
    if case == 1.0 || case == 14.0
        f = SA_F64[ispositive(WEST_face(itp)), ispositive(SOUTH_face(itp))]
    elseif case == 2.0 || case == 13.0
        f = SA_F64[ispositive(SOUTH_face(itp)), ispositive(EAST_face(itp))]
    elseif case == 3.0 || case == 12.0
        f = SA_F64[ispositive(WEST_face(itp)), ispositive(EAST_face(itp))]
    elseif case == 4.0 || case == 11.0
        f = SA_F64[ispositive(EAST_face(itp)), ispositive(NORTH_face(itp))]
    elseif case == 6.0 || case == 9.0
        f = SA_F64[ispositive(SOUTH_face(itp)), ispositive(NORTH_face(itp))]
    elseif case == 7.0 || case == 8.0
        f = SA_F64[ispositive(WEST_face(itp)), ispositive(NORTH_face(itp))]
    else
        f = @SVector zeros(2)
    end
    return float(f)
end

function face_capacities(a, itp, case, II)
    if case == 1.0 || case == 14.0
        a[II, 1] = ispositive(WEST_face(itp))
        a[II, 2] = ispositive(SOUTH_face(itp))
    elseif case == 2.0 || case == 13.0
        a[II, 2] = ispositive(SOUTH_face(itp))
        a[II, 3] = ispositive(EAST_face(itp))
    elseif case == 3.0 || case == 12.0
        a[II, 1] = ispositive(WEST_face(itp))
        a[II, 3] = ispositive(EAST_face(itp))
    elseif case == 4.0 || case == 11.0
        a[II, 3] = ispositive(EAST_face(itp))
        a[II, 4] = ispositive(NORTH_face(itp))
    elseif case == 6.0 || case == 9.0
        a[II, 2] = ispositive(SOUTH_face(itp))
        a[II, 4] = ispositive(NORTH_face(itp))
    elseif case == 7.0 || case == 8.0
        a[II, 1] = ispositive(WEST_face(itp))
        a[II, 4] = ispositive(NORTH_face(itp))
    end
end

function average_face_capacities(a, case, cases, II)
    if case == 1.0
        f = SA_F64[WE(a, cases, II, 0.0, 1.0), SN(a, cases, II, 0.0, 1.0)]
    elseif case == 14.0
        f = SA_F64[WE(a, cases, II, 1.0, 0.0), SN(a, cases, II, 1.0, 0.0)]
    elseif case == 2.0
        f = SA_F64[SN(a, cases, II, 1.0, 0.0), EW(a, cases,II, 0.0, 1.0)]
    elseif case == 13.0
        f = SA_F64[SN(a, cases, II, 0.0, 1.0), EW(a, cases,II, 1.0, 0.0)]
    elseif case == 3.0
        f = SA_F64[WE(a, cases, II, 0.0, 1.0), EW(a, cases,II, 0.0, 1.0)]
    elseif case == 12.0
        f = SA_F64[WE(a, cases, II, 1.0, 0.0), EW(a, cases,II, 1.0, 0.0)]
    elseif case == 4.0
        f = SA_F64[EW(a, cases, II, 1.0, 0.0), NS(a, cases, II, 1.0, 0.0)]
    elseif case == 11.0
        f = SA_F64[EW(a, cases, II, 0.0, 1.0), NS(a, cases, II, 0.0, 1.0)]
    elseif case == 6.0
        f = SA_F64[SN(a, cases, II, 1.0, 0.0), NS(a, cases, II, 1.0, 0.0)]
    elseif case == 9.0
        f = SA_F64[SN(a, cases, II, 0.0, 1.0), NS(a, cases, II, 0.0, 1.0)]
    elseif case == 7.0
        f = SA_F64[WE(a, cases, II, 0.0, 1.0), NS(a, cases, II, 1.0, 0.0)]
    elseif case == 8.0
        f = SA_F64[WE(a, cases, II, 1.0, 0.0), NS(a, cases, II, 0.0, 1.0)]
    else
        f = @SVector zeros(2)
    end
    return float(f)
end

function average_face_capacities!(a)
    # average face capacities due to mismatch in iso 
    n = size(a,1)
    @inbounds @threads for i = 1:n
        @inbounds tmp1 = a[i,1:end-1,3]
        @inbounds tmp2 = a[i,2:end,1]
        @inbounds avg_cap = 0.5*(tmp1 .+ tmp2)
        @inbounds a[i,1:end-1,3] .= avg_cap
        @inbounds a[i,2:end,1] .= avg_cap
    end
    n = size(a,2)
    @inbounds @threads for i = 1:n
        @inbounds tmp1 = a[1:end-1,i,4]
        @inbounds tmp2 = a[2:end,i,2]
        @inbounds avg_cap = 0.5*(tmp1 .+ tmp2)
        @inbounds a[1:end-1,i,4] .= avg_cap
        @inbounds a[2:end,i,2] .= avg_cap
    end
end

# function average_face_capacities(a, case, cases, II)
#     if case == 1.0 || case == 14.0
#         f = SA_F64[WE(a, II), SN(a, II)]
#     elseif case == 2.0 || case == 13.0
#         f = SA_F64[SN(a, II), EW(a,II)]
#     elseif case == 3.0 || case == 12.0
#         f = SA_F64[WE(a, II), EW(a,II)]
#     elseif case == 4.0 || case == 11.0
#         f = SA_F64[EW(a,II), NS(a, II)]
#     elseif case == 6.0 || case == 9.0
#         f = SA_F64[SN(a, II), NS(a, II)]
#     elseif case == 7.0 || case == 8.0
#         f = SA_F64[WE(a, II), NS(a, II)]
#     else
#         f = @SVector zeros(2)
#     end
#     return float(f)
# end


function get_cells_indices(iso, inside)
    local M = 1
    local S = 1
    local L = 1
    MIXED = Vector{CartesianIndex{2}}(undef, length(inside))
    SOLID = Vector{CartesianIndex{2}}(undef, length(inside))
    LIQUID = Vector{CartesianIndex{2}}(undef, length(inside))
    @simd for II in inside
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

function get_cells_indices(iso, u, inside, boundaries)
    local M = 1
    local S = 1
    local L = 1
    MIXED = Vector{CartesianIndex{2}}(undef, length(inside)+length(boundaries))
    SOLID = Vector{CartesianIndex{2}}(undef, length(inside)+length(boundaries))
    LIQUID = Vector{CartesianIndex{2}}(undef, length(inside)+length(boundaries))
    @simd for II in inside
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
    @simd for II in boundaries
        @inbounds if u[II] <= 0.
            SOLID[S] = II
            S+=1
        else
            LIQUID[L] = II
            L+=1
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

function projection_2points(absolute_position, mid, α, L0, h)
    β = opposite(α)
    tmpS = Point(mid.x + cos(α), mid.y + sin(α))
    tmpL = Point(mid.x + cos(β), mid.y + sin(β))
    S2 = Point(NaN, NaN)
    L2 = Point(NaN, NaN)
    Sflag = false
    Lflag = false
    if -π/4 < α <= π/4
        S1 = affine(mid, tmpS, Point(-1., NaN))
        S2 = affine(mid, tmpS, Point(-2., NaN))
        L1 = affine(mid, tmpL, Point(1., NaN))
        L2 = affine(mid, tmpL, Point(2., NaN))
    elseif π/4 < α <= 3π/4
        S1 = affine(mid, tmpS, Point(NaN, -1.))
        S2 = affine(mid, tmpS, Point(NaN, -2.))
        L1 = affine(mid, tmpL, Point(NaN, 1.))
        L2 = affine(mid, tmpL, Point(NaN, 2.))
    elseif α > 3π/4 || α <= -3π/4
        S1 = affine(mid, tmpS, Point(1., NaN))
        S2 = affine(mid, tmpS, Point(2., NaN))
        L1 = affine(mid, tmpL, Point(-1., NaN))
        L2 = affine(mid, tmpL, Point(-2., NaN))
    elseif -3π/4 < α <= -π/4
        S1 = affine(mid, tmpS, Point(NaN, 1.))
        S2 = affine(mid, tmpS, Point(NaN, 2.))
        L1 = affine(mid, tmpL, Point(NaN, -1.))
        L2 = affine(mid, tmpL, Point(NaN, -2.))
    end
    if S2 != Point(NaN, NaN) && L2 != Point(NaN, NaN)
        Sflag = absolute_position + h*S2 < L0
        Lflag = absolute_position + h*L2 < L0
    end
    pos = absolute_position + h*mid
    return Gradient(Sflag, β, mid, S1, S2, distance(mid, S1), distance(mid, S2), pos), Gradient(Lflag, α, mid, L1, L2, distance(mid, L1), distance(mid, L2), pos)
end

function kill_dead_cells!(T, EMPTY)
    @inbounds @threads for II in EMPTY
        T[II] = 0.
    end
end

function init_fresh_cells!(T, projection, FRESH)
    @inbounds @threads for II in FRESH
        if projection[II].flag
            T_1, T_2 = interpolated_temperature(projection[II].angle, projection[II].point1, projection[II].point2, T, II)
            if π/4 <= projection[II].angle <= 3π/4 || -π/4 >= projection[II].angle >= -3π/4
                T[II] = y_extrapolation(T_1, T_2, projection[II].point1, projection[II].point2, projection[II].mid_point)
            else
                T[II] = x_extrapolation(T_1, T_2, projection[II].point1, projection[II].point2, projection[II].mid_point)
            end
        end
    end
end

function init_fresh_cells!(u, V, projection, FRESH)
    @inbounds @threads for II in FRESH
        if projection[II].flag
            u_1, u_2 = interpolated_temperature(projection[II].angle, projection[II].point1, projection[II].point2, V, II)
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
