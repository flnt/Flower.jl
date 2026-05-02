# front_tracking.jl
# ==================
# Front-Tracking for Flower.jl using the interface geometry already computed
# by update_all_ls_data (projections + connectivity from convert_interfacial_D_to_segments).
#
# Three marker velocity modes:
#   :exact     →  V_n = num.mdot_rho  (constant — shrinking bubble test)
#   :bilinear    →  V_n bilinearly interpolated from grid_p.V
#   :peskin    →  V_n via Peskin regularised δ kernel from grid_p.V
#
# Two LS options (recompute_ls kwarg):
#   true   →  φ rebuilt from scratch (exact signed distance, no drift)
#   false  →  grid_p.LS[1].u untouched (test pure marker motion)
#
# Integration in run.jl — three additions:
#   1. Before time loop, after update_all_ls_data:
#         ft = FrontTracker(num, grid_p)
#   2. Inside time loop, advection block:
#         elseif num.advection_LS_mode_symb === :front_tracking
#             ft_step!(ft, num, grid_p, grid_u, grid_v;
#                      vel_mode=:exact, recompute_ls=true)
#             ft_pdi!(ft, num, grid_p, grid_u, grid_v)
#   3. update_all_ls_data called right after (unchanged).

# ─────────────────────────────────────────────────────────────────────────────
# Data structure
# ─────────────────────────────────────────────────────────────────────────────

mutable struct FrontTracker
    x              :: Vector{Float64}   # marker x  (n_markers,)
    y              :: Vector{Float64}   # marker y
    normal         :: Matrix{Float64}   # marker y
    connectivities :: Matrix{Int64}     # (n_seg, 2)  vertex index pairs
    n_markers      :: Int64
    n_seg          :: Int64
    vn_markers     :: Vector{Float64}   # workspace: normal speed at each marker
    # Grid metadata
    nx :: Int;  ny :: Int;  Δ :: Float64
    xc1 :: Float64;  yc1 :: Float64    # centre of cell (j=1, i=1)
end

"""
    FrontTracker(num, grid_p)

Initialise from the interface geometry already in grid_p.LS[1].
Uses convert_interfacial_D_to_segments for vertex positions and connectivity,
exactly matching what Flower writes to PDI.
"""
function FrontTracker(num, grid_p)


    init = :circle

    if init === :circle
        radius=num.current_radius; center=(num.intfc_x, num.intfc_y)
        interface_length = 2*pi*radius

        num_vertices=round(Int64,interface_length / grid_p.dx[2,2])
        n_m = num_vertices
        n_s = n_m
        print("\n num_vertices ",num_vertices)
        # Generate vertices for a circle
        intfc_vtx_x, intfc_vtx_y = generate_circle_vertices(radius, center, num_vertices)

        # Initialize connectivities for a closed circle
        # connectivities_flat, num_segments = generate_circle_connectivities(num_vertices)
        # connectivities = reshape(connectivities_flat, (n_s, 2))

        connectivities, num_segments = generate_circle_connectivities(num_vertices)

        # print("\n connectivities_flat ",connectivities_flat)
        # Compute normals for each vertex (radial direction)
        normals = compute_circle_normals(intfc_vtx_x, intfc_vtx_y, center)
        

        # print("\n connectivities ",connectivities)

    else 
        intfc_vtx_x, intfc_vtx_y, intfc_vtx_connectivities,
            intfc_vtx_num, intfc_seg_num =
                convert_interfacial_D_to_segments(num, grid_p,1)

        n_m = Int64(intfc_vtx_num)
        n_s = Int64(intfc_seg_num)

        print("\n FrontTracker")
        print("\n intfc_vtx_connectivities",intfc_vtx_connectivities)

        # intfc_vtx_connectivities is a flat vector of length 2*n_s (v1, v2, v1, v2, ...)
        # reshape into (n_s, 2)
        connectivities = reshape(copy(intfc_vtx_connectivities[1:2*n_s]), 2, n_s)'   # (n_s, 2)

        print("\n conn",conn)

    end


    Δ   = num.Δ
    nx  = grid_p.nx
    ny  = grid_p.ny
    xc1 = grid_p.x[1, 1]
    yc1 = grid_p.y[1, 1]

    return FrontTracker(
        copy(intfc_vtx_x[1:n_m]),
        copy(intfc_vtx_y[1:n_m]),
        normals,
        connectivities, n_m, n_s,
        zeros(n_m),
        nx, ny, Δ, xc1, yc1
    )
end



# """
#     FrontTracker(num, grid_p; radius=1.0, center=(0.0, 0.0), num_vertices=32)

# Initialize a FrontTracker for a circle with the given radius, center, and number of vertices.
# """
# function FrontTracker(num, grid_p; radius=1.0, center=(0.0, 0.0), num_vertices=32)
#     # Generate vertices for a circle
#     x, y = generate_circle_vertices(radius, center, num_vertices)

#     # Initialize connectivities for a closed circle
#     connectivities, num_segments = generate_circle_connectivities(num_vertices)

#     # Compute normals for each vertex (radial direction)
#     normals = compute_circle_normals(x, y, center)

#     # Extract grid parameters
#     Δ   = num.Δ
#     nx  = grid_p.nx
#     ny  = grid_p.ny
#     xc1 = grid_p.x[1, 1]
#     yc1 = grid_p.y[1, 1]

#     return FrontTracker(
#         x,
#         y,
#         connectivities,
#         num_vertices,
#         num_segments,
#         normals,
#         nx,
#         ny,
#         Δ,
#         xc1,
#         yc1
#     )
# end

"""
    generate_circle_vertices(radius, center, num_vertices)

Generate vertices for a circle with the given radius, center, and number of vertices.
"""
function generate_circle_vertices(radius, center, num_vertices)
    x = zeros(num_vertices)
    y = zeros(num_vertices)
    cx, cy = center

    for i in 1:num_vertices
        θ = 2π * (i - 1) / num_vertices  # Angle for vertex i
        x[i] = cx + radius * cos(θ)
        y[i] = cy + radius * sin(θ)
    end

    return x, y
end

"""
    generate_circle_connectivities(num_vertices)

Generate connectivities for a closed circle with the given number of vertices.
"""
function generate_circle_connectivities(num_vertices)
    connectivities = Matrix{Int64}(undef, num_vertices, 2)

    for i in 1:num_vertices
        j = i % num_vertices + 1
        connectivities[i, 1] = i
        connectivities[i, 2] = j
    end

    return connectivities, num_vertices
end


"""
    compute_circle_normals(x, y, center)

Compute normals for each vertex of a circle (radial direction).
"""
function compute_circle_normals(x, y, center)
    cx, cy = center
    num_vertices = length(x)
    normals = zeros(num_vertices, 2)

    for i in 1:num_vertices
        # Vector from center to vertex
        dx = -x[i] + cx
        dy = -y[i] + cy

        # Normalize to get the normal (radial direction)
        norm = sqrt(dx^2 + dy^2)
        if norm > 0
            normals[i, 1] = dx / norm
            normals[i, 2] = dy / norm
        else
            normals[i, 1] = 0.0
            normals[i, 2] = 0.0
        end
    end

    return normals
end


"""
    returns x,y,f,connectivities,vtx_index,num_seg which can be used to plot interfacial value
    supports on interface 
    field_index is the index of the interfacial value: ex 2 for LS 1 since field_index = 1 is the index of the bulk value
"""
function convert_interfacial_D_to_segments(num,gp,iLS)

    # x = zeros(0)
    # y = zeros(0)
    # f = zeros(0)
    x = Float64[]
    y = Float64[]
    # f = Float64[]

    connectivities=Int64[]
    ij_index=-ones(Int64,(gp.ny, gp.nx))
    vtx_index = 0 #Integer(0) 

 


    for II in gp.LS[iLS].MIXED
        push!( x, gp.LS[1].mid_point[II].x * gp.dx[II] + gp.x[II] )
        push!( y, gp.LS[1].mid_point[II].y * gp.dy[II] + gp.y[II] )

        # push!( f, field[II] )

        # a[g.ny*g.nx*(p-1)+1:g.ny*g.nx*p]
        
        pII = lexicographic(II, gp.ny)

        # push!( f, field_D[gp.ny*gp.nx*(field_index-1) + pII] )



        ij_index[II] = vtx_index

        # print("\n vtx ", vtx_index," ", II, " ",gp.x[II]-gp.dx[II]/2," ",gp.x[II]+gp.dx[II]/2," ",x," ",gp.y[II]-gp.dy[II]/2,
        # " ",gp.y[II]+gp.dy[II]/2,y," ",gp.x[II]," ",gp.y[II])


        #TODO store connvectivities or similar so that after: only need to filter and store non-null values ? 
        #no need to do the work n times for each scalar

        vtx_index +=1
    end

    num_seg = 0
    for j in 1:gp.ny-1
        for i in 1:gp.nx-1
            # II = CartesianIndex(i,j)
            if ij_index[j,i] >=0

                if ij_index[j+1,i]>=0  
                    push!(connectivities, ij_index[j,i])
                    push!(connectivities,ij_index[j+1,i])
                    num_seg+=1
                    # print("\n seg ", i," ",j," ",ij_index[j,i], " ", ij_index[j+1,i])
                end

                if ij_index[j,i+1]>=0  
                    push!(connectivities, ij_index[j,i])
                    push!(connectivities,ij_index[j,i+1])
                    num_seg+=1
                    # print("\n seg ", i," ",j," ",ij_index[j,i], " ", ij_index[j,i+1])
                end

            end

        end
    end

 
    return x,y,connectivities,vtx_index,num_seg

end

# ─────────────────────────────────────────────────────────────────────────────
# Outward unit normal at each marker via central difference along arc
# ─────────────────────────────────────────────────────────────────────────────
function _marker_normals(x, y, n)
    nx_ = Vector{Float64}(undef, n)
    ny_ = Vector{Float64}(undef, n)
    @inbounds for k in 1:n
        kp = mod1(k+1, n);  km = mod1(k-1, n)
        tx  = x[kp] - x[km]
        ty  = y[kp] - y[km]
        mag = sqrt(tx^2 + ty^2) + 1e-14
        # nx_[k] =  ty / mag   # outward: rotate CCW
        # ny_[k] = -tx / mag
        nx_[k] = - ty / mag   
        ny_[k] = + tx / mag
    end
    return nx_, ny_
end

# ─────────────────────────────────────────────────────────────────────────────
# Bilinear interpolation of a (ny×nx) cell-centred field at (px, py)
# ─────────────────────────────────────────────────────────────────────────────
@inline function _bilinear(field, px, py, Δ, xc1, yc1, nx, ny)
    fi = (px - xc1) / Δ + 1.0
    fj = (py - yc1) / Δ + 1.0
    i0 = clamp(floor(Int, fi), 1, nx-1)
    j0 = clamp(floor(Int, fj), 1, ny-1)
    tx = clamp(fi - i0, 0.0, 1.0)
    ty = clamp(fj - j0, 0.0, 1.0)
    return ((1-tx)*(1-ty)*field[j0,   i0  ] +
                tx*(1-ty)*field[j0,   i0+1] +
            (1-tx)*   ty *field[j0+1, i0  ] +
                tx*   ty *field[j0+1, i0+1])
end

@inline function interp_u(grid_u, px, py)
    @unpack V, x, y, nx, ny = grid_u
    # u is shifted in x by -Δ/2
    dx = grid_u.dx[2,2]
    return _bilinear(V, px - 0.5*dx, py, dx, x, y, nx, ny)
end

@inline function interp_v(grid_v, px, py)
    @unpack V, x, y, nx, ny = grid_v
    # v is shifted in y by -Δ/2
    dy = grid_v.dy[2,2]
    return _bilinear(V, px, py - 0.5*dy, dy, x, y, nx, ny)
end

function compute_normal_velocity!(ft, grid_u, grid_v)
    n_m = ft.n_markers

    for k in 1:n_m
        px = ft.x[k]
        py = ft.y[k]

        u = interp_u(grid_u, px, py)
        v = interp_v(grid_v, px, py)

        nx = ft.normal[k,1]
        ny = ft.normal[k,2]

        ft.un[k] = u*nx + v*ny
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Peskin 4-point δ kernel (Peskin 2002, eq. 2.3)
# ─────────────────────────────────────────────────────────────────────────────
@inline function _peskin_φ(r)
    ra = abs(r)
    ra <= 1.0 && return 0.125*(3.0 - 2ra + sqrt(1.0 + 4ra - 4ra^2))
    ra <= 2.0 && return 0.125*(5.0 - 2ra - sqrt(-7.0 + 12ra - 4ra^2))
    return 0.0
end

function _peskin_interp(field, px, py, Δ, xc1, yc1, nx, ny)
    fi = (px - xc1) / Δ + 1.0
    fj = (py - yc1) / Δ + 1.0
    ic = round(Int, fi);  jc = round(Int, fj)
    val = 0.0
    @inbounds for dj in -1:2, di in -1:2
        i = clamp(ic+di, 1, nx);  j = clamp(jc+dj, 1, ny)
        rx = (px - (xc1 + (i-1)*Δ)) / Δ
        ry = (py - (yc1 + (j-1)*Δ)) / Δ
        val += _peskin_φ(rx) * _peskin_φ(ry) * field[j, i]
    end
    return val
end

# ─────────────────────────────────────────────────────────────────────────────
# Compute scalar normal velocity at each marker
# ─────────────────────────────────────────────────────────────────────────────
function _marker_velocities!(vn, x, y, n, grid_p, grid_u,grid_v,num, vel_mode::Symbol,ft)
    Δ = num.Δ;  nx = grid_p.nx;  ny = grid_p.ny
    xc1 = grid_p.x[1,1];  yc1 = grid_p.y[1,1]

    if vel_mode === :exact
        # Constant V_n = ṁ/ρ  (add num.mdot_rho to your Numerical struct,
        # or replace with the appropriate field name you already have)
        Vn = - num.mass_transfer_rate_imposed_value / num.rho1
        fill!(vn, Vn)

    elseif vel_mode === :bilinear
        @inbounds for k in 1:n
            # vn[k] = _bilinear(grid_p.V, x[k], y[k], Δ, xc1, yc1, nx, ny)
            
            px = x[k]
            py = y[k]

            u = interp_u(grid_u, px, py)
            v = interp_v(grid_v, px, py)

            nx = ft.normal[k,1]
            ny = ft.normal[k,2]

            vn[k] = u*nx + v*ny

            vn[k] -=  num.mass_transfer_rate_imposed_value / num.rho1
        end

    elseif vel_mode === :peskin
        @inbounds for k in 1:n
            vn[k] = _peskin_interp(grid_p.V, x[k], y[k], Δ, xc1, yc1, nx, ny)
        end

    else
        error("ft_step!: unknown vel_mode=$vel_mode. Use :exact, :bilinear, or :peskin.")
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Advect markers — Heun RK2 in normal direction
# ─────────────────────────────────────────────────────────────────────────────
function _advect_markers!(ft::FrontTracker, grid_p, grid_u, grid_v, num, vel_mode::Symbol)
    dt = num.timestep_n
    n  = ft.n_markers

    nx1, ny1 = _marker_normals(ft.x, ft.y, n)
    _marker_velocities!(ft.vn_markers, ft.x, ft.y, n, grid_p, grid_u,grid_v,num, vel_mode,ft)
    vx1 = ft.vn_markers .* nx1
    vy1 = ft.vn_markers .* ny1

    # Predictor positions
    x_s = ft.x .+ dt .* vx1
    y_s = ft.y .+ dt .* vy1

    # Stage 2 — normals and velocity at predicted positions
    nx2, ny2 = _marker_normals(x_s, y_s, n)
    vn2 = similar(ft.vn_markers)

    # Build a lightweight temp struct sharing grid metadata
    # ft_s = FrontTracker(x_s, y_s, ft.connectivities, n, ft.n_seg,
    #                     vn2, ft.nx, ft.ny, ft.Δ, ft.xc1, ft.yc1)
    _marker_velocities!(vn2, x_s, y_s, n, grid_p, grid_u,grid_v, num, vel_mode,ft)

    vx2 = vn2 .* nx2
    vy2 = vn2 .* ny2

    @. ft.x = ft.x + 0.5*dt*(vx1 + vx2)
    @. ft.y = ft.y + 0.5*dt*(vy1 + vy2)
end

# ─────────────────────────────────────────────────────────────────────────────
# Redistribute markers to equal arc-length spacing
# ─────────────────────────────────────────────────────────────────────────────
function _redistribute!(ft::FrontTracker)
    n  = ft.n_markers
    ds = [sqrt((ft.x[mod1(k+1,n)] - ft.x[k])^2 +
               (ft.y[mod1(k+1,n)] - ft.y[k])^2) for k in 1:n]
    L  = sum(ds);  ds_t = L / n
    s_cum = [0.0; cumsum(ds)]
    x_new = similar(ft.x);  y_new = similar(ft.y)
    j = 1
    for k in 1:n
        s_k = (k-1)*ds_t
        while j < n && s_cum[j+1] < s_k;  j += 1;  end
        t  = ds[j] > 1e-14 ? (s_k - s_cum[j]) / ds[j] : 0.0
        t  = clamp(t, 0.0, 1.0)
        jp = mod1(j+1, n)
        x_new[k] = (1-t)*ft.x[j] + t*ft.x[jp]
        y_new[k] = (1-t)*ft.y[j] + t*ft.y[jp]
    end
    copyto!(ft.x, x_new);  copyto!(ft.y, y_new)
    # Rebuild sequential connectivity
    for k in 1:ft.n_seg
        ft.connectivities[k,1] = k
        ft.connectivities[k,2] = mod1(k+1, ft.n_markers)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Recompute signed-distance φ from polygonal front (exact, no drift)
# ─────────────────────────────────────────────────────────────────────────────
function _recompute_ls!(grid_p, ft::FrontTracker)
    n = ft.n_markers;  Δ = ft.Δ
    xc1 = ft.xc1;  yc1 = ft.yc1
    φ = grid_p.LS[1].u

    @inbounds for j in 1:ft.ny, i in 1:ft.nx
        px = xc1 + (i-1)*Δ
        py = yc1 + (j-1)*Δ
        d_min = Inf
        for k in 1:n
            kp = mod1(k+1, n)
            ax, ay = ft.x[k],  ft.y[k]
            bx, by = ft.x[kp], ft.y[kp]
            abx = bx-ax;  aby = by-ay
            t   = clamp((( px-ax)*abx + (py-ay)*aby) / (abx^2+aby^2+1e-14), 0.0, 1.0)
            d   = sqrt((px - ax - t*abx)^2 + (py - ay - t*aby)^2)
            d < d_min && (d_min = d)
        end
        # Sign via ray casting
        inside = false
        for k in 1:n
            kp = mod1(k+1, n)
            y1, y2 = ft.y[k], ft.y[kp]
            x1, x2 = ft.x[k], ft.x[kp]
            if (y1 <= py < y2) || (y2 <= py < y1)
                x_cross = x1 + (py-y1)/(y2-y1+1e-14)*(x2-x1)
                px < x_cross && (inside = !inside)
            end
        end
        # φ[j,i] = inside ? -d_min : d_min
        φ[j,i] = inside ? d_min : -d_min
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Public API
# ─────────────────────────────────────────────────────────────────────────────

"""
    ft_step!(ft, num, grid_p, grid_u, grid_v;
             vel_mode=:exact, recompute_ls=true, redistribute=true)

One front-tracking step.
After return, grid_p.LS[1].u holds updated φ (if recompute_ls=true).
Call update_all_ls_data immediately after.
"""
function ft_step!(ft::FrontTracker, num, grid_p, grid_u, grid_v;
    vel_mode::Symbol   = :exact,
    recompute_ls::Bool = true,
    redistribute::Bool = true)
    _advect_markers!(ft, grid_p, grid_u, grid_v, num, vel_mode)
    redistribute && _redistribute!(ft)
    recompute_ls && _recompute_ls!(grid_p, ft)
end

"""
    ft_radius(ft) → Float64

Mean distance of markers from their centroid.
"""
function ft_radius(ft::FrontTracker)
    cx = sum(ft.x) / ft.n_markers
    cy = sum(ft.y) / ft.n_markers
    return sum(sqrt.((ft.x .- cx).^2 .+ (ft.y .- cy).^2)) / ft.n_markers
end

"""
    ft_pdi!(ft, num, grid_p, grid_u, grid_v; event="update_levelset")

Write marker positions, connectivity, barycentre, and level-set to PDI.
Matches the exact field names used by convert_interfacial_D_to_segments output.

"""
function ft_pdi!(ft::FrontTracker, num, grid_p, grid_u, grid_v;
                 event::String = "write_front_tracking")
    n_m = Int64(ft.n_markers)
    n_s = Int64(ft.n_seg)
    bx  = sum(ft.x) / n_m
    by  = sum(ft.y) / n_m

    radius = 0.0
    @inbounds for i in 1:n_m
        dx = ft.x[i] - bx
        dy = ft.y[i] - by
        radius += sqrt(dx*dx + dy*dy)
    end
    radius /= n_m



    # Flatten connectivity row-major for PDI: [v1_seg1, v2_seg1, v1_seg2, ...]
    conn_flat = Vector{Int64}(vec(ft.connectivities'))

    # PDI_status = @ccall "libpdi".PDI_multi_expose(
    #     event::Cstring,
    #     "nstep"::Cstring,              num.current_iter::Ref{Clonglong},        PDI_OUT::Cint,
    #     "time"::Cstring,               num.time::Ref{Cdouble},                  PDI_OUT::Cint,
    #     # "levelset_p"::Cstring,         grid_p.LS[num.iLSpdi].u::Ptr{Cdouble},  PDI_OUT::Cint,
    #     # "levelset_u"::Cstring,         grid_u.LS[num.iLSpdi].u::Ptr{Cdouble},  PDI_OUT::Cint,
    #     # "levelset_v"::Cstring,         grid_v.LS[num.iLSpdi].u::Ptr{Cdouble},  PDI_OUT::Cint,
    #     "ft_n_markers"::Cstring,       n_m::Ref{Clonglong},                     PDI_OUT::Cint,
    #     "ft_n_seg"::Cstring,           n_s::Ref{Clonglong},                     PDI_OUT::Cint,
    #     "ft_marker_x"::Cstring,        ft.x::Ptr{Cdouble},                      PDI_OUT::Cint,
    #     "ft_marker_y"::Cstring,        ft.y::Ptr{Cdouble},                      PDI_OUT::Cint,
    #     "ft_connectivities"::Cstring,  conn_flat::Ptr{Clonglong},               PDI_OUT::Cint,
    #     # "barycenter_x_coord"::Cstring, bx::Ref{Cdouble},                        PDI_OUT::Cint,
    #     # "barycenter_y_coord"::Cstring, by::Ref{Cdouble},                        PDI_OUT::Cint,
    #     C_NULL::Ptr{Cvoid})::Cint
    ft.normal[:, 1], ft.normal[:, 2] = _marker_normals(ft.x, ft.y, ft.n_markers)

    PDI_status = @ccall "libpdi".PDI_multi_expose("write_front_tracking"::Cstring,
    "nstep"::Cstring, num.current_iter ::Ref{Clonglong}, PDI_OUT::Cint,
    "time"::Cstring, num.time::Ref{Cdouble}, PDI_OUT::Cint,                                 
    "intfc_vtx_num"::Cstring, ft.n_markers::Ref{Clonglong}, PDI_OUT::Cint, 
    "intfc_seg_num"::Cstring, ft.n_markers::Ref{Clonglong}, PDI_OUT::Cint, 
    "intfc_vtx_x"::Cstring, ft.x::Ptr{Cdouble}, PDI_OUT::Cint,
    "intfc_vtx_y"::Cstring, ft.y::Ptr{Cdouble}, PDI_OUT::Cint,
    "intfc_normal"::Cstring, ft.normal::Ptr{Cdouble}, PDI_OUT::Cint,
    # "intfc_vtx_field"::Cstring, intfc_vtx_field::Ptr{Cdouble}, PDI_OUT::Cint,
    "intfc_vtx_connectivities"::Cstring, conn_flat::Ptr{Clonglong}, PDI_OUT::Cint,
    # "barycenter_x_coord"::Cstring, barycenter_x_coord::Ref{Cdouble}, PDI_OUT::Cint,
    "radius"::Cstring, radius::Ref{Cdouble}, PDI_OUT::Cint,
    C_NULL::Ptr{Cvoid})::Cint

    return PDI_status
end
