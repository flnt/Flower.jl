# -----------------------------
# u-momentum x-flux (east/west)
# u is at (i-1/2, j) when p is at (i,j)
# -----------------------------
@inline function Fx_u_upwind(u, j, i)
    ax = 0.5 * (u[j, i] + u[j, i+1])   # average velocity at face #TODO variable spacing dx...
    return ax ≥ 0 ? u[j, i] * ax : u[j, i+1] * ax
end

# -----------------------------
# u-momentum y-flux (north/south)
# u at (i-1/2, j)
# v at (i, j-1/2)
# -----------------------------
@inline function Fy_u_upwind(u, v, j, i)
    vy = 0.5 * (v[j+1, i-1] + v[j+1, i]) # interpolate v to u-face
    return vy ≥ 0 ? u[j, i] * vy : u[j+1, i] * vy
end

# -----------------------------
# v-momentum x-flux (east/west)
# u at (i-1/2, j)
# v at (i, j-1/2)
# -----------------------------
@inline function Fx_v_upwind(u, v, j, i)
    u_face = 0.5 * ( u[j,   i+1] + u[j-1, i+1] )  # interpolate u to v-face
    return u_face >= 0  ? u_face * v[j, i] : u_face * v[j, i+1]
end

# -----------------------------
# v-momentum y-flux (north/south)
# v at (i, j-1/2) when p is at (i,j)
# -----------------------------
@inline function Fy_v_upwind(v, j, i)
    v_face = 0.5 * (v[j, i] + v[j+1, i])   # average at north face
    return v_face ≥ 0 ? v_face*v[j, i] : v_face*v[j+1, i]
end


"""
    compute_fluxes_upwind(num,u, v, rho, dx_p, dy_p)

Compute the convective fluxes for a staggered grid system with variable grid spacing
using an upwind scheme.

# Arguments
- `u::Matrix{Float64}`: x-velocity component (staggered in x)
- `v::Matrix{Float64}`: y-velocity component (staggered in y)
- `rho::Matrix{Float64}`: Density field (defined at cell centers)
- `dx_p::Vector{Float64}`: x-direction grid spacing for the p grid (cell-centered)
- `dy_p::Vector{Float64}`: y-direction grid spacing for the p grid (cell-centered)

# Returns
- `conv_x::Matrix{Float64}`: x-direction flux (staggered in x)
- `conv_y::Matrix{Float64}`: y-direction flux (staggered in y)
"""
function compute_fluxes_upwind(num,u, v, rho_u, rho_v, dx_u, dy_u ,dx_v, dy_v,grid_u,grid_v)
  
    # Initialize flux arrays with proper staggering
    conv_x = zeros(grid_u)  # conv_x is staggered in x (one less column than u)
    conv_y = zeros(grid_v)   # conv_y is staggered in y (one less row than v)

    

    # Compute convective flux of u in x-direction: ∂(u u)/∂x + ∂(v u)/∂y
    # Compute x-direction flux (u at (i-1/2, j))
    for j in 1:grid_u.ny
        for i in 1:grid_u.nx
            # Calculate u velocities at cell faces
            # u_east = 0.5 * (u[j, i+1] + u[j, i])
            # u_west = 0.5 * (u[j, i]   + u[j, i-1])

            # u_north = 0.5 * (u[j+1, i] + u[j, i]) # u at y = y_{j+1/2}, x = x_{i+1/2}
            # u_south = 0.5 * (u[j, i  ] + u[j-1, i]) # u at y = y_{j-1/2}, x = x_{i+1/2}


            # # Calculate v velocities at cell faces
            # v_north = 0.5 * (v[j+1, i-1] + v[j+1, i])
            # v_south = 0.5 * (v[j,   i-1] + v[j,   i])

            # Compute flux in x-direction (integral approximation)
            # Using the density at the cell center rho[j, i]
            rho_c = rho_u[j, i]
            
            # try
            #     dy = dy_u[j,i-1]
            # catch e
            #     error("\n dy ",j,i,size(dy_u))
            # end

            # try
            #     dx = dx_u[j,i]
            # catch e
            #     error("\n dy ",j,i,size(dy_u))
            # end


            # Fx_plus  = Fx_u_upwind(u, i,   j)
            # Fx_minus = Fx_u_upwind(u, i-1, j)
            # Fy_plus  = Fy_u_upwind(u, v, i, j)
            # Fy_minus = Fy_u_upwind(u, v, i, j-1)

            Fx_plus  = Fx_u_upwind(u, j, i)      # east face: u^2 at (i+1/2, j)
            Fx_minus = Fx_u_upwind(u, j, i-1)    # west face: u^2 at (i-1/2, j)

            Fy_plus  = Fy_u_upwind(u, v, j, i)   # north face: u*v at (i-1/2, j+1/2)
            Fy_minus = Fy_u_upwind(u, v, j-1, i) # south face: u*v at (i-1/2, j-1/2)

            # print("\n i j",i," ",j)
            if num.non_dimensionalize == 2 
                conv_x[j,i] =rho_c * (
                    (Fx_plus - Fx_minus ) / dx_u[j,i]+
                    (Fy_plus - Fy_minus ) / dy_u[j,i]
                ) 

            elseif num.non_dimensionalize == 0
                conv_x[j,i] = (Fx_plus - Fx_minus ) * dy_u[j,i] + (Fy_plus - Fy_minus ) * dx_u[j,i]
                

            else
                conv_x[j,i] =rho_c * (
                    (Fx_plus * dy_u[j,i] - Fx_minus * dy_u[j,i-1]) +
                    (Fy_plus * dx_u[j,i] - Fy_minus * dx_u[j-1,i])
                ) 
            end

            # conv_x[j, i] = rho_c * (
            #     ( u_east * u_east * dy_u[j,i] - u_west * u_west * dy_u[j,i-1] ) +
            #     ( v_north * u_north * dx_u[j,i] - v_south * u_south * dx_u[j-1,i] )
            # )

        end
    end

    # Compute convective flux of v in y-direction: ∂(uv)/∂x + ∂(vv)/∂y
    for j in 1:grid_v.ny
        for i in 1:grid_v.nx


            # Calculate v velocities at cell faces
            v_north = 0.5 * (v[j+1, i] + v[j, i])  # y = y_{j+1/2}, x = x_i
            v_south = 0.5 * (v[j,   i] + v[j-1, i])  # y = y_{j-1/2}, x = x_i

            v_east  = 0.5 * (v[j, i+1] + v[j, i])  # x = x_{i+1/2}, y = y_{j-1/2}
            v_west  = 0.5 * (v[j, i]   + v[j, i-1])  # x = x_{i-1/2}, y = y_{j-1/2}

            # u interpolated to east/west faces of v control volume
            u_east = 0.5 * ( u[j,   i+1] + u[j-1, i+1] )  # y = y_{j-1/2}, x = x_{i+1/2}
            u_west = 0.5 * ( u[j,   i]   + u[j-1, i]   )  # y = y_{j-1/2}, x = x_{i-1/2}

            rho_c = rho_v[j, i]  # density at v control volume center (or interpolated)

            # convective flux
            # Fx_plus  = Fx_v_upwind(u, v, i,   j)   # east face: uv at (i+1/2, j)
            # Fx_minus = Fx_v_upwind(u, v, i-1, j)   # west face: uv at (i-1/2, j)

            # Fy_plus  = Fy_v_upwind(v, i, j)        # north face: v^2 at (i, j+1/2)
            # Fy_minus = Fy_v_upwind(v, i, j-1)      # south face: v^2 at (i, j-1/2)
            Fx_plus  = Fx_v_upwind(u, v, j, i)      # east face: u*v at (i+1/2, j)
            Fx_minus = Fx_v_upwind(u, v, j, i-1)    # west face: u*v at (i-1/2, j)

            Fy_plus  = Fy_v_upwind(v, j, i)         # north face: v^2 at (i, j+1/2)
            Fy_minus = Fy_v_upwind(v, j-1, i)       # south face: v^2 at (i, j-1/2)

            if num.non_dimensionalize == 2 
               conv_y[j,i] = rho_c * (
                    (Fx_plus - Fx_minus) / dx_v[j,i] +
                    (Fy_plus - Fy_minus ) / dy_v[j,i] 
                )

            elseif num.non_dimensionalize == 0
                conv_y[j,i] = (Fx_plus - Fx_minus) * dy_v[j,i] + (Fy_plus - Fy_minus ) * dx_v[j,i] 
                

            else
                conv_y[j,i] = rho_c * (
                    (Fx_plus  * dy_v[j,i] - Fx_minus * dy_v[j,i-1]) +
                    (Fy_plus  * dx_v[j,i] - Fy_minus * dx_v[j-1,i])
                )
            end

            # conv_y[j, i] = rho_c * (
            #     ( u_east * v_east * dy_v[j,i] - u_west * v_west  * dy_v[j,i-1]) +   # x-flux
            #     ( v_north * v_north * dx_v[j,i]- v_south * v_south * dx_v[j-1,i])  # y-flux
            # )
            
            # if j == 1
            #     # print("\n check flux ",j," ",i," ",u_east," ",v_east," ",  dy_v[j,i]," ",u_west ," ", v_west  ," ", dy_v[j,i-1] ," ", v_north ," ", dx_v[j,i] ," ", v_south ," ", dx_v[j-1,i] )
            #     printstyled(color=:green, @sprintf "\n i %.5i j %.5i flux %.2e u_east %.2e v_east %.2e u_west %.2e v_west %.2e v_north %.2e v_south %.2e dy_v[j,i-1] %.2e dx_v[j-1,i] %.2e \n" i j conv_y[j, i] u_east v_east u_west v_west v_north v_south dy_v[j,i-1] dx_v[j-1,i])
            # end

        end
    end

   

    return conv_x, conv_y
end



#TODO
"""

Compute the convective fluxes for a staggered grid system with variable grid spacing
using a CUI scheme.

# Arguments
- `u::Matrix{Float64}`: x-velocity component (staggered in x)
- `v::Matrix{Float64}`: y-velocity component (staggered in y)
- `rho::Matrix{Float64}`: Density field (defined at cell centers)
- `dx_p::Vector{Float64}`: x-direction grid spacing for the p grid (cell-centered)
- `dy_p::Vector{Float64}`: y-direction grid spacing for the p grid (cell-centered)

# Returns
- `conv_x::Matrix{Float64}`: x-direction flux (staggered in x)
- `conv_y::Matrix{Float64}`: y-direction flux (staggered in y)
"""
function compute_fluxes_CUI(u, v, rho_u, rho_v, dx_u, dy_u ,dx_v, dy_v,grid_u,grid_v)
  
    # Initialize flux arrays with proper staggering
    conv_x = zeros(grid_u)  # conv_x is staggered in x (one less column than u)
    conv_y = zeros(grid_v)   # conv_y is staggered in y (one less row than v)

    # Compute convective flux of u in x-direction: ∂(u u)/∂x + ∂(v u)/∂y
    # Compute x-direction flux (u at (i-1/2, j))
    for j in 1:grid_u.ny
        for i in 1:grid_u.nx
           
            rho_c = rho_u[j, i]
            
            Fx_plus  = Fx_u_upwind(u, j, i)      # east face: u^2 at (i+1/2, j)
            Fx_minus = Fx_u_upwind(u, j, i-1)    # west face: u^2 at (i-1/2, j)

            Fy_plus  = Fy_u_upwind(u, v, j, i)   # north face: u*v at (i-1/2, j+1/2)
            Fy_minus = Fy_u_upwind(u, v, j-1, i) # south face: u*v at (i-1/2, j-1/2)

            conv_x[j,i] =rho_c * (
                (Fx_plus * dy_u[j,i] - Fx_minus * dy_u[j,i-1]) +
                (Fy_plus * dx_u[j,i] - Fy_minus * dx_u[j-1,i])
            ) 

        end
    end

    # Compute convective flux of v in y-direction: ∂(uv)/∂x + ∂(vv)/∂y
    for j in 1:grid_v.ny
        for i in 1:grid_v.nx
            rho_c = rho_v[j, i]  # density at v control volume center (or interpolated)

            Fx_plus  = Fx_v_upwind(u, v, j, i)      # east face: u*v at (i+1/2, j)
            Fx_minus = Fx_v_upwind(u, v, j, i-1)    # west face: u*v at (i-1/2, j)

            Fy_plus  = Fy_v_upwind(v, j, i)         # north face: v^2 at (i, j+1/2)
            Fy_minus = Fy_v_upwind(v, j-1, i)       # south face: v^2 at (i, j-1/2)

            conv_y[j,i] = rho_c * (
                (Fx_plus  * dy_v[j,i] - Fx_minus * dy_v[j,i-1]) +
                (Fy_plus  * dx_v[j,i] - Fy_minus * dx_v[j-1,i])
            )
            
        end
    end

    return conv_x, conv_y
end