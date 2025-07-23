"""

"""
function update_one_fluid_density_viscosity(num,grid_p,grid_u,grid_v,volume_fraction,levelset_one_fluid,rho_one_fluid,
    rho_one_fluid_u,rho_one_fluid_v,)

    volume_fraction .= grid_p.LS[end].geoL.cap[:,:,5]
    levelset_one_fluid .= grid_p.LS[end].u


    # print("\nvolume fraction update")
    # display(volume_fraction)
    # display(grid_p.LS[end].geoL.cap[:,:,5])


    if num.rho_one_fluid_average == 0 #arithmetic average
        rho_one_fluid .= num.rho1 * grid_p.LS[end].geoL.cap[:,:,5] .+ num.rho2 * (1.0 .- grid_p.LS[end].geoL.cap[:,:,5] )
        rho_one_fluid_u .= num.rho1 * grid_u.LS[end].geoL.cap[:,:,5] .+ num.rho2 * (1.0 .- grid_u.LS[end].geoL.cap[:,:,5] )
        rho_one_fluid_v .= num.rho1 * grid_v.LS[end].geoL.cap[:,:,5] .+ num.rho2 * (1.0 .- grid_v.LS[end].geoL.cap[:,:,5] )
    end

    # if num.mu_one_fluid_average == 0 #arithmetic average
    #     mu_one_fluid  .= (num.mu1 - num.mu2) * grid_p.LS[end].geoL.cap[:,:,5]  .+ num.mu2
    #     mu_one_fluid_u  .= (num.mu1 - num.mu2) * grid_u.LS[end].geoL.cap[:,:,5] .+ num.mu2  
    #     mu_one_fluid_v  .= (num.mu1 - num.mu2) * grid_v.LS[end].geoL.cap[:,:,5] .+ num.mu2  
    # elseif num.mu_one_fluid_average == 1 #harmonic average
    #     mu_one_fluid   .= harmonic_average_one_fluid.(num.mu1,num.mu2,grid_p.LS[end].geoL.cap[:,:,5])
    #     mu_one_fluid_u .= harmonic_average_one_fluid.(num.mu1,num.mu2,grid_u.LS[end].geoL.cap[:,:,5])
    #     mu_one_fluid_v .= harmonic_average_one_fluid.(num.mu1,num.mu2,grid_v.LS[end].geoL.cap[:,:,5])
    # end

    # PDI_status = @ccall "libpdi".PDI_multi_expose("print_one_fluid"::Cstring,
    # "rho_one_fluid"::Cstring, rho_one_fluid::Ptr{Cdouble}, PDI_OUT::Cint,
    # "mu_one_fluid"::Cstring, mu_one_fluid::Ptr{Cdouble}, PDI_OUT::Cint,
    # C_NULL::Ptr{Cvoid})::Cint

    PDI_status = @ccall "libpdi".PDI_multi_expose("write_one_fluid"::Cstring,
    "nstep"::Cstring, num.current_i ::Ref{Clonglong}, PDI_OUT::Cint,
    "rho_one_fluid"::Cstring, rho_one_fluid::Ptr{Cdouble}, PDI_OUT::Cint,
    # "mu_one_fluid"::Cstring, mu_one_fluid::Ptr{Cdouble}, PDI_OUT::Cint,
    "volume_fraction"::Cstring, volume_fraction::Ptr{Cdouble}, PDI_OUT::Cint,
    C_NULL::Ptr{Cvoid})::Cint

end

function harmonic_average_one_fluid(mu1,mu2, volume_fraction)
    return (mu1*mu2) / (mu2 * volume_fraction + (1.0 - volume_fraction) * mu1)
end


"""
interpolate based on four node values (rectangle)
"""
function bilinear_interpolation(x, y, x1, y1, x2, y2, Q11, Q12, Q21, Q22)
        # Calculate the intermediate terms
        term1 = (x2 - x) * (y2 - y) * Q11 / ((x2 - x1) * (y2 - y1))
        term2 = (x - x1) * (y2 - y) * Q21 / ((x2 - x1) * (y2 - y1))
        term3 = (x2 - x) * (y - y1) * Q12 / ((x2 - x1) * (y2 - y1))
        term4 = (x - x1) * (y - y1) * Q22 / ((x2 - x1) * (y2 - y1))

        # print("\n term1 ", term1," ", term2, " ", term3," ",term4 , " ",((x2 - x1) * (y2 - y1)))
        # Sum the terms to get the interpolated value
        return term1 + term2 + term3 + term4
    end


"""
bilinear interpolation, based on grid_p, assuming constant mesh spacing
"""
function bilinear_interpolation(grid_p, x, y,values)
    dx = grid_p.dx[2,2] #constant dx
    dy = grid_p.dy[2,2] #constant dx
    # print("\n dx dy ",dx," dy ",dy)
    # Calculate the indices and weights for interpolation
    i0 = floor(Int, x / dx) #TODO
    j0 = floor(Int, y / dy)
    i1 = i0 + 1
    j1 = j0 + 1
    
    print("\nindices "," i0 ",i0," i1 ",i1," j0 ",j0," j1 ",j1)
    print("\ngrid "," i0 j0 ",grid_p.x[j0,i0]," i1 j0 ",grid_p.x[j0,i1]," i0 j1 ",grid_p.x[j1,i0]," i1 j1 ",grid_p.x[j1,i1])

    # Calculate the weights
    wx = (x - grid_p.x[j0,i0]) / dx
    wy = (y - grid_p.y[j0,i0]) / dy

    # Perform bilinear interpolation
    value = (1 - wx) * (1 - wy) * values[j0,i0] +
            wx * (1 - wy) * values[j0,i1] +
            (1 - wx) * wy * values[j1,i0] +
            wx * wy * values[j1, i1]

    return value
end


"""
store every nodes (border included) in a 2D matrix for interpolations
uses geoL, not geoS
"""
function create_2D_grid_x(grid_p,add_x=true,add_y=true)

    if add_x
        nx = grid_p.nx+2
        rangex = 2:grid_p.nx+1
    else
        nx = grid_p.nx
        rangex = 1:grid_p.nx
    end

    if add_y
        ny = grid_p.ny+2
        rangey = 2:grid_p.ny+1
    else
        ny = grid_p.ny
        rangey= 1:grid_p.ny
    end

    #all_grid_v_nodes_2D_x_for_dv_dx_interp = create_2D_grid_x(grid_v,true,false)
    # print("\n create_2D_grid_x ", nx," ",ny)

    create_2D_grid = zeros(ny,nx)
    

    x_centroid = grid_p.x .+ getproperty.(grid_p.LS[1].geoL.centroid, :x) .* grid_p.dx #geoS
    # y_centroid = grid_p.y .+ getproperty.(grid_p.LS[1].geoS.centroid, :y) .* grid_p.dy



    # print("\n x_centroid")
    # display(x_centroid)

    create_2D_grid[rangey,rangex] = x_centroid #grid_p.x

    x_bc_left = grid_p.x[:,1] .- grid_p.dx[:,1] ./ 2.0

    # y_bc_bottom = grid_p.y[1,:] .- grid_p.dy[1,:] ./ 2.0

    # y_bc_top = grid_p.y[end,:] .+ grid_p.dy[end,:] ./ 2.0

    x_bc_right = grid_p.x[:,end] .+ grid_p.dx[:,end] ./ 2.0

    # create_2D_grid[1,2:grid_p.nx] = create_2D_grid[2,2:grid_p.nx]

    # create_2D_grid[end,2:grid_p.nx] = create_2D_grid[end-1,2:grid_p.nx]

    # display(create_2D_grid)

    if add_x
        # create_2D_grid[2:grid_p.ny+1,1] = x_bc_left
        # create_2D_grid[2:grid_p.ny+1,end] = x_bc_right
        create_2D_grid[rangey,1] = x_bc_left
        create_2D_grid[rangey,end] = x_bc_right
    end

    if add_y
        create_2D_grid[1,:] = create_2D_grid[2,:]
        create_2D_grid[end,:] = create_2D_grid[end-1,:]
    end

    return create_2D_grid
end


"""
store every nodes (border included) in a 2D matrix for interpolations
uses geoL, not geoS
"""
function create_2D_grid_y(grid_p,add_x=true,add_y=true)
    
    if add_x
        nx = grid_p.nx+2
        rangex = 2:grid_p.nx+1
    else
        nx = grid_p.nx
        rangex = 1:grid_p.nx
    end

    if add_y
        ny = grid_p.ny+2
        rangey = 2:grid_p.ny+1
    else
        ny = grid_p.ny
        rangey= 1:grid_p.ny
    end

    
    create_2D_grid = zeros(ny,nx)
    
    # x_centroid = grid_p.x .+ getproperty.(grid_p.LS[1].geoS.centroid, :x) .* grid_p.dx
    y_centroid = grid_p.y .+ getproperty.(grid_p.LS[1].geoL.centroid, :y) .* grid_p.dy #geoS

    # print("\n y centroid")

    # display(y_centroid)

    create_2D_grid[rangey,rangex] = y_centroid #grid_p.x

    # x_bc_left = grid_p.x[:,1] .- grid_p.dx[:,1] ./ 2.0

    y_bc_bottom = grid_p.y[1,:] .- grid_p.dy[1,:] ./ 2.0

    y_bc_top = grid_p.y[end,:] .+ grid_p.dy[end,:] ./ 2.0

    # x_bc_right = grid_p.x[:,end] .+ grid_p.dx[:,end] ./ 2.0

    # create_2D_grid[1,2:grid_p.nx] = create_2D_grid[2,2:grid_p.nx]

    # create_2D_grid[end,2:grid_p.nx] = create_2D_grid[end-1,2:grid_p.nx]

    # display(create_2D_grid)

    # create_2D_grid[2:grid_p.ny+1,1] = x_bc_left

    # create_2D_grid[2:grid_p.ny+1,end] = x_bc_right

    # create_2D_grid[1,:] = create_2D_grid[2,:]

    # create_2D_grid[end,:] = create_2D_grid[end-1,:]

    # create_2D_grid[2:grid_p.ny+1,1] = create_2D_grid[2:grid_p.ny+1,2]
    # create_2D_grid[2:grid_p.ny+1,end] = create_2D_grid[2:grid_p.ny+1,end-1]

    # create_2D_grid[1,2:grid_p.nx+1] = y_bc_bottom

    # create_2D_grid[end,2:grid_p.nx+1] = y_bc_top

    # create_2D_grid[end,1] = create_2D_grid[end,2]
    # create_2D_grid[end,end] = create_2D_grid[end,end-1]


    if add_x
        # create_2D_grid[2:grid_p.ny+1,1] = create_2D_grid[2:grid_p.ny+1,2]
        # create_2D_grid[2:grid_p.ny+1,end] = create_2D_grid[2:grid_p.ny+1,end-1]
        create_2D_grid[rangey,1] = create_2D_grid[rangey,2]
        create_2D_grid[rangey,end] = create_2D_grid[rangey,end-1]
    end

    if add_y
        create_2D_grid[1,rangex] = y_bc_bottom
        create_2D_grid[end,rangex] = y_bc_top

        if add_x
            #top corner
            create_2D_grid[end,1] = create_2D_grid[end,2]
            create_2D_grid[end,end] = create_2D_grid[end,end-1]
        end
    end

    return create_2D_grid
end


"""
store every nodes (border included) in a 2D matrix for interpolations
"""
function create_2D_grid_volume_fraction(grid_p,volume_fraction)
    #assuming no contact angle method
    create_2D_grid = zeros(grid_p.ny+2,grid_p.nx+2)

    # x_centroid = grid_p.x .+ getproperty.(grid_p.LS[1].geoS.centroid, :x) .* grid_p.dx
    # y_centroid = grid_p.y .+ getproperty.(grid_p.LS[1].geoS.centroid, :y) .* grid_p.dy

    create_2D_grid[2:grid_p.ny+1,2:grid_p.nx+1] = volume_fraction


    # x_bc_left = grid_p.x[:,1] .- grid_p.dx[:,1] ./ 2.0

    # y_bc_bottom = grid_p.y[1,:] .- grid_p.dy[1,:] ./ 2.0

    # y_bc_top = grid_p.y[end,:] .+ grid_p.dy[end,:] ./ 2.0

    # x_bc_right = grid_p.x[:,end] .+ grid_p.dx[:,end] ./ 2.0

    # # create_2D_grid[1,2:grid_p.nx] = create_2D_grid[2,2:grid_p.nx]

    # # create_2D_grid[end,2:grid_p.nx] = create_2D_grid[end-1,2:grid_p.nx]

    # # display(create_2D_grid)

    # create_2D_grid[2:grid_p.ny+1,1] = x_bc_left

    # create_2D_grid[2:grid_p.ny+1,end] = x_bc_right

    create_2D_grid[1,:] = create_2D_grid[2,:]

    create_2D_grid[end,:] = create_2D_grid[end-1,:]

    create_2D_grid[:,1] = create_2D_grid[:,2]

    create_2D_grid[:,end] = create_2D_grid[:,end-1]

    create_2D_grid[1,1] = (create_2D_grid[1,2] + create_2D_grid[2,1])/2
    create_2D_grid[1,end] = (create_2D_grid[1,end-1] + create_2D_grid[2,end])/2
    create_2D_grid[end,1] = (create_2D_grid[end,2] + create_2D_grid[end-1,2])/2
    create_2D_grid[end,end] = (create_2D_grid[end,end-1] + create_2D_grid[end-1,end])/2

    #TODO if contact angle
    # Q11 = volume_fraction[j-1,i-1]
    # Q12 = volume_fraction[j,i-1]
    # Q21 = vecb_R(volume_fraction_1D,grid_p)[j-1]
    # Q22 = vecb_R(volume_fraction_1D,grid_p)[j]


    return create_2D_grid
end



"""
Compute the a smooth heavyside from levelset

# Arguments
- `phi`: A scalar value representing the level set function value.
- `epsilon`: A small positive number defining the width of the transition region.

"""
function levelset_heavyside(phi, epsilon)
    if phi < -epsilon
        return 0
    elseif phi > epsilon
        return 1
    else
        return 0.5 * (1 + phi / epsilon + sin(pi * phi / epsilon) / pi)
    end
end




"""
Smooth the volume fraction in the hope of improving the computation of the curvature
"""
function smooth_vof_2d!(grid_p,vof_field, num_smoothings,smoothed_vof)
    # Get the dimensions of the VOF field
    @unpack nx, ny = grid_p

    # Create a copy of the VOF field to store the smoothed values

    # kernel_size
    # Define the kernel radius
    # r = kernel_size ÷ 2

    start_x = 2
    end_x = nx -1

    start_y = 2
    end_y = ny -1

    smoothed_vof = copy(vof_field)
    smoothed_vof_prev = similar(smoothed_vof)



    for iter in 1:num_smoothings
    
    
        # print("\n nx ny ",nx, " ", ny, " ",size(vof_field), " r ",r)
        # Iterate over each cell in the VOF field
        smoothed_vof_prev .= smoothed_vof

    


        for j in start_y:end_y
            for i in start_x:end_x

                # # Initialize the sum and count for the averaging
                # sum = 0.0
                # count = 0

                # # Iterate over the neighboring cells within the kernel
                # for di in -r:r
                #     for dj in -r:r
                #         # Calculate the neighboring cell indices
                #         ni, nj = i + di, j + dj
                        
                #         # print("\n Smoothing VOF ",i," ",j," ",ni," ",nj)
                #         # Check if the neighboring cell is within the bounds
                #         if 1 ≤ ni ≤ nx && 1 ≤ nj ≤ ny
                #             # Add the VOF value of the neighboring cell to the sum
                #             sum += vof_field[nj, ni]
                #             count += 1
                #         end
                #     end
                # end
                # # Calculate the smoothed VOF value
                # smoothed_vof[j, i] = sum / count
                # print("\n smoothed ",sum/count,"vof_field ",vof_field[j,i]," ",count,)

                smoothed_vof[j,i] = smoothed_vof_prev[j,i]/2.0 + (smoothed_vof_prev[j,i-1]+smoothed_vof_prev[j-1,i]+smoothed_vof_prev[j+1,i]+smoothed_vof_prev[j,i+1])/8.0
                
                print("\n smoothed ","vof_field ",vof_field[j,i]," ",smoothed_vof[j,i])

            end
        end

        print("\n smoothed ")
        display(smoothed_vof_prev)
        display(smoothed_vof)
        display(vof_field)
    end
    # return smoothed_vof
    print("\n smoothed ")
    display(smoothed_vof)
end


"""
solves Navier-Stokes equations with a pressure projection method. 


#### Variables and Data Structures
- `vec1(ucorrD, grid_u)`: Velocity correction for the horizontal grid_p.
- `vec1(vcorrD, grid_v)`: Velocity correction for the vertical grid_p.
- `vec1(rhs_ϕ, grid_p)`: Right-hand side of the Poisson equation.
- `vec1(pD, grid_p)`: Pressure correction.
- `vec1(uD, grid_u)`: Updated horizontal velocity.
- `vec1(vD, grid_v)`: Updated vertical velocity.
- `vec1(result_p, grid_p)`: Pressure correction potential.
- `ϕ`: Pressure correction potential.
- `u`: Updated horizontal velocity.
- `v`: Updated vertical velocity.
- `p`: Pressure.
- `opC_p`, `opC_u`, `opC_v`: Operator matrices for pressure, horizontal velocity, and vertical velocity, respectively.
- `geo`, `geo_u`, `geo_v`: Geometric data for the grid_p.
- `bc_int`: Interfacial boundary conditions.
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
   - `vecb(vcorrD, grid_v) .= uvD[border_v_velocity]`: Updates the vertical velocity correction.
   - `kill_dead_cells!(vec1(vcorrD,grid_v), grid_v, geo_v[end])`: Removes dead cells from the vertical velocity correction grid_p.
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
function pressure_projection_one_fluid!(
    time_scheme, bc_int,
    num, grid_p, geo, grid_u, geo_u, grid_v, geo_v, ph,
    BC_u, BC_v, BC_p,
    opC_p, opC_u, opC_v, op_conv,
    Au, Bu, Av, Bv, Aϕ, Auv, Buv,rhs_uv,
    Lpm1, bc_Lpm1, bc_Lpm1_b, Lum1, bc_Lum1, bc_Lum1_b, Lvm1, bc_Lvm1, bc_Lvm1_b,
    Cum1, Cvm1, Mum1, Mvm1,
    periodic_x, periodic_y, advection, ls_advection, current_i, Ra, navier,
    volume_fraction,
    levelset_one_fluid,
    rho_one_fluid, 
    # mu_one_fluid,
    rho_one_fluid_u, 
    # mu_one_fluid_u,
    rho_one_fluid_v, 
    # mu_one_fluid_v,
    tmp_vec_p,
    tmp_vec_p0,
    pres_free_suface,jump_mass_flux,mass_flux
    )
    @unpack Re, τ, σ, g, β, nLS, nNavier = num
    @unpack p, pD, ϕ, u, v, ucorrD, vcorrD, uD, vD, ucorr, vcorr, uT = ph
    @unpack Cu, Cv, CUTCu, CUTCv = op_conv

    iτ = 1.0 / τ
    # irho1 = 1.0 ./ rho_one_fluid
    # mu1_over_rho1 = mu_one_fluid ./ rho_one_fluid
    # mu1_over_rho1 = num.mu1 / num.rho1 


    # II = CartesianIndex(div(grid_v.ny,2),1)
    # pII = lexicographic(II,grid_v.ny)
    # print("\nLv[pII,:] ",Lvm1[pII,:])
    # print("\bc_Lvm1_b[pII,:] ",bc_Lvm1_b[pII,:])

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

        ∇ϕ_x = opC_u.AxT * opC_u.Rx * vec1(pD,grid_p) .+ opC_u.Gx_b * vecb(pD,grid_p)
        ∇ϕ_y = opC_v.AyT * opC_v.Ry * vec1(pD,grid_p) .+ opC_v.Gy_b * vecb(pD,grid_p)
        #region cut-cell
        # for iLS in 1:nLS
        #     ∇ϕ_x .+= opC_u.Gx[iLS] * veci(pD,grid_p,iLS+1)
        #     ∇ϕ_y .+= opC_v.Gy[iLS] * veci(pD,grid_p,iLS+1)
        # end
        #endregion cut-cell

        ph.Gxm1 .= 0.0 #TODO
        ph.Gym1 .= 0.0
        
        # gradient \times volume of cell (cells at border: volume = dx*dy/2)
        ph.Gxm1 .= copy(∇ϕ_x) #∇ϕ_x
        ph.Gym1 .= copy(∇ϕ_y) #∇ϕ_y

        # gradient check
        PDI_status = @ccall "libpdi".PDI_multi_expose("print_pressure_gradient_in_prediction"::Cstring,
        "grad_x_1D"::Cstring, ph.Gxm1::Ptr{Cdouble}, PDI_OUT::Cint,
        "grad_y_1D"::Cstring, ph.Gym1::Ptr{Cdouble}, PDI_OUT::Cint,
        "p_1D"::Cstring, ph.pD::Ptr{Cdouble}, PDI_OUT::Cint,
        C_NULL::Ptr{Cvoid})::Cint

        ∇ϕ_x .= 0.0
        ∇ϕ_y .= 0.0

        grad_x = zeros(grid_u)
        grad_y = zeros(grid_v)
        compute_grad_T_x_T_y_array_u_v_capacities!(num, grid_p, grid_u, grid_v, opC_u, opC_v, grad_x, grad_y, ph.pD)
    
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


    nip = grid_p.nx * grid_p.ny

    niu = grid_u.nx * grid_u.ny
    nbu = 2 * grid_u.nx + 2 * grid_u.ny
    # ntu = (nLS - nNavier + 1) * niu + nbu
    ntu = niu + nbu

    niv = grid_v.nx * grid_v.ny
    nbv = 2 * grid_v.nx + 2 * grid_v.ny
    # ntv = (nLS - nNavier + 1) * niv + nbv
    ntv = niv + nbv

    ntNavier = num.nNavier * nip


    # Array indices 
    bulk_u_velocity = 1:niu
    bulk_v_velocity =  ntu+1:ntu+niv
    # bulk_tangential_velocity = 

    border_u_velocity = ntu-nbu+1:ntu
    border_v_velocity = ntu+ntv-nbv+1:ntu+ntv


    if num.prediction == "PmIIimposedpressure" || num.prediction == "PmIIimposedpressureBCincrement" 
        BC_Poisson = Boundaries() #Neumann everywhere
    else
        BC_Poisson = copy(BC_p) 
    end

    if is_Forward_Euler(time_scheme)
        rhs_u, rhs_v, rhs_ϕ, rhs_uv, Lp, bc_Lp, bc_Lp_b, Lu, diffusion_LS_u, diffusion_border_u, Lv, diffusion_LS_v, diffusion_border_v = set_Forward_Euler_one_fluid!(
            bc_int, num, grid_p, geo, grid_u, geo_u, grid_v, geo_v,
            opC_p, opC_u, opC_v, BC_Poisson,BC_u, BC_v,
            Au, Bu, Av, Bv, Aϕ, Auv, Buv,
            volume_fraction,rho_one_fluid_u,rho_one_fluid_v,
            Lpm1, bc_Lpm1, bc_Lpm1_b, Lum1, bc_Lum1, bc_Lum1_b, Lvm1, bc_Lvm1, bc_Lvm1_b,
            Mum1, Mvm1, op_conv, ph,
            periodic_x, periodic_y, advection, ls_advection, navier,rhs_uv,
        )
    elseif is_Crank_Nicolson(time_scheme)
        @error("\n Crank_Nicolson not implemented")
        # rhs_u, rhs_v, rhs_ϕ, Lp, bc_Lp, bc_Lp_b, Lu, diffusion_LS_u, diffusion_border_u, Lv, diffusion_LS_v, diffusion_border_v = set_Crank_Nicolson!(
        #     bc_int, num, grid_p, geo, grid_u, geo_u, grid_v, geo_v,
        #     opC_p, opC_u, opC_v, BC_Poisson, BC_u, BC_v,
        #     Au, Bu, Av, Bv, Aϕ,
        #     Lpm1, bc_Lpm1, bc_Lpm1_b, Lum1, bc_Lum1, bc_Lum1_b, Lvm1, bc_Lvm1, bc_Lvm1_b,
        #     Mum1, Mvm1, mu1_over_rho1, op_conv, ph,
        #     periodic_x, periodic_y, advection, ls_advection
        # )
    end

    # ra_x = Ra .* sin(β) .* opC_u.M * vec(hcat(zeros(grid_u.ny), ph.T))
    # ra_y = Ra .* cos(β) .* opC_v.M * vec(vcat(zeros(1,grid_v.nx), ph.T))

   

    Convu = fzeros(grid_u)
    Convv = fzeros(grid_v)

    #region convection
    if num.non_dimensionalize == 0
        Cui = Cu * vec(u) .+ CUTCu
        Cvi = Cv * vec(v) .+ CUTCv
    elseif num.non_dimensionalize == -1 
        # Element-wise multiplication
        Cui = Cu * vec(rho_one_fluid_u .* u) .+ CUTCu #_rho should be included  #TODO 
        Cvi = Cv * vec(rho_one_fluid_v .* v) .+ CUTCv #_rho
    end
    #endregion convection


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


    if num.one_fluid_model == 1 
        volumic_surface_tension_u = zeros(grid_u)
        volumic_surface_tension_v = zeros(grid_v)

        if num.surface_tension == 0
            compute_surface_tension_VOF!(num,grid_p, grid_u, grid_v, opC_p, opC_u, opC_v, 
            volume_fraction,levelset_one_fluid,volumic_surface_tension_u,volumic_surface_tension_v,tmp_vec_p,tmp_vec_p0)
        elseif num.surface_tension == 1
            compute_surface_tension_LS!(num,grid_p, grid_u, grid_v, opC_p, opC_u, opC_v, 
            volume_fraction,levelset_one_fluid,volumic_surface_tension_u,volumic_surface_tension_v,tmp_vec_p,tmp_vec_p0)
        end
    end
    
    PDI_status = @ccall "libpdi".PDI_multi_expose("write_one_fluid_surface_tension_concise"::Cstring,
    "nstep"::Cstring, num.current_i ::Ref{Clonglong}, PDI_OUT::Cint,
    # "rho_one_fluid"::Cstring, rho_one_fluid::Ptr{Cdouble}, PDI_OUT::Cint,
    # "mu_one_fluid"::Cstring, mu_one_fluid::Ptr{Cdouble}, PDI_OUT::Cint,
    # "volume_fraction"::Cstring, volume_fraction::Ptr{Cdouble}, PDI_OUT::Cint,
    # "grad_u"::Cstring, normal_and_dirac_u::Ptr{Cdouble}, PDI_OUT::Cint,
    # "grad_v"::Cstring, normal_and_dirac_v::Ptr{Cdouble}, PDI_OUT::Cint,
    # "curvature_p"::Cstring, curvature_p::Ptr{Cdouble}, PDI_OUT::Cint,
    # "curvature_u"::Cstring, curvature_u::Ptr{Cdouble}, PDI_OUT::Cint,
    # "curvature_v"::Cstring, curvature_v::Ptr{Cdouble}, PDI_OUT::Cint,
    "volumic_surface_tension_u"::Cstring, volumic_surface_tension_u::Ptr{Cdouble}, PDI_OUT::Cint,
    "volumic_surface_tension_v"::Cstring, volumic_surface_tension_v::Ptr{Cdouble}, PDI_OUT::Cint,
    # "normal_angle"::Cstring, grid_p.LS[iLSpdi].α::Ptr{Cdouble}, PDI_OUT::Cint,
    # "normal_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
    # "normal_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,  
    C_NULL::Ptr{Cvoid})::Cint


    # TODO PDI_multi_expose() #Cu u CUTCu


   

    #region Navier

    #fill u part at bulk_u_velocity
    uvm1 = zeros(ntu + ntv + nNavier * nip)
    uvm1[bulk_u_velocity] .= vec1(uD,grid_u)
    uvm1[bulk_v_velocity] .= vec1(vD,grid_v)
    uvm1[border_u_velocity] .= vecb(uD,grid_u)
    uvm1[border_v_velocity] .= vecb(vD,grid_v)

    #region cut-cell
    # _iLS = 1
    # for iLS in 1:num.nLS
    #     if !is_navier(bc_int[iLS]) && !is_navier_cl(bc_int[iLS])
    #         uvm1[_iLS*niu+1:(_iLS+1)*niu] .= veci(uD,grid_u,iLS+1)
    #         uvm1[ntu+_iLS*niv+1:ntu+(_iLS+1)*niv] .= veci(vD,grid_v,iLS+1)
    #         _iLS += 1
    #     end
    # end
    #endregion cut-cell

    velocity_block = 1:ntu + ntv + nNavier * nip

    printstyled(color=:red, @sprintf "\n temporal\n")

    PDI_status = @ccall "libpdi".PDI_multi_expose("rhs_uv"::Cstring,
    "rhs_uv_len"::Cstring, length(rhs_uv)::Ref{Clonglong}, PDI_OUT::Cint,
    "rhs_uv_1D"::Cstring, rhs_uv::Ptr{Cdouble}, PDI_OUT::Cint,
    C_NULL::Ptr{Cvoid})::Cint

    @views mul!(rhs_uv[velocity_block], Buv[velocity_block,velocity_block], uvm1, 1.0, 1.0)
    

    PDI_status = @ccall "libpdi".PDI_multi_expose("rhs_uv"::Cstring,
    "rhs_uv_len"::Cstring, length(rhs_uv)::Ref{Clonglong}, PDI_OUT::Cint,
    "rhs_uv_1D"::Cstring, rhs_uv::Ptr{Cdouble}, PDI_OUT::Cint,
    C_NULL::Ptr{Cvoid})::Cint

    printstyled(color=:red, @sprintf "\n temporal\n")

    # print("\nAuv")
    # print(Auv)
    # print("\nAuv")
    
    # print("\n  size(Auv.colptr) ", size(Auv.colptr)," ",length(Auv.colptr))
    # print("\n  size(Auv.rowval) ", size(Auv.rowval))
    # print("\n  size(Auv.nzval) ", size(Auv.nzval))
    # print("\n  Auv.n ", Auv.n)
    # print("\n  Auv.m ", Auv.m)



    # PDI_status = @ccall "libpdi".PDI_multi_expose("print_matrix"::Cstring,
    # "Auv_n"::Cstring, Auv.n::Ref{Clonglong}, PDI_OUT::Cint,
    # "Auv_m"::Cstring, Auv.m::Ref{Clonglong}, PDI_OUT::Cint,
    # "Auv_colptr_len"::Cstring, length(Auv.colptr)::Ref{Clonglong}, PDI_OUT::Cint,
    # "Auv_rowval_len"::Cstring, length(Auv.rowval)::Ref{Clonglong}, PDI_OUT::Cint,
    # "Auv_nzval_len"::Cstring, length(Auv.nzval)::Ref{Clonglong}, PDI_OUT::Cint,
    # "Auv_colptr_1D"::Cstring, Auv.colptr::Ptr{Clonglong}, PDI_OUT::Cint,
    # "Auv_rowval_1D"::Cstring, Auv.rowval::Ptr{Clonglong}, PDI_OUT::Cint,
    # "Auv_nzval_1D"::Cstring, Auv.nzval::Ptr{Cdouble}, PDI_OUT::Cint,
    # C_NULL::Ptr{Cvoid})::Cint



    if num.non_dimensionalize == 0
        grav_x = g .* sin(β) .* opC_u.M * fones(grid_u)
        grav_y = g .* cos(β) .* opC_v.M * fones(grid_v)
    elseif num.non_dimensionalize == -1
        diag_rho_u = Diagonal(vec(rho_one_fluid_u))
        diag_rho_v = Diagonal(vec(rho_one_fluid_v))

        # print("\n size diag_rho_u ",size(diag_rho_u))
        # print("\n size grav_x ",size(grav_x)," opC_u.M ", size(opC_u.M))

        grav_x = g .* sin(β) .* opC_u.M * diag_rho_u * fones(grid_u)
        grav_y = g .* cos(β) .* opC_v.M * diag_rho_v * fones(grid_v)
    else
        grav_x = g .* sin(β) .* opC_u.M * fones(grid_u)
        grav_y = g .* cos(β) .* opC_v.M * fones(grid_v)
    end

    # printstyled(color=:red, @sprintf "\n gravity \n")

    # display(grav_y)

    rhs_uv[bulk_u_velocity] .+= τ .* grav_x #τ * rho_one_fluid_u .* grav_x

    rhs_uv[bulk_u_velocity] .-= τ .* Convu #rho in Convu

     PDI_status = @ccall "libpdi".PDI_multi_expose("rhs_uv"::Cstring,
    "rhs_uv_len"::Cstring, length(rhs_uv)::Ref{Clonglong}, PDI_OUT::Cint,
    "rhs_uv_1D"::Cstring, rhs_uv::Ptr{Cdouble}, PDI_OUT::Cint,
    C_NULL::Ptr{Cvoid})::Cint

    # print("\n volumic_surface_tension_u")
    # display(volumic_surface_tension_u)
    # display(volumic_surface_tension_v)

    # print("\n rhs_uv",rhs_uv)

    # rhs_uv[bulk_u_velocity] .+= τ .* ra_x
    if num.pressure_velocity_coupling == 0
        if num.non_dimensionalize == 0
            rhs_uv[bulk_u_velocity] .-= τ .* ph.Gxm1 ./ vec(rho_one_fluid_u)
            #Surface tension
            rhs_uv[bulk_u_velocity] .+= τ .* vec(volumic_surface_tension_u) ./ vec(rho_one_fluid_u)

            print("\n surface tension u")
            PDI_status = @ccall "libpdi".PDI_multi_expose("rhs_uv"::Cstring,
            "rhs_uv_len"::Cstring, length(rhs_uv)::Ref{Clonglong}, PDI_OUT::Cint,
            "rhs_uv_1D"::Cstring, rhs_uv::Ptr{Cdouble}, PDI_OUT::Cint,
            C_NULL::Ptr{Cvoid})::Cint

        else
            rhs_uv[bulk_u_velocity] .-= τ .* ph.Gxm1 
            #Surface tension
            rhs_uv[bulk_u_velocity] .+= τ .* vec(volumic_surface_tension_u)

        end

    end 

   
    rhs_uv[bulk_v_velocity] .+= τ .* grav_y
    rhs_uv[bulk_v_velocity] .-= τ .* Convv
    # rhs_uv[bulk_v_velocity] .+= τ .* ra_y
    if num.pressure_velocity_coupling == 0
         
        if num.non_dimensionalize == 0
            rhs_uv[bulk_v_velocity] .-= τ .* ph.Gym1 ./ vec(rho_one_fluid_v)
            #Surface tension
            rhs_uv[bulk_v_velocity] .+= τ .* vec(volumic_surface_tension_v) ./  vec(rho_one_fluid_v)
        else
            rhs_uv[bulk_v_velocity] .-= τ .* ph.Gym1
            #Surface tension
            rhs_uv[bulk_v_velocity] .+= τ .* vec(volumic_surface_tension_v)
        end
    end

   
    PDI_status = @ccall "libpdi".PDI_multi_expose("rhs_uv"::Cstring,
    "rhs_uv_len"::Cstring, length(rhs_uv)::Ref{Clonglong}, PDI_OUT::Cint,
    "rhs_uv_1D"::Cstring, rhs_uv::Ptr{Cdouble}, PDI_OUT::Cint,
    C_NULL::Ptr{Cvoid})::Cint

    #region cut-cell
    # @views kill_dead_cells!(rhs_uv[bulk_u_velocity], grid_u, geo_u[end])
    # @views kill_dead_cells!(rhs_uv[bulk_v_velocity], grid_v, geo_v[end])
    # _iLS = 1
    # for iLS in 1:nLS
    #     sbu = _iLS*niu+1:(_iLS+1)*niu
    #     sbv = ntu+_iLS*niv+1:ntu+(_iLS+1)*niv
    #     if !is_navier(bc_int[iLS]) && !is_navier_cl(bc_int[iLS])
    #         @views kill_dead_cells!(rhs_uv[sbu], grid_u, geo_u[end])
    #         @views kill_dead_cells!(rhs_uv[sbv], grid_v, geo_v[end])
    #         _iLS += 1
    #     end
    # end
    #endregion cut-cell

    if num.pressure_velocity_coupling == 3
        # uvD = ones(ntu + ntv + nNavier * nip + (num.nLS + 1) * nip + nbp)
        uvD = zeros(ntu + ntv + nNavier * nip + nip)
    elseif num.pressure_velocity_coupling == 4
        # uvD = ones(ntu + ntv + nNavier * nip + (num.nLS + 1) * nip + nbp)

        n_phase = 2
        # nip = grid_p.nx * grid_p.ny
        # nbp =  2 * grid_p.nx + 2 * grid_p.ny

        # niu = grid_u.nx * grid_u.ny
        # nbu = 2 * grid_u.nx + 2 * grid_u.ny
        # ntu = (nLS - nNavier + 1) * niu + nbu

        # niv = grid_v.nx * grid_v.ny
        # nbv = 2 * grid_v.nx + 2 * grid_v.ny
        # ntv = (nLS - nNavier + 1) * niv + nbv

        # ntNavier = num.nNavier * nip

        ntu1 = (nLS - nNavier + 1) * niu 
        ntv1 = (nLS - nNavier + 1) * niv 


        uvD = zeros(nphase * (ntu1 + ntv1 + nNavier * nip + (num.nLS + 1) * nip) + nbu + nbv)


    elseif num.pressure_velocity_coupling == 0
        uvD = ones(ntu + ntv + nNavier * nip)
    else
        # uvD = ones(ntu + ntv + nNavier * nip + (num.nLS + 1) * nip + nbp)
        uvD = zeros(ntu + ntv + nNavier * nip + (num.nLS + 1) * nip + nbp)
    end

    # print("\n before solving")
    PDI_status = @ccall "libpdi".PDI_multi_expose("rhs_uv"::Cstring,
    "rhs_uv_len"::Cstring, length(rhs_uv)::Ref{Clonglong}, PDI_OUT::Cint,
    "rhs_uv_1D"::Cstring, rhs_uv::Ptr{Cdouble}, PDI_OUT::Cint,
    C_NULL::Ptr{Cvoid})::Cint

    # print("\n diag Auv ", Auv.nzval)
    # display(Auv)
    # print("\n size Auv ",size(Auv))
    
    try
        @time uvD .= Auv \ rhs_uv
    catch e
        uvD .= Inf
        println(e)
    end


    vec1(ucorrD, grid_u) .= uvD[bulk_u_velocity]
    vecb(ucorrD, grid_u) .= uvD[border_u_velocity]
    #region cut-cell
    # kill_dead_cells!(vec1(ucorrD,grid_u), grid_u, geo_u[end])
    #endregion cut-cell
    
    ucorr .= reshape(vec1(ucorrD,grid_u), grid_u)

    vec1(vcorrD, grid_v) .= uvD[bulk_v_velocity]
    vecb(vcorrD, grid_v) .= uvD[border_v_velocity]
    #region cut-cell
    # kill_dead_cells!(vec1(vcorrD,grid_v), grid_v, geo_v[end])
    #endregion cut-cell
    
    vcorr .= reshape(vec1(vcorrD,grid_v), grid_v)

    #region cut-cell
    # nNav = 0
    # _iLS = 1
    # for iLS in 1:nLS
    #     if !is_navier(bc_int[iLS]) && !is_navier_cl(bc_int[iLS])
    #         veci(ucorrD,grid_u,iLS+1) .= uvD[_iLS*niu+1:(_iLS+1)*niu]
    #         kill_dead_cells!(veci(ucorrD,grid_u,iLS+1), grid_u, geo_u[end])

    #         veci(vcorrD,grid_v,iLS+1) .= uvD[ntu+_iLS*niv+1:ntu+(_iLS+1)*niv]
    #         kill_dead_cells!(veci(vcorrD,grid_v,iLS+1), grid_v, geo_v[end])
    #         _iLS += 1
    #     else
    #         @inbounds uT[nNav+1,:] .= vec(uvD[ntu+ntv+1+nNav*nip:ntu+ntv+(nNav+1)*nip])
    #         nNav += 1
    #     end
    # end
    #endregion cut-cell

    #endregion Navier


    PDI_status = @ccall "libpdi".PDI_multi_expose("print_velocity_prediction"::Cstring,
    "u_1D"::Cstring, ucorrD::Ptr{Cdouble}, PDI_OUT::Cint,
    "v_1D"::Cstring, vcorrD::Ptr{Cdouble}, PDI_OUT::Cint,
    "p_1D"::Cstring, ph.pD::Ptr{Cdouble}, PDI_OUT::Cint,
    C_NULL::Ptr{Cvoid})::Cint

    II = CartesianIndex(div(grid_v.ny,2),1)
    pII = lexicographic(II,grid_v.ny)
    # print("\n A coeff ",Av[pII,:])
    # print("\nM ",opC_p.iMy.diag[pII])

    # print("\nAv[pII,:]./mu1_over_rho1 ",Av[pII,:]./mu1_over_rho1)

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
    #region cut-cell
    # for iLS in 1:nLS
    #     if !is_navier(bc_int[iLS]) && !is_navier_cl(bc_int[iLS])
    #         Duv .+= opC_p.Gx[iLS] * veci(ucorrD,grid_u,iLS+1) .+ 
    #                 opC_p.Gy[iLS] * veci(vcorrD,grid_v,iLS+1)
    #     end
    # end
    #endregion cut-cell

    #TODO replace

    # Poisson equation: source term
    # divergence of velocity / dt
    vec1(rhs_ϕ,grid_p) .= iτ .* Duv

    #region needs to be corrected/documented for the signs, free surface pressure BC 
    
    # pres_free_suface = 0.0
    #TODO Marangoni
    #TODO phase change
    # diff_inv_rho = 1.0/rho1 - 1.0/rho2
    # jump_mass_flux = 0.0 #TODO

    #region cut-cell
    # if jump_mass_flux
    #     for iLS in 1:nLS
    #         if is_fs(bc_int[iLS])
    #             Smat = strain_rate(iLS, opC_u, opC_v, opC_p)
    #             S = Smat[1,1] * vec1(ucorrD,grid_u) .+ Smat[1,2] * veci(ucorrD,grid_u,iLS+1) .+
    #                 Smat[2,1] * vec1(vcorrD,grid_v) .+ Smat[2,2] * veci(vcorrD,grid_v,iLS+1)
    
    #             fs_mat = opC_p.HxT[iLS] * opC_p.Hx[iLS] .+ opC_p.HyT[iLS] * opC_p.Hy[iLS]
    #             veci(rhs_ϕ,grid_p,iLS+1) .= -2.0 .* mu1_over_rho1 .* S .+ Diagonal(diag(fs_mat)) * ( σ .* vec(grid_p.LS[iLS].κ) .- pres_free_suface .- diff_inv_rho * mass_flux ^ 2)
    #         end
    #     end
    # else
    #     for iLS in 1:nLS
    #         if is_fs(bc_int[iLS])
    #             Smat = strain_rate(iLS, opC_u, opC_v, opC_p)
    #             S = Smat[1,1] * vec1(ucorrD,grid_u) .+ Smat[1,2] * veci(ucorrD,grid_u,iLS+1) .+
    #                 Smat[2,1] * vec1(vcorrD,grid_v) .+ Smat[2,2] * veci(vcorrD,grid_v,iLS+1)

    #             fs_mat = opC_p.HxT[iLS] * opC_p.Hx[iLS] .+ opC_p.HyT[iLS] * opC_p.Hy[iLS]
    #             veci(rhs_ϕ,grid_p,iLS+1) .= -2.0 .* mu1_over_rho1 .* S .+ Diagonal(diag(fs_mat)) * ( σ .* vec(grid_p.LS[iLS].κ) .- pres_free_suface )
    #         end
    #     end
    # end
    #endregion cut-cell

    # Remove nullspace by adding small quantity to main diagonal
    if num.null_space == 0
        @inbounds @threads for i in 1:Aϕ.m
            @inbounds Aϕ[i,i] += 1e-10
        end
    end
    kill_dead_cells!(vec1(rhs_ϕ,grid_p), grid_p, geo[end])
    for iLS in 1:nLS
        kill_dead_cells!(veci(rhs_ϕ,grid_p,iLS+1), grid_p, geo[end])
    end
    # @time bicgstabl!(result_p, Aϕ, rhs_ϕ, Pl = Diagonal(Aϕ), log = true)

    #endregion needs to be corrected/documented for the signs, free surface pressure BC 


    # print("size pressure ",size(rhs_ϕ)," ",size(result_p)," ",size(Aϕ))
    # Solve Poisson equation
    # \phi^{n+1}: result_p
    # @time result_p .= Aϕ \ rhs_ϕ
    result_p = Aϕ \ rhs_ϕ

    

    #region cut-cell
    # kill_dead_cells!(vec1(result_p,grid_p), grid_p, geo[end])
    # for iLS in 1:nLS
    #     kill_dead_cells!(veci(result_p,grid_p,iLS+1), grid_p, geo[end])
    # end
    #endregion cut-cell

    ϕ .= reshape(vec1(result_p,grid_p), grid_p)

    iMu = Diagonal(inv_weight_eps2.(num.epsilon_mode,num.epsilon_vol,opC_u.M.diag))
    iMv = Diagonal(inv_weight_eps2.(num.epsilon_mode,num.epsilon_vol,opC_v.M.diag))
    # Gradient of pressure, eq. 17 in 
    #"A Conservative Cartesian Cut-Cell Method for Mixed Boundary Conditions and the Incompressible Navier-Stokes Equations on Staggered Meshes"
    ∇ϕ_x = opC_u.AxT * opC_u.Rx * vec(ϕ) .+ opC_u.Gx_b * vecb(result_p,grid_p)
    ∇ϕ_y = opC_v.AyT * opC_v.Ry * vec(ϕ) .+ opC_v.Gy_b * vecb(result_p,grid_p)
  
    #region cut-cell
    # for iLS in 1:nLS
    #     ∇ϕ_x .+= opC_u.Gx[iLS] * veci(result_p,grid_p,iLS+1)
    #     ∇ϕ_y .+= opC_v.Gy[iLS] * veci(result_p,grid_p,iLS+1)
    # end
    #endregion cut-cell

    # ∇ϕ_x = irho1 .* opC_u.AxT * opC_u.Rx * vec(ϕ) .+ opC_u.Gx_b * vecb(result_p,grid_p)
    # ∇ϕ_y = irho1 .* opC_v.AyT * opC_v.Ry * vec(ϕ) .+ opC_v.Gy_b * vecb(result_p,grid_p)
    # for iLS in 1:nLS
    #     ∇ϕ_x .+= irho1 .* opC_u.Gx[iLS] * veci(result_p,grid_p,iLS+1)
    #     ∇ϕ_y .+= irho1 .* opC_v.Gy[iLS] * veci(result_p,grid_p,iLS+1)
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
        vec1(pD,grid_p) .+= vec(ϕ) #no τ  since div u not rho1

    elseif num.prediction == "PmII" || num.prediction == "PmIIimposedpressure" || num.prediction == "PmIIimposedpressureBCincrement"
        # \nabla_h p^{n+1/2} = \nabla_h p^{n-1/2} + \nabla_h \phi^{n+1} - 
        #nu dt/2
        # \Delta t \nabla_h^2 \phi^{n+1} = \nabla_h \cdot \mathbf{u}^{*} \quad \text{in } \Omega
        

        print("\n max p vec1 ",maximum(vec1(pD,grid_p)))
        print("\n min p",minimum(ϕ)," max ",maximum(ϕ))
        print("\n max p",minimum(num.mu_cin1./2 .* reshape(iM * Duv,grid_p))," ",maximum(num.mu_cin1./2 .* reshape(iM * Duv,grid_p)))

        # todo zero neumann bc !!!!
        # vec1(pD,grid_p) .+= vec(ϕ .- num.mu_cin1./2 .* reshape(iM * Duv,grid_p))
        # or,
        # for better readability
        vec1(pD,grid_p) .= vec1(pD,grid_p) .+ vec(ϕ .- num.mu_cin1./2 .* reshape(iM * Duv,grid_p)) 
        # vec1(pD,grid_p) .+= vec(ϕ .- num.mu_cin1./2 .* reshape(iM * Duv,grid_p)) 

        print("\n max p vec1 ",maximum(vec1(pD,grid_p)))
        print("\n min p",minimum(ϕ)," max ",maximum(ϕ))
        print("\n max p",minimum(num.mu_cin1./2 .* reshape(iM * Duv,grid_p))," ",maximum(num.mu_cin1./2 .* reshape(iM * Duv,grid_p)))



    elseif num.prediction == "PmIII"
        # \Delta t \nabla_h^2 \phi^{n+1} = \nabla_h \cdot \mathbf{u}^{*} \quad \text{in } \Omega
        vec1(pD,grid_p) .= vec(ϕ .- num.mu_cin1./2 .* reshape(iM * Duv,grid_p)) #no contribution from p^{n-1/2}
    
    elseif num.prediction == "Flower" #occursin("Flower",num.prediction)
        # \nabla_h p^{n+1/2} = \nabla_h \phi^{n+1} # TODO: does not correspond to any formula in Brown 2001 ?
        vec1(pD,grid_p) .= vec(ϕ) #.- mu1_over_rho1 .* reshape(iM * Duv, grid_p))
    
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
        tmp_vec_p = zeros(grid_p) 
        # tmp_vec_p .= 0.0
        get_height!(grid_p.LS[1],grid_p.ind,grid_p.dx,grid_p.dy,grid_p.LS[end].geoS,tmp_vec_p) #here tmp_vec_p solid #TODO geoS ??

        init_fields_multiple_levelsets!(num,ph.pD,ph.p,tmp_vec_p,BC_p,grid_p,num.pres_intfc,"pL")

    else #update pressure at boundaries
        #region cut-cell
        # for iLS in 1:nLS
        #     veci(pD,grid_p,iLS+1) .= veci(result_p,grid_p,iLS+1)
        # end
        #endregion cut-cell
        vecb(pD,grid_p) .= vecb(result_p,grid_p) 
    end


    #endregion interfacial pressure is overwritten

    #TODO reapply Neumann BC for pressure boundary values ?


    p .= reshape(vec1(pD,grid_p), grid_p)

    #TODO
    # compute_grad_p!(num,grid_p, grid_u, grid_v, pD, opC_p, opC_u, opC_v)


    # else
    #     vec1(pD,grid_p) .= vec(p) .+ vec(ϕ) #.- mu1_over_rho1 .* iM * Duv
    #     vec2(pD,grid_p) .+= vec2(result_p,grid_p)
    #     vecb(pD,grid_p) .+= vecb(result_p,grid_p)
    #     p .= reshape(vec1(pD,grid_p), grid_p)
    # end

    # vec1(∇ϕ_x,grid_p) .*= irho1 
    # vec1(∇ϕ_y,grid_p) .*= irho1

    
    # u .= ucorr .- τ .* reshape(iMu * ∇ϕ_x, grid_u)
    # v .= vcorr .- τ .* reshape(iMv * ∇ϕ_y, grid_v)

    u .= ucorr .- τ .* reshape(iMu * ∇ϕ_x, grid_u) ./ rho_one_fluid_u
    v .= vcorr .- τ .* reshape(iMv * ∇ϕ_y, grid_v) ./ rho_one_fluid_v
    #region cut-cell
    # kill_dead_cells!(u, grid_u, geo_u[end])
    # kill_dead_cells!(v, grid_v, geo_v[end])
    #endregion cut-cell

    vec1(uD,grid_u) .= vec(u)
    vecb(uD,grid_u) .= vecb(ucorrD,grid_u)
    vec1(vD,grid_v) .= vec(v)
    vecb(vD,grid_v) .= vecb(vcorrD,grid_v)

    #region cut-cell
    # for iLS in 1:nLS
    #     if !is_navier(bc_int[iLS]) && !is_navier_cl(bc_int[iLS])
    #         veci(uD,grid_u,iLS+1) .= veci(ucorrD,grid_u,iLS+1)
    #         veci(vD,grid_v,iLS+1) .= veci(vcorrD,grid_v,iLS+1)
    #     end
    #     # if is_fs(bc_int[iLS])
    #     #     @inbounds for II in grid_u.ind.all_indices
    #     #         pII = lexicographic(II, grid_u.ny)
    #     #         if abs(veci(ucorrD,grid_u,iLS+1)[pII]) > 1e-12
    #     #             veci(ucorrD,grid_u,iLS+1)[pII] -= (τ .* iMu * ∇ϕ_x)[pII]
    #     #         end
    #     #     end
    #     #     @inbounds for II in grid_v.ind.all_indices
    #     #         pII = lexicographic(II, grid_v.ny)
    #     #         if abs(veci(vcorrD,grid_v,iLS+1)[pII]) > 1e-12
    #     #             veci(vcorrD,grid_v,iLS+1)[pII] -= (τ .* iMv * ∇ϕ_y)[pII]
    #     #         end
    #     #     end
    #     # end
    # end
    #endregion cut-cell


    grad_x = zeros(grid_u)
    grad_y = zeros(grid_v)

    #region cut-cell
    # compute_grad_T_x_T_y_array_u_v_capacities!(num, grid_p, grid_u, grid_v, opC_u, opC_v, grad_x, grad_y, ph.pD)
    #endregion cut-cell

    compute_grad_T_x_T_y_array_u_v_capacities!(num, grid_p, grid_u, grid_v, opC_u, opC_v, grad_x, grad_y, ph.pD)
    # compute_grad_T_x_T_y_array_u_v_capacities_one_fluid!(num, grid_p, grid_u, grid_v, opC_u, opC_v, grad_x, grad_y, ph.pD)

    # TODO give another set of capacities (allocate twice or special cases everywhere?) 

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

    #endregion correction 

    #end pressure velocity

    return Lp, bc_Lp, bc_Lp_b, Lu, diffusion_LS_u, diffusion_border_u, Lv, diffusion_LS_v, diffusion_border_v, opC_p.M, opC_u.M, opC_v.M, Cui, Cvi
end


"""
    set_Forward_Euler_one_fluid!(
        bc_int, num, grid_p, geo, grid_u, geo_u, grid_v, geo_v,
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
- `grid_p`: Grid structure.
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
- `Lu`, `diffusion_LS_u`, `diffusion_border_u`: Laplacian matrix and boundary conditions for the u-component.
- `Lv`, `diffusion_LS_v`, `diffusion_border_v`: Laplacian matrix and boundary conditions for the v-component.

### Description

1. **Advection Setup**: If advection is enabled, the convection terms are set up using the `set_convection!` function.
2. **Laplacian Matrices**: If advection is enabled, the Laplacian matrices are updated using the `set_matrices!` function. Otherwise, the Laplacian matrices are used.
3. **Right-Hand Side Vectors**:
   - For the Stokes system (`navier` is false), the right-hand side vectors for the u and v components are computed using the `FE_set_momentum` function.
   - For the Navier-Stokes system (`navier` is true), the right-hand side vector for the coupled system is computed using the `FE_set_momentum_coupled` function.
4. **Pressure Poisson Equation**: The right-hand side vector for the pressure Poisson equation is computed using the `set_poisson` function.

"""
function set_Forward_Euler_one_fluid!(
    bc_int, num, grid_p, geo, grid_u, geo_u, grid_v, geo_v,
    opC_p, opC_u, opC_v, BC_p, BC_u, BC_v,
    Au, Bu, Av, Bv, Aϕ, Auv, Buv,
    volume_fraction,rho_one_fluid_u,rho_one_fluid_v,
    Lpm1, bc_Lpm1, bc_Lpm1_b, Lum1, bc_Lum1, bc_Lum1_b, Lvm1, bc_Lvm1, bc_Lvm1_b,
    Mum1, Mvm1, op_conv, ph,
    periodic_x, periodic_y, advection, ls_advection, navier,rhs_uv = nothing)

    if advection
        set_convection!(num, grid_p, geo[end], grid_u, grid_u.LS, grid_v, grid_v.LS, ph.u, ph.v, op_conv, ph, BC_u, BC_v,opC_p, opC_u, opC_v)
    end

    if ls_advection
        update_all_ls_data(num, grid_p, grid_u, grid_v, bc_int, periodic_x, periodic_y, false)

        laps = set_matrices!(
            num, grid_p, geo, grid_u, geo_u, grid_v, geo_v,
            opC_p, opC_u, opC_v,
            periodic_x, periodic_y
        )

     

        # print("\n iMu",opC_u.M.diag)
        # print("\n iMv",opC_v.M.diag)

        # # display(reshape(vec1(iMu,grid_u),grid_u))
        # for j in 1:grid_u.ny
        #     for i in 1:grid_u.nx
        #         pII = lexicographic(CartesianIndex(j,i),grid_u.ny)
        #         print("\n iMu ",j," ",i," ",opC_u.M.diag[pII])
        #     end
        # end

        # if any(opC_u.M.diag  == 0 )
        #     @error("\n One-fluid model error M null")
        # end

    else
        laps = Lpm1, bc_Lpm1, bc_Lpm1_b, Lum1, bc_Lum1, bc_Lum1_b, Lvm1, bc_Lvm1, bc_Lvm1_b
    end

    Lp, bc_Lp, bc_Lp_b, Lu, diffusion_LS_u, diffusion_border_u, Lv, diffusion_LS_v, diffusion_border_v = laps



    # print("\n ls_advection ", ls_advection)

    # for j in 1:grid_u.ny
    #     for i in 1:grid_u.nx
    #         pII = lexicographic(CartesianIndex(j,i),grid_u.ny)
    #         print("\n iMu ",j," ",i," ",opC_u.M.diag[pII])
    #     end
    # end

    # if any(opC_u.M.diag  == 0 )
    #     @error("\n One-fluid model error M null")
    # end

    #region if num.one_fluid_model == 0
    # if num.one_fluid_model == 0
        

    #     diffusion_bulk_u = mu_over_rho.*Lu
    #     diffusion_LS_u = mu_over_rho.*diffusion_LS_u
    #     diffusion_border_u = mu_over_rho.*diffusion_border_u

    #     diffusion_bulk_v = mu_over_rho.*Lv
    #     diffusion_LS_v = mu_over_rho.*diffusion_LS_v
    #     diffusion_border_v = mu_over_rho.*diffusion_border_v

    # else

    # cf set_poisson_variable_coeff 


    #region Poisson variable coefficient
    
    # coeffD = viscosity 90 degrees ???


    # Laplacian: bulk 
    # L = BxT * iMx * Bx
    # L is of size nx*ny
    # tmp_x ((nx+1)*ny , nx*ny)
    # BxT (nx*ny, (nx+1)*ny )
    # Bx ((nx+1)*ny, nx*ny)
    # iMx ((nx+1)*ny, (nx+1)*ny )
    # iMx_b (nx+1)*ny,2*nx +2*ny
    # Hx_b  2*nx +2*ny,2*nx +2*ny

    # mat_coeffDx_b should be of size ((nx+1)*ny,(nx+1)*ny)

    if num.nLS>1
        printstyled(color=:red, @sprintf "\n TODO coeff poisson interpolation multiple levelsets:\n")
    end

    #region Viscosity coefficient for \frac{\partial u}{\partial x}
    #cf test in orientation.jl
    viscosity_coeff_for_du_dx = zeros(grid_u.ny,grid_u.nx+1)
    viscosity_coeff_for_du_dx[:,2:grid_u.nx] = volume_fraction #grid_p.LS[end].geoL.cap[:,:,5]

    #TODO contact angle change  viscosity_coeff_for_du_dx[:,1] and at end (not interpolating right now)
    viscosity_coeff_for_du_dx[:,1] = viscosity_coeff_for_du_dx[:,2]
    viscosity_coeff_for_du_dx[:,end] = viscosity_coeff_for_du_dx[:,end-1]

    # print("\n debug volume_fraction")
    # display(volume_fraction)

    # print("\n debug viscosity_coeff_for_du_dx")
    # display(viscosity_coeff_for_du_dx)


    # arithmetic average 
    # factor 2
    viscosity_coeff_for_du_dx .= 2 * ( (num.mu1 - num.mu2) * viscosity_coeff_for_du_dx  .+ num.mu2 )

    # display(viscosity_coeff_for_du_dx)


    # print("\n num.mu1 num.mu2 ",num.mu1 ," ",num.mu2)
    # display(viscosity_coeff_for_du_dx)
    PDI_status = @ccall "libpdi".PDI_multi_expose("viscosity_coeff_for_du_dx"::Cstring,
    "nstep"::Cstring, num.current_i ::Ref{Clonglong}, PDI_OUT::Cint,
    "viscosity_coeff_for_du_dx"::Cstring, viscosity_coeff_for_du_dx::Ptr{Cdouble}, PDI_OUT::Cint,        
    C_NULL::Ptr{Cvoid})::Cint

    diag_viscosity_coeff_for_du_dx = Diagonal(vec(viscosity_coeff_for_du_dx))
    #endregion Viscosity coefficient for \frac{\partial u}{\partial x}


    #region Viscosity coefficient for \frac{\partial u}{\partial y}


    viscosity_coeff_for_du_dy = zeros(grid_u.ny+1,grid_u.nx)
    viscosity_coeff_for_dv_dx = zeros(grid_v.ny,grid_v.nx+1)


    #region interpolate 
    volume_fraction_full = create_2D_grid_volume_fraction(grid_p,volume_fraction)
    grid_x_full_2D = create_2D_grid_x(grid_p,true,true)
    grid_y_full_2D = create_2D_grid_y(grid_p,true,true)

    all_grid_u_nodes_2D_x_for_du_dy_interp = create_2D_grid_x(grid_u,false,true)
    all_grid_u_nodes_2D_y_for_du_dy_interp = create_2D_grid_y(grid_u,false,true)

    # printstyled(color=:magenta, @sprintf "\n grid_x_full_2D \n") 

    # display(grid_x_full_2D)
    
    # printstyled(color=:magenta, @sprintf "\n grid_y_full_2D \n") 

    # display(grid_y_full_2D)

    # printstyled(color=:magenta, @sprintf "\n x for viscosity_coeff_for_du_dy \n") 

    # display(all_grid_u_nodes_2D_x_for_du_dy_interp)    

    # printstyled(color=:magenta, @sprintf "\n all_grid_u_nodes_2D_y_for_du_dy_interp \n") 

    # display(all_grid_u_nodes_2D_y_for_du_dy_interp)


    for j in 1:grid_u.ny+1
        for i in 1:grid_u.nx

            # print("\nvolume_fraction i ",i," j ",j,"\n")

            interp_coord_x = all_grid_u_nodes_2D_x_for_du_dy_interp[j,i] 

            interp_coord_y = (all_grid_u_nodes_2D_y_for_du_dy_interp[j,i] + all_grid_u_nodes_2D_y_for_du_dy_interp[j+1,i])/2 

            x1 = grid_x_full_2D[j,i]
            x2 = grid_x_full_2D[j,i+1]
            y1 = grid_y_full_2D[j,i]
            y2 = grid_y_full_2D[j+1,i+1]

            Q11 = volume_fraction_full[j,i]
            Q12 = volume_fraction_full[j+1,i]
            Q21 = volume_fraction_full[j,i+1]
            Q22 = volume_fraction_full[j+1,i+1]

            #     (x2,y2)
            #  Q12   Q22
            # | .    . |
            # |   x    | interpolate based on four nodes of value Q11, Q12, Q21, Q22
            # | .    . |
            #   Q11  Q21
            #(x1,y1)
            volume_fraction_face = bilinear_interpolation(interp_coord_x, interp_coord_y, x1, y1, x2, y2, Q11, Q12, Q21, Q22)

            # printstyled(color=:green, @sprintf "\n i %.5i j %.5i x %.2e y %.2e x1 %.2e y1 %.2e x2 %.2e y2 %.2e Q11 %.2e Q12 %.2e Q21 %.2e Q22 %.2e\n" i j interp_coord_x interp_coord_y x1 y1 x2 y2 Q11 Q12 Q21 Q22)
            # print("\nvolume_fraction i ",i," j ",j," ",volume_fraction_face)
            
            viscosity_coeff_for_du_dy[j,i] = harmonic_average_one_fluid(num.mu1,num.mu2,volume_fraction_face)

        end
    end
    

    PDI_status = @ccall "libpdi".PDI_multi_expose("viscosity_coeff_for_du_dy"::Cstring,
    "nstep"::Cstring, num.current_i ::Ref{Clonglong}, PDI_OUT::Cint,
    "viscosity_coeff_for_du_dy"::Cstring, viscosity_coeff_for_du_dy::Ptr{Cdouble}, PDI_OUT::Cint,        
    C_NULL::Ptr{Cvoid})::Cint

    # printstyled(color=:cyan, @sprintf "\n viscosity_coeff_for_du_dy min max\n" min_visc max_visc)
    # display(viscosity_coeff_for_du_dy)
    #endregion interpolate




    diag_viscosity_coeff_for_du_dy = Diagonal(vec(viscosity_coeff_for_du_dy))

    #endregion Viscosity coefficient for \frac{\partial u}{\partial x}


    mul!(opC_u.tmp_x, diag_viscosity_coeff_for_du_dx * opC_u.iMx, opC_u.Bx)
    diffusion_bulk_u = opC_u.BxT * opC_u.tmp_x

    # \frac{\partial}{\partial x} \left( \mu \frac{\partial u}{\partial x} \right)
    mul!(opC_u.tmp_y, diag_viscosity_coeff_for_du_dy * opC_u.iMy, opC_u.By)
    diffusion_bulk_u = diffusion_bulk_u .+ opC_u.ByT * opC_u.tmp_y

    # replaces Lu


    #region v diffusion

    #region Viscosity coefficient for \frac{\partial v}{\partial y}
    #cf test in orientation.jl
    viscosity_coeff_for_dv_dy = zeros(grid_v.ny+1,grid_v.nx)

    viscosity_coeff_for_dv_dy[2:grid_v.ny,:] = volume_fraction #grid_p.LS[end].geoL.cap[:,:,5]

    #TODO contact angle change  viscosity_coeff_for_dv_dy[:,1] and at end (not interpolating right now)
    viscosity_coeff_for_dv_dy[1,:] = viscosity_coeff_for_dv_dy[2,:]
    viscosity_coeff_for_dv_dy[end,:] = viscosity_coeff_for_dv_dy[end-1,:]




    # arithmetic average 
    viscosity_coeff_for_dv_dy .= 2 * ( (num.mu1 - num.mu2) * viscosity_coeff_for_dv_dy  .+ num.mu2 )

    # display(viscosity_coeff_for_dv_dy)
    PDI_status = @ccall "libpdi".PDI_multi_expose("viscosity_coeff_for_dv_dy"::Cstring,
    "nstep"::Cstring, num.current_i ::Ref{Clonglong}, PDI_OUT::Cint,
    "viscosity_coeff_for_dv_dy"::Cstring, viscosity_coeff_for_dv_dy::Ptr{Cdouble}, PDI_OUT::Cint,        
    C_NULL::Ptr{Cvoid})::Cint

    diag_viscosity_coeff_for_dv_dy = Diagonal(vec(viscosity_coeff_for_dv_dy))

    #endregion Viscosity coefficient for \frac{\partial v}{\partial y}


    #region interpolate 
    # volume_fraction_full = create_2D_grid_volume_fraction(grid_p,volume_fraction)
    # grid_x_full_2D = create_2D_grid_x(grid_p,true,true)
    # grid_y_full_2D = create_2D_grid_y(grid_p,true,true)

    # print("\n all_grid_v_nodes_2D_x_for_dv_dx_interp")
    all_grid_v_nodes_2D_x_for_dv_dx_interp = create_2D_grid_x(grid_v,true,false)
    # print("\n all_grid_v_nodes_2D_y_for_dv_dx_interp")
    all_grid_v_nodes_2D_y_for_dv_dx_interp = create_2D_grid_y(grid_v,true,false)

    for j in 1:grid_v.ny
        for i in 1:grid_v.nx+1

            # print("\nvolume_fraction i ",i," j ",j,"\n")

            interp_coord_x = (all_grid_v_nodes_2D_x_for_dv_dx_interp[j,i] + all_grid_v_nodes_2D_x_for_dv_dx_interp[j,i+1])/2 

            interp_coord_y = all_grid_v_nodes_2D_y_for_dv_dx_interp[j,i]

            x1 = grid_x_full_2D[j,i]
            x2 = grid_x_full_2D[j,i+1]
            y1 = grid_y_full_2D[j,i]
            y2 = grid_y_full_2D[j+1,i+1]

            Q11 = volume_fraction_full[j,i]
            Q12 = volume_fraction_full[j+1,i]
            Q21 = volume_fraction_full[j,i+1]
            Q22 = volume_fraction_full[j+1,i+1]

            #     (x2,y2)
            #  Q12   Q22
            # | .    . |
            # |   x    | interpolate based on four nodes of value Q11, Q12, Q21, Q22
            # | .    . |
            #   Q11  Q21
            #(x1,y1)
            volume_fraction_face = bilinear_interpolation(interp_coord_x, interp_coord_y, x1, y1, x2, y2, Q11, Q12, Q21, Q22)

            # printstyled(color=:green, @sprintf "\n i %.5i j %.5i x %.2e y %.2e x1 %.2e y1 %.2e x2 %.2e y2 %.2e Q11 %.2e Q12 %.2e Q21 %.2e Q22 %.2e\n" i j interp_coord_x interp_coord_y x1 y1 x2 y2 Q11 Q12 Q21 Q22)
            # print("\nvolume_fraction i ",i," j ",j," ",volume_fraction_face)
            
            viscosity_coeff_for_dv_dx[j,i] = harmonic_average_one_fluid(num.mu1,num.mu2,volume_fraction_face)

        end
    end
    

    PDI_status = @ccall "libpdi".PDI_multi_expose("viscosity_coeff_for_dv_dx"::Cstring,
    "nstep"::Cstring, num.current_i ::Ref{Clonglong}, PDI_OUT::Cint,
    "viscosity_coeff_for_dv_dx"::Cstring, viscosity_coeff_for_dv_dx::Ptr{Cdouble}, PDI_OUT::Cint,        
    C_NULL::Ptr{Cvoid})::Cint

    # printstyled(color=:cyan, @sprintf "\n viscosity_coeff_for_du_dy min max\n" min_visc max_visc)
    # display(viscosity_coeff_for_du_dy)
    #endregion interpolate

    diag_viscosity_coeff_for_dv_dx = Diagonal(vec(viscosity_coeff_for_dv_dx))

    
    mul!(opC_v.tmp_x, diag_viscosity_coeff_for_dv_dx * opC_v.iMx, opC_v.Bx)
    diffusion_bulk_v = opC_v.BxT * opC_v.tmp_x

    mul!(opC_v.tmp_y, diag_viscosity_coeff_for_dv_dy * opC_v.iMy, opC_v.By)
    diffusion_bulk_v = diffusion_bulk_v .+ opC_v.ByT * opC_v.tmp_y

    #endregion v diffusion


    # TODO: check localisation of fields
    
    #TODO wall control volume kappa for bulk vs bulk control volume kappa different  
    
 
    # print("\n opC_u.iMx_b", size(opC_u.iMx_b))

    # printstyled(color=:green, @sprintf "\n diffusion_border_u ")

    # print("\n diffusion_border_u")
    # display(diffusion_border_u)

    # viscosity_coeff_for_u_border_x = zeros(grid_p.ny, grid_p.nx+2) # copy(coeffDu)
    # viscosity_coeff_for_u_border_y = zeros(grid_p.ny+1, grid_p.nx+1) # copy(coeffDv)

    # viscosity_coeff_for_v_border_x = zeros(grid_p.ny+1, grid_p.nx+1) #copy(coeffDu)
    # viscosity_coeff_for_v_border_y = zeros(grid_p.ny+2, grid_p.nx) #copy(coeffDv)



    # @inbounds @threads for II in grid_p.ind.b_left[1] #[2:end-1]
    # # @inbounds @threads for II in grid_u.ind.b_left[1] #[2:end-1]
    #     coeffDu[II] = (vecb_L(coeffD,grid_p)[II[1]]+reshape(veci(coeffD,grid_p),grid_p)[II])/2.0 # veci(coeff,grid_p) or elec_cond
    #     # print("\n left, II ",II, coeffDu[II])
    # end

    # @inbounds @threads for II in grid_p.ind.b_right[1] #[2:end-1]
    #     coeffDu[ δx⁺(II) ] = (vecb_R(coeffD,grid_p)[II[1]]+reshape(veci(coeffD,grid_p),grid_p)[II])/2.0
    #     # print("\n right, II ",II, coeffDu[II])

    # end

    # @inbounds @threads for II in grid_p.ind.b_bottom[1] #[2:end-1]
    #     coeffDv[II] = (vecb_B(coeffD,grid_p)[II[2]]+reshape(veci(coeffD,grid_p),grid_p)[II])/2.0
    # end

    # @inbounds @threads for II in grid_p.ind.b_top[1] #[2:end-1]
    #     coeffDv[δy⁺(II)] = (vecb_T(coeffD,grid_p)[II[2]]+reshape(veci(coeffD,grid_p),grid_p)[II])/2.0
    # end

    # coeffD = fones()

    # # Interpolate at the border
    # interpolate_scalar_to_staggered_u_v_grids_at_border!(num,grid_u,coeffD,viscosity_coeff_for_u_border_x,viscosity_coeff_for_u_border_y)
    # interpolate_scalar_to_staggered_u_v_grids_at_border!(num,grid_v,coeffD,viscosity_coeff_for_v_border_x,viscosity_coeff_for_v_border_y)


    # viscosity_coeff_for_du_dx

    # TODO Here we assume that the viscosity is the same for border control volume associated to the gradient and for ...
    # Hypothesis viscosity_coeff_for_u_border_x = viscosity_coeff_for_du_dx
    
    # viscosity_coeff_for_u_border_x should be of size (ny, nx+2)
    # viscosity_coeff_for_u_border_y should be of size (ny+1, nx+1)

    # viscosity_coeff_for_v_border_x should be of size (ny+1, nx+1)
    # viscosity_coeff_for_v_border_y should be of size (ny+2, nx)

    # diag_viscosity_coeff_for_u_border_x = Diagonal(vec(viscosity_coeff_for_u_border_x)) 
    # diag_viscosity_coeff_for_u_border_y = Diagonal(vec(viscosity_coeff_for_u_border_y)) 

    # diag_viscosity_coeff_for_v_border_x = Diagonal(vec(viscosity_coeff_for_v_border_x)) 
    # diag_viscosity_coeff_for_v_border_y = Diagonal(vec(viscosity_coeff_for_v_border_y)) 

    # diffusion_border_u = (opC_u.BxT * diag_viscosity_coeff_for_u_border_x * opC_u.iMx_b * opC_u.Hx_b .+ opC_u.ByT * diag_viscosity_coeff_for_u_border_y * opC_u.iMy_b * opC_u.Hy_b)
    # diffusion_border_v = (opC_v.BxT * diag_viscosity_coeff_for_v_border_x * opC_v.iMx_b * opC_v.Hx_b .+ opC_v.ByT * diag_viscosity_coeff_for_v_border_y * opC_v.iMy_b * opC_v.Hy_b)

    diffusion_border_u = (opC_u.BxT * diag_viscosity_coeff_for_du_dx * opC_u.iMx_b * opC_u.Hx_b .+ 
                          opC_u.ByT * diag_viscosity_coeff_for_du_dy * opC_u.iMy_b * opC_u.Hy_b)
    diffusion_border_v = (opC_v.BxT * diag_viscosity_coeff_for_dv_dx * opC_v.iMx_b * opC_v.Hx_b .+ 
                          opC_v.ByT * diag_viscosity_coeff_for_dv_dy * opC_v.iMy_b * opC_v.Hy_b)


    # # cf  for Poisson
    #  coeffDu_border = copy(coeffDu)
    # coeffDv_border = copy(coeffDv)

    # # Interpolate conductivity at center of control volumes for potential gradient at the border
    # interpolate_scalar_to_staggered_u_v_grids_at_border!(num,grid_p,coeffD,coeffDu_border,coeffDv_border)

    # coeffDx_border = veci(coeffDu_border,grid_u)
    # coeffDy_border = veci(coeffDv_border,grid_v)
    # mat_coeffDx_b = Diagonal(vec(coeffDx_border)) 
    # mat_coeffDy_b = Diagonal(vec(coeffDy_border))

    # # bc_L_b = (BxT * mat_coeffDx_b * iMx_b * Hx_b .+ ByT * mat_coeffDy_b * iMy_b  * Hy_b)



    # diffusion_bulk_u   = mu_over_rho.*Lu
    # diffusion_LS_u     = 0.0 * iRe.*diffusion_LS_u
    # diffusion_border_u = iRe.*diffusion_border_u


    # diffusion_bulk_v   = mu_over_rho.*Lv
    # diffusion_LS_v     = 0.0 * iRe.*diffusion_LS_u
    # diffusion_border_u = iRe.*diffusion_border_u
    
    #endregion Poisson variable coefficient

    # end
    


    #region cross-terms
    
    # function laplacian(
    #mul!(tmp_x, iMx, Bx)
    # L = BxT * tmp_x
    # mul!(tmp_y, iMy, By)
    # L = L .+ ByT * tmp_y

  

    mul!(opC_v.tmp_x, diag_viscosity_coeff_for_dv_dx * opC_v.iMx, opC_v.Bx)
    # diffusion_bulk_v = opC_v.BxT * opC_v.tmp_x

    # Test
    # cross_term_diffusion_bulk_d_dv_dx_dy =

    print("\n size opC_v.iMx, opC_v.Bx",size(opC_v.tmp_x))
    print("\n size opC_u.ByT",size(opC_u.ByT))
    print("\n size opC_v.BxT",size(opC_v.BxT))


    cross_term_diffusion_bulk_d_dv_dx_dy = opC_u.ByT * opC_v.tmp_x

    # pII = lexicographic(CartesianIndex(div(grid_u.ny,2),div(grid_u.nx,2)),grid_u.ny)
    pII = lexicographic(CartesianIndex(5,5),grid_u.ny)

    print("\n cross_term_diffusion_bulk_d_dv_dx_dy ",pII)

    print("\n cross_term_diffusion_bulk_d_dv_dx_dy ",cross_term_diffusion_bulk_d_dv_dx_dy[pII,:])

    nip = grid_p.nx * grid_p.ny
    nbp = 2 * grid_p.nx + 2 * grid_p.ny

    niu = grid_u.nx * grid_u.ny
    nbu = 2 * grid_u.nx + 2 * grid_u.ny
    # ntu = (nLS - nNavier + 1) * niu + nbu
    ntu = niu + nbu

    niv = grid_v.nx * grid_v.ny
    nbv = 2 * grid_v.nx + 2 * grid_v.ny
    # ntv = (nLS - nNavier + 1) * niv + nbv
    ntv = niv + nbv

    bulk_u_velocity = 1:niu
    bulk_v_velocity =  ntu+1:ntu+niv
    # bulk_tangential_velocity = 

    border_u_velocity = ntu-nbu+1:ntu
    border_v_velocity = ntu+ntv-nbv+1:ntu+ntv
  

    print("\n indices ",bulk_u_velocity," ",bulk_v_velocity," ",border_u_velocity," ",border_v_velocity)

    print("\n indices ",size(bulk_u_velocity)," ",size(bulk_v_velocity)," ",size(border_u_velocity)," ",size(border_v_velocity))


    #TODO shift stencil

    mul!(opC_u.tmp_y, diag_viscosity_coeff_for_du_dy * opC_u.iMy, opC_u.By)
    # diffusion_bulk_u = diffusion_bulk_u .+ opC_u.ByT * opC_u.tmp_y
    cross_term_diffusion_bulk_d_du_dy_dx =  opC_v.BxT * opC_u.tmp_y


    # diffusion_border_u = (opC_u.BxT * diag_viscosity_coeff_for_du_dx * opC_u.iMx_b * opC_u.Hx_b 
    #         .+ opC_u.ByT * diag_viscosity_coeff_for_du_dy * opC_u.iMy_b * opC_u.Hy_b)

    # diffusion_border_v = (opC_v.BxT * diag_viscosity_coeff_for_dv_dx * opC_v.iMx_b * opC_v.Hx_b 
    #         .+ opC_v.ByT * diag_viscosity_coeff_for_dv_dy * opC_v.iMy_b * opC_v.Hy_b)


    cross_term_diffusion_bulk_d_dv_dx_dy_border = opC_u.ByT * diag_viscosity_coeff_for_dv_dx * opC_v.iMx_b * opC_v.Hx_b 
    
    cross_term_diffusion_bulk_d_du_dy_dx_border = opC_v.BxT * diag_viscosity_coeff_for_du_dy * opC_u.iMy_b * opC_u.Hy_b


    #endregion cross-terms

    #endregion if num.one_fluid_model == 0





    if num.pressure_velocity_coupling == 0

        if num.one_fluid_model == 1

            rhs_u = nothing
            rhs_v = nothing
            # rhs_ϕ = nothing
            diffusion_LS_u = nothing
            diffusion_LS_v = nothing
            rhs_uv = FE_set_momentum_coupled2_one_fluid(
            bc_int, num, grid_p, grid_u, grid_v,
            opC_p, opC_u, opC_v,
            Auv, Buv,
            rhs_uv,
            diffusion_bulk_u, diffusion_LS_u, diffusion_border_u, Mum1, BC_u,
            diffusion_bulk_v, diffusion_LS_v, diffusion_border_v, Mvm1, BC_v,
            cross_term_diffusion_bulk_d_dv_dx_dy,cross_term_diffusion_bulk_d_du_dy_dx,
            cross_term_diffusion_bulk_d_dv_dx_dy_border,cross_term_diffusion_bulk_d_du_dy_dx_border,rho_one_fluid_u,rho_one_fluid_v,
            ls_advection,BC_p,ph
        )

        

        elseif !navier
            rhs_u = FE_set_momentum(
                num, grid_u, opC_u,
                Au, Bu,
                diffusion_bulk_u, diffusion_LS_u, diffusion_border_u, Mum1, BC_u,
                ls_advection
            )
            rhs_v = FE_set_momentum(
                num, grid_v, opC_v,
                Av, Bv,
                diffusion_bulk_v, diffusion_LS_v, diffusion_border_v, Mvm1, BC_v,
                ls_advection
            )
            rhs_uv = nothing
        else
            rhs_u = nothing
            rhs_v = nothing
            rhs_uv = FE_set_momentum_coupled(
                bc_int, num, grid_p, grid_u, grid_v,
                opC_p, opC_u, opC_v,
                Auv, Buv,
                diffusion_bulk_u, diffusion_LS_u, diffusion_border_u, Mum1, BC_u,
                diffusion_bulk_v, diffusion_LS_v, diffusion_border_v, Mvm1, BC_v,
                ls_advection
            )
        end

        a0_p = []
        for i in 1:num.nLS
            push!(a0_p, zeros(grid_p))
        end
        # rhs_ϕ = set_poisson(
        #     bc_int, num, grid_p, a0_p, opC_p, opC_u, opC_v,
        #     Aϕ, Lp, bc_Lp, bc_Lp_b, BC_p,
        #     ls_advection
        # )

        # vecb(rhs,grid_p) .= +χ_b * vec(a0_b) #was - in set_poisson , -a0, a1 = -1,...

        # In solve_poisson, the equation a \frac{\partial p}{\partial n} + bp = g , 
        # if inhomogeneous Neumann: -1 at bottom and left
        # +1 sign at top and right
        
        rhs_ϕ = solve_poisson(
            bc_int, num, grid_p, a0_p, opC_p, opC_u, opC_v,
            Aϕ, Lp, bc_Lp, bc_Lp_b, BC_p,
            ls_advection
        )

        print("\n rhs phi ",size(rhs_ϕ))
        

    elseif num.pressure_velocity_coupling > 1
        rhs_u = nothing
        rhs_v = nothing
        rhs_ϕ = nothing
        rhs_uv = FE_set_momentum_coupled2(
            bc_int, num, grid_p, grid_u, grid_v,
            opC_p, opC_u, opC_v,
            Auv, Buv,
            rhs_uv,
            diffusion_bulk_u, diffusion_LS_u, diffusion_border_u, Mum1, BC_u,
            diffusion_bulk_v, diffusion_LS_v, diffusion_border_v, Mvm1, BC_v,
            ls_advection,BC_p,ph
        )

    end
    
    return rhs_u, rhs_v, rhs_ϕ, rhs_uv, Lp, bc_Lp, bc_Lp_b, Lu, diffusion_LS_u, diffusion_border_u, Lv, diffusion_LS_v, diffusion_border_v #TODO
end



# compute_grad_...
# mass_flux...
# -sigma kappa noru (noru grad c contains dirac)


"""
Set the system matrix for Forward-Euler scheme

* bc_interface: interface boundary condition
* B contains Mum1 and Mvm1 (cell volumes)

## Modifications to `rhs`

1. **For each interface `iLS`:**
   - If the boundary condition is not Navier and not Navier-CL:
     ```julia
     @inbounds rhs[sbu] .= opu.χ[iLS] * vec(a0u)
     @inbounds rhs[sbv] .= opv.χ[iLS] * vec(a0v)
     ```
   - Otherwise (if the boundary condition is Navier or Navier-CL):
     ```julia
     @inbounds rhs[ntu+ntv+1+nNav1*nip:ntu+ntv+(nNav1+1)*nip] .= opp.χ[iLS] * vec(a0p)
     nNav1 += 1
     ```

2. **For the outer boundaries/borders:**
   ```julia
   @inbounds rhs[border_u_velocity] .= opu.χ_b * vec(a0_bu) #BC for u component on borders
   @inbounds rhs[border_v_velocity] .= opv.χ_b * vec(a0_bv) #BC for v component on borders

robin BC : source term a0
                
At the moment, the Levelset is not computed at borders/interfaces ? Only bulk
"""
function FE_set_momentum_coupled2_one_fluid(
    bc_interface, num, grid_p, grid_u, grid_v,
    opp, opu, opv,
    A, B,
    rhs,
    diffusion_bulk_u, diffusion_LS_u, diffusion_border_u, Mum1, BCu,
    diffusion_bulk_v, diffusion_LS_v, diffusion_border_v, Mvm1, BCv,
    cross_term_diffusion_bulk_d_dv_dx_dy,cross_term_diffusion_bulk_d_du_dy_dx,
    cross_term_diffusion_bulk_d_dv_dx_dy_border,cross_term_diffusion_bulk_d_du_dy_dx_border,
    rho_one_fluid_u,rho_one_fluid_v,
    ls_advection::Bool,
    BCp,ph=nothing
    )
    @unpack τ, Re, nLS, nNavier = num


    printstyled(color=:red, @sprintf "\n coupled pressure-velocity FE_set_momentum_coupled2\n")


    #region init
    iRe = num.visc_coeff

    nip = grid_p.nx * grid_p.ny
    nbp = 2 * grid_p.nx + 2 * grid_p.ny

    niu = grid_u.nx * grid_u.ny
    nbu = 2 * grid_u.nx + 2 * grid_u.ny
    # ntu = (nLS - nNavier + 1) * niu + nbu
    ntu = niu + nbu

    niv = grid_v.nx * grid_v.ny
    nbv = 2 * grid_v.nx + 2 * grid_v.ny
    # ntv = (nLS - nNavier + 1) * niv + nbv
    ntv = niv + nbv

    #Reset to zero
    rhs .= 0.0 

    #region BC borders u
    a0_bu = zeros(nbu)
    _a1_bu = zeros(nbu)
    _b_bu = zeros(nbu)
    for iLS in 1:num.nLS
        set_borders!(grid_u, grid_u.LS[iLS].cl, grid_u.LS[iLS].u, a0_bu, _a1_bu, _b_bu, BCu, num.n_ext_cl)
    end
    a1_bu = Diagonal(vec(_a1_bu))
    b_bu = Diagonal(vec(_b_bu))
    #endregion BC borders u

    #region BC borders v
    a0_bv = zeros(nbv)
    _a1_bv = zeros(nbv)
    _b_bv = zeros(nbv)
    for iLS in 1:num.nLS
        set_borders!(grid_v, grid_v.LS[iLS].cl, grid_v.LS[iLS].u, a0_bv, _a1_bv, _b_bv, BCv, num.n_ext_cl)
    end
    a1_bv = Diagonal(vec(_a1_bv))
    b_bv = Diagonal(vec(_b_bv))
    #endregion BC borders v


    #region BC borders p
    if num.pressure_velocity_coupling == 2
        # rhs = fnzeros(grid_p, num)

        a0_bp = zeros(nbp)
        _a1_bp = zeros(nbp)
        _b_bp = zeros(nbp)
        for iLS in 1:num.nLS
            set_borders_poisson!(grid_p, grid_p.LS[iLS].cl, grid_p.LS[iLS].u, a0_bp, _a1_bp, _b_bp, BCp, num.n_ext_cl)
        end
        a1_bp = Diagonal(vec(_a1_bp))
        b_bp = Diagonal(vec(_b_bp))
    end

    # if ls_advection
    #     # Poisson equation
    #     # A[1:ni,1:ni] = pad(L, -4.0)
    #     # A[1:ni,end-nb+1:end] = bc_L_b

    #     # Boundary conditions for outer boundaries
    #     A[end-nb+1:end,1:ni] = b_bp * (HxT_b * iMx_b' * Bx .+ HyT_b * iMy_b' * By) 
    #     A[end-nb+1:end,end-nb+1:end] = pad(b_bp * (HxT_b * iMx_bd * Hx_b .+ HyT_b * iMy_bd * Hy_b) .+ χ_b * a1_bp,- 4.0)
    # end

    # for iLS in 1:num.nLS
    #     if ls_advection
    #         if is_dirichlet(bc_type[iLS])
    #             __a1 = 1.0
    #             __a2 = 0.0
    #             __b = 0.0
    #         elseif is_neumann(bc_type[iLS])
    #             __a1 = 0.0
    #             __a2 = 0.0
    #             __b = 1.0
    #         elseif is_robin(bc_type[iLS])
    #             __a1 = 1.0
    #             __a2 = 0.0
    #             __b = 1.0
    #         elseif is_fs(bc_type[iLS])
    #             __a1 = 0.0
    #             __a2 = 1.0
    #             __b = 0.0
    #         elseif is_wall_no_slip(bc_type[iLS])
    #             __a1 = 0.0
    #             __a2 = 0.0
    #             __b = 1.0
    #         elseif is_navier(bc_type[iLS])
    #             __a1 = 0.0
    #             __a2 = 0.0
    #             __b = 1.0
    #         elseif is_navier_cl(bc_type[iLS])
    #             __a1 = 0.0
    #             __a2 = 0.0
    #             __b = 1.0
    #         else
    #             __a1 = 0.0
    #             __a2 = 0.0
    #             __b = 1.0
    #         end
    
    #         _a1 = ones(grid_p) .* __a1
    #         a1 = Diagonal(vec(_a1))
    #         _a2 = ones(grid_p) .* __a2
    #         a2 = Diagonal(vec(_a2))
    #         _b = ones(grid_p) .* __b
    #         b = Diagonal(vec(_b))

    #         fs_mat = HxT[iLS] * Hx[iLS] .+ HyT[iLS] * Hy[iLS]

    #         sb = iLS*ni+1:(iLS+1)*ni
            
    #         # Poisson equation
    #         A[1:ni,sb] = bc_L[iLS]
    #         # Boundary conditions for inner boundaries
    #         A[sb,1:ni] = b * (HxT[iLS] * iMx * Bx .+ HyT[iLS] * iMy * By)
    #         # Contribution to Neumann BC from other boundaries
    #         for i in 1:num.nLS
    #             if i != iLS
    #                 A[sb,i*ni+1:(i+1)*ni] = b * (HxT[iLS] * iMx * Hx[i] .+ HyT[iLS] * iMy * Hy[i])
    #             end
    #         end
    #         A[sb,sb] = pad(
    #             b *(HxT[iLS] * iMx * Hx[iLS] .+ HyT[iLS] * iMy * Hy[iLS]) .+ χ[iLS] * a1 .-
    #             a2 * Diagonal(diag(fs_mat)), -4.0
    #         )
    #         #TODO was +b... here in set_poisson 
    #         A[sb,end-nb+1:end] = b * (HxT[iLS] * iMx_b * Hx_b .+ HyT[iLS] * iMy_b * Hy_b) #why was it different in set_poisson ?
    #         # Boundary conditions for outer boundaries
    #         A[end-nb+1:end,sb] = b_b * (HxT_b * iMx_b' * Hx[iLS] .+ HyT_b * iMy_b' * Hy[iLS])
    #     end

    #     veci(rhs,grid_p,iLS+1) .= +χ[iLS] * vec(a0[iLS]) #was - in set_poisson
    # end

    # vecb(rhs,grid_p) .= +χ_b * vec(a0_b) #was - in set_poisson
    #endregion BC borders p


    #endregion init

    # Array indices 
    bulk_u_velocity = 1:niu
    bulk_v_velocity =  ntu+1:ntu+niv
    # bulk_tangential_velocity = 

    border_u_velocity = ntu-nbu+1:ntu
    border_v_velocity = ntu+ntv-nbv+1:ntu+ntv

    # nt = (num.nLS - num.nNavier + 1) * ni_uv + num.nNavier * nip + nb_uv + (num.nLS + 1) * nip + nbp
    ntNavier = num.nNavier * nip
    # bulk_pressure = ntu+ntv+ntNavier+1:ntu+ntv+ntNavier+nip
    # #ntu+ntv+1:ntu+ntv+nip
    # border_pressure = ntu+ntv+ntNavier+(num.nLS + 1)*nip+1:ntu+ntv+ntNavier+(num.nLS + 1)*nip+nbp

    # print("\n len rhs_uv ",size(rhs))

    # print("bulk_pressure ",bulk_pressure)
    # print("border_pressure ",border_pressure)




    if ls_advection
        A.nzval .= 0.0

        pII = lexicographic(CartesianIndex(5,5),grid_u.ny)

        print("\n A[bulk_u_velocity,bulk_u_velocity] before ",pII)

        print("\n dt ",τ )
        print("\n diffusion_bulk_u ",diffusion_bulk_u[pII,:] )
        # example with mu1=mu2=1 and dt =1 : factor 2 for x, so 2-4 2 and 1 -2 1
        # diffusion_bulk_u   [101]  =  2.0
        # [132]  =  1.0
        # [133]  =  -6.0
        # [134]  =  1.0
        # [165]  =  2.0

        pII = lexicographic(CartesianIndex(1,5),grid_u.ny)
        pIIv = lexicographic(CartesianIndex(1,5),grid_v.ny)

        print("\n pIIv ",ntu + pIIv," pII ",pII)

        print("\n diffusion_bulk_u ",diffusion_bulk_u[pII,:] )

        pII = lexicographic(CartesianIndex(5,1),grid_u.ny)
        pIIv = lexicographic(CartesianIndex(5,1),grid_v.ny)

        print("\n pIIv ",ntu + pIIv," pII ",pII)

        print("\n diffusion_bulk_u ",diffusion_bulk_u[pII,:] )

        pII = lexicographic(CartesianIndex(1,1),grid_u.ny)
        pIIv = lexicographic(CartesianIndex(1,1),grid_v.ny)

        print("\n pIIv ",ntu + pIIv," pII ",pII)

        print("\n diffusion_bulk_u ",diffusion_bulk_u[pII,:] )


        # Implicit part of viscous term
        if num.non_dimensionalize == 0
            A[bulk_u_velocity,bulk_u_velocity] = pad_crank_nicolson(opu.M .- τ .* diffusion_bulk_u, grid_u, τ)
        else
            A[bulk_u_velocity,bulk_u_velocity] = pad_crank_nicolson(rho_one_fluid_u*opu.M .- τ .* diffusion_bulk_u, grid_u, τ)
        end

        # example 
        # A[bulk_u_velocity,bulk_u_velocity]   [101]  =  -2.0
        # [132]  =  -1.0
        # [133]  =  6.00391
        # [134]  =  -1.0
        # [165]  =  -2.0
        #same with opu.M


        print("\n A[bulk_u_velocity,bulk_u_velocity] ",A[pII,:])

        print("\n size(A) ",size(A))

        print("\n size(bulk_u_velocity) ",size(bulk_u_velocity))
        print("\n bulk_u_velocity ",bulk_u_velocity)
        print("\n bulk_v_velocity ",bulk_v_velocity)

        pIIv = lexicographic(CartesianIndex(5,5),grid_v.ny)

        # A[bulk_u_velocity,ntu + pIIv] .= 333
        
        print("\n pIIv ",ntu + pIIv)

        print("\n A[bulk_u_velocity,bulk_u_velocity] ",A[pII,:])
        
        print("\n grid_p ",grid_p.dx[5,5]," ", grid_u.dx[5,5] ," ",grid_v.dx[5,5]," ")
        print("\n grid_p ",grid_p.dy[5,5]," ", grid_u.dy[5,5] ," ",grid_v.dy[5,5]," ")


        A[bulk_u_velocity,bulk_v_velocity] = - τ * cross_term_diffusion_bulk_d_dv_dx_dy
        

        print("\n A[bulk_u_velocity,bulk_u_velocity] ",A[pII,:])


        # Contribution to implicit part of viscous term from outer boundaries
        A[bulk_u_velocity,border_u_velocity] = - τ .* diffusion_border_u 

        print("\n A[bulk_u_velocity,bulk_u_velocity] after border ",A[pII,:])

       
        A[bulk_u_velocity,border_v_velocity] = - τ .* cross_term_diffusion_bulk_d_dv_dx_dy_border

        # Boundary conditions for outer boundaries
        A[border_u_velocity,bulk_u_velocity] = b_bu * (opu.HxT_b * opu.iMx_b' * opu.Bx .+ opu.HyT_b * opu.iMy_b' * opu.By)
        A[border_u_velocity,border_u_velocity] = pad(b_bu * (
            opu.HxT_b * opu.iMx_bd * opu.Hx_b .+ 
            opu.HyT_b * opu.iMy_bd * opu.Hy_b
        ) .- opu.χ_b * a1_bu)

        # Implicit part of viscous term
        if num.non_dimensionalize == 0
            A[bulk_v_velocity,bulk_v_velocity] = pad_crank_nicolson(opv.M .- τ .* diffusion_bulk_v, grid_v, τ)
        else
            A[bulk_v_velocity,bulk_v_velocity] = pad_crank_nicolson(rho_one_fluid_v * opv.M .- τ .* diffusion_bulk_v, grid_v, τ)
        end

        A[bulk_v_velocity,bulk_u_velocity] = - τ .* cross_term_diffusion_bulk_d_du_dy_dx


        # Contribution to implicit part of viscous term from outer boundaries
        A[bulk_v_velocity,border_v_velocity] = - τ .* diffusion_border_v 
        A[bulk_v_velocity,border_u_velocity] = - τ .* cross_term_diffusion_bulk_d_du_dy_dx_border

        
        # Boundary conditions for outer boundaries
        A[border_v_velocity,bulk_v_velocity] = b_bv * (opv.HxT_b * opv.iMx_b' * opv.Bx .+ opv.HyT_b * opv.iMy_b' * opv.By)
        A[border_v_velocity,border_v_velocity] = pad(b_bv * (
            opv.HxT_b * opv.iMx_bd * opv.Hx_b .+ 
            opv.HyT_b * opv.iMy_bd * opv.Hy_b
        ) .- opv.χ_b * a1_bv)

        # TODO pad 1 or -4
        #TODO sign divergence not same u v and p

        pII = lexicographic(CartesianIndex(5,5),grid_u.ny)
        pIIv = lexicographic(CartesianIndex(5,5),grid_v.ny)

        print("\n pIIv ",ntu + pIIv," pII ",pII, " 5 5")

        print("\n A[bulk_u_velocity,bulk_u_velocity] ",A[pII,:])

        pII = lexicographic(CartesianIndex(1,5),grid_u.ny)
        pIIv = lexicographic(CartesianIndex(1,5),grid_v.ny)

        print("\n pIIv ",ntu + pIIv," pII ",pII, " 1 5 ")

        print("\n A[bulk_u_velocity,bulk_u_velocity] ",A[pII,:])

        pII = lexicographic(CartesianIndex(5,1),grid_u.ny)
        pIIv = lexicographic(CartesianIndex(5,1),grid_v.ny)

        print("\n pIIv ",ntu + pIIv," pII ",pII, " 5 1")

        print("\n A[bulk_u_velocity,bulk_u_velocity] ",A[pII,:])

        pII = lexicographic(CartesianIndex(1,1),grid_u.ny)
        pIIv = lexicographic(CartesianIndex(1,1),grid_v.ny)

        print("\n pIIv ",ntu + pIIv," pII ",pII)

        print("\n A[bulk_u_velocity,bulk_u_velocity] ",A[pII,:])

        #region coupled pression
        
        if num.pressure_velocity_coupling == 2

            # Poisson equation
            # A[1:ni,1:ni] = pad(L, -4.0)
            # A[1:ni,end-nb+1:end] = bc_L_b
    
            # Boundary conditions for outer boundaries
            A[border_pressure,bulk_pressure] = b_bp * (opp.HxT_b * opp.iMx_b' * opp.Bx .+ opp.HyT_b * opp.iMy_b' * opp.By) 
            A[border_pressure,border_pressure] = pad(b_bp * (
                opp.HxT_b * opp.iMx_bd * opp.Hx_b .+ 
                opp.HyT_b * opp.iMy_bd * opp.Hy_b) .+ opp.χ_b * a1_bp) #-4.0
            
        end
        
        #region Implicit gradient of pressure (volume integrated)

        # cf. Explicit gradient of pressure
        # ∇ϕ_x = opu.AxT * opu.Rx * vec1(pD,grid_p) .+ opu.Gx_b * vecb(pD,grid_p)
        # ∇ϕ_y = opv.AyT * opv.Ry * vec1(pD,grid_p) .+ opv.Gy_b * vecb(pD,grid_p)
        # for iLS in 1:nLS
        #     ∇ϕ_x .+= opu.Gx[iLS] * veci(pD,grid_p,iLS+1)
        #     ∇ϕ_y .+= opv.Gy[iLS] * veci(pD,grid_p,iLS+1)
        # end
        if num.pressure_velocity_coupling > 0


            irho1 = 1.0/num.rho1
            factor = num.τ * irho1
            A[bulk_u_velocity,bulk_pressure] = factor * opu.AxT * opu.Rx #TODO + or multiply by cell volume ? + not required, only component of matrix
            A[bulk_v_velocity,bulk_pressure] = factor * opv.AyT * opv.Ry 
        
            if num.pressure_velocity_coupling !=3
                #Outer boundaries
                A[bulk_u_velocity,border_pressure] = factor * opu.Gx_b
                A[bulk_v_velocity,border_pressure] = factor * opv.Gy_b
            

                for iLS in 1:nLS
                    interfacial_nb_iLS_pressure = ntu+ntv+nNavier*nip+nip*iLS+1:ntu+ntv+nNavier*nip+(iLS+1)*nip
                    A[bulk_u_velocity,interfacial_nb_iLS_pressure] .+= factor * opu.Gx[iLS] 
                    A[bulk_v_velocity,interfacial_nb_iLS_pressure] .+= factor * opv.Gy[iLS] 
                end
            
            end #num.pressure_velocity_coupling !=3
        end
        #endregion Implicit gradient of pressure (volume integrated)

        #region check coupled matrix
        # check_coupled_matrix()
        #endregion check coupled matrix

        #endregion coupled pression


        #TODO sbu sbv interfacial_nb_iLS_pressure

        #region divergence of velocity: -div U for symmetry
        if num.pressure_velocity_coupling > 0

            #equation written on lines ntu+ntv+1:ntu+ntv+nip 

            range_divergence = ntu+ntv+ntNavier+1:ntu+ntv+ntNavier+nip
            # u bulk
            A[range_divergence,bulk_u_velocity] = -opp.AxT 
            # u border
            A[range_divergence,border_u_velocity] = -opp.Gx_b
            # v bulk
            A[range_divergence,bulk_v_velocity] = -opp.AyT 
            # v border
            A[range_divergence,border_v_velocity] = -opp.Gy_b

            for iLS in 1:nLS
                sbu = iLS*niu+1:(iLS+1)*niu
                sbv = ntu+iLS*niv+1:ntu+(iLS+1)*niv
                # hypothesis no blowing so normal velocity null if Navier
                if !is_navier(bc_interface[iLS]) && !is_navier_cl(bc_interface[iLS]) 
                    A[range_divergence,sbu] = -opp.Gx[iLS] 
                    A[range_divergence,sbv] = -opp.Gy[iLS] 
                end
            end

            # divergence of velocity explicit
            # Duv = opp.AxT * vec1(ucorrD,grid_u) .+ opp.Gx_b * vecb(ucorrD,grid_u) .+
            #       opp.AyT * vec1(vcorrD,grid_v) .+ opp.Gy_b * vecb(vcorrD,grid_v)
            # for iLS in 1:nLS
            #     if !is_navier(bc_int[iLS]) && !is_navier_cl(bc_int[iLS])
            #         Duv .+= opp.Gx[iLS] * veci(ucorrD,grid_u,iLS+1) .+ 
            #                 opp.Gy[iLS] * veci(vcorrD,grid_v,iLS+1)
            #     end
            # end
        end
        #endregion divergence of velocity: -div U for symmetry
        print("size Mum1 " ,size(Mum1), " rho_one_fluid_u ", size(rho_one_fluid_u))

        if num.non_dimensionalize == 0
            B[bulk_u_velocity,bulk_u_velocity] = Mum1 #TODO rho_u ???
            B[bulk_v_velocity,bulk_v_velocity] = Mvm1
        elseif num.non_dimensionalize == 1 #TODO
            B[bulk_u_velocity,bulk_u_velocity] = rho_one_fluid_u * Mum1 #TODO rho_u ???
            B[bulk_v_velocity,bulk_v_velocity] = rho_one_fluid_v * Mvm1
        end
    end #ls_advection


    index_Navier = 0
    _iLS = 1
    for iLS in 1:num.nLS
        #TODO can be improved for readability / compilation:
        # allocates ones...
        # better tot do mapping function if cte... else if matrix ...

        range_nb_iLS_Navier = ntu+ntv+1+index_Navier*nip:ntu+ntv+(index_Navier+1)*nip

        #region BC iLS
        #for Navier: tangential velocity = \lambda grad tangential vel 
        a0u, a1u, bu, a0v, a1v, bv, a0p, bp = set_velocity_boundary_conditions(bc_interface, iLS, grid_u, grid_v, grid_p, num)
        #endregion BC iLS


        interfacial_nb_1_u_velocity = _iLS*niu+1:(_iLS+1)*niu
        interfacial_nb_1_v_velocity = ntu+_iLS*niv+1:ntu+(_iLS+1)*niv

        if ls_advection
            if !is_navier_cl(bc_interface[iLS]) && !is_navier(bc_interface[iLS])
                #region not Navier
                # Contribution to implicit part of viscous term from inner boundaries
                A[bulk_u_velocity,interfacial_nb_1_u_velocity] = - τ .* diffusion_LS_u[iLS]
                A[bulk_v_velocity,interfacial_nb_1_v_velocity] = - τ .* diffusion_LS_v[iLS] 
                # Boundary conditions for inner boundaries
                A[interfacial_nb_1_u_velocity,bulk_u_velocity] = bu * (opu.HxT[iLS] * opu.iMx * opu.Bx .+ opu.HyT[iLS] * opu.iMy * opu.By)
                A[interfacial_nb_1_v_velocity,bulk_v_velocity] = bv * (opv.HxT[iLS] * opv.iMx * opv.Bx .+ opv.HyT[iLS] * opv.iMy * opv.By)
                # Contribution to Neumann BC from other boundaries
                nNav2 = 0
                for i in 1:num.nLS
                    if i != iLS && (!is_navier_cl(bc_interface[i]) && !is_navier(bc_interface[i]))
                        A[interfacial_nb_1_u_velocity,i*niu+1:(i+1)*niu] = bu * (
                            opu.HxT[iLS] * opu.iMx * opu.Hx[i] .+
                            opu.HyT[iLS] * opu.iMy * opu.Hy[i]
                        )
                        A[interfacial_nb_1_v_velocity,ntu+i*niv+1:ntu+(i+1)*niv] = bv * (
                            opv.HxT[iLS] * opv.iMx * opv.Hx[i] .+
                            opv.HyT[iLS] * opv.iMy * opv.Hy[i]
                        )
                    elseif i != iLS
                        sinα = Diagonal(vec(sin.(grid_p.LS[i].α)))
                        # replace!(sinα.diag, NaN=>0.0)
                        cosα = Diagonal(vec(cos.(grid_p.LS[i].α)))
                        # replace!(cosα.diag, NaN=>0.0)

                        if any(isnan, sinα.diag) || any(isnan, cosα.diag)
                            @error("NaN FE_set_momentum_coupled")
                            replace!(sinα.diag, NaN=>0.0)
                            replace!(cosα.diag, NaN=>0.0)
                        end

                        # not Navier, averaging coefficients differently computed 
                        interpolate_x = interpolating_coefficient_Navier(grid_u,grid_p,i,bc_interface[iLS])
                        interpolate_y = interpolating_coefficient_Navier(grid_v,grid_p,i,bc_interface[iLS])

                        A[interfacial_nb_1_u_velocity,ntu+ntv+1+nNav2*nip:ntu+ntv+(nNav2+1)*nip] = bu * (
                            opu.HxT[iLS] * opu.iMx * opu.Hx[i] .+
                            opu.HyT[iLS] * opu.iMy * opu.Hy[i]
                        ) * interpolate_x * sinα
                        A[interfacial_nb_1_v_velocity,ntu+ntv+1+nNav2*nip:ntu+ntv+(nNav2+1)*nip] = bv * (
                            opv.HxT[iLS] * opv.iMx * opv.Hx[i] .+
                            opv.HyT[iLS] * opv.iMy * opv.Hy[i]
                        ) * interpolate_y * (-cosα)

                        nNav2 += 1
                    end
                end
                A[interfacial_nb_1_u_velocity,interfacial_nb_1_u_velocity] = pad(bu * (
                    opu.HxT[iLS] * opu.iMx * opu.Hx[iLS] .+
                    opu.HyT[iLS] * opu.iMy * opu.Hy[iLS]
                ) .- opu.χ[iLS] * a1u)
                A[interfacial_nb_1_v_velocity,interfacial_nb_1_v_velocity] = pad(bv * (
                    opv.HxT[iLS] * opv.iMx * opv.Hx[iLS] .+
                    opv.HyT[iLS] * opv.iMy * opv.Hy[iLS]
                ) .- opv.χ[iLS] * a1v)

                A[interfacial_nb_1_u_velocity,border_u_velocity] = bu * (
                    opu.HxT[iLS] * opu.iMx_b * opu.Hx_b .+ opu.HyT[iLS] * opu.iMy_b * opu.Hy_b
                )
                A[interfacial_nb_1_v_velocity,border_v_velocity] = bv * (
                    opv.HxT[iLS] * opv.iMx_b * opv.Hx_b .+ opv.HyT[iLS] * opv.iMy_b * opv.Hy_b
                )
                # Boundary conditions for outer boundaries
                A[border_u_velocity,interfacial_nb_1_u_velocity] = b_bu * (
                    opu.HxT_b * opu.iMx_b' * opu.Hx[iLS] .+ opu.HyT_b * opu.iMy_b' * opu.Hy[iLS]
                )
                A[border_v_velocity,interfacial_nb_1_v_velocity] = b_bv * (
                    opv.HxT_b * opv.iMx_b' * opv.Hx[iLS] .+ opv.HyT_b * opv.iMy_b' * opv.Hy[iLS]
                )

                _iLS += 1
                #endregion not Navier

            else 
                #region Navier
                # Tangential component of velocity if Navier BC #if !is_navier_cl(bc_interface[iLS]) && !is_navier(bc_interface[iLS])
                sinα_p = Diagonal(vec(sin.(grid_p.LS[iLS].α)))
                cosα_p = Diagonal(vec(cos.(grid_p.LS[iLS].α)))
                sinα_u = Diagonal(vec(sin.(grid_u.LS[iLS].α)))
                cosα_v = Diagonal(vec(cos.(grid_v.LS[iLS].α)))

                if any(isnan, sinα_p.diag) || any(isnan, cosα_p.diag) || any(isnan, sinα_u.diag) || any(isnan, cosα_v.diag)
                    @error("NaN FE_set_momentum_coupled")
                    replace!(sinα.diag, NaN=>0.0)
                    replace!(cosα.diag, NaN=>0.0)
                    replace!(sinα_u.diag, NaN=>0.0)
                    replace!(cosα_v.diag, NaN=>0.0)
                end

                # Contribution to implicit part of viscous term from inner boundaries
           
                # Navier BC, averaging coefficient computed differently
                interpolate_x = interpolating_coefficient_Navier(grid_u,grid_p,bc_interface[iLS])
                interpolate_y = interpolating_coefficient_Navier(grid_v,grid_p,bc_interface[iLS])
                # implicit part of viscous stress, from tangential velocity at wall with Navier, 
                # interpolated from scalar to u, v grids
                A[bulk_u_velocity,range_nb_iLS_Navier] = - iRe * τ .* (
                    opu.BxT * opu.iMx * opu.Hx[iLS] .+
                    opu.ByT * opu.iMy * opu.Hy[iLS]
                ) * sinα_u * interpolate_x
                A[bulk_v_velocity,range_nb_iLS_Navier] = - iRe * τ .* (
                    opv.BxT * opv.iMx * opv.Hx[iLS] .+
                    opv.ByT * opv.iMy * opv.Hy[iLS]
                ) * (-cosα_v) * interpolate_y

                # Boundary conditions for inner boundaries
                
                interpolate_u_to_p, interpolate_v_to_p = interpolating_coefficient_Navier_uv_grids_to_p_grid_volume(num,grid_p,grid_u,grid_v,iLS)

                A[range_nb_iLS_Navier,bulk_u_velocity] = bp * (
                    opp.HxT[iLS] * opp.iMx * opp.Bx .+
                    opp.HyT[iLS] * opp.iMy * opp.By
                ) * sinα_p * interpolate_u_to_p
                A[range_nb_iLS_Navier,bulk_v_velocity] = bp * (
                    opp.HxT[iLS] * opp.iMx * opp.Bx .+
                    opp.HyT[iLS] * opp.iMy * opp.By
                ) * (-cosα_p) * interpolate_v_to_p

                for i in 1:num.nLS
                    if i != iLS && (!is_navier_cl(bc_interface[i]) && !is_navier(bc_interface[i]))
                        #TODO why isn't it with volume like in interpolating_coefficient_Navier_uv_grids_to_p_grid_volume ?
                        interpolate_u_to_p, interpolate_v_to_p = interpolating_coefficient_Navier_uv_grids_to_p_grid_height(num,grid_p,grid_u,grid_v,i)

                        A[range_nb_iLS_Navier,i*niu+1:(i+1)*niu] = bp * (
                            opp.HxT[iLS] * opp.iMx * opp.Hx[i] .+
                            opp.HyT[iLS] * opp.iMy * opp.Hy[i]
                        ) * sinα_p * interpolate_u_to_p
                        A[range_nb_iLS_Navier,ntu+i*niv+1:ntu+(i+1)*niv] = bp * (
                            opp.HxT[iLS] * opp.iMx * opp.Hx[i] .+
                            opp.HyT[iLS] * opp.iMy * opp.Hy[i]
                        ) * (-cosα_p) * interpolate_v_to_p
                    end
                end
                A[range_nb_iLS_Navier,range_nb_iLS_Navier] = pad(bp * (
                    opp.HxT[iLS] * opp.iMx * opp.Hx[iLS] .+
                    opp.HyT[iLS] * opp.iMy * opp.Hy[iLS]
                ) .+ opp.χ[iLS])

                # coefficients (u,v) frame to (normal, tangent)
                sin_alpha_border = Diagonal(zeros(nbp))
                sin_alpha_border.diag[1:grid_p.ny] .= sin.(grid_p.LS[iLS].α[:,1])
                sin_alpha_border.diag[grid_p.ny+1:grid_p.ny+grid_p.nx] .= sin.(grid_p.LS[iLS].α[1,:])
                sin_alpha_border.diag[grid_p.ny+grid_p.nx+1:2gp.ny+grid_p.nx] .= sin.(grid_p.LS[iLS].α[:,end])
                sin_alpha_border.diag[2gp.ny+grid_p.nx+1:end] .= sin.(grid_p.LS[iLS].α[end,:])
                cos_alpha_border = Diagonal(zeros(nbp))
                cos_alpha_border.diag[1:grid_p.ny] .= cos.(grid_p.LS[iLS].α[:,1])
                cos_alpha_border.diag[grid_p.ny+1:grid_p.ny+grid_p.nx] .= cos.(grid_p.LS[iLS].α[1,:])
                cos_alpha_border.diag[grid_p.ny+grid_p.nx+1:2gp.ny+grid_p.nx] .= cos.(grid_p.LS[iLS].α[:,end])
                cos_alpha_border.diag[2gp.ny+grid_p.nx+1:end] .= cos.(grid_p.LS[iLS].α[end,:])

                #region interpolation coefficients
                interpolate_u_to_p = spdiagm(nbp, nbu, 0 => zeros(nbp), 1 => zeros(nbp-1))
                for ii in 1:grid_p.ny
                    interpolate_u_to_p[ii,ii] = 1.0
                    interpolate_u_to_p[grid_p.ny+grid_p.nx+ii,grid_u.ny+grid_u.nx+ii] = 1.0
                end
                for ii in 1:grid_p.nx
                    interpolate_u_to_p[ii+grid_p.ny,ii+grid_u.ny] = 0.5
                    interpolate_u_to_p[ii+grid_p.ny,ii+grid_u.ny+1] = 0.5
                    interpolate_u_to_p[ii+2gp.ny+grid_p.nx,ii+2gu.ny+grid_u.nx] = 0.5
                    interpolate_u_to_p[ii+2gp.ny+grid_p.nx,ii+2gu.ny+grid_u.nx+1] = 0.5
                end
                interpolate_v_to_p = spdiagm(nbp, nbv, 0 => zeros(nbp), 1 => zeros(nbp-1))
                for ii in 1:grid_p.ny
                    interpolate_v_to_p[ii,ii] = 0.5
                    interpolate_v_to_p[ii,ii+1] = 0.5
                    interpolate_v_to_p[grid_p.ny+grid_p.nx+ii,grid_v.ny+grid_v.nx+ii] = 0.5
                    interpolate_v_to_p[grid_p.ny+grid_p.nx+ii,grid_v.ny+grid_v.nx+ii+1] = 0.5
                end
                for ii in 1:grid_p.nx
                    interpolate_v_to_p[ii+grid_p.ny,ii+grid_v.ny] = 1.0
                    interpolate_v_to_p[ii+2gp.ny+grid_p.nx,ii+2gv.ny+grid_v.nx] = 1.0
                end
                #endregion interpolation coefficients

                #region Navier slip 
                A[range_nb_iLS_Navier,border_u_velocity] = bp * (
                    opp.HxT[iLS] * opp.iMx_b * opp.Hx_b .+ 
                    opp.HyT[iLS] * opp.iMy_b * opp.Hy_b
                ) * sin_alpha_border * interpolate_u_to_p
                A[range_nb_iLS_Navier,border_v_velocity] = bp * (
                    opp.HxT[iLS] * opp.iMx_b * opp.Hx_b .+ 
                    opp.HyT[iLS] * opp.iMy_b * opp.Hy_b
                ) * (-cos_alpha_border) * interpolate_v_to_p
                #endregion Navier slip
                
                # Boundary conditions for outer boundaries
                #region interpolation coefficients
                interpolate_u_to_p = spdiagm(nbu, nbp, 0 => zeros(nbp), 1 => zeros(nbp-1))
                for ii in 1:grid_u.ny
                    interpolate_u_to_p[ii,ii] = 1.0
                    interpolate_u_to_p[grid_u.ny+grid_u.nx+ii,grid_p.ny+grid_p.nx+ii] = 1.0
                end
                interpolate_u_to_p[grid_u.ny+1,grid_p.ny+1]  = 1.0
                interpolate_u_to_p[grid_u.ny+grid_u.nx,grid_p.ny+grid_p.nx]  = 1.0
                interpolate_u_to_p[2gu.ny+grid_u.nx+1,2gp.ny+grid_p.nx+1]  = 1.0
                interpolate_u_to_p[end,end]  = 1.0
                for ii in 2:(grid_u.nx-1)
                    interpolate_u_to_p[ii+grid_u.ny,ii+grid_p.ny-1] = 0.5
                    interpolate_u_to_p[ii+grid_u.ny,ii+grid_p.ny] = 0.5
                    interpolate_u_to_p[ii+2gu.ny+grid_u.nx,ii+2gp.ny+grid_p.nx-1] = 0.5
                    interpolate_u_to_p[ii+2gu.ny+grid_u.nx,ii+2gp.ny+grid_p.nx] = 0.5
                end
                interpolate_v_to_p = spdiagm(nbv, nbp, 0 => zeros(nbp), 1 => zeros(nbp-1))
                interpolate_v_to_p[1,1]  = 1.0
                interpolate_v_to_p[grid_v.ny,grid_p.ny]  = 1.0
                interpolate_v_to_p[grid_v.ny+grid_v.nx+1,grid_p.ny+grid_p.nx+1]  = 1.0
                interpolate_v_to_p[2gv.ny+grid_v.nx,2gp.ny+grid_p.nx]  = 1.0
                for ii in 2:(grid_v.ny-1)
                    interpolate_v_to_p[ii,ii-1] = 0.5
                    interpolate_v_to_p[ii,ii] = 0.5
                    interpolate_v_to_p[grid_v.ny+grid_v.nx+ii,grid_p.ny+grid_p.nx+ii] = 0.5
                    interpolate_v_to_p[grid_v.ny+grid_v.nx+ii,grid_p.ny+grid_p.nx+ii] = 0.5
                end
                for ii in 1:grid_v.nx
                    interpolate_v_to_p[ii+grid_v.ny,ii+grid_p.ny] = 1.0
                    interpolate_v_to_p[ii+2gv.ny+grid_v.nx,ii+2gp.ny+grid_p.nx] = 1.0
                end
                #endregion interpolation coefficients

                A[border_u_velocity,range_nb_iLS_Navier] = b_bu * (
                    opu.HxT_b * opu.iMx_b' * opu.Hx[iLS] .+ opu.HyT_b * opu.iMy_b' * opu.Hy[iLS]
                ) * interpolate_x * sinα_p
                A[border_v_velocity,range_nb_iLS_Navier] = b_bv * (
                    opv.HxT_b * opv.iMx_b' * opv.Hx[iLS] .+ opv.HyT_b * opv.iMy_b' * opv.Hy[iLS]
                ) * interpolate_y * (-cosα_p)
                
                #endregion Navier

            end # if !is_navier_cl(bc_interface[iLS]) && !is_navier(bc_interface[iLS])
        end #advection

        if !is_navier_cl(bc_interface[iLS]) && !is_navier(bc_interface[iLS])
            @inbounds rhs[interfacial_nb_1_u_velocity] .= opu.χ[iLS] * vec(a0u)
            @inbounds rhs[interfacial_nb_1_v_velocity] .= opv.χ[iLS] * vec(a0v)
        else
            @inbounds rhs[range_nb_iLS_Navier] .= opp.χ[iLS] * vec(a0p)
            index_Navier += 1
        end #!is_navier_cl(bc_interface[iLS]) && !is_navier(bc_interface[iLS])

    end #for iLS in 1:num.nLS

    
    @inbounds rhs[border_u_velocity] .= opu.χ_b * vec(a0_bu) #BC for u component on borders
    @inbounds rhs[border_v_velocity] .= opv.χ_b * vec(a0_bv) #BC for v component on borders

    if num.pressure_velocity_coupling == 2
        @inbounds rhs[border_pressure] .= opp.χ_b * vec(a0_bp) #BC for u component on borders
    end


    #region set first cells to boundary velocity
    if num.pressure_velocity_coupling == 3
        set_first_cells!(A,rhs,grid_u,ntu-nbu,0,true,false,true,false)
        set_first_cells!(A,rhs,grid_v,ntu+ntv-nbv,ntu,false,true,false,true)
    end
    #endregion set first cells to boundary velocity



    #rhs div u : 0 for Incompressible

    return rhs
end #end FE_set_momentum_coupled2


function check_coupled_matrix()

    #Test pressure part
    # Adummy = copy(A)
    # Adummy.nzval .= 0.0

    # nt = ntu+ntv + nNavier * nip + (num.nLS + 1) * nip + nbp
    # ncol_A = ntu+ntv + nNavier * nip + nip

    ni_p = grid_p.nx * grid_p.ny
    nb_p = 2 * grid_p.nx + 2 * grid_p.ny

    ni_u = grid_u.nx * grid_u.ny
    nb_u = 2 * grid_u.nx + 2 * grid_u.ny

    ni_v = grid_v.nx * grid_v.ny
    nb_v = 2 * grid_v.nx + 2 * grid_v.ny

    ni_uv = ni_u + ni_v
    nb_uv = nb_u + nb_v

    nt = (num.nLS - num.nNavier + 1) * ni_uv + num.nNavier * ni_p + nb_uv + (num.nLS + 1) * ni_p + nb_p
    ncol_A = (num.nLS - num.nNavier + 1) * ni_uv + num.nNavier * ni_p + nb_uv + ni_p

    Adummy = spzeros(ncol_A, nt)


    # length_colptr = length(Adummy.colptr)
    # length_rowval = length(Adummy.rowval)
    # length_nzval = length(Adummy.nzval)

    # PDI_status = @ccall "libpdi".PDI_multi_expose("print_matrix"::Cstring,
    # "Auv_n"::Cstring, Adummy.n::Ref{Clonglong}, PDI_OUT::Cint,
    # "Auv_m"::Cstring, Adummy.m::Ref{Clonglong}, PDI_OUT::Cint,
    # "Auv_colptr_len"::Cstring, length_colptr::Ref{Clonglong}, PDI_OUT::Cint,
    # "Auv_rowval_len"::Cstring, length_rowval::Ref{Clonglong}, PDI_OUT::Cint,
    # "Auv_nzval_len"::Cstring, length_nzval::Ref{Clonglong}, PDI_OUT::Cint,
    # # "Auv_colptr_len"::Cstring, length(Adummy.colptr)::Ref{Clonglong}, PDI_OUT::Cint,
    # # "Auv_rowval_len"::Cstring, length(Adummy.rowval)::Ref{Clonglong}, PDI_OUT::Cint,
    # # "Auv_nzval_len"::Cstring, length(Adummy.nzval)::Ref{Clonglong}, PDI_OUT::Cint,
    # "Auv_colptr_1D"::Cstring, Adummy.colptr::Ptr{Clonglong}, PDI_OUT::Cint,
    # "Auv_rowval_1D"::Cstring, Adummy.rowval::Ptr{Clonglong}, PDI_OUT::Cint,
    # "Auv_nzval_1D"::Cstring, Adummy.nzval::Ptr{Cdouble}, PDI_OUT::Cint,
    # C_NULL::Ptr{Cvoid})::Cint

    irho1 = 1.0/num.rho1
    factor = num.τ * irho1
    Adummy[bulk_u_velocity,bulk_pressure] = factor * opu.AxT * opu.Rx #TODO + or multiply by cell volume ? + not required, only component of matrix
    Adummy[bulk_v_velocity,bulk_pressure] = factor * opv.AyT * opv.Ry 
    #Outer boundaries
    Adummy[bulk_u_velocity,border_pressure] = factor * opu.Gx_b
    Adummy[bulk_v_velocity,border_pressure] = factor * opv.Gy_b


    for iLS in 1:nLS
        interfacial_nb_iLS_pressure = ntu+ntv+nNavier*nip+nip*iLS+1:ntu+ntv+nNavier*nip+(iLS+1)*nip
        Adummy[bulk_u_velocity,interfacial_nb_iLS_pressure] .+= factor * opu.Gx[iLS] 
        Adummy[bulk_v_velocity,interfacial_nb_iLS_pressure] .+= factor * opv.Gy[iLS] 
    end

    # print("\n A before ")

    # # print("\n  Auv.n ", A.n)
    # # print("\n  Auv.m ", A.m)

    # # print("\n  length(A.colptr) ", length(A.colptr))
    # # print("\n  length(A.rowval) ", length(A.rowval))
    # # print("\n  length(A.nzval) ", length(A.nzval))

    # # print("\n  Auv.colptr ", A.colptr)
    # # print("\n  Auv.rowval ", A.rowval)
    # # print("\n  Auv.nzval ", A.nzval)

    # length_colptr = length(A.colptr)
    # length_rowval = length(A.rowval)
    # length_nzval = length(A.nzval)


    # # nx: $nx
    # # ny: $ny
    # # nNavier: $nb_Navier_slip_BC
    # # nb_levelsets: $nb_levelsets

    # PDI_status = @ccall "libpdi".PDI_multi_expose("print_matrix_test"::Cstring,
    # "Auv_n"::Cstring, A.n::Ref{Clonglong}, PDI_OUT::Cint,
    # "Auv_m"::Cstring, A.m::Ref{Clonglong}, PDI_OUT::Cint,
    # "nx"::Cstring, grid_p.nx::Ref{Clonglong}, PDI_OUT::Cint,
    # "ny"::Cstring, grid_p.ny::Ref{Clonglong}, PDI_OUT::Cint,
    # "nb_Navier_slip_BC"::Cstring, num.nNavier::Ref{Clonglong}, PDI_OUT::Cint,
    # # "Auv_colptr_len"::Cstring, length_colptr::Ref{Clonglong}, PDI_OUT::Cint,
    # # "Auv_rowval_len"::Cstring, length_rowval::Ref{Clonglong}, PDI_OUT::Cint,
    # # "Auv_nzval_len"::Cstring, length_nzval::Ref{Clonglong}, PDI_OUT::Cint,
    # # # "Auv_colptr_len"::Cstring, length(A.colptr)::Ref{Clonglong}, PDI_OUT::Cint,
    # # # "Auv_rowval_len"::Cstring, length(A.rowval)::Ref{Clonglong}, PDI_OUT::Cint,
    # # # "Auv_nzval_len"::Cstring, length(A.nzval)::Ref{Clonglong}, PDI_OUT::Cint,
    # # "Auv_colptr_1D"::Cstring, A.colptr::Ptr{Clonglong}, PDI_OUT::Cint,
    # # "Auv_rowval_1D"::Cstring, A.rowval::Ptr{Clonglong}, PDI_OUT::Cint,
    # # "Auv_nzval_1D"::Cstring, A.nzval::Ptr{Cdouble}, PDI_OUT::Cint,
    # C_NULL::Ptr{Cvoid})::Cint

    # PDI_status = @ccall "libpdi".PDI_multi_expose("print_matrix"::Cstring,
    # "Auv_n"::Cstring, A.n::Ref{Clonglong}, PDI_OUT::Cint,
    # "Auv_m"::Cstring, A.m::Ref{Clonglong}, PDI_OUT::Cint,
    # "Auv_colptr_len"::Cstring, length_colptr::Ref{Clonglong}, PDI_OUT::Cint,
    # "Auv_rowval_len"::Cstring, length_rowval::Ref{Clonglong}, PDI_OUT::Cint,
    # "Auv_nzval_len"::Cstring, length_nzval::Ref{Clonglong}, PDI_OUT::Cint,
    # # "Auv_colptr_len"::Cstring, length(A.colptr)::Ref{Clonglong}, PDI_OUT::Cint,
    # # "Auv_rowval_len"::Cstring, length(A.rowval)::Ref{Clonglong}, PDI_OUT::Cint,
    # # "Auv_nzval_len"::Cstring, length(A.nzval)::Ref{Clonglong}, PDI_OUT::Cint,
    # "Auv_colptr_1D"::Cstring, A.colptr::Ptr{Clonglong}, PDI_OUT::Cint,
    # "Auv_rowval_1D"::Cstring, A.rowval::Ptr{Clonglong}, PDI_OUT::Cint,
    # "Auv_nzval_1D"::Cstring, A.nzval::Ptr{Cdouble}, PDI_OUT::Cint,
    # C_NULL::Ptr{Cvoid})::Cint

    # print("\n Adummy ",Adummy)
    print("\n Adummy ")

    # print("\n  Auv.n ", Adummy.n)
    # print("\n  Auv.m ", Adummy.m)

    # print("\n  length(Adummy.colptr) ", length(Adummy.colptr))
    # print("\n  length(Adummy.rowval) ", length(Adummy.rowval))
    # print("\n  length(Adummy.nzval) ", length(Adummy.nzval))

    # print("\n  Auv.colptr ", Adummy.colptr)
    # print("\n  Auv.rowval ", Adummy.rowval)
    # print("\n  Auv.nzval ", Adummy.nzval)

    
    #region print grad p part of matrix
    PDI_status = @ccall "libpdi".PDI_multi_expose("print_matrix"::Cstring,
    "Auv_n"::Cstring, Adummy.n::Ref{Clonglong}, PDI_OUT::Cint,
    "Auv_m"::Cstring, Adummy.m::Ref{Clonglong}, PDI_OUT::Cint,
    # "Auv_colptr_len"::Cstring, length_colptr::Ref{Clonglong}, PDI_OUT::Cint,
    # "Auv_rowval_len"::Cstring, length_rowval::Ref{Clonglong}, PDI_OUT::Cint,
    # "Auv_nzval_len"::Cstring, length_nzval::Ref{Clonglong}, PDI_OUT::Cint,
    "Auv_colptr_len"::Cstring, length(Adummy.colptr)::Ref{Clonglong}, PDI_OUT::Cint,
    "Auv_rowval_len"::Cstring, length(Adummy.rowval)::Ref{Clonglong}, PDI_OUT::Cint,
    "Auv_nzval_len"::Cstring, length(Adummy.nzval)::Ref{Clonglong}, PDI_OUT::Cint,
    "Auv_colptr_1D"::Cstring, Adummy.colptr::Ptr{Clonglong}, PDI_OUT::Cint,
    "Auv_rowval_1D"::Cstring, Adummy.rowval::Ptr{Clonglong}, PDI_OUT::Cint,
    "Auv_nzval_1D"::Cstring, Adummy.nzval::Ptr{Cdouble}, PDI_OUT::Cint,
    C_NULL::Ptr{Cvoid})::Cint
    #endregion print grad p part of matrix



    uvD = zeros(ntu + ntv + nNavier * nip + (num.nLS + 1) * nip + nbp)

    uvD[1:ntu] .= ph.uD
    uvD[ntu+1:ntu+ntv] .= ph.vD 
    uvD[ntu+ntv+ntNavier+1:ntu+ntv+ntNavier+(num.nLS+1)*nip+nbp] .= ph.pD

    PDI_status = @ccall "libpdi".PDI_multi_expose("rhs_uv"::Cstring,
    "rhs_uv_len"::Cstring, length(uvD)::Ref{Clonglong}, PDI_OUT::Cint,
    "rhs_uv_1D"::Cstring, uvD::Ptr{Cdouble}, PDI_OUT::Cint,
    C_NULL::Ptr{Cvoid})::Cint

    # print("\n test Adummy\n",Adummy*uvD)

    # print("\n test Adummy\n",Adummy*uvD/grid_p.dx[1,1]^2)

    # print("\n test grid_p ",grid_p.dx[1,1])
    # print("\n factor ", factor )

    control_volumes = ones(ncol_A)

    # control_volumes[1:ntu] .= vec(grid_u.LS[1].geoL.dcap[:,:,5])
    # control_volumes[[ntu+1:ntu+ntv]] .= vec(grid_v.LS[1].geoL.dcap[:,:,5])
    # control_volumes[ntu+ntv+ntNavier+1:ntu+ntv+ntNavier+(num.nLS+1)*nip+nbp] .= vec(grid_p.LS[1].geoL.dcap[:,:,5])
    
    control_volumes[1:ni_u] .= vec(grid_u.LS[1].geoL.dcap[:,:,5]) #non-dimensionalise bulk rhs

    control_volumes[ntu+1:ntu+ni_v] .= vec(grid_v.LS[1].geoL.dcap[:,:,5]) #non-dimensionalise bulk rhs
    # control_volumes[ntu+ntv+ntNavier+1:ntu+ntv+ntNavier+(num.nLS+1)*nip+nbp] .= vec(grid_p.LS[1].geoL.dcap[:,:,5])

    #niu 

    # print("\n test Adummy\n",Adummy*uvD/factor/grid_p.dx[1,1]^2)

    # print("\n test volume\n",vec(grid_u.LS[1].geoL.dcap[:,:,5]))

    # print("\n test Adummy\n",Adummy*uvD/factor ./ control_volumes)

    

    test_matrix = zeros(ntu + ntv + nNavier * nip + (num.nLS + 1) * nip + nbp)

    test_matrix = Adummy*uvD/factor ./ control_volumes

    PDI_status = @ccall "libpdi".PDI_multi_expose("check_coupled_matrix"::Cstring,
    "vec_1D_len"::Cstring, length(test_matrix)::Ref{Clonglong}, PDI_OUT::Cint,
    "vec_1D"::Cstring, test_matrix::Ptr{Cdouble}, PDI_OUT::Cint,
    C_NULL::Ptr{Cvoid})::Cint
end


"""
uses vector_convection

set BC
"""
function set_convection_with_rho!(
    num, grid_p, geo, grid_u, LS_u, grid_v, LS_v,
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

    if num.prediction == "PmIII" #cf Brown 2001, BC are to be corrected (intermediate step)
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
        compute_grad_T_x_T_y_array_u_v_capacities!(num, grid_p, grid_u, grid_v, opC_u, opC_v, grad_x, grad_y, ph.pD)
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


    #compare with:
    # scalar_convection!
    # inside : B[pII] += -0.5 * Dx[II] * ((A3 - B1) * Du[δx⁺(II)] + (B1 - A1) * Du[II])


    # ...
    # @inbounds O[pII,pII] += -0.25 * (A3_2 - B1_2) * Du[δx⁺(II)]
    # ...
    # @inbounds B[pII] += -0.25 * u[II] * (A3_2 - B1_2) * Du[δx⁺(II)]
    # ...

    # Compute convection Cu, CUTCu, Cv, CUTCv
    # vector_convection_with_rho!
    vector_convection!(dir, GridFCx, Cu, CUTCu, u, v, Du_x, Du_y, Dv_x, Dv_y,
            geo.dcap, grid_p.nx, grid_p.ny, BC_u, grid_u.ind.inside,
            grid_u.ind.b_left[1], grid_u.ind.b_bottom[1], grid_u.ind.b_right[1], grid_u.ind.b_top[1])
    vector_convection!(dir, GridFCy, Cv, CUTCv, u, v, Du_x, Du_y, Dv_x, Dv_y,
            geo.dcap, grid_p.nx, grid_p.ny, BC_v, grid_v.ind.inside,
            grid_v.ind.b_left[1], grid_v.ind.b_bottom[1], grid_v.ind.b_right[1], grid_v.ind.b_top[1])
    
    return nothing
end



"""
    fills O and B
"""
function vector_convection_with_rho!(::Dirichlet, ::Type{GridFCx}, O, B, u, v, Du_x, Du_y, Dv_x, Dv_y, cap, n, ny, BC, inside, b_left, b_bottom, b_right, b_top)
    B .= 0.0
    @inbounds @threads for II in inside
        fill_inside_conv!(GridFCx, O, B, u, v, Du_x, Dv_y, cap, ny, II)
    end

    @inbounds @threads for II in vcat(b_left, b_bottom[2:end-1], b_right, b_top[2:end-1])
        pII = lexicographic(II, ny)
        @inbounds O[pII,pII] = 0.0
    end
    bnds = (b_left, b_bottom[2:end-1], b_top[2:end-1])
    bc = ((Du_x, Dv_y), (Du_x, Dv_y), (Du_x, Dv_y))
    for (bnd, (Du, Dv)) in zip(bnds, bc)
        @inbounds @threads for II in bnd
            vec_convx_1!(II, O, B, u, Du, Dv, cap, ny)
        end
    end
    bnds = (b_bottom[2:end-1], b_right, b_top[2:end-1])
    bc = ((Du_x, Dv_y), (Du_x, Dv_y), (Du_x, Dv_y))
    for (bnd, (Du, Dv)) in zip(bnds, bc)
        @inbounds @threads for II in bnd
            vec_convx_2!(II, O, B, u, Du, Dv, cap, ny)
        end
    end
    @inbounds @threads for II in b_bottom[2:end-1]
        vec_convx_3!(II, O, v, cap, ny)
    end
    @inbounds @threads for II in b_top[2:end-1]
        vec_convx_4!(II, O, v, cap, ny)
    end
    @inbounds @threads for II in b_left[2:end]
        vec_convx_5!(II, O, v, cap, n, ny, BC)
    end
    @inbounds @threads for II in b_left[1:end-1]
        vec_convx_6!(II, O, v, cap, n, ny, BC)
    end
    @inbounds @threads for II in b_right[2:end]
        vec_convx_7!(II, O, v, cap, n, ny, BC)
    end
    @inbounds @threads for II in b_right[1:end-1]
        vec_convx_8!(II, O, v, cap, n, ny, BC)
    end

    if is_periodic(BC.left) && is_periodic(BC.right)
        @inbounds for (II, JJ) in zip(b_left, b_right)
            pII = lexicographic(II, ny)
            pJJ = lexicographic(JJ, ny)
            A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(JJ))
            
            Auim1, Aui = A1_1 * u[JJ], A3_1 * u[II]

            Au1 = 0.5 * (Auim1 + Aui)

            @inbounds O[pII,pII] += -0.5 * Au1
            @inbounds O[pII,pJJ] = -0.5 * Au1

            @inbounds O[pII,pII] += -0.25 * (A3_1 - B1_1) * Du_x[II]
            @inbounds O[pII,pII] += -0.25 * (B1_1 - A1_1) * Du_x[JJ]

            @inbounds O[pII,pII] += -0.25 * (A4_1 - B2_1) * Dv_y[δy⁺(δx⁻(JJ))]
            @inbounds O[pII,pII] += -0.25 * (B2_1 - A2_1) * Dv_y[δx⁻(JJ)]

            @inbounds B[pII] += -0.25 * u[II] * (A3_1 - B1_1) * Du_x[II]
            @inbounds B[pII] += -0.25 * u[II] * (B1_1 - A1_1) * Du_x[JJ]

            @inbounds B[pII] += -0.25 * u[II] * (A4_1 - B2_1) * Dv_y[δy⁺(δx⁻(JJ))]
            @inbounds B[pII] += -0.25 * u[II] * (B2_1 - A2_1) * Dv_y[δx⁻(JJ)]
        end
        @inbounds for (II, JJ) in zip(b_right, b_left)
            pII = lexicographic(II, ny)
            pJJ = lexicographic(JJ, ny)
            A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ)
            
            Aui, Auip1 = A1_2 * u[II], A3_2 * u[JJ]

            Au3 = 0.5 * (Aui + Auip1)

            @inbounds O[pII,pII] += 0.5 * Au3
            @inbounds O[pII,pJJ] = 0.5 * Au3

            @inbounds O[pII,pII] += -0.25 * (A3_2 - B1_2) * Du_x[JJ]
            @inbounds O[pII,pII] += -0.25 * (B1_2 - A1_2) * Du_x[II]

            @inbounds O[pII,pII] += -0.25 * (A4_2 - B2_2) * Dv_y[δy⁺(δx⁻(II))]
            @inbounds O[pII,pII] += -0.25 * (B2_2 - A2_2) * Dv_y[δx⁻(II)]

            @inbounds B[pII] += -0.25 * u[II] * (A3_2 - B1_2) * Du_x[JJ]
            @inbounds B[pII] += -0.25 * u[II] * (B1_2 - A1_2) * Du_x[II]

            @inbounds B[pII] += -0.25 * u[II] * (A4_2 - B2_2) * Dv_y[δy⁺(δx⁻(II))]
            @inbounds B[pII] += -0.25 * u[II] * (B2_2 - A2_2) * Dv_y[δx⁻(II)]
        end
    end
    if is_periodic(BC.bottom) && is_periodic(BC.top)
        @inbounds for (II,JJ) in zip(b_bottom[2:end-1], b_top[2:end-1])
            pII = lexicographic(II, ny)
            pJJ = lexicographic(JJ, ny)
            A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(II))
            A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
            
            Avim1jm1 = A2_1 * v[δx⁻(II)]
            Avip1jm1 = A2_2 * v[II]
    
            Au2 = 0.5 * (Avim1jm1 + Avip1jm1)
    
            @inbounds O[pII,pII] += -0.5 * Au2
            @inbounds O[pII,pJJ] = -0.5 * Au2
        end
        @inbounds for (II,JJ) in zip(b_top[2:end-1], b_bottom[2:end-1])
            pII = lexicographic(II, ny)
            pJJ = lexicographic(JJ, ny)
            A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(II))
            A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
            
            Avim1jp1 = A4_1 * v[δx⁻(δy⁺(JJ))]
            Avip1jp1 = A4_2 * v[δy⁺(JJ)]

            Au4 = 0.5 * (Avim1jp1 + Avip1jp1)
    
            @inbounds O[pII,pII] += 0.5 * Au4
            @inbounds O[pII,pJJ] = 0.5 * Au4
        end

        ii = b_left[1]
        pii = lexicographic(ii, ny)
        A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, ii)
        
        Avip1jm1 = A2_2 * v[ii]

        Au2 = 0.5 * Avip1jm1

        if is_periodic(BC.left)
            JJ = ii + CartesianIndex(0, n)
            A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, JJ)
            Avim1jm1 = A2_1 * v[δx⁻(JJ)]
            Au2 += 0.5 * Avim1jm1
        end

        JJ = ii + CartesianIndex(ny-1, 0)
        pJJ = lexicographic(JJ, ny)
        @inbounds O[pii,pii] += -0.5 * Au2
        @inbounds O[pii,pJJ] = -0.5 * Au2
        
        ii = b_left[end]
        pii = lexicographic(ii, ny)
        A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, ii)
        
        Avip1jp1 = A4_2 * v[δy⁺(ii)]

        Au4 = 0.5 * Avip1jp1

        if is_periodic(BC.left)
            JJ = ii + CartesianIndex(0, n)
            A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, JJ)
            Avim1jp1 = A4_1 * v[δy⁺(δx⁻(JJ))]
            Au4 += 0.5 * Avim1jp1
        end

        JJ = ii + CartesianIndex(-ny+1, 0)
        pJJ = lexicographic(JJ, ny)
        @inbounds O[pii,pii] += 0.5 * Au4
        @inbounds O[pii,pJJ] = 0.5 * Au4

        ii = b_right[1]
        pii = lexicographic(ii, ny)
        A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(ii))
        
        Avim1jm1 = A2_1 * v[δx⁻(ii)]

        Au2 = 0.5 * Avim1jm1

        if is_periodic(BC.right)
            JJ = ii + CartesianIndex(0, -n)
            A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ)
            Avip1jm1 = A2_2 * v[JJ]
            Au2 += 0.5 * Avip1jm1
        end

        JJ = ii + CartesianIndex(ny-1, 0)
        pJJ = lexicographic(JJ, ny)
        @inbounds O[pii,pii] += -0.5 * Au2
        @inbounds O[pii,pJJ] = -0.5 * Au2
        
        ii = b_right[end]
        pii = lexicographic(ii, ny)
        A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(ii))
        
        Avim1jp1 = A4_1 * v[δy⁺(δx⁻(ii))]

        Au4 = 0.5 * Avim1jp1

        if is_periodic(BC.right)
            JJ = ii + CartesianIndex(0, -n)
            A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ)
            Avip1jp1 = A4_2 * v[δy⁺(JJ)]
            Au4 += 0.5 * Avip1jp1
        end
        
        JJ = ii + CartesianIndex(-ny+1, 0)
        pJJ = lexicographic(JJ, ny)
        @inbounds O[pii,pii] += 0.5 * Au4
        @inbounds O[pii,pJJ] = 0.5 * Au4
    end

    return nothing
end



"""
### Variables

- `GridFCy`: v-grid_p.
- `O`: A matrix used to store intermediate results during the computation.
- `B`: A vector used to store boundary conditions.
- `u`: Velocity field in the x-direction.
- `v`: Velocity field in the y-direction.
- `Du_x`: 
- `Dv_y`: 
- `cap`: Capacities related to the convection terms.
- `ny`: Number of grid_p points in the y-direction.
- `inside`: Indices representing the interior grid_p points.
- `b_left`, `b_bottom`, `b_right`, `b_top`: Indices representing the boundary grid_p points.
- `BC`: Boundary conditions structure.


### Functions

- `fill_inside_conv!`: Updates the interior grid_p points based on convection terms.
- `vec_convy_1!`, `vec_convy_2!`, `vec_convy_3!`, `vec_convy_4!`, `vec_convy_5!`, `vec_convy_6!`, `vec_convy_7!`, `vec_convy_8!`: Functions to handle different boundary conditions and convection terms.
- `get_capacities_convection`: Retrieves capacities related to convection terms for a given grid_p point.
- `lexicographic`: Converts a multi-dimensional index to a linear index.

"""
function vector_convection_with_rho!(::Dirichlet, ::Type{GridFCy}, O, B, u, v, Du_x, Du_y, Dv_x, Dv_y, cap, n, ny, BC, inside, b_left, b_bottom, b_right, b_top)
    B .= 0.0
    @inbounds @threads for II in inside
        fill_inside_conv!(GridFCy, O, B, u, v, Du_x, Dv_y, cap, ny, II)
    end

    @inbounds @threads for II in vcat(b_left, b_bottom[2:end-1], b_right, b_top[2:end-1])
        pII = lexicographic(II, ny+1)
        @inbounds O[pII,pII] = 0.0
    end
    bnds = (b_left[2:end-1], b_bottom, b_right[2:end-1])
    bc = ((Du_x, Dv_y), (Du_x, Dv_y), (Du_x, Dv_y))
    for (bnd, (Du, Dv)) in zip(bnds, bc)
        @inbounds @threads for II in bnd
            vec_convy_1!(II, O, B, v, Du, Dv, cap, ny)
        end
    end
    bnds = (b_left[2:end-1], b_right[2:end-1], b_top)
    bc = ((Du_x, Dv_y), (Du_x, Dv_y), (Du_x, Dv_y))
    for (bnd, (Du, Dv)) in zip(bnds, bc)
        @inbounds @threads for II in bnd
            vec_convy_2!(II, O, B, v, Du, Dv, cap, ny)
        end
    end
    @inbounds @threads for II in b_left[2:end-1]
        vec_convy_3!(II, O, u, cap, ny)
    end
    @inbounds @threads for II in b_right[2:end-1]
        vec_convy_4!(II, O, u, cap, ny)
    end
    @inbounds @threads for II in b_bottom[2:end]
        vec_convy_5!(II, O, u, cap, ny, BC)
    end
    @inbounds @threads for II in b_bottom[1:end-1]
        vec_convy_6!(II, O, u, cap, ny, BC)
    end
    @inbounds @threads for II in b_top[2:end]
        vec_convy_7!(II, O, u, cap, ny, BC)
    end
    @inbounds @threads for II in b_top[1:end-1]
        vec_convy_8!(II, O, u, cap, ny, BC)
    end

    if is_periodic(BC.bottom) && is_periodic(BC.top)
        @inbounds for (II, JJ) in zip(b_bottom, b_top)
            pII = lexicographic(II, ny+1)
            pJJ = lexicographic(JJ, ny+1)
            A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(JJ))
            
            Avim1, Avi = A2_1 * v[JJ], A4_1 * v[II]

            Au2 = 0.5 * (Avim1 + Avi)

            @inbounds O[pII,pII] += -0.5 * Au2
            @inbounds O[pII,pJJ] = -0.5 * Au2

            @inbounds O[pII,pII] += -0.25 * (A4_1 - B2_1) * Dv_y[II]
            @inbounds O[pII,pII] += -0.25 * (B2_1 - A2_1) * Dv_y[JJ]

            @inbounds O[pII,pII] += -0.25 * (A3_1 - B1_1) * Du_x[δx⁺(δy⁻(JJ))]
            @inbounds O[pII,pII] += -0.25 * (B1_1 - A1_1) * Du_x[δy⁻(JJ)]

            @inbounds B[pII] += -0.25 * v[II] * (A4_1 - B2_1) * Dv_y[II]
            @inbounds B[pII] += -0.25 * v[II] * (B2_1 - A2_1) * Dv_y[JJ]

            @inbounds B[pII] += -0.25 * v[II] * (A3_1 - B1_1) * Du_x[δx⁺(δy⁻(JJ))]
            @inbounds B[pII] += -0.25 * v[II] * (B1_1 - A1_1) * Du_x[δy⁻(JJ)]
        end
        @inbounds for (II, JJ) in zip(b_top, b_bottom)
            pII = lexicographic(II, ny+1)
            pJJ = lexicographic(JJ, ny+1)
            A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ)
            
            Avi, Avip1 = A2_2 * v[II], A4_2 * v[JJ]

            Au4 = 0.5 * (Avi + Avip1)

            @inbounds O[pII,pII] += 0.5 * Au4
            @inbounds O[pII,pJJ] = 0.5 * Au4

            @inbounds O[pII,pII] += -0.25 * (A4_2 - B2_2) * Dv_y[JJ]
            @inbounds O[pII,pII] += -0.25 * (B2_2 - A2_2) * Dv_y[II]

            @inbounds O[pII,pII] += -0.25 * (A3_2 - B1_2) * Du_x[δx⁺(δy⁻(II))]
            @inbounds O[pII,pII] += -0.25 * (B1_2 - A1_2) * Du_x[δy⁻(II)]

            @inbounds B[pII] += -0.25 * v[II] * (A4_2 - B2_2) * Dv_y[JJ]
            @inbounds B[pII] += -0.25 * v[II] * (B2_2 - A2_2) * Dv_y[II]

            @inbounds B[pII] += -0.25 * v[II] * (A3_2 - B1_2) * Du_x[δx⁺(δy⁻(II))]
            @inbounds B[pII] += -0.25 * v[II] * (B1_2 - A1_2) * Du_x[δy⁻(II)]
        end
    end
    if is_periodic(BC.left) && is_periodic(BC.right)
        @inbounds for (II,JJ) in zip(b_left[2:end-1], b_right[2:end-1])
            pII = lexicographic(II, ny+1)
            pJJ = lexicographic(JJ, ny+1)
            A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(II))
            A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
            
            Auim1jm1 = A1_1 * u[δy⁻(II)]
            Auim1jp1 = A1_2 * u[II]
    
            Au1 = 0.5 * (Auim1jm1 + Auim1jp1)
    
            @inbounds O[pII,pII] += -0.5 * Au1
            @inbounds O[pII,pJJ] = -0.5 * Au1
        end
        @inbounds for (II,JJ) in zip(b_right[2:end-1], b_left[2:end-1])
            pII = lexicographic(II, ny+1)
            pJJ = lexicographic(JJ, ny+1)
            A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(II))
            A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
            
            Auip1jm1 = A3_1 * u[δx⁺(δy⁻(II))]
            Auip1jp1 = A3_2 * u[δx⁺(II)]
    
            Au3 = 0.5 * (Auip1jm1 + Auip1jp1)
    
            @inbounds O[pII,pII] += 0.5 * Au3
            @inbounds O[pII,pJJ] = 0.5 * Au3
        end

        ii = b_bottom[1]
        pii = lexicographic(ii, ny+1)
        A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, ii)
        
        Auim1jp1 = A1_2 * u[ii]

        Au1 = 0.5 * Auim1jp1

        if is_periodic(BC.bottom)
            JJ = ii + CartesianIndex(ny-1, 0)
            A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, JJ)
            Auim1jm1 = A1_1 * u[JJ]
            Au1 += 0.5 * Auim1jm1
        end

        JJ = ii + CartesianIndex(0, n-1)
        pJJ = lexicographic(JJ, ny+1)
        @inbounds O[pii,pii] += -0.5 * Au1
        @inbounds O[pii,pJJ] = -0.5 * Au1
        
        ii = b_bottom[end]
        pii = lexicographic(ii, ny+1)
        A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, ii)
        
        Auip1jp1 = A3_2 * u[δx⁺(ii)]

        Au3 = 0.5 * Auip1jp1

        if is_periodic(BC.bottom)
            JJ = ii + CartesianIndex(ny-1, 0)
            A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, JJ)
            Auip1jm1 = A3_1 * u[δx⁺(JJ)]
            Au3 += 0.5 * Auip1jm1
        end

        JJ = ii + CartesianIndex(0, -n+1)
        pJJ = lexicographic(JJ, ny+1)
        @inbounds O[pii,pii] += 0.5 * Au3
        @inbounds O[pii,pJJ] = 0.5 * Au3
        
        ii = b_top[1]
        pii = lexicographic(ii, ny+1)
        A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(ii))
        
        Auim1jm1 = A1_1 * u[δy⁻(ii)]

        Au1 = 0.5 * Auim1jm1

        if is_periodic(BC.top)
            JJ = ii + CartesianIndex(-ny, 0)
            A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ)
            Auim1jp1 = A1_2 * u[JJ]
            Au1 += 0.5 * Auim1jp1
        end

        JJ = ii + CartesianIndex(0, n-1)
        pJJ = lexicographic(JJ, ny+1)
        @inbounds O[pii,pii] += -0.5 * Au1
        @inbounds O[pii,pJJ] = -0.5 * Au1
        
        ii = b_top[end]
        pii = lexicographic(ii, ny+1)
        A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(ii))
        
        Auip1jm1 = A3_1 * u[δx⁺(δy⁻(ii))]

        Au3 = 0.5 * Auip1jm1

        if is_periodic(BC.top)
            JJ = ii + CartesianIndex(-ny, 0)
            A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ)
            Auip1jp1 = A3_2 * u[δx⁺(JJ)]
            Au3 += 0.5 * Auip1jp1
        end

        JJ = ii + CartesianIndex(0, -n+1)
        pJJ = lexicographic(JJ, ny+1)
        @inbounds O[pii,pii] += 0.5 * Au3
        @inbounds O[pii,pJJ] = 0.5 * Au3
    end

    return nothing
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
function set_poisson_one_fluid(
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
    # for iLS in 1:num.nLS
    #     set_borders!(grid, grid.LS[iLS].cl, grid.LS[iLS].u, a0_b, _a1_b, _b_b, BC, num.n_ext_cl)
    # end
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

    # for iLS in 1:num.nLS
    #     if ls_advection
    #         if is_dirichlet(bc_type[iLS])
    #             __a1 = -1.0
    #             __a2 = 0.0
    #             __b = 0.0
    #         elseif is_neumann(bc_type[iLS])
    #             __a1 = 0.0
    #             __a2 = 0.0
    #             __b = 1.0
    #         elseif is_robin(bc_type[iLS])
    #             __a1 = -1.0
    #             __a2 = 0.0
    #             __b = 1.0
    #         elseif is_fs(bc_type[iLS])
    #             __a1 = 0.0
    #             __a2 = 1.0
    #             __b = 0.0
    #         elseif is_wall_no_slip(bc_type[iLS])
    #             __a1 = 0.0
    #             __a2 = 0.0
    #             __b = 1.0
    #         elseif is_navier(bc_type[iLS])
    #             __a1 = 0.0
    #             __a2 = 0.0
    #             __b = 1.0
    #         elseif is_navier_cl(bc_type[iLS])
    #             __a1 = 0.0
    #             __a2 = 0.0
    #             __b = 1.0
    #         else
    #             __a1 = 0.0
    #             __a2 = 0.0
    #             __b = 1.0
    #         end
    
    #         _a1 = ones(grid) .* __a1
    #         a1 = Diagonal(vec(_a1))
    #         _a2 = ones(grid) .* __a2
    #         a2 = Diagonal(vec(_a2))
    #         _b = ones(grid) .* __b
    #         b = Diagonal(vec(_b))

    #         fs_mat = HxT[iLS] * Hx[iLS] .+ HyT[iLS] * Hy[iLS]

    #         sb = iLS*ni+1:(iLS+1)*ni
            
    #         # Poisson equation
    #         A[1:ni,sb] = bc_L[iLS]
    #         # Boundary conditions for inner boundaries
    #         A[sb,1:ni] = -b * (HxT[iLS] * iMx * Bx .+ HyT[iLS] * iMy * By)
    #         # Contribution to Neumann BC from other boundaries
    #         for i in 1:num.nLS
    #             if i != iLS
    #                 A[sb,i*ni+1:(i+1)*ni] = -b * (HxT[iLS] * iMx * Hx[i] .+ HyT[iLS] * iMy * Hy[i])
    #             end
    #         end
    #         A[sb,sb] = -pad(
    #             b * (HxT[iLS] * iMx * Hx[iLS] .+ HyT[iLS] * iMy * Hy[iLS]) .- χ[iLS] * a1 .+
    #             a2 * Diagonal(diag(fs_mat)), 4.0
    #         )
    #         A[sb,end-nb+1:end] = b * (HxT[iLS] * iMx_b * Hx_b .+ HyT[iLS] * iMy_b * Hy_b)
    #         # Boundary conditions for outer boundaries
    #         A[end-nb+1:end,sb] = -b_b * (HxT_b * iMx_b' * Hx[iLS] .+ HyT_b * iMy_b' * Hy[iLS])
    #     end

    #     veci(rhs,grid,iLS+1) .= -χ[iLS] * vec(a0[iLS])
    # end

    vecb(rhs,grid) .= -χ_b * vec(a0_b)
    
    return rhs
end