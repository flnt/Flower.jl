"""

"""
function update_one_fluid_density_viscosity(num,gp,gu,gv,volume_fraction,rho_one_fluid,mu_one_fluid,
    rho_one_fluid_u,mu_one_fluid_u,rho_one_fluid_v,mu_one_fluid_v)

    volume_fraction .= gp.LS[end].geoL.cap[:,:,5]

    if num.rho_one_fluid_average == 0 #arithmetic average
        rho_one_fluid .= num.rho1 * gp.LS[end].geoL.cap[:,:,5] .+ num.rho2 * (1.0 .- gp.LS[end].geoL.cap[:,:,5] )
        rho_one_fluid_u .= num.rho1 * gu.LS[end].geoL.cap[:,:,5] .+ num.rho2 * (1.0 .- gu.LS[end].geoL.cap[:,:,5] )
        rho_one_fluid_v .= num.rho1 * gv.LS[end].geoL.cap[:,:,5] .+ num.rho2 * (1.0 .- gv.LS[end].geoL.cap[:,:,5] )
    end

    if num.mu_one_fluid_average == 0 #arithmetic average
        mu_one_fluid  .= num.mu1  * gp.LS[end].geoL.cap[:,:,5]  .+ num.mu2  * (1.0 .- gp.LS[end].geoL.cap[:,:,5] )
        mu_one_fluid_u  .= num.mu1  * gu.LS[end].geoL.cap[:,:,5] .+ num.mu2  * (1.0 .- gu.LS[end].geoL.cap[:,:,5] )
        mu_one_fluid_v  .= num.mu1  * gv.LS[end].geoL.cap[:,:,5] .+ num.mu2  * (1.0 .- gv.LS[end].geoL.cap[:,:,5] )
    elseif num.mu_one_fluid_average == 1 #harmonic average
        mu_one_fluid   .= harmonic_average_one_fluid.(num.mu1,num.mu2,gp.LS[end].geoL.cap[:,:,5])
        mu_one_fluid_u .= harmonic_average_one_fluid.(num.mu1,num.mu2,gu.LS[end].geoL.cap[:,:,5])
        mu_one_fluid_v .= harmonic_average_one_fluid.(num.mu1,num.mu2,gv.LS[end].geoL.cap[:,:,5])
    end

    # PDI_status = @ccall "libpdi".PDI_multi_expose("print_one_fluid"::Cstring,
    # "rho_one_fluid"::Cstring, rho_one_fluid::Ptr{Cdouble}, PDI_OUT::Cint,
    # "mu_one_fluid"::Cstring, mu_one_fluid::Ptr{Cdouble}, PDI_OUT::Cint,
    # C_NULL::Ptr{Cvoid})::Cint

    PDI_status = @ccall "libpdi".PDI_multi_expose("write_one_fluid"::Cstring,
    "nstep"::Cstring, num.current_i ::Ref{Clonglong}, PDI_OUT::Cint,
    "rho_one_fluid"::Cstring, rho_one_fluid::Ptr{Cdouble}, PDI_OUT::Cint,
    "mu_one_fluid"::Cstring, mu_one_fluid::Ptr{Cdouble}, PDI_OUT::Cint,
    # "volume_fraction"::Cstring, volume_fraction::Ptr{Cdouble}, PDI_OUT::Cint,
    C_NULL::Ptr{Cvoid})::Cint

end

function harmonic_average_one_fluid(mu1,mu2, volume_fraction)
    return (mu1*mu2) / (mu2 * volume_fraction + (1.0 - volume_fraction) * mu1)
end

"""
store every nodes (border included) in a 2D matrix for interpolations
"""
function create_2D_grid_x(gp,add_x=true,add_y=true)

    if add_x
        nx = gp.nx+2
        rangex = 2:gp.nx+1
    else
        nx = gp.nx
        rangex = 1:gp.nx
    end

    if add_y
        ny = gp.ny+2
        rangey = 2:gp.ny+1
    else
        ny = gp.ny
        rangey= 1:gp.ny
    end

    
    create_2D_grid = zeros(ny,nx)
    

    x_centroid = gp.x .+ getproperty.(gp.LS[1].geoS.centroid, :x) .* gp.dx
    # y_centroid = gp.y .+ getproperty.(gp.LS[1].geoS.centroid, :y) .* gp.dy

    create_2D_grid[rangey,rangex] = x_centroid #grid.x

    x_bc_left = gp.x[:,1] .- gp.dx[:,1] ./ 2.0

    # y_bc_bottom = gp.y[1,:] .- gp.dy[1,:] ./ 2.0

    # y_bc_top = gp.y[end,:] .+ gp.dy[end,:] ./ 2.0

    x_bc_right = gp.x[:,end] .+ gp.dx[:,end] ./ 2.0

    # create_2D_grid[1,2:gp.nx] = create_2D_grid[2,2:gp.nx]

    # create_2D_grid[end,2:gp.nx] = create_2D_grid[end-1,2:gp.nx]

    # display(create_2D_grid)

    if add_x
        create_2D_grid[2:gp.ny+1,1] = x_bc_left
        create_2D_grid[2:gp.ny+1,end] = x_bc_right
    end

    if add_y
        create_2D_grid[1,:] = create_2D_grid[2,:]
        create_2D_grid[end,:] = create_2D_grid[end-1,:]
    end

    return create_2D_grid
end


"""
store every nodes (border included) in a 2D matrix for interpolations
"""
function create_2D_grid_y(gp,add_x=true,add_y=true)
    
    if add_x
        nx = gp.nx+2
        rangex = 2:gp.nx+1
    else
        nx = gp.nx
        rangex = 1:gp.nx
    end

    if add_y
        ny = gp.ny+2
        rangey = 2:gp.ny+1
    else
        ny = gp.ny
        rangey= 1:gp.ny
    end

    
    create_2D_grid = zeros(ny,nx)
    
    # x_centroid = gp.x .+ getproperty.(gp.LS[1].geoS.centroid, :x) .* gp.dx
    y_centroid = gp.y .+ getproperty.(gp.LS[1].geoS.centroid, :y) .* gp.dy

    create_2D_grid[rangey,rangex] = y_centroid #grid.x

    # x_bc_left = gp.x[:,1] .- gp.dx[:,1] ./ 2.0

    y_bc_bottom = gp.y[1,:] .- gp.dy[1,:] ./ 2.0

    y_bc_top = gp.y[end,:] .+ gp.dy[end,:] ./ 2.0

    # x_bc_right = gp.x[:,end] .+ gp.dx[:,end] ./ 2.0

    # create_2D_grid[1,2:gp.nx] = create_2D_grid[2,2:gp.nx]

    # create_2D_grid[end,2:gp.nx] = create_2D_grid[end-1,2:gp.nx]

    # display(create_2D_grid)

    # create_2D_grid[2:gp.ny+1,1] = x_bc_left

    # create_2D_grid[2:gp.ny+1,end] = x_bc_right

    # create_2D_grid[1,:] = create_2D_grid[2,:]

    # create_2D_grid[end,:] = create_2D_grid[end-1,:]

    # create_2D_grid[2:gp.ny+1,1] = create_2D_grid[2:gp.ny+1,2]
    # create_2D_grid[2:gp.ny+1,end] = create_2D_grid[2:gp.ny+1,end-1]

    # create_2D_grid[1,2:gp.nx+1] = y_bc_bottom

    # create_2D_grid[end,2:gp.nx+1] = y_bc_top

    # create_2D_grid[end,1] = create_2D_grid[end,2]
    # create_2D_grid[end,end] = create_2D_grid[end,end-1]


    if add_x
        create_2D_grid[2:gp.ny+1,1] = create_2D_grid[2:gp.ny+1,2]
        create_2D_grid[2:gp.ny+1,end] = create_2D_grid[2:gp.ny+1,end-1]
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
function create_2D_grid_volume_fraction(gp)
    #assuming no contact angle method
    create_2D_grid = zeros(gp.ny+2,gp.nx+2)

    # x_centroid = gp.x .+ getproperty.(gp.LS[1].geoS.centroid, :x) .* gp.dx
    # y_centroid = gp.y .+ getproperty.(gp.LS[1].geoS.centroid, :y) .* gp.dy

    create_2D_grid[2:gp.ny+1,2:gp.nx+1] = gp.LS[end].geoL.cap[:,:,5]


    # x_bc_left = gp.x[:,1] .- gp.dx[:,1] ./ 2.0

    # y_bc_bottom = gp.y[1,:] .- gp.dy[1,:] ./ 2.0

    # y_bc_top = gp.y[end,:] .+ gp.dy[end,:] ./ 2.0

    # x_bc_right = gp.x[:,end] .+ gp.dx[:,end] ./ 2.0

    # # create_2D_grid[1,2:gp.nx] = create_2D_grid[2,2:gp.nx]

    # # create_2D_grid[end,2:gp.nx] = create_2D_grid[end-1,2:gp.nx]

    # # display(create_2D_grid)

    # create_2D_grid[2:gp.ny+1,1] = x_bc_left

    # create_2D_grid[2:gp.ny+1,end] = x_bc_right

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
    # Q21 = vecb_R(volume_fraction_1D,grid)[j-1]
    # Q22 = vecb_R(volume_fraction_1D,grid)[j]


    return create_2D_grid
end


function compute_surface_tension(num,volumic_surface_tension_u,volumic_surface_tension_v)


    volume_fraction_1D = zeros(grid)

    vec1(volume_fraction_1D,grid) .= volume_fraction 

    compute_grad_T_x_T_y_array_u_v_capacities!(num, grid, grid_u, grid_v, opC_u, opC_v, normal_and_dirac_u, normal_and_dirac_v, volume_fraction_1D)

    # TODO border ???

    # divergence of velocity explicit
    curvature_p = opp.AxT * vec(normal_and_dirac_u) .+ opp.AyT * vec(normal_and_dirac_v) 
    
    compute_curvature_border = 0

    if compute_curvature_border >0
        # normal_and_dirac_u_1D = ...
        # normal_and_dirac_v_1D = ...
        curvature_p .+= opp.Gx_b * vecb(normal_and_dirac_u_1D,grid_u) .+ opp.Gy_b * vecb(normal_and_dirac_v_1D,grid_v)
    end

    #no interface: one-fluid model      
    # for iLS in 1:nLS
    #     if !is_navier(bc_int[iLS]) && !is_navier_cl(bc_int[iLS])
    #         curvature_p .+= opp.Gx[iLS] * veci(ucorrD,grid_u,iLS+1) .+ 
    #                 opp.Gy[iLS] * veci(vcorrD,grid_v,iLS+1)
    #     end
    # end


    # divergence of velocity explicit
    # Duv = opp.AxT * vec1(ucorrD,grid_u) .+ opp.Gx_b * vecb(ucorrD,grid_u) .+
    #       opp.AyT * vec1(vcorrD,grid_v) .+ opp.Gy_b * vecb(vcorrD,grid_v)
    # for iLS in 1:nLS
    #     if !is_navier(bc_int[iLS]) && !is_navier_cl(bc_int[iLS])
    #         Duv .+= opp.Gx[iLS] * veci(ucorrD,grid_u,iLS+1) .+ 
    #                 opp.Gy[iLS] * veci(vcorrD,grid_v,iLS+1)
    #     end
    # end

    interpolate_scalar!(grid, grid_u, grid_v, curvature, curvature_u, curvature_v)

    volumic_surface_tension_u = - num.sigma .* curvature_u .* normal_and_dirac_u
    volumic_surface_tension_v = - num.sigma .* curvature_v .* normal_and_dirac_v

    PDI_status = @ccall "libpdi".PDI_multi_expose("write_one_fluid_surface_tension"::Cstring,
    "nstep"::Cstring, num.current_i ::Ref{Clonglong}, PDI_OUT::Cint,
    "rho_one_fluid"::Cstring, rho_one_fluid::Ptr{Cdouble}, PDI_OUT::Cint,
    "mu_one_fluid"::Cstring, mu_one_fluid::Ptr{Cdouble}, PDI_OUT::Cint,
    "volume_fraction"::Cstring, volume_fraction::Ptr{Cdouble}, PDI_OUT::Cint,
    "curvature_p"::Cstring, curvature_p::Ptr{Cdouble}, PDI_OUT::Cint,
    "curvature_u"::Cstring, curvature_u::Ptr{Cdouble}, PDI_OUT::Cint,
    "curvature_v"::Cstring, curvature_v::Ptr{Cdouble}, PDI_OUT::Cint,
    "volumic_surface_tension_u"::Cstring, volumic_surface_tension_u::Ptr{Cdouble}, PDI_OUT::Cint,
    "volumic_surface_tension_v"::Cstring, volumic_surface_tension_v::Ptr{Cdouble}, PDI_OUT::Cint,
    C_NULL::Ptr{Cvoid})::Cint

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
function pressure_projection_one_fluid!(
    time_scheme, bc_int,
    num, grid, geo, grid_u, geo_u, grid_v, geo_v, ph,
    BC_u, BC_v, BC_p,
    opC_p, opC_u, opC_v, op_conv,
    Au, Bu, Av, Bv, Aϕ, Auv, Buv,
    Lpm1, bc_Lpm1, bc_Lpm1_b, Lum1, bc_Lum1, bc_Lum1_b, Lvm1, bc_Lvm1, bc_Lvm1_b,
    Cum1, Cvm1, Mum1, Mvm1,
    periodic_x, periodic_y, advection, ls_advection, current_i, Ra, navier,
    rho_one_fluid, mu_one_fluid,
    rho_one_fluid_u, mu_one_fluid_u,
    rho_one_fluid_v, mu_one_fluid_v,
    pres_free_suface,jump_mass_flux,mass_flux
    )
    @unpack Re, τ, σ, g, β, nLS, nNavier = num
    @unpack p, pD, ϕ, ϕD, u, v, ucorrD, vcorrD, uD, vD, ucorr, vcorr, uT = ph
    @unpack Cu, Cv, CUTCu, CUTCv = op_conv

    iτ = 1.0 / τ

    irho1 = 1.0 ./ rho_one_fluid
    mu1_over_rho1 = mu_one_fluid ./ rho_one_fluid
    
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

        ∇ϕ_x = opC_u.AxT * opC_u.Rx * vec1(pD,grid) .+ opC_u.Gx_b * vecb(pD,grid)
        ∇ϕ_y = opC_v.AyT * opC_v.Ry * vec1(pD,grid) .+ opC_v.Gy_b * vecb(pD,grid)
        #region cut-cell
        # for iLS in 1:nLS
        #     ∇ϕ_x .+= opC_u.Gx[iLS] * veci(pD,grid,iLS+1)
        #     ∇ϕ_y .+= opC_v.Gy[iLS] * veci(pD,grid,iLS+1)
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
        rhs_u, rhs_v, rhs_ϕ, rhs_uv, Lp, bc_Lp, bc_Lp_b, Lu, bc_Lu, bc_Lu_b, Lv, bc_Lv, bc_Lv_b = set_Forward_Euler_one_fluid!(
            bc_int, num, grid, geo, grid_u, geo_u, grid_v, geo_v,
            opC_p, opC_u, opC_v, BC_Poisson,BC_u, BC_v,
            Au, Bu, Av, Bv, Aϕ, Auv, Buv,
            Lpm1, bc_Lpm1, bc_Lpm1_b, Lum1, bc_Lum1, bc_Lum1_b, Lvm1, bc_Lvm1, bc_Lvm1_b,
            Mum1, Mvm1, mu1_over_rho1, op_conv, ph,
            periodic_x, periodic_y, advection, ls_advection, navier,rhs_uv = nothing,
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


    if num.one_fluid == 1 
        compute_surface_tension(num,)
    end
    
    # TODO PDI_multi_expose() #Cu u CUTCu

    # u and v are coupled if a Navier slip BC is employed inside, otherwise they are uncoupled
    #region not Navier
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
        vec1(rhs_u,grid_u) .-= τ .* ph.Gxm1 ./ rho_one_fluid_u

        #Surface tension
        vec1(rhs_u,grid_u) .+= τ .* volumic_surface_tension_u ./ rho_one_fluid_u

        
        
        #region cut-cell
        # kill_dead_cells!(vec1(rhs_u,grid_u), grid_u, geo_u[end])
        # for iLS in 1:nLS
        #     kill_dead_cells!(veci(rhs_u,grid_u,iLS+1), grid_u, geo_u[end])
        # end
        #endregion cut-cell

        # @time bicgstabl!(ucorrD, Au, rhs_u, log=true)
        try
            # @time bicgstabl!(ucorrD, Au, rhs_u, Pl=Diagonal(Au), log=true)
            @time ucorrD .= Au \ rhs_u
        catch e
            ucorrD .= Inf
            println(e)
        end

        #region cut-cell
        # kill_dead_cells!(vec1(ucorrD,grid_u), grid_u, geo_u[end])
        # for iLS in 1:nLS
        #     kill_dead_cells!(veci(ucorrD,grid_u,iLS+1), grid_u, geo_u[end])
        # end
        #endregion cut-cell

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

        vec1(rhs_v,grid_v) .-= τ .* ph.Gym1 ./ rho_one_fluid_v

        PDI_status = @ccall "libpdi".PDI_multi_expose("rhs_v"::Cstring,
        "v_1D"::Cstring, rhs_v::Ptr{Cdouble}, PDI_OUT::Cint,
        C_NULL::Ptr{Cvoid})::Cint

        #Surface tension
        vec1(rhs_v,grid_v) .+= τ .* volumic_surface_tension_v ./ rho_one_fluid_v

        #region cut-cell
        # kill_dead_cells!(vec1(rhs_v,grid_v), grid_v, geo_v[end])
        # for iLS in 1:nLS
        #     kill_dead_cells!(veci(rhs_v,grid_v,iLS+1), grid_v, geo_v[end])
        # end
        #endregion cut-cell

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

        #region cut-cell
        # kill_dead_cells!(vec1(vcorrD,grid_v), grid_v, geo_v[end])
        # for iLS in 1:nLS
        #     kill_dead_cells!(veci(vcorrD,grid_v,iLS+1), grid_v, geo_v[end])
        # end
        #endregion cut-cell

        vcorr .= reshape(vec1(vcorrD,grid_v), grid_v)
    #endregion not Navier

    #region Navier

    else #navier
        #fill u part at 1:niu
        uvm1 = zeros(ntu + ntv + nNavier * nip)
        uvm1[1:niu] .= vec1(uD,grid_u)
        uvm1[ntu+1:ntu+niv] .= vec1(vD,grid_v)
        uvm1[ntu-nbu+1:ntu] .= vecb(uD,grid_u)
        uvm1[ntu+ntv-nbv+1:ntu+ntv] .= vecb(vD,grid_v)

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


        rhs_uv .+=  Buv * uvm1

        rhs_uv[1:niu] .+= τ .* grav_x
        rhs_uv[1:niu] .-= τ .* Convu
        rhs_uv[1:niu] .+= τ .* ra_x
        rhs_uv[1:niu] .-= τ .* ph.Gxm1 ./ rho_one_fluid_u

        rhs_uv[ntu+1:ntu+niv] .+= τ .* grav_y
        rhs_uv[ntu+1:ntu+niv] .-= τ .* Convv
        rhs_uv[ntu+1:ntu+niv] .+= τ .* ra_y
        rhs_uv[ntu+1:ntu+niv] .-= τ .* ph.Gym1 ./ rho_one_fluid_v

        #region cut-cell
        # @views kill_dead_cells!(rhs_uv[1:niu], grid_u, geo_u[end])
        # @views kill_dead_cells!(rhs_uv[ntu+1:ntu+niv], grid_v, geo_v[end])
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

        uvD = ones(ntu + ntv + nNavier * nip)
        try
            @time uvD .= Auv \ rhs_uv
        catch e
            uvD .= Inf
            println(e)
        end

        vec1(ucorrD, grid_u) .= uvD[1:niu]
        vecb(ucorrD, grid_u) .= uvD[ntu-nbu+1:ntu]
        #region cut-cell
        # kill_dead_cells!(vec1(ucorrD,grid_u), grid_u, geo_u[end])
        #endregion cut-cell
        
        ucorr .= reshape(vec1(ucorrD,grid_u), grid_u)

        vec1(vcorrD, grid_v) .= uvD[ntu+1:ntu+niv]
        vecb(vcorrD, grid_v) .= uvD[ntu+ntv-nbv+1:ntu+ntv]
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

    end #navier
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
    vec1(rhs_ϕ,grid) .= iτ .* Duv

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
    #             veci(rhs_ϕ,grid,iLS+1) .= -2.0 .* mu1_over_rho1 .* S .+ Diagonal(diag(fs_mat)) * ( σ .* vec(grid.LS[iLS].κ) .- pres_free_suface .- diff_inv_rho * mass_flux ^ 2)
    #         end
    #     end
    # else
    #     for iLS in 1:nLS
    #         if is_fs(bc_int[iLS])
    #             Smat = strain_rate(iLS, opC_u, opC_v, opC_p)
    #             S = Smat[1,1] * vec1(ucorrD,grid_u) .+ Smat[1,2] * veci(ucorrD,grid_u,iLS+1) .+
    #                 Smat[2,1] * vec1(vcorrD,grid_v) .+ Smat[2,2] * veci(vcorrD,grid_v,iLS+1)

    #             fs_mat = opC_p.HxT[iLS] * opC_p.Hx[iLS] .+ opC_p.HyT[iLS] * opC_p.Hy[iLS]
    #             veci(rhs_ϕ,grid,iLS+1) .= -2.0 .* mu1_over_rho1 .* S .+ Diagonal(diag(fs_mat)) * ( σ .* vec(grid.LS[iLS].κ) .- pres_free_suface )
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
    kill_dead_cells!(vec1(rhs_ϕ,grid), grid, geo[end])
    for iLS in 1:nLS
        kill_dead_cells!(veci(rhs_ϕ,grid,iLS+1), grid, geo[end])
    end
    # @time bicgstabl!(ϕD, Aϕ, rhs_ϕ, Pl = Diagonal(Aϕ), log = true)

    #endregion needs to be corrected/documented for the signs, free surface pressure BC 

    # Solve Poisson equation
    # \phi^{n+1}: ϕD
    @time ϕD .= Aϕ \ rhs_ϕ

    #region cut-cell
    # kill_dead_cells!(vec1(ϕD,grid), grid, geo[end])
    # for iLS in 1:nLS
    #     kill_dead_cells!(veci(ϕD,grid,iLS+1), grid, geo[end])
    # end
    #endregion cut-cell

    ϕ .= reshape(vec1(ϕD,grid), grid)

    iMu = Diagonal(inv_weight_eps2.(num.epsilon_mode,num.epsilon_vol,opC_u.M.diag))
    iMv = Diagonal(inv_weight_eps2.(num.epsilon_mode,num.epsilon_vol,opC_v.M.diag))
    # Gradient of pressure, eq. 17 in 
    #"A Conservative Cartesian Cut-Cell Method for Mixed Boundary Conditions and the Incompressible Navier-Stokes Equations on Staggered Meshes"
    ∇ϕ_x = opC_u.AxT * opC_u.Rx * vec(ϕ) .+ opC_u.Gx_b * vecb(ϕD,grid)
    ∇ϕ_y = opC_v.AyT * opC_v.Ry * vec(ϕ) .+ opC_v.Gy_b * vecb(ϕD,grid)
  
    #region cut-cell
    # for iLS in 1:nLS
    #     ∇ϕ_x .+= opC_u.Gx[iLS] * veci(ϕD,grid,iLS+1)
    #     ∇ϕ_y .+= opC_v.Gy[iLS] * veci(ϕD,grid,iLS+1)
    # end
    #endregion cut-cell

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
        #region cut-cell
        # for iLS in 1:nLS
        #     veci(pD,grid,iLS+1) .= veci(ϕD,grid,iLS+1)
        # end
        #endregion cut-cell
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
    # compute_grad_T_x_T_y_array_u_v_capacities!(num, grid, grid_u, grid_v, opC_u, opC_v, grad_x, grad_y, ph.pD)
    #endregion cut-cell

    compute_grad_T_x_T_y_array_u_v_capacities_one_fluid!(num, grid, grid_u, grid_v, opC_u, opC_v, grad_x, grad_y, ph.pD)

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


    return Lp, bc_Lp, bc_Lp_b, Lu, bc_Lu, bc_Lu_b, Lv, bc_Lv, bc_Lv_b, opC_p.M, opC_u.M, opC_v.M, Cui, Cvi
end


"""
    set_Forward_Euler_one_fluid!(
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
function set_Forward_Euler_one_fluid!(
    bc_int, num, grid, geo, grid_u, geo_u, grid_v, geo_v,
    opC_p, opC_u, opC_v, BC_p, BC_u, BC_v,
    Au, Bu, Av, Bv, Aϕ, Auv, Buv,
    Lpm1, bc_Lpm1, bc_Lpm1_b, Lum1, bc_Lum1, bc_Lum1_b, Lvm1, bc_Lvm1, bc_Lvm1_b,
    Mum1, Mvm1, mu_over_rho, op_conv, ph,
    periodic_x, periodic_y, advection, ls_advection, navier,rhs_uv = nothing,
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


    if num.one_fluid == 0
        
        Lp, bc_Lp, bc_Lp_b, Lu, bc_Lu, bc_Lu_b, Lv, bc_Lv, bc_Lv_b = laps

        diffusion_bulk_u = mu_over_rho.*Lu
        diffusion_LS_u = mu_over_rho.*bc_Lu
        diffusion_border_u = mu_over_rho.*bc_Lu_b

        diffusion_bulk_v = mu_over_rho.*Lv
        diffusion_LS_v = mu_over_rho.*bc_Lv
        diffusion_border_v = mu_over_rho.*bc_Lv_b

    else
        # cf set_poisson_variable_coeff 


        #region Poisson variable coefficient
        
        # coeffD = viscosity 90 degrees ???


        coeffD = fzeros(grid)
        vec1(coeffD,grid) = mu_one_fluid

        recopy_bulk_to_border(coeffD,mu_one_fluid)  


        #interpolate coefficient
        coeffD_borders = vecb(coeffD,grid)
        interpolate_scalar!(grid, grid_u, grid_v, reshape(veci(coeffD,grid,1), grid), coeffDu, coeffDv)

        # coeffDx_bulk = veci(coeffDu,grid_u)
        # coeffDy_bulk = veci(coeffDv,grid_v)

        mat_coeffDx = Diagonal(vec(coeffDu)) 
        mat_coeffDy = Diagonal(vec(coeffDv)) 

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
        viscosity_coeff_for_du_dx[:,2:grid_u.nx] = grid.LS[end].geoL.cap[:,:,5]

        #TODO contact angle change  viscosity_coeff_for_du_dx[:,1] and at end (not interpolating right now)
        viscosity_coeff_for_du_dx[:,1] = viscosity_coeff_for_du_dx[:,2]
        viscosity_coeff_for_du_dx[:,end] = viscosity_coeff_for_du_dx[:,end-1]

        # arithmetic average 
        viscosity_coeff_for_du_dx .= 2 * (num.mu1  * viscosity_coeff_for_du_dx  .+ num.mu2  * viscosity_coeff_for_du_dx)

        # display(viscosity_coeff_for_du_dx)
        # PDI_mult_expose()

        diag_viscosity_coeff_for_du_dx = Diagonal(vec(viscosity_coeff_for_du_dx))
        #endregion Viscosity coefficient for \frac{\partial u}{\partial x}


        #region Viscosity coefficient for \frac{\partial u}{\partial y}


        viscosity_coeff_for_du_dy = zeros(grid_u.ny+1,grid_u.nx)

  
        #region interpolate 
        volume_fraction_full = create_2D_grid_volume_fraction(grid)
        grid_x_full_2D = create_2D_grid_x(grid,true,true)
        grid_y_full_2D = create_2D_grid_y(grid,true,true)

        all_grid_u_nodes_2D_x_for_du_dy_interp = create_2D_grid_x(gu,false,true)
        all_grid_u_nodes_2D_y_for_du_dy_interp = create_2D_grid_y(gu,false,true)

        printstyled(color=:magenta, @sprintf "\n grid_x_full_2D \n") 

        display(grid_x_full_2D)
        
        printstyled(color=:magenta, @sprintf "\n grid_y_full_2D \n") 

        display(grid_y_full_2D)

        printstyled(color=:magenta, @sprintf "\n x for viscosity_coeff_for_du_dy \n") 

        display(all_grid_u_nodes_2D_x_for_du_dy_interp)    
        # display(grid_u.x)    
        # display(x_centroid_u)    

        
        printstyled(color=:magenta, @sprintf "\n all_grid_u_nodes_2D_y_for_du_dy_interp \n") 

        display(all_grid_u_nodes_2D_y_for_du_dy_interp)


        for j in 1:grid_u.ny+1
            for i in 1:grid_u.nx

                print("\nvolume_fraction i ",i," j ",j,"\n")

                du_dy_coord_x = all_grid_u_nodes_2D_x_for_du_dy_interp[j,i] 

                du_dy_coord_y = (all_grid_u_nodes_2D_y_for_du_dy_interp[j,i] + all_grid_u_nodes_2D_y_for_du_dy_interp[j+1,i])/2 

                x1 = grid_x_full_2D[j,i]
                x2 = grid_x_full_2D[j,i+1]
                y1 = grid_y_full_2D[j,i]
                y2 = grid_y_full_2D[j+1,i+1]

                Q11 = volume_fraction_full[j,i]
                Q12 = volume_fraction_full[j+1,i]
                Q21 = volume_fraction_full[j,i+1]
                Q22 = volume_fraction_full[j+1,i+1]

                volume_fraction_face = bilinear_interpolation(du_dy_coord_x, du_dy_coord_y, x1, y1, x2, y2, Q11, Q12, Q21, Q22)

                # printstyled(color=:green, @sprintf "\n i %.5i j %.5i x %.2e y %.2e x1 %.2e y1 %.2e x2 %.2e y2 %.2e Q11 %.2e Q12 %.2e Q21 %.2e Q22 %.2e\n" i j du_dy_coord_x du_dy_coord_y x1 y1 x2 y2 Q11 Q12 Q21 Q22)
                # print("\nvolume_fraction i ",i," j ",j," ",volume_fraction_face)
               
                viscosity_coeff_for_du_dy[j,i] = harmonic_average_one_fluid(num.mu1,num.mu2,volume_fraction_face)

            end
        end
        
        printstyled(color=:cyan, @sprintf "\n viscosity_coeff_for_du_dy\n")

        display(viscosity_coeff_for_du_dy)

        #endregion interpolate




        diag_viscosity_coeff_for_du_dy = Diagonal(vec(viscosity_coeff_for_du_dy))

        #endregion Viscosity coefficient for \frac{\partial u}{\partial x}


        mul!(opC_u.tmp_x, diag_viscosity_coeff_for_du_dx * opC_u.iMx, opC_u.Bx)
        diffusion_bulk_u = opC_u.BxT * opC_u.tmp_x

        # \frac{\partial}{\partial x} \left( \mu \frac{\partial u}{\partial x} \right)
        mul!(opC_u.tmp_y, diag_viscosity_coeff_for_du_dy * opC_u.iMy, opC_u.By)
        diffusion_bulk_u = diffusion_bulk_u .+ opC_u.ByT * tmp_y



        #region v diffusion

        #region Viscosity coefficient for \frac{\partial v}{\partial y}
        #cf test in orientation.jl
        viscosity_coeff_for_dv_dy = zeros(grid_v.ny+1,grid_u.nx)

        viscosity_coeff_for_dv_dy[2:grid_u.ny,:] = grid.LS[end].geoL.cap[:,:,5]

        #TODO contact angle change  viscosity_coeff_for_dv_dy[:,1] and at end (not interpolating right now)
        viscosity_coeff_for_dv_dy[1,:] = viscosity_coeff_for_dv_dy[2,:]
        viscosity_coeff_for_dv_dy[end,:] = viscosity_coeff_for_dv_dy[end-1,:]

        # arithmetic average 
        viscosity_coeff_for_dv_dy .= 2 * (num.mu1  * viscosity_coeff_for_dv_dy  .+ num.mu2  * viscosity_coeff_for_dv_dy)

        # display(viscosity_coeff_for_dv_dy)
        # PDI_mult_expose()

        diag_viscosity_coeff_for_dv_dy = Diagonal(vec(viscosity_coeff_for_dv_dy))

        #endregion Viscosity coefficient for \frac{\partial v}{\partial y}

      
        mul!(opC_v.tmp_x, diag_viscosity_coeff_for_dv_dx * opC_v.iMx, opC_v.Bx)
        diffusion_bulk_v = opC_v.BxT * opC_v.tmp_x

        mul!(opC_v.tmp_y, diag_viscosity_coeff_for_dv_dy * opC_v.iMy, opC_v.By)
        diffusion_bulk_v = diffusion_bulk_u .+ opC_v.ByT * tmp_y

        #endregion v diffusion


        # TODO: check localisation of fields
        
        #TODO wall control volume kappa for bulk vs bulk control volume kappa different  
        # mul!(tmp_x, mat_coeffDx * iMx, Bx)
        # L = BxT * tmp_x
        # mul!(tmp_y, mat_coeffDy * iMy, By)
        # L = L .+ ByT * tmp_y
        
        
        #TODO divide by rho


        # diffusion_bulk_u   = mu_over_rho.*Lu
        diffusion_LS_u     = 0.0 * iRe.*bc_Lu
        diffusion_border_u = iRe.*bc_Lu_b

        # diffusion_bulk_v   = mu_over_rho.*Lv
        diffusion_LS_v     = 0.0 * iRe.*bc_Lu
        diffusion_border_u = iRe.*bc_Lu_b

        #TODO wall
        #TODO border ?
        
      

        coeffDu_border = copy(coeffDu)
        coeffDv_border = copy(coeffDv)

        # Interpolate conductivity at center of control volumes for potential gradient at the border
        interpolate_scalar_to_staggered_u_v_grids_at_border!(num,grid,coeffD,coeffDu_border,coeffDv_border)

        coeffDx_border = veci(coeffDu_border,grid_u)
        coeffDy_border = veci(coeffDv_border,grid_v)

        mat_coeffDx_b = Diagonal(vec(coeffDx_border)) 
        mat_coeffDy_b = Diagonal(vec(coeffDy_border))

        #Boundary for Laplacian
        # bc_L_b = (BxT * mat_coeffDx_b * iMx_b * Hx_b .+ ByT * mat_coeffDy_b * iMy_b  * Hy_b)


        bc_Lu_b = (opC_u.BxT * diag_viscosity_coeff_for_u_border_x * opC_u.iMx_b * opC_u.Hx_b .+ opC_u.ByT * diag_viscosity_coeff_for_u_border_y * opC_u.iMy_b * opC_u.Hy_b)
        
        bc_Lv_b = (opC_v.BxT * diag_viscosity_coeff_for_v_border_x * opC_v.iMx_b * opC_v.Hx_b .+ opC_v.ByT * diag_viscosity_coeff_for_v_border_y * opC_v.iMy_b * opC_v.Hy_b)

        #endregion Poisson variable coefficient

  

    end

    if num.pressure_velocity_coupling == 0
        if !navier
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
                bc_int, num, grid, grid_u, grid_v,
                opC_p, opC_u, opC_v,
                Auv, Buv,
                diffusion_bulk_u, diffusion_LS_u, diffusion_border_u, Mum1, BC_u,
                diffusion_bulk_v, diffusion_LS_v, diffusion_border_v, Mvm1, BC_v,
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

    elseif num.pressure_velocity_coupling > 1
        rhs_u = nothing
        rhs_v = nothing
        rhs_ϕ = nothing
        rhs_uv = FE_set_momentum_coupled2(
            bc_int, num, grid, grid_u, grid_v,
            opC_p, opC_u, opC_v,
            Auv, Buv,
            rhs_uv,
            diffusion_bulk_u, diffusion_LS_u, diffusion_border_u, Mum1, BC_u,
            diffusion_bulk_v, diffusion_LS_v, diffusion_border_v, Mvm1, BC_v,
            ls_advection,BC_p,ph
        )

    end
    
    return rhs_u, rhs_v, rhs_ϕ, rhs_uv, Lp, bc_Lp, bc_Lp_b, Lu, bc_Lu, bc_Lu_b, Lv, bc_Lv, bc_Lv_b
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
   @inbounds rhs[ntu-nbu+1:ntu] .= opu.χ_b * vec(a0_bu) #BC for u component on borders
   @inbounds rhs[ntu+ntv-nbv+1:ntu+ntv] .= opv.χ_b * vec(a0_bv) #BC for v component on borders

robin BC : source term a0
                
At the moment, the Levelset is not computed at borders/interfaces ? Only bulk
"""
function FE_set_momentum_coupled2_one_fluid(
    bc_interface, num, gp, gu, gv,
    opp, opu, opv,
    A, B,
    rhs,
    Lu, bc_Lu, bc_Lu_b, Mum1, BCu,
    Lv, bc_Lv, bc_Lv_b, Mvm1, BCv,
    ls_advection,BCp,ph=nothing
    )
    @unpack τ, Re, nLS, nNavier = num


    printstyled(color=:red, @sprintf "\n coupled pressure-velocity FE_set_momentum_coupled2\n")


    #region init
    iRe = num.visc_coeff

    nip = gp.nx * gp.ny
    nbp = 2 * gp.nx + 2 * gp.ny

    niu = gu.nx * gu.ny
    nbu = 2 * gu.nx + 2 * gu.ny
    ntu = (nLS - nNavier + 1) * niu + nbu

    niv = gv.nx * gv.ny
    nbv = 2 * gv.nx + 2 * gv.ny
    ntv = (nLS - nNavier + 1) * niv + nbv

    #Reset to zero
    rhs .= 0.0 

    #region BC borders u
    a0_bu = zeros(nbu)
    _a1_bu = zeros(nbu)
    _b_bu = zeros(nbu)
    for iLS in 1:num.nLS
        set_borders!(gu, gu.LS[iLS].cl, gu.LS[iLS].u, a0_bu, _a1_bu, _b_bu, BCu, num.n_ext_cl)
    end
    a1_bu = Diagonal(vec(_a1_bu))
    b_bu = Diagonal(vec(_b_bu))
    #endregion BC borders u

    #region BC borders v
    a0_bv = zeros(nbv)
    _a1_bv = zeros(nbv)
    _b_bv = zeros(nbv)
    for iLS in 1:num.nLS
        set_borders!(gv, gv.LS[iLS].cl, gv.LS[iLS].u, a0_bv, _a1_bv, _b_bv, BCv, num.n_ext_cl)
    end
    a1_bv = Diagonal(vec(_a1_bv))
    b_bv = Diagonal(vec(_b_bv))
    #endregion BC borders v


    #region BC borders p
    if num.pressure_velocity_coupling == 2
        # rhs = fnzeros(gp, num)

        a0_bp = zeros(nbp)
        _a1_bp = zeros(nbp)
        _b_bp = zeros(nbp)
        for iLS in 1:num.nLS
            set_borders_poisson!(gp, gp.LS[iLS].cl, gp.LS[iLS].u, a0_bp, _a1_bp, _b_bp, BCp, num.n_ext_cl)
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
    
    #         _a1 = ones(gp) .* __a1
    #         a1 = Diagonal(vec(_a1))
    #         _a2 = ones(gp) .* __a2
    #         a2 = Diagonal(vec(_a2))
    #         _b = ones(gp) .* __b
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

    #     veci(rhs,gp,iLS+1) .= +χ[iLS] * vec(a0[iLS]) #was - in set_poisson
    # end

    # vecb(rhs,gp) .= +χ_b * vec(a0_b) #was - in set_poisson
    #endregion BC borders p


    #endregion init

    # Array indices 
    bulk_u_velocity = 1:niu
    bulk_v_velocity = ntu+1:ntu+niv
    # bulk_tangential_velocity = 

    border_u_velocity = ntu-nbu+1:ntu
    border_v_velocity = ntu+ntv-nbv+1:ntu+ntv

    # nt = (num.nLS - num.nNavier + 1) * ni_uv + num.nNavier * nip + nb_uv + (num.nLS + 1) * nip + nbp
    ntNavier = num.nNavier * nip
    bulk_pressure = ntu+ntv+ntNavier+1:ntu+ntv+ntNavier+nip
    #ntu+ntv+1:ntu+ntv+nip
    border_pressure = ntu+ntv+ntNavier+(num.nLS + 1)*nip+1:ntu+ntv+ntNavier+(num.nLS + 1)*nip+nbp

    print("\n len rhs_uv ",size(rhs))

    print("bulk_pressure ",bulk_pressure)
    print("border_pressure ",border_pressure)


    if ls_advection
        A.nzval .= 0.0
        # Implicit part of viscous term
        A[bulk_u_velocity,bulk_u_velocity] = pad_crank_nicolson(opu.M .- τ .* Lu, gu, τ)
        # Contribution to implicit part of viscous term from outer boundaries
        A[bulk_u_velocity,border_u_velocity] = - τ .* bc_Lu_b

        # Boundary conditions for outer boundaries
        A[border_u_velocity,bulk_u_velocity] = b_bu * (opu.HxT_b * opu.iMx_b' * opu.Bx .+ opu.HyT_b * opu.iMy_b' * opu.By)
        A[border_u_velocity,border_u_velocity] = pad(b_bu * (
            opu.HxT_b * opu.iMx_bd * opu.Hx_b .+ 
            opu.HyT_b * opu.iMy_bd * opu.Hy_b
        ) .- opu.χ_b * a1_bu)

        # Implicit part of viscous term
        A[bulk_v_velocity,bulk_v_velocity] = pad_crank_nicolson(opv.M .- τ .* Lv, gv, τ)
        # Contribution to implicit part of viscous term from outer boundaries
        A[bulk_v_velocity,border_v_velocity] = - τ .* bc_Lv_b

        
        # Boundary conditions for outer boundaries
        A[border_v_velocity,bulk_v_velocity] = b_bv * (opv.HxT_b * opv.iMx_b' * opv.Bx .+ opv.HyT_b * opv.iMy_b' * opv.By)
        A[border_v_velocity,border_v_velocity] = pad(b_bv * (
            opv.HxT_b * opv.iMx_bd * opv.Hx_b .+ 
            opv.HyT_b * opv.iMy_bd * opv.Hy_b
        ) .- opv.χ_b * a1_bv)

        # TODO pad 1 or -4
        #TODO sign divergence not same u v and p

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
        # ∇ϕ_x = opu.AxT * opu.Rx * vec1(pD,grid) .+ opu.Gx_b * vecb(pD,grid)
        # ∇ϕ_y = opv.AyT * opv.Ry * vec1(pD,grid) .+ opv.Gy_b * vecb(pD,grid)
        # for iLS in 1:nLS
        #     ∇ϕ_x .+= opu.Gx[iLS] * veci(pD,grid,iLS+1)
        #     ∇ϕ_y .+= opv.Gy[iLS] * veci(pD,grid,iLS+1)
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
        
        B[bulk_u_velocity,bulk_u_velocity] = Mum1
        B[bulk_v_velocity,bulk_v_velocity] = Mvm1
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
        a0u, a1u, bu, a0v, a1v, bv, a0p, bp = set_velocity_boundary_conditions(bc_interface, iLS, gu, gv, gp, num)
        #endregion BC iLS


        interfacial_nb_1_u_velocity = _iLS*niu+1:(_iLS+1)*niu
        interfacial_nb_1_v_velocity = ntu+_iLS*niv+1:ntu+(_iLS+1)*niv

        if ls_advection
            if !is_navier_cl(bc_interface[iLS]) && !is_navier(bc_interface[iLS])
                #region not Navier
                # Contribution to implicit part of viscous term from inner boundaries
                A[bulk_u_velocity,interfacial_nb_1_u_velocity] = - τ .* bc_Lu[iLS]
                A[bulk_v_velocity,interfacial_nb_1_v_velocity] = - τ .* bc_Lv[iLS] 
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
                        sinα = Diagonal(vec(sin.(gp.LS[i].α)))
                        # replace!(sinα.diag, NaN=>0.0)
                        cosα = Diagonal(vec(cos.(gp.LS[i].α)))
                        # replace!(cosα.diag, NaN=>0.0)

                        if any(isnan, sinα.diag) || any(isnan, cosα.diag)
                            @error("NaN FE_set_momentum_coupled")
                            replace!(sinα.diag, NaN=>0.0)
                            replace!(cosα.diag, NaN=>0.0)
                        end

                        # not Navier, averaging coefficients differently computed 
                        interpolate_x = interpolating_coefficient_Navier(gu,gp,i,bc_interface[iLS])
                        interpolate_y = interpolating_coefficient_Navier(gv,gp,i,bc_interface[iLS])

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
                sinα_p = Diagonal(vec(sin.(gp.LS[iLS].α)))
                cosα_p = Diagonal(vec(cos.(gp.LS[iLS].α)))
                sinα_u = Diagonal(vec(sin.(gu.LS[iLS].α)))
                cosα_v = Diagonal(vec(cos.(gv.LS[iLS].α)))

                if any(isnan, sinα_p.diag) || any(isnan, cosα_p.diag) || any(isnan, sinα_u.diag) || any(isnan, cosα_v.diag)
                    @error("NaN FE_set_momentum_coupled")
                    replace!(sinα.diag, NaN=>0.0)
                    replace!(cosα.diag, NaN=>0.0)
                    replace!(sinα_u.diag, NaN=>0.0)
                    replace!(cosα_v.diag, NaN=>0.0)
                end

                # Contribution to implicit part of viscous term from inner boundaries
           
                # Navier BC, averaging coefficient computed differently
                interpolate_x = interpolating_coefficient_Navier(gu,gp,bc_interface[iLS])
                interpolate_y = interpolating_coefficient_Navier(gv,gp,bc_interface[iLS])
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
                
                interpolate_u_to_p, interpolate_v_to_p = interpolating_coefficient_Navier_uv_grids_to_p_grid_volume(num,gp,gu,gv,iLS)

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
                        interpolate_u_to_p, interpolate_v_to_p = interpolating_coefficient_Navier_uv_grids_to_p_grid_height(num,gp,gu,gv,i)

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
                sin_alpha_border.diag[1:gp.ny] .= sin.(gp.LS[iLS].α[:,1])
                sin_alpha_border.diag[gp.ny+1:gp.ny+gp.nx] .= sin.(gp.LS[iLS].α[1,:])
                sin_alpha_border.diag[gp.ny+gp.nx+1:2gp.ny+gp.nx] .= sin.(gp.LS[iLS].α[:,end])
                sin_alpha_border.diag[2gp.ny+gp.nx+1:end] .= sin.(gp.LS[iLS].α[end,:])
                cos_alpha_border = Diagonal(zeros(nbp))
                cos_alpha_border.diag[1:gp.ny] .= cos.(gp.LS[iLS].α[:,1])
                cos_alpha_border.diag[gp.ny+1:gp.ny+gp.nx] .= cos.(gp.LS[iLS].α[1,:])
                cos_alpha_border.diag[gp.ny+gp.nx+1:2gp.ny+gp.nx] .= cos.(gp.LS[iLS].α[:,end])
                cos_alpha_border.diag[2gp.ny+gp.nx+1:end] .= cos.(gp.LS[iLS].α[end,:])

                #region interpolation coefficients
                interpolate_u_to_p = spdiagm(nbp, nbu, 0 => zeros(nbp), 1 => zeros(nbp-1))
                for ii in 1:gp.ny
                    interpolate_u_to_p[ii,ii] = 1.0
                    interpolate_u_to_p[gp.ny+gp.nx+ii,gu.ny+gu.nx+ii] = 1.0
                end
                for ii in 1:gp.nx
                    interpolate_u_to_p[ii+gp.ny,ii+gu.ny] = 0.5
                    interpolate_u_to_p[ii+gp.ny,ii+gu.ny+1] = 0.5
                    interpolate_u_to_p[ii+2gp.ny+gp.nx,ii+2gu.ny+gu.nx] = 0.5
                    interpolate_u_to_p[ii+2gp.ny+gp.nx,ii+2gu.ny+gu.nx+1] = 0.5
                end
                interpolate_v_to_p = spdiagm(nbp, nbv, 0 => zeros(nbp), 1 => zeros(nbp-1))
                for ii in 1:gp.ny
                    interpolate_v_to_p[ii,ii] = 0.5
                    interpolate_v_to_p[ii,ii+1] = 0.5
                    interpolate_v_to_p[gp.ny+gp.nx+ii,gv.ny+gv.nx+ii] = 0.5
                    interpolate_v_to_p[gp.ny+gp.nx+ii,gv.ny+gv.nx+ii+1] = 0.5
                end
                for ii in 1:gp.nx
                    interpolate_v_to_p[ii+gp.ny,ii+gv.ny] = 1.0
                    interpolate_v_to_p[ii+2gp.ny+gp.nx,ii+2gv.ny+gv.nx] = 1.0
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
                for ii in 1:gu.ny
                    interpolate_u_to_p[ii,ii] = 1.0
                    interpolate_u_to_p[gu.ny+gu.nx+ii,gp.ny+gp.nx+ii] = 1.0
                end
                interpolate_u_to_p[gu.ny+1,gp.ny+1]  = 1.0
                interpolate_u_to_p[gu.ny+gu.nx,gp.ny+gp.nx]  = 1.0
                interpolate_u_to_p[2gu.ny+gu.nx+1,2gp.ny+gp.nx+1]  = 1.0
                interpolate_u_to_p[end,end]  = 1.0
                for ii in 2:(gu.nx-1)
                    interpolate_u_to_p[ii+gu.ny,ii+gp.ny-1] = 0.5
                    interpolate_u_to_p[ii+gu.ny,ii+gp.ny] = 0.5
                    interpolate_u_to_p[ii+2gu.ny+gu.nx,ii+2gp.ny+gp.nx-1] = 0.5
                    interpolate_u_to_p[ii+2gu.ny+gu.nx,ii+2gp.ny+gp.nx] = 0.5
                end
                interpolate_v_to_p = spdiagm(nbv, nbp, 0 => zeros(nbp), 1 => zeros(nbp-1))
                interpolate_v_to_p[1,1]  = 1.0
                interpolate_v_to_p[gv.ny,gp.ny]  = 1.0
                interpolate_v_to_p[gv.ny+gv.nx+1,gp.ny+gp.nx+1]  = 1.0
                interpolate_v_to_p[2gv.ny+gv.nx,2gp.ny+gp.nx]  = 1.0
                for ii in 2:(gv.ny-1)
                    interpolate_v_to_p[ii,ii-1] = 0.5
                    interpolate_v_to_p[ii,ii] = 0.5
                    interpolate_v_to_p[gv.ny+gv.nx+ii,gp.ny+gp.nx+ii] = 0.5
                    interpolate_v_to_p[gv.ny+gv.nx+ii,gp.ny+gp.nx+ii] = 0.5
                end
                for ii in 1:gv.nx
                    interpolate_v_to_p[ii+gv.ny,ii+gp.ny] = 1.0
                    interpolate_v_to_p[ii+2gv.ny+gv.nx,ii+2gp.ny+gp.nx] = 1.0
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
        set_first_cells!(A,rhs,gu,ntu-nbu,0,true,false,true,false)
        set_first_cells!(A,rhs,gv,ntu+ntv-nbv,ntu,false,true,false,true)
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

    ni_p = gp.nx * gp.ny
    nb_p = 2 * gp.nx + 2 * gp.ny

    ni_u = gu.nx * gu.ny
    nb_u = 2 * gu.nx + 2 * gu.ny

    ni_v = gv.nx * gv.ny
    nb_v = 2 * gv.nx + 2 * gv.ny

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
    # "nx"::Cstring, gp.nx::Ref{Clonglong}, PDI_OUT::Cint,
    # "ny"::Cstring, gp.ny::Ref{Clonglong}, PDI_OUT::Cint,
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

    # print("\n test Adummy\n",Adummy*uvD/gp.dx[1,1]^2)

    # print("\n test grid ",gp.dx[1,1])
    # print("\n factor ", factor )

    control_volumes = ones(ncol_A)

    # control_volumes[1:ntu] .= vec(gu.LS[1].geoL.dcap[:,:,5])
    # control_volumes[[ntu+1:ntu+ntv]] .= vec(gv.LS[1].geoL.dcap[:,:,5])
    # control_volumes[ntu+ntv+ntNavier+1:ntu+ntv+ntNavier+(num.nLS+1)*nip+nbp] .= vec(gp.LS[1].geoL.dcap[:,:,5])
    
    control_volumes[1:ni_u] .= vec(gu.LS[1].geoL.dcap[:,:,5]) #non-dimensionalise bulk rhs

    control_volumes[ntu+1:ntu+ni_v] .= vec(gv.LS[1].geoL.dcap[:,:,5]) #non-dimensionalise bulk rhs
    # control_volumes[ntu+ntv+ntNavier+1:ntu+ntv+ntNavier+(num.nLS+1)*nip+nbp] .= vec(gp.LS[1].geoL.dcap[:,:,5])

    #niu 

    # print("\n test Adummy\n",Adummy*uvD/factor/gp.dx[1,1]^2)

    # print("\n test volume\n",vec(gu.LS[1].geoL.dcap[:,:,5]))

    # print("\n test Adummy\n",Adummy*uvD/factor ./ control_volumes)

    

    test_matrix = zeros(ntu + ntv + nNavier * nip + (num.nLS + 1) * nip + nbp)

    test_matrix = Adummy*uvD/factor ./ control_volumes

    PDI_status = @ccall "libpdi".PDI_multi_expose("check_coupled_matrix"::Cstring,
    "vec_1D_len"::Cstring, length(test_matrix)::Ref{Clonglong}, PDI_OUT::Cint,
    "vec_1D"::Cstring, test_matrix::Ptr{Cdouble}, PDI_OUT::Cint,
    C_NULL::Ptr{Cvoid})::Cint
end