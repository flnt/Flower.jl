
"""
solves Navier-Stokes equations with a pressure projection method. 


#### Variables and Data Structures
- `vec1(u_predictionD, grid_u)`: Velocity correction for the horizontal grid_p.
- `vec1(v_predictionD, grid_v)`: Velocity correction for the vertical grid_p.
- `vec1(rhs_phi, grid_p)`: Right-hand side of the Poisson equation.
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
- `mass_transfer_rate`: Mass transfer rate.
- `pres_free_suface`: Free surface pressure.
- `diff_inv_rho`: Difference in inverse densities.
- `jump_mass_transfer_rate`: Flag for mass flux jump.
- `timestep_n`: Time step.
- `A_phi`: Matrix for the Poisson equation.
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
   - `vecb(v_predictionD, grid_v) .= uvD[border_v_velocity]`: Updates the vertical velocity correction.
   - `kill_dead_cells!(vec1(v_predictionD,grid_v), grid_v, geo_v[end])`: Removes dead cells from the vertical velocity correction grid_p.
   - `v_prediction .= reshape(vec1(v_predictionD,grid_v), grid_v)`: Reshapes the vertical velocity correction.

2. **Navier and Non-Navier Boundary Conditions**
   - Loop through linear solvers (`iLS`) to apply boundary conditions:
     - If not Navier or Navier-CL, update and apply boundary conditions for horizontal and vertical velocities.
     - If Navier or Navier-CL, update the Navier matrix.

3. **Divergence Calculation**
   - Calculate the divergence of the velocity corrections (`velocity_divergence`).
   - Add contributions from internal boundary conditions.

4. **Poisson Equation**
   - Set the right-hand side of the Poisson equation (`rhs_phi`).
   - Handle free surface conditions and Marangoni effects if `jump_mass_transfer_rate` is true.
   - Remove nullspace from the matrix `A_phi`.
   - Apply boundary conditions and solve the Poisson equation using `A_phi / rhs_phi`.

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
- The Poisson equation is solved using a linear solver (`A_phi / rhs_phi`).
- The pressure correction is applied to update the velocities.
- Dead cells are removed from the grids to maintain numerical stability.

---

This documentation provides a high-level overview of the code's functionality and the key operations performed. For detailed implementation of specific functions or operations, refer to the corresponding sections of the code.
"""
function solve_one_fluid_NS_no_phase!(
    time_scheme, bc_int,
    num, grid_p, geo, grid_u, geo_u, grid_v, geo_v,
    p, pD, ϕ, u, v, u_predictionD, v_predictionD, uD, vD, u_prediction, v_prediction, uT,
    pres_grad_x, pres_grad_y,
    phase_change_currently_activated,
    BC_u, BC_v, BC_p,
    opC_p, opC_u, opC_v, op_conv,
    Au, Bu, Av, Bv, A_phi, Auv, Buv,rhs_uv,
    Lpm1, bc_Lpm1, bc_Lpm1_b, Lum1, bc_Lum1, bc_Lum1_b, Lvm1, bc_Lvm1, bc_Lvm1_b,
    Cum1, Cvm1, Mum1, Mvm1,
    periodic_x, periodic_y, advection, ls_advection, current_iter, Ra, navier,
    volume_fraction,
    levelset_one_fluid,
    rho_one_fluid, 
    mu_one_fluid,
    rho_one_fluid_u, 
    # mu_one_fluid_u,
    rho_one_fluid_v, 
    # mu_one_fluid_v,
    volumic_surface_tension_u,
    volumic_surface_tension_v,
    convection_u,convection_v,
    viscosity_coeff_for_du_dx ,
    viscosity_coeff_for_du_dy ,
    viscosity_coeff_for_dv_dx ,
    viscosity_coeff_for_dv_dy,
    tmp_vec_p,
    tmp_vec_p0,
    rhs_phi,
    pres_free_suface,jump_mass_transfer_rate,mass_transfer_rate
    )

    ph=nothing

    @unpack Re, timestep_n, σ, g, β, nLS, nNavier = num
    # @unpack p, pD, ϕ, u, v, u_predictionD, v_predictionD, uD, vD, u_prediction, v_prediction, uT = ph
    @unpack Cu, Cv, CUTCu, CUTCv = op_conv

    u0 = copy(u)
    v0 = copy(v)

    idt = 1.0 / timestep_n
    # irho1 = 1.0 ./ rho_one_fluid
    # mu1_over_rho1 = mu_one_fluid ./ rho_one_fluid
    # mu1_over_rho1 = num.mu1 / num.rho1 


    # II = CartesianIndex(div(grid_v.ny,2),1)
    # pII = lexicographic(II,grid_v.ny)
    # print("\nLv[pII,:] ",Lvm1[pII,:])
    # print("\bc_Lvm1_b[pII,:] ",bc_Lvm1_b[pII,:])

  


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

    if num.pressure_velocity_coupling == 0
        if num.prediction == "PmIIimposedpressure" || 
            num.prediction == "PmIIimposedpressureBCincrement" || 
            num.prediction == "PmIIimposedpressure_nodiv" ||
            num.prediction == "testpressure"
            BC_Poisson = Boundaries() #Neumann everywhere
        elseif num.prediction == "PmIIimposedpressure_nodiv_2" || num.prediction == "PmIIimposedpressure_nodiv_4"
            BC_Poisson = Boundaries(top=Dirichlet()) 
        elseif num.prediction == "PmIIimposedpressure_nodiv_3"
            BC_Poisson = Boundaries(top=Dirichlet(),
                                    bottom=Dirichlet()) 
        else
            BC_Poisson = copy(BC_p) 
        end
    else 
        BC_Poisson = nothing 
    end
    
    if is_Forward_Euler(time_scheme)
        rhs_u, rhs_v, rhs_phi, rhs_uv, Lp, bc_Lp, bc_Lp_b, Lu, diffusion_LS_u, diffusion_border_u, Lv, diffusion_LS_v, diffusion_border_v = set_Forward_Euler_one_fluid!(
            bc_int, num, grid_p, geo, grid_u, geo_u, grid_v, geo_v,
            opC_p, opC_u, opC_v, BC_Poisson,BC_u, BC_v,
            Au, Bu, Av, Bv, A_phi, rhs_phi,Auv, Buv,
            volume_fraction,rho_one_fluid_u,rho_one_fluid_v,
            mass_transfer_rate,
            viscosity_coeff_for_du_dx ,
            viscosity_coeff_for_du_dy ,
            viscosity_coeff_for_dv_dx ,
            viscosity_coeff_for_dv_dy,
            Lpm1, bc_Lpm1, bc_Lpm1_b, Lum1, bc_Lum1, bc_Lum1_b, Lvm1, bc_Lvm1, bc_Lvm1_b,
            Mum1, Mvm1, op_conv, ph,
            periodic_x, periodic_y, advection, ls_advection, navier,rhs_uv,
        )
    elseif is_Crank_Nicolson(time_scheme)
        @error("\n Crank_Nicolson not implemented")
        # rhs_u, rhs_v, rhs_phi, Lp, bc_Lp, bc_Lp_b, Lu, diffusion_LS_u, diffusion_border_u, Lv, diffusion_LS_v, diffusion_border_v = set_Crank_Nicolson!(
        #     bc_int, num, grid_p, geo, grid_u, geo_u, grid_v, geo_v,
        #     opC_p, opC_u, opC_v, BC_Poisson, BC_u, BC_v,
        #     Au, Bu, Av, Bv, A_phi,
        #     Lpm1, bc_Lpm1, bc_Lpm1_b, Lum1, bc_Lum1, bc_Lum1_b, Lvm1, bc_Lvm1, bc_Lvm1_b,
        #     Mum1, Mvm1, mu1_over_rho1, op_conv, ph,
        #     periodic_x, periodic_y, advection, ls_advection
        # )
    end

    # ra_x = Ra .* sin(β) .* opC_u.M * vec(hcat(zeros(grid_u.ny), T))
    # ra_y = Ra .* cos(β) .* opC_v.M * vec(vcat(zeros(1,grid_v.nx), T))

   
      PDI_status = @ccall "libpdi".PDI_multi_expose("print_before_prediction"::Cstring,
    "u_1D"::Cstring, uD::Ptr{Cdouble}, PDI_OUT::Cint,
    "v_1D"::Cstring, vD::Ptr{Cdouble}, PDI_OUT::Cint,
    "p_1D"::Cstring, pD::Ptr{Cdouble}, PDI_OUT::Cint,
    C_NULL::Ptr{Cvoid})::Cint

    #region prediction 
   
    #region add gradient of pressure to prediction
    # Compute gradient of pressure localized on u and v grids , times volume
    # $ \nabla p^{n-1/2} $
    if num.pressure_velocity_coupling == 0
        if num.prediction == "PmI" || 
        num.prediction == "PmII" || 
        num.prediction == "PmIIimposedpressure" || 
        num.prediction == "PmIIimposedpressureBCincrement" || 
        num.prediction == "PmIIimposedpressure_nodiv" ||
        num.prediction == "PmIIimposedpressure_nodiv_2" ||
        num.prediction == "PmIIimposedpressure_nodiv_3" ||
        num.prediction == "PmIIimposedpressure_nodiv_4" ||
        num.prediction == "testpressure"

            #cf Brown 2001

            ∇ϕ_x = opC_u.AxT * opC_u.Rx * vec1(pD,grid_p) .+ opC_u.Gx_b * vecb(pD,grid_p)
            ∇ϕ_y = opC_v.AyT * opC_v.Ry * vec1(pD,grid_p) .+ opC_v.Gy_b * vecb(pD,grid_p)
            #region cut-cell
            # for iLS in 1:nLS
            #     ∇ϕ_x .+= opC_u.Gx[iLS] * veci(pD,grid_p,iLS+1)
            #     ∇ϕ_y .+= opC_v.Gy[iLS] * veci(pD,grid_p,iLS+1)
            # end
            #endregion cut-cell

            pres_grad_x .= 0.0 #TODO
            pres_grad_y .= 0.0
            
            # gradient \times volume of cell (cells at border: volume = dx*dy/2)
            pres_grad_x .= copy(∇ϕ_x) #∇ϕ_x
            pres_grad_y .= copy(∇ϕ_y) #∇ϕ_y

            # gradient check
            PDI_status = @ccall "libpdi".PDI_multi_expose("print_pressure_gradient_in_prediction"::Cstring,
            "grad_x_1D"::Cstring, pres_grad_x::Ptr{Cdouble}, PDI_OUT::Cint,
            "grad_y_1D"::Cstring, pres_grad_y::Ptr{Cdouble}, PDI_OUT::Cint,
            "p_1D"::Cstring, pD::Ptr{Cdouble}, PDI_OUT::Cint,
            C_NULL::Ptr{Cvoid})::Cint

            ∇ϕ_x .= 0.0
            ∇ϕ_y .= 0.0

            grad_x = zeros(grid_u)
            grad_y = zeros(grid_v)
            compute_grad_T_x_T_y_array_u_v_capacities!(num, grid_p, grid_u, grid_v, opC_u, opC_v, grad_x, grad_y, pD)
        
            #TODO divergence level
            
        
            PDI_status = @ccall "libpdi".PDI_multi_expose("check_pressure_velocity_end"::Cstring,
            # "grad_x"::Cstring,grad_x::Ptr{Cdouble}, PDI_OUT::Cint,
            # "grad_y"::Cstring, grad_y::Ptr{Cdouble}, PDI_OUT::Cint,
            "grad_u"::Cstring,grad_x::Ptr{Cdouble}, PDI_OUT::Cint,
            "grad_v"::Cstring, grad_y::Ptr{Cdouble}, PDI_OUT::Cint,
            "u_1D"::Cstring, u_predictionD::Ptr{Cdouble}, PDI_OUT::Cint,
            "v_1D"::Cstring, v_predictionD::Ptr{Cdouble}, PDI_OUT::Cint,
            "p_1D"::Cstring, pD::Ptr{Cdouble}, PDI_OUT::Cint,
            C_NULL::Ptr{Cvoid})::Cint

            PDI_status = @ccall "libpdi".PDI_multi_expose("grad_pres_y"::Cstring,
            # "grad_x"::Cstring,grad_x::Ptr{Cdouble}, PDI_OUT::Cint,
            # "grad_y"::Cstring, grad_y::Ptr{Cdouble}, PDI_OUT::Cint,
            # "grad_u"::Cstring,grad_x::Ptr{Cdouble}, PDI_OUT::Cint,
            "grad_pres_y"::Cstring, grad_y::Ptr{Cdouble}, PDI_OUT::Cint,
            # "grad_pres_coupled_y"::Cstring, grad_y[2:end-1,:]::Ptr{Cdouble}, PDI_OUT::Cint,
            # "u_1D"::Cstring, u_predictionD::Ptr{Cdouble}, PDI_OUT::Cint,
            # "v_1D"::Cstring, v_predictionD::Ptr{Cdouble}, PDI_OUT::Cint,
            # "p_1D"::Cstring, pD::Ptr{Cdouble}, PDI_OUT::Cint,
            C_NULL::Ptr{Cvoid})::Cint
        

        end
    end
    #endregion add gradient to prediction



    #region convection
    if num.non_dimensionalize == 0
        # Cui = Cu * vec(u) .+ CUTCu
        # Cvi = Cv * vec(v) .+ CUTCv
        Cui = Cu * vec(u) 
        Cvi = Cv * vec(v) 

        # print("\n max CUTCv ",maximum(CUTCv))


    elseif num.non_dimensionalize == -1 
        # Element-wise multiplication
        Cui = Cu * vec(rho_one_fluid_u .* u) .+ CUTCu #_rho should be included  #TODO 
        Cvi = Cv * vec(rho_one_fluid_v .* v) .+ CUTCv #_rho
    end

    if num.convection == 0
       
        if advection
            # scheme
            if current_iter == 1
                convection_u .+= Cui
                convection_v .+= Cvi
            else
                convection_u .+= 1.5 .* Cui .- 0.5 .* Cum1 #Cui returned at the end of function to Cum1
                convection_v .+= 1.5 .* Cvi .- 0.5 .* Cvm1
            end
        end
    else
        convection_u .= 0.0
        convection_v .= 0.0

    end


    #endregion convection


    



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



    PDI_status = @ccall "libpdi".PDI_multi_expose("print_matrix"::Cstring,
    "Auv_n"::Cstring, Auv.n::Ref{Clonglong}, PDI_OUT::Cint,
    "Auv_m"::Cstring, Auv.m::Ref{Clonglong}, PDI_OUT::Cint,
    "Auv_colptr_len"::Cstring, length(Auv.colptr)::Ref{Clonglong}, PDI_OUT::Cint,
    "Auv_rowval_len"::Cstring, length(Auv.rowval)::Ref{Clonglong}, PDI_OUT::Cint,
    "Auv_nzval_len"::Cstring, length(Auv.nzval)::Ref{Clonglong}, PDI_OUT::Cint,
    "Auv_colptr_1D"::Cstring, Auv.colptr::Ptr{Clonglong}, PDI_OUT::Cint,
    "Auv_rowval_1D"::Cstring, Auv.rowval::Ptr{Clonglong}, PDI_OUT::Cint,
    "Auv_nzval_1D"::Cstring, Auv.nzval::Ptr{Cdouble}, PDI_OUT::Cint,
    C_NULL::Ptr{Cvoid})::Cint


    #TODO doc beta 

    if num.prediction == "Flowergravproj"
        grav_x = fzeros(grid_u)
        grav_y = fzeros(grid_v)
    else

        if num.non_dimensionalize == 0
            # grav_x = g .* sin(β) .* opC_u.M * fones(grid_u)
            # grav_y = g .* cos(β) .* opC_v.M * fones(grid_v)
            # print("\n test gravity")

            grav_x = fzeros(grid_u)
            grav_y = g .* opC_v.M * fones(grid_v)

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
    end

    # grav_x .= 0.0 
    # grav_y .= 0.0
    # print("\n test grav...")

    # print("\n opC_u.M ",opC_u.M)
    # print("\n opC_v.M ",opC_v.M)

    
    # printstyled(color=:red, @sprintf "\n gravity \n")

    # display(grav_y)

    rhs_uv[bulk_u_velocity] .-= timestep_n .* grav_x #timestep_n * rho_one_fluid_u .* grav_x

    rhs_uv[bulk_u_velocity] .-= timestep_n .* convection_u #rho in convection_u

     PDI_status = @ccall "libpdi".PDI_multi_expose("rhs_uv"::Cstring,
    "rhs_uv_len"::Cstring, length(rhs_uv)::Ref{Clonglong}, PDI_OUT::Cint,
    "rhs_uv_1D"::Cstring, rhs_uv::Ptr{Cdouble}, PDI_OUT::Cint,
    C_NULL::Ptr{Cvoid})::Cint

    # print("\n volumic_surface_tension_u")
    # display(volumic_surface_tension_u)
    # display(volumic_surface_tension_v)

    # print("\n rhs_uv",rhs_uv)
    # print("\n ph GGGGG",rhs_uv)

    # display(pres_grad_x)

    # display(pres_grad_y)

    # pres_grad_x .= 0.0
    # pres_grad_y .= 0.0



    # rhs_uv[bulk_u_velocity] .+= timestep_n .* ra_x
    if num.pressure_velocity_coupling == 0
        if num.non_dimensionalize == 0
            rhs_uv[bulk_u_velocity] .-= timestep_n .* pres_grad_x ./ vec(rho_one_fluid_u)

            print("\n grad pressure u")
            PDI_status = @ccall "libpdi".PDI_multi_expose("rhs_uv"::Cstring,
            "rhs_uv_len"::Cstring, length(rhs_uv)::Ref{Clonglong}, PDI_OUT::Cint,
            "rhs_uv_1D"::Cstring, rhs_uv::Ptr{Cdouble}, PDI_OUT::Cint,
            C_NULL::Ptr{Cvoid})::Cint
        else
            rhs_uv[bulk_u_velocity] .-= timestep_n .* pres_grad_x 
        end
    end 

    # rhs_uv[bulk_u_velocity] .+= timestep_n .* ra_x
    if num.non_dimensionalize == 0
        
        #Surface tension
        rhs_uv[bulk_u_velocity] .+= timestep_n .* vec(volumic_surface_tension_u) ./ vec(rho_one_fluid_u)

        print("\n surface tension u")
        PDI_status = @ccall "libpdi".PDI_multi_expose("rhs_uv"::Cstring,
        "rhs_uv_len"::Cstring, length(rhs_uv)::Ref{Clonglong}, PDI_OUT::Cint,
        "rhs_uv_1D"::Cstring, rhs_uv::Ptr{Cdouble}, PDI_OUT::Cint,
        C_NULL::Ptr{Cvoid})::Cint

    else
        #Surface tension
        rhs_uv[bulk_u_velocity] .+= timestep_n .* vec(volumic_surface_tension_u)

    end

    

   
    # print("\n grav_y ",grav_y)

    rhs_uv[bulk_v_velocity] .-= timestep_n .* grav_y
    
    print("\n grav y")

    PDI_status = @ccall "libpdi".PDI_multi_expose("rhs_uv"::Cstring,
    "rhs_uv_len"::Cstring, length(rhs_uv)::Ref{Clonglong}, PDI_OUT::Cint,
    "rhs_uv_1D"::Cstring, rhs_uv::Ptr{Cdouble}, PDI_OUT::Cint,
    C_NULL::Ptr{Cvoid})::Cint

    rhs_uv[bulk_v_velocity] .-= timestep_n .* convection_v

    conv_y = reshape(convection_v,grid_v)

    grav_y_2D = reshape(grav_y,grid_v)

    PDI_status = @ccall "libpdi".PDI_multi_expose("conv_y"::Cstring,
    # "grad_x"::Cstring,grad_x::Ptr{Cdouble}, PDI_OUT::Cint,
    # "grad_y"::Cstring, grad_y::Ptr{Cdouble}, PDI_OUT::Cint,
    # "grad_u"::Cstring,grad_x::Ptr{Cdouble}, PDI_OUT::Cint,
    "conv_y"::Cstring, conv_y::Ptr{Cdouble}, PDI_OUT::Cint,
    # "u_1D"::Cstring, u_predictionD::Ptr{Cdouble}, PDI_OUT::Cint,
    # "v_1D"::Cstring, v_predictionD::Ptr{Cdouble}, PDI_OUT::Cint,
    # "p_1D"::Cstring, pD::Ptr{Cdouble}, PDI_OUT::Cint,
    C_NULL::Ptr{Cvoid})::Cint

    PDI_status = @ccall "libpdi".PDI_multi_expose("grav_y"::Cstring,
    # "grad_x"::Cstring,grad_x::Ptr{Cdouble}, PDI_OUT::Cint,
    # "grad_y"::Cstring, grad_y::Ptr{Cdouble}, PDI_OUT::Cint,
    # "grad_u"::Cstring,grad_x::Ptr{Cdouble}, PDI_OUT::Cint,
    "grav_y"::Cstring, grav_y_2D::Ptr{Cdouble}, PDI_OUT::Cint,
    # "u_1D"::Cstring, u_predictionD::Ptr{Cdouble}, PDI_OUT::Cint,
    # "v_1D"::Cstring, v_predictionD::Ptr{Cdouble}, PDI_OUT::Cint,
    # "p_1D"::Cstring, pD::Ptr{Cdouble}, PDI_OUT::Cint,
    C_NULL::Ptr{Cvoid})::Cint

    print("\n conv y")

    PDI_status = @ccall "libpdi".PDI_multi_expose("rhs_uv"::Cstring,
    "rhs_uv_len"::Cstring, length(rhs_uv)::Ref{Clonglong}, PDI_OUT::Cint,
    "rhs_uv_1D"::Cstring, rhs_uv::Ptr{Cdouble}, PDI_OUT::Cint,
    C_NULL::Ptr{Cvoid})::Cint

    # print("\n test rhs 0 ")
    # rhs_uv .= 0.0

    # rhs_uv[bulk_v_velocity] .+= timestep_n .* ra_y
    if num.pressure_velocity_coupling == 0
         
        if num.non_dimensionalize == 0
            rhs_uv[bulk_v_velocity] .-= timestep_n .* pres_grad_y ./ vec(rho_one_fluid_v)

            print("\n grad pressure y")

            PDI_status = @ccall "libpdi".PDI_multi_expose("rhs_uv"::Cstring,
            "rhs_uv_len"::Cstring, length(rhs_uv)::Ref{Clonglong}, PDI_OUT::Cint,
            "rhs_uv_1D"::Cstring, rhs_uv::Ptr{Cdouble}, PDI_OUT::Cint,
            C_NULL::Ptr{Cvoid})::Cint
        else
            rhs_uv[bulk_v_velocity] .-= timestep_n .* pres_grad_y
            
        end
    end


          
    if num.non_dimensionalize == 0

        #Surface tension
        rhs_uv[bulk_v_velocity] .+= timestep_n .* vec(volumic_surface_tension_v) ./  vec(rho_one_fluid_v)

        print("\n surface tension y")

        PDI_status = @ccall "libpdi".PDI_multi_expose("rhs_uv"::Cstring,
        "rhs_uv_len"::Cstring, length(rhs_uv)::Ref{Clonglong}, PDI_OUT::Cint,
        "rhs_uv_1D"::Cstring, rhs_uv::Ptr{Cdouble}, PDI_OUT::Cint,
        C_NULL::Ptr{Cvoid})::Cint

    else
        #Surface tension
        rhs_uv[bulk_v_velocity] .+= timestep_n .* vec(volumic_surface_tension_v)
    end


    #region phase change
    range_divergence = ntu+ntv+ntNavier+1:ntu+ntv+ntNavier+nip

    # if num.phase_change_method == 2 && num.phase_change_currently_activated == 1
    #     rhs_uv[range_divergence] = vec(mass_transfer_rate * ( 1.0/num.rho1 - 1.0/num.rho2 ) )
    #     # print("\n TODO sign factor divergence and Dirac and one sided, not a problem ?", range_divergence)
    # elseif num.phase_change_method == 4 && num.phase_change_currently_activated == 1 #TODO
    #     # iM = Diagonal(inv_weight_eps2.(num.epsilon_mode,num.epsilon_vol,vec(geo[end].dcap[:,:,5])))

    #     # rhs_uv[range_divergence] = vec(geo[end].dcap[:,:,5] * mass_transfer_rate * ( 1.0/num.rho1 - 1.0/num.rho2 ) )
    #     # Dirac so intfc length
    #     # print("\n opC_p.χ[1] * mass_transfer_rate ",opC_p.χ[1] * mass_transfer_rate)
    #     # Dirac : opC_p.χ[1]/geo[end].dcap[:,:,5] so volume integrated intfc len * \dot m (1/rho-...)
    #     # rhs_uv[range_divergence] = opC_p.χ[1] * vec( mass_transfer_rate * ( 1.0/num.rho1 - 1.0/num.rho2 ) )
    #     rhs_uv[range_divergence] =  vec( mass_transfer_rate * ( 1.0/num.rho1 - 1.0/num.rho2 ) ) #opC_p.χ[1] in the redistribution

    #     # II = CartesianIndex(3,37)
    #     # pII =lexicographic(II,grid_p.ny)
    #     # print("\n op.χ[1] NS",opC_p.χ[1].diag[pII])
    # end

    if phase_change_currently_activated == 1 && num.pressure_velocity_coupling == 3
        rhs_uv[range_divergence] =  vec( mass_transfer_rate * ( 1.0/num.rho1 - 1.0/num.rho2 ) )

        PDI_status = @ccall "libpdi".PDI_multi_expose("rhs_uv_divergence"::Cstring,
        "rhs_uv_divergence"::Cstring, rhs_uv[range_divergence]::Ptr{Cdouble}, PDI_OUT::Cint,
        C_NULL::Ptr{Cvoid})::Cint

        #TODO BC u v mass transfer outflow , so Neumann
    end

    #endregion phase change

    NS_force_y = reshape(-grav_y .-pres_grad_y ./ vec(rho_one_fluid_v),grid_v)

    PDI_status = @ccall "libpdi".PDI_multi_expose("NS_force_y"::Cstring,
    # "grad_x"::Cstring,grad_x::Ptr{Cdouble}, PDI_OUT::Cint,
    # "grad_y"::Cstring, grad_y::Ptr{Cdouble}, PDI_OUT::Cint,
    # "grad_u"::Cstring,grad_x::Ptr{Cdouble}, PDI_OUT::Cint,
    "NS_force_y"::Cstring, NS_force_y::Ptr{Cdouble}, PDI_OUT::Cint,
    # "u_1D"::Cstring, u_predictionD::Ptr{Cdouble}, PDI_OUT::Cint,
    # "v_1D"::Cstring, v_predictionD::Ptr{Cdouble}, PDI_OUT::Cint,
    # "p_1D"::Cstring, pD::Ptr{Cdouble}, PDI_OUT::Cint,
    C_NULL::Ptr{Cvoid})::Cint
   
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

    PDI_status = @ccall "libpdi".PDI_multi_expose("write_rhs_uv_v"::Cstring,
    "rhs_uv_v"::Cstring, rhs_uv[bulk_v_velocity]::Ptr{Cdouble}, PDI_OUT::Cint,
    C_NULL::Ptr{Cvoid})::Cint

    #region check_one_fluid
    # check_one_fluid = false
    # if check_one_fluid
    #     II = CartesianIndex(div(grid_u.ny,4),div(grid_u.nx,2)) #center
    #     pII = lexicographic(II, grid_u.ny)
    #     print("\n test A ",Auv[pII,:])

    #     print("\n test rhs ",rhs_uv[pII])


    #     II = CartesianIndex(div(grid_u.ny,2),div(grid_u.nx,2)) #center
    #     pII = lexicographic(II, grid_u.ny)
    #     print("\n test A ",Auv[pII,:])

    #     print("\n test rhs ",rhs_uv[pII])

    #     # print("\n size Auv ",size(Auv))

    #     # print("\n test rhs 0 ")
    #     # rhs_uv .= 0.0
    # end
    #endregion check_one_fluid


    print("\n test A")
    II = CartesianIndex(1,1)
    pII = lexicographic(II, grid_v.ny)
    print("\n test A v",Auv[pII+ntu,:])
    print("\n end test A \n")


    if num.pressure_velocity_coupling == 3
        print("\n Setting first cells")
        if BC_u.left isa Dirichlet
            print("\n Setting first cells: Dirichlet")
            set_first_cells!(Auv,rhs_uv,grid_u,ntu-nbu,0,true,false,true,false)
        else
            print("\n Setting first cells: Neumann")
            set_first_cells_Neumann!(Auv,rhs_uv,grid_u,ntu-nbu,0,true,false,true,false)
        end

        if BC_v.bottom isa Dirichlet
            print("\n Setting first cells: Dirichlet")
            set_first_cells!(Auv,rhs_uv,grid_v,ntu+ntv-nbv,ntu,false,true,false,true)
        else
            print("\n Setting first cells: Neumann")
            set_first_cells_Neumann!(Auv,rhs_uv,grid_v,ntu+ntv-nbv,ntu,false,true,false,true)
        end
        print("\n End setting first cells \n")
    end

    print("\n test A")
    II = CartesianIndex(1,1)
    pII = lexicographic(II, grid_v.ny)
    print("\n test A v",Auv[pII+ntu,:])
    print("\n end test A \n")

    
    #region solver
    try
        @time uvD .= Auv \ rhs_uv
    catch e
        uvD .= Inf
        printstyled(color=:red, @sprintf "\n --------------------------------------------------------------\n")
        printstyled(color=:red, @sprintf "\n NS solver error solve_one_fluid_NS_no_phase!\n")
        println(e)
        @error("NS solver error solve_one_fluid_NS_no_phase!")

        num.status = 1
        print("\n rhs_uv ", rhs_uv)
        display(Auv)
        printstyled(color=:red, @sprintf "\n --------------------------------------------------------------\n")

    end
    #endregion solver

    
    #region check_one_fluid
    # check_one_fluid = false
    # if check_one_fluid
    #     for j in 1:grid_v.ny
    #         II = CartesianIndex(j,div(grid_v.nx,2)) #center
    #         pII = lexicographic(II, grid_v.ny)
    #         IIp = CartesianIndex(j,div(grid_p.nx,2)) #center
    #         pIIp = lexicographic(IIp, grid_p.ny)
    #         print("\n test A v",II," ",pII," v ",uvD[pII+ntu]," rhs ",rhs_uv[pII+ntu]," ",uvD[ntu+ntv+pIIp-1]," ",uvD[ntu+ntv+pIIp]," ",uvD[ntu+ntv+pIIp+1]," grad ",(uvD[ntu+ntv+pIIp+1]-uvD[ntu+ntv+pIIp])*40," ",(uvD[ntu+ntv+pIIp]-uvD[ntu+ntv+pIIp-1])*40," ",Auv[pII+ntu,:])
    #     end

    #     II = CartesianIndex(div(grid_v.ny,4),div(grid_v.nx,2)) #center
    #     pII = lexicographic(II, grid_v.ny)

    #     print("\n test A v",Auv[pII+ntu,:])
    #     print("\n test A v",uvD[pII+ntu])
    #     print("\n test rhs ",rhs_uv[pII+ntu])

    #     print("\n pII +ntu ",pII+ntu)
    #     print("\n test A v",uvD[pII+ntu+1])
    #     print("\n test A v",uvD[pII+ntu-1])
    #     print("\n test A v",uvD[pII+ntu+grid_v.ny])
    #     print("\n test A v",uvD[pII+ntu-grid_v.ny])

    #     IIp = CartesianIndex(div(grid_p.ny,4),div(grid_p.nx,2)) #center
    #     pIIp = lexicographic(IIp, grid_p.ny)
    #     print("\n test A p",uvD[ntu+ntv+pIIp])
    #     print("\n test A p",uvD[ntu+ntv+pIIp+1])
    #     print("\n test A p",uvD[ntu+ntv+pIIp-1])


    #     II = CartesianIndex(div(grid_u.ny,4),div(grid_u.nx,2)) #center
    #     pII = lexicographic(II, grid_u.ny)
    #     print("\n test A ",Auv[pII,:])

    

    #     PDI_status = @ccall "libpdi".PDI_multi_expose("check_pressure_rising"::Cstring,
    #     # "u_1D"::Cstring, u_predictionD::Ptr{Cdouble}, PDI_OUT::Cint,
    #     # "v_1D"::Cstring, v_predictionD::Ptr{Cdouble}, PDI_OUT::Cint,
    #     "p_1D"::Cstring, pD::Ptr{Cdouble}, PDI_OUT::Cint,
    #     C_NULL::Ptr{Cvoid})::Cint

    #     # print("\n test rhs ",rhs_uv[pII])
    #     # print("\n pII ",pII)
    #     # print("\n pII ",uvD[pII])
    #     # print("\n pII ",uvD[1539])
    #     # print("\n pII ",uvD[pII+1])
    #     # print("\n pII ",uvD[pII+2])
        
    #     # print("\n pII ",uvD[1620])
    #     # print("\n pII ",uvD[5000])
    #     # print("\n pII ",uvD[5000+grid_p.ny])
    #     # print("\n pII ",uvD[5000])
    #     # print("\n pII ",uvD[5000])

    #     # print("\n pII ",uvD[8464])
    #     # print("\n pII ",uvD[8544])

        
    #     # vec1(uD,grid_u) .= uvD[1:niu]
    #     # vecb(uD,grid_u) .= uvD[niu+1:ntu]

    #     # vec1(vD,grid_v) .= uvD[ntu+1:ntu+niv]
    #     # vecb(vD,grid_v) .= uvD[ntu+1+niv:ntu+ntv]




    # #     [1539]  =  -0.0001
    # #   [1540]  =  0.001225
    # #   [1541]  =  -0.0001
    # #   [1620]  =  -0.0002
    # #   [5000]  =  -0.0001
    # #   [5001]  =  0.0001
    # #   [5081]  =  0.0001
    # #   [5082]  =  -0.0001
    # #   [8464]  =  -2.5e-6
    # #   [8544]  =  2.5e-6



    #     II = CartesianIndex(div(grid_u.ny,2),div(grid_u.nx,2)) #center
    #     pII = lexicographic(II, grid_u.ny)
    #     print("\n test A ",Auv[pII,:])

    #     print("\n test rhs ",rhs_uv[pII])
    
    # end
    #endregion check_one_fluid

    vec1(u_predictionD, grid_u) .= uvD[bulk_u_velocity]
    vecb(u_predictionD, grid_u) .= uvD[border_u_velocity]
    #region cut-cell
    # kill_dead_cells!(vec1(u_predictionD,grid_u), grid_u, geo_u[end])
    #endregion cut-cell
    
    u_prediction .= reshape(vec1(u_predictionD,grid_u), grid_u)

    vec1(v_predictionD, grid_v) .= uvD[bulk_v_velocity]
    vecb(v_predictionD, grid_v) .= uvD[border_v_velocity]
    #region cut-cell
    # kill_dead_cells!(vec1(v_predictionD,grid_v), grid_v, geo_v[end])
    #endregion cut-cell
    
    v_prediction .= reshape(vec1(v_predictionD,grid_v), grid_v)


    II = CartesianIndex(div(grid_v.ny,4),div(grid_v.nx,2)) #center
    pII = lexicographic(II, grid_v.ny)

    print("\n test A v",uvD[pII+ntu])

    print("\n v_prediction ",v_prediction[II]," ",v_predictionD[pII])

    #region cut-cell
    # nNav = 0
    # _iLS = 1
    # for iLS in 1:nLS
    #     if !is_navier(bc_int[iLS]) && !is_navier_cl(bc_int[iLS])
    #         veci(u_predictionD,grid_u,iLS+1) .= uvD[_iLS*niu+1:(_iLS+1)*niu]
    #         kill_dead_cells!(veci(u_predictionD,grid_u,iLS+1), grid_u, geo_u[end])

    #         veci(v_predictionD,grid_v,iLS+1) .= uvD[ntu+_iLS*niv+1:ntu+(_iLS+1)*niv]
    #         kill_dead_cells!(veci(v_predictionD,grid_v,iLS+1), grid_v, geo_v[end])
    #         _iLS += 1
    #     else
    #         @inbounds uT[nNav+1,:] .= vec(uvD[ntu+ntv+1+nNav*nip:ntu+ntv+(nNav+1)*nip])
    #         nNav += 1
    #     end
    # end
    #endregion cut-cell

    #endregion Navier


    PDI_status = @ccall "libpdi".PDI_multi_expose("print_velocity_prediction"::Cstring,
    "u_1D"::Cstring, u_predictionD::Ptr{Cdouble}, PDI_OUT::Cint,
    "v_1D"::Cstring, v_predictionD::Ptr{Cdouble}, PDI_OUT::Cint,
    "p_1D"::Cstring, pD::Ptr{Cdouble}, PDI_OUT::Cint,
    C_NULL::Ptr{Cvoid})::Cint

    u_prediction = reshape(vec1(u_predictionD,grid_u),grid_u)
    v_prediction = reshape(vec1(v_predictionD,grid_v),grid_v)

    PDI_status = @ccall "libpdi".PDI_multi_expose("write_velocity_prediction"::Cstring,
    "u_prediction"::Cstring, u_prediction::Ptr{Cdouble}, PDI_OUT::Cint,
    "v_prediction"::Cstring, v_prediction::Ptr{Cdouble}, PDI_OUT::Cint,
    # "p_1D"::Cstring, pD::Ptr{Cdouble}, PDI_OUT::Cint,
    C_NULL::Ptr{Cvoid})::Cint

    # II = CartesianIndex(div(grid_v.ny,2),1)
    # pII = lexicographic(II,grid_v.ny)
    # print("\n A coeff ",Av[pII,:])
    # print("\nM ",opC_p.iMy.diag[pII])

    # print("\nAv[pII,:]./mu1_over_rho1 ",Av[pII,:]./mu1_over_rho1)

    # print("\nmu1_over_rho1 ",mu1_over_rho1)


    # # Test analytical vel
    # u_predictionD = copy(uD)
    # v_predictionD = copy(vD)
    # print("\n test analytical vel ")

    # PDI_status = @ccall "libpdi".PDI_multi_expose("print_velocity_prediction"::Cstring,
    # "u_1D"::Cstring, u_predictionD::Ptr{Cdouble}, PDI_OUT::Cint,
    # "v_1D"::Cstring, v_predictionD::Ptr{Cdouble}, PDI_OUT::Cint,
    # "p_1D"::Cstring, pD::Ptr{Cdouble}, PDI_OUT::Cint,
    # C_NULL::Ptr{Cvoid})::Cint


    #endregion prediction 

    #region correction 

    # Compute divergence of velocity
    velocity_divergence = opC_p.AxT * vec1(u_predictionD,grid_u) .+ opC_p.Gx_b * vecb(u_predictionD,grid_u) .+
          opC_p.AyT * vec1(v_predictionD,grid_v) .+ opC_p.Gy_b * vecb(v_predictionD,grid_v)
    #region cut-cell
    # for iLS in 1:nLS
    #     if !is_navier(bc_int[iLS]) && !is_navier_cl(bc_int[iLS])
    #         velocity_divergence .+= opC_p.Gx[iLS] * veci(u_predictionD,grid_u,iLS+1) .+ 
    #                 opC_p.Gy[iLS] * veci(v_predictionD,grid_v,iLS+1)
    #     end
    # end
    #endregion cut-cell

    #TODO replace




    #region check divergence
    #TODO function
    normalise_velocity_divergence = abs.(opC_p.AxT * vec1(u_predictionD,grid_u)) .+ abs.(opC_p.Gx_b * vecb(u_predictionD,grid_u)) .+
                                    abs.(opC_p.AyT * vec1(v_predictionD,grid_v)) .+ abs.(opC_p.Gy_b * vecb(v_predictionD,grid_v))
    # for iLS in 1:nLS
    #     if !is_navier(bc_int[iLS]) && !is_navier_cl(bc_int[iLS])
    #         normalise_velocity_divergence .+= abs.(opC_p.Gx[iLS] * veci(uD,grid_u,iLS+1)) .+ 
    #                 abs.(opC_p.Gy[iLS] * veci(vD,grid_v,iLS+1))
    #     end
    # end

  

    # PDI_status = @ccall "libpdi".PDI_multi_expose("rhs_uv"::Cstring,
    # "rhs_uv_len"::Cstring, length(rhs_uv)::Ref{Clonglong}, PDI_OUT::Cint,
    # "rhs_uv_1D"::Cstring, rhs_uv::Ptr{Cdouble}, PDI_OUT::Cint,
    # C_NULL::Ptr{Cvoid})::Cint

    # PDI_status = @ccall "libpdi".PDI_multi_expose("rhs_uv"::Cstring,
    # "rhs_uv_len"::Cstring, length(F_residual)::Ref{Clonglong}, PDI_OUT::Cint,
    # "rhs_uv_1D"::Cstring, F_residual::Ptr{Cdouble}, PDI_OUT::Cint,
    # C_NULL::Ptr{Cvoid})::Cint


    # PDI_status = @ccall "libpdi".PDI_multi_expose("check_divergence"::Cstring,
    # "nstep"::Cstring, num.current_iter::Ref{Clonglong}, PDI_OUT::Cint,
    # # "grad_x"::Cstring,grad_x::Ptr{Cdouble}, PDI_OUT::Cint,
    # # "grad_y"::Cstring, grad_y::Ptr{Cdouble}, PDI_OUT::Cint,
    # # "grad_u"::Cstring,grad_x::Ptr{Cdouble}, PDI_OUT::Cint,
    # # "grad_v"::Cstring, grad_y::Ptr{Cdouble}, PDI_OUT::Cint,
    # # "u_1D"::Cstring,  uD::Ptr{Cdouble}, PDI_OUT::Cint,
    # # "v_1D"::Cstring,  vD::Ptr{Cdouble}, PDI_OUT::Cint,
    # # "p_1D"::Cstring, pD::Ptr{Cdouble}, PDI_OUT::Cint,
    # "velocity_divergence"::Cstring, velocity_divergence::Ptr{Cdouble}, PDI_OUT::Cint,
    # "normalise_velocity_divergence"::Cstring, normalise_velocity_divergence::Ptr{Cdouble}, PDI_OUT::Cint,
    # # "max_abs_residual"::Cstring, max_abs_residual::Ref{Cdouble}, PDI_OUT::Cint,
    # # "max_abs_rhs"::Cstring, max_abs_rhs::Ref{Cdouble}, PDI_OUT::Cint,
    # C_NULL::Ptr{Cvoid})::Cint

    #endregion check divergence


    #region correction

    #region velocity correction, compute pressure with Poisson
    if num.pressure_velocity_coupling == 0
        # Poisson equation: source term
        # divergence of velocity / dt
        # vec1(rhs_phi,grid_p) .= idt .* velocity_divergence
        vec1(rhs_phi,grid_p) .= velocity_divergence

        if phase_change_currently_activated == 1
            vec1(rhs_phi,grid_p) .+=  vec( mass_transfer_rate * ( 1.0/num.rho1 - 1.0/num.rho2 ) )

            # PDI_status = @ccall "libpdi".PDI_multi_expose("rhs_uv_divergence"::Cstring,
            # "rhs_uv_divergence"::Cstring, vec1(rhs_phi,grid_p)::Ptr{Cdouble}, PDI_OUT::Cint,
            # C_NULL::Ptr{Cvoid})::Cint

            #TODO BC u v mass transfer outflow , so Neumann
        end



        #region needs to be corrected/documented for the signs, free surface pressure BC 
        
        # pres_free_suface = 0.0
        #TODO Marangoni
        #TODO phase change
        # diff_inv_rho = 1.0/rho1 - 1.0/rho2
        # jump_mass_transfer_rate = 0.0 #TODO

        #region cut-cell
        # if jump_mass_transfer_rate
        #     for iLS in 1:nLS
        #         if is_fs(bc_int[iLS])
        #             Smat = strain_rate(iLS, opC_u, opC_v, opC_p)
        #             S = Smat[1,1] * vec1(u_predictionD,grid_u) .+ Smat[1,2] * veci(u_predictionD,grid_u,iLS+1) .+
        #                 Smat[2,1] * vec1(v_predictionD,grid_v) .+ Smat[2,2] * veci(v_predictionD,grid_v,iLS+1)
        
        #             fs_mat = opC_p.HxT[iLS] * opC_p.Hx[iLS] .+ opC_p.HyT[iLS] * opC_p.Hy[iLS]
        #             veci(rhs_phi,grid_p,iLS+1) .= -2.0 .* mu1_over_rho1 .* S .+ Diagonal(diag(fs_mat)) * ( σ .* vec(grid_p.LS[iLS].κ) .- pres_free_suface .- diff_inv_rho * mass_transfer_rate ^ 2)
        #         end
        #     end
        # else
        #     for iLS in 1:nLS
        #         if is_fs(bc_int[iLS])
        #             Smat = strain_rate(iLS, opC_u, opC_v, opC_p)
        #             S = Smat[1,1] * vec1(u_predictionD,grid_u) .+ Smat[1,2] * veci(u_predictionD,grid_u,iLS+1) .+
        #                 Smat[2,1] * vec1(v_predictionD,grid_v) .+ Smat[2,2] * veci(v_predictionD,grid_v,iLS+1)

        #             fs_mat = opC_p.HxT[iLS] * opC_p.Hx[iLS] .+ opC_p.HyT[iLS] * opC_p.Hy[iLS]
        #             veci(rhs_phi,grid_p,iLS+1) .= -2.0 .* mu1_over_rho1 .* S .+ Diagonal(diag(fs_mat)) * ( σ .* vec(grid_p.LS[iLS].κ) .- pres_free_suface )
        #         end
        #     end
        # end
        #endregion cut-cell

        # Remove nullspace by adding small quantity to main diagonal
        if num.null_space == 0
            @inbounds @threads for i in 1:A_phi.m
                @inbounds A_phi[i,i] += 1e-10
            end
        end
        kill_dead_cells!(vec1(rhs_phi,grid_p), grid_p, geo[end])
        for iLS in 1:nLS
            kill_dead_cells!(veci(rhs_phi,grid_p,iLS+1), grid_p, geo[end])
        end
        # @time bicgstabl!(result_p, A_phi, rhs_phi, Pl = Diagonal(A_phi), log = true)

        #endregion needs to be corrected/documented for the signs, free surface pressure BC 


        # print("size pressure ",size(rhs_phi)," ",size(result_p)," ",size(A_phi))
        # Solve Poisson equation
        # \phi^{n+1}: result_p
        # @time result_p .= A_phi \ rhs_phi
        result_p = A_phi \ rhs_phi

        

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
        #     pres_grad_x .+= ∇ϕ_x
        #     pres_grad_y .+= ∇ϕ_y
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
            vec1(pD,grid_p) .+= vec(ϕ) #no timestep_n  since div u not rho1

        elseif num.prediction == "PmII" || num.prediction == "PmIIimposedpressure" || num.prediction == "PmIIimposedpressureBCincrement"
            # \nabla_h p^{n+1/2} = \nabla_h p^{n-1/2} + \nabla_h \phi^{n+1} - 
            #nu dt/2
            # \Delta t \nabla_h^2 \phi^{n+1} = \nabla_h \cdot \mathbf{u}^{*} \quad \text{in } \Omega
            

            print("\n max p vec1 ",maximum(vec1(pD,grid_p)))
            print("\n min p",minimum(ϕ)," max ",maximum(ϕ))
            # print("\n max p",minimum(num.mu_cin1./2 .* reshape(iM * velocity_divergence,grid_p))," ",maximum(num.mu_cin1./2 .* reshape(iM * velocity_divergence,grid_p)))

            # todo zero neumann bc !!!!
            # vec1(pD,grid_p) .+= vec(ϕ .- num.mu_cin1./2 .* reshape(iM * velocity_divergence,grid_p))
            # or,
            # for better readability
            # vec1(pD,grid_p) .= vec1(pD,grid_p) .+ vec(ϕ .- num.mu_cin1./2 .* reshape(iM * velocity_divergence,grid_p)) 
            vec1(pD,grid_p) .= vec1(pD,grid_p) .+ vec(ϕ .- mu_one_fluid./rho_one_fluid./2 .* reshape(iM * velocity_divergence,grid_p)) 
            
            # vec1(pD,grid_p) .+= vec(ϕ .- num.mu_cin1./2 .* reshape(iM * velocity_divergence,grid_p)) 

            print("\n max p vec1 ",maximum(vec1(pD,grid_p)))
            print("\n min p",minimum(ϕ)," max ",maximum(ϕ))
            # print("\n max p",minimum(num.mu_cin1./2 .* reshape(iM * velocity_divergence,grid_p))," ",maximum(num.mu_cin1./2 .* reshape(iM * velocity_divergence,grid_p)))


        elseif num.prediction == "PmIIimposedpressure_nodiv" ||
        num.prediction == "PmIIimposedpressure_nodiv_2" ||
        num.prediction == "PmIIimposedpressure_nodiv_3" ||
        num.prediction == "PmIIimposedpressure_nodiv_4"

            vec1(pD,grid_p) .= vec1(pD,grid_p) .+ vec(ϕ)

        elseif num.prediction == "PmIII"
            #TODO which zverage better here for higher order term mu_one_fluid./rho_one_fluid ?
            # \Delta t \nabla_h^2 \phi^{n+1} = \nabla_h \cdot \mathbf{u}^{*} \quad \text{in } \Omega
            vec1(pD,grid_p) .= vec(ϕ .- mu_one_fluid./rho_one_fluid ./2 .* reshape(iM * velocity_divergence,grid_p)) #no contribution from p^{n-1/2}
        
        elseif num.prediction == "Flower"  #occursin("Flower",num.prediction)
            # \nabla_h p^{n+1/2} = \nabla_h \phi^{n+1} # TODO: does not correspond to any formula in Brown 2001 ?
            vec1(pD,grid_p) .= vec(ϕ) #.- mu1_over_rho1 .* reshape(iM * velocity_divergence, grid_p))
        
        elseif  num.prediction == "Flowergravproj" 

            vec1(pD,grid_p) .= vec(ϕ) .- vec(rho_one_fluid .* g .* grid_p.y)

        elseif num.prediction == "testpressure"
            print("test pressure")

        else
            @error("wrong prediction method, does not exist")
        end

        # PDI_status = @ccall "libpdi".PDI_multi_expose("print_pressure_projection"::Cstring,
        # "u_1D"::Cstring, u_predictionD::Ptr{Cdouble}, PDI_OUT::Cint,
        # "v_1D"::Cstring, v_predictionD::Ptr{Cdouble}, PDI_OUT::Cint,
        # "p_1D"::Cstring, pD::Ptr{Cdouble}, PDI_OUT::Cint,
        # C_NULL::Ptr{Cvoid})::Cint

        #region interfacial pressure is overwritten
        # TODO check and document: for pressure-imposed Poiseuille, not sure
        #  
        if num.prediction == "PmIIimposedpressure" || num.prediction == "PmIIimposedpressure_nodiv" ||
        num.prediction == "PmIIimposedpressure_nodiv_2" || 
        num.prediction == "PmIIimposedpressure_nodiv_3" || 
        num.prediction == "PmIIimposedpressure_nodiv_4" || 

        num.prediction == "testpressure"


        elseif num.prediction == "PmIIimposedpressureBCincrement"
            #TODO reapply BC 
            #TODO
            #init p for Poiseuille : grad_x = 0
            #increment : grad_x=0
            # so grad_x still zero ? even with div ?
            tmp_vec_p = zeros(grid_p) 
            # tmp_vec_p .= 0.0
            get_height!(grid_p.LS[1],grid_p.ind,grid_p.dx,grid_p.dy,grid_p.LS[end].geoS,tmp_vec_p) #here tmp_vec_p solid #TODO geoS ??

            init_fields_multiple_levelsets!(num,pD,p,tmp_vec_p,BC_p,grid_p,num.pres_intfc,"pL")

        elseif  num.prediction == "Flowergravproj" 
            vecb(pD,grid_p) .= vecb(result_p,grid_p)  #TODO.+ ... vec(rho_one_fluid*g*grid_p.y)

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
        #     vec1(pD,grid_p) .= vec(p) .+ vec(ϕ) #.- mu1_over_rho1 .* iM * velocity_divergence
        #     vec2(pD,grid_p) .+= vec2(result_p,grid_p)
        #     vecb(pD,grid_p) .+= vecb(result_p,grid_p)
        #     p .= reshape(vec1(pD,grid_p), grid_p)
        # end

        # vec1(∇ϕ_x,grid_p) .*= irho1 
        # vec1(∇ϕ_y,grid_p) .*= irho1

        
        # u .= u_prediction .- timestep_n .* reshape(iMu * ∇ϕ_x, grid_u)
        # v .= v_prediction .- timestep_n .* reshape(iMv * ∇ϕ_y, grid_v)

        u .= u_prediction .- timestep_n .* reshape(iMu * ∇ϕ_x, grid_u) ./ rho_one_fluid_u
        v .= v_prediction .- timestep_n .* reshape(iMv * ∇ϕ_y, grid_v) ./ rho_one_fluid_v
        #region cut-cell
        # kill_dead_cells!(u, grid_u, geo_u[end])
        # kill_dead_cells!(v, grid_v, geo_v[end])
        #endregion cut-cell

        vec1(uD,grid_u) .= vec(u)
        vecb(uD,grid_u) .= vecb(u_predictionD,grid_u)
        vec1(vD,grid_v) .= vec(v)
        vecb(vD,grid_v) .= vecb(v_predictionD,grid_v)

        if num.prediction == "PmIIimposedpressure_nodiv_4"

            PDI_status = @ccall "libpdi".PDI_multi_expose("check_pressure_velocity_BC"::Cstring,
            "nstep"::Cstring, num.current_iter ::Ref{Clonglong}, PDI_OUT::Cint,
            "time"::Cstring, num.time::Ref{Cdouble}, PDI_OUT::Cint,
            "u_1D"::Cstring, uD::Ptr{Cdouble}, PDI_OUT::Cint,
            "v_1D"::Cstring, vD::Ptr{Cdouble}, PDI_OUT::Cint,
            "p_1D"::Cstring, pD::Ptr{Cdouble}, PDI_OUT::Cint,
            C_NULL::Ptr{Cvoid})::Cint


            # apply Neumann BC for u and v since the gradient is not available in the wall
            # vecb(uD,grid_u) .= vecb(u_predictionD,grid_u)
            # vecb(vD,grid_v) .= vecb(v_predictionD,grid_v)

            tmp_vec_u = zeros(grid_u) 
            # tmp_vec_u .= 0.0
            get_height!(grid_u.LS[1],grid_u.ind,grid_u.dx,grid_u.dy,grid_u.LS[end].geoS,tmp_vec_u) #here tmp_vec_p solid #TODO geoS ??

            init_fields_multiple_levelsets!(num,uD,u,tmp_vec_u,BC_u,grid_u,num.pres_intfc,"uL")

            tmp_vec_v = zeros(grid_v) 
            # tmp_vec_u .= 0.0
            get_height!(grid_v.LS[1],grid_v.ind,grid_v.dx,grid_v.dy,grid_v.LS[end].geoS,tmp_vec_v) #here tmp_vec_p solid #TODO geoS ??

            init_fields_multiple_levelsets!(num,vD,v,tmp_vec_v,BC_v,grid_v,num.pres_intfc,"vL")

            PDI_status = @ccall "libpdi".PDI_multi_expose("check_pressure_velocity_BC"::Cstring,
            "nstep"::Cstring, num.current_iter ::Ref{Clonglong}, PDI_OUT::Cint,
            "time"::Cstring, num.time::Ref{Cdouble}, PDI_OUT::Cint,
            "u_1D"::Cstring, uD::Ptr{Cdouble}, PDI_OUT::Cint,
            "v_1D"::Cstring, vD::Ptr{Cdouble}, PDI_OUT::Cint,
            "p_1D"::Cstring, pD::Ptr{Cdouble}, PDI_OUT::Cint,
            C_NULL::Ptr{Cvoid})::Cint

        end

        #region cut-cell
        # for iLS in 1:nLS
        #     if !is_navier(bc_int[iLS]) && !is_navier_cl(bc_int[iLS])
        #         veci(uD,grid_u,iLS+1) .= veci(u_predictionD,grid_u,iLS+1)
        #         veci(vD,grid_v,iLS+1) .= veci(v_predictionD,grid_v,iLS+1)
        #     end
        #     # if is_fs(bc_int[iLS])
        #     #     @inbounds for II in grid_u.ind.all_indices
        #     #         pII = lexicographic(II, grid_u.ny)
        #     #         if abs(veci(u_predictionD,grid_u,iLS+1)[pII]) > 1e-12
        #     #             veci(u_predictionD,grid_u,iLS+1)[pII] -= (timestep_n .* iMu * ∇ϕ_x)[pII]
        #     #         end
        #     #     end
        #     #     @inbounds for II in grid_v.ind.all_indices
        #     #         pII = lexicographic(II, grid_v.ny)
        #     #         if abs(veci(v_predictionD,grid_v,iLS+1)[pII]) > 1e-12
        #     #             veci(v_predictionD,grid_v,iLS+1)[pII] -= (timestep_n .* iMv * ∇ϕ_y)[pII]
        #     #         end
        #     #     end
        #     # end
        # end
        #endregion cut-cell


        grad_x = zeros(grid_u)
        grad_y = zeros(grid_v)

        #region cut-cell
        # compute_grad_T_x_T_y_array_u_v_capacities!(num, grid_p, grid_u, grid_v, opC_u, opC_v, grad_x, grad_y, pD)
        #endregion cut-cell

        compute_grad_T_x_T_y_array_u_v_capacities!(num, grid_p, grid_u, grid_v, opC_u, opC_v, grad_x, grad_y, pD)
        # compute_grad_T_x_T_y_array_u_v_capacities_one_fluid!(num, grid_p, grid_u, grid_v, opC_u, opC_v, grad_x, grad_y, pD)

        # TODO give another set of capacities (allocate twice or special cases everywhere?) 

        #TODO divergence level
    
        PDI_status = @ccall "libpdi".PDI_multi_expose("check_pressure_velocity_end"::Cstring,
        # "grad_x"::Cstring,grad_x::Ptr{Cdouble}, PDI_OUT::Cint,
        # "grad_y"::Cstring, grad_y::Ptr{Cdouble}, PDI_OUT::Cint,
        "grad_u"::Cstring,grad_x::Ptr{Cdouble}, PDI_OUT::Cint,
        "grad_v"::Cstring, grad_y::Ptr{Cdouble}, PDI_OUT::Cint,
        "u_1D"::Cstring, u_predictionD::Ptr{Cdouble}, PDI_OUT::Cint,
        "v_1D"::Cstring, v_predictionD::Ptr{Cdouble}, PDI_OUT::Cint,
        "p_1D"::Cstring, pD::Ptr{Cdouble}, PDI_OUT::Cint,
        C_NULL::Ptr{Cvoid})::Cint

        # print("\n num.mu_one_fluid_average " , num.mu_one_fluid_average)
        # display(mu_one_fluid)
    
    end #if num.pressure_velocity_coupling == 0
    # #region velocity correction, compute pressure with Poisson
    
    #endregion correction 



    #region end coupled

    if num.pressure_velocity_coupling == 3
        F_residual = similar(rhs_uv)

        # F_residual .= Auv*uvD_dummy .-rhs_uv
        F_residual .= Auv*uvD 

        printstyled(color=:red, @sprintf "\n A u\n")

        # print("\n Auv*uvD  ", F_residual)

        PDI_status = @ccall "libpdi".PDI_multi_expose("rhs_uv"::Cstring,
        "rhs_uv_len"::Cstring, length(F_residual)::Ref{Clonglong}, PDI_OUT::Cint,
        "rhs_uv_1D"::Cstring, F_residual::Ptr{Cdouble}, PDI_OUT::Cint,
        C_NULL::Ptr{Cvoid})::Cint

        printstyled(color=:red, @sprintf "\n A u\n")


        F_residual .-= rhs_uv

        # print("\n F_residual ", F_residual)


        # try
        #     if num.pressure_velocity_solver == 1
        #         printstyled(color=:red, @sprintf "\nBICGSTAB(2)\n")

        #         # bicgstabl!(x, A, b, l; kwargs...)
        #         bicgstabl!(uvD, Auv, rhs_uv, 2; tol=1.0e-6,maxiter=100) #residual normalised in Julia ? 
        #         # tol: (relative) stopping tolerance of the method;
        #         # verbose: print information during the iterations;
        #         # maxiter: maximum number of allowed iterations;
        #         # Pl and Pr: left and right precon
        #     else
        #         @time uvD .= Auv \ rhs_uv
        #     end

        # # try
        # #     @time uvD .= Auv \ rhs_uv
        # catch e
        #     printstyled(color=:red, @sprintf "\nError coupled pressure-velocity\n")
        #     # println(e)
        #     print(e)
        #     uvD .= Inf
        # printstyled(color=:green, @sprintf "\n NS solver error solve_one_fluid_NS_no_phase!\n")
        # println(e)
        # @error("NS solver error solve_one_fluid_NS_no_phase!")
        # num.status = 1
        # end
        # print("\n size uvD ",ntu, " ",size(uD)," ",size(vec1(uD,grid_u))," " ,size(uvD[1:nbu]))
        vec1(uD,grid_u) .= uvD[1:niu]
        vecb(uD,grid_u) .= uvD[niu+1:ntu]

        vec1(vD,grid_v) .= uvD[ntu+1:ntu+niv]
        vecb(vD,grid_v) .= uvD[ntu+1+niv:ntu+ntv]

        u .= reshape(vec1(uD,grid_u), grid_u)
        v .= reshape(vec1(vD,grid_v), grid_v)
    end

    # print("\n velocity")
    # display(v)

    # uD .= uvD[1:ntu]
    # vD .= uvD[ntu+1:ntu+ntv]

    if num.pressure_velocity_coupling ==3

       
        
        vec1(pD,grid_p) .= uvD[ntu+ntv+ntNavier+1:ntu+ntv+ntNavier+nip]

        #region correct pressure with top corner value (pressure known at a constant from NS)

        p .= reshape(vec1(pD,grid_p), grid_p)

        corr_p = p[grid_p.ny,1]
        print("\n corr_p ",corr_p)
        p .-= corr_p
        vec1(pD,grid_p) .= vec(p)
        #endregion correct pressure with top corner value (pressure known at a constant from NS)

    # else
    #     pD .= uvD[ntu+ntv+ntNavier+1:ntu+ntv+ntNavier+(num.nLS+1)*nip+nbp]
    end
    # print("uvD",uvD)


    grad_x = zeros(grid_u)
    grad_y = zeros(grid_v)
    compute_grad_T_x_T_y_array_u_v_capacities!(num, grid_p, grid_u, grid_v, opC_u, opC_v, grad_x, grad_y, pD)

    #TODO divergence level

    # Compute divergence of velocity
    velocity_divergence = opC_p.AxT * vec1(uD,grid_u) .+ opC_p.Gx_b * vecb(uD,grid_u) .+
                          opC_p.AyT * vec1(vD,grid_v) .+ opC_p.Gy_b * vecb(vD,grid_v)
    for iLS in 1:nLS
        if !is_navier(bc_int[iLS]) && !is_navier_cl(bc_int[iLS])
            velocity_divergence .+= opC_p.Gx[iLS] * veci(uD,grid_u,iLS+1) .+ 
                    opC_p.Gy[iLS] * veci(vD,grid_v,iLS+1)
        end
    end

    normalise_velocity_divergence = abs.(opC_p.AxT * vec1(uD,grid_u)) .+ abs.(opC_p.Gx_b * vecb(uD,grid_u)) .+
                                    abs.(opC_p.AyT * vec1(vD,grid_v)) .+ abs.(opC_p.Gy_b * vecb(vD,grid_v))
    for iLS in 1:nLS
        if !is_navier(bc_int[iLS]) && !is_navier_cl(bc_int[iLS])
            normalise_velocity_divergence .+= abs.(opC_p.Gx[iLS] * veci(uD,grid_u,iLS+1)) .+ 
                    abs.(opC_p.Gy[iLS] * veci(vD,grid_v,iLS+1))
        end
    end

    if num.pressure_velocity_coupling == 3
        F_residual .= Auv*uvD .-rhs_uv

        max_abs_residual = maximum(abs.(F_residual))
        max_abs_rhs = maximum(abs.(rhs_uv))

        # PDI_status = @ccall "libpdi".PDI_multi_expose("rhs_uv"::Cstring,
        # "rhs_uv_len"::Cstring, length(rhs_uv)::Ref{Clonglong}, PDI_OUT::Cint,
        # "rhs_uv_1D"::Cstring, rhs_uv::Ptr{Cdouble}, PDI_OUT::Cint,
        # C_NULL::Ptr{Cvoid})::Cint

        # PDI_status = @ccall "libpdi".PDI_multi_expose("rhs_uv"::Cstring,
        # "rhs_uv_len"::Cstring, length(F_residual)::Ref{Clonglong}, PDI_OUT::Cint,
        # "rhs_uv_1D"::Cstring, F_residual::Ptr{Cdouble}, PDI_OUT::Cint,
        # C_NULL::Ptr{Cvoid})::Cint

    
        PDI_status = @ccall "libpdi".PDI_multi_expose("check_coupled_solver_iteration"::Cstring,
        # "grad_x"::Cstring,grad_x::Ptr{Cdouble}, PDI_OUT::Cint,
        # "grad_y"::Cstring, grad_y::Ptr{Cdouble}, PDI_OUT::Cint,
        "grad_u"::Cstring,grad_x::Ptr{Cdouble}, PDI_OUT::Cint,
        "grad_v"::Cstring, grad_y::Ptr{Cdouble}, PDI_OUT::Cint,
        "u_1D"::Cstring,  uD::Ptr{Cdouble}, PDI_OUT::Cint,
        "v_1D"::Cstring,  vD::Ptr{Cdouble}, PDI_OUT::Cint,
        "p_1D"::Cstring, pD::Ptr{Cdouble}, PDI_OUT::Cint,
        "velocity_divergence"::Cstring, velocity_divergence::Ptr{Cdouble}, PDI_OUT::Cint,
        "normalise_velocity_divergence"::Cstring, normalise_velocity_divergence::Ptr{Cdouble}, PDI_OUT::Cint,
        "max_abs_residual"::Cstring, max_abs_residual::Ref{Cdouble}, PDI_OUT::Cint,
        "max_abs_rhs"::Cstring, max_abs_rhs::Ref{Cdouble}, PDI_OUT::Cint,
        C_NULL::Ptr{Cvoid})::Cint

        PDI_status = @ccall "libpdi".PDI_multi_expose("check_divergence"::Cstring,
        "nstep"::Cstring, num.current_iter::Ref{Clonglong}, PDI_OUT::Cint,
        "velocity_divergence"::Cstring, velocity_divergence::Ptr{Cdouble}, PDI_OUT::Cint,
        "normalise_velocity_divergence"::Cstring, normalise_velocity_divergence::Ptr{Cdouble}, PDI_OUT::Cint,
        C_NULL::Ptr{Cvoid})::Cint
        
    end


    kill_dead_cells!(vec1(uD,grid_u), grid_u, geo_u[end])
    u .= reshape(vec1(uD,grid_u), grid_u)
    kill_dead_cells!(vec1(vD,grid_v), grid_v, geo_v[end])
    v .= reshape(vec1(vD,grid_v), grid_v)
    
    PDI_status = @ccall "libpdi".PDI_multi_expose("print_velocity_prediction"::Cstring,
    "u_1D"::Cstring, uD::Ptr{Cdouble}, PDI_OUT::Cint,
    "v_1D"::Cstring, vD::Ptr{Cdouble}, PDI_OUT::Cint,
    "p_1D"::Cstring, pD::Ptr{Cdouble}, PDI_OUT::Cint,
    C_NULL::Ptr{Cvoid})::Cint


    PDI_status = @ccall "libpdi".PDI_multi_expose("grad_pres_y"::Cstring,
    # "grad_x"::Cstring,grad_x::Ptr{Cdouble}, PDI_OUT::Cint,
    # "grad_y"::Cstring, grad_y::Ptr{Cdouble}, PDI_OUT::Cint,
    # "grad_u"::Cstring,grad_x::Ptr{Cdouble}, PDI_OUT::Cint,
    "grad_pres_y"::Cstring, grad_y::Ptr{Cdouble}, PDI_OUT::Cint,
    # "u_1D"::Cstring, u_predictionD::Ptr{Cdouble}, PDI_OUT::Cint,
    # "v_1D"::Cstring, v_predictionD::Ptr{Cdouble}, PDI_OUT::Cint,
    # "p_1D"::Cstring, pD::Ptr{Cdouble}, PDI_OUT::Cint,
    C_NULL::Ptr{Cvoid})::Cint
    
    # NS_force_y = reshape(-grav_y .-pres_grad_y ./ vec(rho_one_fluid_v),grid_v)
   
   
       
   compute_grad_T_x_T_y_array_u_v_capacities_cell_integrated!(num, grid_p, grid_u, grid_v, opC_u, opC_v, grad_x, grad_y, pD)

    NS_force_y = reshape(-grav_y,grid_v)
    NS_force_y .-= grad_y ./ rho_one_fluid_v


    PDI_status = @ccall "libpdi".PDI_multi_expose("NS_force_y"::Cstring,
    # "grad_x"::Cstring,grad_x::Ptr{Cdouble}, PDI_OUT::Cint,
    # "grad_y"::Cstring, grad_y::Ptr{Cdouble}, PDI_OUT::Cint,
    # "grad_u"::Cstring,grad_x::Ptr{Cdouble}, PDI_OUT::Cint,
    "NS_force_y"::Cstring, NS_force_y::Ptr{Cdouble}, PDI_OUT::Cint,
    # "u_1D"::Cstring, u_predictionD::Ptr{Cdouble}, PDI_OUT::Cint,
    # "v_1D"::Cstring, v_predictionD::Ptr{Cdouble}, PDI_OUT::Cint,
    # "p_1D"::Cstring, pD::Ptr{Cdouble}, PDI_OUT::Cint,
    C_NULL::Ptr{Cvoid})::Cint


    # vec1(u_predictionD, grid_u) .= uvD[1:niu]
    # vecb(u_predictionD, grid_u) .= uvD[ntu-nbu+1:ntu]
    # kill_dead_cells!(vec1(u_predictionD,grid_u), grid_u, geo_u[end])
    # u_prediction .= reshape(vec1(u_predictionD,grid_u), grid_u)

    # vec1(v_predictionD, grid_v) .= uvD[ntu+1:ntu+niv]
    # vecb(v_predictionD, grid_v) .= uvD[ntu+ntv-nbv+1:ntu+ntv]
    # kill_dead_cells!(vec1(v_predictionD,grid_v), grid_v, geo_v[end])
    # v_prediction .= reshape(vec1(v_predictionD,grid_v), grid_v)

    nNav = 0
    _iLS = 1
    for iLS in 1:nLS
        if !is_navier(bc_int[iLS]) && !is_navier_cl(bc_int[iLS])
            veci(uD,grid_u,iLS+1) .= uvD[_iLS*niu+1:(_iLS+1)*niu]
            kill_dead_cells!(veci(uD,grid_u,iLS+1), grid_u, geo_u[end])

            veci(vD,grid_v,iLS+1) .= uvD[ntu+_iLS*niv+1:ntu+(_iLS+1)*niv]
            kill_dead_cells!(veci(vD,grid_v,iLS+1), grid_v, geo_v[end])
            _iLS += 1
        else
            @inbounds uT[nNav+1,:] .= vec(uvD[ntu+ntv+1+nNav*nip:ntu+ntv+(nNav+1)*nip])
            nNav += 1
        end
    end
    #endregion Navier


    PDI_status = @ccall "libpdi".PDI_multi_expose("print_velocity_prediction"::Cstring,
    "u_1D"::Cstring, uD::Ptr{Cdouble}, PDI_OUT::Cint,
    "v_1D"::Cstring, vD::Ptr{Cdouble}, PDI_OUT::Cint,
    "p_1D"::Cstring, pD::Ptr{Cdouble}, PDI_OUT::Cint,
    C_NULL::Ptr{Cvoid})::Cint

    II = CartesianIndex(div(grid_v.ny,2),1)
    pII = lexicographic(II,grid_v.ny)
    # print("\n A coeff ",Av[pII,:])
    # print("\nM ",opC_p.iMy.diag[pII])

    # print("\nAv[pII,:]./mu1_over_rho1 ",Av[pII,:]./mu1_over_rho1)

    grad_x = zeros(grid_u)
    grad_y = zeros(grid_v)
    compute_grad_T_x_T_y_array_u_v_capacities!(num, grid_p, grid_u, grid_v, opC_u, opC_v, grad_x, grad_y, pD)

    #TODO divergence level
  



    #endregion end coupled

    #region check_acceleration
    # printstyled(color=:red, @sprintf "\n Check acceleration \n")

    # rho_l = 1000.0
    # rho_g = 100.0
    # g = 9.81e-1
    # print("\n accel ", 2*(rho_l-rho_g)/(rho_l+2*rho_g)*g)

    # #v only, TODO interp and u
    # # acceleration = (v.-v0)/num.timestep_n
    # # display(acceleration)
    
    # printstyled(color=:red, @sprintf "\n Check acceleration \n")

    # print("\n opv.M ",opC_v.M[1,:])

    # # print("\n opv.M ",opC_v.M[pII,:])

    # print("\n opv.M ",opC_v.M[pII,:])



    #endregion check_acceleration


    PDI_status = @ccall "libpdi".PDI_multi_expose("check_NS_end"::Cstring,
    "u_1D"::Cstring, uD::Ptr{Cdouble}, PDI_OUT::Cint,
    "v_1D"::Cstring, vD::Ptr{Cdouble}, PDI_OUT::Cint,
    "p_1D"::Cstring, pD::Ptr{Cdouble}, PDI_OUT::Cint,
    "stop_simulation"::Cstring, num.stop_simulation::Ref{Clonglong}, PDI_OUT::Cint,
    C_NULL::Ptr{Cvoid})::Cint


     # Compute divergence of velocity
    velocity_divergence = opC_p.AxT * vec1(u_predictionD,grid_u) .+ opC_p.Gx_b * vecb(u_predictionD,grid_u) .+
          opC_p.AyT * vec1(v_predictionD,grid_v) .+ opC_p.Gy_b * vecb(v_predictionD,grid_v)
    #region cut-cell
    # for iLS in 1:nLS
    #     if !is_navier(bc_int[iLS]) && !is_navier_cl(bc_int[iLS])
    #         velocity_divergence .+= opC_p.Gx[iLS] * veci(u_predictionD,grid_u,iLS+1) .+ 
    #                 opC_p.Gy[iLS] * veci(v_predictionD,grid_v,iLS+1)
    #     end
    # end
    #endregion cut-cell

    #region check divergence
    #TODO function
    normalise_velocity_divergence = abs.(opC_p.AxT * vec1(u_predictionD,grid_u)) .+ abs.(opC_p.Gx_b * vecb(u_predictionD,grid_u)) .+
                                    abs.(opC_p.AyT * vec1(v_predictionD,grid_v)) .+ abs.(opC_p.Gy_b * vecb(v_predictionD,grid_v))
    # for iLS in 1:nLS
    #     if !is_navier(bc_int[iLS]) && !is_navier_cl(bc_int[iLS])
    #         normalise_velocity_divergence .+= abs.(opC_p.Gx[iLS] * veci(uD,grid_u,iLS+1)) .+ 
    #                 abs.(opC_p.Gy[iLS] * veci(vD,grid_v,iLS+1))
    #     end
    # end

    # PDI_status = @ccall "libpdi".PDI_multi_expose("check_divergence"::Cstring,
    # "nstep"::Cstring, num.current_iter::Ref{Clonglong}, PDI_OUT::Cint,
    # "velocity_divergence"::Cstring, velocity_divergence::Ptr{Cdouble}, PDI_OUT::Cint,
    # "normalise_velocity_divergence"::Cstring, normalise_velocity_divergence::Ptr{Cdouble}, PDI_OUT::Cint,
    # C_NULL::Ptr{Cvoid})::Cint

    #endregion check divergence


    return Lp, bc_Lp, bc_Lp_b, Lu, diffusion_LS_u, diffusion_border_u, Lv, diffusion_LS_v, diffusion_border_v, opC_p.M, opC_u.M, opC_v.M, Cui, Cvi
end
