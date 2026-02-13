function compute_mass_transfer_rate_main!(num, grid_p, grid_u, grid_v, op, phL, phS, BC_int, electrolysis, electrolysis_phase_change_case, 
    periodic_x, periodic_y, λ, Vmean, iLSpdi, mode_2d, show_every, 
    mass_transfer_rate, 
    mass_transfer_rate_vec1,
    mass_transfer_rate_vecb,mass_transfer_rate_veci, mass_transfer_rate_redistributed, tmp_vec_p, tmp_vec_p0, tmp_vec_p1,
    nb_gaz_acceptors, volume_fraction, 
    interface_length,total_interface_length)

    # print("\n electrolysis_phase_change_case ",electrolysis_phase_change_case)
    #TODO print case, quantity, ...

    if electrolysis && electrolysis_phase_change_case != "None"
        printstyled(color=:magenta, @sprintf "\n integrate_mass_transfer_rate_over_interface\n")
        
        # print("\n total_interface_length ",total_interface_length)

        # total_interface_length = compute_interface_length!(num, grid_p, 1, interface_length)

        print("\n total_interface_length ",total_interface_length)

        if total_interface_length == 0.0
            @error("\n error total_interface_length")
        end

        # @views integrate_mass_transfer_rate_over_interface(num,grid_p,op.opC_pL,phL.trans_scalD[:,1],mass_transfer_rate_vec1,mass_transfer_rate_vecb,mass_transfer_rate_veci,mass_transfer_rate)
        # @views integrate_mass_transfer_rate_over_interface_2(num,grid_p,op.opC_pL,phL.trans_scalD[:,1],mass_transfer_rate_vec1,mass_transfer_rate_vecb,mass_transfer_rate_veci,mass_transfer_rate)


        # check_JC()

        @views integrate_mass_transfer_rate_over_interface(num,grid_p,op.opC_pL,phL.trans_scalD[:,1],mass_transfer_rate_vec1,
        mass_transfer_rate_vecb,mass_transfer_rate_veci, tmp_vec_p, tmp_vec_p0, tmp_vec_p1, 
        mass_transfer_rate,num.index_phase_change) #1
        #here mass_transfer_rate is integrated on cell interface part
        PDI_status = @ccall "libpdi".PDI_multi_expose("check_mass_transfer_rate_NS"::Cstring,
        "mass_transfer_rate"::Cstring, mass_transfer_rate::Ptr{Cdouble}, PDI_OUT::Cint,
        C_NULL::Ptr{Cvoid})::Cint

        # print("\n sum mass flux all levelsets (walls and interfaces alike) ", sum(mass_transfer_rate),"\n ")
    end

    #    grid_p.LS[i].α  which is the angle of the outward point normal with respect to the horizontal axis

    for iLS in 1:num.nLS
        if is_stefan(BC_int[iLS])
            update_stefan_velocity(num, grid_p, iLS, grid_p.LS[iLS].u, phS.T, phL.T, periodic_x, periodic_y, λ, Vmean)
        elseif is_fs(BC_int[iLS]) || (occursin("levelset",electrolysis_phase_change_case) && iLS == num.iLSbubble)
            printstyled(color=:green, @sprintf "\n grid_p.V %.2e max abs(u) : %.2e max abs(v)%.2e\n" maximum(abs.(grid_p.V)) maximum(abs.(phL.u)) maximum(abs.(phL.v)))

            if electrolysis_phase_change_case!="none"    
                if occursin("levelset",electrolysis_phase_change_case)

                    printstyled(color=:magenta, @sprintf "\n phase-change for LS %.2i " iLS)

                    # plot_electrolysis_velocity!(num, grid_p, grid_p.LS, grid_p.V, TL, MIXED, periodic_x, periodic_y, concentration_scal_intfc)

                    # TODO send to PDI points and velocity for phase change like in plot_electrolysis_velocity!
                    
                    # Minus sign because normal points toward bubble and num.varnH2 for gaz, not liquid phase 

                    
                    if num.advection_LS_mode !=10    


                        PDI_status = @ccall "libpdi".PDI_multi_expose("check_mass_transfer_rate_NS"::Cstring,
                        "mass_transfer_rate"::Cstring, mass_transfer_rate::Ptr{Cdouble}, PDI_OUT::Cint,
                        "mass_transfer_rate_redistributed"::Cstring, mass_transfer_rate_redistributed::Ptr{Cdouble}, PDI_OUT::Cint,
                        "nb_gaz_acceptors"::Cstring, nb_gaz_acceptors::Ptr{Cdouble}, PDI_OUT::Cint,
                        C_NULL::Ptr{Cvoid})::Cint
                        
                        # display(nb_gaz_acceptors)

                        #use interface_length or temp_vec_p0

                        flower_status = compute_mass_transfer_rate!(num, grid_p, grid_u, grid_v, iLS, phL.uD, phL.vD, 
                        periodic_x, periodic_y, num.average_velocity, phL.trans_scalD[:,num.index_phase_change],phL.trans_scal[:,:,num.index_phase_change],
                        num.diffusion_coeff[num.index_phase_change],num.concentration0[num.index_phase_change],
                        electrolysis_phase_change_case,mass_transfer_rate, mass_transfer_rate_redistributed,
                        nb_gaz_acceptors,volume_fraction,interface_length,total_interface_length)

                        PDI_status = @ccall "libpdi".PDI_multi_expose("check_mass_transfer_rate_NS"::Cstring,
                        "mass_transfer_rate"::Cstring, mass_transfer_rate::Ptr{Cdouble}, PDI_OUT::Cint,
                        "mass_transfer_rate_redistributed"::Cstring, mass_transfer_rate_redistributed::Ptr{Cdouble}, PDI_OUT::Cint,
                        "nb_gaz_acceptors"::Cstring, nb_gaz_acceptors::Ptr{Cdouble}, PDI_OUT::Cint,
                        C_NULL::Ptr{Cvoid})::Cint
                        
                        
                        @ccall "libpdi".PDI_multi_expose("write_mass_transfer_rate_redistributed"::Cstring,
                        "mass_transfer_rate"::Cstring, mass_transfer_rate_redistributed::Ptr{Cdouble}, PDI_OUT::Cint,
                        "mass_transfer_rate_before_redistribution"::Cstring, mass_transfer_rate::Ptr{Cdouble}, PDI_OUT::Cint,   
                        "nb_gaz_acceptors"::Cstring, nb_gaz_acceptors::Ptr{Cdouble}, PDI_OUT::Cint,                               
                        C_NULL::Ptr{Cvoid})::Cvoid

                        #not here, need to use mass transfer rate without redistribution for contribution to advection
                        # if num.phase_change_method == 5
                        #     print("\n no redistribution")
                        # else
                        #     mass_transfer_rate .= mass_transfer_rate_redistributed
                        # end

                        flower_status = compute_phase_change_velocity_electrolysis!(num, grid_p, grid_u, grid_v, iLS, phL.uD, phL.vD, 
                        periodic_x, periodic_y, num.average_velocity, phL.trans_scalD[:,num.index_phase_change],phL.trans_scal[:,:,num.index_phase_change],
                        num.diffusion_coeff[num.index_phase_change],num.concentration0[num.index_phase_change],
                        electrolysis_phase_change_case,mass_transfer_rate, mass_transfer_rate_redistributed,
                        nb_gaz_acceptors,volume_fraction,interface_length)



                  

                        
                        PDI_status = @ccall "libpdi".PDI_multi_expose("check_advection"::Cstring,
                        "levelset_p"::Cstring, grid_p.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                        "advection_velocity_p"::Cstring, grid_p.V::Ptr{Cdouble}, PDI_OUT::Cint,
                        "advection_velocity_u"::Cstring, grid_u.V::Ptr{Cdouble}, PDI_OUT::Cint,
                        "advection_velocity_v"::Cstring, grid_v.V::Ptr{Cdouble}, PDI_OUT::Cint,                        
                        C_NULL::Ptr{Cvoid})::Cint

                        PDI_status = @ccall "libpdi".PDI_multi_expose("write_normal_phase_change_velocity"::Cstring,
                        "normal_phase_change_velocity"::Cstring, grid_p.V::Ptr{Cdouble}, PDI_OUT::Cint,                                     
                        C_NULL::Ptr{Cvoid})::Cint

                        #apply redistribution
                        # if num.phase_change_method == 5
                        #     print("\n no redistribution")

                        # else
                        #     mass_transfer_rate .= mass_transfer_rate_redistributed
                        # end
                        
                        mass_transfer_rate .= mass_transfer_rate_redistributed

                        @ccall "libpdi".PDI_multi_expose("write_mass_transfer_rate_only"::Cstring,
                        "mass_transfer_rate"::Cstring, mass_transfer_rate::Ptr{Cdouble}, PDI_OUT::Cint,
                        # "mass_transfer_rate_bulk"::Cstring, mass_transfer_rate_vec1_2::Ptr{Cdouble}, PDI_OUT::Cint,
                        # "mass_transfer_rate_border"::Cstring, mass_transfer_rate_vecb_2::Ptr{Cdouble}, PDI_OUT::Cint,
                        # "mass_transfer_rate_intfc"::Cstring, mass_transfer_rate_veci_2::Ptr{Cdouble}, PDI_OUT::Cint,
                        C_NULL::Ptr{Cvoid})::Cvoid


                        PDI_status = @ccall "libpdi".PDI_multi_expose("check_mass_transfer_rate_NS"::Cstring,
                        "mass_transfer_rate"::Cstring, mass_transfer_rate::Ptr{Cdouble}, PDI_OUT::Cint,
                        C_NULL::Ptr{Cvoid})::Cint



                        if flower_status !=0
                            printstyled(color=:red, @sprintf "\n Stopping simulation %.3i " flower_status)
                            # return status
                        end
                    end

                        # # iLS = 1
                        # # intfc_length = 0.0
                        # # @inbounds @threads for II in grid_p.LS[iLS].MIXED
                        # #     intfc_length += 
                        # # end


                        # printstyled(color=:green, @sprintf "\n pi*R %.2e len : %.2e \n" π*num.R intfc_length)

                        # #TODO u-vphase change

                        # #TODO check velocity
                        # @inbounds @threads for II in grid_p.LS[iLS].MIXED
                        #     grid_p.V[II] = sum(mass_transfer_rate) * num.diffusion_coeff[num.index_phase_change] *(1.0/num.rho2-1.0/num.rho1).*num.diffusion_coeff[num.index_phase_change].*num.MWH2
                        # end


                    if num.mass_transfer_rate == 0
                        num.varnH2 = num.sum_mass_transfer_rate * num.diffusion_coeff[num.index_phase_change] 

                        num.new_nH2 = num.nH2 + num.varnH2 * num.timestep_n

                        print("\n varn ",num.varnH2 ," dt ", num.timestep_n," dn ",num.varnH2 * num.timestep_n, " sum ", num.sum_mass_transfer_rate)
                        printstyled(color=:green, @sprintf "\n it %.5i Mole: %.2e dn %.2e new num.nH2 %.2e \n" num.current_iter num.nH2 num.varnH2*num.timestep_n num.new_nH2)

                        if num.varnH2 < 0.0 
                            # print(@sprintf "error num.nH2 %.2e dnH2 %.2e new num.nH2 %.2e\n" num.nH2-num.varnH2*num.timestep_n num.varnH2*num.timestep_n num.nH2 )
                            @error ("error num.nH2")
                            crashed = true
                            num.new_nH2 = num.nH2
                            print("wrong num.nH2 ")
                            # println(@sprintf "\n CRASHED after %d iterations \n" num.current_iter)
                            # return status
                        else
                            num.nH2 = num.new_nH2
                        end

                    end #num.mass_transfer_rate == 0


                end
                # update_free_surface_velocity(num, grid_u, grid_v, iLS, phL.uD, phL.vD, periodic_x, periodic_y)
                printstyled(color=:green, @sprintf "\n grid_p.V %.2e max abs(u) : %.2e max abs(v)%.2e\n" maximum(abs.(grid_p.V)) maximum(abs.(phL.u)) maximum(abs.(phL.v)))
                
                printstyled(color=:green, @sprintf "\n grid_p.V %.2e dx : %.2e CFL %.2e\n" maximum(abs.(grid_p.V)) grid_p.dx[1,1] maximum(abs.(grid_p.V))*num.timestep_n/grid_p.dx[1,1])


            else
                printstyled(color=:magenta, @sprintf "\n update_free_surface_velocity")
                
                update_free_surface_velocity(num, grid_u, grid_v, iLS, phL.uD, phL.vD, periodic_x, periodic_y)
            end

            # printstyled(color=:magenta, @sprintf "\n update_free_surface_velocity")
            # #TODO
            # update_free_surface_velocity(num, grid_u, grid_v, 1, phL.uD, phL.vD, periodic_x, periodic_y)


        
        elseif (electrolysis && occursin("Khalighi",electrolysis_phase_change_case))

            if ((num.current_iter-1)%show_every == 0) 
                # print_electrolysis_statistics(num,grid_p,phL)
                PDI_status = @ccall "libpdi".PDI_multi_expose("print_variables"::Cstring,
                "nstep"::Cstring, num.current_iter::Ref{Clonglong}, PDI_OUT::Cint,
                "time"::Cstring, num.current_time::Ref{Cdouble}, PDI_OUT::Cint,
                "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
                "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
                "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
                "levelset_p"::Cstring, grid_p.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                "levelset_u"::Cstring, grid_u.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                "levelset_v"::Cstring, grid_v.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                "trans_scal_1DT"::Cstring, phL.trans_scalD'::Ptr{Cdouble}, PDI_OUT::Cint,
                "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,                       
                C_NULL::Ptr{Cvoid})::Cint
            end

            num.previous_radius = num.current_radius

            # Minus sign because normal points toward bubble and num.varnH2 for gaz, not liquid phase 
            num.varnH2 =  sum(mass_transfer_rate) * num.diffusion_coeff[num.index_phase_change] 

            #TODO mode_2d==0 flux corresponds to cylinder of length 1
            #2D cylinder reference length
            if mode_2d == 1
                num.varnH2 .*= num.ref_thickness_2d
            end

            #Pliquid is the average value of p over the bubble interface plus the ambient operating pressure (P).
            p_liq= num.pres0 + mean(veci(phL.pD,grid_p,2)) #TODO here one bubble
            # p_g=p_liq + 2 * num.σ / num.current_radius #3D
            p_g=p_liq + num.σ / num.current_radius #2D

            num.new_nH2 = num.nH2 + num.varnH2 * num.timestep_n

            
            printstyled(color=:green, @sprintf "\n it %.5i Mole: %.2e dn %.2e new num.nH2 %.2e \n" num.current_iter num.nH2 num.varnH2*num.timestep_n num.new_nH2)

            if num.varnH2 < 0.0 
                # print(@sprintf "error num.nH2 %.2e dnH2 %.2e new num.nH2 %.2e\n" num.nH2-num.varnH2*num.timestep_n num.varnH2*num.timestep_n num.nH2 )
                @error ("error num.nH2")
                crashed = true
                num.new_nH2 = num.nH2
                print("wrong num.nH2 ")
                # println(@sprintf "\n CRASHED after %d iterations \n" num.current_iter)
                # return status
            end

            if occursin("Khalighi_no_update",electrolysis_phase_change_case)

            else
                num.nH2 = num.new_nH2
            end
            
            #TODO using num.temperature0
            if mode_2d == 0
                num.current_radius = cbrt(3.0 * num.nH2 * num.Ru * num.temperature0/( 4.0 * pi * p_g) )
            elseif mode_2d == 1
                num.current_radius = sqrt(num.nH2 * num.Ru * num.temperature0/( pi * p_g * num.ref_thickness_2d) )
            elseif mode_2d == 2
                num.current_radius = sqrt(num.nH2/(num.concentration0[num.index_phase_change] * pi))
            elseif mode_2d == 3
                num.current_radius = sqrt(2*num.nH2/(num.concentration0[num.index_phase_change] * pi))
            elseif mode_2d == 4 #TODO
                num.current_radius = sqrt(num.nH2 * num.Ru * num.temperature0/( pi * p_g) )
            end

            printstyled(color=:green, @sprintf "\n radius num.CFL: %.2e \n" (num.current_radius-num.previous_radius)/(num.L0/grid_p.nx))

            if (num.current_radius-num.previous_radius)/(num.L0/grid_p.nx) > 0.5
                printstyled(color=:red, @sprintf "\n radius num.CFL: %.2e \n" (num.current_radius-num.previous_radius)/(num.L0/grid_p.nx))
                @error ("num.CFL radius")
                crashed = true
                # return status
            end

            
            printstyled(color=:cyan, @sprintf "\n div(0,grad): %.5i %.2e %.2e %.2e %.2e\n" grid_p.nx num.timestep_n num.L0/grid_p.nx (num.current_radius-num.previous_radius)/(num.L0/grid_p.nx) sum(mass_transfer_rate))
            
            printstyled(color=:green, @sprintf "\n num.n(H2): %.2e added %.2e old R %.2e new R %.2e \n" num.nH2 num.varnH2*num.timestep_n num.previous_radius num.current_radius)
            printstyled(color=:green, @sprintf "\n p0: %.2e p_liq %.2e p_lapl %.2e \n" num.pres0 p_liq p_g)

            if mode_2d == 3
                grid_p.LS[1].u .= sqrt.((grid_p.x.- num.xcoord).^ 2 + (grid_p.y .- num.ycoord) .^ 2) - (num.current_radius) * ones(grid_p.ny, grid_p.nx)                  
            else
                grid_p.LS[1].u .= sqrt.((grid_p.x .- num.xcoord .- num.current_radius .+ num.R ).^ 2 + (grid_p.y .- num.ycoord) .^ 2) - (num.current_radius) * ones(grid_p.ny, grid_p.nx)
            end
            # init_franck!(grid_p, TL, R, num.T_inf, 0)
            # u

        elseif (electrolysis && electrolysis_phase_change_case == "imposed_radius")

            #num.CFL 0.5
            num.current_radius = num.current_radius + grid_p.dx[1,1]/2

            grid_p.LS[1].u .= sqrt.((grid_p.x.- num.xcoord).^ 2 + (grid_p.y .- num.ycoord) .^ 2) - (num.current_radius) * ones(grid_p.ny, grid_p.nx)                  

        elseif (electrolysis && electrolysis_phase_change_case == "imposed_radius4")

            #num.CFL 0.5
            num.current_radius = num.current_radius + grid_p.dx[1,1]/4

            grid_p.LS[1].u .= sqrt.((grid_p.x.- num.xcoord).^ 2 + (grid_p.y .- num.ycoord) .^ 2) - (num.current_radius) * ones(grid_p.ny, grid_p.nx)                  


        end #phase change

    end #iLS

    # return status

end


# """
# to check mass conservation
# """
# function compute_conservation_mass(num,phL, grid_u, grid_v, op)

#     #TODO mutliple levelsets 1 or end
#     # cf bc_matrix_borders!(grid_p, opC_u.Gx_b, opC_v.Gy_b, opC_p.Gx_b, opC_p.Gy_b, geo.dcap)

#     conservation = 0.0

#     # conservation += -dot(vecb_L(phL.uD, grid_u), grid_u.LS[1].geoL.dcap[:, 1, 1])
#     # conservation += dot(vecb_R(phL.uD, grid_u), grid_u.LS[1].geoL.dcap[:, grid_u.nx, 3])  # right capacity: u
#     # conservation += -dot(vecb_B(phL.vD, grid_v), grid_v.LS[1].geoL.dcap[1, :, 2])  # bottom capacity: v
#     # conservation += dot(vecb_T(phL.vD, grid_v), grid_v.LS[1].geoL.dcap[grid_v.ny, :, 4])  # top capacity: v

#     # print("\n dummy")
#     # vecb_L(phL.uD, grid_u) .= - 1.0

#     # vecb_B(phL.vD, grid_v) .= -1.0

#     # vecb_R(phL.uD, grid_u) .= 1.0
#     # vecb_T(phL.vD, grid_v) .=1.0


#     conservation += -dot(left_border_view(op.opC_pL.Gx_b,grid_u) , vecb_L(phL.uD, grid_u))
#     conservation += -dot(bottom_border_view(op.opC_pL.Gy_b,grid_v) , vecb_B(phL.vD, grid_v))
#     conservation += dot(right_border_view(op.opC_pL.Gx_b,grid_u) , vecb_R(phL.uD, grid_u))
#     conservation += dot(top_border_view(op.opC_pL.Gy_b ,grid_v), vecb_T(phL.vD, grid_v))

#     # conservation *= num.rho1 * num.timestep_n

#     return conservation
# end


  function compute_conservation_mass(num,phL, grid_p, grid_u, grid_v,rho_one_fluid)

        #TODO mutliple levelsets 1 or end
        # cf bc_matrix_borders!(grid_p, opC_u.Gx_b, opC_v.Gy_b, opC_p.Gx_b, opC_p.Gy_b, geo.dcap)

        conservation = 0.0
        # print("\n compute_conservation_mass")
        # print("\n conservation ",conservation)

        # border_values = extract_border_capacities_1D(grid_p, dcap)  grid_p.LS[1].geoL.dcap
        # border_values = extract_border_capacities_1D_one_fluid(grid_p) 
        
        # vecx = grid_p.dx * rho_one_fluid
        # vecy = grid_p.dy * rho_one_fluid
        
        #TODO finer info loc face for rho
        vecx = grid_p.dx .* rho_one_fluid 
        vecy = grid_p.dy .* rho_one_fluid

        # vecx = grid_p.dx 
        # vecy = grid_p.dy 


        # print("\n vecx ",vecx)
        # print("\n vecy ",vecy)

        border_values = extract_vec_1D_one_fluid(grid_p,vecx,vecy)

        # TODO more precise loc face ...
        # rho_one_fluid

        # print("\n",border_values)
        # print("\n",vecb_L(phL.uD, grid_u))
        # print("\n",vecb_B(phL.vD, grid_v))
        # print("\n",vecb_R(phL.uD, grid_u))
        # print("\n",vecb_T(phL.vD, grid_v))

        conservation += dot(left_border_view(border_values,grid_p), vecb_L(phL.uD, grid_u))
        # print("\n conservation ",conservation)
        conservation += dot(bottom_border_view(border_values,grid_p) , vecb_B(phL.vD, grid_v))
        # print("\n conservation ",conservation)
        conservation += dot(right_border_view(border_values,grid_p) , vecb_R(phL.uD, grid_u))
        # print("\n conservation ",conservation)
        conservation += dot(top_border_view(border_values ,grid_p), vecb_T(phL.vD, grid_v))
        # print("\n conservation ",conservation)

        # conservation *= num.rho1 * num.timestep_n
        conservation *=  num.timestep_n
        # conservation *= num.rho2 * num.timestep_n


        # print("\n conservation ",conservation)

        # print("\n compute_conservation_mass")


    return conservation
    end

    

    function extract_border_capacities_1D(grid, dcap)
        @unpack nx, ny, ind = grid
        @unpack b_left, b_bottom, b_right, b_top = ind

        # Initialize a 1D array to store the border values
        border_values = zeros(2 * (nx + ny))

        @inbounds @threads for i in 1:ny
            II = CartesianIndex(i, 1)
            @inbounds A1 = dcap[II, 1]
            @inbounds border_values[i] = -A1  # Left border
        end

        @inbounds @threads for i in 1:nx
            II = CartesianIndex(1, i)
            @inbounds A2 = dcap[II, 2]
            @inbounds border_values[ny + i] = -A2  # Bottom border
        end

        @inbounds @threads for i in 1:ny
            II = CartesianIndex(i, nx)
            @inbounds A3 = dcap[II, 3]
            @inbounds border_values[ny + nx + i] = A3  # Right border
        end

        @inbounds @threads for i in 1:nx
            II = CartesianIndex(ny, i)
            @inbounds A4 = dcap[II, 4]
            @inbounds border_values[ny + nx + ny + i] = A4  # Top border
        end

        return border_values
    end

    function extract_border_capacities_1D_one_fluid(grid_p)
        @unpack nx, ny = grid_p
        @unpack dx, dy = grid_p

        # Initialize a 1D array to store the border values
        border_values = zeros(2 * (nx + ny))

        # Compute border values for left, bottom, right, top
        @inbounds @threads for i in 1:ny
            border_values[i] = -dx[i, 1]  # Left border
        end

        @inbounds @threads for i in 1:nx
            border_values[ny + i] = -dy[1, i]  # Bottom border
        end

        @inbounds @threads for i in 1:ny
            border_values[ny + nx + i] = dx[i, nx]  # Right border
        end

        @inbounds @threads for i in 1:nx
            border_values[ny + nx + ny + i] = dy[ny, i]  # Top border
        end

        return border_values
    end

    function extract_vec_1D_one_fluid(grid_p,vecx,vecy)
        @unpack nx, ny = grid_p
        # @unpack dx, dy = grid_p

        # Initialize a 1D array to store the border values
        border_values = zeros(2 * (nx + ny))

        # Compute border values for left, bottom, right, top
        @inbounds @threads for i in 1:ny
            border_values[i] = -vecx[i, 1]  # Left border
        end

        @inbounds @threads for i in 1:nx
            border_values[ny + i] = -vecy[1, i]  # Bottom border
        end

        @inbounds @threads for i in 1:ny
            border_values[ny + nx + i] = vecx[i, nx]  # Right border
        end

        @inbounds @threads for i in 1:nx
            border_values[ny + nx + ny + i] = vecy[ny, i]  # Top border
        end

        return border_values
    end