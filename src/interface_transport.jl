
"""
advect with grid_u.V and grid_v.V
"""
function advection_u_and_v(grid_p, grid_u, grid_v, iLS, θ_out, num, BC_int, BC_u, rhs_LS, periodic_x, periodic_y)
    rhs_LS .= 0.0
    grid_p.LS[iLS].A.nzval .= 0.0
    grid_p.LS[iLS].B.nzval .= 0.0
    IIOE!(grid_p, grid_u, grid_v, grid_p.LS[iLS].A, grid_p.LS[iLS].B, θ_out, num.timestep_n, periodic_x, periodic_y)
    BC_LS_interior!(num, grid_p, grid_u, grid_v, iLS, grid_p.LS[iLS].A, grid_p.LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
    BC_LS!(grid_p, grid_p.LS[iLS].u, grid_p.LS[iLS].A, grid_p.LS[iLS].B, rhs_LS, BC_u)
    utmp .= reshape(gmres(grid_p.LS[iLS].A, grid_p.LS[iLS].B * vec(grid_p.LS[iLS].u) .+ rhs_LS), grid_p)

    rhs_LS .= 0.0
    S2IIOE!(grid_p, grid_u, grid_v, grid_p.LS[iLS].A, grid_p.LS[iLS].B, utmp, grid_p.LS[iLS].u, θ_out, num.timestep_n, periodic_x, periodic_y)
    BC_LS_interior!(num, grid_p, grid_u, grid_v, iLS, grid_p.LS[iLS].A, grid_p.LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
    BC_LS!(grid_p, grid_p.LS[iLS].u, grid_p.LS[iLS].A, grid_p.LS[iLS].B, rhs_LS, BC_u)
    grid_p.LS[iLS].u .= reshape(gmres(grid_p.LS[iLS].A, grid_p.LS[iLS].B * vec(grid_p.LS[iLS].u) .+ rhs_LS), grid_p)
end



"""
compute normal component of velocity 
"""
function compute_normal_component_of_velocity(num,grid_u, grid_v, u, v, grid_p, iLS, II)
    cap1 = grid_u.LS[iLS].geoL.cap[II,5]
    cap3 = grid_u.LS[iLS].geoL.cap[δx⁺(II),5]
    # interp_velocity_p_grid_x = (u[II] * cap1 + u[δx⁺(II)] * cap3) / (cap1 + cap3 + num.epsilon_vol)
     denominator_x = cap1 + cap3
    if denominator_x > 0.0
        interp_velocity_p_grid_x = (u[II] * cap1 + u[δx⁺(II)] * cap3) / denominator_x
    else
        interp_velocity_p_grid_x = 0.0 #TODO
        print("MIXED cell with zero denominator")
    end

    cap2 = grid_v.LS[iLS].geoL.cap[II,5]
    cap4 = grid_v.LS[iLS].geoL.cap[δy⁺(II),5]
    # interp_velocity_p_grid_y = (v[II] * cap2 + v[δy⁺(II)] * cap4) / (cap2 + cap4 + num.epsilon_vol)
    denominator_y = cap2 + cap4
    if denominator_y > 0.0
        interp_velocity_p_grid_y = (v[II] * cap2 + v[δy⁺(II)] * cap4) / denominator_y
    else
        interp_velocity_p_grid_y = 0.0 #TODO
    end

    # Calculate magnitude and angle
    interp_velocity_p_grid = sqrt(interp_velocity_p_grid_x^2 + interp_velocity_p_grid_y^2)
    β = atan(interp_velocity_p_grid_y, interp_velocity_p_grid_x)
    if grid_p.LS[iLS].α[II] > 0.0 && β < 0.0
        β += 2π
    end
    if grid_p.LS[iLS].α[II] < 0.0 && β > 0.0
        β -= 2π
    end

    normal_comp = interp_velocity_p_grid * cos(β - grid_p.LS[iLS].α[II])

    return normal_comp

end


# """
# Project velocities to the normal and use advection scheme for advection just in the normal direction
#     not efficient, unecessary allocation
# """
# function compute_normal_velocity(grid_p, grid_u, grid_v, iLS, tmpVx, tmpVy)
#    tmpVx .=0.0
#    tmpVy .=0.0

#     # grid_p.V .= 0.0

#     @inbounds @threads for II in grid_p.LS[iLS].MIXED
#         cap1 = grid_u.LS[iLS].geoL.cap[II,5]
#         cap3 = grid_u.LS[iLS].geoL.cap[δx⁺(II),5]
#         # tmpVx[II] = (grid_u.V[II] * cap1 + grid_u.V[δx⁺(II)] * cap3) / (cap1 + cap3 + eps(0.01))
#         denominator_x = cap1 + cap3
#         if denominator_x > 0.0
#             tmpVx[II] = (grid_u.V[II] * cap1 + grid_u.V[δx⁺(II)] * cap3) / denominator_x
#         else
#             tmpVx[II] = 0.0 #TODO
#         end

#         cap2 = grid_v.LS[iLS].geoL.cap[II,5]
#         cap4 = grid_v.LS[iLS].geoL.cap[δy⁺(II),5]
#         # tmpVy[II] = (grid_v.V[II] * cap2 + grid_v.V[δy⁺(II)] * cap4) / (cap2 + cap4 + eps(0.01))
#         denominator_y = cap2 + cap4
#         if denominator_y > 0.0
#             tmpVy[II] = (grid_v.V[II] * cap2 + grid_v.V[δy⁺(II)] * cap4) / denominator_y
#         else
#             tmpVy[II] = 0.0 #TODO
#         end
        
#         # Calculate magnitude and angle
#         tmpV = sqrt(tmpVx[II]^2 + tmpVy[II]^2)
#         β = atan(tmpVy[II], tmpVx[II])

#         # Adjust angle based on interface orientation
#         if grid_p.LS[iLS].α[II] > 0.0 && β < 0.0
#             β += 2π
#         elseif grid_p.LS[iLS].α[II] < 0.0 && β > 0.0
#             β -= 2π
#         end

#         # Project velocity onto the normal direction and add it to grid_p.V (added because may have phase change contribution too)
#         grid_p.V[II] += tmpV * cos(β - grid_p.LS[iLS].α[II])
#     end

# end


"""
select method to compute the interfacial velocity for interface transport
advection_LS_mode:
* 12 : u and v
* 13 : 

"""
function select_advection!(num, grid_p, BC_int, BC_u, grid_u, grid_v, CFL_sc, periodic_x, periodic_y, 
    θ_out, rhs_LS, utmp, electrolysis_phase_change_case, mass_transfer_rate, levelset_1D,volume_fraction,
    tmp_vec_p,tmp_vec_p0,
    tmp_vec_u,tmp_vec_v,tmp_vec_u0,tmp_vec_v0,tmp_vec_u1,tmp_vec_v1,
    op,
    phL::Phase{Float64},u_extended=nothing, v_extended=nothing)

    print("\n select_advection")

   #In the current implementations, the first cell is not solved: levelset not advected
    for (iLS, bc) in enumerate(BC_int)
        print("\n advecting LS number ",iLS," num.advection_LS_mode ",num.advection_LS_mode )

        # if no mixed cells, do not advect
        if isempty(grid_p.LS[iLS].MIXED)
            print("\nNo advection")
            continue
        end

        if is_stefan(bc) #normal advection
            IIOE_normal!(grid_p, grid_p.LS[iLS].A, grid_p.LS[iLS].B, grid_p.LS[iLS].u, grid_p.V, CFL_sc, periodic_x, periodic_y)
            grid_p.LS[iLS].u .= reshape(gmres(grid_p.LS[iLS].A, grid_p.LS[iLS].B * vec(grid_p.LS[iLS].u)), grid_p)

        # elseif is_fs(bc) || (occursin("levelset",electrolysis_phase_change_case) && iLS == num.iLSbubble)
        else 
            if num.advection_LS_mode == 0
                advection_u_and_v(grid_p, grid_u, grid_v, iLS, θ_out, num, BC_int, BC_u, rhs_LS, periodic_x, periodic_y)

            elseif num.advection_LS_mode == 1

                # Project velocities to the normal and use advection scheme for advection just
                # in the normal direction
                grid_p.V .= 0.0
                @inbounds @threads for II in grid_p.LS[iLS].MIXED
                    grid_p.V[II] += compute_normal_component_of_velocity(num,grid_u, grid_v, grid_u.V, grid_v.V, grid_p, iLS, II)
                end

                i_ext, l_ext, b_ext, r_ext, t_ext = indices_extension(grid_p, grid_p.LS[iLS], grid_p.ind.inside, periodic_x, periodic_y)
                field_extension!(grid_p, grid_p.LS[iLS].u, grid_p.V, i_ext, l_ext, b_ext, r_ext, t_ext, num.NB, periodic_x, periodic_y)

                rhs_LS .= 0.0
                IIOE_normal!(grid_p, grid_p.LS[iLS].A, grid_p.LS[iLS].B, grid_p.LS[iLS].u, grid_p.V, CFL_sc, periodic_x, periodic_y)
                BC_LS!(grid_p, grid_p.LS[iLS].u, grid_p.LS[iLS].A, grid_p.LS[iLS].B, rhs_LS, BC_u)
                BC_LS_interior!(num, grid_p, iLS, grid_p.LS[iLS].A, grid_p.LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
                grid_p.LS[iLS].u .= reshape(gmres(grid_p.LS[iLS].A, grid_p.LS[iLS].B * vec(grid_p.LS[iLS].u) .+ rhs_LS), grid_p)

                # Impose contact angle if a wall is present
                rhs_LS .= 0.0
                grid_p.LS[iLS].A.nzval .= 0.0
                grid_p.LS[iLS].B.nzval .= 0.0
                for II in grid_p.ind.all_indices
                    pII = lexicographic(II, grid_p.ny)
                    grid_p.LS[iLS].A[pII,pII] = 1.0
                    grid_p.LS[iLS].B[pII,pII] = 1.0
                end
                BC_LS_interior!(num, grid_p, iLS, grid_p.LS[iLS].A, grid_p.LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
                grid_p.LS[iLS].u .= reshape(gmres(grid_p.LS[iLS].A, grid_p.LS[iLS].B * vec(grid_p.LS[iLS].u) .+ rhs_LS), grid_p)

            elseif num.advection_LS_mode == 2

                update_radius_from_contact_line(num,grid_p, grid_p.LS[iLS].u, BC_u)


                printstyled(color=:green, @sprintf "\n grid_p p u v max : %.2e %.2e %.2e\n" maximum(abs.(grid_p.V[grid_p.LS[iLS].MIXED])) maximum(abs.(grid_u.V[grid_p.LS[iLS].MIXED])) maximum(abs.(grid_v.V[grid_v.LS[iLS].MIXED])))

                IIOE_normal!(grid_p, grid_p.LS[iLS].A, grid_p.LS[iLS].B, grid_p.LS[iLS].u, grid_p.V, CFL_sc, periodic_x, periodic_y)
                
                grid_p.LS[iLS].u .= reshape(gmres(grid_p.LS[iLS].A, grid_p.LS[iLS].B * vec(grid_p.LS[iLS].u)), grid_p)

                update_radius_from_contact_line(num,grid_p, grid_p.LS[iLS].u, BC_u)


            elseif num.advection_LS_mode == 3

                # from scalar grid_p, normal to u and v
                #TODO 
                print("\n setting 0 u and v velocities \n")
                grid_u.V .= 0
                grid_v.V .= 0

                interpolate_scalar!(grid_p, grid_u, grid_v, grid_p.V, grid_u.V, grid_v.V)

                normalx = cos.(grid_u.LS[iLS].α)
                normaly = sin.(grid_v.LS[iLS].α)
            
                grid_u.V .*= normalx
                grid_v.V .*= normaly

                print("\n num.advection_LS_mode == 3 iLS", iLS)


                update_radius_from_contact_line(num,grid_p, grid_p.LS[iLS].u, BC_u)


                printstyled(color=:green, @sprintf "\n grid_p p u v max : %.2e %.2e %.2e\n" maximum(abs.(grid_p.V[grid_p.LS[iLS].MIXED])) maximum(abs.(grid_u.V[grid_p.LS[iLS].MIXED])) maximum(abs.(grid_v.V[grid_v.LS[iLS].MIXED])))

                advection_u_and_v(grid_p, grid_u, grid_v, iLS, θ_out, num, BC_int, BC_u, rhs_LS, periodic_x, periodic_y)
                rhs_LS .= 0.0
                grid_p.LS[iLS].A.nzval .= 0.0
                grid_p.LS[iLS].B.nzval .= 0.0
                IIOE!(grid_p, grid_u, grid_v, grid_p.LS[iLS].A, grid_p.LS[iLS].B, θ_out, num.timestep_n, periodic_x, periodic_y)
                BC_LS_interior!(num, grid_p, grid_u, grid_v, iLS, grid_p.LS[iLS].A, grid_p.LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
                BC_LS!(grid_p, grid_p.LS[iLS].u, grid_p.LS[iLS].A, grid_p.LS[iLS].B, rhs_LS, BC_u)
                utmp .= reshape(gmres(grid_p.LS[iLS].A, grid_p.LS[iLS].B * vec(grid_p.LS[iLS].u) .+ rhs_LS), grid_p)

                rhs_LS .= 0.0
                S2IIOE!(grid_p, grid_u, grid_v, grid_p.LS[iLS].A, grid_p.LS[iLS].B, utmp, grid_p.LS[iLS].u, θ_out, num.timestep_n, periodic_x, periodic_y)
                BC_LS_interior!(num, grid_p, grid_u, grid_v, iLS, grid_p.LS[iLS].A, grid_p.LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
                BC_LS!(grid_p, grid_p.LS[iLS].u, grid_p.LS[iLS].A, grid_p.LS[iLS].B, rhs_LS, BC_u)
                grid_p.LS[iLS].u .= reshape(gmres(grid_p.LS[iLS].A, grid_p.LS[iLS].B * vec(grid_p.LS[iLS].u) .+ rhs_LS), grid_p)


                update_radius_from_contact_line(num,grid_p, grid_p.LS[iLS].u, BC_u)

            # elseif num.advection_LS_mode == 4

            #     previous_radius = num.current_radius

            #     sign_mass_transfer_rate = 
            #     # Minus sign because normal points toward bubble and varnH2 for gaz, not liquid phase 
            #     varnH2 = sign_mass_transfer_rate * sum(mass_transfer_rate) * num.diffusion_coeff[num.index_phase_change] 

            #     #TODO mode_2d==0 flux corresponds to cylinder of length 1
            #     #2D cylinder reference length
            #     if mode_2d == 1
            #         varnH2 .*= num.ref_thickness_2d
            #     end

                
            #     nH2 = nH2 + varnH2 * num.timestep_n

            #     # printstyled(color=:green, @sprintf "\n it %.5i Mole: %.2e dn %.2e new nH2 %.2e \n" num.current_iter nH2 varnH2*num.timestep_n new_nH2)

            #     # if varnH2 < 0.0 
            #     #     # print(@sprintf "error nH2 %.2e dnH2 %.2e new nH2 %.2e\n" nH2-varnH2*num.timestep_n varnH2*num.timestep_n nH2 )
            #     #     @error ("error nH2")
            #     #     crashed = true
            #     #     new_nH2 = nH2
            #     #     print("wrong nH2 ")
            #     #     # println(@sprintf "\n CRASHED after %d iterations \n" num.current_iter)
            #     #     return
            #     # end

            #     # if occursin("Khalighi_no_update",electrolysis_phase_change_case)

            #     # else
            #     #     nH2 = new_nH2
            #     # end
                
            #     #TODO using num.temperature0
            #     if mode_2d == 0
            #         num.current_radius = cbrt(3.0 * nH2 * num.Ru * num.temperature0/( 4.0 * pi * p_g) )
            #     elseif mode_2d == 1
            #         num.current_radius = sqrt(nH2 * num.Ru * num.temperature0/( pi * p_g * num.ref_thickness_2d) )
            #     elseif mode_2d == 2
            #         num.current_radius = sqrt(nH2/(num.concentration0[num.index_phase_change] * pi))
            #     elseif mode_2d == 3
            #         num.current_radius = sqrt(2*nH2/(num.concentration0[num.index_phase_change] * pi))
            #     end

            #     printstyled(color=:green, @sprintf "\n radius num.CFL: %.2e \n" (num.current_radius-previous_radius)/(num.L0/grid_p.nx))

            #     if (num.current_radius-previous_radius)/(num.L0/grid_p.nx) > 0.5
            #         printstyled(color=:red, @sprintf "\n radius num.CFL: %.2e \n" (num.current_radius-previous_radius)/(num.L0/grid_p.nx))
            #         @error ("num.CFL radius")
            #         crashed = true
            #         return
            #     end

                
            #     printstyled(color=:cyan, @sprintf "\n div(0,grad): %.5i %.2e %.2e %.2e \n" grid_p.nx num.timestep_n num.L0/grid_p.nx (num.current_radius-previous_radius)/(num.L0/grid_p.nx)) 
            #     # sum(mass_transfer_rate))
                
            #     printstyled(color=:green, @sprintf "\n num.n(H2): %.2e added %.2e old R %.2e new R %.2e \n" nH2 varnH2*num.timestep_n previous_radius num.current_radius)
            #     printstyled(color=:green, @sprintf "\n p0: %.2e p_liq %.2e p_lapl %.2e \n" num.pres0 p_liq p_g)

            #     if mode_2d == 3
            #         grid_p.LS[1].u .= sqrt.((grid_p.x.- num.xcoord).^ 2 + (grid_p.y .- num.ycoord) .^ 2) - (num.current_radius) * ones(grid_p.ny, grid_p.nx)                  
            #     else
            #         grid_p.LS[1].u .= sqrt.((grid_p.x .- num.xcoord .- num.current_radius .+ num.R ).^ 2 + (grid_p.y .- num.ycoord) .^ 2) - (num.current_radius) * ones(grid_p.ny, grid_p.nx)
            #     end


            elseif num.advection_LS_mode == 5
                print("\n num.advection_LS_mode == 5 iLS", iLS)


                update_radius_from_contact_line(num,grid_p, grid_p.LS[iLS].u, BC_u)

                grid_p.V .=0.25*grid_p.dx[1,1]/num.timestep_n  

                printstyled(color=:green, @sprintf "\n grid_p p u v max : %.2e %.2e %.2e\n" maximum(abs.(grid_p.V[grid_p.LS[iLS].MIXED])) maximum(abs.(grid_u.V[grid_p.LS[iLS].MIXED])) maximum(abs.(grid_v.V[grid_v.LS[iLS].MIXED])))

                IIOE_normal!(grid_p, grid_p.LS[iLS].A, grid_p.LS[iLS].B, grid_p.LS[iLS].u, grid_p.V, CFL_sc, periodic_x, periodic_y)
                grid_p.LS[iLS].u .= reshape(gmres(grid_p.LS[iLS].A, grid_p.LS[iLS].B * vec(grid_p.LS[iLS].u)), grid_p)



                update_radius_from_contact_line(num,grid_p, grid_p.LS[iLS].u, BC_u)

            elseif (num.advection_LS_mode == 6) || (num.advection_LS_mode == 7)

                rhs_LS .= 0.0
                # grid_p.LS[iLS].A.nzval .= 0.0
                # grid_p.LS[iLS].B.nzval .= 0.0

                print("\n num.advection_LS_mode == 2 iLS", iLS)


                update_radius_from_contact_line(num,grid_p, grid_p.LS[iLS].u, BC_u)


                printstyled(color=:green, @sprintf "\n grid_p p u v max : %.2e %.2e %.2e\n" maximum(abs.(grid_p.V[grid_p.LS[iLS].MIXED])) maximum(abs.(grid_u.V[grid_p.LS[iLS].MIXED])) maximum(abs.(grid_v.V[grid_v.LS[iLS].MIXED])))

                IIOE_normal!(grid_p, grid_p.LS[iLS].A, grid_p.LS[iLS].B, grid_p.LS[iLS].u, grid_p.V, CFL_sc, periodic_x, periodic_y)

                # IIOE_normal_indices!(grid_p, grid_p.LS[iLS].A, grid_p.LS[iLS].B, grid_p.LS[iLS].u, grid_p.V, CFL_sc, periodic_x, periodic_y,grid_p.ind.all_indices)
                # grid_p.LS[iLS].u .= reshape(gmres(grid_p.LS[iLS].A, grid_p.LS[iLS].B * vec(grid_p.LS[iLS].u)), grid_p)


                # rhs_LS .= 0.0
                # grid_p.LS[iLS].A.nzval .= 0.0
                # grid_p.LS[iLS].B.nzval .= 0.0
                # IIOE!(grid_p, grid_u, grid_v, grid_p.LS[iLS].A, grid_p.LS[iLS].B, θ_out, num.timestep_n, periodic_x, periodic_y)
                # BC_LS_interior!(num, grid_p, grid_u, grid_v, iLS, grid_p.LS[iLS].A, grid_p.LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
                # BC_LS!(grid_p, grid_p.LS[iLS].u, grid_p.LS[iLS].A, grid_p.LS[iLS].B, rhs_LS, BC_u)
                # utmp .= reshape(gmres(grid_p.LS[iLS].A, grid_p.LS[iLS].B * vec(grid_p.LS[iLS].u) .+ rhs_LS), grid_p)

                # rhs_LS .= 0.0
                # S2IIOE!(grid_p, grid_u, grid_v, grid_p.LS[iLS].A, grid_p.LS[iLS].B, utmp, grid_p.LS[iLS].u, θ_out, num.timestep_n, periodic_x, periodic_y)
                # BC_LS_interior!(num, grid_p, grid_u, grid_v, iLS, grid_p.LS[iLS].A, grid_p.LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
                # BC_LS!(grid_p, grid_p.LS[iLS].u, grid_p.LS[iLS].A, grid_p.LS[iLS].B, rhs_LS, BC_u)
                # grid_p.LS[iLS].u .= reshape(gmres(grid_p.LS[iLS].A, grid_p.LS[iLS].B * vec(grid_p.LS[iLS].u) .+ rhs_LS), grid_p)

                

                if num.advection_LS_mode == 6

                    # IIOE_normal!(grid_p, grid_p.LS[iLS].A, grid_p.LS[iLS].B, grid_p.LS[iLS].u, grid_p.V, CFL_sc, periodic_x, periodic_y)
                    BC_LS!(grid_p, grid_p.LS[iLS].u, grid_p.LS[iLS].A, grid_p.LS[iLS].B, rhs_LS, BC_u)
                    BC_LS_interior!(num, grid_p, grid_u, grid_v, iLS, grid_p.LS[iLS].A, grid_p.LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
                
                end
                
                grid_p.LS[iLS].u .= reshape(gmres(grid_p.LS[iLS].A, grid_p.LS[iLS].B * vec(grid_p.LS[iLS].u) .+ rhs_LS), grid_p)

                # # Impose contact angle if a wall is present
                # rhs_LS .= 0.0
                # grid_p.LS[iLS].A.nzval .= 0.0
                # grid_p.LS[iLS].B.nzval .= 0.0
                # for II in grid_p.ind.all_indices
                #     pII = lexicographic(II, grid_p.ny)
                #     grid_p.LS[iLS].A[pII,pII] = 1.0
                #     grid_p.LS[iLS].B[pII,pII] = 1.0
                # end
                # BC_LS_interior!(num, grid_p, iLS, grid_p.LS[iLS].A, grid_p.LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
                # grid_p.LS[iLS].u .= reshape(gmres(grid_p.LS[iLS].A, grid_p.LS[iLS].B * vec(grid_p.LS[iLS].u) .+ rhs_LS), grid_p)





                update_radius_from_contact_line(num,grid_p, grid_p.LS[iLS].u, BC_u)


            elseif num.advection_LS_mode == 8
                print("\n num.advection_LS_mode == 8 iLS", iLS)
                print("\n deprecated", iLS)

                # nghost = 1

                # Aghost, Bghost = allocate_ghost_matrices(grid_p.nx,grid_p.ny,nghost)

                # print("\n periodic ",periodic_x," y ",periodic_y)


                # update_radius_from_contact_line(num,grid_p, grid_p.LS[iLS].u, BC_u)


                # printstyled(color=:green, @sprintf "\n grid_p p u v max : %.2e %.2e %.2e\n" maximum(abs.(grid_p.V[grid_p.LS[iLS].MIXED])) maximum(abs.(grid_u.V[grid_p.LS[iLS].MIXED])) maximum(abs.(grid_v.V[grid_v.LS[iLS].MIXED])))

                # # IIOE_normal!(grid_p, grid_p.LS[iLS].A, grid_p.LS[iLS].B, grid_p.LS[iLS].u, grid_p.V, CFL_sc, periodic_x, periodic_y)
                # # IIOE_normal_indices!(grid_p, grid_p.LS[iLS].A, grid_p.LS[iLS].B, grid_p.LS[iLS].u, grid_p.V, CFL_sc, periodic_x, periodic_y,grid_p.ind.all_indices)

                # # print("\n sizes LS A ",size(grid_p.LS[iLS].A), " B ", size(grid_p.LS[iLS].B)," LS ",size(grid_p.LS[iLS].u)," V ",size(grid_p.V),"\n")

                # #recopy value in ghost cell
                # LSghost = init_ghost_neumann(grid_p.LS[iLS].u,grid_p.nx,grid_p.ny,nghost)
                
                # Vghost = init_ghost_neumann(grid_p.V,grid_p.nx,grid_p.ny,nghost)

                # # print("\n LSghost \n")
                # # print("\n ",LSghost[0,:])
                # # print("\n ",LSghost[1,:])
                # # print("\n ",LSghost[:,0])
                # # print("\n ",LSghost[:,1])

                # IIOE_normal_indices!(grid_p, Aghost, Bghost, grid_p.LS[iLS].u, LSghost, 
                # grid_p.V, CFL_sc, periodic_x, periodic_y,grid_p.ind.all_indices)

                # # IIOE_normal_indices!(grid_p, grid_p.LS[iLS].A, grid_p.LS[iLS].B, grid_p.LS[iLS].u, LSghost, 
                # # grid_p.V, CFL_sc, periodic_x, periodic_y,grid_p.ind.all_indices)

                # # IIOE_normal_indices!(grid_p, Aghost, Bghost, LSghost, Vghost, CFL_sc, periodic_x, periodic_y,grid_p.ind.all_indices)

                # # print("\n sizes LS A ",size(grid_p.LS[iLS].A), " B ", size(grid_p.LS[iLS].B)," LS ",size(grid_p.LS[iLS].u)," V ",size(grid_p.V),"\n")
                # # print("\n sizes LS A ",size(OffsetArrays.no_offset_view(Aghost)), " B ", size(OffsetArrays.no_offset_view(Bghost))," LS ",size(LSghost)," V ",size(grid_p.V),"\n")

                # OffsetArrays.no_offset_view(LSghost) .= reshape(gmres(OffsetArrays.no_offset_view(Aghost), OffsetArrays.no_offset_view(Bghost) * vec(OffsetArrays.no_offset_view(LSghost))), (grid_p.ny+2,grid_p.nx+2))

                # grid_p.LS[iLS].u .= LSghost[1:grid_p.ny,1:grid_p.nx]

                # # grid_p.LS[iLS].u .= reshape(gmres(OffsetArrays.no_offset_view(Aghost), OffsetArrays.no_offset_view(Bghost) * vec(grid_p.LS[iLS].u)), grid_p)



                # update_radius_from_contact_line(num,grid_p, grid_p.LS[iLS].u, BC_u)

            elseif ((num.advection_LS_mode == 9) || (num.advection_LS_mode == 10))
                print("\n num.advection_LS_mode == 9 or 10 iLS", iLS)

                if num.advection_LS_mode == 10
                    grid_p.V .=0.25*grid_p.dx[2,2]/num.timestep_n  
                    
                    print("\n dummy adv", 0.25*grid_p.dx[2,2]/num.timestep_n , " ",grid_p.dx[2,2]," ",num.timestep_n)

                end

                nghost = 1

                Aghost, Bghost = allocate_ghost_matrices_2(grid_p.nx,grid_p.ny,nghost)


                update_radius_from_contact_line(num,grid_p, grid_p.LS[iLS].u, BC_u)


                printstyled(color=:green, @sprintf "\n grid_p p u v max : %.2e %.2e %.2e\n" maximum(abs.(grid_p.V[grid_p.LS[iLS].MIXED])) maximum(abs.(grid_u.V[grid_p.LS[iLS].MIXED])) maximum(abs.(grid_v.V[grid_v.LS[iLS].MIXED])))

                # IIOE_normal!(grid_p, grid_p.LS[iLS].A, grid_p.LS[iLS].B, grid_p.LS[iLS].u, grid_p.V, CFL_sc, periodic_x, periodic_y)
                # IIOE_normal_indices!(grid_p, grid_p.LS[iLS].A, grid_p.LS[iLS].B, grid_p.LS[iLS].u, grid_p.V, CFL_sc, periodic_x, periodic_y,grid_p.ind.all_indices)

                # print("\n sizes LS A ",size(grid_p.LS[iLS].A), " B ", size(grid_p.LS[iLS].B)," LS ",size(grid_p.LS[iLS].u)," V ",size(grid_p.V),"\n")


                LSghost = init_ghost_neumann_2(grid_p.LS[iLS].u,grid_p.nx,grid_p.ny,nghost)

                # print("\n LSghost \n")
                # print("\n ",LSghost[0,:])
                # print("\n ",LSghost[1,:])
                # print("\n ",LSghost[:,0])
                # print("\n ",LSghost[:,1])
                                        
                Vghost = init_ghost_neumann_2(grid_p.V,grid_p.nx,grid_p.ny,nghost)

                IIOE_normal_indices_2!(grid_p, Aghost, Bghost, grid_p.LS[iLS].u, LSghost, 
                Vghost, CFL_sc, periodic_x, periodic_y,nghost)

                # IIOE_normal_indices!(grid_p, grid_p.LS[iLS].A, grid_p.LS[iLS].B, grid_p.LS[iLS].u, LSghost, 
                # grid_p.V, CFL_sc, periodic_x, periodic_y,grid_p.ind.all_indices)

                # IIOE_normal_indices!(grid_p, Aghost, Bghost, LSghost, Vghost, CFL_sc, periodic_x, periodic_y,grid_p.ind.all_indices)

                # print("\n sizes LS A ",size(grid_p.LS[iLS].A), " B ", size(grid_p.LS[iLS].B)," LS ",size(grid_p.LS[iLS].u)," V ",size(grid_p.V),"\n")
                # print("\n sizes LS A ",size(OffsetArrays.no_offset_view(Aghost)), " B ", size(OffsetArrays.no_offset_view(Bghost))," LS ",size(LSghost)," V ",size(grid_p.V),"\n")

                LSghost .= reshape(gmres(Aghost, Bghost * vec(LSghost)), (grid_p.ny+2*nghost,grid_p.nx+2*nghost))


                # grid_p.LS[iLS].u .= LSghost[1:grid_p.ny,1:grid_p.nx]

                #Store result of LS advection without ghost cells
                for j=1:grid_p.ny
                    for i=1:grid_p.nx
                        grid_p.LS[iLS].u[j,i] = LSghost[j+1,i+1]
                    end
                end

                # grid_p.LS[iLS].u .= reshape(gmres(OffsetArrays.no_offset_view(Aghost), OffsetArrays.no_offset_view(Bghost) * vec(grid_p.LS[iLS].u)), grid_p)



                update_radius_from_contact_line(num,grid_p, grid_p.LS[iLS].u, BC_u)

            elseif num.advection_LS_mode == 11 || num.advection_LS_mode == 12 

                printstyled(color=:magenta, @sprintf "\n update_free_surface_velocity")
                #TODO
                update_free_surface_velocity(num, grid_u, grid_v, 1, phL.uD, phL.vD, periodic_x, periodic_y)

                # display(grid_v.V)

                # "write_"

                i_u_ext, l_u_ext, b_u_ext, r_u_ext, t_u_ext = indices_extension(grid_u, grid_u.LS[1], grid_u.ind.inside, periodic_x, periodic_y)
                i_v_ext, l_v_ext, b_v_ext, r_v_ext, t_v_ext = indices_extension(grid_v, grid_v.LS[1], grid_v.ind.inside, periodic_x, periodic_y)

                field_extension!(grid_u, grid_u.LS[1].u, grid_u.V, i_u_ext, l_u_ext, b_u_ext, r_u_ext, t_u_ext, num.NB, periodic_x, periodic_y)
                field_extension!(grid_v, grid_v.LS[1].u, grid_v.V, i_v_ext, l_v_ext, b_v_ext, r_v_ext, t_v_ext, num.NB, periodic_x, periodic_y)
                
                # printstyled(color=:magenta, @sprintf "\n extended")

                # display(grid_v.V)



                if num.advection_LS_mode == 11 
                    print("\n dummy grid_v.V ")

                    grid_v.V .=0.25*grid_p.dx[2,2]/num.timestep_n  
                end

                print("\n grid_v.V ")
                printstyled(color=:green, @sprintf "\n grid_p p u v max : %.2e %.2e %.2e\n" maximum(abs.(grid_p.V[grid_p.LS[iLS].MIXED])) maximum(abs.(grid_u.V[grid_u.LS[iLS].MIXED])) maximum(abs.(grid_v.V[grid_v.LS[iLS].MIXED])))


                rhs_LS .= 0.0
                grid_p.LS[iLS].A.nzval .= 0.0
                grid_p.LS[iLS].B.nzval .= 0.0
                IIOE!(grid_p, grid_u, grid_v, grid_p.LS[iLS].A, grid_p.LS[iLS].B, θ_out, num.timestep_n, periodic_x, periodic_y)
                # BC_LS_interior!(num, grid_p, grid_u, grid_v, iLS, grid_p.LS[iLS].A, grid_p.LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
                BC_LS!(grid_p, grid_p.LS[iLS].u, grid_p.LS[iLS].A, grid_p.LS[iLS].B, rhs_LS, BC_u)
                utmp .= reshape(gmres(grid_p.LS[iLS].A, grid_p.LS[iLS].B * vec(grid_p.LS[iLS].u) .+ rhs_LS), grid_p)

                rhs_LS .= 0.0
                S2IIOE!(grid_p, grid_u, grid_v, grid_p.LS[iLS].A, grid_p.LS[iLS].B, utmp, grid_p.LS[iLS].u, θ_out, num.timestep_n, periodic_x, periodic_y)
                # BC_LS_interior!(num, grid_p, grid_u, grid_v, iLS, grid_p.LS[iLS].A, grid_p.LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
                BC_LS!(grid_p, grid_p.LS[iLS].u, grid_p.LS[iLS].A, grid_p.LS[iLS].B, rhs_LS, BC_u)
                grid_p.LS[iLS].u .= reshape(gmres(grid_p.LS[iLS].A, grid_p.LS[iLS].B * vec(grid_p.LS[iLS].u) .+ rhs_LS), grid_p)

            #region bulk +phase-change velocity    
                elseif num.advection_LS_mode == 13 || num.advection_LS_mode == 14 || num.advection_LS_mode == 15 || num.advection_LS_mode == 16 || num.advection_LS_mode == 17

                if num.time > num.nucleation_time #TODO more precisely no mass transfer but velocity 
                    
                    num.phase_change_currently_activated = 1 #TODO before after
                    
                    #region bulk velocity
                    
                    if num.advection_LS_mode == 16
                        grid_u.V .= u_extended
                        grid_v.V .= v_extended
                    else
                        grid_u.V .= phL.u #reshape(vec1(phL.uD,grid_u), grid_u)
                        grid_v.V .= phL.v #reshape(vec1(phL.vD,grid_v), grid_v)
                    end
                    
                    PDI_status = @ccall "libpdi".PDI_multi_expose("write_advection_velocity_bulk"::Cstring,
                    "advection_velocity_bulk_u"::Cstring, grid_u.V::Ptr{Cdouble}, PDI_OUT::Cint,
                    "advection_velocity_bulk_v"::Cstring, grid_v.V::Ptr{Cdouble}, PDI_OUT::Cint,                                       
                    C_NULL::Ptr{Cvoid})::Cint
                    
                    #pbm localisation interp vs normal

                    interpolate_staggered_u_v_to_scalar_grid_one_fluid_or_one_phase!(num,grid_p,grid_u,grid_v,phL.u,phL.v,tmp_vec_p,tmp_vec_p0)

                    PDI_status = @ccall "libpdi".PDI_multi_expose("write_interpolated_advection_velocity_bulk"::Cstring,                               
                    "advection_velocity_bulk_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
                    "advection_velocity_bulk_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,                                  
                    C_NULL::Ptr{Cvoid})::Cint

                    #endregion bulk velocity


                    #region add phase-change contribution (localised)
                    tmp_vec_u0 .= 0.0 
                    tmp_vec_v0 .= 0.0
                    # project normal contribution to u and v : tmp_vec_u0 and tmp_vec_v0
                    interpolate_scalar_Dirac_to_u_v!(grid_p, grid_u, grid_v, grid_p.V, tmp_vec_u0, tmp_vec_v0)

                    PDI_status = @ccall "libpdi".PDI_multi_expose("write_mass_transfer_rate_uv"::Cstring,
                    "mass_transfer_rate_u"::Cstring, tmp_vec_u0::Ptr{Cdouble}, PDI_OUT::Cint,
                    "mass_transfer_rate_v"::Cstring, tmp_vec_v0::Ptr{Cdouble}, PDI_OUT::Cint,                                       
                    C_NULL::Ptr{Cvoid})::Cint

                    tmp_vec_u .= 0.0
                    tmp_vec_v .= 0.0
                    tmp_vec_u1 .= 0.0
                    tmp_vec_v1 .= 0.0
                    compute_unit_normal(num,grid_p, grid_u, grid_v, 
                    op.opC_uL, op.opC_vL,levelset_1D,
                    volume_fraction,
                    tmp_vec_p,tmp_vec_p0,
                    tmp_vec_u,tmp_vec_v, #normal_and_dirac_u, normal_and_dirac_v,
                    tmp_vec_u1,tmp_vec_v1,#normal_u, normal_v, 
                    )

                    PDI_status = @ccall "libpdi".PDI_multi_expose("write_normal_uv"::Cstring,
                    "normal_u"::Cstring, tmp_vec_u1::Ptr{Cdouble}, PDI_OUT::Cint,
                    "normal_v"::Cstring, tmp_vec_v1::Ptr{Cdouble}, PDI_OUT::Cint,                                       
                    C_NULL::Ptr{Cvoid})::Cint

                    tmp_vec_u0 .*= tmp_vec_u1
                    tmp_vec_v0 .*= tmp_vec_v1

                    PDI_status = @ccall "libpdi".PDI_multi_expose("write_advection_velocity_phase_change"::Cstring,
                    "advection_velocity_phase_change_u"::Cstring, tmp_vec_u0::Ptr{Cdouble}, PDI_OUT::Cint,
                    "advection_velocity_phase_change_v"::Cstring, tmp_vec_v0::Ptr{Cdouble}, PDI_OUT::Cint,                                       
                    C_NULL::Ptr{Cvoid})::Cint

                    interpolate_staggered_u_v_to_scalar_grid_one_fluid_or_one_phase!(num,grid_p,grid_u,grid_v,
                    tmp_vec_u0,tmp_vec_v0,tmp_vec_p,tmp_vec_p0)

                    PDI_status = @ccall "libpdi".PDI_multi_expose("write_interpolated_advection_velocity_phase_change"::Cstring,                               
                    "advection_velocity_phase_change_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
                    "advection_velocity_phase_change_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,                                  
                    C_NULL::Ptr{Cvoid})::Cint

                    grid_u.V .+= tmp_vec_u0
                    grid_v.V .+= tmp_vec_v0
                    
                    #endregion add phase-change contribution (localised)

                    PDI_status = @ccall "libpdi".PDI_multi_expose("write_advection_velocity_before_extension"::Cstring,
                    "advection_velocity_before_extension_u"::Cstring, grid_u.V::Ptr{Cdouble}, PDI_OUT::Cint,
                    "advection_velocity_before_extension_v"::Cstring, grid_v.V::Ptr{Cdouble}, PDI_OUT::Cint,                                       
                    C_NULL::Ptr{Cvoid})::Cint

                    interpolate_staggered_u_v_to_scalar_grid_one_fluid_or_one_phase!(num,grid_p,grid_u,grid_v,
                    grid_u.V,grid_v.V,tmp_vec_p,tmp_vec_p0)

                    PDI_status = @ccall "libpdi".PDI_multi_expose("write_interpolated_advection_velocity_before_extension"::Cstring,                               
                    "advection_velocity_before_extension_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
                    "advection_velocity_before_extension_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,                                  
                    C_NULL::Ptr{Cvoid})::Cint

        
                    if num.extend_field == 0
                        i_u_ext, l_u_ext, b_u_ext, r_u_ext, t_u_ext = indices_extension(grid_u, grid_u.LS[iLS], grid_u.ind.inside, periodic_x, periodic_y)
                        i_v_ext, l_v_ext, b_v_ext, r_v_ext, t_v_ext = indices_extension(grid_v, grid_v.LS[iLS], grid_v.ind.inside, periodic_x, periodic_y)

                        field_extension!(grid_u, grid_u.LS[iLS].u, grid_u.V, i_u_ext, l_u_ext, b_u_ext, r_u_ext, t_u_ext, num.NB, periodic_x, periodic_y)
                        field_extension!(grid_v, grid_v.LS[iLS].u, grid_v.V, i_v_ext, l_v_ext, b_v_ext, r_v_ext, t_v_ext, num.NB, periodic_x, periodic_y)
                    end


                    PDI_status = @ccall "libpdi".PDI_multi_expose("check_advection"::Cstring,
                    "levelset_p"::Cstring, grid_p.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                    "advection_velocity_p"::Cstring, grid_p.V::Ptr{Cdouble}, PDI_OUT::Cint,
                    "advection_velocity_u"::Cstring, grid_u.V::Ptr{Cdouble}, PDI_OUT::Cint,
                    "advection_velocity_v"::Cstring, grid_v.V::Ptr{Cdouble}, PDI_OUT::Cint,                        
                    C_NULL::Ptr{Cvoid})::Cint

                    PDI_status = @ccall "libpdi".PDI_multi_expose("write_advection"::Cstring,
                    "levelset_p"::Cstring, grid_p.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                    "advection_velocity_p"::Cstring, grid_p.V::Ptr{Cdouble}, PDI_OUT::Cint,
                    "advection_velocity_u"::Cstring, grid_u.V::Ptr{Cdouble}, PDI_OUT::Cint,
                    "advection_velocity_v"::Cstring, grid_v.V::Ptr{Cdouble}, PDI_OUT::Cint,                        
                    C_NULL::Ptr{Cvoid})::Cint

                    PDI_status = @ccall "libpdi".PDI_multi_expose("write_advection_velocity"::Cstring,
                    "advection_velocity_u"::Cstring, grid_u.V::Ptr{Cdouble}, PDI_OUT::Cint,
                    "advection_velocity_v"::Cstring, grid_v.V::Ptr{Cdouble}, PDI_OUT::Cint,                                       
                    C_NULL::Ptr{Cvoid})::Cint

        

                    if num.advection_LS_mode == 13
                        rhs_LS .= 0.0
                        grid_p.LS[iLS].A.nzval .= 0.0
                        grid_p.LS[iLS].B.nzval .= 0.0
                        IIOE!(grid_p, grid_u, grid_v, grid_p.LS[iLS].A, grid_p.LS[iLS].B, θ_out, num.timestep_n, periodic_x, periodic_y)
                        BC_LS_interior!(num, grid_p, grid_u, grid_v, iLS, grid_p.LS[iLS].A, grid_p.LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
                        BC_LS!(grid_p, grid_p.LS[iLS].u, grid_p.LS[iLS].A, grid_p.LS[iLS].B, rhs_LS, BC_u)
                        utmp .= reshape(gmres(grid_p.LS[iLS].A, grid_p.LS[iLS].B * vec(grid_p.LS[iLS].u) .+ rhs_LS), grid_p)

                        rhs_LS .= 0.0
                        S2IIOE!(grid_p, grid_u, grid_v, grid_p.LS[iLS].A, grid_p.LS[iLS].B, utmp, grid_p.LS[iLS].u, θ_out, num.timestep_n, periodic_x, periodic_y)
                        BC_LS_interior!(num, grid_p, grid_u, grid_v, iLS, grid_p.LS[iLS].A, grid_p.LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
                        BC_LS!(grid_p, grid_p.LS[iLS].u, grid_p.LS[iLS].A, grid_p.LS[iLS].B, rhs_LS, BC_u)
                        grid_p.LS[iLS].u .= reshape(gmres(grid_p.LS[iLS].A, grid_p.LS[iLS].B * vec(grid_p.LS[iLS].u) .+ rhs_LS), grid_p)
                    elseif num.advection_LS_mode == 14
                        #region ghost cell adv in normal direction
                        if num.extend_field == 0
                            i_ext, l_ext, b_ext, r_ext, t_ext = indices_extension(grid_p, grid_p.LS[iLS], grid_p.ind.inside, periodic_x, periodic_y)
                            field_extension!(grid_p, grid_p.LS[iLS].u, grid_p.V, i_ext, l_ext, b_ext, r_ext, t_ext, num.NB, periodic_x, periodic_y)
                        end

                        nghost = 1
                        Aghost, Bghost = allocate_ghost_matrices_2(grid_p.nx,grid_p.ny,nghost)
                        # update_radius_from_contact_line(num,grid_p, grid_p.LS[iLS].u, BC_u)
                        printstyled(color=:green, @sprintf "\n grid_p p u v max : %.2e %.2e %.2e\n" maximum(abs.(grid_p.V[grid_p.LS[iLS].MIXED])) maximum(abs.(grid_u.V[grid_p.LS[iLS].MIXED])) maximum(abs.(grid_v.V[grid_v.LS[iLS].MIXED])))
                        LSghost = init_ghost_neumann_2(grid_p.LS[iLS].u,grid_p.nx,grid_p.ny,nghost)
                        Vghost = init_ghost_neumann_2(grid_p.V,grid_p.nx,grid_p.ny,nghost)
                        IIOE_normal_indices_2!(grid_p, Aghost, Bghost, grid_p.LS[iLS].u, LSghost, 
                        Vghost, CFL_sc, periodic_x, periodic_y,nghost)
                        LSghost .= reshape(gmres(Aghost, Bghost * vec(LSghost)), (grid_p.ny+2*nghost,grid_p.nx+2*nghost))

                        #Store result of LS advection without ghost cells
                        for j=1:grid_p.ny
                            for i=1:grid_p.nx
                                grid_p.LS[iLS].u[j,i] = LSghost[j+1,i+1]
                            end
                        end

                        #endregion ghost cell adv in normal direction

                    elseif num.advection_LS_mode == 15

                        #region ghost cell adv in normal direction
                        # Project velocities to the normal and use advection scheme for advection just
                        # in the normal direction
                        
                        # tmpVx = zeros(grid_p)
                        # tmpVy = zeros(grid_p)
                        # # grid_p.V .= 0.0 #keep phase change contribution
                        # @inbounds @threads for II in grid_p.LS[iLS].MIXED
                        #     cap1 = grid_u.LS[iLS].geoL.cap[II,5]
                        #     cap3 = grid_u.LS[iLS].geoL.cap[δx⁺(II),5]
                        #     tmpVx[II] = (grid_u.V[II] * cap1 + grid_u.V[δx⁺(II)] * cap3) / (cap1 + cap3 + eps(0.01)) 

                        #     cap2 = grid_v.LS[iLS].geoL.cap[II,5]
                        #     cap4 = grid_v.LS[iLS].geoL.cap[δy⁺(II),5]
                        #     tmpVy[II] = (grid_v.V[II] * cap2 + grid_v.V[δy⁺(II)] * cap4) / (cap2 + cap4 + eps(0.01))

                        #     tmpV = sqrt(tmpVx[II]^2 + tmpVy[II]^2)
                        #     β = atan(tmpVy[II], tmpVx[II])
                        #     if grid_p.LS[iLS].α[II] > 0.0 && β < 0.0
                        #         β += 2π
                        #     end
                        #     if grid_p.LS[iLS].α[II] < 0.0 && β > 0.0
                        #         β -= 2π
                        #     end

                        #     grid_p.V[II] += tmpV * cos(β - grid_p.LS[iLS].α[II])
                        # end
                        
                        # Project velocities to the normal and use advection scheme for advection just
                        # in the normal direction
                        @inbounds @threads for II in grid_p.LS[iLS].MIXED
                            grid_p.V[II] += compute_normal_component_of_velocity(num,grid_u, grid_v, grid_u.V, grid_v.V, grid_p, iLS, II)
                        end


                        # gridp.V .+= interp

                        if num.extend_field == 0
                            i_ext, l_ext, b_ext, r_ext, t_ext = indices_extension(grid_p, grid_p.LS[iLS], grid_p.ind.inside, periodic_x, periodic_y)
                            field_extension!(grid_p, grid_p.LS[iLS].u, grid_p.V, i_ext, l_ext, b_ext, r_ext, t_ext, num.NB, periodic_x, periodic_y)
                        end

                        nghost = 1
                        Aghost, Bghost = allocate_ghost_matrices_2(grid_p.nx,grid_p.ny,nghost)
                        # update_radius_from_contact_line(num,grid_p, grid_p.LS[iLS].u, BC_u)
                        printstyled(color=:green, @sprintf "\n grid_p p u v max : %.2e %.2e %.2e\n" maximum(abs.(grid_p.V[grid_p.LS[iLS].MIXED])) maximum(abs.(grid_u.V[grid_p.LS[iLS].MIXED])) maximum(abs.(grid_v.V[grid_v.LS[iLS].MIXED])))
                        LSghost = init_ghost_neumann_2(grid_p.LS[iLS].u,grid_p.nx,grid_p.ny,nghost)
                        Vghost = init_ghost_neumann_2(grid_p.V,grid_p.nx,grid_p.ny,nghost)
                        IIOE_normal_indices_2!(grid_p, Aghost, Bghost, grid_p.LS[iLS].u, LSghost, 
                        Vghost, CFL_sc, periodic_x, periodic_y,nghost)
                        LSghost .= reshape(gmres(Aghost, Bghost * vec(LSghost)), (grid_p.ny+2*nghost,grid_p.nx+2*nghost))

                        #Store result of LS advection without ghost cells
                        for j=1:grid_p.ny
                            for i=1:grid_p.nx
                                grid_p.LS[iLS].u[j,i] = LSghost[j+1,i+1]
                            end
                        end

                        #endregion ghost cell adv in normal direction


                    elseif num.advection_LS_mode == 16 || num.advection_LS_mode == 17

                        #region ghost cell adv in normal direction
                        # Project velocities to the normal and use advection scheme for advection just
                        # in the normal direction

                        if num.advection_LS_mode == 16
                            # grid_p.V not reset to zero to keep phase change contribution
                            @inbounds @threads for II in grid_p.LS[iLS].MIXED
                                grid_p.V[II] += compute_normal_component_of_velocity(num,grid_u, grid_v, u_extended, v_extended, grid_p, iLS, II)
                            end
                        end

                        # gridp.V .+= interp

                        if num.extend_field == 0
                            i_ext, l_ext, b_ext, r_ext, t_ext = indices_extension(grid_p, grid_p.LS[iLS], grid_p.ind.inside, periodic_x, periodic_y)
                            field_extension!(grid_p, grid_p.LS[iLS].u, grid_p.V, i_ext, l_ext, b_ext, r_ext, t_ext, num.NB, periodic_x, periodic_y)
                        end

                        

                        PDI_status = @ccall "libpdi".PDI_multi_expose("write_normal_velocity_intfc_LS_ext"::Cstring,
                        "normal_velocity_intfc_LS_ext"::Cstring, grid_p.V::Ptr{Cdouble}, PDI_OUT::Cint,
                        C_NULL::Ptr{Cvoid})::Cint

                        nghost = 1
                        Aghost, Bghost = allocate_ghost_matrices_2(grid_p.nx,grid_p.ny,nghost)
                        # update_radius_from_contact_line(num,grid_p, grid_p.LS[iLS].u, BC_u)
                        printstyled(color=:green, @sprintf "\n grid_p p u v max : %.2e %.2e %.2e\n" maximum(abs.(grid_p.V[grid_p.LS[iLS].MIXED])) maximum(abs.(grid_u.V[grid_p.LS[iLS].MIXED])) maximum(abs.(grid_v.V[grid_v.LS[iLS].MIXED])))
                        LSghost = init_ghost_neumann_2(grid_p.LS[iLS].u,grid_p.nx,grid_p.ny,nghost)
                        Vghost = init_ghost_neumann_2(grid_p.V,grid_p.nx,grid_p.ny,nghost)
                        IIOE_normal_indices_2!(grid_p, Aghost, Bghost, grid_p.LS[iLS].u, LSghost, 
                        Vghost, CFL_sc, periodic_x, periodic_y,nghost)
                        LSghost .= reshape(gmres(Aghost, Bghost * vec(LSghost)), (grid_p.ny+2*nghost,grid_p.nx+2*nghost))

                        #Store result of LS advection without ghost cells
                        for j=1:grid_p.ny
                            for i=1:grid_p.nx
                                grid_p.LS[iLS].u[j,i] = LSghost[j+1,i+1]
                            end
                        end

                        #endregion ghost cell adv in normal direction



                    end
                    
                    update_radius_from_contact_line(num,grid_p, grid_p.LS[iLS].u, BC_u)

                    radius_pdi = [0.0]

                    PDI_status = @ccall "libpdi".PDI_multi_expose("compute_radius"::Cstring,
                    "levelset_p"::Cstring, grid_p.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                    "mesh_p_x"::Cstring, grid_p.x::Ptr{Cdouble}, PDI_OUT::Cint,
                    "mesh_p_y"::Cstring, grid_p.y::Ptr{Cdouble}, PDI_OUT::Cint,
                    # "radius"::Cstring, test_radius::Ref{Cdouble}, PDI_INOUT::Cint,  
                    "radius_vec"::Cstring, radius_pdi::Ptr{Cdouble}, PDI_INOUT::Cint,                             
                    C_NULL::Ptr{Cvoid})::Cint

                    

                    num.current_radius = radius_pdi[1]

        

                else
                    printstyled(color=:red, @sprintf "\n no levelset advection before nucleation \n" )

                end

                printstyled(color=:red, @sprintf "\n after advection 13 radius: %.2e \n" num.current_radius)


            #endregion bulk +phase-change velocity    


            end #num.advection_LS_mode == 
    end # if is_stefan(bc)  or ... #fs or phase change levelse
    end

    # printstyled(color=:red, @sprintf "\n after advection_LS_mode radius: %.2e \n" num.current_radius)

end



