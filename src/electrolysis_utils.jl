"""
Compute average of values when geo.cap[II,5] > eps 
"""
function average!(T::Matrix, grid, geo,num)
    @unpack ind = grid
    @unpack eps, = num
    average= 0.0
    numcells=0
    @inbounds @threads for II in ind.all_indices
        if geo.cap[II,5] >= eps
            average +=T[II]
            numcells +=1 
        end
    end
    return average/numcells
end


"""
center(r, θ)

Returns the 
"""
@inline center(r, θ) = r * cos(π - θ)


"""
get distance between one LS and centroid of cell (defined by all LS) for intialisation
    
"""
function get_height!(LS,ind,dx,dy,geo,H)

    # printstyled(color=:green, @sprintf "\n get height \n")
    # print("\n ind left",ind.b_left[1], " bottom "  ,ind.b_bottom[1]," right ", ind.b_right[1]," top ", ind.b_top[1] )
    # print("\n H ",size(H)," \n")
    # print("\n LS.mid_point ",size(LS.mid_point)," \n")
    # print("\n dx ",size(dx)," \n")
    # print("\n dy ",size(dy)," \n")
    # print("\n geo.centroid ",size(geo.centroid)," \n")

    @inbounds @threads for II in vcat(ind.b_left[1], ind.b_bottom[1], ind.b_right[1], ind.b_top[1])
        H[II] = distance(LS.mid_point[II], geo.centroid[II], dx[II], dy[II])
    end   
end


# """init bulk interfacial and border values of field"""
# function init_fields_2!(TD,T,H,BC,grid,dir_val_intfc)

#     vec1(TD,grid) .= vec(T)
#     vec2(TD,grid) .= dir_val_intfc

#     if is_neumann(BC.left)
#         # printstyled(color=:green, @sprintf "\n init_fields_2! sizes : %.5i : %.5i : %.5i: %.5i \n" size(vecb_L(TD,grid)) size(T[:,1]) size(H[:,1]) size(BC.left.val))
        
#         vecb_L(TD,grid) .= T[:,1] .+ H[:,1] .* BC.left.val
#     else
#         vecb_L(TD,grid) .= BC.left.val #.* ones(grid.ny)
#     end

#     if is_neumann(BC.bottom)
#         vecb_B(TD,grid) .= T[1,:] .+ H[1,:] .* BC.bottom.val
#     else
#         vecb_B(TD,grid) .= BC.bottom.val #.* ones(grid.nx)
#     end

#     if is_neumann(BC.right)
#         vecb_R(TD,grid) .= T[:,end] .+ H[:,end] .* BC.right.val 
#     else
#         vecb_R(TD,grid) .= BC.right.val #.* ones(grid.ny)
#     end

#     if is_neumann(BC.top)
#         vecb_T(TD,grid) .= T[end,:] .+ H[end,:] .* BC.top.val
#     else
#         vecb_T(TD,grid) .= BC.top.val #.* ones(grid.nx)
#     end

# end


# function _DBG_type_and_size(_s_var_name::AbstractString, _value::Any, _is_DBG::Bool)
#     if _is_DBG
#         if ~isempty(size(_value))
#             println(@sprintf("%-60s, type: ", _s_var_name), @sprintf("%20s", typeof(_value)), ",\t value: ", size(_value))
#         elseif isreal(_value)
#             println(@sprintf("%-60s, type: ", _s_var_name), @sprintf("%20s", typeof(_value)), ",\t value: ", @sprintf("%.3f", _value))
#         else
#             println(@sprintf("%-60s, type: ", _s_var_name), @sprintf("%20s", typeof(_value)), ",\t size: ",  _value)
#         end
#         flush(stdout)
#     end
# end

# macro _DBG_investigate_variable(var, enabled)
#     name = string(var)
#     return esc(:(_DBG_type_and_size($name, $var, $enabled)))
# end


"""
init_Neumann_iLS

dist for geoL, iLS (not LS[end])
"""
function init_Neumann_iLS(num,TD,BC,grid,dir_val_intfc,iLS)

    for II in grid.LS[iLS].MIXED
                            

        # a0[II] .= butler_volmer_no_concentration_potential_Neumann.(num,
        # reshape(veci(ph.phi_eleD, grid,iLS+1),grid),
        # reshape(veci(ph.trans_scalD[:,2],grid,iLS+1),grid),
        # num.temperature0)

        pII = lexicographic(II, grid.ny)

        # H[II] = distance(grid.LS[1].mid_point[II], geo.centroid[II], dx[II], dy[II])
        #end for centroid
        dist = distance(grid.LS[iLS].mid_point[II], grid.LS[end].geoL.centroid[II], grid.dx[II], grid.dy[II]) #grid.LS[iLS].geoL.centroid ? or end?

        veci(TD,grid,iLS+1)[pII] = dir_val_intfc + dist * BC.LS[iLS].val

        if dist == 0
            printstyled(color=:red, @sprintf "\n Neumann %.2e %.2e %.2e %.2e" dir_val_intfc dist BC.LS[iLS].val veci(TD,grid,iLS+1)[pII])
        else
            printstyled(color=:green, @sprintf "\n Neumann %.2e %.2e %.2e %.2e" dir_val_intfc dist BC.LS[iLS].val veci(TD,grid,iLS+1)[pII])
        end
        # # if grid.LS[iLS].geoL.cap[II,5] < num.ϵ
    
        # if grid.LS[end].geoL.cap[II,5] > num.ϵ #TODO clearer eps

        #     a0[II] = butler_volmer_no_concentration_potential_Neumann.(num,
        #     veci(ph.phi_eleD, grid,iLS+1)[pII],
        #     veci(ph.trans_scalD[:,2],grid,iLS+1)[pII],
        #     num.temperature0)

        #     if veci(ph.trans_scalD[:,2],grid,iLS+1)[pII] < num.ϵ
        #         a0[II] = butler_volmer_no_concentration_potential_Neumann.(num,
        #         reshape(veci(ph.phi_eleD, grid,iLS+1),grid),
        #         ph.trans_scal[II,2],
        #         num.temperature0)
        #     end

        #     # print("\n II",II,"BC ", BC.LS[iLS].val)
        #     # printstyled(color=:red, @sprintf "\n Butler %.2e %.2e \n" a0[II] reshape(veci(ph.trans_scalD[:,2],grid,iLS+1),grid)[II])


        #     # a0[II] = butler_volmer_no_concentration_potential_Neumann.(num,
        #     # reshape(veci(ph.phi_eleD, grid,iLS+1),grid)[II],
        #     # reshape(veci(ph.trans_scalD[:,2],grid,iLS+1),grid)[II],
        #     # num.temperature0)


        #     if grid.LS[iLS].geoL.cap[II,5] < num.ϵ #volume
        #         #use bulk conductivity of mixed cell
        #         # butler_volmer_no_concentration_potential_Neumann!.(num,
        #         # reshape(veci(ph.phi_eleD, grid,iLS+1),grid),
        #         # ph.trans_scal[II,2],
        #         # num.temperature0,
        #         # a0[II]) #TODO if temperature solved temperature[II]

        #         a0[II] = butler_volmer_no_concentration_potential_Neumann.(num,
        #         reshape(veci(ph.phi_eleD, grid,iLS+1),grid),
        #         ph.trans_scal[II,2],
        #         num.temperature0) #TODO if temperature solved temperature[II]

        #         # a0[II]  = butler_volmer_no_concentration.(num.alpha_a,num.alpha_c,num.Faraday,num.i0,veci(ph.phi_eleD, grid,iLS+1),
        #         # num.phi_ele1,num.Ru,num.temperature0)./ph.elec_cond[II]
        #     end

        #     if grid.LS[iLS].geoL.cap[II,1] < num.ϵ #here length so volume approx num.ϵ^2
        #         a0[II] = 1.0
        #     end

        # end #grid.LS[end].geoL.cap[II,5] > num.ϵ: liquid cell
        # # TODO
        # #Remove Nan when dividing by conductivity which may be null
        # # kill_dead_bc_left_wall!(vecb(elec_condD,grid), grid, iLS,1.0)
        #     #Remove Nan when dividing by conductivity which may be null

    end   

end


"""
    init bulk interfacial and border values of field
    
    For Neumann BC, one may extrapolate spatially a value

    

    For a radial flow, the velocity components u and v are calculated using the boundary conditions
    and the angle of the normal to the interface computed from the Levelset.


"""
function init_fields_multiple_levelsets!(num,TD,T,H,BC,grid,dir_val_intfc,str)

    if BC.init_mode == "None" || BC.init_mode == "False" #no init
        return
    end

    vec1(TD,grid) .= vec(T)

    if str == "uL"

        #TODO multiple LS: init grid.V
        iLS = 1
        # For radial flow 
        if BC.LS[iLS].val != 0.0

            printstyled(color=:green, @sprintf "\n Initialising velocity from BC\n")

            for II in grid.ind.inside
              
                grid.V[II] = BC.LS[iLS].val * cos(grid.LS[iLS].α[II]+π) # u
            end

            print("\n BC velocity int ", BC)

        end

    end

    if str == "vL"

        #TODO multiple LS: init grid.V
        iLS = 1

        
        if BC.LS[iLS].val != 0.0

            printstyled(color=:green, @sprintf "\n Initialising velocity from BC\n")

            for II in grid.ind.inside
            
                grid.V[II] = BC.LS[iLS].val * sin(grid.LS[iLS].α[II]+π) # v
            end
            
            print("\n BC velocity int ", BC)

        end

    end


    if str =="scalL" && num.nLS>1
        for iLS in 1:num.nLS
            try
                print(BC.LS[iLS])
                if is_dirichlet(BC.LS[iLS])
                    veci(TD,grid,iLS+1) .= BC.LS[iLS].val
                    printstyled(color=:green, @sprintf "\n Dirichlet iLS %.2i %.2e\n" iLS BC.LS[iLS].val)

                elseif is_neumann(BC.LS[iLS]) #TODO all init  in main_current_folder.jl or here?
                    init_Neumann_iLS(num,TD,BC,grid,dir_val_intfc,iLS)
                else 
                    print("\n BC TODO",BC.LS[iLS])
                end

                # veci(TD,grid,iLS+1) .= BC.LS[iLS].val
                printstyled(color=:green, @sprintf "\n init iLS %.2i %.2e\n" iLS BC.LS[iLS].val)
            catch e
                print(e)
                @error("BC not defined for var iLS")
                # @_DBG_investigate_variable(TD, true)
            end

            # if is_dirichlet(BC.LS[iLS])
            #     veci(TD,grid,iLS+1) .= BC.LS[iLS].val
            # else #TODO all init  in main_current_folder.jl or here?
            #     veci(TD,grid,iLS+1) .= BC.LS[iLS].val
            # end
        end

    else
        if num.nLS>1
            vec2(TD,grid) .= dir_val_intfc
        end
    end

    if is_neumann(BC.left)
        # printstyled(color=:green, @sprintf "\n init_fields_2! sizes : %.5i : %.5i : %.5i: %.5i \n" size(vecb_L(TD,grid)) size(T[:,1]) size(H[:,1]) size(BC.left.val))
        
        vecb_L(TD,grid) .= T[:,1] .+ H[:,1] .* BC.left.val
    else
        vecb_L(TD,grid) .= BC.left.val #.* ones(grid.ny)
    end

    if is_neumann(BC.bottom)
        vecb_B(TD,grid) .= T[1,:] .+ H[1,:] .* BC.bottom.val
    else
        vecb_B(TD,grid) .= BC.bottom.val #.* ones(grid.nx)
    end

    if is_neumann(BC.right)
        vecb_R(TD,grid) .= T[:,end] .+ H[:,end] .* BC.right.val 
    else
        vecb_R(TD,grid) .= BC.right.val #.* ones(grid.ny)
    end

    if is_neumann(BC.top)
        vecb_T(TD,grid) .= T[end,:] .+ H[end,:] .* BC.top.val
    else
        vecb_T(TD,grid) .= BC.top.val #.* ones(grid.nx)
    end

end


"""
Computes average value at interface for scalar
"""
function mean_intfc_non_null(scalD,iscal,grid,iLS)

    index = iLS+1
    num=0
    nonzero = 0.0

    # cf veci @view a[g.ny*g.nx*(p-1)+1:g.ny*g.nx*p]

    for i in grid.ny*grid.nx*(index-1)+1:grid.ny*grid.nx*index
        if abs(scalD[i,iscal]) .> 0.0
            nonzero += scalD[i,iscal]
            num += 1
        end
    end

    if num == 0
        print("\n no intfc in mean_intfc_non_null")
        return 0
    else
        nonzero /= num
        return nonzero
    end
    
end


"""
Computes average value at interface for scalar
"""
function mean_intfc_non_null_v2(scalD,grid,iLS)

    index = iLS+1
    num=0
    nonzero = 0.0

    # cf veci @view a[g.ny*g.nx*(p-1)+1:g.ny*g.nx*p]

    for i in grid.ny*grid.nx*(index-1)+1:grid.ny*grid.nx*index
        if abs(scalD[i]) .> 0.0
            nonzero += scalD[i]
            num += 1
        end
    end

    if num == 0
        print("\n no intfc in mean_intfc_non_null")
        return 0
    else
        nonzero /= num
        return nonzero
    end
    
end

"""
Computes average value at interface for scalar
    index: 1 gives bulk
    2 gives 1st interface
"""
function mean_intfc_non_null_v3(scalD,grid,index)

    num=0
    nonzero = 0.0

    # cf veci @view a[g.ny*g.nx*(p-1)+1:g.ny*g.nx*p]

    for i in grid.ny*grid.nx*(index-1)+1:grid.ny*grid.nx*(index)
        if abs(scalD[i]) .> 0.0
            nonzero += scalD[i]
            num += 1
        end
    end

    if num == 0
        print("\n no intfc in mean_intfc_non_null")
        return 0
    else
        nonzero /= num
        return nonzero
    end
    
end


"""
Rf(θ, V)
Returns the 
"""
@inline Rf(θ, V) = sqrt(V / (θ - sin(θ) * cos(θ)))


"""
RR0(θ)
Returns the 
"""
@inline RR0(θ) = sqrt(π / (2 * (θ - sin(θ) * cos(θ))))


"""
From kinetic_energy

!!! TODO

"""
function scal_magnitude(phL, phS, gp, gu, gv)
    #TODO eps, den eps
    LS_u =gu.LS[1]
    LS_v = gv.LS[1]
    phL.p .= (
        (phL.u[:,2:end].^2.0 .* LS_u.geoL.dcap[:,2:end,6] .+ 
        phL.u[:,1:end-1].^2.0 .* LS_u.geoL.dcap[:,1:end-1,6]) ./ 
        (LS_u.geoL.dcap[:,1:end-1,6] .+ LS_u.geoL.dcap[:,2:end,6] .+ 1e-8 )
    )
    phL.p .+= (
        (phL.v[2:end,:].^2.0 .* LS_v.geoL.dcap[2:end,:,7] .+ 
        phL.v[1:end-1,:].^2.0 .* LS_v.geoL.dcap[1:end-1,:,7]) ./
        (LS_v.geoL.dcap[1:end-1,:,7] .+ LS_v.geoL.dcap[2:end,:,7] .+ 1e-8 )
    )
    phL.p .+= (
        (phS.u[:,2:end].^2.0 .* LS_u.geoS.dcap[:,2:end,6] .+ 
        phS.u[:,1:end-1].^2.0 .* LS_u.geoS.dcap[:,1:end-1,6]) ./ 
        (LS_u.geoS.dcap[:,1:end-1,6] .+ LS_u.geoS.dcap[:,2:end,6] .+ 1e-8 )
    )
    phL.p .+= (
        (phS.v[2:end,:].^2.0 .* LS_v.geoS.dcap[2:end,:,7] .+ 
        phS.v[1:end-1,:].^2.0 .* LS_v.geoS.dcap[1:end-1,:,7]) ./
        (LS_v.geoS.dcap[1:end-1,:,7] .+ LS_v.geoS.dcap[2:end,:,7] .+ 1e-8 )
    )

    #####################################################################
    # phL.p .= (
    #     (phL.u[:,2:end].^2.0 .* LS_u.geoL.dcap[:,2:end,6] .+ 
    #     phL.u[:,1:end-1].^2.0 .* LS_u.geoL.dcap[:,1:end-1,6]) ./ 
    #     (LS_u.geoL.dcap[:,1:end-1,6] .+ LS_u.geoL.dcap[:,2:end,6] )
    # )
    # phL.p .+= (
    #     (phL.v[2:end,:].^2.0 .* LS_v.geoL.dcap[2:end,:,7] .+ 
    #     phL.v[1:end-1,:].^2.0 .* LS_v.geoL.dcap[1:end-1,:,7]) ./
    #     (LS_v.geoL.dcap[1:end-1,:,7] .+ LS_v.geoL.dcap[2:end,:,7] )
    # )
    # phL.p .+= (
    #     (phS.u[:,2:end].^2.0 .* LS_u.geoS.dcap[:,2:end,6] .+ 
    #     phS.u[:,1:end-1].^2.0 .* LS_u.geoS.dcap[:,1:end-1,6]) ./ 
    #     (LS_u.geoS.dcap[:,1:end-1,6] .+ LS_u.geoS.dcap[:,2:end,6] )
    # )
    # phL.p .+= (
    #     (phS.v[2:end,:].^2.0 .* LS_v.geoS.dcap[2:end,:,7] .+ 
    #     phS.v[1:end-1,:].^2.0 .* LS_v.geoS.dcap[1:end-1,:,7]) ./
    #     (LS_v.geoS.dcap[1:end-1,:,7] .+ LS_v.geoS.dcap[2:end,:,7] )
    # )

    phL.p .= sqrt.(phL.p)
end

"""
From kinetic_energy
!!! TODO

"""
function scal_magnitude_L(ph, gp, gu, gv)

    LS_u =gu.LS[1]
    LS_v = gv.LS[1]
    # LS =gp.LS[1]

    ph.p .= (
        (ph.u[:,2:end].^2.0 .* LS_u.geoL.dcap[:,2:end,6] .+ 
        ph.u[:,1:end-1].^2.0 .* LS_u.geoL.dcap[:,1:end-1,6]) ./ 
        (LS_u.geoL.dcap[:,1:end-1,6] .+ LS_u.geoL.dcap[:,2:end,6] 
        # .+ 1e-8
        )
    )
    ph.p .+= (
        (ph.v[2:end,:].^2.0 .* LS_v.geoL.dcap[2:end,:,7] .+ 
        ph.v[1:end-1,:].^2.0 .* LS_v.geoL.dcap[1:end-1,:,7]) ./
        (LS_v.geoL.dcap[1:end-1,:,7] .+ LS_v.geoL.dcap[2:end,:,7] 
        # .+ 1e-8
        )
    )
    ph.p .= sqrt.(ph.p)
    # ph.p .= sqrt.(ph.p .* LS.geoL.dcap[:,:,5])

end


"""
    computes relative errors of bulk variable for convergence study
"""
function relative_errors(T, Tanalytical, pos, cap, h)

    l1_rel_error = 0.0
    l2_rel_error = 0.0
    linfty_rel_error = 0.0
    l1_rel_error_den = 0.0 
    l2_rel_error_den = 0.0 
    linfty_rel_error_den = 0.0 
    
    volume = 0.0
    max_diff = 0.0
 

    if size(Tanalytical) != size(T) #if Tanalytical a slice (assuming it is an x slice)

        # print("\n size ",size(Tanalytical)," ",size(T))

        @inbounds for ii in pos

            volume = cap[ii]*h^2
            if volume > 0.0

                i_slice = ii[2]
                # print("\n ii ",ii," ",i_slice," ",Tanalytical[i_slice]," ",T[ii])
                abs_diff = abs(Tanalytical[i_slice] .- T[ii])
                abs_val = abs(Tanalytical[i_slice])

                l1_rel_error += volume * abs_diff
                l1_rel_error_den += volume * abs_val

                l2_rel_error += volume * abs_diff^2
                l2_rel_error_den += volume * (Tanalytical[i_slice])^2

                
                if (abs_diff > linfty_rel_error) linfty_rel_error = abs_diff end
                if (abs_val > linfty_rel_error_den) linfty_rel_error_den = abs_val end

            end
        end

    else

        @inbounds for ii in pos

            volume = cap[ii]*h^2
            if volume > 0.0
                abs_diff = abs(Tanalytical[ii] .- T[ii])
                abs_val = abs(Tanalytical[ii])

                l1_rel_error += volume * abs_diff
                l1_rel_error_den += volume * abs_val

                l2_rel_error += volume * abs_diff^2
                l2_rel_error_den += volume * (Tanalytical[ii])^2

                
                if (abs_diff > linfty_rel_error) linfty_rel_error = abs_diff end
                if (abs_val > linfty_rel_error_den) linfty_rel_error_den = abs_val end

            end
        end
    end

    


    l1_rel_error = l1_rel_error / l1_rel_error_den
    l2_rel_error = sqrt(l2_rel_error / l2_rel_error_den )
    linfty_rel_error = linfty_rel_error / linfty_rel_error_den

    print("\n l1    ",l1_rel_error ," ", l1_rel_error," ",l1_rel_error_den)
    print("\n l2    ",l2_rel_error ," ", l2_rel_error," ",l2_rel_error_den)
    print("\n linfty",linfty_rel_error ," ", linfty_rel_error," ",linfty_rel_error_den)

    # if linfty_rel_error < l1_rel_error

    return l1_rel_error, l2_rel_error, linfty_rel_error
end


"""
    computes relative errors of interfacial variable for convergence study
"""
function relative_errors_interface(T, Tanalytical, pos, cap, h)


    # l1_rel_error = 0.0
    # l2_rel_error = 0.0
    linfty_rel_error = 0.0
    # l1_rel_error_den = 0.0 
    # l2_rel_error_den = 0.0 
    linfty_rel_error_den = 0.0 
    
    volume = 0.0
    max_diff = 0.0

    @inbounds for ii in pos

        volume = cap[ii]*h^2
        if volume > 0.0
            abs_diff = abs(Tanalytical[ii] .- T[ii])
            abs_val = abs(Tanalytical[ii])

            # l1_rel_error += volume * abs_diff
            # l1_rel_error_den += volume * abs_val

            # l2_rel_error += volume * abs_diff^2
            # l2_rel_error_den += volume * (Tanalytical[ii])^2

            if (abs_diff > linfty_rel_error) linfty_rel_error = abs_diff end
            if (abs_val > linfty_rel_error_den) linfty_rel_error_den = abs_val end

        end
    end

    # l1_rel_error = l1_rel_error / l1_rel_error_den
    # l2_rel_error = sqrt(l2_rel_error / l2_rel_error_den )
    linfty_rel_error = linfty_rel_error / linfty_rel_error_den

    # return l1_rel_error, l2_rel_error, linfty_rel_error
    return linfty_rel_error
end



"""
compute the average in specified cells 
"""
function compute_interface_average(num,scalar_1D_vec, grid, iLS)
    min_scal = 0.0
    max_scal = 0.0
    average = 0.0
    count = 0 
    
    # volume = 0.0
    index = iLS+1


        # @inbounds for II in grid.LS[iLS].MIXED
        # # print("\n II update ",II, grid.LS[end].u[II], " iso end ",grid.LS[end].iso[II]," iso 1 ",grid.LS[1].iso[II])
        # if grid.LS[end].iso[II] < 14.5 #15.0 -0.5 # check if inside domain defined by other LS 
        # # if grid.LS[end].u[II]>0.0 # check if inside domain defined by other LS 
        # # if grid.LS[2].u[II]>0.0 #second wall
        #     # print("\n cells for free surface", II," x ",grid.x[II]," LS[iLS] ",grid.LS[iLS].u[II]," LS[end] ",grid.LS[end].u[II]," LS[2] ",grid.LS[2].u[II])
        #     # grid.V[II] = mass_transfer_rate[II] * factor_velocity

        #     num_mixed_cells += 1
        #     # intfc_length_cell !=0 since mixed cell

        

        #     #compute interface length
        #     χx = (grid.LS[iLS].geoL.dcap[II,3] .- grid.LS[iLS].geoL.dcap[II,1]) .^ 2
        #     χy = (grid.LS[iLS].geoL.dcap[II,4] .- grid.LS[iLS].geoL.dcap[II,2]) .^ 2
        #     intfc_length_cell = sqrt(χx + χy)

        #     if intfc_length_cell > num.epsilon_dist_mass_transfer_rate



    @inbounds for II in grid.LS[iLS].MIXED

        if grid.LS[iLS].iso[II] < 14.5 #not solid (15) (i.e. liquid or mixed)
            χx = (grid.LS[iLS].geoL.dcap[II,3] .- grid.LS[iLS].geoL.dcap[II,1]) .^ 2
            χy = (grid.LS[iLS].geoL.dcap[II,4] .- grid.LS[iLS].geoL.dcap[II,2]) .^ 2
            intfc_length_cell = sqrt(χx + χy)

            if intfc_length_cell > num.epsilon_dist_mass_transfer_rate
            
                pII = lexicographic(II, grid.ny)
                index_1D = grid.ny*grid.nx*(iLS) + pII

                min_scal = scalar_1D_vec[index_1D]
                max_scal = scalar_1D_vec[index_1D]
                break
            end

        end 
    end

    @inbounds for II in grid.LS[iLS].MIXED
        if grid.LS[iLS].iso[II] < 14.5 #not solid (15) (i.e. liquid or mixed)
            # if grid.LS[iLS].iso[II] < 14.5 #not solid (15) (i.e. liquid or mixed)
            χx = (grid.LS[iLS].geoL.dcap[II,3] .- grid.LS[iLS].geoL.dcap[II,1]) .^ 2
            χy = (grid.LS[iLS].geoL.dcap[II,4] .- grid.LS[iLS].geoL.dcap[II,2]) .^ 2
            intfc_length_cell = sqrt(χx + χy)

            if intfc_length_cell > num.epsilon_dist_mass_transfer_rate

                # cf veci @view a[g.ny*g.nx*(p-1)+1:g.ny*g.nx*p]
                
                # II = CartesianIndex(jplot, iplot) #(id_y, id_x)
                pII = lexicographic(II, grid.ny)

                # index_1D = grid.ny*grid.nx*(index-1) + pII
                index_1D = grid.ny*grid.nx*(iLS) + pII

                min_scal = min(min_scal,scalar_1D_vec[index_1D])
                max_scal = max(max_scal,scalar_1D_vec[index_1D])

                average += scalar_1D_vec[index_1D]
                count += 1
            end

        end 
    end

    if count == 0 
        print("\n no interface, no average")
    else
        average = average / count
    end

    return min_scal,max_scal,average

end


"""
compute the average in specified cells 
"""
function compute_bulk_or_interface_average(scalar_1D_vec, grid, iLS)
    min_scal = 0.0
    max_scal = 0.0
    average = 0.0
    count = 0 
    
    # volume = 0.0
    index = iLS+1


    @inbounds for II in grid.ind.all_indices
        if grid.LS[iLS].iso[II] < 14.5 #not solid (15) (i.e. liquid or mixed)
            
            pII = lexicographic(II, grid.ny)
            index_1D = grid.ny*grid.nx*(iLS) + pII

            min_scal = scalar_1D_vec[index_1D]
            max_scal = scalar_1D_vec[index_1D]
            break

        end 
    end #loop liquid + mixed

    @inbounds for II in grid.ind.all_indices
        if grid.LS[iLS].iso[II] < 14.5 #not solid (15) (i.e. liquid or mixed)
            

            # cf veci @view a[g.ny*g.nx*(p-1)+1:g.ny*g.nx*p]
            
            # II = CartesianIndex(jplot, iplot) #(id_y, id_x)
            pII = lexicographic(II, grid.ny)

            # index_1D = grid.ny*grid.nx*(index-1) + pII
            index_1D = grid.ny*grid.nx*(iLS) + pII

            min_scal = min(min_scal,scalar_1D_vec[index_1D])
            max_scal = max(max_scal,scalar_1D_vec[index_1D])

            average += scalar_1D_vec[index_1D]
            count += 1

        end 
    end #loop liquid + mixed

    if count == 0 
        print("\n no interface, no average")
    else
        average = average / count
    end

    return min_scal,max_scal,average

end

"""
find sign changes used to compute radius
"""
function find_sign_changes(slice::AbstractVector)
    # Ensure the slice has at least two elements
    if length(slice) < 2
        throw(ArgumentError("Slice must have at least two elements to detect sign changes."))
    end

    # Iterate through the slice to find the first sign change
    for i in 1:(length(slice)-1)
        if slice[i] * slice[i+1] < 0
            return (i, i+1)  # Return indices as a tuple
        end
    end

    # If no sign change is found, return nothing 
    return nothing,nothing
end

# function find_sign_changes(slice)
#     # print('len',len(slice),slice)
#     # min_dist = np.min(abs(slice))
#     # min_dist_tmp = np.max(abs(slice))
#     for i = 1:size(slice)
#         if (slice[i] * slice[i+1]) < 0
#             i1 = i
#             i2 = i+1
#             break
#         end

#     end
#     return(i1,i2)

# end

"""

"""
function compute_bubble_drop_radius(num, grid_p)
    
    volume_cell = grid_p.LS[num.iLSpdi].geoL.cap[:, :, 5]
   
    # Calculate center of mass
    center_of_mass_x, center_of_mass_y = calculate_centroid(
        grid_p.x, grid_p.y, volume_cell
    )

    # Find indices for bubble mass center
    indices_bubble_mass_center = find_slice_coord_bubble_mass_center(
        center_of_mass_x, center_of_mass_y, num, grid_p
    )

    # # Compute horizontal radius
    # radius_horizontal = compute_radius_from_levelset_slice(
    #     grid_p.LS[num.iLSpdi].u[indices_bubble_mass_center[1], :],
    #     grid_p.y[indices_bubble_mass_center[1], :]
    # )

    # Compute vertical radius
    # radius_vertical = compute_radius_from_levelset_slice(
    #     grid_p.LS[num.iLSpdi].u[:, div(grid_p.nx, 2)],
    #     grid_p.y[:, div(grid_p.nx, 2)]
    # )
    print("\n indices_bubble_mass_center ",indices_bubble_mass_center)

    slice_indices_list = [
    (1, :),  # Horizontal slice at bottom wall (for bubble at wall)
    (indices_bubble_mass_center[1], :),  # Horizontal slice at bubble mass center
    (: , div(grid_p.nx, 2)), # Vertical slice at middle of domain
    (: , indices_bubble_mass_center[2])  # Vertical slice at bubble mass center
    ]

    radii = Vector{Union{Float64, Nothing}}(undef, length(slice_indices_list))
    for (i, slice_indices) in enumerate(slice_indices_list)
        print("\n i slice ",i, " ",slice_indices)
        radii[i] = compute_radius_from_levelset_slice(
            grid_p.LS[num.iLSpdi].u, grid_p.x, grid_p.y, slice_indices)
    end

    print("\nradii ",radii)
    # print("\nradii ",skipmissing(radii))
    
    # num.current_radius = maximum(radii) 
    # num.current_radius = maximum(skipmissing(radii))

    # if isempty(radii)
    #     print("\n empty radii")
    # end

    if all(x -> x === nothing, radii)
        # error("All radii are nothing")
        print("\n All radii are nothing")
    else
        num.current_radius = maximum(x for x in radii if x !== nothing)
    end


   

    # # Handle errors and set current radius
    # if isnothing(radius_vertical) && isnothing(radius_horizontal)
    #     @error "Error: Both radius_vertical and radius_horizontal are nothing."
    # elseif isnothing(radius_vertical)
    #     @warn "Warning: radius_vertical is nothing. Using radius_horizontal."
    #     num.current_radius = radius_horizontal
    # elseif isnothing(radius_horizontal)
    #     @warn "Warning: radius_horizontal is nothing. Using radius_vertical."
    #     num.current_radius = radius_vertical
    # else
    #     num.current_radius = max(radius_horizontal, radius_vertical)
    # end

end

# """
#     compute_radii_from_slices(
#         u::AbstractMatrix,
#         x_grid::AbstractMatrix,
#         y_grid::AbstractMatrix,
#         slices::Vector{Tuple{Union{Int, Colon}, Union{Int, Colon}}};
#         center_x::Real=0.0,
#         center_y::Real=0.0
#     )

# Compute the radius for each slice in `slices` from a level set function.

# # Arguments
# - `u`: Level set function (2D matrix).
# - `x_grid`, `y_grid`: Grid coordinates (2D matrices).
# - `slices`: Vector of slice indices, e.g., `[(nx, :), (: , ny), ...]`.
# - `center_x`, `center_y`: Reference point (default: `(0.0, 0.0)`).

# # Returns
# - Vector of radii (one for each slice) or `nothing` for slices with no interface.
# """
# function compute_radii_from_slices(
#     u::AbstractMatrix,
#     x_grid::AbstractMatrix,
#     y_grid::AbstractMatrix,
#     slices::Vector{Tuple{Union{Int, Colon}, Union{Int, Colon}}};
#     center_x::Real=0.0,
#     center_y::Real=0.0
# )
#     radii = Vector{Union{Float64, Nothing}}(undef, length(slices))
#     for (i, slice) in enumerate(slices)
#         radii[i] = compute_radius_from_levelset_slice(
#             u, x_grid, y_grid, slice)
#     end
#     return radii
# end


"""

"""
function compute_radius_from_levelset_slice(
    u::AbstractMatrix,
    x_grid::AbstractMatrix,
    y_grid::AbstractMatrix,
    slice_indices)     
   

    # Determine if the slice is vertical or horizontal
    is_vertical = slice_indices[2] isa Int  
    is_horizontal = slice_indices[1] isa Int 
    print("\n slice_indices ",slice_indices)
    # Select the correct coordinate and center
    if is_vertical
        coord_slice = y_grid[:, slice_indices[2]] 
        # center = center_y
    elseif is_horizontal
        coord_slice = x_grid[slice_indices[1], :]  
        # center = center_x
    else
        error("Slice indices must be of the form `(nx, :)` or `(: , ny)`.")
    end

    slice = u[slice_indices...]
    # x_slice = x_grid[slice_indices...]
    # y_slice = y_grid[slice_indices...]

    # slice = u[slice_indices]

    dx = coord_slice[2]-coord_slice[1]

    # # print(colored('first','red'))
    # i1,i2 = find_one_minimum(slice,coord_slice,eps)
    # print('i1 i2',i1,i2)

    # # print(colored('second','red'))
    # itmp = max(i1,i2)
    # # print('itmp',itmp)
    # slice2 = slice[itmp+1:] 
    # i3,i4 = find_one_minimum(slice2,coord_slice,eps)
    # i3+= itmp+1
    # i4+= itmp+1
    # print('i3 i4',i3,i4)

    # print("\n slice ",slice)

    i1,i2 = find_sign_changes(slice)
    if (isnothing(i1) || isnothing(i2)) 
        return nothing
    end
    # print('i1 i2',i1,i2)
    itmp = max(i1,i2)
    slice2 = slice[itmp+1:end]
    i3,i4 = find_sign_changes(slice2)
    
    if (isnothing(i3) || isnothing(i4)) 
        return nothing
    end

    i3+= itmp+1
    i4+= itmp+1
    # print('i3 i4',i3,i4)


    a = (slice[i1]-slice[i2])/((coord_slice[i1]-coord_slice[i2]))
    interp1 = coord_slice[i1]-slice[i1]/a
    # print('x1',coord_slice[i1],coord_slice[i2],interp1)

    a = (slice[i3]-slice[i4])/((coord_slice[i3]-coord_slice[i4]))
    interp2 = coord_slice[i3]-slice[i3]/a
    # print('x1',coord_slice[i3],coord_slice[i4],interp2)

    radius = abs(interp2-interp1)/2

    return radius
end


"""

"""
function find_slice_coord_bubble_mass_center(center_of_mass_x,center_of_mass_y,num,grid_p)
    dx = grid_p.dx[2] - grid_p.x[1]
    dy = grid_p.dy[2] - grid_p.y[1]
   
    xmin = grid_p.x[1,1] #with regards to first scalar node (inner) (1,1), 
    #hence we can take dx, and not dx/2 or dx in the case we took num.x[1,1] the domain corner
    ymin = grid_p.y[1,1] #with regards to first scalar node (inner) (1,1)


    
    print("\n mass center ", center_of_mass_x, " ",center_of_mass_y)

    i = floor(Int,(center_of_mass_x-xmin)/dx)+1
    j = floor(Int,(center_of_mass_y-ymin)/dy)+1
    #trunc

    print("\n mass center coord i ", i, " j ", j, " xmin ",xmin, " ymin ",ymin)

    return(j,i)
end
# """
# To read BC from dict
# """
# function read_BC(dict)

#     try
#         BC = BoundariesInt(
#             left   = eval(Meta.parseall(dict.left)),
#             right   = eval(Meta.parseall(dict.right)),
#             bottom   = eval(Meta.parseall(dict.bottom)),
#             top   = eval(Meta.parseall(dict.top)),
#             LS   = eval(Meta.parseall(dict.LS)),
#           )

#         # print(BC)
#         return BC

#     catch error
#         printstyled(color=:red, @sprintf "\n Initialization error \n")
#         print(error)
#         BC = BoundariesInt()
#         # print(BC)
#         return BC
#     end

# end


# macro read_BC(dict)
#     # Generate the code for the macro
#     quote
       
#         try
#             BC = BoundariesInt(
#                 left   = eval(Meta.parseall($(dict).left)),
#                 right  = eval(Meta.parseall($(dict).right)),
#                 bottom = eval(Meta.parseall($(dict).bottom)),
#                 top    = eval(Meta.parseall($(dict).top)),
#                 LS     = eval(Meta.parseall($(dict).LS)),
#             )

#             # print(BC)
#             return BC

#         catch error
#             printstyled(color=:red, @sprintf "\n Initialization error \n")
#             print(error)
#             BC = BoundariesInt()
#             # print(BC)
#             return BC
#         end
#     end
# end
