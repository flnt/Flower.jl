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
    
    TODO Neumann...

"""
function init_fields_multiple_levelsets!(num,TD,T,H,BC,grid,dir_val_intfc,str)

    vec1(TD,grid) .= vec(T)

    if str == "uL"

        #TODO multiple LS: init grid.V
        iLS = 1

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
        vec2(TD,grid) .= dir_val_intfc
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

    l1_rel_error = l1_rel_error / l1_rel_error_den
    l2_rel_error = sqrt(l2_rel_error / l2_rel_error_den )
    linfty_rel_error = linfty_rel_error / linfty_rel_error_den

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
function compute_interface_average(scalar_1D_vec, grid, iLS)
    min_scal = 0.0
    max_scal = 0.0
    average = 0.0
    count = 0 
    
    # volume = 0.0
    index = iLS+1


    @inbounds for II in grid.LS[iLS].MIXED
        # if grid.LS[iLS].iso[II] < 14.5 #not solid (15) (i.e. liquid or mixed)
            
            pII = lexicographic(II, grid.ny)
            index_1D = grid.ny*grid.nx*(iLS) + pII

            min_scal = scalar_1D_vec[index_1D]
            max_scal = scalar_1D_vec[index_1D]
            break

        # end 
    end

    @inbounds for II in grid.LS[iLS].MIXED
        # if grid.LS[iLS].iso[II] < 14.5 #not solid (15) (i.e. liquid or mixed)
            

            # cf veci @view a[g.ny*g.nx*(p-1)+1:g.ny*g.nx*p]
            
            # II = CartesianIndex(jplot, iplot) #(id_y, id_x)
            pII = lexicographic(II, grid.ny)

            # index_1D = grid.ny*grid.nx*(index-1) + pII
            index_1D = grid.ny*grid.nx*(iLS) + pII

            min_scal = min(min_scal,scalar_1D_vec[index_1D])
            max_scal = max(max_scal,scalar_1D_vec[index_1D])

            average += scalar_1D_vec[index_1D]
            count += 1

        # end 
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