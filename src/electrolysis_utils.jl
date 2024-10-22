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


""""height for intialisation"""
function get_height!(grid,ind,dx,dy,geo,H)
    @inbounds @threads for II in vcat(ind.b_left[1], ind.b_bottom[1], ind.b_right[1], ind.b_top[1])
        H[II] = distance(grid.LS[1].mid_point[II], geo.centroid[II], dx[II], dy[II])
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
        dist = distance(grid.LS[iLS].mid_point[II], grid.LS[iLS].geoL.centroid[II], grid.dx[II], grid.dy[II]) #grid.LS[iLS].geoL.centroid ? or end?

        veci(TD,grid,iLS+1)[pII] = dir_val_intfc + dist * BC.LS[iLS].val

        printstyled(color=:green, @sprintf "\n Neumann %.2e %.2e %.2e %.2e" dir_val_intfc dist BC.LS[iLS].val veci(TD,grid,iLS+1)[pII])

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

    if str =="scalL"
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
Rf(θ, V)
Returns the 
"""
@inline Rf(θ, V) = sqrt(V / (θ - sin(θ) * cos(θ)))


"""
RR0(θ)
Returns the 
"""
@inline RR0(θ) = sqrt(π / (2 * (θ - sin(θ) * cos(θ))))

