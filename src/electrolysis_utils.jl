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


"""init bulk interfacial and border values of field"""
function init_fields_2!(TD,T,H,BC,grid,dir_val_intfc)

    vec1(TD,grid) .= vec(T)
    vec2(TD,grid) .= dir_val_intfc

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

