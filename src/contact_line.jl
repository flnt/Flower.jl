## This file is part of the Flower project.
## It contains the routines and procedures to compute the Young stresses.
## It gives values to gp.Young as function of the BCs and the levelset (gp.u) on the domain.
## @author: R. A. S. Frantz (2023)

function prtnn(r) # for debug 
    println("Values: ", [x for x in vec(r) if x != 0 && !isnan(x)])
end

function update_levelset_contact_angle(periodic_x, periodic_y, BC, gp, tracer, LSRHS)

    function compute_theta_ext(gp, tracer)
        if is_navier(BC.bottom.t)
            for II in gp.ind.b_bottom[1]
                if !isnothing(gp.α[II]) && gp.α[II] > 0
                    θapp = gp.α[δy⁺(II)] # apparent is located on grid above the boundary
                    # compute first order (forward) derivative of the tracer in the normal direction
                    ∂n = (tracer[δy⁺(II)] - tracer[II]) / gp.dy[II]
                    # compute extrapolated angle using apparent angle 
                    gp.θext[II] = θapp + 1.5 * gp.dy[II] * ((gp.κ[II] * sqrt(1.0 + ∂n)) / sin(θapp))
                end
            end
        end 
        if is_navier(BC.top.t)
            for II in gp.ind.b_top[1]
                if !isnothing(gp.α[II]) && gp.α[II] > 0
                    θapp = gp.α[δy⁻(II)]
                    ∂n = (tracer[II] - tracer[δy⁻(II)]) / gp.dy[II]
                    gp.θext[II] = θapp + 1.5 * gp.dy[II] * ((gp.κ[II] * sqrt(1.0 + ∂n)) / sin(θapp))
                end
            end
        end
    end

    compute_theta_ext(gp, tracer)

    if !periodic_y

        if is_navier(BC.top.t)
            for (II) in gp.ind.b_top[1]
                i = lexicographic(II, gp.ny)
                #LSRHS[i] = cos(gp.θext[II])
            end
        end # is_navier(BC.top.t)
        if is_navier(BC.bottom.t)
            for (II) in gp.ind.b_bottom[1]
                i = lexicographic(II, gp.ny)
                #LSRHS[i] = cos(gp.θext[II])
            end
        end # is_navier(BC.bottom.t)
    end # periodic_y
    return LSRHS
end

function update_levelset_matrices(periodic_x, periodic_y, BC, grid, LSA, LSB)
    # Imposing Newman boundary condition to the levelset function
    # LSA contains the product of the mass matrix and the Laplacian for the advection scheme
    # LSAΦ^{n+1} = LSBΦ^{n} + LSC * u_target 
    # The boundary conditions are Φ^{n+1}_{i,1} = Φ^{n+1}_{i,2} and Φ^{n+1}_{i,end} = Φ^{n+1}_{i,end-1} 

    if !periodic_y

        if is_navier(BC.bottom.t)
            for (II, JJ) in zip(grid.ind.b_bottom[1], grid.ind.b_bottom[2]) # first and second rows
                i = lexicographic(II, grid.ny)
                j = lexicographic(JJ, grid.ny)
                LSA[i, i] = 1.0 # set 1st row of LSA to 1
                LSB[i, i] = 0.0 # set 1st row of LSB to 0
                LSA[i, j] = -1.0 # set 2nd row of LSA to -1
            end
        end

        if is_navier(BC.top.t)
            for (II, JJ) in zip(grid.ind.b_top[1], grid.ind.b_top[2]) # first and second rows
                i = lexicographic(II, grid.ny)
                j = lexicographic(JJ, grid.ny)
                LSA[i, i] = 1.0 # set 1st row of LSA to 1
                LSB[i, i] = 0.0 # set 1st row of LSB to 0
                LSA[i, j] = -1.0 # set 2nd row of LSA to -1
            end
        end

    end
    return LSA, LSB
end

function update_Young_stress(mixed, bc, gu, gv, num)

    #unfinished function
    #Young_stress_internal(mixed, gu, gv, num)

    if is_navier(bc.top.t)
        gu.Young[end, :] .= Young_stress_boundary(gu, num, gu.ind.b_top[1])
    end
    if is_navier(bc.bottom.t)
        gu.Young[1, :] .= Young_stress_boundary(gu, num, gu.ind.b_bottom[1])
    end

    # if is_navier(bc.left.t)
    #     gu.Young[:, 1] .= Young_stress_boundary(gu, num, gu.ind.b_left[1])
    # end
    # if is_navier(bc.right.t)
    #     gu.Young[:, end] .= Young_stress_boundary(gu, num, gu.ind.b_right[1])
    # end

end

function Young_stress_boundary(gp, num, indices)

    function locate_index_boundary(gp, ind)
        index_CL = falses(length(ind)) # Initialize boolean array of same length as `ind`    
        
        @inbounds for i in eachindex(ind) # loop over indices in `ind`

            # Check if iso value is nonzero and if index is on a boundary
            iso_value = gp.iso[ind][i]
            
            if iso_value != 0.0 && ind in (gp.ind.b_bottom[1], gp.ind.b_top[1], gp.ind.b_left[1], gp.ind.b_right[1])
                # Set corresponding value in `index_CL` based on iso value
                index_CL[i] = (ind == gp.ind.b_bottom[1] && iso_value in (1.0, 2.0, 6.0, 9.0, 13.0, 14.0)) ||
                              (ind == gp.ind.b_top[1] && iso_value in (4.0, 6.0, 7.0, 8.0, 9.0, 11.0)) ||
                              (ind == gp.ind.b_left[1] && iso_value in (1.0, 3.0, 7.0, 8.0, 12.0, 14.0)) ||
                              (ind == gp.ind.b_right[1] && iso_value in (2.0, 3.0, 4.0, 11.0, 12.0, 13.0))
            end

        end
        return index_CL # Return boolean array indicating whether each index satisfies condition
    end

    function get_contact_line_points_boundary(gp, indices, index_CL)
        # initialize x_CL_vec as a zero vector of size indices
        x_CL_vec = zeros(size(indices))

        for i in findall(index_CL) # loop over all non-zero elements of index_CL
        
            if indices == gp.ind.b_bottom[1]
        
                x = gp.cut_points[1, i] # retrieve the cut points along the bottom boundary
                if x[1].y == -0.5 # check which of the two points is on the boundary
                    x_CL_vec[i] = x[1].x
                elseif x[2].y == -0.5
                    x_CL_vec[i] = x[2].x
                end
                x_CL_vec[i] *= gp.dx[1, i]
                x_CL_vec[i] += gp.x[1, i] # add the x-coordinate of the left boundary
        
            elseif indices == gp.ind.b_top[1]
        
                x = gp.cut_points[end, i] # retrieve the cut points along the top boundary
                if x[1].y == 0.5 # check which of the two points is on the boundary
                    x_CL_vec[i] = x[1].x
                elseif x[2].y == 0.5
                    x_CL_vec[i] = x[2].x
                end
                x_CL_vec[i] *= gp.dx[end, i]
                x_CL_vec[i] += gp.x[end, i] # add the x-coordinate of the right boundary
        
            elseif indices == gp.ind.b_left[1]
        
                x = gp.cut_points[i, 1] # retrieve the cut points along the left boundary
                if x[1].x == -0.5 # check which of the two points is on the boundary
                    x_CL_vec[i] = x[1].y
                elseif x[2].x == -0.5
                    x_CL_vec[i] = x[2].y
                end
                x_CL_vec[i] *= gp.dy[i, 1]
                x_CL_vec[i] += gp.y[i, 1] # add the y-coordinate of the bottom boundary
        
            elseif indices == gp.ind.b_right[1]
        
                x = gp.cut_points[i, end] # retrieve the cut points along the right boundary
                if x[1].x == 0.5 # check which of the two points is on the boundary
                    x_CL_vec[i] = x[1].y
                elseif x[2].x == 0.5
                    x_CL_vec[i] = x[2].y
                end
                x_CL_vec[i] *= gp.dy[i, end]
                x_CL_vec[i] += gp.y[i, end] # add the y-coordinate of the top boundary
        
            end
        end
        return x_CL_vec # return the vector of cut point coordinates
    end

    # find the indices of the domain boundary where the contact line is located
    index_CL = locate_index_boundary(gp, indices)

    # get the x coordinates of the contact line at the boundary
    x_CL_vec = get_contact_line_points_boundary(gp, indices, index_CL)
    
    # initialize stress vector and create a temporary array for relative x coordinates
    stress_vector = zeros(size(indices))
    rel_x = similar(stress_vector)
    
    # create a temporary array for the bell function
    bell_function = zeros(size(indices))
    
    for i in findall(index_CL) # loop over all non-zeros

        if indices == gp.ind.b_bottom[1]
        
            rel_x .= gp.x[1, :] .- x_CL_vec[i]
            contact_angle = gp.α[1, i]
            target_angle = gp.θe[1, i]
        
        elseif indices == gp.ind.b_top[1]
        
            rel_x .= gp.x[end, i] .- x_CL_vec[i]
            contact_angle = gp.α[end, i]
            target_angle = gp.θe[end, i]
        
        elseif indices == gp.ind.b_left[1]
        
            rel_x .= gp.y[i, 1] .- x_CL_vec[i]
            contact_angle = gp.α[i, 1]
            target_angle = gp.θe[i, 1]
        
        elseif indices == gp.ind.b_right[1]
        
            rel_x .= gp.y[i, end] .- x_CL_vec[i]
            contact_angle = gp.α[i, end]
            target_angle = gp.θe[i, end]
        
        end

        # compute the cosine value of the contact angle and apply stress only if it is less than the threshold
        cos_value = (cos(contact_angle) - cos(target_angle)) / num.Ca

        #cos_boolean = isless(contact_angle, target_angle)
        #cos_value *= cos_boolean ? 1.0 : 0.0
        # compute the bell function for the current x coordinates and add to the stress vector
        bell_function .= (1.0 .- tanh.(rel_x[:] ./ num.εCA) .^ 2) ./ num.εCA
        
        stress_vector .+= bell_function .* cos_value

    end
    return stress_vector
end


# function Young_stress_internal(mixed, grid_u, grid_v, num)

#     function locate_CL_coordinates(gp, mixed)
#         mixed_CL = CartesianIndex{2}[]
#         for i in mixed
#             iso_value = gp.iso[i]
#             if iso_value != 0 && iso_value != 15
#                 push!(mixed_CL, i)
#             end
#         end
#         return mixed_CL
#     end

#     mixed_CL = locate_CL_coordinates(grid_u, mixed)
#     #@show mixed_CL

#     # how to compute the stress with curvature ???

# end