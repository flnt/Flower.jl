## This file is part of the Flower project.
## It contains the routines and procedures to compute the Young stresses.
## It gives values to gp.Young as function of the BCs and the levelset (gp.u) on the domain.
## @author: R. A. S. Frantz (2023)

function prtnn(arr) # for debug 
    flattened = vec(arr)
    non_zero_non_nan_values = [x for x in flattened if x != 0 && !isnan(x)]
    println("Non-zero, non-NaN values: ", non_zero_non_nan_values)
end

function update_levelset_contact_angle(periodic_x,periodic_y,gp,tracer,LSRHS)
    
    function compute_theta_ext(gp,tracer)
        for (II) in gp.ind.b_bottom[1]

            # compute first order derivative of the tracer in the normal direction
            dh = (tracer[δy⁺(δy⁺(II))]-tracer[II])/(gp.dy[δy⁺(II)]+gp.dy[II])

            θapp = gp.α[δy⁺(II)] # apparent is located on grid above the boundary

            # compute extrapolated angle using apparent angle 
            gp.θext[II] = θapp + 1.5*gp.dy[II]*((gp.κ[II]*sqrt(1.0+dh))/sin(θapp))
            
            #gp.θext[II] = 90. * π / 180 # forcing angle normal to the boundary

        end
    end

    if ! periodic_y
        # for (II) in gp.ind.b_top[1]
        #     i = lexicographic(II, gp.ny)
        #     LSRHS[i] = cos(gp.θext[II])
        # end
        for (II) in gp.ind.b_bottom[1]
            i = lexicographic(II, gp.ny)
            LSRHS[i] = cos(gp.θext[II])
        end
    end 
    return LSRHS
end

function update_levelset_matrices(periodic_x,periodic_y,grid,LSA,LSB)
    # Imposing Newman boundary condition to the levelset function
    # LSA contains the product of the mass matrix and the Laplacian for the advection scheme
    # LSAΦ^{n+1} = LSBΦ^{n} + LSC * u_target 
    # The boundary conditions are Φ^{n+1}_{i,1} = Φ^{n+1}_{i,2} and Φ^{n+1}_{i,end} = Φ^{n+1}_{i,end-1} 

    if ! periodic_y

        # Handle bottom boundary
        for (II,JJ) in zip(grid.ind.b_bottom[1], grid.ind.b_bottom[2]) # first and second rows
            i = lexicographic(II, grid.ny)
            j = lexicographic(JJ, grid.ny)
            LSA[i,i] = 1. # set 1st row of LSA to 1
            LSB[i,i] = 0. # set 1st row of LSB to 0
            LSA[i,j] = -1. # set 2nd row of LSA to -1
        end

        # Handle top boundary
        # for (II,JJ) in zip(grid.ind.b_top[1], grid.ind.b_top[2]) # first and second rows
        #     i = lexicographic(II, grid.ny)
        #     j = lexicographic(JJ, grid.ny)
        #     LSA[i,i] = 1. # set 1st row of LSA to 1
        #     LSB[i,i] = 0. # set 1st row of LSB to 0
        #     LSA[i,j] = -1. # set 2nd row of LSA to -1
        # end

    end
    return LSA, LSB
end 

function update_Young_stress(mixed, grid_u, grid_v, num)

    #compute_young_stress_internal(grid_u,num,mixed)
    grid_u.Young[1,:] .= compute_young_stress(grid_u,num,grid_u.ind.b_bottom[1])
    #grid_v.Young[1,:] .= compute_young_stress(grid_v,num,grid_v.ind.b_bottom[1])
    
end 

function compute_young_stress(gp, num, indices)

    function locate_index_domain_boundary(gp, ind)
        index_CL = falses(length(ind)) # Initialize boolean array of same length as `ind`    
        # Loop over indices in `ind`
        @inbounds for i in eachindex(ind)
            # Check if iso value is nonzero and if index is on a boundary
            iso_value = gp.iso[ind][i]
            if iso_value != 0.0 && ind in (gp.ind.b_bottom[1], gp.ind.b_top[1], gp.ind.b_left[1], gp.ind.b_right[1])
                # Set corresponding value in `index_CL` based on iso value
                index_CL[i] = (ind == gp.ind.b_bottom[1] && iso_value in (1., 2., 6., 9., 13., 14.)) ||
                              (ind == gp.ind.b_top[1] && iso_value in (4., 6., 7., 8., 9., 11.)) ||
                              (ind == gp.ind.b_left[1] && iso_value in (1., 3., 7., 8., 12., 14.)) ||
                              (ind == gp.ind.b_right[1] && iso_value in (2., 3., 4., 11., 12., 13.))
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
    index_CL = locate_index_domain_boundary(gp,indices)
    # get the x coordinates of the contact line at the boundary
    x_CL_vec = get_contact_line_points_boundary(gp, indices, index_CL)
    # initialize stress vector and create a temporary array for relative x coordinates
    stress_vector = zeros(size(indices)); rel_x = similar(stress_vector)
    # create a temporary array for the bell function
    bell_function = zeros(size(indices))
    for i in findall(index_CL) # loop over all non-zeros
        if indices == gp.ind.b_bottom[1]
            rel_x .= gp.x[1, :] .- x_CL_vec[i]
            contact_angle = gp.α[1, i]
        elseif indices == gp.ind.b_top[1]
            rel_x .= gp.x[end, i] .- x_CL_vec[i]
            contact_angle = gp.α[end, i]
        elseif indices == gp.ind.b_left[1]
            rel_x .= gp.y[i, 1] .- x_CL_vec[i]
            contact_angle = gp.α[i, 1]
        elseif indices == gp.ind.b_right[1]
            rel_x .= gp.y[i, end] .- x_CL_vec[i]
            contact_angle = gp.α[i, end]
        end
        # compute the cosine value of the contact angle and apply stress only if it is less than the threshold
        cos_value = (cos(contact_angle) - cos(num.θe * π / 180))/num.Ca
        cos_boolean = isless(contact_angle, num.θe * π / 180)
        cos_value *= cos_boolean ? 1.0 : 0.0 
        # compute the bell function for the current x coordinates and add to the stress vector
        bell_function .= (1.0 .- tanh.(rel_x[:] ./ num.εCA) .^ 2) ./ num.εCA
        stress_vector .+= bell_function .* cos_value
    end
    return stress_vector
end


function compute_young_stress_internal(gp, num, mixed)

    function locate_CL_coordinates(gp, mixed)
        mixed_CL = CartesianIndex{2}[]
        for i in mixed
            iso_value = gp.iso[i]
            if iso_value != 0 && iso_value != 15
                push!(mixed_CL, i)
            end
        end
        return mixed_CL
    end

    mixed_CL = locate_CL_coordinates(gp, mixed)
    #@show mixed_CL

    # how to compute the stress with curvature ???

end

