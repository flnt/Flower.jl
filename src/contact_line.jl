## This file is part of the Flower project.
## It contains the routines and procedures to compute the Young stresses.
## It gives values to gp.Young as function of the BCs and the levelset (gp.u) on the domain.
## @author: R. A. S. Frantz (2023)

function update_Young_stress(gp,num)

    gp.Young[1,:] = compute_young_stress(gp,num,gp.ind.b_bottom[1])

end 

function locate_index_domain_boundary(gp, ind)
    index_CL = zeros(Bool, size(ind))
    if ind == gp.ind.b_bottom[1]
        @inbounds for i in eachindex(ind)
            iso_value = gp.iso[ind[i]]
            index_CL[i] = iso_value != 0.0 && iso_value in (1., 2., 6., 9., 13., 14.) ? true : index_CL[i]
        end
    elseif ind == gp.ind.b_top[1]
        @inbounds for i in eachindex(ind)
            iso_value = gp.iso[ind[i]]
            index_CL[i] = iso_value != 0.0 && iso_value in (4., 6., 8., 9., 13., 14.) ? true : index_CL[i]
        end
    elseif ind == gp.ind.b_left[1]
        @inbounds for i in eachindex(ind)
            iso_value = gp.iso[ind[i]]
            index_CL[i] = iso_value != 0.0 && iso_value in (1., 2., 6., 9., 13., 14.) ? true : index_CL[i]
        end
    elseif ind == gp.ind.b_right[1]
        @inbounds for i in eachindex(ind)
            iso_value = gp.iso[ind[i]]
            index_CL[i] = iso_value != 0.0 && iso_value in (1., 2., 6., 9., 13., 14.) ? true : index_CL[i]
        end
    end
    return index_CL
end


# function locate_index_levelset(mx)
#     for ii in mx # loop over mx points 
#         if mx[ii] != 0 or mx[ii] != 15
#           get the x, y coordinate
#           set point to true 
#         end
#     end
# end 

function get_cut_points(gp, indices, index_CL)
    x_Cl_vec = zeros(size(indices)) # zero vector of size indices
    for i in findall(index_CL) # loop over true values of index_CL
    if indices == gp.ind.b_bottom[1]
        x = gp.cut_points[1,i] # n,n matrix
        if x[1].y == -0.5
            x_Cl_vec[i] = x[1].x
        elseif x[2].y == -0.5
            x_Cl_vec[i] = x[2].x
        end
        x_Cl_vec[i] *= gp.dx[1,i]
        x_Cl_vec[i] += gp.x[1,i]
    elseif indices == gp.ind.b_top[1]
        x = gp.cut_points[end,i]
        if x[1].y == 0.5
            x_Cl_vec[i] = x[1].x
        elseif x[2].y == 0.5
            x_Cl_vec[i] = x[2].x
        end
        x_Cl_vec[i] *= gp.dx[end,i]
        x_Cl_vec[i] += gp.x[end,i]
    elseif indices == gp.ind.b_left[1]
        x = gp.cut_points[i,1]
        if x[1].x == -0.5
            x_Cl_vec[i] = x[1].y
        elseif x[2].x == -0.5
            x_Cl_vec[i] = x[2].y
        end
        x_Cl_vec[i] *= gp.dy[i,1]
        x_Cl_vec[i] += gp.y[i,1]
    elseif indices == gp.ind.b_right[1]
        x = gp.cut_points[i,end]
        if x[1].x == 0.5
            x_Cl_vec[i] = x[1].y
        elseif x[2].x == 0.5
            x_Cl_vec[i] = x[2].y
        end
        x_Cl_vec[i] *= gp.dy[i,end]
        x_Cl_vec[i] += gp.y[i,end]
    end
end
    return x_Cl_vec
end

# function compute_bell_function(gp,num,x_Cl_vec,index_CL,indices)
#     bell_function, rel_x = zeros(size(x_Cl_vec)), similar(x_Cl_vec)
#     for i in findall(index_CL) # loop over all non-zeros
#         @show i,x_Cl_vec[i]
#         if indices == gp.ind.b_bottom[1]
#             rel_x .= gp.x[1,:] .- x_Cl_vec[i]
#         elseif indices == gp.ind.b_top[1]
#             rel_x .= gp.x[end,i] .- x_Cl_vec[i]
#         elseif indices == gp.ind.b_left[1]
#             rel_x .= gp.y[i,1] .- x_Cl_vec[i]
#         elseif indices == gp.ind.b_right[1]
#             rel_x .= gp.y[i,end] .- x_Cl_vec[i]
#         end
#         bell_function .+= (1.0 .- tanh.(rel_x ./ num.εCA).^2)/num.εCA
#     end 
#     return bell_function
# end

# function test_bell(gp,num)
#     indices = gp.ind.b_bottom[1]
#     index_CL = locate_index_domain_boundary(gp,indices)
#     x_Cl_vec = get_cut_points(gp, indices, index_CL)
#     return compute_bell_function(gp,num,x_Cl_vec,index_CL,indices)
# end

# value to fill a0 
# function compute_young_stress(gp,num,indices)
#     index_CL = locate_index_domain_boundary(gp,indices)
#     x_Cl_vec = get_cut_points(gp, indices, index_CL)
#     bellf = compute_bell_function(gp,num,x_Cl_vec,index_CL,indices)
#     cos_values = cos.(gp.α[1, :]) .- cos(num.θe * π / 180)
#     divisor = 1.0 / num.Ca
#     stress = bellf[:] .* divisor .* cos_values[:]
#     return replace!(stress, NaN=>0.0)
# end

function compute_young_stress(gp, num, indices)
    index_CL = locate_index_domain_boundary(gp,indices)
    x_Cl_vec = get_cut_points(gp, indices, index_CL)
    stress = zeros(size(indices))
    rel_x = similar(stress)
    bell_function = zeros(size(indices))
    for i in findall(index_CL)  # loop over all non-zeros
        if indices == gp.ind.b_bottom[1]
            rel_x .= gp.x[1, :] .- x_Cl_vec[i]
            contact_angle = gp.α[1, i]
        elseif indices == gp.ind.b_top[1]
            rel_x .= gp.x[end, i] .- x_Cl_vec[i]
            contact_angle = gp.α[end, i]
        elseif indices == gp.ind.b_left[1]
            rel_x .= gp.y[i, 1] .- x_Cl_vec[i]
            contact_angle = gp.α[i, 1]
        elseif indices == gp.ind.b_right[1]
            rel_x .= gp.y[i, end] .- x_Cl_vec[i]
            contact_angle = gp.α[i, end]
        end
        cos_value = (cos(contact_angle) - cos(num.θe * π / 180))/num.Ca
        bell_function .= (1.0 .- tanh.(rel_x[:] ./ num.εCA) .^ 2) ./ num.εCA
        stress .+= bell_function .* cos_value
    end
    return stress
end