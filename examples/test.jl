
function compute_young_stress3(gp, num, indices)
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

compute_young_stress3(gp,num,gp.ind.b_bottom[1])