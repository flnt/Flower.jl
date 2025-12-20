
"""
    compute surface tension (CSF, Brackbill) with LS
"""
function compute_surface_tension_LS!(num,grid, grid_u, grid_v, opC_p, opC_u, opC_v,
    volume_fraction,levelset_one_fluid,volumic_surface_tension_u,volumic_surface_tension_v,tmp_vec_p,tmp_vec_p0,
    levelset_1D, levelset_heavyside_2D, normal_and_dirac_u, normal_and_dirac_v,
    normal_u, normal_v, curvature_u, curvature_v
    )


    #HACK reset centroids 
    # grid.LS[1].geoL.centroid .= Point(0.0, 0.0)
    # grid_u.LS[1].geoL.centroid .= Point(0.0, 0.0)
    # grid_v.LS[1].geoL.centroid .= Point(0.0, 0.0)

    # grid.LS[1].geoL.centroid.x .= 0.0
    # grid_u.LS[1].geoL.centroid.x .= 0.0
    # grid_v.LS[1].geoL.centroid.x .= 0.0

    # grid.LS[1].geoL.centroid.y .= 0.0
    # grid_u.LS[1].geoL.centroid.y .= 0.0
    # grid_v.LS[1].geoL.centroid.y .= 0.0


    # for grid_iter in [grid,grid_u,grid_v]
    #     for II in grid_iter.ind.inside
    #         grid_iter.LS[1].geoL.centroid[II] = Point(0.0,0.0)
    #     end
    # end


    #TODO   LS.mid_point[II] = ?

    iLS = 1

    # levelset_one_fluid = grid.LS[iLS].u saved or end
    # display(levelset_one_fluid)

    # print("\n size LS ",size(levelset_one_fluid))


    volumic_surface_tension_u .= 0.0
    volumic_surface_tension_v .= 0.0
    

    
    #region heavyside
    heavyside_epsilon = grid.dx[2,2]

    # levelset_1D = fnzeros(grid,num)

    if num.one_fluid_normal == 0

        for j in 1:grid.ny
            for i in 1:grid.nx
                pII = lexicographic(CartesianIndex(j,i),grid.ny)
                levelset_1D[pII] = levelset_heavyside(levelset_one_fluid[j,i],heavyside_epsilon)
            end
        end
    elseif num.one_fluid_normal == 1
        
        for j in 1:grid.ny
            for i in 1:grid.nx
                pII = lexicographic(CartesianIndex(j,i),grid.ny)
                levelset_1D[pII] = levelset_to_binary(levelset_one_fluid[j,i])
            end
        end

    end

    # vec1(levelset_1D,grid) .= vec(levelset_one_fluid)
    # levelset_1D

    #endregion heavyside

    # levelset_heavyside_2D = zeros(grid)

    levelset_heavyside_2D .= reshape(vec1(levelset_1D,grid), grid)

    # print("\n num ",num.nLS)

    PDI_status = @ccall "libpdi".PDI_multi_expose("write_one_fluid_levelset"::Cstring,
    "nstep"::Cstring, num.current_iter ::Ref{Clonglong}, PDI_OUT::Cint,
    "levelset_surface_tension"::Cstring, levelset_one_fluid::Ptr{Cdouble}, PDI_OUT::Cint,
    "levelset_heavyside"::Cstring, levelset_heavyside_2D::Ptr{Cdouble}, PDI_OUT::Cint,
    C_NULL::Ptr{Cvoid})::Cint



    # vecb_L(levelset_1D,grid) .= grid.LS[iLS].u[:,1] #90 degrees contact angle
    # vecb_B(levelset_1D,grid) .= grid.LS[iLS].u[1,:]
    # vecb_R(levelset_1D,grid) .= grid.LS[iLS].u[:,end]
    # vecb_T(levelset_1D,grid) .= grid.LS[iLS].u[end,:]


    vecb_L(levelset_1D,grid) .= levelset_heavyside_2D[:,1] #90 degrees contact angle
    vecb_B(levelset_1D,grid) .= levelset_heavyside_2D[1,:]
    vecb_R(levelset_1D,grid) .= levelset_heavyside_2D[:,end]
    vecb_T(levelset_1D,grid) .= levelset_heavyside_2D[end,:]
    
    # normal_and_dirac_u = zeros(grid_u)
    # normal_and_dirac_v = zeros(grid_v)
    # normal_u = zeros(grid_u)
    # normal_v = zeros(grid_v)

    compute_unit_normal(num,grid, grid_u, grid_v, 
    # opC_p, 
    opC_u, opC_v,levelset_1D,
    volume_fraction,
    # levelset_one_fluid,volumic_surface_tension_u,volumic_surface_tension_v,
    tmp_vec_p,tmp_vec_p0,
    normal_and_dirac_u,normal_and_dirac_v,
    normal_u,normal_v,
    )


    #region compute curvature Flower.jl
   
    # curvature_p = compute_curvature_cutcell_operator(opC_p,normal_and_dirac_u,normal_and_dirac_v)

    #endregion compute curvature Flower.jl


    # printstyled(color=:red, @sprintf "\n curvature_p with LS method")

    curvature_p = compute_curvature_levelset(levelset_one_fluid,grid.dx[2,2],grid.dy[2,2])

    # display(curvature_p)

  
    #region interpolate curvature from scalar to u and v grids
    # curvature_u = zeros(grid_u)
    # curvature_v = zeros(grid_v)

    # interpolate_scalar!(grid, grid_u, grid_v, curvature_p, curvature_u, curvature_v)

    interpolate_scalar_one_fluid_or_one_phase!(num,grid, grid_u, grid_v, curvature_p, curvature_u, curvature_v)


    #endregion interpolate curvature from scalar to u and v grids

   
    if num.constant_curvature >0
        curvature_u .= 4.0
        curvature_v .= 4.0
    end 

    # integrate the volumic surface tension
    volumic_surface_tension_u .= - num.sigma .* curvature_u .* normal_and_dirac_u
    volumic_surface_tension_v .= - num.sigma .* curvature_v .* normal_and_dirac_v

    # volumic_surface_tension_u .= - num.sigma .* curvature_u .* normal_u * grid.dx[2,2]
    # volumic_surface_tension_v .= - num.sigma .* curvature_v .* normal_v * grid.dx[2,2]

    # volumic_surface_tension_u .= - num.sigma .* curvature_u .* normal_u /grid.dx[2,2]
    # volumic_surface_tension_v .= - num.sigma .* curvature_v .* normal_v /grid.dx[2,2]

    # volumic_surface_tension_u .= - num.sigma .* curvature_u .* normal_u 
    # volumic_surface_tension_v .= - num.sigma .* curvature_v .* normal_v 

    

    iLSpdi = 1

    # norm
    # tmp_vec_p
    # tmp_vec_p0

    PDI_status = @ccall "libpdi".PDI_multi_expose("write_one_fluid_surface_tension"::Cstring,
    "nstep"::Cstring, num.current_iter ::Ref{Clonglong}, PDI_OUT::Cint,
    # "rho_one_fluid"::Cstring, rho_one_fluid::Ptr{Cdouble}, PDI_OUT::Cint,
    # "mu_one_fluid"::Cstring, mu_one_fluid::Ptr{Cdouble}, PDI_OUT::Cint,
    "volume_fraction"::Cstring, volume_fraction::Ptr{Cdouble}, PDI_OUT::Cint,
    "grad_u"::Cstring, normal_and_dirac_u::Ptr{Cdouble}, PDI_OUT::Cint,
    "grad_v"::Cstring, normal_and_dirac_v::Ptr{Cdouble}, PDI_OUT::Cint,
    "curvature_p"::Cstring, curvature_p::Ptr{Cdouble}, PDI_OUT::Cint,
    "curvature_u"::Cstring, curvature_u::Ptr{Cdouble}, PDI_OUT::Cint,
    "curvature_v"::Cstring, curvature_v::Ptr{Cdouble}, PDI_OUT::Cint,
    "volumic_surface_tension_u"::Cstring, volumic_surface_tension_u::Ptr{Cdouble}, PDI_OUT::Cint,
    "volumic_surface_tension_v"::Cstring, volumic_surface_tension_v::Ptr{Cdouble}, PDI_OUT::Cint,
    "normal_angle"::Cstring, grid.LS[iLSpdi].α::Ptr{Cdouble}, PDI_OUT::Cint,
    "normal_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
    "normal_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,  
    "mean_curvature"::Cstring, num.mean_curvature::Ref{Cdouble}, PDI_OUT::Cint,  
    C_NULL::Ptr{Cvoid})::Cint

    # print("\n volumic_surface_tension_u")
    # display(volumic_surface_tension_u)
    # display(volumic_surface_tension_v)

end


"""
Compute the curvature using the level set method
"""
 function compute_curvature_levelset(phi, Δx, Δy)

    
    # display(phi)


    # Compute the derivatives using central differences
    phi_x = (circshift(phi, (-1, 0)) - circshift(phi, (1, 0))) / (2 * Δx)
    phi_y = (circshift(phi, (0, -1)) - circshift(phi, (0, 1))) / (2 * Δy)

    # display(phi_x)

    # display(phi_y)


    phi_xx = (circshift(phi, (-1, 0)) - 2 * phi + circshift(phi, (1, 0))) / (Δx^2)
    phi_yy = (circshift(phi, (0, -1)) - 2 * phi + circshift(phi, (0, 1))) / (Δy^2)
    phi_xy = (circshift(phi, (-1, -1)) + circshift(phi, (1, 1)) - circshift(phi, (-1, 1)) - circshift(phi, (1, -1))) / (4 * Δx * Δy)

    # Compute the curvature
    numerator = phi_x.^2 .* phi_yy + phi_y.^2 .* phi_xx - 2 .* phi_x .* phi_y .* phi_xy
    denominator = (phi_x.^2 + phi_y.^2 ).^(3/2)
    curvature = numerator ./ denominator

    return curvature
end




"""
interpolates normal components to normalize
"""
function compute_unit_normal(num,grid, grid_u, grid_v, 
    # opC_p, 
    opC_u, opC_v,levelset_1D,
    volume_fraction,
    # levelset_one_fluid,volumic_surface_tension_u,volumic_surface_tension_v,
    tmp_vec_p,tmp_vec_p0,
    normal_and_dirac_u,
    normal_and_dirac_v,
    normal_u,
    normal_v,
    )

    # cell-averaged gradient of levelset_1D
    compute_grad_T_x_T_y_array_u_v_capacities!(num, grid, grid_u, grid_v, opC_u, opC_v, normal_and_dirac_u, normal_and_dirac_v, levelset_1D)


    # or
    # ∇ϕ_x = opC_u.AxT * opC_u.Rx * vec1(TD,grid) .+ opC_u.Gx_b * vecb(TD,grid)
    # ∇ϕ_y = opC_v.AyT * opC_v.Ry * vec1(TD,grid) .+ opC_v.Gy_b * vecb(TD,grid)

    #TODO check interpolation
    interpolate_grid_liquid_2!(num, grid, 
    grid_u.LS[end], grid_v.LS[end], 
    normal_and_dirac_u, normal_and_dirac_v, 
    tmp_vec_p, tmp_vec_p0) #compute normal x and y on scalar grid

    #compute normal x and y on scalar grid 

# interpolate_scalar_one_fluid_or_one_phase!

    # normalize 
    # normal_and_dirac_u .= normal_and_dirac_u / sqrt(normal_and_dirac_u**2+interpolate_v_to_u())

    # print("\n size normal u v ",size(normal_and_dirac_u)," ",size(normal_and_dirac_v))


    grid_u_x_full_2D = create_2D_grid_x(grid_u,false,true)
    grid_u_y_full_2D = create_2D_grid_y(grid_u,false,true)

       
    # print("\n grid_u_x_full_2D")
    # display(grid_u_x_full_2D)

    # print("\n grid_u_y_full_2D")
    # display(grid_u_y_full_2D)

    #TODO centroid wrong


    grid_v_x_full_2D = create_2D_grid_x(grid_v,true,false)
    

    grid_v_y_full_2D = create_2D_grid_y(grid_v,true,false)

    # print("\n grid_v_x_full_2D")
    # display(grid_v_x_full_2D)
    # print("\n grid_v_y_full_2D")
    # display(grid_v_y_full_2D)

    # print("\n grid_v_y_full_2D")
    # display(grid_v_y_full_2D)

    # print("\n test grid_v")
    # display(grid_v.y)

    # print("\n test grid_u")
    # display(grid_u.x)

    grad_v_full = similar(grid_v_x_full_2D)
    grad_v_full .= 0.0
    grad_v_full[:,2:end-1] = normal_and_dirac_v

    #TODO for now, assuming no contact, so n assumed to be 0

    # print("\n grad_v_full")
    # display(grad_v_full)

    grad_u_full = similar(grid_u_x_full_2D)
    grad_u_full .= 0.0
    grad_u_full[2:end-1,:] = normal_and_dirac_u

    # print("\n grad_u_full")
    # display(grad_u_full)
    
    x_centroid_u = grid_u.x .+ getproperty.(grid_u.LS[1].geoL.centroid, :x) .* grid_u.dx
    y_centroid_u = grid_u.y .+ getproperty.(grid_u.LS[1].geoL.centroid, :y) .* grid_u.dy

    # print("\nx_centroid_u ",x_centroid_u)
    # print("\ny_centroid_u ",y_centroid_u)

    x_centroid_v = grid_v.x .+ getproperty.(grid_v.LS[1].geoL.centroid, :x) .* grid_v.dx
    y_centroid_v = grid_v.y .+ getproperty.(grid_v.LS[1].geoL.centroid, :y) .* grid_v.dy


    #region normalize x component of grad  
    for j in 1:grid_u.ny
        for i in 1:grid_u.nx
            # print("\n normal u v ",j," ",i)

            interp_coord_x = x_centroid_u[j,i] #grid_u.x[j,i] #x_centroid_u[j,i]  but no interface

            interp_coord_y = y_centroid_u[j,i] # grid_u.y[j,i]

            x1 = grid_v_x_full_2D[j,i]
            x2 = grid_v_x_full_2D[j,i+1]
            y1 = grid_v_y_full_2D[j,i]
            y2 = grid_v_y_full_2D[j+1,i+1]

            Q11 = grad_v_full[j,i]
            Q12 = grad_v_full[j+1,i]
            Q21 = grad_v_full[j,i+1]
            Q22 = grad_v_full[j+1,i+1]

            #     (x2,y2)
            #  Q12   Q22
            # | .    . |
            # |   x    | interpolate based on four nodes of value Q11, Q12, Q21, Q22
            # | .    . |
            #   Q11  Q21
            #(x1,y1)
            interpolate_v_to_u = bilinear_interpolation(interp_coord_x, interp_coord_y, x1, y1, x2, y2, Q11, Q12, Q21, Q22)

            # printstyled(color=:green, @sprintf "\n i %.5i j %.5i x %.2e y %.2e x1 %.2e y1 %.2e x2 %.2e y2 %.2e Q11 %.2e Q12 %.2e Q21 %.2e Q22 %.2e\n" i j interp_coord_x interp_coord_y x1 y1 x2 y2 Q11 Q12 Q21 Q22)
            # print("\nvolume_fraction i ",i," j ",j," ",volume_fraction_face)
            
            norm = sqrt(normal_and_dirac_u[j,i]^2+interpolate_v_to_u^2)
            if norm != 0.0
                normal_u[j,i] = normal_and_dirac_u[j,i] / norm
            end
        end
    end

    #endregion normalize x component of grad 


    #region normalize y component of grad 

    # grad_u_full = similar(grid_u_x_full_2D)
    # grad_u_full .= 0.0
    # grad_u_full[2:end-1,:] = normal_and_dirac_u

    # print("\n grad_u_full")
    # display(grad_u_full)

    for j in 1:grid_v.ny
        for i in 1:grid_v.nx
            # print("\n normal u v ",j," ",i)

            interp_coord_x = x_centroid_v[j,i] #grid_v.x[j,i] #x_centroid_u[j,i]  but no interface

            interp_coord_y = y_centroid_v[j,i] #grid_v.y[j,i]

            x1 = grid_u_x_full_2D[j,i]
            x2 = grid_u_x_full_2D[j,i+1]
            y1 = grid_u_y_full_2D[j,i]
            y2 = grid_u_y_full_2D[j+1,i+1]

            Q11 = grad_u_full[j,i]
            Q12 = grad_u_full[j+1,i]
            Q21 = grad_u_full[j,i+1]
            Q22 = grad_u_full[j+1,i+1]

            #     (x2,y2)
            #  Q12   Q22
            # | .    . |
            # |   x    | interpolate based on four nodes of value Q11, Q12, Q21, Q22
            # | .    . |
            #   Q11  Q21
            #(x1,y1)
            interpolate_u_to_v = bilinear_interpolation(interp_coord_x, interp_coord_y, x1, y1, x2, y2, Q11, Q12, Q21, Q22)
            
            # printstyled(color=:green, @sprintf "\n i %.5i j %.5i x %.2e y %.2e x1 %.2e y1 %.2e x2 %.2e y2 %.2e Q11 %.2e Q12 %.2e Q21 %.2e Q22 %.2e\n" i j interp_coord_x interp_coord_y x1 y1 x2 y2 Q11 Q12 Q21 Q22)

            norm = sqrt(normal_and_dirac_v[j,i]^2+interpolate_u_to_v^2)
            if norm !=0.0
                normal_v[j,i] = normal_and_dirac_v[j,i] / norm
            end
        end
    end

    #endregion normalize y component of grad 

    # printstyled(color=:green, @sprintf "\n normal_and_dirac_u")
    # display(normal_and_dirac_u)
    # display(normal_and_dirac_v)

    #region integrate gradient of levelset_1D
    #fill normal_and_dirac_u and normal_and_dirac_v with gradient of levelset_1D
    compute_grad_T_x_T_y_array_u_v_capacities_cell_integrated!(num, grid, grid_u, grid_v, opC_u, opC_v, 
    normal_and_dirac_u, normal_and_dirac_v, 
    levelset_1D) 
    #endregion integrate gradient of levelset_1D

end

"""
interpolates normal components to normalize
"""
function compute_unit_normal_debug(num,grid, grid_u, grid_v, 
    # opC_p, 
    opC_u, opC_v,levelset_1D,
    volume_fraction,
    # levelset_one_fluid,volumic_surface_tension_u,volumic_surface_tension_v,
    tmp_vec_p,tmp_vec_p0,
    normal_and_dirac_u,
    normal_and_dirac_v
    )



    compute_grad_T_x_T_y_array_u_v_capacities!(num, grid, grid_u, grid_v, opC_u, opC_v, normal_and_dirac_u, normal_and_dirac_v, levelset_1D)

    #TODO check interpolation
    interpolate_grid_liquid_2!(num, grid, grid_u.LS[end], grid_v.LS[end], normal_and_dirac_u, normal_and_dirac_v, tmp_vec_p, tmp_vec_p0) #compute normal x and y on scalar grid
    #compute normal x and y on scalar grid 

    # normalize 
    # normal_and_dirac_u .= normal_and_dirac_u / sqrt(normal_and_dirac_u**2+interpolate_v_to_u())

    print("\n size normal u v ",size(normal_and_dirac_u)," ",size(normal_and_dirac_v))


    grid_u_x_full_2D = create_2D_grid_x(grid_u,false,true)
    grid_u_y_full_2D = create_2D_grid_y(grid_u,false,true)

       
    print("\n grid_u_x_full_2D")
    display(grid_u_x_full_2D)

    print("\n grid_u_y_full_2D")
    display(grid_u_y_full_2D)

    #TODO centroid wrong


    grid_v_x_full_2D = create_2D_grid_x(grid_v,true,false)
    

    grid_v_y_full_2D = create_2D_grid_y(grid_v,true,false)

    print("\n grid_v_x_full_2D")
    display(grid_v_x_full_2D)
    print("\n grid_v_y_full_2D")
    display(grid_v_y_full_2D)

    # print("\n grid_v_y_full_2D")
    # display(grid_v_y_full_2D)

    # print("\n test grid_v")
    # display(grid_v.y)

    # print("\n test grid_u")
    # display(grid_u.x)

    grad_v_full = similar(grid_v_x_full_2D)
    grad_v_full .= 0.0
    grad_v_full[:,2:end-1] = normal_and_dirac_v

    #TODO for now, assuming no contact, so n assumed to be 0

    print("\n grad_v_full")
    display(grad_v_full)

    grad_u_full = similar(grid_u_x_full_2D)
    grad_u_full .= 0.0
    grad_u_full[2:end-1,:] = normal_and_dirac_u

    print("\n grad_u_full")
    display(grad_u_full)
    
    x_centroid_u = grid_u.x .+ getproperty.(grid_u.LS[1].geoL.centroid, :x) .* grid_u.dx
    y_centroid_u = grid_u.y .+ getproperty.(grid_u.LS[1].geoL.centroid, :y) .* grid_u.dy

    # print("\nx_centroid_u ",x_centroid_u)
    # print("\ny_centroid_u ",y_centroid_u)

    x_centroid_v = grid_v.x .+ getproperty.(grid_v.LS[1].geoL.centroid, :x) .* grid_v.dx
    y_centroid_v = grid_v.y .+ getproperty.(grid_v.LS[1].geoL.centroid, :y) .* grid_v.dy


    #region normalize x component of grad  
    for j in 1:grid_u.ny
        for i in 1:grid_u.nx
            print("\n normal u v ",j," ",i)

            interp_coord_x = x_centroid_u[j,i] #grid_u.x[j,i] #x_centroid_u[j,i]  but no interface

            interp_coord_y = y_centroid_u[j,i] # grid_u.y[j,i]

            x1 = grid_v_x_full_2D[j,i]
            x2 = grid_v_x_full_2D[j,i+1]
            y1 = grid_v_y_full_2D[j,i]
            y2 = grid_v_y_full_2D[j+1,i+1]

            Q11 = grad_v_full[j,i]
            Q12 = grad_v_full[j+1,i]
            Q21 = grad_v_full[j,i+1]
            Q22 = grad_v_full[j+1,i+1]

            #     (x2,y2)
            #  Q12   Q22
            # | .    . |
            # |   x    | interpolate based on four nodes of value Q11, Q12, Q21, Q22
            # | .    . |
            #   Q11  Q21
            #(x1,y1)
            interpolate_v_to_u = bilinear_interpolation(interp_coord_x, interp_coord_y, x1, y1, x2, y2, Q11, Q12, Q21, Q22)

            printstyled(color=:green, @sprintf "\n i %.5i j %.5i x %.2e y %.2e x1 %.2e y1 %.2e x2 %.2e y2 %.2e Q11 %.2e Q12 %.2e Q21 %.2e Q22 %.2e\n" i j interp_coord_x interp_coord_y x1 y1 x2 y2 Q11 Q12 Q21 Q22)
            # print("\nvolume_fraction i ",i," j ",j," ",volume_fraction_face)
            
            norm = sqrt(normal_and_dirac_u[j,i]^2+interpolate_v_to_u^2)
            if norm != 0.0
                normal_and_dirac_u[j,i] = normal_and_dirac_u[j,i] / norm
            end
        end
    end

    #endregion normalize x component of grad 


    #region normalize y component of grad 

    # grad_u_full = similar(grid_u_x_full_2D)
    # grad_u_full .= 0.0
    # grad_u_full[2:end-1,:] = normal_and_dirac_u

    # print("\n grad_u_full")
    # display(grad_u_full)

    for j in 1:grid_v.ny
        for i in 1:grid_v.nx
            print("\n normal u v ",j," ",i)

            interp_coord_x = x_centroid_v[j,i] #grid_v.x[j,i] #x_centroid_u[j,i]  but no interface

            interp_coord_y = y_centroid_v[j,i] #grid_v.y[j,i]

            x1 = grid_u_x_full_2D[j,i]
            x2 = grid_u_x_full_2D[j,i+1]
            y1 = grid_u_y_full_2D[j,i]
            y2 = grid_u_y_full_2D[j+1,i+1]

            Q11 = grad_u_full[j,i]
            Q12 = grad_u_full[j+1,i]
            Q21 = grad_u_full[j,i+1]
            Q22 = grad_u_full[j+1,i+1]

            #     (x2,y2)
            #  Q12   Q22
            # | .    . |
            # |   x    | interpolate based on four nodes of value Q11, Q12, Q21, Q22
            # | .    . |
            #   Q11  Q21
            #(x1,y1)
            interpolate_u_to_v = bilinear_interpolation(interp_coord_x, interp_coord_y, x1, y1, x2, y2, Q11, Q12, Q21, Q22)
            
            printstyled(color=:green, @sprintf "\n i %.5i j %.5i x %.2e y %.2e x1 %.2e y1 %.2e x2 %.2e y2 %.2e Q11 %.2e Q12 %.2e Q21 %.2e Q22 %.2e\n" i j interp_coord_x interp_coord_y x1 y1 x2 y2 Q11 Q12 Q21 Q22)

            norm = sqrt(normal_and_dirac_v[j,i]^2+interpolate_u_to_v^2)
            if norm !=0.0
                normal_and_dirac_v[j,i] = normal_and_dirac_v[j,i] / norm
            end
        end
    end

    #endregion normalize y component of grad 

    printstyled(color=:green, @sprintf "\n normal_and_dirac_u")
    display(normal_and_dirac_u)
    display(normal_and_dirac_v)

end


function compute_curvature_cutcell_operator(opC_p,normal_and_dirac_u,normal_and_dirac_v)
     # TODO border ???

    # divergence of velocity explicit
    curvature_p_1D = opC_p.AxT * vec(normal_and_dirac_u) .+ opC_p.AyT * vec(normal_and_dirac_v) 


    compute_curvature_border = 0

    if compute_curvature_border > 0
        # normal_and_dirac_u_1D = ...
        # normal_and_dirac_v_1D = ...
        curvature_p_1D .+= opC_p.Gx_b * vecb(normal_and_dirac_u_1D,grid_u) .+ opC_p.Gy_b * vecb(normal_and_dirac_v_1D,grid_v)
    end

    curvature_p = reshape(veci(curvature_p_1D,grid,1), grid)


    # printstyled(color=:green, @sprintf "\n curvature_p")

    # display(curvature_p)


    # ∇ϕ_x = iMu * ∇ϕ_x
    # ∇ϕ_y = iMv * ∇ϕ_y
    iM = Diagonal(inv_weight_eps2.(num.epsilon_mode,num.epsilon_vol,vec(grid.LS[end].geoL.dcap[:,:,5])))
    # iM = Diagonal(inv_weight_eps2.(num.epsilon_mode,num.epsilon_vol,opC_p.M.diag))

   

    # printstyled(color=:green, @sprintf "\n volume")
    # # print("\n vol ",(grid.dx[2,2])^2)
    # # display(grid.LS[end].geoL.dcap[:,:,5])
    # display(grid.LS[1].geoL.dcap[:,:,5])


    # curvature_p = iM * curvature_p
    # iM = Diagonal(inv_weight_eps2.(num.epsilon_mode,num.epsilon_vol,vec(grid.LS[end].geoL.dcap[:,:,5])))
    
    curvature_p = curvature_p ./ grid.LS[1].geoL.dcap[:,:,5]
    
    # printstyled(color=:green, @sprintf "\n curvature_p")
    # display(curvature_p)

    #no interface: one-fluid model      
    # for iLS in 1:nLS
    #     if !is_navier(bc_int[iLS]) && !is_navier_cl(bc_int[iLS])
    #         curvature_p .+= opC_p.Gx[iLS] * veci(u_predictionD,grid_u,iLS+1) .+ 
    #                 opC_p.Gy[iLS] * veci(v_predictionD,grid_v,iLS+1)
    #     end
    # end


    # divergence of velocity explicit
    # Duv = opC_p.AxT * vec1(u_predictionD,grid_u) .+ opC_p.Gx_b * vecb(u_predictionD,grid_u) .+
    #       opC_p.AyT * vec1(v_predictionD,grid_v) .+ opC_p.Gy_b * vecb(v_predictionD,grid_v)
    # for iLS in 1:nLS
    #     if !is_navier(bc_int[iLS]) && !is_navier_cl(bc_int[iLS])
    #         Duv .+= opC_p.Gx[iLS] * veci(u_predictionD,grid_u,iLS+1) .+ 
    #                 opC_p.Gy[iLS] * veci(v_predictionD,grid_v,iLS+1)
    #     end
    # end
end