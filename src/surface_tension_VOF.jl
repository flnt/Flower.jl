
"""
    compute surface tension (Brackbill) 
"""
function compute_surface_tension_VOF!(num,grid, grid_u, grid_v, opC_p, opC_u, opC_v,volume_fraction,levelset_one_fluid,volumic_surface_tension_u,volumic_surface_tension_v,tmp_vec_p,tmp_vec_p0)

    volumic_surface_tension_u .= 0.0
    volumic_surface_tension_v .= 0.0
    

    volume_fraction_1D = zeros(grid)

    # print("\n num ",num.nLS)
    # print("\n size ",size(volume_fraction_1D))

    if num.smooth_VOF >0
        # kernel_size = 3 # Size of the smoothing kernel
        
        print("\n smooth VOF")

        smoothed_volume_fraction = similar(volume_fraction)

        smooth_vof_2d!(grid, volume_fraction, num.smooth_VOF,smoothed_volume_fraction)

        #we do not have ghost cells, so we recopy the values at 2 to 1
        smoothed_volume_fraction[:,1] = smoothed_volume_fraction[:,2]
        smoothed_volume_fraction[:,end] = smoothed_volume_fraction[:,end-1]
        smoothed_volume_fraction[1,:] = smoothed_volume_fraction[2,:]
        smoothed_volume_fraction[end,:] = smoothed_volume_fraction[end-1,:]

        display(volume_fraction)

        display(smoothed_volume_fraction)

           
        vec1(volume_fraction_1D,grid) .= vec(smoothed_volume_fraction) 

        PDI_status = @ccall "libpdi".PDI_multi_expose("write_one_fluid_smoothed_volume_fraction"::Cstring,
        "nstep"::Cstring, num.current_iter ::Ref{Clonglong}, PDI_OUT::Cint,
        # "rho_one_fluid"::Cstring, rho_one_fluid::Ptr{Cdouble}, PDI_OUT::Cint,
        # "mu_one_fluid"::Cstring, mu_one_fluid::Ptr{Cdouble}, PDI_OUT::Cint,
        "smoothed_volume_fraction"::Cstring, smoothed_volume_fraction::Ptr{Cdouble}, PDI_OUT::Cint,
        # "grad_u"::Cstring, normal_and_dirac_u::Ptr{Cdouble}, PDI_OUT::Cint,
        # "grad_v"::Cstring, normal_and_dirac_v::Ptr{Cdouble}, PDI_OUT::Cint,
        # "curvature_p"::Cstring, curvature_p::Ptr{Cdouble}, PDI_OUT::Cint,
        # "curvature_u"::Cstring, curvature_u::Ptr{Cdouble}, PDI_OUT::Cint,
        # "curvature_v"::Cstring, curvature_v::Ptr{Cdouble}, PDI_OUT::Cint,
        # "volumic_surface_tension_u"::Cstring, volumic_surface_tension_u::Ptr{Cdouble}, PDI_OUT::Cint,
        # "volumic_surface_tension_v"::Cstring, volumic_surface_tension_v::Ptr{Cdouble}, PDI_OUT::Cint,
        # "normal_angle"::Cstring, grid.LS[iLSpdi].α::Ptr{Cdouble}, PDI_OUT::Cint,
        # "normal_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
        # "normal_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,  
        C_NULL::Ptr{Cvoid})::Cint

        vecb_L(volume_fraction_1D,grid) .= smoothed_volume_fraction[:,1] #90 degrees contact angle
        vecb_B(volume_fraction_1D,grid) .= smoothed_volume_fraction[1,:]
        vecb_R(volume_fraction_1D,grid) .= smoothed_volume_fraction[:,end]
        vecb_T(volume_fraction_1D,grid) .= smoothed_volume_fraction[end,:]

    else
        vec1(volume_fraction_1D,grid) .= vec(volume_fraction) 
        vecb_L(volume_fraction_1D,grid) .= volume_fraction[:,1] #90 degrees contact angle
        vecb_B(volume_fraction_1D,grid) .= volume_fraction[1,:]
        vecb_R(volume_fraction_1D,grid) .= volume_fraction[:,end]
        vecb_T(volume_fraction_1D,grid) .= volume_fraction[end,:]
    end





    # vecb_L(volume_fraction_1D,grid) .= volume_fraction[:,1] #90 degrees contact angle
    # vecb_B(volume_fraction_1D,grid) .= volume_fraction[1,:]
    # vecb_R(volume_fraction_1D,grid) .= volume_fraction[:,end]
    # vecb_T(volume_fraction_1D,grid) .= volume_fraction[end,:]

    # vecb_L(volume_fraction_1D,grid) .= vec1(volume_fraction_1D,grid)[:,1] #90 degrees contact angle
    # vecb_B(volume_fraction_1D,grid) .= vec1(volume_fraction_1D,grid)[1,:]
    # vecb_R(volume_fraction_1D,grid) .= vec1(volume_fraction_1D,grid)[:,end]
    # vecb_T(volume_fraction_1D,grid) .= vec1(volume_fraction_1D,grid)[end,:]


    # print("\n volume_fraction_1D ",minimum(volume_fraction_1D), " ", maximum(volume_fraction_1D))


    
    normal_and_dirac_u = zeros(grid_u)
    normal_and_dirac_v = zeros(grid_v)

    compute_grad_T_x_T_y_array_u_v_capacities!(num, grid, grid_u, grid_v, opC_u, opC_v, normal_and_dirac_u, normal_and_dirac_v, volume_fraction_1D)

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


    # TODO border ???

    # divergence of velocity explicit
    curvature_p = opC_p.AxT * vec(normal_and_dirac_u) .+ opC_p.AyT * vec(normal_and_dirac_v) 
    
    compute_curvature_border = 0

    if compute_curvature_border > 0
        # normal_and_dirac_u_1D = ...
        # normal_and_dirac_v_1D = ...
        curvature_p .+= opC_p.Gx_b * vecb(normal_and_dirac_u_1D,grid_u) .+ opC_p.Gy_b * vecb(normal_and_dirac_v_1D,grid_v)
    end

    #no interface: one-fluid model      
    # for iLS in 1:nLS
    #     if !is_navier(bc_int[iLS]) && !is_navier_cl(bc_int[iLS])
    #         curvature_p .+= opC_p.Gx[iLS] * veci(ucorrD,grid_u,iLS+1) .+ 
    #                 opC_p.Gy[iLS] * veci(v_predictionD,grid_v,iLS+1)
    #     end
    # end


    # divergence of velocity explicit
    # Duv = opC_p.AxT * vec1(ucorrD,grid_u) .+ opC_p.Gx_b * vecb(ucorrD,grid_u) .+
    #       opC_p.AyT * vec1(v_predictionD,grid_v) .+ opC_p.Gy_b * vecb(v_predictionD,grid_v)
    # for iLS in 1:nLS
    #     if !is_navier(bc_int[iLS]) && !is_navier_cl(bc_int[iLS])
    #         Duv .+= opC_p.Gx[iLS] * veci(ucorrD,grid_u,iLS+1) .+ 
    #                 opC_p.Gy[iLS] * veci(v_predictionD,grid_v,iLS+1)
    #     end
    # end


    curvature_u = zeros(grid_u)
    curvature_v = zeros(grid_v)

    interpolate_scalar!(grid, grid_u, grid_v, curvature_p, curvature_u, curvature_v)

    volumic_surface_tension_u = - num.sigma .* curvature_u .* normal_and_dirac_u
    volumic_surface_tension_v = - num.sigma .* curvature_v .* normal_and_dirac_v
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
    C_NULL::Ptr{Cvoid})::Cint

end
