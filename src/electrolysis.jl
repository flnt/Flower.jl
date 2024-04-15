#From set_heat in heat_coupled.jl, poisson.jl
"""    
set convection and diffusion for scalar
"""
function set_scalar_transport!(bc_type, num, grid, op, geo, ph, θd, BC_T, MIXED, projection,
    A, B,
    op_conv, grid_u, geo_u, grid_v, geo_v,
    periodic_x, periodic_y, convection, ls_advection, BC_int, diffusion_coeff_scal)
    @unpack τ, aniso = num
    @unpack nx, ny, dx, dy, ind  = grid
    @unpack all_indices, inside, b_left, b_bottom, b_right, b_top = ind
    @unpack Bx, By, BxT, ByT, Hx, Hy, HxT, HyT, M, iMx, iMy, χ = op
    @unpack CT, CUTCT = op_conv
    @unpack u, v, uD, vD = ph

    ni = nx * ny
    nb = 2 * nx + 2 * ny
    nt = 2 * ni + nb

    if is_dirichlet(bc_type)
        # printstyled(color=:green, @sprintf "\n Dirichlet : %.2e\n" bc_type.val )
        __a0 = bc_type.val
        __a1 = -1.0
        __b = 0.0
    elseif is_neumann(bc_type)
        __a0 = bc_type.val
        __a1 = 0.0
        __b = 1.0
    elseif is_robin(bc_type)
        __a0 = bc_type.val
        __a1 = -1.0
        __b = 1.0
    elseif is_stefan(bc_type)
        __a0 = θd
        __a1 = -1.0
        __b = 0.0
    elseif is_wall(bc_type)
        __a0 = bc_type.val
        __a1 = -1.0
        __b = 0.0
    else
        __a0 = bc_type.val
        __a1 = -1.0
        __b = 0.0
    end

    # Flags with BCs
    a0 = ones(grid) .* __a0
    # if aniso
    #     apply_anisotropy(num, grid, grid.LS[1].κ, a0, MIXED, projection)
    # else
    #     apply_curvature(num, grid, grid.LS[1].κ, a0, all_indices)
    # end
    _a1 = ones(grid) .* __a1
    a1 = Diagonal(vec(_a1))
    _b = ones(grid) .* __b
    b = Diagonal(vec(_b))

    a0_b = zeros(nb)
    _a1_b = zeros(nb)
    _b_b = zeros(nb)
    set_borders!(grid, grid.LS[1].cl, grid.LS[1].u, a0_b, _a1_b, _b_b, BC_T, num.n_ext_cl)
    a1_b = Diagonal(vec(_a1_b))
    b_b = Diagonal(vec(_b_b))

    if convection
        HT = zeros(grid)
        @inbounds @threads for II in vcat(b_left[1], b_bottom[1], b_right[1], b_top[1])
            HT[II] = distance(grid.LS[1].mid_point[II], geo.centroid[II], dx[II], dy[II])
        end    
        bcTx, bcTy = set_bc_bnds(dir, a0, HT, BC_T)
    
        bnds_u = [grid_u.ind.b_left[1], grid_u.ind.b_bottom[1], grid_u.ind.b_right[1], grid_u.ind.b_top[1]]
        bnds_v = [grid_v.ind.b_left[1], grid_v.ind.b_bottom[1], grid_v.ind.b_right[1], grid_v.ind.b_top[1]]
        Δu = [grid_u.dx[1,1], grid_u.dy[1,1], grid_u.dx[end,end], grid_u.dy[end,end]] .* 0.5
        Δv = [grid_v.dx[1,1], grid_v.dy[1,1], grid_v.dx[end,end], grid_v.dy[end,end]] .* 0.5
        
        Hu = zeros(grid_u)
        for i in eachindex(bnds_u)
            for II in bnds_u[i]
                Hu[II] = Δu[i]
            end
        end

        Hv = zeros(grid_v)
        for i in eachindex(bnds_v)
            for II in bnds_v[i]
                Hv[II] = Δv[i]
            end
        end
    
        bcU = zeros(grid_u)
        bcU .= reshape(vec2(uD,grid_u), grid_u)
        bcU[:,1] .= vecb_L(uD, grid_u)
        bcU[1,:] .= vecb_B(uD, grid_u)
        bcU[:,end] .= vecb_R(uD, grid_u)
        bcU[end,:] .= vecb_T(uD, grid_u)
        
        bcV = zeros(grid_v)
        bcV .= reshape(vec2(vD,grid_v), grid_v)
        bcV[:,1] .= vecb_L(vD, grid_v)
        bcV[1,:] .= vecb_B(vD, grid_v)
        bcV[:,end] .= vecb_R(vD, grid_v)
        bcV[end,:] .= vecb_T(vD, grid_v)

        scalar_convection!(dir, CT, CUTCT, u, v, bcTx, bcTy, bcU, bcV, geo.dcap, ny, 
            BC_T, inside, b_left[1], b_bottom[1], b_right[1], b_top[1]
        )
    end

    if ls_advection
        update_all_ls_data(num, grid, grid_u, grid_v, BC_int, periodic_x, periodic_y, false)

        # Mass matrices
        M.diag .= vec(geo.dcap[:,:,5])
        Mx = zeros(ny,nx+1)
        for II in ind.all_indices
            Mx[II] = geo.dcap[II,8]
        end
        for II in ind.b_right[1]
            Mx[δx⁺(II)] = geo.dcap[II,10]
        end
        My = zeros(ny+1,nx)
        for II in ind.all_indices
            My[II] = geo.dcap[II,9]
        end
        for II in ind.b_top[1]
            My[δy⁺(II)] = geo.dcap[II,11]
        end
        iMx.diag .= 1. ./ (vec(Mx) .+ eps(0.01))
        iMy.diag .= 1. ./ (vec(My) .+ eps(0.01))

        # Discrete gradient and divergence operators
        divergence_B!(BxT, ByT, geo.dcap, ny, ind.all_indices)
        mat_assign!(Bx, sparse(-BxT'))
        mat_assign!(By, sparse(-ByT'))

        # Matrices for interior BCs
        for iLS in 1:num.nLS
            bc_matrix!(grid, Hx[iLS], Hy[iLS], geo.dcap, geo.dcap, ny, ind.all_indices)

            mat_assign_T!(HxT[iLS], sparse(Hx[iLS]'))
            mat_assign_T!(HyT[iLS], sparse(Hy[iLS]'))

            periodic_bcs!(grid, Bx, By, Hx[iLS], Hy[iLS], periodic_x, periodic_y)

            χx = (geo.dcap[:,:,3] .- geo.dcap[:,:,1]) .^ 2
            χy = (geo.dcap[:,:,4] .- geo.dcap[:,:,2]) .^ 2
            χ[iLS].diag .= sqrt.(vec(χx .+ χy))
        end
        mat_assign!(BxT, sparse(-Bx'))
        mat_assign!(ByT, sparse(-By'))

        # Matrices for borders BCs
        set_boundary_indicator!(grid, geo, geo, op)
        mass_matrix_borders!(ind, op.iMx_b, op.iMy_b, op.iMx_bd, op.iMy_bd, geo.dcap, ny)
        bc_matrix_borders!(grid, ind, op.Hx_b, op.Hy_b, geo.dcap)
        mat_assign_T!(op.HxT_b, sparse(op.Hx_b'))
        mat_assign_T!(op.HyT_b, sparse(op.Hy_b'))
        periodic_bcs_borders!(grid, op.Hx_b, op.Hy_b, periodic_x, periodic_y)
    end

    LT = BxT * iMx * Bx .+ ByT * iMy * By
    LD = BxT * iMx * Hx[1] .+ ByT * iMy * Hy[1]
    LD_b = BxT * op.iMx_b * op.Hx_b .+ ByT * op.iMy_b * op.Hy_b

    # Implicit part of heat equation
    A[1:ni,1:ni] = pad_crank_nicolson(M .- 0.5 .* τ .* diffusion_coeff_scal .* LT, grid, τ)
    A[1:ni,ni+1:2*ni] = - 0.5 .* τ .* diffusion_coeff_scal .* LD
    A[1:ni,end-nb+1:end] = - 0.5 .* τ .* diffusion_coeff_scal .* LD_b

    # Interior BC
    A[ni+1:2*ni,1:ni] = b * (HxT[1] * iMx * Bx .+ HyT[1] * iMy * By)
    A[ni+1:2*ni,ni+1:2*ni] = pad(b * (HxT[1] * iMx * Hx[1] .+ HyT[1] * iMy * Hy[1]) .- χ[1] * a1)

    # Border BCs
    A[end-nb+1:end,1:ni] = b_b * (op.HxT_b * op.iMx_b' * Bx .+ op.HyT_b * op.iMy_b' * By)
    A[end-nb+1:end,ni+1:2*ni] = b_b * (op.HxT_b * op.iMx_b' * Hx[1] .+ op.HyT_b * op.iMy_b' * op.Hy[1])
    A[end-nb+1:end,end-nb+1:end] = pad(b_b * (op.HxT_b * op.iMx_bd * op.Hx_b .+ op.HyT_b * op.iMy_bd * op.Hy_b) .- op.χ_b * a1_b, 4.0)

    # Explicit part of heat equation
    B[1:ni,1:ni] = M .+ 0.5 .* τ .* diffusion_coeff_scal .* LT
    B[1:ni,ni+1:2*ni] = 0.5 .* τ .* diffusion_coeff_scal .* LD
    B[1:ni,end-nb+1:end] = 0.5 .* τ .* diffusion_coeff_scal .* LD_b

    rhs = fnzeros(grid, num)
    if convection
        vec1(rhs,grid) .-= τ .* CUTCT
    end
    vec2(rhs,grid) .+= χ[1] * vec(a0)
    vecb(rhs,grid) .+= op.χ_b * vec(a0_b)

    return rhs
end


"""
Butler-Volmer model, relation between exchande current and potential at electrode, with effects of concentration
c0_H2: reference concentration
"""
function butler_volmer_concentration(alphaa,alphac,c_H2,c0_H2,c_H2O,c0_H2O,c_KOH,c0_KOH,Faraday,i0,phi_ele,phi_ele1,Ru,temperature0)
    eta = phi_ele1 - phi_ele
    i_current = i0*(sqrt(c_H2/c0_H2)*(c_KOH/c0_KOH)*exp(alphaa*Faraday*eta/(Ru*temperature0))-(c_H2O/c0_H2O)*exp(-alphac*Faraday*eta/(Ru*temperature0)))
    return i_current
end

"""
Butler-Volmer model, relation between exchande current and potential at electrode without concentration
"""
function butler_volmer_no_concentration(alphaa,alphac,Faraday,i0,phi_ele,phi_ele1,Ru,temperature0)
    eta = phi_ele1 - phi_ele
    i_current = i0*(exp(alphaa*Faraday*eta/(Ru*temperature0))-exp(-alphac*Faraday*eta/(Ru*temperature0)))
    return i_current
end



"""
From Stefan_velocity!
"""
function electrolysis_velocity!(num, grid, LS, V, TL, MIXED, periodic_x, periodic_y, concentration_scal)
    # @unpack θd, ϵ_κ, ϵ_V, m, θ₀, aniso = num
    @unpack geoS, geoL, κ = LS

    V .= 0
    @inbounds @threads for II in MIXED
        # ϵ_c = ifelse(aniso, anisotropy(ϵ_κ, m, geoS.projection[II].angle, θ₀), ϵ_κ)
        # ϵ_v = ifelse(aniso, anisotropy(ϵ_V, m, geoS.projection[II].angle, θ₀), ϵ_V)
        # θ_d = (θd - ϵ_c*κ[II] - ϵ_v*V[II])

        θ_d = concentration_scal

        # dTS = 0.
        dTL = 0.
        # if geoS.projection[II].flag
        #     T_1, T_2 = interpolated_temperature(grid, geoS.projection[II].angle, geoS.projection[II].point1, geoS.projection[II].point2, TS, II, periodic_x, periodic_y)
        #     dTS = normal_gradient(geoS.projection[II].d1, geoS.projection[II].d2, T_1, T_2, θ_d)
        # else
        #     T_1 = interpolated_temperature(grid, geoS.projection[II].angle, geoS.projection[II].point1, TS, II, periodic_x, periodic_y)
        #     dTS = normal_gradient(geoS.projection[II].d1, T_1, θ_d)
        # end
        if geoL.projection[II].flag
            T_1, T_2 = interpolated_temperature(grid, geoL.projection[II].angle, geoL.projection[II].point1, geoL.projection[II].point2, TL, II, periodic_x, periodic_y)
            dTL = normal_gradient(geoL.projection[II].d1, geoL.projection[II].d2, T_1, T_2, θ_d)
        else
            T_1 = interpolated_temperature(grid, geoL.projection[II].angle, geoL.projection[II].point1, TL, II, periodic_x, periodic_y)
            dTL = normal_gradient(geoL.projection[II].d1, T_1, θ_d)
        end
        V[II] = dTL + dTS
    end
    return nothing
end

"""
From update_free_surface_velocity and update_stefan_velocity
"""
function update_free_surface_velocity_electrolysis(num, grid, grid_u, grid_v, iLS, uD, vD, periodic_x, periodic_y, Vmean, c_L,diffusion_coeff_scal,concentration_scal,opC)
    @unpack MWH2,rho1,rho2=num
    @unpack χ,=opC
    #TODO check sign 
    electrolysis_velocity!(num, grid, grid.LS[iLS], grid.V, c_L, grid.LS[iLS].MIXED, periodic_x, periodic_y, concentration_scal)
    # grid.V[grid.LS[iLS].MIXED] .*= 1. ./ λ

    #H2
    #TODO *surface: *χ or χ_b ?

    grid.V[grid.LS[iLS].MIXED] .*=-(1/rho1-1/rho2)*diffusion_coeff_scal[1]*MWH2*χ  


    if Vmean
        a = mean(grid.V[grid.LS[iLS].MIXED])
        grid.V[grid.LS[iLS].MIXED] .= a
    end

    grid_u.V .= reshape(veci(uD,grid_u,iLS+1), (grid_u.ny, grid_u.nx))
    grid_v.V .= reshape(veci(vD,grid_v,iLS+1), (grid_v.ny, grid_v.nx))

    i_u_ext, l_u_ext, b_u_ext, r_u_ext, t_u_ext = indices_extension(grid_u, grid_u.LS[iLS], grid_u.ind.inside, periodic_x, periodic_y)
    i_v_ext, l_v_ext, b_v_ext, r_v_ext, t_v_ext = indices_extension(grid_v, grid_v.LS[iLS], grid_v.ind.inside, periodic_x, periodic_y)

    field_extension!(grid_u, grid_u.LS[iLS].u, grid_u.V, i_u_ext, l_u_ext, b_u_ext, r_u_ext, t_u_ext, num.NB, periodic_x, periodic_y)
    field_extension!(grid_v, grid_v.LS[iLS].u, grid_v.V, i_v_ext, l_v_ext, b_v_ext, r_v_ext, t_v_ext, num.NB, periodic_x, periodic_y)
end


"""
Interpolate velocity on scalar grid for regular grids for vizualisation
"""
function interpolate_grid_liquid(gp,gu,gv,u,v)
    
    # us = p .*0
    # vs = p .*0

    us=zeros(gp)
    vs=zeros(gp)

    LS_u =gu.LS[1]
    LS_v = gv.LS[1]
    us .= (
        (u[:,2:end].^2.0 .* LS_u.geoL.dcap[:,2:end,6] .+ 
        u[:,1:end-1].^2.0 .* LS_u.geoL.dcap[:,1:end-1,6]) ./ 
        (LS_u.geoL.dcap[:,1:end-1,6] .+ LS_u.geoL.dcap[:,2:end,6] .+ 1e-8 )
    )
    vs .= (
        (v[2:end,:].^2.0 .* LS_v.geoL.dcap[2:end,:,7] .+ 
        v[1:end-1,:].^2.0 .* LS_v.geoL.dcap[1:end-1,:,7]) ./
        (LS_v.geoL.dcap[1:end-1,:,7] .+ LS_v.geoL.dcap[2:end,:,7] .+ 1e-8 )
    )

    # phL.p .+= (
    #     (phS.u[:,2:end].^2.0 .* LS_u.geoS.dcap[:,2:end,6] .+ 
    #     phS.u[:,1:end-1].^2.0 .* LS_u.geoS.dcap[:,1:end-1,6]) ./ 
    #     (LS_u.geoS.dcap[:,1:end-1,6] .+ LS_u.geoS.dcap[:,2:end,6] .+ 1e-8 )
    # )
    # phL.p .+= (
    #     (phS.v[2:end,:].^2.0 .* LS_v.geoS.dcap[2:end,:,7] .+ 
    #     phS.v[1:end-1,:].^2.0 .* LS_v.geoS.dcap[1:end-1,:,7]) ./
    #     (LS_v.geoS.dcap[1:end-1,:,7] .+ LS_v.geoS.dcap[2:end,:,7] .+ 1e-8 )
    # )


    # for j = 1:gp.ny
    # for i = 1:gp.nx
    #     us[:,j,i]=(u[:,j,i]+u[:,j,i+1])/2
    #     vs[:,j,i]=(v[:,j,i]+v[:,j+1,i])/2
    # end
    # end
    
    return us,vs
end


"""
Interpolate velocity on scalar grid for regular grids for vizualisation
"""
function interpolate_regular_grid(grid,fwdL)
    
    us = fwdL.p .*0
    vs = fwdL.p .*0

    for j = 1:grid.ny
    for i = 1:grid.nx
        us[:,j,i]=(fwdL.u[:,j,i]+fwdL.u[:,j,i+1])/2
        vs[:,j,i]=(fwdL.v[:,j,i]+fwdL.v[:,j+1,i])/2
    end
    end
    
    return us,vs
end


"""
Interpolate velocity on scalar grid for regular grids for vizualisation
"""
function interpolate_regular_grid(grid,fwdL,u,v)
    
    us = fwdL.p .*0
    vs = fwdL.p .*0

    for j = 1:grid.ny
    for i = 1:grid.nx
        us[:,j,i]=(u[:,j,i]+u[:,j,i+1])/2
        vs[:,j,i]=(v[:,j,i]+v[:,j+1,i])/2
    end
    end
    
    return us,vs
end

# """
# From kinetic_energy
# """
# function scal_magnitude(ph, gp, gu, gv)
#     ph.p = zeros(gp)

#     ph.p .= (
#         (fwdL.u[i,:,2:end].^2.0 .* gu.geoL.dcap[:,2:end,6] .+ 
#         fwdL.u[i,:,1:end-1].^2.0 .* gu.geoL.dcap[:,1:end-1,6]) ./ 
#         (gu.geoL.dcap[:,1:end-1,6] .+ gu.geoL.dcap[:,2:end,6] .+ 1e-8)
#     )
#     ph.p .+= (
#         (fwdL.v[i,2:end,:].^2.0 .* gv.geoL.dcap[2:end,:,7] .+ 
#         fwdL.v[i,1:end-1,:].^2.0 .* gv.geoL.dcap[1:end-1,:,7]) ./
#         (gv.geoL.dcap[1:end-1,:,7] .+ gv.geoL.dcap[2:end,:,7] .+ 1e-8)
#     )
#     ph.p .= sqrt.(ph.p .* gp.geoL.dcap[:,:,5])

# end



"""
From plot_grid
"""
function plot_grid_fig!(fig,ax,
    num, grid;
    linewidth = 0.5,
    limitsx = false,
    limitsy = false,
    hide = false,
    skipx = 0,
    skipy = 0,
    xscale =1.0,
    yscale = 1.0,
    )

    x = grid.x_nodes ./xscale
    y = grid.y_nodes ./yscale
  
    for i = 1:skipx+1:grid.nx+1
        lines!(ones(grid.ny+1) .* x[i], y, linewidth = linewidth, color = :black)
    end
    for i = 1:skipy+1:grid.ny+1
        lines!(x, ones(grid.nx+1) .* y[i], linewidth = linewidth, color = :black)
    end

end


"""
    make_video(
        grid, field_u, field = nothing;
        title_prefix = "video",
        title_suffix = "",
        xlabel = L"x",
        ylabel = L"y",
        colormap = :viridis,
        minv = 0.0,
        maxv = 0.0,
        limitsx = false,
        limitsy = false,
        var = 1,
        framerate = 24,
        step = 1,
        step0 = 1,
        stepf = size(field_u, 1)
    )

Generates and saves a video displaying a `field` and an interface `field_u`.
"""
function make_video_vec(
    num, grid, field_u, field = nothing,
     vecu = nothing, vecv=nothing;
    title_prefix = "video",
    title_suffix = "",
    xlabel = L"x",
    ylabel = L"y",
    colormap = :viridis,
    sz = (1600, 1000),
    minv = 0.0,
    maxv = 0.0,
    limitsx = false,
    limitsy = false,
    var = 1,
    framerate = 24,
    step = 1,
    step0 = 1,
    stepf = size(field_u, 2),
    xscale = 1, 
    yscale = 1, 
    xticks = nothing,
    yticks = nothing,
    scalscale = 1, 
    scalelabel = nothing,
    scalticks = nothing,
    )

    x = grid.x[1,:] ./xscale
    y = grid.y[:,1] ./yscale

    u = field_u[:,step0:stepf,:,:]
    plot_hmap = true
    if isnothing(field)
        plot_hmap = false
    else
        if length(size(field)) == 2
            if var == 1
                z = reshape(
                    field[step0:stepf, 1:grid.ny*grid.nx], 
                    (stepf-step0+1, grid.ny,grid.nx)
                )
            else
                z = reshape(
                    field[step0:stepf, grid.ny*grid.nx+1:2*grid.ny*grid.nx],
                    (stepf-step0+1, grid.ny,grid.nx)
                )
            end
        else
            z = field[step0:stepf,:,:]
        end

        z = z ./scalscale
    end

    if minv == maxv == 0.0
        var_colorrange = true
    else
        var_colorrange = false
    end

    if isa(limitsx, Tuple{Float64, Float64}) || isa(limitsx, Tuple{Int, Int})
        lx = limitsx
    else
        lx = (min(x...) - grid.dx[1,1] / 2, max(x...) + grid.dx[1,end] / 2)
    end
    if isa(limitsy, Tuple{Float64, Float64}) || isa(limitsy, Tuple{Int, Int})
        ly = limitsy
    else
        ly = (min(y...) - grid.dy[1,1] / 2, max(y...) + grid.dy[end,1] / 2)
    end

    obs = Observable{Int32}(1)
    iterator = range(0, size(u, 2) - 1, step=step)

    # fig = Figure(size = sz)
    fig = Figure()
    ax  = Axis(fig[1,1], 
    # aspect=DataAspect(), 
    aspect=1,
    xlabel = xlabel, ylabel = ylabel,
    xticks=xticks, yticks=yticks,
        xtickalign = 0,  ytickalign = 0)
    if plot_hmap
        print("scalticks", scalticks, "xticks", xticks, "yticks", yticks,"\n")

        if !var_colorrange
            hmap = heatmap!(ax, x, y, @lift(z[$obs,:,:]'), colormap = colormap, 
                colorrange = (minv, maxv))
        else
            # if !isnothing(scalticks)
            #     print("scalticks", scalticks)

            #     hmap = heatmap!(x, y, @lift(z[$obs,:,:]'), colormap = colormap, ticks=scalticks)
            # else
            hmap = heatmap!(ax, x, y, @lift(z[$obs,:,:]'), colormap = colormap)
            # end
        end
    end
    if !plot_hmap
        contour!(x, y, @lift(u[1,$obs,:,:]'), levels = -10:grid.dx[1,1]:10, linewidth = 2.0)
    end
    for iLS in 1:num.nLS
        contour!(x, y, @lift(u[iLS,$obs,:,:]'), levels = [0.0], color = :red, linewidth = 3.0)
    end

    #TODO
    if !isnothing(vecu)
        arrows!(grid.x[1,:], grid.y[:,1], vecu[1,:,:], vecv[1,:,:], 
        arrowsize = 10, 
        lengthscale = 0.3,
        # arrowcolor = strength, 
        # linecolor = strength,
        )
    end

    # contourf!(x, y, u[2,1,:,:]', levels = 0:0, color = :red, linewidth = 3, extendlow = :gray);
    if plot_hmap
        if !isnothing(scalticks)
            print("scalticks", scalticks)
            cbar = fig[1,2] = Colorbar(fig, hmap, labelpadding = 0, ticks=scalticks)
            
            # Colorbar(fig[1, 2], hmap, 
            # label = "Reverse sequential colormap",
            # ticks=scalticks)
            # cbar.ticks = (scalticks)
            #https://stackoverflow.com/questions/73555826/modifying-label-and-tick-on-color-bar-of-plots-jl-plot

            Colorbar(fig[1, 2], hmap, ticks = -1:0.25:1)

        else
            cbar = fig[1,2] = Colorbar(fig, hmap, labelpadding = 0)
        end

    end

    if !isnothing(scalelabel)
        # Label(fig[1, 1, Top()], halign = :right, scalelabel)
        Label(fig[1, 2, Top()], halign = :center, scalelabel)
    end

    # limits!(ax, lx[1], lx[2], ly[1], ly[2])
    # colgap!(fig.layout, 10)
    # rowgap!(fig.layout, 10)
    # colsize!(fig.layout, 1, widths(ax.scene.viewport[])[1])
    # rowsize!(fig.layout, 1, widths(ax.scene.viewport[])[2])
    # resize_to_layout!(fig)

    vid = record(
            fig, title_prefix*title_suffix*".mp4", iterator;
            framerate = framerate
        ) do it
        obs[] = it+1
    end

    return vid
end

function print_electrolysis_statistics(nb_transported_scalars,phL)

    # normscal1L = norm(phL.trans_scal[:,:,1])
    # normscal2L = norm(phL.trans_scal[:,:,2])
    # normscal3L = norm(phL.trans_scal[:,:,3])
    # normphi_eleL = norm(phL.phi_ele)
    # normTL = norm(phL.T)
    # normiL=0.0 #TODO

    

    minscal1L=minscal2L=minscal3L=minphi_eleL=minTL=miniL=0.0
    maxscal1L=maxscal2L=maxscal3L=maxphi_eleL=maxTL=maxiL=0.0
    moyscal1L=moyscal2L=moyscal3L=moyphi_eleL=moyTL=moyiL=0.0



    minscal1L = minimum(phL.trans_scal[:,:,1])
    maxscal1L = maximum(phL.trans_scal[:,:,1])
    moyscal1L = mean(phL.trans_scal[:,:,1])

    if nb_transported_scalars>1
        minscal2L = minimum(phL.trans_scal[:,:,2])
        maxscal2L = maximum(phL.trans_scal[:,:,2])
        moyscal2L = mean(phL.trans_scal[:,:,2])
    
        if nb_transported_scalars>2
            minscal3L = minimum(phL.trans_scal[:,:,3])
            maxscal3L = maximum(phL.trans_scal[:,:,3])
            moyscal3L = mean(phL.trans_scal[:,:,3])
        end

    end

    # if heat
    #     minTL = minimum(phL.T)
    #     maxTL = maximum(phL.T)
    #     moyTL = mean(phL.T)
    # end

    minTL = minimum(phL.T)
    maxTL = maximum(phL.T)
    moyTL = mean(phL.T)

    minphi_eleL = minimum(phL.phi_ele)
    miniL=minimum(phL.i_current_mag)

    maxphi_eleL = maximum(phL.phi_ele)
    maxiL=maximum(phL.i_current_mag)

    moyphi_eleL = mean(phL.phi_ele)
    moyiL=mean(phL.i_current_mag)

    # print("$(@sprintf("norm(cH2) %.6e", normscal1L))\t$(@sprintf("norm(KOH) %.6e", normscal2L))\t$(@sprintf("norm(H2O) %.6e", normscal3L))\n")
    # print("$(@sprintf("norm(phi_ele) %.6e", normphi_eleL))\t$(@sprintf("norm(T) %.6e", normTL))\t$(@sprintf("norm(i) %.6e", normiL))\n")

    print("$(@sprintf("min(cH2) %.6e", minscal1L))\t$(@sprintf("min(KOH) %.6e", minscal2L))\t$(@sprintf("min(H2O) %.6e", minscal3L))\n")
    print("$(@sprintf("max(cH2) %.6e", maxscal1L))\t$(@sprintf("max(KOH) %.6e", maxscal2L))\t$(@sprintf("max(H2O) %.6e", maxscal3L))\n")
    print("$(@sprintf("moy(cH2) %.6e", moyscal1L))\t$(@sprintf("moy(KOH) %.6e", moyscal2L))\t$(@sprintf("moy(H2O) %.6e", moyscal3L))\n")

    print("$(@sprintf("min(phi_ele) %.6e", minphi_eleL))\t$(@sprintf("min(T) %.6e", minTL))\t$(@sprintf("min(i) %.6e", miniL))\n")
    print("$(@sprintf("max(phi_ele) %.6e", maxphi_eleL))\t$(@sprintf("max(T) %.6e", maxTL))\t$(@sprintf("max(i) %.6e", maxiL))\n")
    print("$(@sprintf("moy(phi_ele) %.6e", moyphi_eleL))\t$(@sprintf("moy(T) %.6e", moyTL))\t$(@sprintf("moy(i) %.6e", moyiL))\n")

end

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
  Compute norm of gradient for exchange current
"""
# Gradient of pressure, eq. 17 in 
#"A Conservative Cartesian Cut-Cell Method for Mixed Boundary Conditions and the Incompressible Navier-Stokes Equations on Staggered Meshes"
#From navier_stokes_coupled.jl
# ∇ϕ_x = opC_u.AxT * opC_u.Rx * vec(phi_ele) .+ opC_u.Gx_b * vecb(phi_eleD,grid)
# ∇ϕ_y = opC_v.AyT * opC_v.Ry * vec(phi_ele) .+ opC_v.Gy_b * vecb(phi_eleD,grid)
# for iLS in 1:nLS
#     ∇ϕ_x .+= opC_u.Gx[iLS] * veci(phi_eleD,grid,iLS+1)
#     ∇ϕ_y .+= opC_v.Gy[iLS] * veci(phi_eleD,grid,iLS+1)
# end
function compute_grad_phi_ele!(num,grid, grid_u, grid_v, phL, phS,  opC_pL, opC_pS)
    
    @unpack nLS = num

    LS_u =grid_u.LS[1]
    LS_v = grid_v.LS[1]
    # LS =gp.LS[1]

    eps_den = 1e-8
    #TODO different way to do it?

    #Liquid phase
    @unpack phi_eleD = phL
    opC_p = opC_pL

    ∇ϕ_x = opC_p.iMx * opC_p.Bx * vec1(phi_eleD,grid) .+ opC_p.iMx_b * opC_p.Hx_b * vecb(phi_eleD,grid)
    ∇ϕ_y = opC_p.iMy * opC_p.By * vec1(phi_eleD,grid) .+ opC_p.iMy_b * opC_p.Hy_b * vecb(phi_eleD,grid)

    for iLS in 1:nLS
        ∇ϕ_x .+= opC_p.iMx * opC_p.Hx[iLS] * veci(phi_eleD,grid,iLS+1)
        ∇ϕ_y .+= opC_p.iMy * opC_p.Hy[iLS] * veci(phi_eleD,grid,iLS+1)
    end

    grd_x = reshape(veci(∇ϕ_x,grid_u,1), grid_u)
    grd_y = reshape(veci(∇ϕ_y,grid_v,1), grid_v)

    phL.i_current_mag .= (
        (grd_x[:,2:end].^2.0 .* LS_u.geoL.dcap[:,2:end,6] .+ 
        grd_x[:,1:end-1].^2.0 .* LS_u.geoL.dcap[:,1:end-1,6]) ./ 
        (LS_u.geoL.dcap[:,1:end-1,6] .+ LS_u.geoL.dcap[:,2:end,6] .+ eps_den )
    )
    phL.i_current_mag .+= (
        (grd_y[2:end,:].^2.0 .* LS_v.geoL.dcap[2:end,:,7] .+ 
        grd_y[1:end-1,:].^2.0 .* LS_v.geoL.dcap[1:end-1,:,7]) ./
        (LS_v.geoL.dcap[1:end-1,:,7] .+ LS_v.geoL.dcap[2:end,:,7] .+ eps_den )
    )

    phL.Eu .= grd_x
    phL.Ev .= grd_y

    #Solid phase
    @unpack phi_eleD = phS

    opC_p = opC_pS

    ∇ϕ_x = opC_p.iMx * opC_p.Bx * vec1(phi_eleD,grid) .+ opC_p.iMx_b * opC_p.Hx_b * vecb(phi_eleD,grid)
    ∇ϕ_y = opC_p.iMy * opC_p.By * vec1(phi_eleD,grid) .+ opC_p.iMy_b * opC_p.Hy_b * vecb(phi_eleD,grid)

    for iLS in 1:nLS
        ∇ϕ_x .+= opC_p.iMx * opC_p.Hx[iLS] * veci(phi_eleD,grid,iLS+1)
        ∇ϕ_y .+= opC_p.iMy * opC_p.Hy[iLS] * veci(phi_eleD,grid,iLS+1)
    end

    grd_x = reshape(veci(∇ϕ_x,grid_u,1), grid_u)
    grd_y = reshape(veci(∇ϕ_y,grid_v,1), grid_v)

    # phS.i_current_mag .= (
    #     (grd_x[:,2:end].^2.0 .* LS_u.geoS.dcap[:,2:end,6] .+ 
    #     grd_x[:,1:end-1].^2.0 .* LS_u.geoS.dcap[:,1:end-1,6]) ./ 
    #     (LS_u.geoS.dcap[:,1:end-1,6] .+ LS_u.geoS.dcap[:,2:end,6] + eps_den )
    # )
    # phS.i_current_mag .+= (
    #     (grd_y[2:end,:].^2.0 .* LS_v.geoS.dcap[2:end,:,7] .+ 
    #     ph.v[1:end-1,:].^2.0 .* LS_v.geoS.dcap[1:end-1,:,7]) ./
    #     (LS_v.geoS.dcap[1:end-1,:,7] .+ LS_v.geoS.dcap[2:end,:,7] + eps_den )
    # )

    phL.i_current_mag .+= (
        (grd_x[:,2:end].^2.0 .* LS_u.geoS.dcap[:,2:end,6] .+ 
        grd_x[:,1:end-1].^2.0 .* LS_u.geoS.dcap[:,1:end-1,6]) ./ 
        (LS_u.geoS.dcap[:,1:end-1,6] .+ LS_u.geoS.dcap[:,2:end,6] .+ eps_den )
    )
    phL.i_current_mag .+= (
        (grd_y[2:end,:].^2.0 .* LS_v.geoS.dcap[2:end,:,7] .+ 
        grd_y[1:end-1,:].^2.0 .* LS_v.geoS.dcap[1:end-1,:,7]) ./
        (LS_v.geoS.dcap[1:end-1,:,7] .+ LS_v.geoS.dcap[2:end,:,7] .+ eps_den )
    )

    phL.i_current_mag .= sqrt.(phL.i_current_mag)
    # phS.i_current_mag .= sqrt.(phS.i_current_mag)

    phS.Eu .= grd_x
    phS.Ev .= grd_y


end

"""
Adapt timestep based on velocity, viscosity, ...
# [Kang 2000, “A Boundary Condition Capturing Method for Multiphase Incompressible Flow”](https://doi.org/10.1023/A:1011178417620) 
"""
function adapt_timestep!(num, phL, phS, grid_u, grid_v,adapt_timestep_mode)
    @unpack rho1,rho2,mu1,mu2,grav_x,grav_y,CFL=num

    #With grid_u
    min_spacing_x = minimum(grid_u.dx)
    min_spacing_y = minimum(grid_u.dy)

    # min_spacing_xyz = min(min_spacing_x,min_spacing_y)

    # τ = min(CFL*Δ^2*Re, CFL*Δ/max(abs.(V)..., abs.(phL.u)..., abs.(phL.v)..., abs.(phS.u)..., abs.(phS.v)...))

    #TODO V ?
    # vel = max(abs.(V)..., abs.(phL.u)..., abs.(phL.v)..., abs.(phS.u)..., abs.(phS.v)...)

    if adapt_timestep_mode==1
        c_conv = max(abs.(phL.u)./grid_u.dx..., abs.(phL.v)./grid_v.dy..., abs.(phS.u)./grid_u.dx..., abs.(phS.v)./grid_v.dy...)

        c_visc = max(mu1/rho1, mu2/rho2) * (2.0/min_spacing_x^2 + 2.0/min_spacing_y^2 ) 
        
        c_grav = sqrt(max(abs(grav_x)/min_spacing_x, abs(grav_y)/min_spacing_y))

        #TODO capillary timestep
        c_surf = 0.0

        # c_surf = sigma*kappa / (min(rho1,rho2)*min_spacing_xyz^2)

        c_tot = (c_conv+c_visc)+ sqrt((c_conv+c_visc)^2+4*c_grav^2+4*c_surf^2 )

        return 2*CFL/c_tot
    elseif adapt_timestep_mode==2
        c_conv = max(abs.(phL.u)./grid_u.dx..., abs.(phL.v)./grid_v.dy..., abs.(phS.u)./grid_u.dx..., abs.(phS.v)./grid_v.dy...)

        return CFL/c_conv
    end

    # printstyled(color=:green, @sprintf "\n rho1 %.2e rho2 %.2e mu1 %.2e mu2 %.2e\n" rho1 rho2 mu1 mu2) 
    # printstyled(color=:green, @sprintf "\n dx %.2e dy %.2e \n" min_spacing_x min_spacing_x) 

    
    # max(mu1/rho1, mu2/rho2) * (2.0/min_spacing_x^2 + 2.0/min_spacing_y^2 )

    #  printstyled(color=:green, @sprintf "\n c_conv %.2e c_visc %.2e c_grad %.2e\n" c_conv c_visc c_grav)
    #  printstyled(color=:green, @sprintf "\n CFL : %.2e dt : %.2e\n" CFL num.τ)
end

function init_fields_2(TD,T,H,BC,grid,dir_val_intfc)

    vec1(TD,grid) .= vec(T)
    vec2(TD,grid) .= dir_val_intfc

    if is_neumann(BC.left)
        vecb_L(TD,grid) .= T[:,1] .+ H[:,1] .* BC.left.val
    else
        vecb_L(TD,grid) .= BC.left.val #.* ones(grid.ny)
    end
    if is_neumann(BC.bottom)
        vecb_B(TD,grid) .= T[1,:] .+ H[1,:] .* BC.bottom.val
    else
        vecb_B(TD,grid) .= BC.bottom.val .* ones(grid.nx)
    end
    if is_neumann(BC.right)
        vecb_R(TD,grid) .= T[:,end] .+ H[:,end] .* BC.right.val 
    else
        vecb_R(TD,grid) .= BC.right.val .* ones(grid.ny)
    end
    if is_neumann(BC.top)
        vecb_T(TD,grid) .= T[end,:] .+ Hu[end,:] .* BC.top.val
    else
        vecb_T(TD,grid) .= BC.top.val .* ones(grid.nx)
    end
    
end