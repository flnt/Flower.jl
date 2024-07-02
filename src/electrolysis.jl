#From set_heat in heat_coupled.jl, poisson.jl

function Poiseuille_fmax(x,v_inlet_max,L0)
    return 4*v_inlet_max*x/L0*(1-x/L0)
end

function Poiseuille_favg(x,v_inlet_moy,L0)
    return 6*v_inlet_moy*x/L0*(1-x/L0)
end



function set_convection_2!(
    num, grid, geo, grid_u, LS_u, grid_v, LS_v,
    u, v, op, ph, BC_u, BC_v
    )
    @unpack Cu, CUTCu, Cv, CUTCv = op
    @unpack uD, vD = ph

    Du_x = zeros(grid_u)
    Du_y = zeros(grid_u)
    for iLS in 1:num.nLS
        Du_x[LS_u[iLS].MIXED] .= reshape(veci(uD,grid_u,iLS+1), grid_u)[LS_u[iLS].MIXED]
        Du_y[LS_u[iLS].MIXED] .= reshape(veci(uD,grid_u,iLS+1), grid_u)[LS_u[iLS].MIXED]
    end
    Du_x[:,1] .= vecb_L(uD,grid_u)
    Du_y[1,:] .= vecb_B(uD,grid_u)
    Du_x[:,end] .= vecb_R(uD,grid_u)
    Du_y[end,:] .= vecb_T(uD,grid_u)

    # bnds_u = [grid_u.ind.b_left[1], grid_u.ind.b_bottom[1], grid_u.ind.b_right[1], grid_u.ind.b_top[1]]
    # bnds_v = [grid_v.ind.b_left[1], grid_v.ind.b_bottom[1], grid_v.ind.b_right[1], grid_v.ind.b_top[1]]
    # Δu = [grid_u.dx[1,1] * 0.25, grid_u.dy[1,1] * 0.5, grid_u.dx[end,end] * 0.25, grid_u.dy[end,end] * 0.5]
    # Δv = [grid_v.dx[1,1] * 0.5, grid_v.dy[1,1] * 0.25, grid_v.dx[end,end] * 0.5, grid_v.dy[end,end] * 0.25]

    # Hu = zeros(grid_u)
    # for i in eachindex(bnds_u)
    #     for II in bnds_u[i]
    #         Hu[II] = Δu[i]
    #     end
    # end

    # Hv = zeros(grid_v)
    # for i in eachindex(bnds_v)
    #     for II in bnds_v[i]
    #         Hv[II] = Δv[i]
    #     end
    # end


    vector_convection!(dir, GridFCx, Cu, CUTCu, u, v, Du_x, Du_y, Dv_x, Dv_y,
            geo.dcap, grid.nx, grid.ny, BC_u, grid_u.ind.inside,
            grid_u.ind.b_left[1], grid_u.ind.b_bottom[1], grid_u.ind.b_right[1], grid_u.ind.b_top[1])
    
    return nothing
end

"""    
set convection and diffusion for scalar
"""
function set_scalar_transport_2!(bc_type, num, grid, op, geo, ph, θd, BC_T, MIXED, projection,
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

    ######################################################################################################
    #Interface boundary condition
    ######################################################################################################
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

    ######################################################################################################
    #Wall BC
    ######################################################################################################
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
        # @inbounds @threads for II in vcat(b_left[1], b_bottom[1], b_right[1], b_top[1])
        @inbounds for II in vcat(b_left[1], b_bottom[1], b_right[1], b_top[1])
            HT[II] = distance(grid.LS[1].mid_point[II], geo.centroid[II], dx[II], dy[II])
        end    
        bcTx, bcTy = set_bc_bnds(dir, a0, HT, BC_T)


        # set_convection_2!(num, grid, geo[end], grid_u, grid_u.LS, grid_v, grid_v.LS, ph.u, ph.v, op_conv, ph, BC_u, BC_v)

    
        # bnds_u = [grid_u.ind.b_left[1], grid_u.ind.b_bottom[1], grid_u.ind.b_right[1], grid_u.ind.b_top[1]]
        # bnds_v = [grid_v.ind.b_left[1], grid_v.ind.b_bottom[1], grid_v.ind.b_right[1], grid_v.ind.b_top[1]]
        # Δu = [grid_u.dx[1,1], grid_u.dy[1,1], grid_u.dx[end,end], grid_u.dy[end,end]] .* 0.5
        # Δv = [grid_v.dx[1,1], grid_v.dy[1,1], grid_v.dx[end,end], grid_v.dy[end,end]] .* 0.5
        
        # Hu = zeros(grid_u)
        # for i in eachindex(bnds_u)
        #     for II in bnds_u[i]
        #         Hu[II] = Δu[i]
        #     end
        # end

        # Hv = zeros(grid_v)
        # for i in eachindex(bnds_v)
        #     for II in bnds_v[i]
        #         Hv[II] = Δv[i]
        #     end
        # end
    
        bcU = zeros(grid_u)
        bcU .= reshape(vec2(uD,grid_u), grid_u)
        bcU[1,:] .= vecb_B(uD, grid_u)
        bcU[end,:] .= vecb_T(uD, grid_u)
        bcU[:,1] .= vecb_L(uD, grid_u)
        bcU[:,end] .= vecb_R(uD, grid_u)
        
        bcV = zeros(grid_v)
        bcV .= reshape(vec2(vD,grid_v), grid_v)
        bcV[:,1] .= vecb_L(vD, grid_v)
        bcV[:,end] .= vecb_R(vD, grid_v)
        bcV[1,:] .= vecb_B(vD, grid_v)
        bcV[end,:] .= vecb_T(vD, grid_v)

        scalar_convection!(dir, CT, CUTCT, u, v, bcTx, bcTy, bcU, bcV, geo.dcap, ny, 
            BC_T, inside, b_left[1], b_bottom[1], b_right[1], b_top[1]
        )

 
        # vector_convection!(dir, GridFCx, Cu, CUTCu, u, v, Du_x, Du_y, Dv_x, Dv_y,
        # geo.dcap, grid.nx, grid.ny, BC_u, grid_u.ind.inside,
        # grid_u.ind.b_left[1], grid_u.ind.b_bottom[1], grid_u.ind.b_right[1], grid_u.ind.b_top[1])
        
        # vector_convection!(dir, GridFCy, Cv, CUTCv, u, v, Du_x, Du_y, Dv_x, Dv_y,
        # geo.dcap, grid.nx, grid.ny, BC_v, grid_v.ind.inside,
        # grid_v.ind.b_left[1], grid_v.ind.b_bottom[1], grid_v.ind.b_right[1], grid_v.ind.b_top[1])



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
    B[1:ni,1:ni] = M .+ 0.5 .* τ .* diffusion_coeff_scal .* LT .- τ .* CT
    B[1:ni,ni+1:2*ni] = 0.5 .* τ .* diffusion_coeff_scal .* LD
    B[1:ni,end-nb+1:end] = 0.5 .* τ .* diffusion_coeff_scal .* LD_b

    rhs = fnzeros(grid, num)
    if convection
        vec1(rhs,grid) .-= τ .* CUTCT
    end

    #IItest = CartesianIndex(66, 4) #(id_y, id_x)
    # IItest = CartesianIndex(69, 3) #(id_y, id_x)
    # pIItest = lexicographic(IItest, grid.ny)

    # print("\n rhs",vec2(rhs,grid)[pIItest])

    vec2(rhs,grid) .+= χ[1] * vec(a0)

    # print("\n rhs",vec2(rhs,grid)[pIItest])

    vecb(rhs,grid) .+= op.χ_b * vec(a0_b)
    
    # print("\n rhs",vec2(rhs,grid)[pIItest])
    # print("\n rhs",bc_type.val*grid.dx[1,1])


    return rhs
end


    function set_bc_bnds_total(::Dirichlet, D, H, BC)
        Dx = copy(D)
        Dy = copy(D)

        if is_neumann(BC.left)
            @inbounds Dx[:,1] .= H[:,1] .* BC.left.val
        elseif is_dirichlet(BC.left)
            @inbounds Dx[:,1] .= BC.left.val
        end
        @inbounds Dy[:,1] .= Dx[:,1]

        if is_neumann(BC.bottom)
            @inbounds Dy[1,:] .= H[1,:] .* BC.bottom.val 
        elseif is_dirichlet(BC.bottom)
            @inbounds Dy[1,:] .= BC.bottom.val
        end
        @inbounds Dx[1,:] .= Dy[1,:]

        if is_neumann(BC.right)
            @inbounds Dx[:,end] .= H[:,end] .* BC.right.val 
        elseif is_dirichlet(BC.right)
            @inbounds Dx[:,end] .= BC.right.val
        end
        @inbounds Dy[:,end] .= Dx[:,end]

        if is_neumann(BC.top)
            @inbounds Dy[end,:] .= H[end,:] .* BC.top.val 
        elseif is_dirichlet(BC.top)
            @inbounds Dy[end,:] .= BC.top.val
        end
        @inbounds Dx[end,:] .= Dy[end,:]


        return Dx, Dy
    end



"""    
scalar transport (convection and diffusion)
"""
function scalar_transport_2!(bc, num, grid, op, geo, ph, concentration0, MIXED, projection,
    op_conv, grid_u, geo_u, grid_v, geo_v,
    periodic_x, periodic_y, convection, ls_advection, BC_int, diffusion_coeff, convection_Cdivu)
    @unpack τ, aniso, nb_transported_scalars,nb_saved_scalars = num
    @unpack nx, ny, dx, dy, ind, LS  = grid
    @unpack all_indices, inside, b_left, b_bottom, b_right, b_top = ind
    @unpack Bx, By, BxT, ByT, Hx, Hy, HxT, HyT, M, iMx, iMy, χ = op
    @unpack CT, CUTCT = op_conv
    @unpack u, v, uD, vD = ph

    printstyled(color=:red, @sprintf "\n levelset: start scalar_transport!\n")
    println(grid.LS[1].geoL.dcap[1,1,:])

    ni = nx * ny
    nb = 2 * nx + 2 * ny
    nt = 2 * ni + nb

    A = spzeros(nt, nt)
    B = spzeros(nt, nt)

    all_CUTCT = zeros(grid.ny * grid.nx, nb_transported_scalars)

    ######################################################################################################
    # Operators
    ######################################################################################################


    if ls_advection
        # update_all_ls_data(num, grid, grid_u, grid_v, BC_int, periodic_x, periodic_y, false) #already false

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
    ######################################################################################################

    update_all_ls_data(num, grid, grid_u, grid_v, BC_int, periodic_x, periodic_y)


    if convection

        bcU = zeros(grid_u)
        bcU .= reshape(vec2(uD,grid_u), grid_u)
        bcU[1,:] .= vecb_B(uD, grid_u)
        bcU[end,:] .= vecb_T(uD, grid_u)
        bcU[:,1] .= vecb_L(uD, grid_u)
        bcU[:,end] .= vecb_R(uD, grid_u)
        
        bcV = zeros(grid_v)
        bcV .= reshape(vec2(vD,grid_v), grid_v)
        bcV[:,1] .= vecb_L(vD, grid_v)
        bcV[:,end] .= vecb_R(vD, grid_v)
        bcV[1,:] .= vecb_B(vD, grid_v)
        bcV[end,:] .= vecb_T(vD, grid_v)

        # print("\n vecb_L ", vecb_L(vD, grid_v))
        # print("\n vecb_B ", vecb_B(vD, grid_v))

        HT = zeros(grid)
        # @inbounds @threads for II in vcat(b_left[1], b_bottom[1], b_right[1], b_top[1])
        @inbounds for II in vcat(b_left[1], b_bottom[1], b_right[1], b_top[1])
            HT[II] = distance(grid.LS[1].mid_point[II], geo.centroid[II], dx[II], dy[II])
        end  

        iscal = 1

        # printstyled(color=:red, @sprintf "\n levelset: before scalar_convection\n")
        # println(grid.LS[1].geoL.dcap[1,1,:])

        ######################################################################################################
        #Interface boundary condition
        ######################################################################################################
        if is_dirichlet(bc[iscal].int)
            __a0 = bc[iscal].int.val        
        elseif is_neumann(bc[iscal].int)
            __a0 = bc[iscal].int.val
        elseif is_robin(bc[iscal].int)
            __a0 = bc_type.val  
        elseif is_stefan(bc[iscal].int)
            __a0 = concentration0[iscal]       
        elseif is_wall(bc[iscal].int)
            __a0 = bc[iscal].int.val         
        else
            __a0 = bc[iscal].int.val
        end

        # Flags with BCs
        a0 = ones(grid) .* __a0

        # The idea there is to have at the corners the contributions from both borders. 
        # bcTx is only used for derivatives in x and bcTy for derivatives in y, so it doesn't matter 
        # that bcTx is not filled at the top and bottom and the same for bcTy
        # Otherwise we would have to select one of the contributions
        bcTx, bcTy = set_bc_bnds(dir, a0, HT, bc[iscal])

        # @views necessary
        @views scalar_convection!(dir, CT, all_CUTCT[:,iscal], u, v, bcTx, bcTy, bcU, bcV, geo.dcap, ny, 
        bc[iscal], inside, b_left[1], b_bottom[1], b_right[1], b_top[1]
        )

        if nb_transported_scalars>1
            for iscal=2:nb_transported_scalars

                # printstyled(color=:red, @sprintf "\n levelset: before scalar_convection\n")
                # println(grid.LS[1].geoL.dcap[1,1,:])

                ######################################################################################################
                #Interface boundary condition
                ######################################################################################################
                if is_dirichlet(bc[iscal].int)
                    __a0 = bc[iscal].int.val        
                elseif is_neumann(bc[iscal].int)
                    __a0 = bc[iscal].int.val
                elseif is_robin(bc[iscal].int)
                    __a0 = bc_type.val  
                elseif is_stefan(bc[iscal].int)
                    __a0 = concentration0[iscal]       
                elseif is_wall(bc[iscal].int)
                    __a0 = bc[iscal].int.val         
                else
                    __a0 = bc[iscal].int.val
                end

                # Flags with BCs
                a0 = ones(grid) .* __a0

                bcTx, bcTy = set_bc_bnds(dir, a0, HT, bc[iscal])
                # @views necessary
                @views scalar_convection_CUTCT!(dir, all_CUTCT[:,iscal], u, v, bcTx, bcTy, bcU, bcV, geo.dcap, ny, 
                bc[iscal], inside, b_left[1], b_bottom[1], b_right[1], b_top[1]
                )
            end
        end       
    end #convection

    # printstyled(color=:red, @sprintf "\n levelset: end scalar_convection\n")
    # println(grid.LS[1].geoL.dcap[1,1,:])


   

    #either update here or compute 
    # update_all_ls_data(num, grid, grid_u, grid_v, BC_int, periodic_x, periodic_y)

    for iscal=1:nb_transported_scalars

        # printstyled(color=:red, @sprintf "\n levelset: start iscal!\n")
        # println(grid.LS[1].geoL.dcap[1,1,:])
    
        diffusion_coeff_scal = diffusion_coeff[iscal]

        # if iscal==3
        #     printstyled(color=:red, @sprintf "\n levelset: forcing diff 0\n")
        #     diffusion_coeff_scal = 0.0
        # end

        ######################################################################################################
        #Interface boundary condition
        ######################################################################################################
        if is_dirichlet(bc[iscal].int)
            # printstyled(color=:green, @sprintf "\n Dirichlet : %.2e\n" bc[iscal].int.val )
            __a0 = bc[iscal].int.val
            __a1 = -1.0
            __b = 0.0
        elseif is_neumann(bc[iscal].int)
            __a0 = bc[iscal].int.val
            __a1 = 0.0
            __b = 1.0
        elseif is_robin(bc[iscal].int)
            __a0 = bc_type.val
            __a1 = -1.0
            __b = 1.0
        elseif is_stefan(bc[iscal].int)
            __a0 = concentration0[iscal]
            __a1 = -1.0
            __b = 0.0
        elseif is_wall(bc[iscal].int)
            __a0 = bc[iscal].int.val
            __a1 = -1.0
            __b = 0.0
        else
            __a0 = bc[iscal].int.val
            __a1 = -1.0
            __b = 0.0
        end

        # Flags with BCs
        a0 = ones(grid) .* __a0

        ######################################################################################################
        #Wall BC
        ######################################################################################################
        _a1 = ones(grid) .* __a1
        a1 = Diagonal(vec(_a1))
        _b = ones(grid) .* __b
        b = Diagonal(vec(_b))

        a0_b = zeros(nb)
        _a1_b = zeros(nb)
        _b_b = zeros(nb)
        set_borders!(grid, grid.LS[1].cl, grid.LS[1].u, a0_b, _a1_b, _b_b, bc[iscal], num.n_ext_cl)
        a1_b = Diagonal(vec(_a1_b))
        b_b = Diagonal(vec(_b_b))

        # if convection
        #     # HT = zeros(grid)
        #     # # @inbounds @threads for II in vcat(b_left[1], b_bottom[1], b_right[1], b_top[1])
        #     # @inbounds for II in vcat(b_left[1], b_bottom[1], b_right[1], b_top[1])
        #     #     HT[II] = distance(grid.LS[1].mid_point[II], geo.centroid[II], dx[II], dy[II])
        #     # end    

        #     # The idea there is to have at the corners the contributions from both borders. 
        #     # bcTx is only used for derivatives in x and bcTy for derivatives in y, so it doesn't matter 
        #     # that bcTx is not filled at the top and bottom and the same for bcTy
        #     # Otherwise we would have to select one of the contributions

        #     bcTx, bcTy = set_bc_bnds(dir, a0, HT, bc[iscal])
        
        #     bcU = zeros(grid_u)
        #     bcU .= reshape(vec2(uD,grid_u), grid_u)
        #     bcU[1,:] .= vecb_B(uD, grid_u)
        #     bcU[end,:] .= vecb_T(uD, grid_u)
        #     bcU[:,1] .= vecb_L(uD, grid_u)
        #     bcU[:,end] .= vecb_R(uD, grid_u)
            
        #     bcV = zeros(grid_v)
        #     bcV .= reshape(vec2(vD,grid_v), grid_v)
        #     bcV[:,1] .= vecb_L(vD, grid_v)
        #     bcV[:,end] .= vecb_R(vD, grid_v)
        #     bcV[1,:] .= vecb_B(vD, grid_v)
        #     bcV[end,:] .= vecb_T(vD, grid_v)

        #     print("\n vecb_L ", vecb_L(vD, grid_v))

        #     print("\n vecb_B ", vecb_B(vD, grid_v))

     

        #     printstyled(color=:red, @sprintf "\n levelset: before scalar_convection\n")
        #     println(grid.LS[1].geoL.dcap[1,1,:])

        #     scalar_convection!(dir, CT, CUTCT, u, v, bcTx, bcTy, bcU, bcV, geo.dcap, ny, 
        #     bc[iscal], inside, b_left[1], b_bottom[1], b_right[1], b_top[1]
        #     )

        # end

        LT = BxT * iMx * Bx .+ ByT * iMy * By
        LD = BxT * iMx * Hx[1] .+ ByT * iMy * Hy[1]
        LD_b = BxT * op.iMx_b * op.Hx_b .+ ByT * op.iMy_b * op.Hy_b

        # A = spzeros(nt, nt)
        # B = spzeros(nt, nt)

        #TODO check no need to reinitialize A and B

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
        B[1:ni,1:ni] = M .+ 0.5 .* τ .* diffusion_coeff_scal .* LT .- τ .* CT
        B[1:ni,ni+1:2*ni] = 0.5 .* τ .* diffusion_coeff_scal .* LD
        B[1:ni,end-nb+1:end] = 0.5 .* τ .* diffusion_coeff_scal .* LD_b

        rhs = fnzeros(grid, num)
        if convection
            vec1(rhs,grid) .-= τ .* all_CUTCT[:,iscal]
        end

        # IItest = CartesianIndex(66, 4) #(id_y, id_x)
        # pIItest = lexicographic(IItest, grid.ny)
        # print("\n rhs",vec2(rhs,grid)[pIItest])

        vec2(rhs,grid) .+= χ[1] * vec(a0)

        # print("\n rhs",vec2(rhs,grid)[pIItest])

        vecb(rhs,grid) .+= op.χ_b * vec(a0_b)
        
        # print("\n rhs",vec2(rhs,grid)[pIItest])
        # print("\n rhs",bc_type.val*grid.dx[1,1])


        if convection_Cdivu
            # Duv = fzeros(grid)
            # Duv = fnzeros(grid,num)

            #TODO BC missing

            Duv = opC_pL.AxT * vec1(ph.uD,grid_u) .+ opC_pL.Gx_b * vecb(ph.uD,grid_u) .+
            opC_pL.AyT * vec1(ph.vD,grid_v) .+ opC_pL.Gy_b * vecb(ph.vD,grid_v)
            for iLS in 1:nLS
                if !is_navier(BC_int[iLS]) && !is_navier_cl(BC_int[iLS])
                    Duv .+= opC_pL.Gx[iLS] * veci(ph.uD,grid_u,iLS+1) .+ 
                            opC_pL.Gy[iLS] * veci(ph.vD,grid_v,iLS+1)
                end
            end
        
            # rhs .+= Duv #.* ph.trans_scalD[:,iscal] multiplied just after

            # vec1(rhs_ϕ,grid) .= rho1 .* iτ .* Duv #TODO
            vec1(rhs,grid) .+= Duv 

            printstyled(color=:green, @sprintf "\n max Duv for C.div(u): %.2e\n" maximum(Duv))


        end
        ##########################################################################################################"
        
        @views mul!(rhs, B, ph.trans_scalD[:,iscal], 1.0, 1.0) #TODO @views not necessary ?


        # if nb_saved_scalars>1
        #     # ph.saved_scal[:,:,5]=reshape(opC_TL.χ[1].diag,grid)
        #     ph.saved_scal[:,:,2]=reshape(veci(rhs,grid,1), grid)
        # end

        # if nb_saved_scalars>4
        #     # ph.saved_scal[:,:,5]=reshape(opC_TL.χ[1].diag,grid)
        #     ph.saved_scal[:,:,5]=reshape(veci(rhs,grid,1), grid)

        #     if nb_saved_scalars>5
        #         ph.saved_scal[:,:,6]=reshape(opC_TL.χ[1].diag,grid)
        #     end
        # end


        # print("\n test left before A/r L", vecb_L(ph.trans_scalD[:,iscal], grid))
        # print("\n test left before A/r T", vecb_T(ph.trans_scalD[:,iscal], grid))


        @views ph.trans_scalD[:,iscal] .= A \ rhs



        # print("\n test left after A/r L", vecb_L(ph.trans_scalD[:,iscal], grid))
        # print("\n test left after A/r T", vecb_T(ph.trans_scalD[:,iscal], grid))

        # @views ph.trans_scal[:,:,iscal] .= reshape(veci(ph.trans_scalD[:,iscal],grid,2), grid)


        # printstyled(color=:green, @sprintf "\n average c %s\n" average!(ph.trans_scal[:,:,iscal], grid, LS[1].geoL, num))


        @views ph.trans_scal[:,:,iscal] .= reshape(veci(ph.trans_scalD[:,iscal],grid,1), grid)

        # printstyled(color=:cyan, @sprintf "\n after resol\n")
        # print("\n",ph.trans_scal[1,:,iscal])
        # print("\n",reshape(veci(rhs,grid,1), grid)[1,:])
        # print("\n",reshape(veci(rhs,grid,1), grid)[2,:])
        # print("\n",reshape(veci(rhs,grid,1), grid)[3,:]) 

        nonzero = veci(ph.trans_scalD[:,iscal],grid,2)[abs.(veci(ph.trans_scalD[:,iscal],grid,2)) .> 0.0]
        # print("nonzero\n")
        # print(nonzero)
        print("\n mean ",mean(nonzero))



        if iscal!=3 #H2O consummed at the electrode, would need to make distinction to make sure the decrease in H2O is physical or not
          
          
            @views kill_dead_cells_val!(ph.trans_scal[:,:,iscal], grid, LS[1].geoL,concentration0[iscal]) 

            # @views veci(ph.trans_scalD[:,iscal],grid,1) .= vec(ph.trans_scal[:,:,iscal])

            if any(ph.trans_scal[:,:,iscal].<concentration0[iscal]*(1-num.concentration_check_factor))
                print("iscal ",iscal)
                printstyled(color=:red, @sprintf "\n concentration: %.2e %.2e \n" minimum(ph.trans_scal[:,:,iscal]) concentration0[iscal]*(1-num.concentration_check_factor))
                # printstyled(color=:red, @sprintf "\n concentration drop: %.2e%% \n" (minimum(ph.trans_scal[:,:,iscal])-concentration0[iscal])/concentration0[iscal]*100)
                @error("concentration too low")
                # return current_i
            end

            printstyled(color=:red, @sprintf "\n concentration variation vs min: %.2e%% \n" (minimum(ph.trans_scal[:,:,iscal])-concentration0[iscal])/concentration0[iscal]*100)

        else 
            @views kill_dead_cells_val!(ph.trans_scal[:,:,iscal], grid, LS[1].geoL,concentration0[iscal]) 

            if any(ph.trans_scal[:,:,iscal].> concentration0[iscal]*(1+num.concentration_check_factor))
                print("iscal ",iscal)
                printstyled(color=:red, @sprintf "\n concentration: %.2e %.2e \n" maximum(ph.trans_scal[:,:,iscal]) concentration0[iscal]*(1-num.concentration_check_factor))
                # printstyled(color=:red, @sprintf "\n concentration increase: %.2e%% \n" (maximum(ph.trans_scal[:,:,iscal])-concentration0[iscal])/concentration0[iscal]*100)
                @error("concentration too high")
                # return current_i

                for jplot in 1:grid.ny
                    for iplot in 1:grid.nx
                        II = CartesianIndex(jplot, iplot) #(id_y, id_x)
                        pII = lexicographic(II, grid.ny)

                        if ph.trans_scalD[pII,iscal] >concentration0[iscal]
                            printstyled(color=:green, @sprintf "\n j: %5i %5i %.2e %.2e %.2e %.2e \n" iplot jplot grid.x[iplot]/num.plot_xscale grid.y[jplot]/num.plot_xscale ph.trans_scalD[pII,iscal] rhs[pII])
                    
                        end
                    end
                end




            end
            printstyled(color=:red, @sprintf "\n concentration variation vs max: %.2e%% \n" (maximum(ph.trans_scal[:,:,iscal])-concentration0[iscal])/concentration0[iscal]*100)
        end

        @views kill_dead_cells_val!(ph.trans_scal[:,:,iscal], grid, LS[1].geoL,0.0)  #reset to zero (for plot)
        #kill_dead_cells_val! can be used for display to avoid displaying a range from 0 to c0 in python


        # print("\n vecb_L(elec_condD, grid) after res 0 \n ", vecb_L(ph.trans_scalD[:,2], grid) )


        # print("all\n")
        # print(ph.trans_scalD[:,iscal])

        # print("\n test", ph.trans_scal[:,1,iscal])
        # print("\n test", ph.trans_scal[1,:,iscal])

        # print("\n test left", vecb_L(ph.trans_scalD[:,iscal], grid))
        # print("\n test right", vecb_R(ph.trans_scalD[:,iscal], grid))
        # print("\n test bottom", vecb_B(ph.trans_scalD[:,iscal], grid))
        # print("\n test top", vecb_T(ph.trans_scalD[:,iscal], grid))

        # print("\n test v bottom", maximum(vecb_B(ph.vD, grid)))
        # print("\n test v top", maximum(vecb_T(ph.vD, grid)))
        # print("\n test v 1 ", maximum(ph.v[1,:]))
        # print("\n test v 2 ", maximum(ph.v[1,:]))

        print("\n test concentration ", minimum(ph.trans_scal[:,:,iscal]), " ", maximum(ph.trans_scal[:,:,iscal]))


        # print("\n test diff", ph.v[1,:] .- vecb_T(ph.vD, grid))



        # for jplot in 1:ny
        #     for iplot in 1:nx

                # II = CartesianIndex(jplot, iplot) #(id_y, id_x)
        #         pII = lexicographic(II, grid.ny)

        #         if ph.trans_scalD[pII,iscal] < 0.0
        #         # if ph.trans_scalD[pII,iscal] < concentration0[iscal]

        #             printstyled(color=:green, @sprintf "\n j: %5i %5i %.2e %.2e %.2e %.2e \n" iplot jplot grid.x[iplot]/num.plot_xscale grid.y[jplot]/num.plot_xscale ph.trans_scalD[pII,iscal] rhs[pII])
            
        #         end
        #     end
        # end

        # printstyled(color=:red, @sprintf "\n levelset: end scal set_scalar_transport!\n")
        # println(grid.LS[1].geoL.dcap[1,1,:])

    end #end loop iscal


    # printstyled(color=:red, @sprintf "\n levelset: end scal set_scalar_transport!\n")
    # println(grid.LS[1].geoL.dcap[1,1,:])

    return nothing
end


"""    
scalar transport (convection and diffusion)
"""
function scalar_transport!(bc, num, grid, op, geo, ph, concentration0, MIXED, projection,
    op_conv, grid_u, geo_u, grid_v, geo_v,
    periodic_x, periodic_y, convection, ls_advection, BC_int, diffusion_coeff, convection_Cdivu)
    @unpack τ, aniso, nb_transported_scalars,nb_saved_scalars = num
    @unpack nx, ny, dx, dy, ind, LS  = grid
    @unpack all_indices, inside, b_left, b_bottom, b_right, b_top = ind
    @unpack Bx, By, BxT, ByT, Hx, Hy, HxT, HyT, M, iMx, iMy, χ = op
    @unpack CT, CUTCT = op_conv
    @unpack u, v, uD, vD = ph

    printstyled(color=:red, @sprintf "\n levelset: start scalar_transport!\n")
    println(grid.LS[1].geoL.dcap[1,1,:])

    ni = nx * ny
    nb = 2 * nx + 2 * ny
    nt = 2 * ni + nb

    A = spzeros(nt, nt)
    B = spzeros(nt, nt)

    all_CUTCT = zeros(grid.ny * grid.nx, nb_transported_scalars)

    ######################################################################################################
    # Operators
    ######################################################################################################

    if convection

        bcU = zeros(grid_u)
        bcU .= reshape(vec2(uD,grid_u), grid_u)
        bcU[1,:] .= vecb_B(uD, grid_u)
        bcU[end,:] .= vecb_T(uD, grid_u)
        bcU[:,1] .= vecb_L(uD, grid_u)
        bcU[:,end] .= vecb_R(uD, grid_u)
        
        bcV = zeros(grid_v)
        bcV .= reshape(vec2(vD,grid_v), grid_v)
        bcV[:,1] .= vecb_L(vD, grid_v)
        bcV[:,end] .= vecb_R(vD, grid_v)
        bcV[1,:] .= vecb_B(vD, grid_v)
        bcV[end,:] .= vecb_T(vD, grid_v)

        # print("\n vecb_L ", vecb_L(vD, grid_v))
        # print("\n vecb_B ", vecb_B(vD, grid_v))

        HT = zeros(grid)
        # @inbounds @threads for II in vcat(b_left[1], b_bottom[1], b_right[1], b_top[1])
        @inbounds for II in vcat(b_left[1], b_bottom[1], b_right[1], b_top[1])
            HT[II] = distance(grid.LS[1].mid_point[II], geo.centroid[II], dx[II], dy[II])
        end  

        iscal = 1

        # printstyled(color=:red, @sprintf "\n levelset: before scalar_convection\n")
        # println(grid.LS[1].geoL.dcap[1,1,:])

        ######################################################################################################
        #Interface boundary condition
        ######################################################################################################
        if is_dirichlet(bc[iscal].int)
            __a0 = bc[iscal].int.val        
        elseif is_neumann(bc[iscal].int)
            __a0 = bc[iscal].int.val
        elseif is_robin(bc[iscal].int)
            __a0 = bc_type.val  
        elseif is_stefan(bc[iscal].int)
            __a0 = concentration0[iscal]       
        elseif is_wall(bc[iscal].int)
            __a0 = bc[iscal].int.val         
        else
            __a0 = bc[iscal].int.val
        end

        # Flags with BCs
        a0 = ones(grid) .* __a0

        # The idea there is to have at the corners the contributions from both borders. 
        # bcTx is only used for derivatives in x and bcTy for derivatives in y, so it doesn't matter 
        # that bcTx is not filled at the top and bottom and the same for bcTy
        # Otherwise we would have to select one of the contributions
        bcTx, bcTy = set_bc_bnds(dir, a0, HT, bc[iscal])

        # @views necessary
        @views scalar_convection!(dir, CT, all_CUTCT[:,iscal], u, v, bcTx, bcTy, bcU, bcV, geo.dcap, ny, 
        bc[iscal], inside, b_left[1], b_bottom[1], b_right[1], b_top[1]
        )

        if nb_transported_scalars>1
            for iscal=2:nb_transported_scalars

                # printstyled(color=:red, @sprintf "\n levelset: before scalar_convection\n")
                # println(grid.LS[1].geoL.dcap[1,1,:])

                ######################################################################################################
                #Interface boundary condition
                ######################################################################################################
                if is_dirichlet(bc[iscal].int)
                    __a0 = bc[iscal].int.val        
                elseif is_neumann(bc[iscal].int)
                    __a0 = bc[iscal].int.val
                elseif is_robin(bc[iscal].int)
                    __a0 = bc_type.val  
                elseif is_stefan(bc[iscal].int)
                    __a0 = concentration0[iscal]       
                elseif is_wall(bc[iscal].int)
                    __a0 = bc[iscal].int.val         
                else
                    __a0 = bc[iscal].int.val
                end

                # Flags with BCs
                a0 = ones(grid) .* __a0

                bcTx, bcTy = set_bc_bnds(dir, a0, HT, bc[iscal])
                # @views necessary
                @views scalar_convection_CUTCT!(dir, all_CUTCT[:,iscal], u, v, bcTx, bcTy, bcU, bcV, geo.dcap, ny, 
                bc[iscal], inside, b_left[1], b_bottom[1], b_right[1], b_top[1]
                )
            end
        end       
    end #convection

    # printstyled(color=:red, @sprintf "\n levelset: end scalar_convection\n")
    # println(grid.LS[1].geoL.dcap[1,1,:])


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
    ######################################################################################################

    #either update here or compute 
    # update_all_ls_data(num, grid, grid_u, grid_v, BC_int, periodic_x, periodic_y)

    for iscal=1:nb_transported_scalars

        # printstyled(color=:red, @sprintf "\n levelset: start iscal!\n")
        # println(grid.LS[1].geoL.dcap[1,1,:])
    
        diffusion_coeff_scal = diffusion_coeff[iscal]

        # if iscal==3
        #     printstyled(color=:red, @sprintf "\n levelset: forcing diff 0\n")
        #     diffusion_coeff_scal = 0.0
        # end

        ######################################################################################################
        #Interface boundary condition
        ######################################################################################################
        if is_dirichlet(bc[iscal].int)
            # printstyled(color=:green, @sprintf "\n Dirichlet : %.2e\n" bc[iscal].int.val )
            __a0 = bc[iscal].int.val
            __a1 = -1.0
            __b = 0.0
        elseif is_neumann(bc[iscal].int)
            __a0 = bc[iscal].int.val
            __a1 = 0.0
            __b = 1.0
        elseif is_robin(bc[iscal].int)
            __a0 = bc_type.val
            __a1 = -1.0
            __b = 1.0
        elseif is_stefan(bc[iscal].int)
            __a0 = concentration0[iscal]
            __a1 = -1.0
            __b = 0.0
        elseif is_wall(bc[iscal].int)
            __a0 = bc[iscal].int.val
            __a1 = -1.0
            __b = 0.0
        else
            __a0 = bc[iscal].int.val
            __a1 = -1.0
            __b = 0.0
        end

        # Flags with BCs
        a0 = ones(grid) .* __a0

        ######################################################################################################
        #Wall BC
        ######################################################################################################
        _a1 = ones(grid) .* __a1
        a1 = Diagonal(vec(_a1))
        _b = ones(grid) .* __b
        b = Diagonal(vec(_b))

        a0_b = zeros(nb)
        _a1_b = zeros(nb)
        _b_b = zeros(nb)
        set_borders!(grid, grid.LS[1].cl, grid.LS[1].u, a0_b, _a1_b, _b_b, bc[iscal], num.n_ext_cl)
        a1_b = Diagonal(vec(_a1_b))
        b_b = Diagonal(vec(_b_b))

        # if convection
        #     # HT = zeros(grid)
        #     # # @inbounds @threads for II in vcat(b_left[1], b_bottom[1], b_right[1], b_top[1])
        #     # @inbounds for II in vcat(b_left[1], b_bottom[1], b_right[1], b_top[1])
        #     #     HT[II] = distance(grid.LS[1].mid_point[II], geo.centroid[II], dx[II], dy[II])
        #     # end    

        #     # The idea there is to have at the corners the contributions from both borders. 
        #     # bcTx is only used for derivatives in x and bcTy for derivatives in y, so it doesn't matter 
        #     # that bcTx is not filled at the top and bottom and the same for bcTy
        #     # Otherwise we would have to select one of the contributions

        #     bcTx, bcTy = set_bc_bnds(dir, a0, HT, bc[iscal])
        
        #     bcU = zeros(grid_u)
        #     bcU .= reshape(vec2(uD,grid_u), grid_u)
        #     bcU[1,:] .= vecb_B(uD, grid_u)
        #     bcU[end,:] .= vecb_T(uD, grid_u)
        #     bcU[:,1] .= vecb_L(uD, grid_u)
        #     bcU[:,end] .= vecb_R(uD, grid_u)
            
        #     bcV = zeros(grid_v)
        #     bcV .= reshape(vec2(vD,grid_v), grid_v)
        #     bcV[:,1] .= vecb_L(vD, grid_v)
        #     bcV[:,end] .= vecb_R(vD, grid_v)
        #     bcV[1,:] .= vecb_B(vD, grid_v)
        #     bcV[end,:] .= vecb_T(vD, grid_v)

        #     print("\n vecb_L ", vecb_L(vD, grid_v))

        #     print("\n vecb_B ", vecb_B(vD, grid_v))

     

        #     printstyled(color=:red, @sprintf "\n levelset: before scalar_convection\n")
        #     println(grid.LS[1].geoL.dcap[1,1,:])

        #     scalar_convection!(dir, CT, CUTCT, u, v, bcTx, bcTy, bcU, bcV, geo.dcap, ny, 
        #     bc[iscal], inside, b_left[1], b_bottom[1], b_right[1], b_top[1]
        #     )

        # end

        LT = BxT * iMx * Bx .+ ByT * iMy * By
        LD = BxT * iMx * Hx[1] .+ ByT * iMy * Hy[1]
        LD_b = BxT * op.iMx_b * op.Hx_b .+ ByT * op.iMy_b * op.Hy_b

        # A = spzeros(nt, nt)
        # B = spzeros(nt, nt)

        #TODO check no need to reinitialize A and B

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
        B[1:ni,1:ni] = M .+ 0.5 .* τ .* diffusion_coeff_scal .* LT .- τ .* CT
        B[1:ni,ni+1:2*ni] = 0.5 .* τ .* diffusion_coeff_scal .* LD
        B[1:ni,end-nb+1:end] = 0.5 .* τ .* diffusion_coeff_scal .* LD_b

        rhs = fnzeros(grid, num)
        if convection
            vec1(rhs,grid) .-= τ .* all_CUTCT[:,iscal]
        end

        # IItest = CartesianIndex(66, 4) #(id_y, id_x)
        # pIItest = lexicographic(IItest, grid.ny)
        # print("\n rhs",vec2(rhs,grid)[pIItest])

        vec2(rhs,grid) .+= χ[1] * vec(a0)

        # print("\n rhs",vec2(rhs,grid)[pIItest])

        vecb(rhs,grid) .+= op.χ_b * vec(a0_b)
        
        # print("\n rhs",vec2(rhs,grid)[pIItest])
        # print("\n rhs",bc_type.val*grid.dx[1,1])


        if convection_Cdivu
            # Duv = fzeros(grid)
            # Duv = fnzeros(grid,num)

            #TODO BC missing

            Duv = opC_pL.AxT * vec1(ph.uD,grid_u) .+ opC_pL.Gx_b * vecb(ph.uD,grid_u) .+
            opC_pL.AyT * vec1(ph.vD,grid_v) .+ opC_pL.Gy_b * vecb(ph.vD,grid_v)
            for iLS in 1:nLS
                if !is_navier(BC_int[iLS]) && !is_navier_cl(BC_int[iLS])
                    Duv .+= opC_pL.Gx[iLS] * veci(ph.uD,grid_u,iLS+1) .+ 
                            opC_pL.Gy[iLS] * veci(ph.vD,grid_v,iLS+1)
                end
            end
        
            # rhs .+= Duv #.* ph.trans_scalD[:,iscal] multiplied just after

            # vec1(rhs_ϕ,grid) .= rho1 .* iτ .* Duv #TODO
            vec1(rhs,grid) .+= Duv 

            printstyled(color=:green, @sprintf "\n max Duv for C.div(u): %.2e\n" maximum(Duv))


        end
        ##########################################################################################################"
        
        @views mul!(rhs, B, ph.trans_scalD[:,iscal], 1.0, 1.0) #TODO @views not necessary ?


        # if nb_saved_scalars>1
        #     # ph.saved_scal[:,:,5]=reshape(opC_TL.χ[1].diag,grid)
        #     ph.saved_scal[:,:,2]=reshape(veci(rhs,grid,1), grid)
        # end

        # if nb_saved_scalars>4
        #     # ph.saved_scal[:,:,5]=reshape(opC_TL.χ[1].diag,grid)
        #     ph.saved_scal[:,:,5]=reshape(veci(rhs,grid,1), grid)

        #     if nb_saved_scalars>5
        #         ph.saved_scal[:,:,6]=reshape(opC_TL.χ[1].diag,grid)
        #     end
        # end


        # print("\n test left before A/r L", vecb_L(ph.trans_scalD[:,iscal], grid))
        # print("\n test left before A/r T", vecb_T(ph.trans_scalD[:,iscal], grid))


        @views ph.trans_scalD[:,iscal] .= A \ rhs



        # print("\n test left after A/r L", vecb_L(ph.trans_scalD[:,iscal], grid))
        # print("\n test left after A/r T", vecb_T(ph.trans_scalD[:,iscal], grid))

        # @views ph.trans_scal[:,:,iscal] .= reshape(veci(ph.trans_scalD[:,iscal],grid,2), grid)


        # printstyled(color=:green, @sprintf "\n average c %s\n" average!(ph.trans_scal[:,:,iscal], grid, LS[1].geoL, num))


        @views ph.trans_scal[:,:,iscal] .= reshape(veci(ph.trans_scalD[:,iscal],grid,1), grid)

        # printstyled(color=:cyan, @sprintf "\n after resol\n")
        # print("\n",ph.trans_scal[1,:,iscal])
        # print("\n",reshape(veci(rhs,grid,1), grid)[1,:])
        # print("\n",reshape(veci(rhs,grid,1), grid)[2,:])
        # print("\n",reshape(veci(rhs,grid,1), grid)[3,:]) 

        nonzero = veci(ph.trans_scalD[:,iscal],grid,2)[abs.(veci(ph.trans_scalD[:,iscal],grid,2)) .> 0.0]
        # print("nonzero\n")
        # print(nonzero)
        print("\n mean ",mean(nonzero))



        if iscal!=3 #H2O consummed at the electrode, would need to make distinction to make sure the decrease in H2O is physical or not
          
          
            @views kill_dead_cells_val!(ph.trans_scal[:,:,iscal], grid, LS[1].geoL,concentration0[iscal]) 

            # @views veci(ph.trans_scalD[:,iscal],grid,1) .= vec(ph.trans_scal[:,:,iscal])

            if any(ph.trans_scal[:,:,iscal].<concentration0[iscal]*(1-num.concentration_check_factor))
                print("iscal ",iscal)
                printstyled(color=:red, @sprintf "\n concentration: %.2e %.2e \n" minimum(ph.trans_scal[:,:,iscal]) concentration0[iscal]*(1-num.concentration_check_factor))
                # printstyled(color=:red, @sprintf "\n concentration drop: %.2e%% \n" (minimum(ph.trans_scal[:,:,iscal])-concentration0[iscal])/concentration0[iscal]*100)
                @error("concentration too low")
                # return current_i
            end

            printstyled(color=:red, @sprintf "\n concentration variation vs min: %.2e%% \n" (minimum(ph.trans_scal[:,:,iscal])-concentration0[iscal])/concentration0[iscal]*100)

        else 
            @views kill_dead_cells_val!(ph.trans_scal[:,:,iscal], grid, LS[1].geoL,concentration0[iscal]) 

            if any(ph.trans_scal[:,:,iscal].> concentration0[iscal]*(1+num.concentration_check_factor))
                print("iscal ",iscal)
                printstyled(color=:red, @sprintf "\n concentration: %.2e %.2e \n" maximum(ph.trans_scal[:,:,iscal]) concentration0[iscal]*(1-num.concentration_check_factor))
                # printstyled(color=:red, @sprintf "\n concentration increase: %.2e%% \n" (maximum(ph.trans_scal[:,:,iscal])-concentration0[iscal])/concentration0[iscal]*100)
                @error("concentration too high")
                # return current_i

                for jplot in 1:grid.ny
                    for iplot in 1:grid.nx
                        II = CartesianIndex(jplot, iplot) #(id_y, id_x)
                        pII = lexicographic(II, grid.ny)

                        if (ph.trans_scalD[pII,iscal] >concentration0[iscal]*(1+num.concentration_check_factor))

                        # if (ph.trans_scalD[pII,iscal] >concentration0[iscal]) 
                            
                            # && (ph.trans_scalD[pII,iscal] == maximum(ph.trans_scalD[:,iscal]) )
                            printstyled(color=:green, @sprintf "\n j: %5i %5i %.2e %.2e %.2e %.2e \n" iplot jplot grid.x[iplot]/num.plot_xscale grid.y[jplot]/num.plot_xscale ph.trans_scalD[pII,iscal] rhs[pII])
                            
                            @views scalar_debug!(dir, CT, all_CUTCT[:,iscal], u, v, bcTx, bcTy, bcU, bcV, geo.dcap, ny, 
                            bc[iscal], inside, b_left[1], b_bottom[1], b_right[1], b_top[1],num,grid,iplot,jplot)
                        end
                    end
                end

            end
            printstyled(color=:red, @sprintf "\n concentration variation vs max: %.2e%% \n" (maximum(ph.trans_scal[:,:,iscal])-concentration0[iscal])/concentration0[iscal]*100)
        end

        @views kill_dead_cells_val!(ph.trans_scal[:,:,iscal], grid, LS[1].geoL,0.0)  #reset to zero (for plot)
        #kill_dead_cells_val! can be used for display to avoid displaying a range from 0 to c0 in python


        # print("\n vecb_L(elec_condD, grid) after res 0 \n ", vecb_L(ph.trans_scalD[:,2], grid) )


        # print("all\n")
        # print(ph.trans_scalD[:,iscal])

        # print("\n test", ph.trans_scal[:,1,iscal])
        # print("\n test", ph.trans_scal[1,:,iscal])

        # print("\n test left", vecb_L(ph.trans_scalD[:,iscal], grid))
        # print("\n test right", vecb_R(ph.trans_scalD[:,iscal], grid))
        # print("\n test bottom", vecb_B(ph.trans_scalD[:,iscal], grid))
        # print("\n test top", vecb_T(ph.trans_scalD[:,iscal], grid))

        # print("\n test v bottom", maximum(vecb_B(ph.vD, grid)))
        # print("\n test v top", maximum(vecb_T(ph.vD, grid)))
        # print("\n test v 1 ", maximum(ph.v[1,:]))
        # print("\n test v 2 ", maximum(ph.v[1,:]))

        print("\n test concentration ", minimum(ph.trans_scal[:,:,iscal]), " ", maximum(ph.trans_scal[:,:,iscal]))


        # print("\n test diff", ph.v[1,:] .- vecb_T(ph.vD, grid))



        # for jplot in 1:ny
        #     for iplot in 1:nx

                # II = CartesianIndex(jplot, iplot) #(id_y, id_x)
        #         pII = lexicographic(II, grid.ny)

        #         if ph.trans_scalD[pII,iscal] < 0.0
        #         # if ph.trans_scalD[pII,iscal] < concentration0[iscal]

        #             printstyled(color=:green, @sprintf "\n j: %5i %5i %.2e %.2e %.2e %.2e \n" iplot jplot grid.x[iplot]/num.plot_xscale grid.y[jplot]/num.plot_xscale ph.trans_scalD[pII,iscal] rhs[pII])
            
        #         end
        #     end
        # end

        # printstyled(color=:red, @sprintf "\n levelset: end scal set_scalar_transport!\n")
        # println(grid.LS[1].geoL.dcap[1,1,:])

    end #end loop iscal


    # printstyled(color=:red, @sprintf "\n levelset: end scal set_scalar_transport!\n")
    # println(grid.LS[1].geoL.dcap[1,1,:])

    return nothing
end

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



    printstyled(color=:red, @sprintf "\n levelset: start set_scalar_transport!\n")
    println(grid.LS[1].geoL.dcap[1,1,:])

    ni = nx * ny
    nb = 2 * nx + 2 * ny
    nt = 2 * ni + nb

    ######################################################################################################
    #Interface boundary condition
    ######################################################################################################
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

    ######################################################################################################
    #Wall BC
    ######################################################################################################
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
        # @inbounds @threads for II in vcat(b_left[1], b_bottom[1], b_right[1], b_top[1])
        @inbounds for II in vcat(b_left[1], b_bottom[1], b_right[1], b_top[1])
            HT[II] = distance(grid.LS[1].mid_point[II], geo.centroid[II], dx[II], dy[II])
        end    
        bcTx, bcTy = set_bc_bnds(dir, a0, HT, BC_T)
    
        # bnds_u = [grid_u.ind.b_left[1], grid_u.ind.b_bottom[1], grid_u.ind.b_right[1], grid_u.ind.b_top[1]]
        # bnds_v = [grid_v.ind.b_left[1], grid_v.ind.b_bottom[1], grid_v.ind.b_right[1], grid_v.ind.b_top[1]]
        # Δu = [grid_u.dx[1,1], grid_u.dy[1,1], grid_u.dx[end,end], grid_u.dy[end,end]] .* 0.5
        # Δv = [grid_v.dx[1,1], grid_v.dy[1,1], grid_v.dx[end,end], grid_v.dy[end,end]] .* 0.5
        
        # Hu = zeros(grid_u)
        # for i in eachindex(bnds_u)
        #     for II in bnds_u[i]
        #         Hu[II] = Δu[i]
        #     end
        # end

        # Hv = zeros(grid_v)
        # for i in eachindex(bnds_v)
        #     for II in bnds_v[i]
        #         Hv[II] = Δv[i]
        #     end
        # end
    
        bcU = zeros(grid_u)
        bcU .= reshape(vec2(uD,grid_u), grid_u)
        bcU[1,:] .= vecb_B(uD, grid_u)
        bcU[end,:] .= vecb_T(uD, grid_u)
        bcU[:,1] .= vecb_L(uD, grid_u)
        bcU[:,end] .= vecb_R(uD, grid_u)
        
        bcV = zeros(grid_v)
        bcV .= reshape(vec2(vD,grid_v), grid_v)
        bcV[:,1] .= vecb_L(vD, grid_v)
        bcV[:,end] .= vecb_R(vD, grid_v)
        bcV[1,:] .= vecb_B(vD, grid_v)
        bcV[end,:] .= vecb_T(vD, grid_v)

        print("\n vecb_L ", vecb_L(vD, grid_v))

        print("\n vecb_B ", vecb_B(vD, grid_v))

        #TODO
        # print("\n test")
        # bcTx, bcTy = set_bc_bnds_total(dir, a0, HT, BC_T)

        printstyled(color=:red, @sprintf "\n levelset: before scalar_convection\n")
        println(grid.LS[1].geoL.dcap[1,1,:])

        scalar_convection!(dir, CT, CUTCT, u, v, bcTx, bcTy, bcU, bcV, geo.dcap, ny, 
            BC_T, inside, b_left[1], b_bottom[1], b_right[1], b_top[1]
        )


        # #####################################################################################################"
        # Du_x = zeros(grid_u)
        # # Du_y = zeros(grid_u)
        # # bcU .= reshape(vec2(uD,grid_u), grid_u)
        
        # for iLS in 1:num.nLS
        #     Du_x[grid_u.LS[iLS].MIXED] .= reshape(veci(uD,grid_u,iLS+1), grid_u)[grid_u.LS[iLS].MIXED]
        #     # Du_y[grid_u.LS[iLS].MIXED] .= reshape(veci(uD,grid_u,iLS+1), grid_u)[grid_u.LS[iLS].MIXED]
        # end
        # Du_x[:,1] .= vecb_L(uD,grid_u)
        # # Du_y[1,:] .= vecb_B(uD,grid_u)
        # Du_x[:,end] .= vecb_R(uD,grid_u)
        # # Du_y[end,:] .= vecb_T(uD,grid_u)
        # # Dv_x = zeros(grid_v)
        # Dv_y = zeros(grid_v)
        # for iLS in 1:num.nLS
        #     # Dv_x[grid_v.LS[iLS].MIXED] .= reshape(veci(vD,grid_v,iLS+1), grid_v)[grid_v.LS[iLS].MIXED]
        #     Dv_y[grid_v.LS[iLS].MIXED] .= reshape(veci(vD,grid_v,iLS+1), grid_v)[grid_v.LS[iLS].MIXED]
        # end
        # # Dv_x[:,1] .= vecb_L(vD,grid_v)
        # Dv_y[1,:] .= vecb_B(vD,grid_v)
        # # Dv_x[:,end] .= vecb_R(vD,grid_v)
        # Dv_y[end,:] .= vecb_T(vD,grid_v)
        # scalar_convection!(dir, CT, CUTCT, u, v, bcTx, bcTy, Du_x, Dv_y, geo.dcap, ny, 
        #     BC_T, inside, b_left[1], b_bottom[1], b_right[1], b_top[1]
        # )
        # #####################################################################################################"

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
    B[1:ni,1:ni] = M .+ 0.5 .* τ .* diffusion_coeff_scal .* LT .- τ .* CT
    B[1:ni,ni+1:2*ni] = 0.5 .* τ .* diffusion_coeff_scal .* LD
    B[1:ni,end-nb+1:end] = 0.5 .* τ .* diffusion_coeff_scal .* LD_b

    rhs = fnzeros(grid, num)
    if convection
        vec1(rhs,grid) .-= τ .* CUTCT
    end

    #IItest = CartesianIndex(66, 4) #(id_y, id_x)
    # IItest = CartesianIndex(69, 3) #(id_y, id_x)
    # pIItest = lexicographic(IItest, grid.ny)

    # print("\n rhs",vec2(rhs,grid)[pIItest])

    vec2(rhs,grid) .+= χ[1] * vec(a0)

    # print("\n rhs",vec2(rhs,grid)[pIItest])

    vecb(rhs,grid) .+= op.χ_b * vec(a0_b)
    
    # print("\n rhs",vec2(rhs,grid)[pIItest])
    # print("\n rhs",bc_type.val*grid.dx[1,1])


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
function electrolysis_velocity!(num, grid, LS, V, TL, MIXED, periodic_x, periodic_y, concentration_scal_intfc)
    @unpack geoS, geoL, κ = LS

    V .= 0
    # @inbounds @threads for II in MIXED
    @inbounds for II in MIXED

        θ_d = concentration_scal_intfc

        # dTS = 0.
        dTL = 0.
        if geoL.projection[II].flag
            T_1, T_2 = interpolated_temperature(grid, geoL.projection[II].angle, geoL.projection[II].point1, geoL.projection[II].point2, TL, II, periodic_x, periodic_y)
            dTL = normal_gradient(geoL.projection[II].d1, geoL.projection[II].d2, T_1, T_2, θ_d)
        else
            T_1 = interpolated_temperature(grid, geoL.projection[II].angle, geoL.projection[II].point1, TL, II, periodic_x, periodic_y)
            dTL = normal_gradient(geoL.projection[II].d1, T_1, θ_d)
        end
        V[II] = dTL #+ dTS

        print("\n grad",II,"val",dTL,"val",geoL.projection[II].flag,"val",geoL.projection[II].angle,"val", geoL.projection[II].point1,"val", geoL.projection[II].point2,"val", T_1,"val",T_2,"val",θ_d)
    end
    return nothing
end

"""
From update_free_surface_velocity and update_stefan_velocity
"""
function update_free_surface_velocity_electrolysis(num, grid, grid_u, grid_v, iLS, uD, vD, periodic_x, periodic_y, Vmean, concentration_scalD, diffusion_coeff_scal,concentration_scal_intfc)
    @unpack MWH2,rho1,rho2=num
    # @unpack χ,=opC
    #TODO check sign 

    printstyled(color=:green, @sprintf "\n grid p u v max : %.2e %.2e %.2e\n" maximum(abs.(grid.V[grid.LS[iLS].MIXED])) maximum(abs.(grid_u.V[grid.LS[iLS].MIXED])) maximum(abs.(grid_v.V[grid_v.LS[iLS].MIXED])))

    plot_electrolysis_velocity!(num, grid, grid.LS[iLS], grid.V, concentration_scalD, grid.LS[iLS].MIXED, periodic_x, periodic_y, concentration_scal_intfc)


    electrolysis_velocity!(num, grid, grid.LS[iLS], grid.V, concentration_scalD, grid.LS[iLS].MIXED, periodic_x, periodic_y, concentration_scal_intfc)

    # concentration_scalu = zeros(grid_u)
    # concentration_scalv = zeros(grid_v)
    # #interpolate
    # interpolate_scalar!(grid, grid_u, grid_v, reshape(veci(concentration_scalD,grid,1), grid), concentration_scalu, concentration_scalv)
    # # Or give scal directly instead of scalD 
    # # interpolate_scalar!(grid, grid_u, grid_v, reshape(veci(concentration_scal,grid,1), grid), concentration_scalu, concentration_scalv)
    # #TODO H: suppose val on interface for u and v same as for p 
    # electrolysis_velocity!(num, grid_u, grid_u.LS[iLS], grid_u.V, concentration_scalu, grid_u.LS[iLS].MIXED, periodic_x, periodic_y, concentration_scal_intfc)
    # electrolysis_velocity!(num, grid_v, grid_v.LS[iLS], grid_v.V, concentration_scalv, grid_v.LS[iLS].MIXED, periodic_x, periodic_y, concentration_scal_intfc)
    # grid.V[grid.LS[iLS].MIXED] .*= 1. ./ λ
    #H2
    #TODO check order rho1 rho2
    #no surface term

    factor = -(1.0/rho2-1.0/rho1).*diffusion_coeff_scal[1].*MWH2
    
    # grid_u.V[grid_u.LS[iLS].MIXED] .*= factor
    # grid_v.V[grid_v.LS[iLS].MIXED] .*= factor

    grid.V[grid.LS[iLS].MIXED] .*= factor

   
    # if Vmean
    #     a = mean(grid.V[grid.LS[iLS].MIXED])
    #     grid.V[grid.LS[iLS].MIXED] .= a
    # end

    grid_u.V .= reshape(vec1(uD,grid_u), grid_u)
    grid_v.V .= reshape(vec1(vD,grid_v), grid_v)

    print("factor",factor)

    printstyled(color=:green, @sprintf "\n grid p u v max : %.2e %.2e %.2e\n" maximum(abs.(grid.V[grid.LS[iLS].MIXED])) maximum(abs.(grid_u.V[grid.LS[iLS].MIXED])) maximum(abs.(grid_v.V[grid_v.LS[iLS].MIXED])))

    #TODO extension


#     grid_u.V .+= reshape(veci(uD,grid_u,iLS+1), (grid_u.ny, grid_u.nx))
#     grid_v.V .+= reshape(veci(vD,grid_v,iLS+1), (grid_v.ny, grid_v.nx))

#     i_u_ext, l_u_ext, b_u_ext, r_u_ext, t_u_ext = indices_extension(grid_u, grid_u.LS[iLS], grid_u.ind.inside, periodic_x, periodic_y)
#     i_v_ext, l_v_ext, b_v_ext, r_v_ext, t_v_ext = indices_extension(grid_v, grid_v.LS[iLS], grid_v.ind.inside, periodic_x, periodic_y)

#     field_extension!(grid_u, grid_u.LS[iLS].u, grid_u.V, i_u_ext, l_u_ext, b_u_ext, r_u_ext, t_u_ext, num.NB, periodic_x, periodic_y)
#     field_extension!(grid_v, grid_v.LS[iLS].u, grid_v.V, i_v_ext, l_v_ext, b_v_ext, r_v_ext, t_v_ext, num.NB, periodic_x, periodic_y)
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



function print_electrolysis_statistics(nb_transported_scalars,grid,phL)

    # normscal1L = norm(phL.trans_scal[:,:,1])
    # normscal2L = norm(phL.trans_scal[:,:,2])
    # normscal3L = norm(phL.trans_scal[:,:,3])
    # normphi_eleL = norm(phL.phi_ele)
    # normTL = norm(phL.T)
    # normiL=0.0 #TODO

    

    minscal1L=minscal2L=minscal3L=minphi_eleL=minTL=miniL=0.0
    maxscal1L=maxscal2L=maxscal3L=maxphi_eleL=maxTL=maxiL=0.0
    moyscal1L=moyscal2L=moyscal3L=moyphi_eleL=moyTL=moyiL=0.0


    if nb_transported_scalars>0
        minscal1L = minimum(phL.trans_scal[:,:,1])
        maxscal1L = maximum(phL.trans_scal[:,:,1])
        moyscal1L = mean(phL.trans_scal[:,:,1])
    end

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
    # print("mean cH2 interface",mean(veci(phL.trans_scalD[:,1],grid,2)))
    # print("mean cH2 interface",average!(reshape(veci(phL.trans_scalD[:,iscal],grid,2), grid), grid, LS[1].geoL, num))

    # printstyled(color=:green, @sprintf "\n average c %s\n" average!(reshape(veci(phL.trans_scalD[:,iscal],grid,2), grid), grid, LS[1].geoL, num))
    
    # print("nonzero\n")
    # print(nonzero)
    # print("\nmean c(H2) interface \n",mean(nonzero))
    if nb_transported_scalars>0
        nonzero = veci(phL.trans_scalD[:,1],grid,2)[abs.(veci(phL.trans_scalD[:,1],grid,2)) .> 0.0]
        printstyled(color=:green, @sprintf "\n mean c(H2) interface : %.2e\n" mean(nonzero))
    end

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

    #Output everything to phL

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
    elseif adapt_timestep_mode==3
        return num.dt0
    end

    # printstyled(color=:green, @sprintf "\n rho1 %.2e rho2 %.2e mu1 %.2e mu2 %.2e\n" rho1 rho2 mu1 mu2) 
    # printstyled(color=:green, @sprintf "\n dx %.2e dy %.2e \n" min_spacing_x min_spacing_x) 

    
    # max(mu1/rho1, mu2/rho2) * (2.0/min_spacing_x^2 + 2.0/min_spacing_y^2 )

    #  printstyled(color=:green, @sprintf "\n c_conv %.2e c_visc %.2e c_grad %.2e\n" c_conv c_visc c_grav)
    #  printstyled(color=:green, @sprintf "\n CFL : %.2e dt : %.2e\n" CFL num.τ)
end

function init_fields_2!(TD,T,H,BC,grid,dir_val_intfc)

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

# function get_S_height(grid,ind,dx,dy)

#     H = zeros(grid)
#     @inbounds @threads for II in vcat(ind.b_left[1], ind.b_bottom[1], ind.b_right[1], ind.b_top[1])
#         H[II] = distance(grid.LS[1].mid_point[II], geoS.centroid[II], dx[II], dy[II])
#     end   

# end

# function get_L_height(grid,ind,dx,dy)

#     H = zeros(grid)
#     @inbounds @threads for II in vcat(ind.b_left[1], ind.b_bottom[1], ind.b_right[1], ind.b_top[1])
#         H[II] = distance(grid.LS[1].mid_point[II], geoL.centroid[II], dx[II], dy[II])
#     end   

# end

function get_height(grid,ind,dx,dy,geo)

    H = zeros(grid)
    @inbounds @threads for II in vcat(ind.b_left[1], ind.b_bottom[1], ind.b_right[1], ind.b_top[1])
        H[II] = distance(grid.LS[1].mid_point[II], geo.centroid[II], dx[II], dy[II])
    end   

    return H

end

function get_uv_height(grid_u,grid_v)

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

    return Hu,Hv

end

"""    
# Arguments
- bc_type: BC for interface, num, grid, 
- a0, 
- opC, 
- opC_u, 
- pC_v,
- A, 
- L, 
- bc_L, 
- bc_L_b, 
- BC: BC for wall
- ls_advection
"""
function set_poisson_variable_coeff(
    bc_type, num, grid, grid_u, grid_v, a0, opC, opC_u, opC_v,
    A, 
    # L, bc_L, bc_L_b, 
    BC,
    ls_advection,coeffD)
    @unpack Bx, By, Hx, Hy, HxT, HyT, χ, M, iMx, iMy, Hx_b, Hy_b, HxT_b, HyT_b, iMx_b, iMy_b, iMx_bd, iMy_bd, χ_b = opC

    ni = grid.nx * grid.ny
    nb = 2 * grid.nx + 2 * grid.ny

    rhs = fnzeros(grid, num)

    a0_b = zeros(nb)
    _a1_b = zeros(nb)
    _b_b = zeros(nb)
    for iLS in 1:num.nLS
        set_borders!(grid, grid.LS[iLS].cl, grid.LS[iLS].u, a0_b, _a1_b, _b_b, BC, num.n_ext_cl)
    end
    a1_b = Diagonal(vec(_a1_b))
    b_b = Diagonal(vec(_b_b))

    # printstyled(color=:green, @sprintf "\n sizes \n")

    # print(size(A[end-nb+1:end,1:ni]),"\n")
    # print(size(b_b),"\n")
    # print(size(HxT_b),"\n")
    # print(size(iMx_b'),"\n")
    # print(size(vec1(coeffD,grid)),"\n")
    # print(size(Bx),"\n")
    # print(size(HyT_b),"\n")
    # print(size(iMy_b'),"\n")
    # print(size(By),"\n")


    # printstyled(color=:green, @sprintf "\n vecb \n")
    # print(size(vecb(coeffD,grid)),"\n")
    # nt = (num.nLS + 1) * ni + nb

    # print(" ni ",ni," nb ",nb," nt ",nt)


    coeffDu = zeros(grid_u)
    coeffDv = zeros(grid_v)
    coeffDx_interface = zeros(grid_u)
    coeffDy_interface = zeros(grid_v)


    coeffD_borders = vecb(coeffD,grid)
    interpolate_scalar!(grid, grid_u, grid_v, reshape(veci(coeffD,grid,1), grid), coeffDu, coeffDv)
    interpolate_scalar!(grid, grid_u, grid_v, reshape(veci(coeffD,grid,2), grid), coeffDx_interface, coeffDy_interface)
    coeffDx_bulk = veci(coeffDu,grid_u)
    coeffDy_bulk = veci(coeffDv,grid_v)

    mat_coeffDx = Diagonal(vec(coeffDx_bulk)) # coeffDx_bulk is a 2d matrix with shape (grid_u.ny, grid_u.nx), multiplies Bx
    mat_coeffDy = Diagonal(vec(coeffDy_bulk)) # coeffDx_bulk is a 2d matrix with shape (grid_v.ny, grid_v.nx), multiplies By
    mat_coeffDx_i = Diagonal(vec(coeffDx_interface)) # coeffDx_interface is a 2d matrix with shape (grid_u.ny, grid_u.nx), multiplies Hx
    mat_coeffDy_i = Diagonal(vec(coeffDy_interface)) # coeffDx_interface is a 2d matrix with shape (grid_v.ny, grid_v.nx), multiplies Hy
    mat_coeffDx_b = Diagonal(vec(coeffD_borders)) # coeffDx_interface is a 1d vector with shape (2grid.ny + 2grid.nx), multiplies Hx_b and Hy_b


    # printstyled(color=:green, @sprintf "\n sizes \n")

    # print(size(A[end-nb+1:end,1:ni]),"\n")
    # print(size(b_b),"\n")
    # print(size(HxT_b),"\n")
    # print(size(iMx_b'),"\n")
    # print(size(vec1(coeffD,grid)),"\n")
    # print(size(Bx),"\n")
    # print(size(HyT_b),"\n")
    # print(size(iMy_b'),"\n")
    # print(size(By),"\n")



    # printstyled(color=:green, @sprintf "\n vecb \n")
    # print(size(vecb(coeffD,grid)),"\n")
    # nt = (num.nLS + 1) * ni + nb

    # print(" ni ",ni," nb ",nb," nt ",nt)

    # print(size(mat_coeffDx),"\n")
    # print(size(mat_coeffDy),"\n")
    # print(size(mat_coeffDx_i),"\n")
    # print(size(mat_coeffDy_i),"\n")
    # print(size(mat_coeffDx_b),"\n")

    # A[sb,1:ni] = -b * (HxT[iLS] * iMx .* mat_coeffDx * Bx .+ HyT[iLS] * iMy * mat_coeffDy * By) #or vec1

#   A  (1024, 65536)
# b_b (1024, 1024)
# HxT_b (1024, 1024)
# iMx_b' (1024, 65792)
# vec1(coeffD,grid) (65536,)
# Bx (65792, 65536)
# HyT_b (1024, 1024)
# iMy_b' (1024, 65792)
# By (65792, 65536)

# mat_coeffDx (65792, 65792)


    # elec_Lpm1_L, elec_bc_Lpm1_L, elec_bc_Lpm1_b_L, elec_Lum1_L, elec_bc_Lum1_L, elec_bc_Lum1_b_L, elec_Lvm1_L, elec_bc_Lvm1_L, elec_bc_Lvm1_b_L = set_matrices!(
        #     num, grid, geoL, grid_u, geo_uL, grid_v, geo_vL,
        #     opC_pL, opC_uL, opC_vL, periodic_x, periodic_y
        # )

    # Update laplacian at each timestep
    # elec_Lpm1_L = laplacian_variable_coeff(opC_pL,elec_condD)
    # elec_bc_Lpm1_L, elec_bc_Lpm1_b_L = laplacian_bc_variable_coeff(opC_pL, num.nLS, elec_condD)

    # L, bc_L, bc_L_b, 

    # L = laplacian_variable_coeff(opC_pL,coeffD)
    # bc_L, bc_L_b = laplacian_bc_variable_coeff(opC_pL, num.nLS, coeffD)

    # mat_coeffDx = Diagonal(vec(coeffDx_bulk)) # coeffDx_bulk is a 2d matrix with shape (grid_u.ny, grid_u.nx), multiplies Bx
    # mat_coeffDy = Diagonal(vec(coeffDy_bulk)) # coeffDx_bulk is a 2d matrix with shape (grid_v.ny, grid_v.nx), multiplies By



    # @unpack Bx, By, BxT, ByT, iMx, iMy, tmp_x, tmp_y = opC
    # @unpack BxT, ByT, Hx, Hy, iMx, iMy, Hx_b, Hy_b, iMx_b, iMy_b = opC

    @unpack BxT, ByT,tmp_x, tmp_y = opC


    #Laplacian
    mul!(tmp_x, iMx, mat_coeffDx * Bx)
    L = BxT * tmp_x
    mul!(tmp_y, iMy, mat_coeffDy * By)
    L = L .+ ByT * tmp_y

    #Boundary for Laplacian
    bc_L = []
    for iLS in 1:num.nLS
        push!(bc_L, BxT * iMx * mat_coeffDx_i *Hx[iLS] .+ ByT * iMy * mat_coeffDy_i *Hy[iLS])
    end

    bc_L_b = (BxT * iMx_b * mat_coeffDx_b *Hx_b .+ ByT * iMy_b * mat_coeffDx_b *Hy_b)

      
    if ls_advection
        # Poisson equation
        A[1:ni,1:ni] = pad(L, -4.0)
        A[1:ni,end-nb+1:end] = bc_L_b

        # Boundary conditions for outer boundaries
        A[end-nb+1:end,1:ni] = -b_b * (HxT_b * iMx_b' * mat_coeffDx * Bx .+ HyT_b * iMy_b' * mat_coeffDy * By)
        A[end-nb+1:end,end-nb+1:end] = -pad(b_b * (HxT_b * iMx_bd * mat_coeffDx_b * Hx_b .+ HyT_b * iMy_bd * mat_coeffDx_b * Hy_b) .- χ_b * a1_b, 4.0)
    end

    for iLS in 1:num.nLS
        if ls_advection
            if is_dirichlet(bc_type[iLS])
                __a0 = bc_type[iLS].val
                __a1 = -1.0
                __a2 = 0.0
                __b = 0.0
            elseif is_neumann(bc_type[iLS])
                __a0 = bc_type[iLS].val
                __a1 = 0.0
                __a2 = 0.0
                __b = 1.0
            elseif is_robin(bc_type[iLS])
                __a0 = bc_type[iLS].val
                __a1 = -1.0
                __a2 = 0.0
                __b = 1.0
            elseif is_fs(bc_type[iLS])
                print("error not implemented set_poisson_variable_coeff",bc_type[iLS])
                @error ("error set_poisson_variable_coeff")

                __a1 = 0.0
                __a2 = 1.0
                __b = 0.0
            elseif is_wall_no_slip(bc_type[iLS])
                print("error not implemented set_poisson_variable_coeff",bc_type[iLS])
                @error ("error set_poisson_variable_coeff")

                __a1 = 0.0
                __a2 = 0.0
                __b = 1.0
            elseif is_flapping(bc_type[iLS])
                print("error not implemented set_poisson_variable_coeff",bc_type[iLS])
                @error ("error set_poisson_variable_coeff")

                __a1 = 0.0
                __a2 = 0.0
                __b = 1.0
            elseif is_navier(bc_type[iLS])
                print("error not implemented set_poisson_variable_coeff",bc_type[iLS])
                @error ("error set_poisson_variable_coeff")

                __a1 = 0.0
                __a2 = 0.0
                __b = 1.0
            elseif is_navier_cl(bc_type[iLS])
                print("error not implemented set_poisson_variable_coeff",bc_type[iLS])
                @error ("error set_poisson_variable_coeff")

                __a1 = 0.0
                __a2 = 0.0
                __b = 1.0
            else
                __a1 = 0.0
                __a2 = 0.0
                __b = 1.0
            end
    

            # Flags with BCs
            a0 = ones(grid) .* __a0


          
    
            _a1 = ones(grid) .* __a1
            a1 = Diagonal(vec(_a1))
            _b = ones(grid) .* __b
            b = Diagonal(vec(_b))

            a0_b = zeros(nb)
            _a1_b = zeros(nb)
            _b_b = zeros(nb)
            set_borders!(grid, grid.LS[1].cl, grid.LS[1].u, a0_b, _a1_b, _b_b, BC, num.n_ext_cl)
            a1_b = Diagonal(vec(_a1_b))
            b_b = Diagonal(vec(_b_b))

           
            _a2 = ones(grid) .* __a2
            a2 = Diagonal(vec(_a2))
         
            fs_mat = HxT[iLS] * Hx[iLS] .+ HyT[iLS] * Hy[iLS]

            sb = iLS*ni+1:(iLS+1)*ni
            
            # Poisson equation
            A[1:ni,sb] = bc_L[iLS]
            # Boundary conditions for inner boundaries
            A[sb,1:ni] = -b * (HxT[iLS] * iMx * mat_coeffDx * Bx .+ HyT[iLS] * iMy * mat_coeffDy * By) #or vec1
            # Contribution to Neumann BC from other boundaries
            for i in 1:num.nLS
                if i != iLS
                    A[sb,i*ni+1:(i+1)*ni] = -b * (HxT[iLS] * iMx * mat_coeffDx_i * Hx[i] .+ HyT[iLS] * iMy * mat_coeffDy_i * Hy[i])
                end
            end
            A[sb,sb] = -pad(
                b * (HxT[iLS] * iMx * mat_coeffDx_i * Hx[iLS] .+ HyT[iLS] * iMy * mat_coeffDy_i * Hy[iLS]) .- χ[iLS] * a1 .+
                a2 * Diagonal(diag(fs_mat)), 4.0
            )
            A[sb,end-nb+1:end] = b * (HxT[iLS] * iMx_b * mat_coeffDx_b * Hx_b .+ HyT[iLS] * iMy_b * mat_coeffDx_b * Hy_b)
            # Boundary conditions for outer boundaries
            A[end-nb+1:end,sb] = -b_b * (HxT_b * iMx_b' * mat_coeffDx_i * Hx[iLS] .+ HyT_b * iMy_b' * mat_coeffDy_i * Hy[iLS])
        end

        # a0_p = []
        # for i in 1:num.nLS
        #     push!(a0_p, zeros(grid))
        # end

        veci(rhs,grid,iLS+1) .= -χ[iLS] * vec(a0) #vec(a0[iLS])
    end

    vecb(rhs,grid) .= -χ_b * vec(a0_b)
    
    return rhs
end

function laplacian_bc_variable_coeff(opC, nLS, grid, coeffD)
    @unpack BxT, ByT, Hx, Hy, iMx, iMy, Hx_b, Hy_b, iMx_b, iMy_b = opC

    coeffD_borders = vecb(coeffD,grid)

    mat_coeffDx_i = Diagonal(vec(coeffDx_interface)) # coeffDx_interface is a 2d matrix with shape (grid_u.ny, grid_u.nx), multiplies Hx
    mat_coeffDy_i = Diagonal(vec(coeffDy_interface)) # coeffDx_interface is a 2d matrix with shape (grid_v.ny, grid_v.nx), multiplies Hy
    mat_coeffDx_b = Diagonal(vec(coeffD_borders)) # coeffDx_interface is a 1d vector with shape (2grid.ny + 2grid.nx), multiplies Hx_b and Hy_b


    bc_L = []
    for iLS in 1:nLS
        push!(bc_L, BxT * iMx * mat_coeffDx_i *Hx[iLS] .+ ByT * iMy * mat_coeffDy_i *Hy[iLS])
    end

    bc_L_b = (BxT * iMx_b * mat_coeffDx_b *Hx_b .+ ByT * iMy_b * mat_coeffDx_b *Hy_b)

    return bc_L, bc_L_b
end

function laplacian_variable_coeff(opC,coeffD)
    @unpack Bx, By, BxT, ByT, iMx, iMy, tmp_x, tmp_y = opC

    mat_coeffDx = Diagonal(vec(coeffDx_bulk)) # coeffDx_bulk is a 2d matrix with shape (grid_u.ny, grid_u.nx), multiplies Bx
    mat_coeffDy = Diagonal(vec(coeffDy_bulk)) # coeffDx_bulk is a 2d matrix with shape (grid_v.ny, grid_v.nx), multiplies By

    mul!(tmp_x, iMx, mat_coeffDx, Bx)
    L = BxT * tmp_x
    mul!(tmp_y, iMy, mat_coeffDy, By)
    L = L .+ ByT * tmp_y

    return L
end
    


"""
  Compute mass flux like in Khalighi 2023
"""
function compute_mass_flux!(num,grid, grid_u, grid_v, phL, phS,  opC_pL, opC_pS,diffusion_coeff,iscal,geo)
    @unpack  χ = opC_pL
    @unpack nLS = num
    #Liquid phase
    @unpack trans_scalD = phL
    if iscal !=0
        scalD= trans_scalD[:,iscal] #H2 species
    else
        @unpack phi_eleD = phL
        scalD = phi_eleD
    end

    #TODO overwriting chi

    opC_p = opC_pL

    iLStmp=1

    # mass_flux_vec1 = fnzeros(grid,num)
    # mass_flux_vecb = fnzeros(grid,num)
    # mass_flux_veci = fnzeros(grid,num)

    mass_flux_vec1 = fzeros(grid)
    mass_flux_vecb = fzeros(grid)
    mass_flux_veci = fzeros(grid)



    # print(size(mass_flux_vec1),"\n")
    # print(size(mass_flux_vecb),"\n")
    # print(size(mass_flux_veci),"\n")

    mass_flux_vec1   = opC_p.HxT[iLStmp] * opC_p.iMx * opC_p.Bx * vec1(scalD,grid) .+ opC_p.HyT[iLStmp] * opC_p.iMy * opC_p.By * vec1(scalD,grid)
    mass_flux_vecb   = opC_p.HxT[iLStmp] * opC_p.iMx_b * opC_p.Hx_b * vecb(scalD,grid) .+ opC_p.HyT[iLStmp] *  opC_p.iMy_b * opC_p.Hy_b * vecb(scalD,grid)

    # printstyled(color=:red, @sprintf "\n vec1 x y %.2e %.2e\n" sum(opC_p.HxT[iLStmp] * opC_p.iMx * opC_p.Bx * vec1(scalD,grid)) sum(opC_p.HyT[iLStmp] * opC_p.iMy * opC_p.By * vec1(scalD,grid)))

    # print(size(mass_flux_vec1),"\n")
    # print(size(mass_flux_vecb),"\n")
    # print(size(mass_flux_veci),"\n")

    # print(size(opC_p.HxT[1] * opC_p.iMx * opC_p.Hx[1] * veci(scalD,grid,1+1)),"\n")


    for iLS in 1:nLS
        mass_flux_veci .+= opC_p.HxT[iLS] * opC_p.iMx * opC_p.Hx[iLS] * veci(scalD,grid,iLS+1)
        mass_flux_veci .+= opC_p.HyT[iLS] * opC_p.iMy * opC_p.Hy[iLS] * veci(scalD,grid,iLS+1)
    end

    # mass_flux = mass_flux_vec1 .+ mass_flux_vecb .+ mass_flux_veci

    # mass_flux_2 = reshape(mass_flux,grid)
    mass_flux_vec1_2 = reshape(mass_flux_vec1,grid)
    mass_flux_vecb_2 = reshape(mass_flux_vecb,grid)
    mass_flux_veci_2 = reshape(mass_flux_veci,grid)

    mass_flux_2 = mass_flux_vec1_2 .+ mass_flux_vecb_2 .+ mass_flux_veci_2



    ######################################################################################


    # Matrices for interior BCs
    # for iLS in 1:num.nLS
    #     χx = (geo.dcap[:,:,3] .- geo.dcap[:,:,1]) .^ 2
    #     χy = (geo.dcap[:,:,4] .- geo.dcap[:,:,2]) .^ 2
    #     χ[iLS].diag .= sqrt.(vec(χx .+ χy))
    # end

    # χx = (geo.dcap[:,:,3] .- geo.dcap[:,:,1]) .^ 2
    # χy = (geo.dcap[:,:,4] .- geo.dcap[:,:,2]) .^ 2
    # #     χ[iLS].diag .= sqrt.(vec(χx .+ χy))
    # radial_flux_surf = mass_flux_2 ./ sqrt.(vec(χx .+ χy))
    # # radial_flux_surf = mass_flux_2 ./ χ[1]

    # printstyled(color=:green, @sprintf "\n Radial flux: %.2e \n" radial_flux_surf)



    ######################################################################################



    return mass_flux_2, mass_flux_vec1_2, mass_flux_vecb_2, mass_flux_veci_2
end

# function compute_mass_flux!(num,grid, grid_u, grid_v, phL, phS,  opC_pL, opC_pS,diffusion_coeff,iscal)
    
#     @unpack nLS = num
#     #Liquid phase
#     @unpack trans_scalD = phL
#     scalD= trans_scalD[:,iscal] #H2 species
#     opC_p = opC_pL

#     iLStmp=1

#     # mass_flux = fnzeros(grid,num)
#     # mass_flux   = opC_p.HxT[iLStmp] * opC_p.iMx * opC_p.Bx * vec1(scalD,grid) .+ opC_p.HxT[iLStmp] * opC_p.iMx_b * opC_p.Hx_b * vecb(scalD,grid)
#     # mass_flux .+= opC_p.HyT[iLStmp] * opC_p.iMy * opC_p.By * vec1(scalD,grid) .+ opC_p.HyT[iLStmp] *  opC_p.iMy_b * opC_p.Hy_b * vecb(scalD,grid)

#     # printstyled(color=:green, @sprintf "\n mass_flux v1 : %.2e \n" sum(mass_flux))

#     # for iLS in 1:nLS
#     #     mass_flux .+= opC_p.HxT[iLS] * opC_p.iMx * opC_p.Hx[iLS] * veci(scalD,grid,iLS+1)
#     #     mass_flux .+= opC_p.HyT[iLS] * opC_p.iMy * opC_p.Hy[iLS] * veci(scalD,grid,iLS+1)

#     #     printstyled(color=:green, @sprintf "\n mass_flux veci : %.2e %.2e\n" sum(opC_p.HxT[iLS] * opC_p.iMx * opC_p.Hx[iLS] * veci(scalD,grid,iLS+1)) sum(opC_p.HyT[iLS] * opC_p.iMy * opC_p.Hy[iLS] * veci(scalD,grid,iLS+1)))
#     # end

#     mass_flux_vec1 = fnzeros(grid,num)
#     mass_flux_vecb = fnzeros(grid,num)
#     mass_flux_veci = fnzeros(grid,num)

#     mass_flux_vec1   = opC_p.HxT[iLStmp] * opC_p.iMx * opC_p.Bx * vec1(scalD,grid) .+ opC_p.HyT[iLStmp] * opC_p.iMy * opC_p.By * vec1(scalD,grid)
#     mass_flux_vecb   = opC_p.HxT[iLStmp] * opC_p.iMx_b * opC_p.Hx_b * vecb(scalD,grid) .+ opC_p.HyT[iLStmp] *  opC_p.iMy_b * opC_p.Hy_b * vecb(scalD,grid)

#     # printstyled(color=:green, @sprintf "\n mass_flux v1 : %.2e \n" sum(mass_flux))
    

#     for iLS in 1:nLS
#         mass_flux_veci .+= opC_p.HxT[iLS] * opC_p.iMx * opC_p.Hx[iLS] * veci(scalD,grid,iLS+1)
#         mass_flux_veci .+= opC_p.HyT[iLS] * opC_p.iMy * opC_p.Hy[iLS] * veci(scalD,grid,iLS+1)

#         # printstyled(color=:green, @sprintf "\n mass_flux veci : %.2e %.2e\n" sum(opC_p.HxT[iLS] * opC_p.iMx * opC_p.Hx[iLS] * veci(scalD,grid,iLS+1)) sum(opC_p.HyT[iLS] * opC_p.iMy * opC_p.Hy[iLS] * veci(scalD,grid,iLS+1)))
#     end

#     mass_flux = mass_flux_vec1 + mass_flux_vecb + mass_flux_veci

#     # mass_flux_sum = sum(mass_flux) * diffusion_coeff[1] 
#     # mass_flux_sum = -sum(mass_flux) * diffusion_coeff[iscal] 
#     #TODO convention normal

#     # printstyled(color=:green, @sprintf "\n mass_flux : %.2e %.2e %.2e\n" max(abs.(mass_flux)...) mass_flux_sum sum(veci(mass_flux,grid,2))*diffusion_coeff[1])
#     return mass_flux, mass_flux_vec1, mass_flux_vecb, mass_flux_veci
# end


# """
#   Compute boundary mass flux like in Khalighi 2023
# """
# function compute_b_mass_flux!(num,grid, grid_u, grid_v, phL, phS,  opC_pL, opC_pS,diffusion_coeff)
    
#     @unpack nLS = num
#     #Liquid phase
#     @unpack trans_scalD = phL
#     scalD= trans_scalD[:,:,1] #H2 species
#     opC_p = opC_pL

#     mass_flux = fnzeros(grid,num)

#     vecb(mass_flux,grid) = ( opC_p.HxT_b * opC_p.iMx_b * opC_p.Hx_b + opC_p.HyT_b * opC_p.iMy_b * opC_p.Hy_b ) * vecb(scalD,grid) * diffusion_coeff[1] 

#     mass_flux_sum = sum(mass_flux) 

#     # printstyled(color=:green, @sprintf "\n mass_flux : %.2e \n" max(abs.(mass_flux)...))
    
#     return mass_flux_sum
# end

"""
    BC_LS_test!(num, cl, grid, A, B, rhs, BC)

Prints the contact angle, from BC_LS!
"""
function BC_LS_test!(grid, u, A, B, rhs, BC)
    @unpack x, y, nx, ny, dx, dy, ind = grid
    @unpack all_indices, b_left, b_bottom, b_right, b_top = ind
    @unpack left, bottom, right, top = BC

    π2 = π / 2.0

    boundaries_idx = [b_left[1], b_bottom[1], b_right[1], b_top[1]]

    left2 = vcat(all_indices[2,2], all_indices[2:end-1,2], all_indices[end-1,2])
    bottom2 = vcat(all_indices[2,2], all_indices[2,2:end-1], all_indices[2,end-1])
    right2 = vcat(all_indices[2,end-1], all_indices[2:end-1,end-1], all_indices[end-1,end-1])
    top2 = vcat(all_indices[end-1,2], all_indices[end-1,2:end-1], all_indices[end-1,end-1])
    boundaries2 = [left2, bottom2, right2, top2]

    left3 = vcat(all_indices[3,3], all_indices[2:end-1,3], all_indices[end-2,3])
    bottom3 = vcat(all_indices[3,3], all_indices[3,2:end-1], all_indices[3,end-2])
    right3 = vcat(all_indices[3,end-2], all_indices[2:end-1,end-2], all_indices[end-2,end-2])
    top3 = vcat(all_indices[end-2,3], all_indices[end-2,2:end-1], all_indices[end-2,end-2])
    boundaries3 = [left3, bottom3, right3, top3]

    boundaries_t = [left, bottom, right, top]

    direction = [y, x, y, x]

    for (i, (idx, idx2, idx3, xy)) in enumerate(zip(boundaries_idx, boundaries2, boundaries3, direction))
        pks, _ = findminima(abs.(u[idx]))
        # if is_neumann(boundaries_t[i])

        #     for (II, JJ) in zip(idx, idx2)
        #         pII = lexicographic(II, grid.ny)
        #         pJJ = lexicographic(JJ, grid.ny)

        #         A[pII,:] .= 0.0
        #         A[pII,pII] = 1.0
        #         A[pII,pJJ] = -1.0
        #         B[pII,:] .= 0.0
        #     end
        if is_neumann_cl(boundaries_t[i]) && maximum(u[idx]) > 0.0 && minimum(u[idx]) < 0.0 && length(pks) >= 2
            pks1 = idx[pks[1]]
            pkse = idx[pks[end]]

            # Gradually update the contact angle
            Δθe = 90.0 * π / 180

        

            # Find current contact angle
            dist = sqrt((x[idx2[pks[1]]] - x[pks1])^2 + (y[idx2[pks[1]]] - y[pks1])^2)
            old = u[pks1] - u[idx2[pks[1]]]
            # Levelset difference between two consecutive points might be bigger
            # than the distance between them if it's not reinitialized often enough
            if abs(old) > dist
                old = sign(old) * dist
            end
            θe_old = acos(old / dist)

            # printstyled(color=:green, @sprintf "\n θe_old : %.2e \n" θe_old)
            printstyled(color=:green, @sprintf "\n θe_old : %.2e °\n" θe_old*180.0/π)


            return x[idx2[pks]],y[idx2[pks]],x[pks1],y[pks1]

            # # Compute new contact angle
            # if abs(boundaries_t[i].θe - θe_old) > Δθe
            #     θe = θe_old + sign(boundaries_t[i].θe - θe_old) * Δθe
            # else
            #     θe = boundaries_t[i].θe
            # end

            # # distance between the center of the drop and the contact line
            # d = abs(xy[pks1] + u[pks1] - (xy[pkse] + u[pkse])) / 2.0

            # for (II, JJ) in zip(idx[2:end-1], idx2[2:end-1])
            #     pII = lexicographic(II, grid.ny)
            #     pJJ = lexicographic(JJ, grid.ny)

            #     A[pII,:] .= 0.0
            #     A[pII,pII] = 1.0
            #     A[pII,pJJ] = -1.0
            #     B[pII,:] .= 0.0

            #     # Compute levelset angle at a distance u[II] from the contact line
            #     if θe < π2
            #         newθ = atan(tan(θe) * (1.0 - u[II] / d))
            #     else
            #         newθ = π - atan(tan(π - θe) * (1.0 - u[II] / d))
            #     end

            #     rhs[pII] = dist * cos(newθ)
            # end
        # elseif is_neumann_cl(boundaries_t[i]) || is_neumann_inh(boundaries_t[i])
        #     for (II, JJ, KK) in zip(idx, idx2, idx3)
        #         pII = lexicographic(II, grid.ny)
        #         pJJ = lexicographic(JJ, grid.ny)

        #         A[pII,:] .= 0.0
        #         A[pII,pII] = 1.0
        #         A[pII,pJJ] = -1.0
        #         B[pII,:] .= 0.0

        #         rhs[pII] = u[JJ] - u[KK]
        #     end
        end
    end

    return nothing
end

function kill_dead_cells_val!(T::Vector, grid, geo,val)
    @unpack ny, ind = grid

    @inbounds @threads for II in ind.all_indices
        pII = lexicographic(II, ny)
     
        if geo.cap[II,5] < 1e-12
            # print("\n replace ",pII," ", geo.cap[II,5], " ",T[pII], " ",val)
            T[pII] = val
        end
    end
end

function kill_dead_cells_val!(S::SubArray{T,N,P,I,L}, grid, geo,val) where {T,N,P<:Vector{T},I,L}
    @unpack ny, ind = grid

    @inbounds @threads for II in ind.all_indices
        pII = lexicographic(II, ny)
        if geo.cap[II,5] < 1e-12
            # print("\n replace ",pII," ", geo.cap[II,5], " ",S[pII], " ",val)
            S[pII] = val
        end
    end
end

function kill_dead_cells_val!(S::SubArray{T,N,P,I,L}, grid, geo,val) where {T,N,P<:Array{T,3},I,L}
    @unpack ind = grid
    # print("kill dead cells mat")
    @inbounds @threads for II in ind.all_indices
        if geo.cap[II,5] < 1e-12
            # print(II, S[II])
            S[II] = val
            # print("v2",S[II])
        end
    end
end



function kill_dead_cells_val_wall!(S::SubArray{T,N,P,I,L}, grid, geo,val) where {T,N,P<:Vector{T},I,L}
    @unpack ny, ind = grid
    # In ind.b_right, the first vector is the boundary and the second one are the cells right next to them
    @inbounds @threads for II in ind.b_left[1]
        pII = lexicographic(II, ny)
        if geo.cap[II,5] < 1e-12
            print("\n replace ",pII," ", geo.cap[II,5], " ",S[pII], " ",val)
            S[pII] = val
        end
    end

    @inbounds @threads for II in ind.b_bottom[1]
        pII = lexicographic(II, ny)
        if geo.cap[II,5] < 1e-12
            print("\n replace ",pII," ", geo.cap[II,5], " ",S[pII], " ",val)
            S[pII] = val
        end
    end

    @inbounds @threads for II in ind.b_right[1]
        pII = lexicographic(II, ny)
        if geo.cap[II,5] < 1e-12
            print("\n replace ",pII," ", geo.cap[II,5], " ",S[pII], " ",val)
            S[pII] = val
        end
    end

    @inbounds @threads for II in ind.b_top[1]
        pII = lexicographic(II, ny)
        if geo.cap[II,5] < 1e-12
            print("\n replace ",pII," ", geo.cap[II,5], " ",S[pII], " ",val)
            S[pII] = val
        end
    end

end

# function kill_dead_bc_left_wall!(BC,grid,iLS,val)
#     for i = 1:grid.ny
#         II = CartesianIndex(i,1)
#         if grid.LS[iLS].geoL.cap[II,5] < 1e-12
#             BC.left.val[i] = val
#         end
#     end
# end

