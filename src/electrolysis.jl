#From set_heat in heat_coupled.jl, poisson.jl

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
function scalar_transport!(bc, num, grid, op, geo, ph, concentration0, MIXED, projection,
    op_conv, grid_u, geo_u, grid_v, geo_v,
    periodic_x, periodic_y, convection, ls_advection, BC_int, diffusion_coeff, convection_Cdivu,A,B,all_CUTCT)
    @unpack τ, aniso, nb_transported_scalars = num
    @unpack nx, ny, dx, dy, ind, LS  = grid
    @unpack all_indices, inside, b_left, b_bottom, b_right, b_top = ind
    @unpack Bx, By, BxT, ByT, Hx, Hy, HxT, HyT, M, iMx, iMy, χ = op
    @unpack CT, CUTCT = op_conv
    @unpack u, v, uD, vD = ph

    printstyled(color=:red, @sprintf "\n levelset: start scalar_transport!\n")
    print("\n nb_transported_scalars ",nb_transported_scalars)
    println("\n grid.LS[1].geoL.dcap[1,1,:]",grid.LS[1].geoL.dcap[1,1,:])

  



    ni = nx * ny
    nb = 2 * nx + 2 * ny
    nt = 2 * ni + nb

    # A = spzeros(nt, nt)
    # B = spzeros(nt, nt)

    # all_CUTCT = zeros(grid.ny * grid.nx, nb_transported_scalars)

    A .= 0.0
    B .= 0.0
    all_CUTCT .=0.0

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

        # # Mass matrices
        # M.diag .= vec(geo.dcap[:,:,5])
        # Mx = zeros(ny,nx+1)
        # for II in ind.all_indices
        #     Mx[II] = geo.dcap[II,8]
        #     inv_weight_diag(num,II,II,ny,iMx,Mx)
        # end
        # for II in ind.b_right[1]
        #     Mx[δx⁺(II)] = geo.dcap[II,10]
        #     inv_weight_diag(num,II,δx⁺(II),ny,iMx,Mx)
        # end
        # My = zeros(ny+1,nx)
        # for II in ind.all_indices
        #     My[II] = geo.dcap[II,9]
        #     inv_weight_diag(num,II,II,ny+1,iMy,My)
        # end
        # for II in ind.b_top[1]
        #     My[δy⁺(II)] = geo.dcap[II,11]
        #     inv_weight_diag(num,II,δy⁺(II),ny+1,iMy,My)
        # end

        # Mass matrices
        M.diag .= vec(geo.dcap[:,:,5])
        Mx = zeros(ny,nx+1)
        for II in ind.all_indices
            Mx[II] = geo.dcap[II,8]
            # inv_weight_diag(num,II,II,ny,iMx,Mx)
        end
        for II in ind.b_right[1]
            Mx[δx⁺(II)] = geo.dcap[II,10]
            # inv_weight_diag(num,II,δx⁺(II),ny,iMx,Mx)
        end
        My = zeros(ny+1,nx)
        for II in ind.all_indices
            My[II] = geo.dcap[II,9]
            # inv_weight_diag(num,II,II,ny+1,iMy,My)
        end
        for II in ind.b_top[1]
            My[δy⁺(II)] = geo.dcap[II,11]
            # inv_weight_diag(num,II,δy⁺(II),ny+1,iMy,My)
        end
     
        # iMx.diag .= 1. ./ (vec(Mx) .+ eps(0.01))
        # iMy.diag .= 1. ./ (vec(My) .+ eps(0.01))
        iMx = Diagonal(inv_weight_eps2.(num.epsilon_mode,num.epsilon_vol,vec(Mx)))
        iMy = Diagonal(inv_weight_eps2.(num.epsilon_mode,num.epsilon_vol,vec(My)))



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
        mass_matrix_borders!(num,ind, op.iMx_b, op.iMy_b, op.iMx_bd, op.iMy_bd, geo.dcap, ny)
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

        # print("\n BC ", bc[iscal].left.val)
        # print("\n ")
        # print("\n a0_b ", a0_b)

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
        A[ni+1:2*ni,end-nb+1:end] = b * (HxT[1] * op.iMx_b * op.Hx_b .+ HyT[1] * op.iMy_b * op.Hy_b)

        # Border BCs
        A[end-nb+1:end,1:ni] = b_b * (op.HxT_b * op.iMx_b' * Bx .+ op.HyT_b * op.iMy_b' * By)
        A[end-nb+1:end,ni+1:2*ni] = b_b * (op.HxT_b * op.iMx_b' * op.Hx[1] .+ op.HyT_b * op.iMy_b' * op.Hy[1])
        A[end-nb+1:end,end-nb+1:end] = pad(b_b * (op.HxT_b * op.iMx_bd * op.Hx_b .+ op.HyT_b * op.iMy_bd * op.Hy_b) .- op.χ_b * a1_b, 4.0)
        # A[end-nb+1:end,end-nb+1:end] = pad(b_b * (op.HxT_b * op.iMx_bd * op.Hx_b .+ op.HyT_b * op.iMy_bd * op.Hy_b) .- op.χ_b * a1_b)


        # printstyled(color=:red, @sprintf "\n test dummy \n")
        # testb = 68
        # testn = ny-testb+1
        # A[end-nb+testn,testn] = -2
        # A[end-nb+testn,ni+testn] = 0
        # A[end-nb+testn,end-nb+testn] = 2

        # Explicit part of heat equation
        B[1:ni,1:ni] = M .+ 0.5 .* τ .* diffusion_coeff_scal .* LT .- τ .* CT
        B[1:ni,ni+1:2*ni] = 0.5 .* τ .* diffusion_coeff_scal .* LD
        B[1:ni,end-nb+1:end] = 0.5 .* τ .* diffusion_coeff_scal .* LD_b

        # for testn in 1:ny
        #     print("\n jmp ", testn," j ",ny-testn+1," chi_b",op.χ_b[end-nb+testn,end-nb+testn])
        # end
        # ##################################################################################################

        # #TODO test
        # print("testBCwall")
        # # Implicit part of heat equation
        # A[1:ni,1:ni] = pad_crank_nicolson(M .- 0.5 .* τ .* diffusion_coeff_scal .* LT, grid, τ)
        # A[1:ni,ni+1:2*ni] = - 0.5 .* τ .* diffusion_coeff_scal .* LD
        # A[1:ni,end-nb+1:end] = - 0.5 .* τ .* diffusion_coeff_scal .* LD_b

        # # Interior BC
        # A[ni+1:2*ni,1:ni] = b * (HxT[1] * iMx * Bx .+ HyT[1] * iMy * By)
        # A[ni+1:2*ni,ni+1:2*ni] = pad((b * χ[1]+ b_b * op.χ_b)/( χ[1]+ op.χ_b)  * (HxT[1] * iMx * Hx[1] .+ HyT[1] * iMy * Hy[1]) .- χ[1] * a1 /( χ[1]+ op.χ_b) - op.χ_b * a1_b /( χ[1]+ op.χ_b))
        # A[ni+1:2*ni,end-nb+1:end] = ((b * χ[1]+ b_b * op.χ_b)/( χ[1]+ op.χ_b))* (HxT[1] * op.iMx_b * op.Hx_b .+ HyT[1] * op.iMy_b * op.Hy_b)

        # # Border BCs
        # A[end-nb+1:end,1:ni] = b_b * (op.HxT_b * op.iMx_b' * Bx .+ op.HyT_b * op.iMy_b' * By)
        # A[end-nb+1:end,ni+1:2*ni] = b_b * (op.HxT_b * op.iMx_b' * Hx[1] .+ op.HyT_b * op.iMy_b' * op.Hy[1])
        # A[end-nb+1:end,end-nb+1:end] = pad(b_b * (op.HxT_b * op.iMx_bd * op.Hx_b .+ op.HyT_b * op.iMy_bd * op.Hy_b) .- op.χ_b * a1_b, 4.0)

        # # Explicit part of heat equation
        # B[1:ni,1:ni] = M .+ 0.5 .* τ .* diffusion_coeff_scal .* LT .- τ .* CT
        # B[1:ni,ni+1:2*ni] = 0.5 .* τ .* diffusion_coeff_scal .* LD
        # B[1:ni,end-nb+1:end] = 0.5 .* τ .* diffusion_coeff_scal .* LD_b


        # ##################################################################################################



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


        if convection_Cdivu>0
            # Duv = fzeros(grid)
            # Duv = fnzeros(grid,num)

            #TODO BC missing

            Duv = op.AxT * vec1(ph.uD,grid_u) .+ op.Gx_b * vecb(ph.uD,grid_u) .+
            op.AyT * vec1(ph.vD,grid_v) .+ op.Gy_b * vecb(ph.vD,grid_v)
            for iLS in 1:num.nLS
                if !is_navier(BC_int[iLS]) && !is_navier_cl(BC_int[iLS])
                    Duv .+= op.Gx[iLS] * veci(ph.uD,grid_u,iLS+1) .+ 
                            op.Gy[iLS] * veci(ph.vD,grid_v,iLS+1)
                end
            end
        
            # rhs .+= Duv #.* ph.trans_scalD[:,iscal] multiplied just after

            # vec1(rhs_ϕ,grid) .= rho1 .* iτ .* Duv #TODO
            vec1(rhs,grid) .+= Duv 

            printstyled(color=:green, @sprintf "\n max Duv for C.div(u): %.2e\n" maximum(Duv))


        end
        ####################################################################################################
        
        @views mul!(rhs, B, ph.trans_scalD[:,iscal], 1.0, 1.0) #TODO @views not necessary ?


        ####################################################################################################      
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
        ####################################################################################################

        @views ph.trans_scalD[:,iscal] .= A \ rhs
        
        #printstyled(color=:magenta, @sprintf "\n Condition number " cond(A,2))
        printstyled(color=:magenta, @sprintf "\n Condition number " cond(Array(A), 2))
        

        ####################################################################################################

        #if iscal == 1

            # print("\n test left after A/r L", vecb_T(ph.trans_scalD[:,iscal], grid))
            # print("\n end ", ph.trans_scal[end,:,iscal])

        #     debug_border_top(grid.nx,grid.ny,A,B,rhs,geo,inside,ind,grid,bc,iscal,num,op,ni,nb,ph,Bx,By)
        #     debug_border_top(grid.nx-1,grid.ny,A,B,rhs,geo,inside,ind,grid,bc,iscal,num,op,ni,nb,ph,Bx,By)
        #     debug_border_top(grid.nx-2,grid.ny,A,B,rhs,geo,inside,ind,grid,bc,iscal,num,op,ni,nb,ph,Bx,By)
        
        #     rhs2 = A * ph.trans_scalD[:,iscal] - rhs

        #     # print("\n test left after A/r L", vecb_T(ph.trans_scalD[:,iscal], grid))
        #     # print("\n rhs2 ", maximum(rhs2))


        #end

        ####################################################################################################
        # if iscal == 1 

        #     # rhs2 = A * ph.trans_scalD[:,iscal] - rhs

            
        #     # print("\n test left after A/r L", vecb_L(ph.trans_scalD[:,iscal], grid))
        #     # for testn in 1:ny
        #     #     printstyled(color=:green, @sprintf "\n jtmp : %.5i j : %.5i chi_b %.2e border %.5e\n" testn ny-testn+1 op.χ_b[end-nb+testn,end-nb+testn] vecb_L(ph.trans_scalD[:,iscal], grid)[testn])
        #     # end

        #     debug_border_left(1,68,A,B,rhs,geo,inside,ind,grid,bc,iscal,num,op,ni,nb,ph,Bx,By)
        #     debug_border_left(1,69,A,B,rhs,geo,inside,ind,grid,bc,iscal,num,op,ni,nb,ph,Bx,By)
        #     debug_border_left(1,70,A,B,rhs,geo,inside,ind,grid,bc,iscal,num,op,ni,nb,ph,Bx,By)


        # end

        # print("\n test left after A/r T", vecb_T(ph.trans_scalD[:,iscal], grid))

        # @views ph.trans_scal[:,:,iscal] .= reshape(veci(ph.trans_scalD[:,iscal],grid,2), grid)


        # printstyled(color=:green, @sprintf "\n average c %s\n" average!(ph.trans_scal[:,:,iscal], grid, LS[1].geoL, num))
        ####################################################################################################


        @views ph.trans_scal[:,:,iscal] .= reshape(veci(ph.trans_scalD[:,iscal],grid,1), grid)


        # printstyled(color=:cyan, @sprintf "\n after resol\n")
        # print("\n",ph.trans_scal[1,:,iscal])
        # print("\n",reshape(veci(rhs,grid,1), grid)[1,:])
        # print("\n",reshape(veci(rhs,grid,1), grid)[2,:])
        # print("\n",reshape(veci(rhs,grid,1), grid)[3,:]) 

        nonzero = veci(ph.trans_scalD[:,iscal],grid,2)[abs.(veci(ph.trans_scalD[:,iscal],grid,2)) .> 0.0]
        # print("nonzero\n")
        # print(nonzero)
        printstyled(color=:green, @sprintf "\n mean  interface : %.2e\n" mean(nonzero))

        print("\n iscal",iscal," BC ",bc[iscal].left.val)

        print("\n test left", vecb_L(ph.trans_scalD[:,iscal], grid))

        print("\n iscal", iscal, "\n ")


        if iscal!=3 #H2O consummed at the electrode, would need to make distinction to make sure the decrease in H2O is physical or not
          
          
            @views kill_dead_cells_val!(ph.trans_scal[:,:,iscal], grid, LS[1].geoL,concentration0[iscal]) 

            # @views veci(ph.trans_scalD[:,iscal],grid,1) .= vec(ph.trans_scal[:,:,iscal])


            # ####################################################################################################
            # #check top 
            # ####################################################################################################

            # not needed, with Neumann 0 at top and non null at left it cannot stay equal to 0.16 the ref value

            # jplot = grid.ny
            # for iplot in 1:grid.nx
            #     II = CartesianIndex(jplot, iplot) #(id_y, id_x)
            #     pII = lexicographic(II, grid.ny)

            #     # printstyled(color=:green, @sprintf "\n i: %5i %.2e %.2e \n" iplot abs(reverse(vecb_T(ph.trans_scalD[:,iscal],grid))[iplot]-concentration0[iscal])/concentration0[iscal] reverse(vecb_T(ph.trans_scalD[:,iscal],grid))[iplot])


            #     if (abs(vecb_T(ph.trans_scalD[:,iscal],grid)[iplot]-concentration0[iscal])/concentration0[iscal]>num.concentration_check_factor)

            #     # if (ph.trans_scalD[pII,iscal] >concentration0[iscal]) 
                    
            #         # && (ph.trans_scalD[pII,iscal] == maximum(ph.trans_scalD[:,iscal]) )
            #         printstyled(color=:green, @sprintf "\n i: %5i j: %5i %.2e %.2e %.2e %.2e \n" iplot jplot grid.x[iplot]/num.plot_xscale grid.y[jplot]/num.plot_xscale ph.trans_scalD[pII,iscal] rhs[pII])
                    
            #         @views scalar_debug!(dir, CT, all_CUTCT[:,iscal], u, v, bcTx, bcTy, bcU, bcV, geo.dcap, ny, 
            #         bc[iscal], inside, b_left[1], b_bottom[1], b_right[1], b_top[1],num,grid,iplot,jplot)

            #         printstyled(color=:red, @sprintf "\n exit 1 error \n" )
            #         @error("exit error")

            #         return

            #     end
            # end

            min_border_bulk = min(minimum(ph.trans_scal[:,:,iscal]),minimum(vecb(ph.trans_scal[:,:,iscal],grid)))

            if min_border_bulk.<concentration0[iscal]*(1-num.concentration_check_factor)
                print("iscal ",iscal)
                printstyled(color=:red, @sprintf "\n concentration: %.10e %.10e \n" min_border_bulk concentration0[iscal]*(1-num.concentration_check_factor))
                # printstyled(color=:red, @sprintf "\n concentration drop: %.2e%% \n" (minimum(ph.trans_scal[:,:,iscal])-concentration0[iscal])/concentration0[iscal]*100)
                @error("concentration too low")
                # return current_i

                ####################################################################################################

                ####################################################################################################
                
                for jplot in 1:grid.ny
                    for iplot in 1:grid.nx
                        II = CartesianIndex(jplot, iplot) #(id_y, id_x)
                        pII = lexicographic(II, grid.ny)

                        if (ph.trans_scalD[pII,iscal] <concentration0[iscal]*(1-num.concentration_check_factor))

                        # if (ph.trans_scalD[pII,iscal] >concentration0[iscal]) 
                            
                            # && (ph.trans_scalD[pII,iscal] == maximum(ph.trans_scalD[:,iscal]) )
                            printstyled(color=:green, @sprintf "\n i: %5i j: %5i %.2e %.2e %.2e %.2e \n" iplot jplot grid.x[iplot]/num.plot_xscale grid.y[jplot]/num.plot_xscale ph.trans_scalD[pII,iscal] rhs[pII])
                            
                            if num.debug== "scalar_debug"
                                @views scalar_debug!(dir, CT, all_CUTCT[:,iscal], u, v, bcTx, bcTy, bcU, bcV, geo.dcap, ny, 
                                bc[iscal], inside, b_left[1], b_bottom[1], b_right[1], b_top[1],num,grid,iplot,jplot)
                            end
                            printstyled(color=:red, @sprintf "\n exit 1 error \n" )
                            @error("TODO exit error")
                            #return 
                            #TODO return and redo iteration,

                        end
                    end
                end
                #############################################################################################################
            end #too low

            printstyled(color=:red, @sprintf "\n concentration variation vs min: %.2e%% \n" (minimum(ph.trans_scal[:,:,iscal])-concentration0[iscal])/concentration0[iscal]*100)

        else 
            @views kill_dead_cells_val!(ph.trans_scal[:,:,iscal], grid, LS[1].geoL,concentration0[iscal]) 

            max_border_bulk = max(maximum(ph.trans_scal[:,:,iscal]),maximum(vecb(ph.trans_scal[:,:,iscal],grid)))
            if max_border_bulk.> concentration0[iscal]*(1+num.concentration_check_factor)
                print("iscal ",iscal)
                printstyled(color=:red, @sprintf "\n concentration: %.10e %.10e \n" max_border_bulk concentration0[iscal]*(1-num.concentration_check_factor))
                # printstyled(color=:red, @sprintf "\n concentration increase: %.2e%% \n" (maximum(ph.trans_scal[:,:,iscal])-concentration0[iscal])/concentration0[iscal]*100)
                @error("concentration too high")
                # return current_i

                #############################################################################################################
                for jplot in 1:grid.ny
                    for iplot in 1:grid.nx
                        II = CartesianIndex(jplot, iplot) #(id_y, id_x)
                        pII = lexicographic(II, grid.ny)

                        if (ph.trans_scalD[pII,iscal] >concentration0[iscal]*(1+num.concentration_check_factor))

                        # if (ph.trans_scalD[pII,iscal] >concentration0[iscal]) 
                            
                            # && (ph.trans_scalD[pII,iscal] == maximum(ph.trans_scalD[:,iscal]) )
                            printstyled(color=:green, @sprintf "\n i: %5i j: %5i %.2e %.2e %.2e %.2e \n" iplot jplot grid.x[iplot]/num.plot_xscale grid.y[jplot]/num.plot_xscale ph.trans_scalD[pII,iscal] rhs[pII])
                            if num.debug== "scalar_debug"
                                @views scalar_debug!(dir, CT, all_CUTCT[:,iscal], u, v, bcTx, bcTy, bcU, bcV, geo.dcap, ny, 
                                bc[iscal], inside, b_left[1], b_bottom[1], b_right[1], b_top[1],num,grid,iplot,jplot)
                            end
                        end
                    end
                end
                #############################################################################################################

            end #too high
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
Butler-Volmer model, relation between exchande current and potential at electrode, with effects of concentration
c0_H2: reference concentration
"""
function butler_volmer_concentration(alpha_a,alpha_c,c_H2,c0_H2,c_H2O,c0_H2O,c_KOH,c0_KOH,Faraday,i0,phi_ele,phi_ele1,Ru,temperature0)
    eta = phi_ele1 - phi_ele
    i_current = i0*(sqrt(c_H2/c0_H2)*(c_KOH/c0_KOH)*exp(alpha_a*Faraday*eta/(Ru*temperature0))-(c_H2O/c0_H2O)*exp(-alpha_c*Faraday*eta/(Ru*temperature0)))
    return i_current
end

"""
Butler-Volmer model, relation between exchande current and potential at electrode without concentration
"""
function butler_volmer_no_concentration(alpha_a,alpha_c,Faraday,i0,phi_ele,phi_ele1,Ru,temperature0)
    eta = phi_ele1 - phi_ele
    i_current = i0*(exp(alpha_a*Faraday*eta/(Ru*temperature0))-exp(-alpha_c*Faraday*eta/(Ru*temperature0)))
    return i_current
end



"""
From Stefan_velocity!
"""
function electrolysis_velocity!(num, grid, LS, V, TL, MIXED, periodic_x, periodic_y, concentration_scal_intfc,electrolysis_phase_change_case, varflux)
    @unpack geoS, geoL, κ = LS

    V .= 0

    if electrolysis_phase_change_case == "levelset"

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

    else

        interfacial_area = 0.0
        @inbounds for II in MIXED
            V[II] = sum(varflux)

            #χx = (geo.dcap[:,:,3] .- geo.dcap[:,:,1]) .^ 2
            #χy = (geo.dcap[:,:,4] .- geo.dcap[:,:,2]) .^ 2
            #χ[iLS].diag .= sqrt.(vec(χx .+ χy))

            χx = (geo.dcap[II,3] .- geo.dcap[II,1]) .^ 2
            χy = (geo.dcap[II,4] .- geo.dcap[II,2]) .^ 2

            interfacial_area += sqrt.(vec(χx .+ χy))
        end 

        V./= interfacial_area

        print("\n phase-change velocity ", sum(varflux)/interfacial_area)

        printstyled(color=:magenta, @sprintf "\n phase-change velocity %.2e area %.2e πR %.2e\n" sum(varflux)/interfacial_area interfacial_area π*num.radius)


        magenta

    end

    return nothing
end

"""
From update_free_surface_velocity and update_stefan_velocity
"""
function update_free_surface_velocity_electrolysis(num, grid, grid_u, grid_v, iLS, uD, vD, periodic_x, periodic_y, Vmean, concentration_scalD, diffusion_coeff_scal,concentration_scal_intfc, electrolysis_phase_change_case)
    @unpack MWH2,rho1,rho2=num
    # @unpack χ,=opC
    #TODO check sign 

    printstyled(color=:green, @sprintf "\n grid p u v max : %.2e %.2e %.2e\n" maximum(abs.(grid.V[grid.LS[iLS].MIXED])) maximum(abs.(grid_u.V[grid.LS[iLS].MIXED])) maximum(abs.(grid_v.V[grid_v.LS[iLS].MIXED])))

    #plot_electrolysis_velocity!(num, grid, grid.LS[iLS], grid.V, concentration_scalD, grid.LS[iLS].MIXED, periodic_x, periodic_y, concentration_scal_intfc)

    electrolysis_velocity!(num, grid, grid.LS[iLS], grid.V, concentration_scalD, grid.LS[iLS].MIXED, periodic_x, periodic_y, concentration_scal_intfc,electrolysis_phase_change_case, varflux)
    

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


"""Interpolate velocity on scalar grid for regular grids for vizualisation

# Arguments
- `gp::Mesh`: scalar grid
"""
function interpolate_grid_liquid!(gp,gu,gv,u,v,us,vs)
    
    
    # us = p .*0
    # vs = p .*0
    LS_u =gu.LS[1]
    LS_v = gv.LS[1]
    us .= (
        (u[:,2:end] .* LS_u.geoL.dcap[:,2:end,6] .+ 
        u[:,1:end-1] .* LS_u.geoL.dcap[:,1:end-1,6]) ./ 
        (LS_u.geoL.dcap[:,1:end-1,6] .+ LS_u.geoL.dcap[:,2:end,6] .+ 1e-8 )
    )
    vs .= (
        (v[2:end,:] .* LS_v.geoL.dcap[2:end,:,7] .+ 
        v[1:end-1,:] .* LS_v.geoL.dcap[1:end-1,:,7]) ./
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

    minuL = minimum(phL.u)
    maxuL = maximum(phL.u)
    moyuL = mean(phL.u)

    minvL = minimum(phL.v)
    maxvL = maximum(phL.v)
    moyvL = mean(phL.v)

    minpL = minimum(phL.p)
    maxpL = maximum(phL.p)
    moypL = mean(phL.p)

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

    print("$(@sprintf("min(u) %.6e", minuL))\t$(@sprintf("min(v) %.6e", minvL))\t$(@sprintf("min(p) %.6e", minpL))\n")
    print("$(@sprintf("max(u) %.6e", maxuL))\t$(@sprintf("max(v) %.6e", maxvL))\t$(@sprintf("max(p) %.6e", maxpL))\n")
    print("$(@sprintf("moy(u) %.6e", moyuL))\t$(@sprintf("moy(v) %.6e", moyvL))\t$(@sprintf("moy(p) %.6e", moypL))\n")

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
function compute_grad_p!(num,grid, grid_u, grid_v, pD, opC_p,opC_u,opC_v)
    
    @unpack nLS = num

    #Liquid phase
    # @unpack pD = phL

    ∇ϕ_x = opC_p.iMx * opC_p.Bx * vec1(pD,grid) .+ opC_p.iMx_b * opC_p.Hx_b * vecb(pD,grid)
    ∇ϕ_y = opC_p.iMy * opC_p.By * vec1(pD,grid) .+ opC_p.iMy_b * opC_p.Hy_b * vecb(pD,grid)

    for iLS in 1:nLS
        ∇ϕ_x .+= opC_p.iMx * opC_p.Hx[iLS] * veci(pD,grid,iLS+1)
        ∇ϕ_y .+= opC_p.iMy * opC_p.Hy[iLS] * veci(pD,grid,iLS+1)
    end

    grd_x = reshape(veci(∇ϕ_x,grid_u,1), grid_u)
    grd_y = reshape(veci(∇ϕ_y,grid_v,1), grid_v)

    printstyled(color=:red, @sprintf "\n grad min max x %.2e %.2e y %.2e %.2e\n" minimum(grd_x) maximum(grd_x) minimum(grd_y) maximum(grd_y))


    ∇ϕ_x = opC_u.AxT * opC_u.Rx * vec1(pD,grid) .+ opC_u.Gx_b * vecb(pD,grid)
    ∇ϕ_y = opC_v.AyT * opC_v.Ry * vec1(pD,grid) .+ opC_v.Gy_b * vecb(pD,grid)
    for iLS in 1:nLS
        ∇ϕ_x .+= opC_u.Gx[iLS] * veci(pD,grid,iLS+1)
        ∇ϕ_y .+= opC_v.Gy[iLS] * veci(pD,grid,iLS+1)
    end

    iMu = Diagonal(inv_weight_eps2.(num.epsilon_mode,num.epsilon_vol,opC_u.M.diag))
    iMv = Diagonal(inv_weight_eps2.(num.epsilon_mode,num.epsilon_vol,opC_v.M.diag))
    ∇ϕ_x = iMu * ∇ϕ_x
    ∇ϕ_y = iMv * ∇ϕ_y
    
    grd_x = reshape(veci(∇ϕ_x,grid_u,1), grid_u)
    grd_y = reshape(veci(∇ϕ_y,grid_v,1), grid_v)

    printstyled(color=:magenta, @sprintf "\n grad min max x %.2e %.2e y %.2e %.2e\n" minimum(grd_x) maximum(grd_x) minimum(grd_y) maximum(grd_y))


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
function compute_mass_flux!(num, grid, grid_u, grid_v, phL, phS, opC_pL, opC_pS, diffusion_coeff, iscal, geo)
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


    # iplot = 1
    # for jplot in 1:grid.ny
    #     #for iplot in 1:grid.nx
    #     II = CartesianIndex(jplot, iplot) #(id_y, id_x)
    #     pII = lexicographic(II, grid.ny)

    #     if mass_flux_vec1_2[II]>0
    #         printstyled(color=:green, @sprintf "\n j %.5i m %.2e HxT %.2e\n" jplot mass_flux_vec1_2[II] opC_p.HxT[iLStmp][II])
    #         printstyled(color=:red, @sprintf "\n iMx %.10e iMy %.10e \n" opC_p.iMy.diag[pII] opC_p.iMy.diag[pII] )
    #         print("\n B ", II," ",opC_p.Bx[pII,pII]," ",opC_p.BxT[pII,pII])
    #     end

    # end


    # iplot=1
    # jplot = 59
    # II = CartesianIndex(jplot, iplot) #(id_y, id_x)
    # pII = lexicographic(II, grid.ny)
    # printstyled(color=:magenta, @sprintf "\n j %.5i m %.2e HxT %.2e HyT %.2e Bx %.2e By %.2e iMx %.10e iMy %.10e\n" jplot mass_flux_vec1_2[II] opC_p.HxT[iLStmp][II] opC_p.HyT[iLStmp][II] opC_p.Bx[pII,pII] opC_p.By[pII,pII] opC_p.iMy.diag[pII] opC_p.iMy.diag[pII])

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

    if num.io_pdi>0
        printstyled(color=:magenta, @sprintf "\n PDI write_mass_flux %.5i \n" num.current_i)

        ######################################################################################
        #TODO careful : nstep needs to be updated beforehand
        @ccall "libpdi".PDI_multi_expose("write_mass_flux"::Cstring,
        # "nstep"::Cstring, nstep::Ref{Clonglong}, PDI_OUT::Cint,
        # "time"::Cstring, time::Ref{Cdouble}, PDI_OUT::Cint,
        "mass_flux"::Cstring, mass_flux_vec1_2::Ptr{Cdouble}, PDI_OUT::Cint,
        "mass_flux_bulk"::Cstring, mass_flux_vec1_2::Ptr{Cdouble}, PDI_OUT::Cint,
        "mass_flux_border"::Cstring, mass_flux_vecb_2::Ptr{Cdouble}, PDI_OUT::Cint,
        "mass_flux_intfc"::Cstring, mass_flux_veci_2::Ptr{Cdouble}, PDI_OUT::Cint,
        # "Bx"::Cstring, opC_p.Bx::Ptr{Cdouble}, PDI_OUT::Cint,
        # "By"::Cstring, opC_p.By::Ptr{Cdouble}, PDI_OUT::Cint,
        # "HxT"::Cstring, opC_p.HxT::Ptr{Cdouble}, PDI_OUT::Cint,
        # "HyT"::Cstring, opC_p.HyT::Ptr{Cdouble}, PDI_OUT::Cint,
        C_NULL::Ptr{Cvoid})::Cvoid

        # @ccall "libpdi".PDI_multi_expose("write_mass_flux"::Cstring,
        # "nstep"::Cstring, nstep::Ref{Clonglong}, PDI_OUT::Cint,
        # "time"::Cstring, time::Ref{Cdouble}, PDI_OUT::Cint,
        # "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
        # "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
        # "levelset_p"::Cstring, LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
        # "levelset_u"::Cstring, grid_u.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
        # "levelset_v"::Cstring, grid_v.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
        # "trans_scal_1D"::Cstring, phL.trans_scalD::Ptr{Cdouble}, PDI_OUT::Cint,
        # "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
        # "i_current_x"::Cstring, Eus::Ptr{Cdouble}, PDI_OUT::Cint,   
        # "i_current_y"::Cstring, Evs::Ptr{Cdouble}, PDI_OUT::Cint,   
        # "velocity_x"::Cstring, us::Ptr{Cdouble}, PDI_OUT::Cint,   
        # "velocity_y"::Cstring, vs::Ptr{Cdouble}, PDI_OUT::Cint,      
        # "radius"::Cstring, current_radius::Ref{Cdouble}, PDI_OUT::Cint, 
        # C_NULL::Ptr{Cvoid})::Cvoid
    end #if num.io_pdi>0

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


function test_LS(grid)
    iLS = 1
    minLS = minimum(grid.LS[iLS].u)
    maxLS = maximum(grid.LS[iLS].u)

    printstyled(color=:red, @sprintf "\n levelset: min %.2e max %.2e\n" minLS maxLS)
    if maximum(grid.LS[iLS].u) < 0.0
        @error("LS initialization")
    end

    if minLS*maxLS > 0.0
        @error("No interface")
    end

end


function FE_set_momentum_debug(
    bc_type, num, grid, opC,
    A, B,
    L, bc_L, bc_L_b, Mm1, BC,
    ls_advection
    )
    @unpack τ = num
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

    if ls_advection
        # Implicit part of viscous term
        A[1:ni,1:ni] = pad_crank_nicolson(M .- τ .* L, grid, τ)
        # Contribution to implicit part of viscous term from outer boundaries
        A[1:ni,end-nb+1:end] = - τ .* bc_L_b
        # Boundary conditions for outer boundaries
        A[end-nb+1:end,1:ni] = b_b * (HxT_b * iMx_b' * Bx .+ HyT_b * iMy_b' * By)
        A[end-nb+1:end,end-nb+1:end] = pad(b_b * (HxT_b * iMx_bd * Hx_b .+ HyT_b * iMy_bd * Hy_b) .- χ_b * a1_b)

        B[1:ni,1:ni] = Mm1
    end

    for iLS in 1:num.nLS
        if is_dirichlet(bc_type[iLS])
            vel = copy(grid.V)
            __a1 = -1.0
            __b = 0.0
        elseif is_neumann(bc_type[iLS])
            vel = 0.0
            __a1 = 0.0
            __b = 1.0
        elseif is_robin(bc_type[iLS])
            vel = 0.0
            __a1 = -1.0
            __b = 1.0
        elseif is_fs(bc_type[iLS])
            vel = 0.0
            __a1 = 0.0
            __b = 1.0
        elseif is_wall_no_slip(bc_type[iLS])
            vel = bc_type[iLS].val
            __a1 = -1.0
            __b = 0.0
                else
            vel = bc_type[iLS].val
            __a1 = -1.0
            __b = 0.0
        end

        a0 = ones(grid) .* vel
        _a1 = ones(grid) .* __a1
        a1 = Diagonal(vec(_a1))
        _b = ones(grid) .* __b
        b = Diagonal(vec(_b))

        sb = iLS*ni+1:(iLS+1)*ni

        if ls_advection
            # Contribution to implicit part of viscous term from inner boundaries
            A[1:ni,sb] = - τ .* bc_L[iLS]
            # Boundary conditions for inner boundaries
            A[sb,1:ni] = b * (HxT[iLS] * iMx * Bx .+ HyT[iLS] * iMy * By)
            # Contribution to Neumann BC from other boundaries
            for i in 1:num.nLS
                if i != iLS
                    A[sb,i*ni+1:(i+1)*ni] = b * (HxT[iLS] * iMx * Hx[i] .+ HyT[iLS] * iMy * Hy[i])
                end
            end
            A[sb,sb] = pad(b * (HxT[iLS] * iMx * Hx[iLS] .+ HyT[iLS] * iMy * Hy[iLS]) .- χ[iLS] * a1)
            A[sb,end-nb+1:end] = b * (HxT[iLS] * iMx_b * Hx_b .+ HyT[iLS] * iMy_b * Hy_b)
            # Boundary conditions for outer boundaries
            A[end-nb+1:end,sb] = b_b * (HxT_b * iMx_b' * Hx[iLS] .+ HyT_b * iMy_b' * Hy[iLS])
        end

        veci(rhs,grid,iLS+1) .= χ[iLS] * vec(a0)
    end

    vecb(rhs,grid) .= χ_b * vec(a0_b)
    
    return rhs
end

# function h5write2(filename, name::AbstractString, data,status; pv...)
    
    # @static if @isdefined(:HDF5)

#         file = h5open(filename, status; pv...)
#         try
#             write(file, name, data)
#         finally
#             close(file)
#         end

#     else
        
#         printstyled(color=:red, @sprintf "\n HDF5 not loaded, not writing with h5write2:\n")


#     end

# end


function pressure_projection_debug!(
    time_scheme, bc_int,
    num, grid, geo, grid_u, geo_u, grid_v, geo_v, ph,
    BC_u, BC_v, BC_p,
    opC_p, opC_u, opC_v, op_conv,
    Au, Bu, Av, Bv, Aϕ, Auv, Buv,
    Lpm1, bc_Lpm1, bc_Lpm1_b, Lum1, bc_Lum1, bc_Lum1_b, Lvm1, bc_Lvm1, bc_Lvm1_b,
    Cum1, Cvm1, Mum1, Mvm1,
    periodic_x, periodic_y, advection, ls_advection, current_i, Ra, navier, pres_free_suface,jump_mass_flux,mass_flux
    )
    @unpack Re, τ, σ, g, β, nLS, nNavier = num
    @unpack p, pD, ϕ, ϕD, u, v, ucorrD, vcorrD, uD, vD, ucorr, vcorr, uT = ph
    @unpack Cu, Cv, CUTCu, CUTCv = op_conv
    @unpack rho1,rho2,visc_coeff = num


    # iRe = 1.0 / Re
    iRe = visc_coeff
    iτ = 1.0 / τ
    irho1 = 1.0/rho1

    printstyled(color=:green, @sprintf "\n rho1 : %.2e irho1 : %.2e iRe : %.2e\n" rho1 irho1 iRe)

    ∇ϕ_x = opC_u.AxT * opC_u.Rx * vec1(pD,grid) .+ opC_u.Gx_b * vecb(pD,grid)
    ∇ϕ_y = opC_v.AyT * opC_v.Ry * vec1(pD,grid) .+ opC_v.Gy_b * vecb(pD,grid)
    for iLS in 1:nLS
        ∇ϕ_x .+= opC_u.Gx[iLS] * veci(pD,grid,iLS+1)
        ∇ϕ_y .+= opC_v.Gy[iLS] * veci(pD,grid,iLS+1)
    end

    grd_x = reshape(veci(∇ϕ_x,grid_u,1), grid_u)
    grd_y = reshape(veci(∇ϕ_y,grid_v,1), grid_v)

    printstyled(color=:cyan, @sprintf "\n B %.2e T %.2e L %.2e R %.2e\n" maximum(abs.(vecb_B(∇ϕ_x,grid_u))) maximum(abs.(vecb_T(∇ϕ_x,grid_u))) maximum(abs.(vecb_L(∇ϕ_y,grid_v))) maximum(abs.(vecb_R(∇ϕ_y,grid_v))))
    printstyled(color=:red, @sprintf "\n B %.2e T %.2e L %.2e R %.2e\n" maximum(abs.(grd_x[end,:])) maximum(abs.(grd_x[1,:])) maximum(abs.(grd_y[:,1])) maximum(abs.(grd_y[:,end])))

    printstyled(color=:red, @sprintf "\n p %.2e p %.2e p %.2e p %.2e \n" maximum(abs.(p[end,:])) maximum(abs.(p[1,:])) maximum(abs.(p[:,end])) maximum(abs.(p[:,1])))

    ptest = reshape(veci(pD,grid,1), grid)
    printstyled(color=:magenta, @sprintf "\n ptest %.2e p %.2e p %.2e p %.2e \n" maximum(abs.(ptest[end,:])) maximum(abs.(ptest[1,:])) maximum(abs.(ptest[:,end])) maximum(abs.(ptest[:,1])))



    ∇ϕ_x .= 0.0
    ∇ϕ_y .= 0.0

    compute_grad_p!(num,grid, grid_u, grid_v, pD, opC_p, opC_u, opC_v)


    if num.prediction == 1 || num.prediction == 2
        printstyled(color=:red, @sprintf "\n pressure_in_prediction \n")

        #TODO
        compute_grad_p!(num,grid, grid_u, grid_v, pD, opC_p, opC_u, opC_v)


        ∇ϕ_x = opC_u.AxT * opC_u.Rx * vec1(pD,grid) .+ opC_u.Gx_b * vecb(pD,grid)
        ∇ϕ_y = opC_v.AyT * opC_v.Ry * vec1(pD,grid) .+ opC_v.Gy_b * vecb(pD,grid)
        for iLS in 1:nLS
            ∇ϕ_x .+= opC_u.Gx[iLS] * veci(pD,grid,iLS+1)
            ∇ϕ_y .+= opC_v.Gy[iLS] * veci(pD,grid,iLS+1)
        end

        grd_x = reshape(veci(∇ϕ_x,grid_u,1), grid_u)
        grd_y = reshape(veci(∇ϕ_y,grid_v,1), grid_v)

        printstyled(color=:cyan, @sprintf "\n B %.2e T %.2e L %.2e R %.2e\n" maximum(abs.(vecb_B(∇ϕ_x,grid_u))) maximum(abs.(vecb_T(∇ϕ_x,grid_u))) maximum(abs.(vecb_L(∇ϕ_y,grid_v))) maximum(abs.(vecb_R(∇ϕ_y,grid_v))))


        grd_xfull = opC_p.iMx * opC_p.Bx * vec1(pD,grid) .+ opC_p.iMx_b * opC_p.Hx_b * vecb(pD,grid)
        grd_yfull = opC_p.iMy * opC_p.By * vec1(pD,grid) .+ opC_p.iMy_b * opC_p.Hy_b * vecb(pD,grid)

        for iLS in 1:num.nLS
            grd_xfull .+= opC_p.iMx * opC_p.Hx[iLS] * veci(pD,grid,iLS+1)
            grd_yfull .+= opC_p.iMy * opC_p.Hy[iLS] * veci(pD,grid,iLS+1)
        end

        grd_x = reshape(veci(grd_xfull,grid_u,1), grid_u)
        grd_y = reshape(veci(grd_yfull,grid_v,1), grid_v)

        # printstyled(color=:red, @sprintf "\n grad min max x %.2e %.2e y %.2e %.2e\n" minimum(grd_x) maximum(grd_x) minimum(grd_y) maximum(grd_y))

        print("\n dt ", num.τ)

        # ph.Gxm1 .+= ∇ϕ_x
        # ph.Gym1 .+= ∇ϕ_y

        ph.Gxm1 .= 0.0
        ph.Gym1 .= 0.0

        # ph.Gxm1 .= grd_xfull
        # ph.Gym1 .= grd_yfull

        ph.Gxm1 .= ∇ϕ_x
        ph.Gym1 .= ∇ϕ_y

        # printstyled(color=:red, @sprintf "\n grad max x %.2e y %.2e %.2e\n" maximum(grd_xfull) maximum(grd_yfull) maximum(ph.Gym1)*irho1)
        printstyled(color=:red, @sprintf "\n full grad min max x %.2e %.2e y %.2e %.2e\n" minimum(grd_xfull) maximum(grd_xfull) minimum(grd_yfull) maximum(grd_yfull))

        printstyled(color=:red, @sprintf "\n full grad min max x %.2e %.2e y %.2e %.2e\n" minimum(ph.Gxm1) maximum(ph.Gxm1) minimum(ph.Gym1) maximum(ph.Gym1))

        printstyled(color=:red, @sprintf "\n full grad min max x %.2e %.2e y %.2e %.2e\n" minimum(τ.*irho1.*ph.Gxm1) maximum(τ.*irho1.*ph.Gxm1) minimum(τ.*irho1.*ph.Gym1) maximum(τ.*irho1.*ph.Gym1))

        printstyled(color=:red, @sprintf "\n full grad min max x %.2e %.2e y %.2e %.2e\n" minimum(∇ϕ_x) maximum(∇ϕ_x) minimum(∇ϕ_y) maximum(∇ϕ_x))


        ∇ϕ_x .= 0.0
        ∇ϕ_y .= 0.0
        

    end

    nip = grid.nx * grid.ny

    niu = grid_u.nx * grid_u.ny
    nbu = 2 * grid_u.nx + 2 * grid_u.ny
    ntu = (nLS - nNavier + 1) * niu + nbu

    niv = grid_v.nx * grid_v.ny
    nbv = 2 * grid_v.nx + 2 * grid_v.ny
    ntv = (nLS - nNavier + 1) * niv + nbv

    if is_FE(time_scheme)
        rhs_u, rhs_v, rhs_ϕ, rhs_uv, Lp, bc_Lp, bc_Lp_b, Lu, bc_Lu, bc_Lu_b, Lv, bc_Lv, bc_Lv_b = set_FE!(
            bc_int, num, grid, geo, grid_u, geo_u, grid_v, geo_v,
            opC_p, opC_u, opC_v, BC_p, BC_u, BC_v,
            Au, Bu, Av, Bv, Aϕ, Auv, Buv,
            Lpm1, bc_Lpm1, bc_Lpm1_b, Lum1, bc_Lum1, bc_Lum1_b, Lvm1, bc_Lvm1, bc_Lvm1_b,
            Mum1, Mvm1, iRe, op_conv, ph,
            periodic_x, periodic_y, advection, ls_advection, navier
        )
    elseif is_CN(time_scheme)
        rhs_u, rhs_v, rhs_ϕ, Lp, bc_Lp, bc_Lp_b, Lu, bc_Lu, bc_Lu_b, Lv, bc_Lv, bc_Lv_b = set_CN!(
            bc_int, num, grid, geo, grid_u, geo_u, grid_v, geo_v,
            opC_p, opC_u, opC_v, BC_p, BC_u, BC_v,
            Au, Bu, Av, Bv, Aϕ,
            Lpm1, bc_Lpm1, bc_Lpm1_b, Lum1, bc_Lum1, bc_Lum1_b, Lvm1, bc_Lvm1, bc_Lvm1_b,
            Mum1, Mvm1, iRe, op_conv, ph,
            periodic_x, periodic_y, advection, ls_advection
        )
    end

    ra_x = Ra .* sin(β) .* opC_u.M * vec(hcat(zeros(grid_u.ny), ph.T))
    ra_y = Ra .* cos(β) .* opC_v.M * vec(vcat(zeros(1,grid_v.nx), ph.T))

    grav_x = g .* sin(β) .* opC_u.M * fones(grid_u)
    grav_y = g .* cos(β) .* opC_v.M * fones(grid_v)

    Convu = fzeros(grid_u)
    Convv = fzeros(grid_v)
    Cui = Cu * vec(u) .+ CUTCu
    Cvi = Cv * vec(v) .+ CUTCv
    if advection
        if current_i == 1
            Convu .+= Cui
            Convv .+= Cvi
        else
            Convu .+= 1.5 .* Cui .- 0.5 .* Cum1 #Cui returned at the end of function to Cum1
            Convv .+= 1.5 .* Cvi .- 0.5 .* Cvm1
        end
    end

    

    # printstyled(color=:green, @sprintf "\n max abs(Cu) : %.2e u: %.2e CUTCu: %.2e \n" maximum(abs.(Cu)) maximum(abs.(u)) maximum(abs.(CUTCu)))

    # u and v are coupled if a Navier slip BC is employed inside, otherwise they are uncoupled
    if !navier
        # if is_wall_no_slip(bc_int)
        #     vec1(uD,grid_u) .= vec(u)
        #     # update_dirichlet_field!(grid_u, uD, u, BC_u)
        #     vec1(rhs_u,grid_u) .+= -τ .* (opC_u.AxT * opC_u.Rx * vec1(pD,grid) .+ opC_u.Gx_b * vecb(pD,grid))
        #     for iLS in 1:nLS
        #         vec1(rhs_u,grid_u) .+= -τ .* (opC_u.Gx[iLS] * veci(pD,grid,iLS+1))
        #     end
        # end
        mul!(rhs_u, Bu, uD, 1.0, 1.0)
        vec1(rhs_u,grid_u) .+= τ .* grav_x
        vec1(rhs_u,grid_u) .-= τ .* Convu
        vec1(rhs_u,grid_u) .+= τ .* ra_x
        printstyled(color=:green, @sprintf "\n rhs u : %.2e uD %.2e Bu %.2e M %.2e \n" maximum(abs.(rhs_u)) maximum(abs.(uD)) maximum(abs.(Bu)) maximum(abs.(Mum1)))

        vec1(rhs_u,grid_u) .-= τ .* irho1 .* ph.Gxm1 
        
        printstyled(color=:green, @sprintf "\n rhs u : %.2e \n" maximum(abs.(rhs_u)))

        kill_dead_cells!(vec1(rhs_u,grid_u), grid_u, geo_u[end])
        for iLS in 1:nLS
            kill_dead_cells!(veci(rhs_u,grid_u,iLS+1), grid_u, geo_u[end])
        end
        # @time bicgstabl!(ucorrD, Au, rhs_u, log=true)
        try
            # @time bicgstabl!(ucorrD, Au, rhs_u, Pl=Diagonal(Au), log=true)
            @time ucorrD .= Au \ rhs_u
        catch e
            ucorrD .= Inf
            println(e)
        end

        # printstyled(color=:green, @sprintf "\n max abs(ucorrD) : %.2e uD: %.2e \n" maximum(abs.(ucorrD)) maximum(abs.(uD)))

        kill_dead_cells!(vec1(ucorrD,grid_u), grid_u, geo_u[end])
        for iLS in 1:nLS
            kill_dead_cells!(veci(ucorrD,grid_u,iLS+1), grid_u, geo_u[end])
        end
        ucorr .= reshape(vec1(ucorrD,grid_u), grid_u)

        # if is_wall_no_slip(bc_int)
        #     vec1(vD,grid_v) .= vec(v)
        #     # update_dirichlet_field!(grid_v, vD, v, BC_v)
        #     vec1(rhs_v,grid_v) .+= -τ .* (opC_v.AyT * opC_v.Ry * vec1(pD,grid) .+opC_v.Gy_b * vecb(pD,grid))
        #     for iLS in 1:nLS
        #         vec1(rhs_v,grid_v) .+= -τ .* (opC_v.Gy[iLS] * veci(pD,grid,iLS+1))
        #     end
        # end
        mul!(rhs_v, Bv, vD, 1.0, 1.0)

        test1 = vec1(rhs_v,grid_v)[1,1]/Poiseuille_fmax(grid_v.x[1,1],num.v_inlet,num.L0)
        test2 = test1 / (grid_v.dx[1,1]^2/2)
        printstyled(color=:red, @sprintf "\n rhs_v vec1 %.10e /pois %.10e /pois %.10e\n" vec1(rhs_v,grid_v)[1,1] test1 test2)

        vec1(rhs_v,grid_v) .+= - τ .* grav_y
        vec1(rhs_v,grid_v) .-= τ .* Convv
        vec1(rhs_v,grid_v) .+= τ .* ra_y

        test1 = vec1(rhs_v,grid_v)[1,1]/Poiseuille_fmax(grid_v.x[1,1],num.v_inlet,num.L0)
        test2 = test1 / (grid_v.dx[1,1]^2/2)
        test3 = vec1(rhs_v,grid_v)[1,1]-Poiseuille_fmax(grid_v.x[1,1],num.v_inlet,num.L0)*(grid_v.dx[1,1]^2/2)
        printstyled(color=:red, @sprintf "\n rhs_v vec1 %.10e /pois %.10e /pois %.10e diff %.10e\n" vec1(rhs_v,grid_v)[1,1] test1 test2 test3)

        # printstyled(color=:green, @sprintf "\n rhs: %.2e vD %.2e \n" maximum(abs.(rhs_v)) maximum(abs.(vD)))
        printstyled(color=:green, @sprintf "\n rhs v : %.2e vD %.2e Bv %.2e M %.2e \n" maximum(abs.(rhs_v)) maximum(abs.(vD)) maximum(abs.(Bv)) maximum(abs.(Mvm1)))


        vec1(rhs_v,grid_v) .-= τ .* irho1 .* ph.Gym1
        printstyled(color=:green, @sprintf "\n rhs: %.2e \n" maximum(abs.(rhs_v)))


        test1 = vec1(rhs_v,grid_v)[1,1]/Poiseuille_fmax(grid_v.x[1,1],num.v_inlet,num.L0)
        test2 = test1 / (grid_v.dx[1,1]^2/2)
        test3 = vec1(rhs_v,grid_v)[1,1]-Poiseuille_fmax(grid_v.x[1,1],num.v_inlet,num.L0)*(grid_v.dx[1,1]^2/2)
        test4 = test3/(τ .* irho1)/ (grid_v.dx[1,1]^2/2)
        printstyled(color=:red, @sprintf "\n rhs_v vec1 %.10e /pois %.10e /pois %.10e diff %.10e diff %.10e\n" vec1(rhs_v,grid_v)[1,1] test1 test2 test3 test4)



        kill_dead_cells!(vec1(rhs_v,grid_v), grid_v, geo_v[end])
        for iLS in 1:nLS
            kill_dead_cells!(veci(rhs_v,grid_v,iLS+1), grid_v, geo_v[end])
        end
        # bicgstabl!(vcorrD, Av, rhs_v, log=true)
        
        
        iplot = 1
        jplot = 1
        II = CartesianIndex(jplot, iplot) #(id_y, id_x)
        # pII = lexicographic(II, grid.ny +1)

        print("\n after kill dead cells ", (grid_v.dx[1,1]^2/2)," full " ,(grid_v.dx[1,1]^2)," test ",geo_v[end].cap[II,5])
        test1 = vec1(rhs_v,grid_v)[1,1]/Poiseuille_fmax(grid_v.x[1,1],num.v_inlet,num.L0)
        test2 = test1 / (grid_v.dx[1,1]^2/2)
        test3 = vec1(rhs_v,grid_v)[1,1]-Poiseuille_fmax(grid_v.x[1,1],num.v_inlet,num.L0)*(grid_v.dx[1,1]^2/2)
        test4 = test3/(τ .* irho1)/ (grid_v.dx[1,1]^2/2)
        printstyled(color=:red, @sprintf "\n rhs_v vec1 %.10e /pois %.10e /pois %.10e diff %.10e diff %.10e\n" vec1(rhs_v,grid_v)[1,1] test1 test2 test3 test4)


        try
            # @time bicgstabl!(vcorrD, Av, rhs_v, Pl=Diagonal(Av), log=true)
            @time vcorrD .= Av \ rhs_v
        catch e
            vcorrD .= Inf
            println(e)
        end

        printstyled(color=:yellow, @sprintf "\n vcorrD \n")

        # iplot = 64
        # jplot = 64
        # II = CartesianIndex(jplot, iplot) #(id_y, id_x)
        # pII = lexicographic(II, grid.ny +1)
        
        # # test = τ .* iRe.*Lv *vcorrD
        # # test =
        # # bc_Lv, bc_Lv_b
        # # print("\n testvisc ",Lv)
        # print("\n ")
        # # print("\n testvisc ",Lv[jplot,iplot])
        # print("\n testvisc ", II," ",Lv[pII,:])
        # printstyled(color=:green, @sprintf "\n Bx: %.10e \n" opC_v.Bx[pII,pII])
        # printstyled(color=:green, @sprintf "\n BxT: %.10e \n" opC_v.BxT[pII,pII])
        # printstyled(color=:green, @sprintf "\n iMx: %.10e \n" opC_v.iMx[pII,pII])
        # printstyled(color=:green, @sprintf "\n Mx: %.10e iMx: %.10e iMx: %.10e\n" geo_v[end].dcap[II,8] 1/geo_v[end].dcap[II,8] 1/(geo_v[end].dcap[II,8]+eps(0.01)))

        

        # @unpack Bx, By, Hx, Hy, HxT, HyT, χ, M, iMx, iMy, Hx_b, Hy_b, HxT_b, HyT_b, iMx_b, iMy_b, iMx_bd, iMy_bd, χ_b = opC
        @unpack  M = opC_v
        print("\n M min ",minimum(M), " max ", maximum(M))


        ni = grid_v.nx * grid_v.ny
        nb = 2 * grid_v.nx + 2 * grid_v.ny
        nt = (num.nLS + 1) * ni + nb

        Avtest = spzeros(nt, nt)
       
        # Implicit part of viscous term
        Avtest[1:ni,1:ni] = iRe .*Lv #pad_crank_nicolson(Lv, grid, τ)
        # Contribution to implicit part of viscous term from outer boundaries
        Avtest[1:ni,end-nb+1:end] = iRe .* bc_Lv_b

        vecv = reshape(vec1(vD,grid_v),grid_v)


        iplot = 1
        jplot = 64
        II = CartesianIndex(jplot, iplot) #(id_y, id_x)
        pII = lexicographic(II, grid.ny +1)
        print("\n ")
        print("\n testvisc ", II," ",Lv[pII,:])
        print("\n testvisc ", II," ",Avtest[pII,:])

        print("\n testvisc ", II," ",bc_Lv_b[pII,:])


        testAv = Avtest * vcorrD .*rho1 
        testAv2 = Avtest * vD .*rho1 

        printstyled(color=:green, @sprintf "\n Avtest * vcorrD/My : %.10e exact %.10e 4/3exact %.10e\n" testAv[pII]*opC_p.iMy.diag[pII] testAv2[pII]*opC_p.iMy.diag[pII] testAv2[pII]*opC_p.iMy.diag[pII]*4/3)

        print("\n op ", rho1*opC_p.iMy.diag[pII]*iRe*(-5*vecv[64,1] +1*vecv[64,2]))
        print("\n op ",vecv[64,1]," op ",vecv[64,2])
        print("\n op ",opC_p.iMy.diag[pII])
        print("\n iRe ", iRe)

        print("\n op ", rho1*iRe)

        print("\n op ", rho1*iRe*(-5*vecv[64,1] +1*vecv[64,2]))

        print("\n testvisc ", II," ",Avtest[pII,pII]," ",Avtest[pII,pII]*opC_p.iMy.diag[pII]," ",Avtest[pII,pII]*opC_p.iMy.diag[pII]*rho1, " ",Avtest[pII,pII]*opC_p.iMy.diag[pII]*rho1*vecv[64,1])




        ####################################################################################################        
        iplot = 2
        jplot = 64
        II = CartesianIndex(jplot, iplot) #(id_y, id_x)
        pII = lexicographic(II, grid.ny +1)
        print("\n ")
        print("\n testvisc ", II," ",Lv[pII,:])

        printstyled(color=:green, @sprintf "\n Avtest * vcorrD/My : %.10e exact %.10e\n" testAv[pII]*opC_p.iMy.diag[pII] testAv2[pII]*opC_p.iMy.diag[pII])
        ####################################################################################################

        ####################################################################################################        
        iplot = 1
        jplot = 1
        II = CartesianIndex(jplot, iplot) #(id_y, id_x)
        pII = lexicographic(II, grid.ny +1)
        print("\n ")
        print("\n testvisc ", II," ",Lv[pII,:])
        print("\n testvisc ", II," ",Avtest[pII,:])

        printstyled(color=:green, @sprintf "\n Avtest * vcorrD/My : %.10e exact %.10e 4/3exact %.10e\n" testAv[pII]*opC_p.iMy.diag[pII] testAv2[pII]*opC_p.iMy.diag[pII] testAv2[pII]*opC_p.iMy.diag[pII]*4/3)
        ####################################################################################################

        #not 
        # printstyled(color=:green, @sprintf "\n Avtest * vcorrD : %.10e Avtest * vcorrD/M : %.10e Avtest * vcorrD/My : %.10e\n" testAv[pII] testAv[pII]*opC_v.iMx_bd[pII,pII] testAv[pII]*opC_v.iMy[pII,pII])
        # ####################################################################################################        
        # iplot = 2
        # jplot = 64
        # II = CartesianIndex(jplot, iplot) #(id_y, id_x)
        # pII = lexicographic(II, grid.ny +1)
        # print("\n ")
        # print("\n testvisc ", II," ",Lv[pII,:])
        # printstyled(color=:green, @sprintf "\n Avtest * vcorrD : %.10e Avtest * vcorrD/M : %.10e Avtest * vcorrD/My : %.10e \n" testAv[pII] testAv[pII]*opC_v.iMx[pII,pII] testAv[pII]*opC_v.iMy[pII,pII])
        # ####################################################################################################


        # 6.103515625000243e-13

        # testAv = Av * vcorrD - 


        # testLv = fnzeros(grid, num)
        # testLv = fnzeros(grid, num)
        # mul!(testLv, Lv, vcorrD, 1.0, 1.0)
        # print("\n testvisc ", II," ",testLv[pII,:])

        # mul!(rhs_v, Bv, vD, 1.0, 1.0)


        # printstyled(color=:green, @sprintf "\n Lv: %.10e \n" Lv[pII,pII])
        # printstyled(color=:green, @sprintf "\n Bx: %.10e \n" opC_v.Bx[pII,pII])
        # printstyled(color=:green, @sprintf "\n BxT: %.10e \n" opC_v.BxT[pII,pII])
        # printstyled(color=:green, @sprintf "\n iMx: %.10e \n" opC_v.iMx[pII,pII])
        # printstyled(color=:green, @sprintf "\n Mx: %.10e iMx: %.10e iMx: %.10e\n" geo_v[end].dcap[II,8] 1/geo_v[end].dcap[II,8] 1/(geo_v[end].dcap[II,8]+eps(0.01)))




        iplot = 2
        jplot = 64
        II = CartesianIndex(jplot, iplot) #(id_y, id_x)
        pII = lexicographic(II, grid.ny +1)
        print("\n ")
        print("\n testvisc ", II," ",Lv[pII,:])
        
        printstyled(color=:green, @sprintf "\n Avtest * vcorrD/My : %.10e exact %.10e\n" testAv[pII]*opC_v.iMy[pII,pII] testAv2[pII]*opC_v.iMy[pII,pII])

        # print("\n testvisc ", II," ",testLv[pII,:])
        # printstyled(color=:green, @sprintf "\n Lv: %.10e \n" Lv[pII,pII])
        # printstyled(color=:green, @sprintf "\n Bx: %.10e \n" opC_v.Bx[pII,pII])
        # printstyled(color=:green, @sprintf "\n BxT: %.10e \n" opC_v.BxT[pII,pII])
        # printstyled(color=:green, @sprintf "\n iMx: %.10e \n" opC_v.iMx[pII,pII])
        # printstyled(color=:green, @sprintf "\n Mx: %.10e iMx: %.10e iMx: %.10e\n" geo_v[end].dcap[II,8] 1/geo_v[end].dcap[II,8] 1/(geo_v[end].dcap[II,8]+eps(0.01)))

        # iplot = 3
        # jplot = 64
        # II = CartesianIndex(jplot, iplot) #(id_y, id_x)
        # pII = lexicographic(II, grid.ny +1)
        # print("\n ")
        # print("\n testvisc ", II," ",Lv[pII,:])
        # printstyled(color=:green, @sprintf "\n Lv: %.10e \n" Lv[pII,pII])
        # printstyled(color=:green, @sprintf "\n Bx: %.10e \n" opC_v.Bx[pII,pII])
        # printstyled(color=:green, @sprintf "\n BxT: %.10e \n" opC_v.BxT[pII,pII])
        # printstyled(color=:green, @sprintf "\n iMx: %.10e \n" opC_v.iMx[pII,pII])
        # printstyled(color=:green, @sprintf "\n Mx: %.10e iMx: %.10e iMx: %.10e\n" geo_v[end].dcap[II,8] 1/geo_v[end].dcap[II,8] 1/(geo_v[end].dcap[II,8]+eps(0.01)))

        # iplot = 4
        # jplot = 64
        # II = CartesianIndex(jplot, iplot) #(id_y, id_x)
        # pII = lexicographic(II, grid.ny +1)
        # print("\n ")
        # print("\n testvisc ", II," ",Lv[pII,:])
        # printstyled(color=:green, @sprintf "\n Lv: %.10e \n" Lv[pII,pII])
        # printstyled(color=:green, @sprintf "\n Bx: %.10e \n" opC_v.Bx[pII,pII])
        # printstyled(color=:green, @sprintf "\n BxT: %.10e \n" opC_v.BxT[pII,pII])
        # printstyled(color=:green, @sprintf "\n iMx: %.10e \n" opC_v.iMx[pII,pII])
        # printstyled(color=:green, @sprintf "\n Mx: %.10e iMx: %.10e iMx: %.10e\n" geo_v[end].dcap[II,8] 1/geo_v[end].dcap[II,8] 1/(geo_v[end].dcap[II,8]+eps(0.01)))


        # ny = grid.ny
    
        # testb = jplot
        # testn = ny-testb+1
        # print("\n test",testn," testb ",testb)
        # # printstyled(color=:green, @sprintf "\n jtmp : %.5i j : %.5i chi_b %.2e  chi_b adim %.2e border %.2e\n" testn testb op.χ_b[end-nb+testn,end-nb+testn] op.χ_b[end-nb+testn,end-nb+testn]/grid.dy[1,1] vecb_L(ph.trans_scalD[:,iscal], grid)[testn])
        # # printstyled(color=:cyan, @sprintf "\n BC %.5e rhs %.5e rhs %.5e \n" bc[iscal].left.val[testn] bc[iscal].left.val[testn]*op.χ_b[end-nb+testn,end-nb+testn] vecb_L(rhs, grid)[testn])
        # # print("\n B ", maximum(B[testb,:])," \n ")
    
        # print("\n A[end-nb+testn,1:ni]", Av[end-nb+testn,1:ni], "\n")
        # print("\n A[end-nb+testn,ni+1:2*ni]", Av[end-nb+testn,ni+1:2*ni], "\n")
        # print("\n A[end-nb+testn,end-nb+1:end]", Av[end-nb+testn,end-nb+1:end], "\n")


        iplot = 1
        jplot = 1
        II = CartesianIndex(jplot, iplot) #(id_y, id_x)
        pII = lexicographic(II, grid.ny +1)
        print("\n ")
        print("\n testvisc ", II," ",Lv[pII,:])
        printstyled(color=:red, @sprintf "\n iMy %.10e %.10e %.10e\n" opC_p.iMy.diag[pII] 1/grid_v.dx[1,1]^2 grid_v.dx[1,1]^2)
        print("\n B ", II," ",opC_p.Bx[pII,pII]," ",opC_p.BxT[pII,pII])

        iplot = 1
        jplot = 64
        II = CartesianIndex(jplot, iplot) #(id_y, id_x)
        pII = lexicographic(II, grid.ny +1)
        print("\n ")
        print("\n testvisc ", II," ",Lv[pII,:])
        printstyled(color=:red, @sprintf "\n iMy %.10e %.10e %.10e\n" opC_p.iMy.diag[pII] 1/grid_v.dx[1,1]^2 grid_v.dx[1,1]^2)
        print("\n B ", II," ",opC_p.Bx[pII,pII]," ",opC_p.BxT[pII,pII])



        iplot = 2
        jplot = 64
        II = CartesianIndex(jplot, iplot) #(id_y, id_x)
        pII = lexicographic(II, grid.ny +1)
        print("\n ")
        print("\n B ", II," ",opC_p.Bx[pII,pII]," ",opC_p.BxT[pII,pII])
       
        # mul!(tmp_x, iMx, Bx)
        # L = BxT * tmp_x
        # mul!(tmp_y, iMy, By)
        # L = L .+ ByT * tmp_y


        # iplot = 1
        # jplot = 1
        # II = CartesianIndex(jplot, iplot) #(id_y, id_x)
        # pII = lexicographic(II, grid_v.ny +1)
        # print("\n ")
        # print("\n testvisc ", II," ",Lv[pII,:])
        # printstyled(color=:red, @sprintf "\n iMy %.10e %.10e %.10e\n" opC_p.iMy.diag[pII] 1/grid_v.dx[1,1]^2 grid_v.dx[1,1]^2)

        # iplot = 1
        # jplot = 64
        # II = CartesianIndex(jplot, iplot) #(id_y, id_x)
        # pII = lexicographic(II, grid_v.ny +1)
        # print("\n ")
        # print("\n testvisc ", II," ",Lv[pII,:])
        # printstyled(color=:red, @sprintf "\n iMy %.10e %.10e %.10e\n" opC_p.iMy.diag[pII] 1/grid_v.dx[1,1]^2 grid_v.dx[1,1]^2)

        

        
        

        #TODO Poiseuille
        test_Poiseuille(num,vcorrD,grid_v)

        printstyled(color=:red, @sprintf "\n vcorrD %.2e %.2e\n" minimum(vcorrD) maximum(vcorrD))

        test_Poiseuille(num,vD,grid_v)


   


        printstyled(color=:red, @sprintf "\n vec1 1\n")
        print(vecv[1,:])

        printstyled(color=:red, @sprintf "\n vecb_B \n" )
        print(vecb_B(vD,grid_v))

        printstyled(color=:red, @sprintf "\n vecb_L vD\n")
        print(vecb_L(vD,grid_v))

        printstyled(color=:red, @sprintf "\n vecb_L vcorrD\n" )
        print(vecb_L(vcorrD,grid_v))


        printstyled(color=:red, @sprintf "\n rhs_v vecb_L \n" )
        print(vecb_L(rhs_v,grid_v))

        kill_dead_cells!(vec1(vcorrD,grid_v), grid_v, geo_v[end])
        for iLS in 1:nLS
            kill_dead_cells!(veci(vcorrD,grid_v,iLS+1), grid_v, geo_v[end])
        end
        vcorr .= reshape(vec1(vcorrD,grid_v), grid_v)
    else
        uvm1 = zeros(ntu + ntv + nNavier * nip)
        uvm1[1:niu] .= vec1(uD,grid_u)
        uvm1[ntu+1:ntu+niv] .= vec1(vD,grid_v)
        uvm1[ntu-nbu+1:ntu] .= vecb(uD,grid_u)
        uvm1[ntu+ntv-nbv+1:ntu+ntv] .= vecb(vD,grid_v)
        _iLS = 1
        for iLS in 1:num.nLS
            if !is_navier(bc_int[iLS]) && !is_navier_cl(bc_int[iLS])
                uvm1[_iLS*niu+1:(_iLS+1)*niu] .= veci(uD,grid_u,iLS+1)
                uvm1[ntu+_iLS*niv+1:ntu+(_iLS+1)*niv] .= veci(vD,grid_v,iLS+1)
                _iLS += 1
            end
        end

        rhs_uv .+=  Buv * uvm1

        rhs_uv[1:niu] .+= τ .* grav_x
        rhs_uv[1:niu] .-= τ .* Convu
        rhs_uv[1:niu] .+= τ .* ra_x
        rhs_uv[1:niu] .-= τ .* irho1 .* ph.Gxm1 

        rhs_uv[ntu+1:ntu+niv] .+= τ .* grav_y
        rhs_uv[ntu+1:ntu+niv] .-= τ .* Convv
        rhs_uv[ntu+1:ntu+niv] .+= τ .* ra_y
        rhs_uv[ntu+1:ntu+niv] .-= τ .* irho1 .* ph.Gym1 

        @views kill_dead_cells!(rhs_uv[1:niu], grid_u, geo_u[end])
        @views kill_dead_cells!(rhs_uv[ntu+1:ntu+niv], grid_v, geo_v[end])
        _iLS = 1
        for iLS in 1:nLS
            sbu = _iLS*niu+1:(_iLS+1)*niu
            sbv = ntu+_iLS*niv+1:ntu+(_iLS+1)*niv
            if !is_navier(bc_int[iLS]) && !is_navier_cl(bc_int[iLS])
                @views kill_dead_cells!(rhs_uv[sbu], grid_u, geo_u[end])
                @views kill_dead_cells!(rhs_uv[sbv], grid_v, geo_v[end])
                _iLS += 1
            end
        end

        uvD = ones(ntu + ntv + nNavier * nip)
        try
            @time uvD .= Auv \ rhs_uv
        catch e
            uvD .= Inf
            println(e)
        end

        vec1(ucorrD, grid_u) .= uvD[1:niu]
        vecb(ucorrD, grid_u) .= uvD[ntu-nbu+1:ntu]
        kill_dead_cells!(vec1(ucorrD,grid_u), grid_u, geo_u[end])
        ucorr .= reshape(vec1(ucorrD,grid_u), grid_u)

        vec1(vcorrD, grid_v) .= uvD[ntu+1:ntu+niv]
        vecb(vcorrD, grid_v) .= uvD[ntu+ntv-nbv+1:ntu+ntv]
        kill_dead_cells!(vec1(vcorrD,grid_v), grid_v, geo_v[end])
        vcorr .= reshape(vec1(vcorrD,grid_v), grid_v)

        nNav = 0
        _iLS = 1
        for iLS in 1:nLS
            if !is_navier(bc_int[iLS]) && !is_navier_cl(bc_int[iLS])
                veci(ucorrD,grid_u,iLS+1) .= uvD[_iLS*niu+1:(_iLS+1)*niu]
                kill_dead_cells!(veci(ucorrD,grid_u,iLS+1), grid_u, geo_u[end])

                veci(vcorrD,grid_v,iLS+1) .= uvD[ntu+_iLS*niv+1:ntu+(_iLS+1)*niv]
                kill_dead_cells!(veci(vcorrD,grid_v,iLS+1), grid_v, geo_v[end])
                _iLS += 1
            else
                @inbounds uT[nNav+1,:] .= vec(uvD[ntu+ntv+1+nNav*nip:ntu+ntv+(nNav+1)*nip])
                nNav += 1
            end
        end
    end

    # printstyled(color=:green, @sprintf "\n max abs(ucorrD) : %.2e vcorrD %.2e \n" maximum(abs.(ucorrD)) maximum(abs.(vcorrD)))


    Duv = opC_p.AxT * vec1(ucorrD,grid_u) .+ opC_p.Gx_b * vecb(ucorrD,grid_u) .+
          opC_p.AyT * vec1(vcorrD,grid_v) .+ opC_p.Gy_b * vecb(vcorrD,grid_v)
    for iLS in 1:nLS
        if !is_navier(bc_int[iLS]) && !is_navier_cl(bc_int[iLS])
            Duv .+= opC_p.Gx[iLS] * veci(ucorrD,grid_u,iLS+1) .+ 
                    opC_p.Gy[iLS] * veci(vcorrD,grid_v,iLS+1)
        end
    end

    #Poisson equation
    # vec1(rhs_ϕ,grid) .= iτ .* Duv
    vec1(rhs_ϕ,grid) .= rho1 .* iτ .* Duv #TODO
    # veci(rhs_ϕ,grid) .*= rho1 #TODO

    # pres_free_suface = 0.0
    #TODO Marangoni
    #TODO phase change
    diff_inv_rho = 1.0/rho1 - 1.0/rho2
    # jump_mass_flux = 0.0 #TODO

    if jump_mass_flux
        for iLS in 1:nLS
            if is_fs(bc_int[iLS])
                Smat = strain_rate(iLS, opC_u, opC_v, opC_p)
                S = Smat[1,1] * vec1(ucorrD,grid_u) .+ Smat[1,2] * veci(ucorrD,grid_u,iLS+1) .+
                    Smat[2,1] * vec1(vcorrD,grid_v) .+ Smat[2,2] * veci(vcorrD,grid_v,iLS+1)
    
                fs_mat = opC_p.HxT[iLS] * opC_p.Hx[iLS] .+ opC_p.HyT[iLS] * opC_p.Hy[iLS]
                veci(rhs_ϕ,grid,iLS+1) .= -2.0 .* iRe .* S .+ Diagonal(diag(fs_mat)) * ( σ .* vec(grid.LS[iLS].κ) .- pres_free_suface .- diff_inv_rho * mass_flux ^ 2)
            end
        end
    else
        for iLS in 1:nLS
            if is_fs(bc_int[iLS])
                Smat = strain_rate(iLS, opC_u, opC_v, opC_p)
                S = Smat[1,1] * vec1(ucorrD,grid_u) .+ Smat[1,2] * veci(ucorrD,grid_u,iLS+1) .+
                    Smat[2,1] * vec1(vcorrD,grid_v) .+ Smat[2,2] * veci(vcorrD,grid_v,iLS+1)

                fs_mat = opC_p.HxT[iLS] * opC_p.Hx[iLS] .+ opC_p.HyT[iLS] * opC_p.Hy[iLS]
                veci(rhs_ϕ,grid,iLS+1) .= -2.0 .* iRe .* S .+ Diagonal(diag(fs_mat)) * ( σ .* vec(grid.LS[iLS].κ) .- pres_free_suface )
            end
        end
    end
    # Remove nullspace by adding small quantity to main diagonal
    if num.null_space == 0
        @inbounds @threads for i in 1:Aϕ.m
            @inbounds Aϕ[i,i] += 1e-10
        end
    end
    kill_dead_cells!(vec1(rhs_ϕ,grid), grid, geo[end])
    for iLS in 1:nLS
        kill_dead_cells!(veci(rhs_ϕ,grid,iLS+1), grid, geo[end])
    end
    # @time bicgstabl!(ϕD, Aϕ, rhs_ϕ, Pl = Diagonal(Aϕ), log = true)

    # rhs_ϕ .*= rho1 #TODO #TODO not BC

    # vec1(rhs_ϕ,grid) .*= rho1 

    # Aϕ .*= irho1

    # vecb(rhs_ϕ,grid) .*= irho1

    @time ϕD .= Aϕ \ rhs_ϕ
    kill_dead_cells!(vec1(ϕD,grid), grid, geo[end])
    for iLS in 1:nLS
        kill_dead_cells!(veci(ϕD,grid,iLS+1), grid, geo[end])
    end
    ϕ .= reshape(vec1(ϕD,grid), grid)

    iMu = Diagonal(inv_weight_eps2.(num.epsilon_mode,num.epsilon_vol,opC_u.M.diag))
    iMv = Diagonal(inv_weight_eps2.(num.epsilon_mode,num.epsilon_vol,opC_v.M.diag))
    # Gradient of pressure, eq. 17 in 
    #"A Conservative Cartesian Cut-Cell Method for Mixed Boundary Conditions and the Incompressible Navier-Stokes Equations on Staggered Meshes"
    ∇ϕ_x = opC_u.AxT * opC_u.Rx * vec(ϕ) .+ opC_u.Gx_b * vecb(ϕD,grid)
    ∇ϕ_y = opC_v.AyT * opC_v.Ry * vec(ϕ) .+ opC_v.Gy_b * vecb(ϕD,grid)
    for iLS in 1:nLS
        ∇ϕ_x .+= opC_u.Gx[iLS] * veci(ϕD,grid,iLS+1)
        ∇ϕ_y .+= opC_v.Gy[iLS] * veci(ϕD,grid,iLS+1)
    end

    # ∇ϕ_x = irho1 .* opC_u.AxT * opC_u.Rx * vec(ϕ) .+ opC_u.Gx_b * vecb(ϕD,grid)
    # ∇ϕ_y = irho1 .* opC_v.AyT * opC_v.Ry * vec(ϕ) .+ opC_v.Gy_b * vecb(ϕD,grid)
    # for iLS in 1:nLS
    #     ∇ϕ_x .+= irho1 .* opC_u.Gx[iLS] * veci(ϕD,grid,iLS+1)
    #     ∇ϕ_y .+= irho1 .* opC_v.Gy[iLS] * veci(ϕD,grid,iLS+1)
    # end

    # if num.prediction == 1 already done
    #     ph.Gxm1 .+= ∇ϕ_x
    #     ph.Gym1 .+= ∇ϕ_y
    # end

    printstyled(color=:magenta, @sprintf "\n full grad min max x %.2e %.2e y %.2e %.2e\n" minimum(∇ϕ_x) maximum(∇ϕ_x) minimum(∇ϕ_y) maximum(∇ϕ_x))


    # iM = Diagonal(1. ./ (vec(geo[end].dcap[:,:,5]) .+ eps(0.01)))

    # iM = Diagonal(inv_weight_eps.(num,geo[end].dcap[:,:,5]))

    iM = Diagonal(inv_weight_eps2.(num.epsilon_mode,num.epsilon_vol,vec(geo[end].dcap[:,:,5])))

    # iM = Diagonal(1. ./ (vec(geo[end].dcap[:,:,5]) ))

    # if is_fs(bc_int)
    if num.prediction == 1
        vec1(pD,grid) .= vec(ϕ .- iRe .* rho1 .* reshape(iM * Duv,grid)) #no τ  since div u not rho1
    elseif num.prediction == 2
        vec1(pD,grid) .+= vec(ϕ .- iRe./2 .* rho1 .* reshape(iM * Duv,grid)) #no τ  since div u not rho1
    else
        vec1(pD,grid) .= vec(ϕ) #.- iRe .* reshape(iM * Duv, grid))
    end
    for iLS in 1:nLS
        veci(pD,grid,iLS+1) .= veci(ϕD,grid,iLS+1)
    end
    vecb(pD,grid) .= vecb(ϕD,grid)
    p .= reshape(vec1(pD,grid), grid)

    #TODO
    compute_grad_p!(num,grid, grid_u, grid_v, pD, opC_p, opC_u, opC_v)


    # else
    #     vec1(pD,grid) .= vec(p) .+ vec(ϕ) #.- iRe .* iM * Duv
    #     vec2(pD,grid) .+= vec2(ϕD,grid)
    #     vecb(pD,grid) .+= vecb(ϕD,grid)
    #     p .= reshape(vec1(pD,grid), grid)
    # end

    # vec1(∇ϕ_x,grid) .*= irho1 
    # vec1(∇ϕ_y,grid) .*= irho1

    
    # u .= ucorr .- τ .* reshape(iMu * ∇ϕ_x, grid_u)
    # v .= vcorr .- τ .* reshape(iMv * ∇ϕ_y, grid_v)

    u .= ucorr .- τ .* irho1 .* reshape(iMu * ∇ϕ_x, grid_u)
    v .= vcorr .- τ .* irho1 .* reshape(iMv * ∇ϕ_y, grid_v)

    kill_dead_cells!(u, grid_u, geo_u[end])
    kill_dead_cells!(v, grid_v, geo_v[end])

    vec1(uD,grid_u) .= vec(u)
    vecb(uD,grid_u) .= vecb(ucorrD,grid_u)
    vec1(vD,grid_v) .= vec(v)
    vecb(vD,grid_v) .= vecb(vcorrD,grid_v)
    for iLS in 1:nLS
        if !is_navier(bc_int[iLS]) && !is_navier_cl(bc_int[iLS])
            veci(uD,grid_u,iLS+1) .= veci(ucorrD,grid_u,iLS+1)
            veci(vD,grid_v,iLS+1) .= veci(vcorrD,grid_v,iLS+1)
        end
        # if is_fs(bc_int[iLS])
        #     @inbounds for II in grid_u.ind.all_indices
        #         pII = lexicographic(II, grid_u.ny)
        #         if abs(veci(ucorrD,grid_u,iLS+1)[pII]) > 1e-12
        #             veci(ucorrD,grid_u,iLS+1)[pII] -= (τ .* iMu * ∇ϕ_x)[pII]
        #         end
        #     end
        #     @inbounds for II in grid_v.ind.all_indices
        #         pII = lexicographic(II, grid_v.ny)
        #         if abs(veci(vcorrD,grid_v,iLS+1)[pII]) > 1e-12
        #             veci(vcorrD,grid_v,iLS+1)[pII] -= (τ .* iMv * ∇ϕ_y)[pII]
        #         end
        #     end
        # end
    end

    # print("\n test u ", vecb_L(uD, grid_u))
    # print("\n test u ", vecb_R(uD, grid_u))
    # print("\n test u ", vecb_B(uD, grid_u))
    # print("\n test u ", vecb_T(uD, grid_u))

    #TODO Poiseuille
    printstyled(color=:yellow, @sprintf "\n before end pressure projection \n")

    test_Poiseuille(num,vD,grid_v)
    #TODO
    compute_grad_p!(num,grid, grid_u, grid_v, pD, opC_p, opC_u, opC_v)


    printstyled(color=:magenta, @sprintf "\n end pressure projection \n")

    return Lp, bc_Lp, bc_Lp_b, Lu, bc_Lu, bc_Lu_b, Lv, bc_Lv, bc_Lv_b, opC_p.M, opC_u.M, opC_v.M, Cui, Cvi
end


function test_laplacian_pressure(num,grid_v,ph, opC_p, Lv, bc_Lv, bc_Lv_b)

    @unpack rho1,rho2,visc_coeff = num
    @unpack vD = ph

    iRe = visc_coeff

    ni = grid_v.nx * grid_v.ny
    nb = 2 * grid_v.nx + 2 * grid_v.ny
    nt = (num.nLS + 1) * ni + nb

    Avtest = spzeros(nt, nt)
   
    # Implicit part of viscous term
    Avtest[1:ni,1:ni] = iRe .*Lv #pad_crank_nicolson(Lv, grid, τ)
    # Contribution to implicit part of viscous term from outer boundaries
    Avtest[1:ni,end-nb+1:end] = iRe .* bc_Lv_b

    # vecv = reshape(vec1(vD,grid_v),grid_v)

    iplot = 1
    jplot = 64
    II = CartesianIndex(jplot, iplot) #(id_y, id_x)
    pII = lexicographic(II, grid_v.ny)
    print("\n ")
    print("\n Laplacian coefficients ", II," ",Lv[pII,:]," ",bc_Lv_b[pII,:])

    # print("\n Laplacian coefficients ", II," ",Lv[pII,:])
    # print("\n testvisc ", II," ",bc_Lv_b[pII,:])
    # print("\n testvisc ", II," ",Avtest[pII,:])

    testAv2 = Avtest * vD .*rho1 
    
    printstyled(color=:green, @sprintf "\n exact %.10e 4/3exact %.10e\n" testAv2[pII]*opC_p.iMy.diag[pII] testAv2[pII]*opC_p.iMy.diag[pII]*4/3)

    return testAv2[pII]*opC_p.iMy.diag[pII]

end


function convert_interfacial_D_to_segments(num,gp,field,iLS)

    # x = zeros(0)
    # y = zeros(0)
    # f = zeros(0)
    x = Float64[]
    y = Float64[]
    f = Float64[]

    connectivities=Int64[]
    ij_index=-ones(Int64,(gp.ny, gp.nx))
    vtx_index = 0 #Integer(0) 

    for II in gp.LS[iLS].MIXED
        push!( x, gp.LS[1].mid_point[II].x * gp.dx[II] + gp.x[II] )
        push!( y, gp.LS[1].mid_point[II].y * gp.dy[II] + gp.y[II] )
        push!( f, field[II] )

        ij_index[II] = vtx_index

        # print("\n vtx ", vtx_index," ", II, " ",gp.x[II]-gp.dx[II]/2," ",gp.x[II]+gp.dx[II]/2," ",x," ",gp.y[II]-gp.dy[II]/2,
        # " ",gp.y[II]+gp.dy[II]/2,y," ",gp.x[II]," ",gp.y[II])


        #TODO store connvectivities or similar so that after: only need to filter and store non-null values ? 
        #no need to do the work n times for each scalar

        vtx_index +=1
    end

    num_seg = 0
    for j in 1:gp.ny-1
        for i in 1:gp.nx-1
            # II = CartesianIndex(i,j)
            if ij_index[j,i] >=0

                if ij_index[j+1,i]>=0  
                    push!(connectivities, ij_index[j,i])
                    push!(connectivities,ij_index[j+1,i])
                    num_seg+=1
                    # print("\n seg ", i," ",j," ",ij_index[j,i], " ", ij_index[j+1,i])
                end

                if ij_index[j,i+1]>=0  
                    push!(connectivities, ij_index[j,i])
                    push!(connectivities,ij_index[j,i+1])
                    num_seg+=1
                    # print("\n seg ", i," ",j," ",ij_index[j,i], " ", ij_index[j,i+1])
                end

            end

        end
    end

    # num_vtx = vtx_index already incremented vtx_index +=1
    
    return x,y,f,connectivities,vtx_index,num_seg

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
"""
center(r, θ)

Returns the 
"""
@inline center(r, θ) = r * cos(π - θ)


function contact_angle_advancing_receding(grid,grid_u, grid_v, iLS, II)

    # From code in run.jl

    # # TODO can be done witu u and v component ?
    # normalx = cos.(grid_u.LS[iLS].α)
    # normaly = sin.(grid_v.LS[iLS].α)

    # # grid_u.V .*= normalx
    # # grid_v.V .*= normaly

    # advancing_receding = normalx[II] * grid_u.V[II] + normaly[II] * grid_v.V[II] #TODO check dim , bounds


    #With interpolation

    normalx = cos.(grid.LS[iLS].α)
    normaly = sin.(grid.LS[iLS].α)

    cap1 = grid_u.LS[iLS].geoL.cap[II,5]
    cap3 = grid_u.LS[iLS].geoL.cap[δx⁺(II),5]
    tmpVx = (grid_u.V[II] * cap1 + grid_u.V[δx⁺(II)] * cap3) * inv_weight_eps(num,cap1 + cap3)
    
    cap2 = grid_v.LS[iLS].geoL.cap[II,5]
    cap4 = grid_v.LS[iLS].geoL.cap[δy⁺(II),5]
    tmpVy = (grid_v.V[II] * cap2 + grid_v.V[δy⁺(II)] * cap4) * inv_weight_eps(num,cap2 + cap4)

    advancing_receding = normalx[II] * tmpVx + normaly[II] * tmpVy

    return advancing_receding



    # V .= 0.0
    # cap1 = grid_u.LS[iLS].geoL.cap[II,5]
    # cap3 = grid_u.LS[iLS].geoL.cap[δx⁺(II),5]
    # tmpVx = (grid_u.V[II] * cap1 + grid_u.V[δx⁺(II)] * cap3) * inv_weight_eps(num,cap1 + cap3)
    
    # cap2 = grid_v.LS[iLS].geoL.cap[II,5]
    # cap4 = grid_v.LS[iLS].geoL.cap[δy⁺(II),5]
    # tmpVy = (grid_v.V[II] * cap2 + grid_v.V[δy⁺(II)] * cap4) * inv_weight_eps(num,cap2 + cap4)
    
    # tmpV = sqrt(tmpVx^2 + tmpVy^2)
    # β = atan(tmpVy, tmpVx)

    # printstyled(color=:magenta, @sprintf "\n beta β: %.2e : %.2e : %.2e" β grid.LS[iLS].α[II] )

    # if grid.LS[iLS].α[II] > 0.0 && β < 0.0
    #     β += 2π
    # end
    # if grid.LS[iLS].α[II] < 0.0 && β > 0.0
    #     β -= 2π
    # end

    # V[II] = tmpV * cos(β - grid.LS[iLS].α[II])

    
    
    
    # Project velocities to the normal and use advecion scheme for advection just
    # in the normal direction
    # tmpVx = zeros(grid)
    # tmpVy = zeros(grid)
    # V .= 0.0
    # @inbounds @threads for II in grid.LS[iLS].MIXED
    #     cap1 = grid_u.LS[iLS].geoL.cap[II,5]
    #     cap3 = grid_u.LS[iLS].geoL.cap[δx⁺(II),5]
    #     # tmpVx[II] = (grid_u.V[II] * cap1 + grid_u.V[δx⁺(II)] * cap3) / (cap1 + cap3 + eps(0.01))
    #     tmpVx[II] = (grid_u.V[II] * cap1 + grid_u.V[δx⁺(II)] * cap3) * inv_weight_eps(num,cap1 + cap3)
        

    #     cap2 = grid_v.LS[iLS].geoL.cap[II,5]
    #     cap4 = grid_v.LS[iLS].geoL.cap[δy⁺(II),5]
    #     # tmpVy[II] = (grid_v.V[II] * cap2 + grid_v.V[δy⁺(II)] * cap4) / (cap2 + cap4 + eps(0.01))
    #     tmpVy[II] = (grid_v.V[II] * cap2 + grid_v.V[δy⁺(II)] * cap4) * inv_weight_eps(num,cap2 + cap4)
        

    #     tmpV = sqrt(tmpVx[II]^2 + tmpVy[II]^2)
    #     β = atan(tmpVy[II], tmpVx[II])
    #     if grid.LS[iLS].α[II] > 0.0 && β < 0.0
    #         β += 2π
    #     end
    #     if grid.LS[iLS].α[II] < 0.0 && β > 0.0
    #         β -= 2π
    #     end

    #     V[II] = tmpV * cos(β - grid.LS[iLS].α[II])
    # end

    # i_ext, l_ext, b_ext, r_ext, t_ext = indices_extension(grid, grid.LS[iLS], grid.ind.inside, periodic_x, periodic_y)
    # field_extension!(grid, grid.LS[iLS].u, V, i_ext, l_ext, b_ext, r_ext, t_ext, num.NB, periodic_x, periodic_y)

    # rhs_LS .= 0.0
    # IIOE_normal!(grid, LS[iLS].A, LS[iLS].B, LS[iLS].u, V, CFL_sc, periodic_x, periodic_y)
    # BC_LS!(grid, LS[iLS].u, LS[iLS].A, LS[iLS].B, rhs_LS, BC_u)
    # BC_LS_interior!(num, grid, iLS, LS[iLS].A, LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
    # LS[iLS].u .= reshape(gmres(LS[iLS].A, LS[iLS].B * vec(LS[iLS].u) .+ rhs_LS), grid)

    # Impose contact angle if a wall is present
    # rhs_LS .= 0.0
    # LS[iLS].A.nzval .= 0.0
    # LS[iLS].B.nzval .= 0.0
    # for II in grid.ind.all_indices
    #     pII = lexicographic(II, grid.ny)
    #     LS[iLS].A[pII,pII] = 1.0
    #     LS[iLS].B[pII,pII] = 1.0
    # end
    # BC_LS_interior!(num, grid, iLS, LS[iLS].A, LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
    # LS[iLS].u .= reshape(gmres(LS[iLS].A, LS[iLS].B * vec(LS[iLS].u) .+ rhs_LS), grid)

end