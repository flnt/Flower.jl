
"""
tests Laplacian for Poiseuille flow
"""
function test_laplacian_pressure(num,grid_v,vD, opC_p, Lv, bc_Lv, bc_Lv_b)

    iRe = num.visc_coeff

    ni = grid_v.nx * grid_v.ny
    nb = 2 * grid_v.nx + 2 * grid_v.ny
    nt = (num.nLS + 1) * ni + nb

    Avtest = spzeros(nt, nt)
   
    # Implicit part of viscous term
    Avtest[1:ni,1:ni] = iRe .*Lv #pad_crank_nicolson(Lv, grid, timestep_n)
    # Contribution to implicit part of viscous term from outer boundaries
    Avtest[1:ni,end-nb+1:end] = iRe .* bc_Lv_b

    # vecv .= reshape(vec1(vD,grid_v),grid_v)

    iplot = 1
    jplot = 64
    II = CartesianIndex(jplot, iplot) #(id_y, id_x)
    pII = lexicographic(II, grid_v.ny)
    print("\n ")
    print("\n Laplacian coefficients ", II," ",Lv[pII,:]," ",bc_Lv_b[pII,:])

    # print("\n Laplacian coefficients ", II," ",Lv[pII,:])
    # print("\n testvisc ", II," ",bc_Lv_b[pII,:])
    # print("\n testvisc ", II," ",Avtest[pII,:])

    # testAv2 = Avtest * vD .* num.rho1 
    # testAv2 = Avtest * vD  
    testAv2 = Avtest * vD .* num.rho1 

    
    printstyled(color=:green, @sprintf "\n exact %.10e 4/3exact %.10e\n" testAv2[pII]*opC_p.iMy.diag[pII] testAv2[pII]*opC_p.iMy.diag[pII]*4/3)

    return testAv2[pII]*opC_p.iMy.diag[pII]

end


"""
test presence of interface
"""
function test_LS(grid)
    iLS = 1
    minLS = minimum(grid.LS[iLS].u)
    maxLS = maximum(grid.LS[iLS].u)

    printstyled(color=:red, @sprintf "\n levelset: min %.2e max %.2e\n" minLS maxLS)
    if maximum(grid.LS[iLS].u) < 0.0
        @error("LS initialization")
    end

    if minLS*maxLS > 0.0
        # @error("No interface")
        printstyled(color=:red, @sprintf "\n WARNING: no interface \n")

    end

end




function FE_set_momentum_debug(
    bc_type, num, grid, opC,
    A, B,
    L, bc_L, bc_L_b, Mm1, BC,
    ls_advection
    )
    @unpack timestep_n = num
    @unpack Bx, By, Hx, Hy, HxT, HyT, χ, M, iMx, iMy, Hx_b, Hy_b, HxT_b, HyT_b, iMx_b, iMy_b, iMx_bd, iMy_bd, χ_b = opC

    ni = grid.nx * grid.ny
    nb = 2 * grid.nx + 2 * grid.ny

    rhs = fnzeros(grid, num) #TODO remove alloc

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
        A[1:ni,1:ni] = pad_crank_nicolson(M .- timestep_n .* L, grid, timestep_n)
        # Contribution to implicit part of viscous term from outer boundaries
        A[1:ni,end-nb+1:end] = - timestep_n .* bc_L_b
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
            A[1:ni,sb] = - timestep_n .* bc_L[iLS]
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
    periodic_x, periodic_y, advection, ls_advection, current_iter, Ra, navier, pres_free_suface,jump_mass_transfer_rate,mass_transfer_rate
    )
    @unpack Re, timestep_n, σ, g, β, nLS, nNavier = num
    @unpack p, pD, ϕ, ϕD, u, v, u_predictionD, v_predictionD, uD, vD, u_prediction, v_prediction, uT = ph
    @unpack Cu, Cv, CUTCu, CUTCv = op_conv
    @unpack rho1,rho2,visc_coeff = num


    # iRe = 1.0 / Re
    iRe = visc_coeff
    itimestep_n = 1.0 / timestep_n
    irho1 = 1.0/rho1

    printstyled(color=:green, @sprintf "\n rho1 : %.2e irho1 : %.2e iRe : %.2e\n" rho1 irho1 iRe)

    ∇ϕ_x = opC_u.AxT * opC_u.Rx * vec1(pD,grid) .+ opC_u.Gx_b * vecb(pD,grid)
    ∇ϕ_y = opC_v.AyT * opC_v.Ry * vec1(pD,grid) .+ opC_v.Gy_b * vecb(pD,grid)
    for iLS in 1:nLS
        ∇ϕ_x .+= opC_u.Gx[iLS] * veci(pD,grid,iLS+1)
        ∇ϕ_y .+= opC_v.Gy[iLS] * veci(pD,grid,iLS+1)
    end

    grd_x .= reshape(veci(∇ϕ_x,grid_u,1), grid_u)
    grd_y .= reshape(veci(∇ϕ_y,grid_v,1), grid_v)

    printstyled(color=:cyan, @sprintf "\n B %.2e T %.2e L %.2e R %.2e\n" maximum(abs.(vecb_B(∇ϕ_x,grid_u))) maximum(abs.(vecb_T(∇ϕ_x,grid_u))) maximum(abs.(vecb_L(∇ϕ_y,grid_v))) maximum(abs.(vecb_R(∇ϕ_y,grid_v))))
    printstyled(color=:red, @sprintf "\n B %.2e T %.2e L %.2e R %.2e\n" maximum(abs.(grd_x[end,:])) maximum(abs.(grd_x[1,:])) maximum(abs.(grd_y[:,1])) maximum(abs.(grd_y[:,end])))

    printstyled(color=:red, @sprintf "\n p %.2e p %.2e p %.2e p %.2e \n" maximum(abs.(p[end,:])) maximum(abs.(p[1,:])) maximum(abs.(p[:,end])) maximum(abs.(p[:,1])))

    ptest .= reshape(veci(pD,grid,1), grid)
    printstyled(color=:magenta, @sprintf "\n ptest %.2e p %.2e p %.2e p %.2e \n" maximum(abs.(ptest[end,:])) maximum(abs.(ptest[1,:])) maximum(abs.(ptest[:,end])) maximum(abs.(ptest[:,1])))



    ∇ϕ_x .= 0.0
    ∇ϕ_y .= 0.0

    compute_grad_p_2!(num,grid, grid_u, grid_v, pD, opC_p, opC_u, opC_v)


    if num.prediction == 1 || num.prediction == 2
        printstyled(color=:red, @sprintf "\n pressure_in_prediction \n")

        #TODO
        compute_grad_p_2!(num,grid, grid_u, grid_v, pD, opC_p, opC_u, opC_v)


        ∇ϕ_x = opC_u.AxT * opC_u.Rx * vec1(pD,grid) .+ opC_u.Gx_b * vecb(pD,grid)
        ∇ϕ_y = opC_v.AyT * opC_v.Ry * vec1(pD,grid) .+ opC_v.Gy_b * vecb(pD,grid)
        for iLS in 1:nLS
            ∇ϕ_x .+= opC_u.Gx[iLS] * veci(pD,grid,iLS+1)
            ∇ϕ_y .+= opC_v.Gy[iLS] * veci(pD,grid,iLS+1)
        end

        grd_x .= reshape(veci(∇ϕ_x,grid_u,1), grid_u)
        grd_y .= reshape(veci(∇ϕ_y,grid_v,1), grid_v)

        printstyled(color=:cyan, @sprintf "\n B %.2e T %.2e L %.2e R %.2e\n" maximum(abs.(vecb_B(∇ϕ_x,grid_u))) maximum(abs.(vecb_T(∇ϕ_x,grid_u))) maximum(abs.(vecb_L(∇ϕ_y,grid_v))) maximum(abs.(vecb_R(∇ϕ_y,grid_v))))


        grd_xfull = opC_p.iMx * opC_p.Bx * vec1(pD,grid) .+ opC_p.iMx_b * opC_p.Hx_b * vecb(pD,grid)
        grd_yfull = opC_p.iMy * opC_p.By * vec1(pD,grid) .+ opC_p.iMy_b * opC_p.Hy_b * vecb(pD,grid)

        for iLS in 1:num.nLS
            grd_xfull .+= opC_p.iMx * opC_p.Hx[iLS] * veci(pD,grid,iLS+1)
            grd_yfull .+= opC_p.iMy * opC_p.Hy[iLS] * veci(pD,grid,iLS+1)
        end

        grd_x .= reshape(veci(grd_xfull,grid_u,1), grid_u)
        grd_y .= reshape(veci(grd_yfull,grid_v,1), grid_v)

        # printstyled(color=:red, @sprintf "\n grad min max x %.2e %.2e y %.2e %.2e\n" minimum(grd_x) maximum(grd_x) minimum(grd_y) maximum(grd_y))

        print("\n dt ", num.timestep_n)

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

        printstyled(color=:red, @sprintf "\n full grad min max x %.2e %.2e y %.2e %.2e\n" minimum(timestep_n.*irho1.*ph.Gxm1) maximum(timestep_n.*irho1.*ph.Gxm1) minimum(timestep_n.*irho1.*ph.Gym1) maximum(timestep_n.*irho1.*ph.Gym1))

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

    if is_Forward_Euler(time_scheme)
        rhs_u, rhs_v, rhs_ϕ, rhs_uv, Lp, bc_Lp, bc_Lp_b, Lu, bc_Lu, bc_Lu_b, Lv, bc_Lv, bc_Lv_b = set_Forward_Euler!(
            bc_int, num, grid, geo, grid_u, geo_u, grid_v, geo_v,
            opC_p, opC_u, opC_v, BC_p, BC_u, BC_v,
            Au, Bu, Av, Bv, Aϕ, Auv, Buv,
            Lpm1, bc_Lpm1, bc_Lpm1_b, Lum1, bc_Lum1, bc_Lum1_b, Lvm1, bc_Lvm1, bc_Lvm1_b,
            Mum1, Mvm1, iRe, op_conv, ph,
            periodic_x, periodic_y, advection, ls_advection, navier
        )
    elseif is_Crank_Nicolson(time_scheme)
        rhs_u, rhs_v, rhs_ϕ, Lp, bc_Lp, bc_Lp_b, Lu, bc_Lu, bc_Lu_b, Lv, bc_Lv, bc_Lv_b = set_Crank_Nicolson!(
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
        if current_iter == 1
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
        #     vec1(rhs_u,grid_u) .+= -timestep_n .* (opC_u.AxT * opC_u.Rx * vec1(pD,grid) .+ opC_u.Gx_b * vecb(pD,grid))
        #     for iLS in 1:nLS
        #         vec1(rhs_u,grid_u) .+= -timestep_n .* (opC_u.Gx[iLS] * veci(pD,grid,iLS+1))
        #     end
        # end
        mul!(rhs_u, Bu, uD, 1.0, 1.0)
        vec1(rhs_u,grid_u) .+= timestep_n .* grav_x
        vec1(rhs_u,grid_u) .-= timestep_n .* Convu
        vec1(rhs_u,grid_u) .+= timestep_n .* ra_x
        printstyled(color=:green, @sprintf "\n rhs u : %.2e uD %.2e Bu %.2e M %.2e \n" maximum(abs.(rhs_u)) maximum(abs.(uD)) maximum(abs.(Bu)) maximum(abs.(Mum1)))

        vec1(rhs_u,grid_u) .-= timestep_n .* irho1 .* ph.Gxm1 
        
        printstyled(color=:green, @sprintf "\n rhs u : %.2e \n" maximum(abs.(rhs_u)))

        kill_dead_cells!(vec1(rhs_u,grid_u), grid_u, geo_u[end])
        for iLS in 1:nLS
            kill_dead_cells!(veci(rhs_u,grid_u,iLS+1), grid_u, geo_u[end])
        end
        # @time bicgstabl!(u_predictionD, Au, rhs_u, log=true)
        try
            # @time bicgstabl!(u_predictionD, Au, rhs_u, Pl=Diagonal(Au), log=true)
            @time u_predictionD .= Au \ rhs_u
        catch e
            u_predictionD .= Inf
            println(e)
        end

        # printstyled(color=:green, @sprintf "\n max abs(u_predictionD) : %.2e uD: %.2e \n" maximum(abs.(u_predictionD)) maximum(abs.(uD)))

        kill_dead_cells!(vec1(u_predictionD,grid_u), grid_u, geo_u[end])
        for iLS in 1:nLS
            kill_dead_cells!(veci(u_predictionD,grid_u,iLS+1), grid_u, geo_u[end])
        end
        u_prediction .= reshape(vec1(u_predictionD,grid_u), grid_u)

        # if is_wall_no_slip(bc_int)
        #     vec1(vD,grid_v) .= vec(v)
        #     # update_dirichlet_field!(grid_v, vD, v, BC_v)
        #     vec1(rhs_v,grid_v) .+= -timestep_n .* (opC_v.AyT * opC_v.Ry * vec1(pD,grid) .+opC_v.Gy_b * vecb(pD,grid))
        #     for iLS in 1:nLS
        #         vec1(rhs_v,grid_v) .+= -timestep_n .* (opC_v.Gy[iLS] * veci(pD,grid,iLS+1))
        #     end
        # end
        mul!(rhs_v, Bv, vD, 1.0, 1.0)

        test1 = vec1(rhs_v,grid_v)[1,1]/Poiseuille_fmax(grid_v.x[1,1],num.v_inlet,num.L0)
        test2 = test1 / (grid_v.dx[1,1]^2/2)
        printstyled(color=:red, @sprintf "\n rhs_v vec1 %.10e /pois %.10e /pois %.10e\n" vec1(rhs_v,grid_v)[1,1] test1 test2)

        vec1(rhs_v,grid_v) .+= - timestep_n .* grav_y
        vec1(rhs_v,grid_v) .-= timestep_n .* Convv
        vec1(rhs_v,grid_v) .+= timestep_n .* ra_y

        test1 = vec1(rhs_v,grid_v)[1,1]/Poiseuille_fmax(grid_v.x[1,1],num.v_inlet,num.L0)
        test2 = test1 / (grid_v.dx[1,1]^2/2)
        test3 = vec1(rhs_v,grid_v)[1,1]-Poiseuille_fmax(grid_v.x[1,1],num.v_inlet,num.L0)*(grid_v.dx[1,1]^2/2)
        printstyled(color=:red, @sprintf "\n rhs_v vec1 %.10e /pois %.10e /pois %.10e diff %.10e\n" vec1(rhs_v,grid_v)[1,1] test1 test2 test3)

        # printstyled(color=:green, @sprintf "\n rhs: %.2e vD %.2e \n" maximum(abs.(rhs_v)) maximum(abs.(vD)))
        printstyled(color=:green, @sprintf "\n rhs v : %.2e vD %.2e Bv %.2e M %.2e \n" maximum(abs.(rhs_v)) maximum(abs.(vD)) maximum(abs.(Bv)) maximum(abs.(Mvm1)))


        vec1(rhs_v,grid_v) .-= timestep_n .* irho1 .* ph.Gym1
        printstyled(color=:green, @sprintf "\n rhs: %.2e \n" maximum(abs.(rhs_v)))


        test1 = vec1(rhs_v,grid_v)[1,1]/Poiseuille_fmax(grid_v.x[1,1],num.v_inlet,num.L0)
        test2 = test1 / (grid_v.dx[1,1]^2/2)
        test3 = vec1(rhs_v,grid_v)[1,1]-Poiseuille_fmax(grid_v.x[1,1],num.v_inlet,num.L0)*(grid_v.dx[1,1]^2/2)
        test4 = test3/(timestep_n .* irho1)/ (grid_v.dx[1,1]^2/2)
        printstyled(color=:red, @sprintf "\n rhs_v vec1 %.10e /pois %.10e /pois %.10e diff %.10e diff %.10e\n" vec1(rhs_v,grid_v)[1,1] test1 test2 test3 test4)



        kill_dead_cells!(vec1(rhs_v,grid_v), grid_v, geo_v[end])
        for iLS in 1:nLS
            kill_dead_cells!(veci(rhs_v,grid_v,iLS+1), grid_v, geo_v[end])
        end
        # bicgstabl!(v_predictionD, Av, rhs_v, log=true)
        
        
        iplot = 1
        jplot = 1
        II = CartesianIndex(jplot, iplot) #(id_y, id_x)
        # pII = lexicographic(II, grid.ny +1)

        print("\n after kill dead cells ", (grid_v.dx[1,1]^2/2)," full " ,(grid_v.dx[1,1]^2)," test ",geo_v[end].cap[II,5])
        test1 = vec1(rhs_v,grid_v)[1,1]/Poiseuille_fmax(grid_v.x[1,1],num.v_inlet,num.L0)
        test2 = test1 / (grid_v.dx[1,1]^2/2)
        test3 = vec1(rhs_v,grid_v)[1,1]-Poiseuille_fmax(grid_v.x[1,1],num.v_inlet,num.L0)*(grid_v.dx[1,1]^2/2)
        test4 = test3/(timestep_n .* irho1)/ (grid_v.dx[1,1]^2/2)
        printstyled(color=:red, @sprintf "\n rhs_v vec1 %.10e /pois %.10e /pois %.10e diff %.10e diff %.10e\n" vec1(rhs_v,grid_v)[1,1] test1 test2 test3 test4)


        try
            # @time bicgstabl!(v_predictionD, Av, rhs_v, Pl=Diagonal(Av), log=true)
            @time v_predictionD .= Av \ rhs_v
        catch e
            v_predictionD .= Inf
            println(e)
        end

        printstyled(color=:yellow, @sprintf "\n v_predictionD \n")

        # iplot = 64
        # jplot = 64
        # II = CartesianIndex(jplot, iplot) #(id_y, id_x)
        # pII = lexicographic(II, grid.ny +1)
        
        # # test = timestep_n .* iRe.*Lv *v_predictionD
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
        Avtest[1:ni,1:ni] = iRe .*Lv #pad_crank_nicolson(Lv, grid, timestep_n)
        # Contribution to implicit part of viscous term from outer boundaries
        Avtest[1:ni,end-nb+1:end] = iRe .* bc_Lv_b

        vecv .= reshape(vec1(vD,grid_v),grid_v)


        iplot = 1
        jplot = 64
        II = CartesianIndex(jplot, iplot) #(id_y, id_x)
        pII = lexicographic(II, grid.ny +1)
        print("\n ")
        print("\n testvisc ", II," ",Lv[pII,:])
        print("\n testvisc ", II," ",Avtest[pII,:])

        print("\n testvisc ", II," ",bc_Lv_b[pII,:])


        testAv = Avtest * v_predictionD .*rho1 
        testAv2 = Avtest * vD .*rho1 

        printstyled(color=:green, @sprintf "\n Avtest * v_predictionD/My : %.10e exact %.10e 4/3exact %.10e\n" testAv[pII]*opC_p.iMy.diag[pII] testAv2[pII]*opC_p.iMy.diag[pII] testAv2[pII]*opC_p.iMy.diag[pII]*4/3)

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

        printstyled(color=:green, @sprintf "\n Avtest * v_predictionD/My : %.10e exact %.10e\n" testAv[pII]*opC_p.iMy.diag[pII] testAv2[pII]*opC_p.iMy.diag[pII])
        ####################################################################################################

        ####################################################################################################        
        iplot = 1
        jplot = 1
        II = CartesianIndex(jplot, iplot) #(id_y, id_x)
        pII = lexicographic(II, grid.ny +1)
        print("\n ")
        print("\n testvisc ", II," ",Lv[pII,:])
        print("\n testvisc ", II," ",Avtest[pII,:])

        printstyled(color=:green, @sprintf "\n Avtest * v_predictionD/My : %.10e exact %.10e 4/3exact %.10e\n" testAv[pII]*opC_p.iMy.diag[pII] testAv2[pII]*opC_p.iMy.diag[pII] testAv2[pII]*opC_p.iMy.diag[pII]*4/3)
        ####################################################################################################

        #not 
        # printstyled(color=:green, @sprintf "\n Avtest * v_predictionD : %.10e Avtest * v_predictionD/M : %.10e Avtest * v_predictionD/My : %.10e\n" testAv[pII] testAv[pII]*opC_v.iMx_bd[pII,pII] testAv[pII]*opC_v.iMy[pII,pII])
        # ####################################################################################################        
        # iplot = 2
        # jplot = 64
        # II = CartesianIndex(jplot, iplot) #(id_y, id_x)
        # pII = lexicographic(II, grid.ny +1)
        # print("\n ")
        # print("\n testvisc ", II," ",Lv[pII,:])
        # printstyled(color=:green, @sprintf "\n Avtest * v_predictionD : %.10e Avtest * v_predictionD/M : %.10e Avtest * v_predictionD/My : %.10e \n" testAv[pII] testAv[pII]*opC_v.iMx[pII,pII] testAv[pII]*opC_v.iMy[pII,pII])
        # ####################################################################################################


        # 6.103515625000243e-13

        # testAv = Av * v_predictionD - 


        # testLv = fnzeros(grid, num)
        # testLv = fnzeros(grid, num)
        # mul!(testLv, Lv, v_predictionD, 1.0, 1.0)
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
        
        printstyled(color=:green, @sprintf "\n Avtest * v_predictionD/My : %.10e exact %.10e\n" testAv[pII]*opC_v.iMy[pII,pII] testAv2[pII]*opC_v.iMy[pII,pII])

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
        test_Poiseuille(num,v_predictionD,grid_v)

        printstyled(color=:red, @sprintf "\n v_predictionD %.2e %.2e\n" minimum(v_predictionD) maximum(v_predictionD))

        test_Poiseuille(num,vD,grid_v)


   


        printstyled(color=:red, @sprintf "\n vec1 1\n")
        print(vecv[1,:])

        printstyled(color=:red, @sprintf "\n vecb_B \n" )
        print(vecb_B(vD,grid_v))

        printstyled(color=:red, @sprintf "\n vecb_L vD\n")
        print(vecb_L(vD,grid_v))

        printstyled(color=:red, @sprintf "\n vecb_L v_predictionD\n" )
        print(vecb_L(v_predictionD,grid_v))


        printstyled(color=:red, @sprintf "\n rhs_v vecb_L \n" )
        print(vecb_L(rhs_v,grid_v))

        kill_dead_cells!(vec1(v_predictionD,grid_v), grid_v, geo_v[end])
        for iLS in 1:nLS
            kill_dead_cells!(veci(v_predictionD,grid_v,iLS+1), grid_v, geo_v[end])
        end
        v_prediction .= reshape(vec1(v_predictionD,grid_v), grid_v)
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

        rhs_uv[1:niu] .+= timestep_n .* grav_x
        rhs_uv[1:niu] .-= timestep_n .* Convu
        rhs_uv[1:niu] .+= timestep_n .* ra_x
        rhs_uv[1:niu] .-= timestep_n .* irho1 .* ph.Gxm1 

        rhs_uv[ntu+1:ntu+niv] .+= timestep_n .* grav_y
        rhs_uv[ntu+1:ntu+niv] .-= timestep_n .* Convv
        rhs_uv[ntu+1:ntu+niv] .+= timestep_n .* ra_y
        rhs_uv[ntu+1:ntu+niv] .-= timestep_n .* irho1 .* ph.Gym1 

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

        vec1(u_predictionD, grid_u) .= uvD[1:niu]
        vecb(u_predictionD, grid_u) .= uvD[ntu-nbu+1:ntu]
        kill_dead_cells!(vec1(u_predictionD,grid_u), grid_u, geo_u[end])
        u_prediction .= reshape(vec1(u_predictionD,grid_u), grid_u)

        vec1(v_predictionD, grid_v) .= uvD[ntu+1:ntu+niv]
        vecb(v_predictionD, grid_v) .= uvD[ntu+ntv-nbv+1:ntu+ntv]
        kill_dead_cells!(vec1(v_predictionD,grid_v), grid_v, geo_v[end])
        v_prediction .= reshape(vec1(v_predictionD,grid_v), grid_v)

        nNav = 0
        _iLS = 1
        for iLS in 1:nLS
            if !is_navier(bc_int[iLS]) && !is_navier_cl(bc_int[iLS])
                veci(u_predictionD,grid_u,iLS+1) .= uvD[_iLS*niu+1:(_iLS+1)*niu]
                kill_dead_cells!(veci(u_predictionD,grid_u,iLS+1), grid_u, geo_u[end])

                veci(v_predictionD,grid_v,iLS+1) .= uvD[ntu+_iLS*niv+1:ntu+(_iLS+1)*niv]
                kill_dead_cells!(veci(v_predictionD,grid_v,iLS+1), grid_v, geo_v[end])
                _iLS += 1
            else
                @inbounds uT[nNav+1,:] .= vec(uvD[ntu+ntv+1+nNav*nip:ntu+ntv+(nNav+1)*nip])
                nNav += 1
            end
        end
    end

    # printstyled(color=:green, @sprintf "\n max abs(u_predictionD) : %.2e v_predictionD %.2e \n" maximum(abs.(u_predictionD)) maximum(abs.(v_predictionD)))


    Duv = opC_p.AxT * vec1(u_predictionD,grid_u) .+ opC_p.Gx_b * vecb(u_predictionD,grid_u) .+
          opC_p.AyT * vec1(v_predictionD,grid_v) .+ opC_p.Gy_b * vecb(v_predictionD,grid_v)
    for iLS in 1:nLS
        if !is_navier(bc_int[iLS]) && !is_navier_cl(bc_int[iLS])
            Duv .+= opC_p.Gx[iLS] * veci(u_predictionD,grid_u,iLS+1) .+ 
                    opC_p.Gy[iLS] * veci(v_predictionD,grid_v,iLS+1)
        end
    end

    #Poisson equation
    # vec1(rhs_ϕ,grid) .= itimestep_n .* Duv
    vec1(rhs_ϕ,grid) .= rho1 .* itimestep_n .* Duv #TODO
    # veci(rhs_ϕ,grid) .*= rho1 #TODO

    # pres_free_suface = 0.0
    #TODO Marangoni
    #TODO phase change
    diff_inv_rho = 1.0/rho1 - 1.0/rho2
    # jump_mass_transfer_rate = 0.0 #TODO

    if jump_mass_transfer_rate
        for iLS in 1:nLS
            if is_fs(bc_int[iLS])
                Smat = strain_rate(iLS, opC_u, opC_v, opC_p)
                S = Smat[1,1] * vec1(u_predictionD,grid_u) .+ Smat[1,2] * veci(u_predictionD,grid_u,iLS+1) .+
                    Smat[2,1] * vec1(v_predictionD,grid_v) .+ Smat[2,2] * veci(v_predictionD,grid_v,iLS+1)
    
                fs_mat = opC_p.HxT[iLS] * opC_p.Hx[iLS] .+ opC_p.HyT[iLS] * opC_p.Hy[iLS]
                veci(rhs_ϕ,grid,iLS+1) .= -2.0 .* iRe .* S .+ Diagonal(diag(fs_mat)) * ( σ .* vec(grid.LS[iLS].κ) .- pres_free_suface .- diff_inv_rho * mass_transfer_rate ^ 2)
            end
        end
    else
        for iLS in 1:nLS
            if is_fs(bc_int[iLS])
                Smat = strain_rate(iLS, opC_u, opC_v, opC_p)
                S = Smat[1,1] * vec1(u_predictionD,grid_u) .+ Smat[1,2] * veci(u_predictionD,grid_u,iLS+1) .+
                    Smat[2,1] * vec1(v_predictionD,grid_v) .+ Smat[2,2] * veci(v_predictionD,grid_v,iLS+1)

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
        vec1(pD,grid) .= vec(ϕ .- iRe .* rho1 .* reshape(iM * Duv,grid)) #no timestep_n  since div u not rho1
    elseif num.prediction == 2
        vec1(pD,grid) .+= vec(ϕ .- iRe./2 .* rho1 .* reshape(iM * Duv,grid)) #no timestep_n  since div u not rho1
    else
        vec1(pD,grid) .= vec(ϕ) #.- iRe .* reshape(iM * Duv, grid))
    end
    for iLS in 1:nLS
        veci(pD,grid,iLS+1) .= veci(ϕD,grid,iLS+1)
    end
    vecb(pD,grid) .= vecb(ϕD,grid)
    p .= reshape(vec1(pD,grid), grid)

    #TODO
    compute_grad_p_2!(num,grid, grid_u, grid_v, pD, opC_p, opC_u, opC_v)


    # else
    #     vec1(pD,grid) .= vec(p) .+ vec(ϕ) #.- iRe .* iM * Duv
    #     vec2(pD,grid) .+= vec2(ϕD,grid)
    #     vecb(pD,grid) .+= vecb(ϕD,grid)
    #     p .= reshape(vec1(pD,grid), grid)
    # end

    # vec1(∇ϕ_x,grid) .*= irho1 
    # vec1(∇ϕ_y,grid) .*= irho1

    
    # u .= u_prediction .- timestep_n .* reshape(iMu * ∇ϕ_x, grid_u)
    # v .= v_prediction .- timestep_n .* reshape(iMv * ∇ϕ_y, grid_v)

    u .= u_prediction .- timestep_n .* irho1 .* reshape(iMu * ∇ϕ_x, grid_u)
    v .= v_prediction .- timestep_n .* irho1 .* reshape(iMv * ∇ϕ_y, grid_v)

    kill_dead_cells!(u, grid_u, geo_u[end])
    kill_dead_cells!(v, grid_v, geo_v[end])

    vec1(uD,grid_u) .= vec(u)
    vecb(uD,grid_u) .= vecb(u_predictionD,grid_u)
    vec1(vD,grid_v) .= vec(v)
    vecb(vD,grid_v) .= vecb(v_predictionD,grid_v)
    for iLS in 1:nLS
        if !is_navier(bc_int[iLS]) && !is_navier_cl(bc_int[iLS])
            veci(uD,grid_u,iLS+1) .= veci(u_predictionD,grid_u,iLS+1)
            veci(vD,grid_v,iLS+1) .= veci(v_predictionD,grid_v,iLS+1)
        end
        # if is_fs(bc_int[iLS])
        #     @inbounds for II in grid_u.ind.all_indices
        #         pII = lexicographic(II, grid_u.ny)
        #         if abs(veci(u_predictionD,grid_u,iLS+1)[pII]) > 1e-12
        #             veci(u_predictionD,grid_u,iLS+1)[pII] -= (timestep_n .* iMu * ∇ϕ_x)[pII]
        #         end
        #     end
        #     @inbounds for II in grid_v.ind.all_indices
        #         pII = lexicographic(II, grid_v.ny)
        #         if abs(veci(v_predictionD,grid_v,iLS+1)[pII]) > 1e-12
        #             veci(v_predictionD,grid_v,iLS+1)[pII] -= (timestep_n .* iMv * ∇ϕ_y)[pII]
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
    compute_grad_p_2!(num,grid, grid_u, grid_v, pD, opC_p, opC_u, opC_v)


    printstyled(color=:magenta, @sprintf "\n end pressure projection \n")

    return Lp, bc_Lp, bc_Lp_b, Lu, bc_Lu, bc_Lu_b, Lv, bc_Lv, bc_Lv_b, opC_p.M, opC_u.M, opC_v.M, Cui, Cvi
end

