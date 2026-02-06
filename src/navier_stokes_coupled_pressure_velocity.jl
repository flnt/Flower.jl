"""
Set the system matrix for Forward-Euler scheme

* bc_interface: interface boundary condition
* B contains Mum1 and Mvm1 (cell volumes)

## Modifications to `rhs`

1. **For each interface `iLS`:**
   - If the boundary condition is not Navier and not Navier-CL:
     ```julia
     @inbounds rhs[sbu] .= opu.χ[iLS] * vec(a0u)
     @inbounds rhs[sbv] .= opv.χ[iLS] * vec(a0v)
     ```
   - Otherwise (if the boundary condition is Navier or Navier-CL):
     ```julia
     @inbounds rhs[ntu+ntv+1+nNav1*nip:ntu+ntv+(nNav1+1)*nip] .= opp.χ[iLS] * vec(a0p)
     nNav1 += 1
     ```

2. **For the outer boundaries/borders:**
    ```julia
    @inbounds rhs[ntu-nbu+1:ntu] .= opu.χ_b * vec(a0_bu) #BC for u component on borders
    @inbounds rhs[ntu+ntv-nbv+1:ntu+ntv] .= opv.χ_b * vec(a0_bv) #BC for v component on borders
    ```
robin BC : source term a0
                
At the moment, the Levelset is not computed at borders/interfaces ? Only bulk
"""
function FE_set_momentum_coupled2(
    bc_interface, num, gp, gu, gv,
    opp, opu, opv,
    A, B,
    rhs,
    Lu, bc_Lu, bc_Lu_b, Mum1, BCu,
    Lv, bc_Lv, bc_Lv_b, Mvm1, BCv,
    ls_advection,BCp,ph=nothing
    )
    @unpack timestep_n, Re, nLS, nNavier = num


    printstyled(color=:red, @sprintf "\n coupled pressure-velocity FE_set_momentum_coupled2\n")


    #region init
    iRe = num.visc_coeff

    nip = gp.nx * gp.ny
    nbp = 2 * gp.nx + 2 * gp.ny

    niu = gu.nx * gu.ny
    nbu = 2 * gu.nx + 2 * gu.ny
    ntu = (nLS - nNavier + 1) * niu + nbu

    niv = gv.nx * gv.ny
    nbv = 2 * gv.nx + 2 * gv.ny
    ntv = (nLS - nNavier + 1) * niv + nbv

    #Reset to zero
    rhs .= 0.0 

    #region BC borders u
    a0_bu = zeros(nbu)
    _a1_bu = zeros(nbu)
    _b_bu = zeros(nbu)
    for iLS in 1:num.nLS
        set_borders!(gu, gu.LS[iLS].cl, gu.LS[iLS].u, a0_bu, _a1_bu, _b_bu, BCu, num.n_ext_cl)
    end
    a1_bu = Diagonal(vec(_a1_bu))
    b_bu = Diagonal(vec(_b_bu))
    #endregion BC borders u

    #region BC borders v
    a0_bv = zeros(nbv)
    _a1_bv = zeros(nbv)
    _b_bv = zeros(nbv)
    for iLS in 1:num.nLS
        set_borders!(gv, gv.LS[iLS].cl, gv.LS[iLS].u, a0_bv, _a1_bv, _b_bv, BCv, num.n_ext_cl)
    end
    a1_bv = Diagonal(vec(_a1_bv))
    b_bv = Diagonal(vec(_b_bv))
    #endregion BC borders v


    #region BC borders p
    if num.pressure_velocity_coupling == 2
        # rhs = fnzeros(gp, num)

        a0_bp = zeros(nbp)
        _a1_bp = zeros(nbp)
        _b_bp = zeros(nbp)
        for iLS in 1:num.nLS
            set_borders_poisson!(gp, gp.LS[iLS].cl, gp.LS[iLS].u, a0_bp, _a1_bp, _b_bp, BCp, num.n_ext_cl)
        end
        a1_bp = Diagonal(vec(_a1_bp))
        b_bp = Diagonal(vec(_b_bp))
    end

    # if ls_advection
    #     # Poisson equation
    #     # A[1:ni,1:ni] = pad(L, -4.0)
    #     # A[1:ni,end-nb+1:end] = bc_L_b

    #     # Boundary conditions for outer boundaries
    #     A[end-nb+1:end,1:ni] = b_bp * (HxT_b * iMx_b' * Bx .+ HyT_b * iMy_b' * By) 
    #     A[end-nb+1:end,end-nb+1:end] = pad(b_bp * (HxT_b * iMx_bd * Hx_b .+ HyT_b * iMy_bd * Hy_b) .+ χ_b * a1_bp,- 4.0)
    # end

    # for iLS in 1:num.nLS
    #     if ls_advection
    #         if is_dirichlet(bc_type[iLS])
    #             __a1 = 1.0
    #             __a2 = 0.0
    #             __b = 0.0
    #         elseif is_neumann(bc_type[iLS])
    #             __a1 = 0.0
    #             __a2 = 0.0
    #             __b = 1.0
    #         elseif is_robin(bc_type[iLS])
    #             __a1 = 1.0
    #             __a2 = 0.0
    #             __b = 1.0
    #         elseif is_fs(bc_type[iLS])
    #             __a1 = 0.0
    #             __a2 = 1.0
    #             __b = 0.0
    #         elseif is_wall_no_slip(bc_type[iLS])
    #             __a1 = 0.0
    #             __a2 = 0.0
    #             __b = 1.0
    #         elseif is_navier(bc_type[iLS])
    #             __a1 = 0.0
    #             __a2 = 0.0
    #             __b = 1.0
    #         elseif is_navier_cl(bc_type[iLS])
    #             __a1 = 0.0
    #             __a2 = 0.0
    #             __b = 1.0
    #         else
    #             __a1 = 0.0
    #             __a2 = 0.0
    #             __b = 1.0
    #         end
    
    #         _a1 = ones(gp) .* __a1
    #         a1 = Diagonal(vec(_a1))
    #         _a2 = ones(gp) .* __a2
    #         a2 = Diagonal(vec(_a2))
    #         _b = ones(gp) .* __b
    #         b = Diagonal(vec(_b))

    #         fs_mat = HxT[iLS] * Hx[iLS] .+ HyT[iLS] * Hy[iLS]

    #         sb = iLS*ni+1:(iLS+1)*ni
            
    #         # Poisson equation
    #         A[1:ni,sb] = bc_L[iLS]
    #         # Boundary conditions for inner boundaries
    #         A[sb,1:ni] = b * (HxT[iLS] * iMx * Bx .+ HyT[iLS] * iMy * By)
    #         # Contribution to Neumann BC from other boundaries
    #         for i in 1:num.nLS
    #             if i != iLS
    #                 A[sb,i*ni+1:(i+1)*ni] = b * (HxT[iLS] * iMx * Hx[i] .+ HyT[iLS] * iMy * Hy[i])
    #             end
    #         end
    #         A[sb,sb] = pad(
    #             b *(HxT[iLS] * iMx * Hx[iLS] .+ HyT[iLS] * iMy * Hy[iLS]) .+ χ[iLS] * a1 .-
    #             a2 * Diagonal(diag(fs_mat)), -4.0
    #         )
    #         #TODO was +b... here in set_poisson 
    #         A[sb,end-nb+1:end] = b * (HxT[iLS] * iMx_b * Hx_b .+ HyT[iLS] * iMy_b * Hy_b) #why was it different in set_poisson ?
    #         # Boundary conditions for outer boundaries
    #         A[end-nb+1:end,sb] = b_b * (HxT_b * iMx_b' * Hx[iLS] .+ HyT_b * iMy_b' * Hy[iLS])
    #     end

    #     veci(rhs,gp,iLS+1) .= +χ[iLS] * vec(a0[iLS]) #was - in set_poisson
    # end

    # vecb(rhs,gp) .= +χ_b * vec(a0_b) #was - in set_poisson
    #endregion BC borders p


    #endregion init

    # Array indices 
    bulk_u_velocity = 1:niu
    bulk_v_velocity = ntu+1:ntu+niv
    # bulk_tangential_velocity = 

    border_u_velocity = ntu-nbu+1:ntu
    border_v_velocity = ntu+ntv-nbv+1:ntu+ntv

    # nt = (num.nLS - num.nNavier + 1) * ni_uv + num.nNavier * nip + nb_uv + (num.nLS + 1) * nip + nbp
    ntNavier = num.nNavier * nip
    bulk_pressure = ntu+ntv+ntNavier+1:ntu+ntv+ntNavier+nip
    #ntu+ntv+1:ntu+ntv+nip
    border_pressure = ntu+ntv+ntNavier+(num.nLS + 1)*nip+1:ntu+ntv+ntNavier+(num.nLS + 1)*nip+nbp

    print("\n len rhs_uv ",size(rhs))

    print("bulk_pressure ",bulk_pressure)
    print("border_pressure ",border_pressure)


    if ls_advection
        A.nzval .= 0.0
        # Implicit part of viscous term
        A[bulk_u_velocity,bulk_u_velocity] = pad_crank_nicolson(opu.M .- timestep_n .* Lu, gu, timestep_n)
        # Contribution to implicit part of viscous term from outer boundaries
        A[bulk_u_velocity,border_u_velocity] = - timestep_n .* bc_Lu_b

        # Boundary conditions for outer boundaries
        A[border_u_velocity,bulk_u_velocity] = b_bu * (opu.HxT_b * opu.iMx_b' * opu.Bx .+ opu.HyT_b * opu.iMy_b' * opu.By)
        A[border_u_velocity,border_u_velocity] = pad(b_bu * (
            opu.HxT_b * opu.iMx_bd * opu.Hx_b .+ 
            opu.HyT_b * opu.iMy_bd * opu.Hy_b
        ) .- opu.χ_b * a1_bu)

        # Implicit part of viscous term
        A[bulk_v_velocity,bulk_v_velocity] = pad_crank_nicolson(opv.M .- timestep_n .* Lv, gv, timestep_n)
        # Contribution to implicit part of viscous term from outer boundaries
        A[bulk_v_velocity,border_v_velocity] = - timestep_n .* bc_Lv_b

        
        # Boundary conditions for outer boundaries
        A[border_v_velocity,bulk_v_velocity] = b_bv * (opv.HxT_b * opv.iMx_b' * opv.Bx .+ opv.HyT_b * opv.iMy_b' * opv.By)
        A[border_v_velocity,border_v_velocity] = pad(b_bv * (
            opv.HxT_b * opv.iMx_bd * opv.Hx_b .+ 
            opv.HyT_b * opv.iMy_bd * opv.Hy_b
        ) .- opv.χ_b * a1_bv)

        # TODO pad 1 or -4
        #TODO sign divergence not same u v and p

        #region coupled pression
        
        if num.pressure_velocity_coupling == 2

            # Poisson equation
            # A[1:ni,1:ni] = pad(L, -4.0)
            # A[1:ni,end-nb+1:end] = bc_L_b
    
            # Boundary conditions for outer boundaries
            A[border_pressure,bulk_pressure] = b_bp * (opp.HxT_b * opp.iMx_b' * opp.Bx .+ opp.HyT_b * opp.iMy_b' * opp.By) 
            A[border_pressure,border_pressure] = pad(b_bp * (
                opp.HxT_b * opp.iMx_bd * opp.Hx_b .+ 
                opp.HyT_b * opp.iMy_bd * opp.Hy_b) .+ opp.χ_b * a1_bp) #-4.0
            
        end
        
        # Implicit gradient of pressure (volume integrated)

        # timestep_n .* irho1 .*
        irho1 = 1.0/num.rho1
        factor = num.timestep_n * irho1
        A[bulk_u_velocity,bulk_pressure] = factor * opu.AxT * opu.Rx #TODO + or multiply by cell volume ? + not required, only component of matrix
        A[bulk_v_velocity,bulk_pressure] = factor * opv.AyT * opv.Ry 
       
        if num.pressure_velocity_coupling !=3
            #Outer boundaries
            A[bulk_u_velocity,border_pressure] = factor * opu.Gx_b
            A[bulk_v_velocity,border_pressure] = factor * opv.Gy_b
        

            for iLS in 1:nLS
                interfacial_nb_iLS_pressure = ntu+ntv+nNavier*nip+nip*iLS+1:ntu+ntv+nNavier*nip+(iLS+1)*nip
                A[bulk_u_velocity,interfacial_nb_iLS_pressure] .+= factor * opu.Gx[iLS] 
                A[bulk_v_velocity,interfacial_nb_iLS_pressure] .+= factor * opv.Gy[iLS] 
            end
        
        end #num.pressure_velocity_coupling !=3

        #region check coupled matrix
        #Test pressure part
        # Adummy = copy(A)
        # Adummy.nzval .= 0.0

        # nt = ntu+ntv + nNavier * nip + (num.nLS + 1) * nip + nbp
        # ncol_A = ntu+ntv + nNavier * nip + nip

        ni_p = gp.nx * gp.ny
        nb_p = 2 * gp.nx + 2 * gp.ny
    
        ni_u = gu.nx * gu.ny
        nb_u = 2 * gu.nx + 2 * gu.ny
    
        ni_v = gv.nx * gv.ny
        nb_v = 2 * gv.nx + 2 * gv.ny

        ni_uv = ni_u + ni_v
        nb_uv = nb_u + nb_v

        nt = (num.nLS - num.nNavier + 1) * ni_uv + num.nNavier * ni_p + nb_uv + (num.nLS + 1) * ni_p + nb_p
        ncol_A = (num.nLS - num.nNavier + 1) * ni_uv + num.nNavier * ni_p + nb_uv + ni_p

        Adummy = spzeros(ncol_A, nt)


        # length_colptr = length(Adummy.colptr)
        # length_rowval = length(Adummy.rowval)
        # length_nzval = length(Adummy.nzval)

        # PDI_status = @ccall "libpdi".PDI_multi_expose("print_matrix"::Cstring,
        # "Auv_n"::Cstring, Adummy.n::Ref{Clonglong}, PDI_OUT::Cint,
        # "Auv_m"::Cstring, Adummy.m::Ref{Clonglong}, PDI_OUT::Cint,
        # "Auv_colptr_len"::Cstring, length_colptr::Ref{Clonglong}, PDI_OUT::Cint,
        # "Auv_rowval_len"::Cstring, length_rowval::Ref{Clonglong}, PDI_OUT::Cint,
        # "Auv_nzval_len"::Cstring, length_nzval::Ref{Clonglong}, PDI_OUT::Cint,
        # # "Auv_colptr_len"::Cstring, length(Adummy.colptr)::Ref{Clonglong}, PDI_OUT::Cint,
        # # "Auv_rowval_len"::Cstring, length(Adummy.rowval)::Ref{Clonglong}, PDI_OUT::Cint,
        # # "Auv_nzval_len"::Cstring, length(Adummy.nzval)::Ref{Clonglong}, PDI_OUT::Cint,
        # "Auv_colptr_1D"::Cstring, Adummy.colptr::Ptr{Clonglong}, PDI_OUT::Cint,
        # "Auv_rowval_1D"::Cstring, Adummy.rowval::Ptr{Clonglong}, PDI_OUT::Cint,
        # "Auv_nzval_1D"::Cstring, Adummy.nzval::Ptr{Cdouble}, PDI_OUT::Cint,
        # C_NULL::Ptr{Cvoid})::Cint

        irho1 = 1.0/num.rho1
        factor = num.timestep_n * irho1
        Adummy[bulk_u_velocity,bulk_pressure] = factor * opu.AxT * opu.Rx #TODO + or multiply by cell volume ? + not required, only component of matrix
        Adummy[bulk_v_velocity,bulk_pressure] = factor * opv.AyT * opv.Ry 
        #Outer boundaries
        Adummy[bulk_u_velocity,border_pressure] = factor * opu.Gx_b
        Adummy[bulk_v_velocity,border_pressure] = factor * opv.Gy_b


        for iLS in 1:nLS
            interfacial_nb_iLS_pressure = ntu+ntv+nNavier*nip+nip*iLS+1:ntu+ntv+nNavier*nip+(iLS+1)*nip
            Adummy[bulk_u_velocity,interfacial_nb_iLS_pressure] .+= factor * opu.Gx[iLS] 
            Adummy[bulk_v_velocity,interfacial_nb_iLS_pressure] .+= factor * opv.Gy[iLS] 
        end

        # print("\n A before ")

        # # print("\n  Auv.n ", A.n)
        # # print("\n  Auv.m ", A.m)

        # # print("\n  length(A.colptr) ", length(A.colptr))
        # # print("\n  length(A.rowval) ", length(A.rowval))
        # # print("\n  length(A.nzval) ", length(A.nzval))

        # # print("\n  Auv.colptr ", A.colptr)
        # # print("\n  Auv.rowval ", A.rowval)
        # # print("\n  Auv.nzval ", A.nzval)

        # length_colptr = length(A.colptr)
        # length_rowval = length(A.rowval)
        # length_nzval = length(A.nzval)


        # # nx: $nx
        # # ny: $ny
        # # nNavier: $nb_Navier_slip_BC
        # # nb_levelsets: $nb_levelsets

        # PDI_status = @ccall "libpdi".PDI_multi_expose("print_matrix_test"::Cstring,
        # "Auv_n"::Cstring, A.n::Ref{Clonglong}, PDI_OUT::Cint,
        # "Auv_m"::Cstring, A.m::Ref{Clonglong}, PDI_OUT::Cint,
        # "nx"::Cstring, gp.nx::Ref{Clonglong}, PDI_OUT::Cint,
        # "ny"::Cstring, gp.ny::Ref{Clonglong}, PDI_OUT::Cint,
        # "nb_Navier_slip_BC"::Cstring, num.nNavier::Ref{Clonglong}, PDI_OUT::Cint,
        # # "Auv_colptr_len"::Cstring, length_colptr::Ref{Clonglong}, PDI_OUT::Cint,
        # # "Auv_rowval_len"::Cstring, length_rowval::Ref{Clonglong}, PDI_OUT::Cint,
        # # "Auv_nzval_len"::Cstring, length_nzval::Ref{Clonglong}, PDI_OUT::Cint,
        # # # "Auv_colptr_len"::Cstring, length(A.colptr)::Ref{Clonglong}, PDI_OUT::Cint,
        # # # "Auv_rowval_len"::Cstring, length(A.rowval)::Ref{Clonglong}, PDI_OUT::Cint,
        # # # "Auv_nzval_len"::Cstring, length(A.nzval)::Ref{Clonglong}, PDI_OUT::Cint,
        # # "Auv_colptr_1D"::Cstring, A.colptr::Ptr{Clonglong}, PDI_OUT::Cint,
        # # "Auv_rowval_1D"::Cstring, A.rowval::Ptr{Clonglong}, PDI_OUT::Cint,
        # # "Auv_nzval_1D"::Cstring, A.nzval::Ptr{Cdouble}, PDI_OUT::Cint,
        # C_NULL::Ptr{Cvoid})::Cint

        # PDI_status = @ccall "libpdi".PDI_multi_expose("print_matrix"::Cstring,
        # "Auv_n"::Cstring, A.n::Ref{Clonglong}, PDI_OUT::Cint,
        # "Auv_m"::Cstring, A.m::Ref{Clonglong}, PDI_OUT::Cint,
        # "Auv_colptr_len"::Cstring, length_colptr::Ref{Clonglong}, PDI_OUT::Cint,
        # "Auv_rowval_len"::Cstring, length_rowval::Ref{Clonglong}, PDI_OUT::Cint,
        # "Auv_nzval_len"::Cstring, length_nzval::Ref{Clonglong}, PDI_OUT::Cint,
        # # "Auv_colptr_len"::Cstring, length(A.colptr)::Ref{Clonglong}, PDI_OUT::Cint,
        # # "Auv_rowval_len"::Cstring, length(A.rowval)::Ref{Clonglong}, PDI_OUT::Cint,
        # # "Auv_nzval_len"::Cstring, length(A.nzval)::Ref{Clonglong}, PDI_OUT::Cint,
        # "Auv_colptr_1D"::Cstring, A.colptr::Ptr{Clonglong}, PDI_OUT::Cint,
        # "Auv_rowval_1D"::Cstring, A.rowval::Ptr{Clonglong}, PDI_OUT::Cint,
        # "Auv_nzval_1D"::Cstring, A.nzval::Ptr{Cdouble}, PDI_OUT::Cint,
        # C_NULL::Ptr{Cvoid})::Cint

        # print("\n Adummy ",Adummy)
        print("\n Adummy ")

        # print("\n  Auv.n ", Adummy.n)
        # print("\n  Auv.m ", Adummy.m)

        # print("\n  length(Adummy.colptr) ", length(Adummy.colptr))
        # print("\n  length(Adummy.rowval) ", length(Adummy.rowval))
        # print("\n  length(Adummy.nzval) ", length(Adummy.nzval))

        # print("\n  Auv.colptr ", Adummy.colptr)
        # print("\n  Auv.rowval ", Adummy.rowval)
        # print("\n  Auv.nzval ", Adummy.nzval)

      
        #region print grad p part of matrix
        PDI_status = @ccall "libpdi".PDI_multi_expose("print_matrix"::Cstring,
        "Auv_n"::Cstring, Adummy.n::Ref{Clonglong}, PDI_OUT::Cint,
        "Auv_m"::Cstring, Adummy.m::Ref{Clonglong}, PDI_OUT::Cint,
        # "Auv_colptr_len"::Cstring, length_colptr::Ref{Clonglong}, PDI_OUT::Cint,
        # "Auv_rowval_len"::Cstring, length_rowval::Ref{Clonglong}, PDI_OUT::Cint,
        # "Auv_nzval_len"::Cstring, length_nzval::Ref{Clonglong}, PDI_OUT::Cint,
        "Auv_colptr_len"::Cstring, length(Adummy.colptr)::Ref{Clonglong}, PDI_OUT::Cint,
        "Auv_rowval_len"::Cstring, length(Adummy.rowval)::Ref{Clonglong}, PDI_OUT::Cint,
        "Auv_nzval_len"::Cstring, length(Adummy.nzval)::Ref{Clonglong}, PDI_OUT::Cint,
        "Auv_colptr_1D"::Cstring, Adummy.colptr::Ptr{Clonglong}, PDI_OUT::Cint,
        "Auv_rowval_1D"::Cstring, Adummy.rowval::Ptr{Clonglong}, PDI_OUT::Cint,
        "Auv_nzval_1D"::Cstring, Adummy.nzval::Ptr{Cdouble}, PDI_OUT::Cint,
        C_NULL::Ptr{Cvoid})::Cint
        #endregion print grad p part of matrix
   


        uvD = zeros(ntu + ntv + nNavier * nip + (num.nLS + 1) * nip + nbp)

        uvD[1:ntu] .= ph.uD
        uvD[ntu+1:ntu+ntv] .= ph.vD 
        uvD[ntu+ntv+ntNavier+1:ntu+ntv+ntNavier+(num.nLS+1)*nip+nbp] .= ph.pD

        PDI_status = @ccall "libpdi".PDI_multi_expose("rhs_uv"::Cstring,
        "rhs_uv_len"::Cstring, length(uvD)::Ref{Clonglong}, PDI_OUT::Cint,
        "rhs_uv_1D"::Cstring, uvD::Ptr{Cdouble}, PDI_OUT::Cint,
        C_NULL::Ptr{Cvoid})::Cint
    
        # print("\n test Adummy\n",Adummy*uvD)

        # print("\n test Adummy\n",Adummy*uvD/gp.dx[1,1]^2)

        # print("\n test grid ",gp.dx[1,1])
        # print("\n factor ", factor )

        control_volumes = ones(ncol_A)

        # control_volumes[1:ntu] .= vec(gu.LS[1].geoL.dcap[:,:,5])
        # control_volumes[[ntu+1:ntu+ntv]] .= vec(gv.LS[1].geoL.dcap[:,:,5])
        # control_volumes[ntu+ntv+ntNavier+1:ntu+ntv+ntNavier+(num.nLS+1)*nip+nbp] .= vec(gp.LS[1].geoL.dcap[:,:,5])
        
        control_volumes[1:ni_u] .= vec(gu.LS[1].geoL.dcap[:,:,5]) #non-dimensionalise bulk rhs

        control_volumes[ntu+1:ntu+ni_v] .= vec(gv.LS[1].geoL.dcap[:,:,5]) #non-dimensionalise bulk rhs
        # control_volumes[ntu+ntv+ntNavier+1:ntu+ntv+ntNavier+(num.nLS+1)*nip+nbp] .= vec(gp.LS[1].geoL.dcap[:,:,5])

        #niu 

        # print("\n test Adummy\n",Adummy*uvD/factor/gp.dx[1,1]^2)

        # print("\n test volume\n",vec(gu.LS[1].geoL.dcap[:,:,5]))

        # print("\n test Adummy\n",Adummy*uvD/factor ./ control_volumes)

        

        test_matrix = zeros(ntu + ntv + nNavier * nip + (num.nLS + 1) * nip + nbp)

        test_matrix = Adummy*uvD/factor ./ control_volumes

        PDI_status = @ccall "libpdi".PDI_multi_expose("check_coupled_matrix"::Cstring,
        "vec_1D_len"::Cstring, length(test_matrix)::Ref{Clonglong}, PDI_OUT::Cint,
        "vec_1D"::Cstring, test_matrix::Ptr{Cdouble}, PDI_OUT::Cint,
        C_NULL::Ptr{Cvoid})::Cint

        #endregion check coupled matrix



        # Explicit gradient of pressure
        # ∇ϕ_x = opu.AxT * opu.Rx * vec1(pD,grid) .+ opu.Gx_b * vecb(pD,grid)
        # ∇ϕ_y = opv.AyT * opv.Ry * vec1(pD,grid) .+ opv.Gy_b * vecb(pD,grid)
        # for iLS in 1:nLS
        #     ∇ϕ_x .+= opu.Gx[iLS] * veci(pD,grid,iLS+1)
        #     ∇ϕ_y .+= opv.Gy[iLS] * veci(pD,grid,iLS+1)
        # end
        #endregion coupled pression

        #TODO sbu sbv interfacial_nb_iLS_pressure

        #region divergence of velocity: -div U for symmetry
        #equation written on lines ntu+ntv+1:ntu+ntv+nip 

        range_divergence = ntu+ntv+ntNavier+1:ntu+ntv+ntNavier+nip
        # u bulk
        A[range_divergence,bulk_u_velocity] = -opp.AxT 
        # u border
        A[range_divergence,border_u_velocity] = -opp.Gx_b
        # v bulk
        A[range_divergence,bulk_v_velocity] = -opp.AyT 
        # v border
        A[range_divergence,border_v_velocity] = -opp.Gy_b

        for iLS in 1:nLS
            sbu = iLS*niu+1:(iLS+1)*niu
            sbv = ntu+iLS*niv+1:ntu+(iLS+1)*niv
            # hypothesis no blowing so normal velocity null if Navier
            if !is_navier(bc_interface[iLS]) && !is_navier_cl(bc_interface[iLS]) 
                A[range_divergence,sbu] = -opp.Gx[iLS] 
                A[range_divergence,sbv] = -opp.Gy[iLS] 
            end
        end

        # divergence of velocity explicit
        # Duv = opp.AxT * vec1(u_predictionD,grid_u) .+ opp.Gx_b * vecb(u_predictionD,grid_u) .+
        #       opp.AyT * vec1(v_predictionD,grid_v) .+ opp.Gy_b * vecb(v_predictionD,grid_v)
        # for iLS in 1:nLS
        #     if !is_navier(bc_int[iLS]) && !is_navier_cl(bc_int[iLS])
        #         Duv .+= opp.Gx[iLS] * veci(u_predictionD,grid_u,iLS+1) .+ 
        #                 opp.Gy[iLS] * veci(v_predictionD,grid_v,iLS+1)
        #     end
        # end
        #endregion divergence of velocity: -div U for symmetry


        B[bulk_u_velocity,bulk_u_velocity] = Mum1
        B[bulk_v_velocity,bulk_v_velocity] = Mvm1
    end

    index_Navier = 0
    _iLS = 1
    for iLS in 1:num.nLS
        #TODO can be improved for readability / compilation:
        # allocates ones...
        # better tot do mapping function if cte... else if matrix ...

        range_nb_iLS_Navier = ntu+ntv+1+index_Navier*nip:ntu+ntv+(index_Navier+1)*nip

        #region BC iLS
        #for Navier: tangential velocity = \lambda grad tangential vel 
        a0u, a1u, bu, a0v, a1v, bv, a0p, bp = set_velocity_boundary_conditions(bc_interface, iLS, gu, gv, gp, num)
        #endregion BC iLS


        interfacial_nb_1_u_velocity = _iLS*niu+1:(_iLS+1)*niu
        interfacial_nb_1_v_velocity = ntu+_iLS*niv+1:ntu+(_iLS+1)*niv

        if ls_advection
            if !is_navier_cl(bc_interface[iLS]) && !is_navier(bc_interface[iLS])
                #region not Navier
                # Contribution to implicit part of viscous term from inner boundaries
                A[bulk_u_velocity,interfacial_nb_1_u_velocity] = - timestep_n .* bc_Lu[iLS]
                A[bulk_v_velocity,interfacial_nb_1_v_velocity] = - timestep_n .* bc_Lv[iLS] 
                # Boundary conditions for inner boundaries
                A[interfacial_nb_1_u_velocity,bulk_u_velocity] = bu * (opu.HxT[iLS] * opu.iMx * opu.Bx .+ opu.HyT[iLS] * opu.iMy * opu.By)
                A[interfacial_nb_1_v_velocity,bulk_v_velocity] = bv * (opv.HxT[iLS] * opv.iMx * opv.Bx .+ opv.HyT[iLS] * opv.iMy * opv.By)
                # Contribution to Neumann BC from other boundaries
                nNav2 = 0
                for i in 1:num.nLS
                    if i != iLS && (!is_navier_cl(bc_interface[i]) && !is_navier(bc_interface[i]))
                        A[interfacial_nb_1_u_velocity,i*niu+1:(i+1)*niu] = bu * (
                            opu.HxT[iLS] * opu.iMx * opu.Hx[i] .+
                            opu.HyT[iLS] * opu.iMy * opu.Hy[i]
                        )
                        A[interfacial_nb_1_v_velocity,ntu+i*niv+1:ntu+(i+1)*niv] = bv * (
                            opv.HxT[iLS] * opv.iMx * opv.Hx[i] .+
                            opv.HyT[iLS] * opv.iMy * opv.Hy[i]
                        )
                    elseif i != iLS
                        sinα = Diagonal(vec(sin.(gp.LS[i].α)))
                        # replace!(sinα.diag, NaN=>0.0)
                        cosα = Diagonal(vec(cos.(gp.LS[i].α)))
                        # replace!(cosα.diag, NaN=>0.0)

                        if any(isnan, sinα.diag) || any(isnan, cosα.diag)
                            @error("NaN FE_set_momentum_coupled")
                            replace!(sinα.diag, NaN=>0.0)
                            replace!(cosα.diag, NaN=>0.0)
                        end

                        # not Navier, averaging coefficients differently computed 
                        interpolate_x = interpolating_coefficient_Navier(gu,gp,i,bc_interface[iLS])
                        interpolate_y = interpolating_coefficient_Navier(gv,gp,i,bc_interface[iLS])

                        A[interfacial_nb_1_u_velocity,ntu+ntv+1+nNav2*nip:ntu+ntv+(nNav2+1)*nip] = bu * (
                            opu.HxT[iLS] * opu.iMx * opu.Hx[i] .+
                            opu.HyT[iLS] * opu.iMy * opu.Hy[i]
                        ) * interpolate_x * sinα
                        A[interfacial_nb_1_v_velocity,ntu+ntv+1+nNav2*nip:ntu+ntv+(nNav2+1)*nip] = bv * (
                            opv.HxT[iLS] * opv.iMx * opv.Hx[i] .+
                            opv.HyT[iLS] * opv.iMy * opv.Hy[i]
                        ) * interpolate_y * (-cosα)

                        nNav2 += 1
                    end
                end
                A[interfacial_nb_1_u_velocity,interfacial_nb_1_u_velocity] = pad(bu * (
                    opu.HxT[iLS] * opu.iMx * opu.Hx[iLS] .+
                    opu.HyT[iLS] * opu.iMy * opu.Hy[iLS]
                ) .- opu.χ[iLS] * a1u)
                A[interfacial_nb_1_v_velocity,interfacial_nb_1_v_velocity] = pad(bv * (
                    opv.HxT[iLS] * opv.iMx * opv.Hx[iLS] .+
                    opv.HyT[iLS] * opv.iMy * opv.Hy[iLS]
                ) .- opv.χ[iLS] * a1v)

                A[interfacial_nb_1_u_velocity,border_u_velocity] = bu * (
                    opu.HxT[iLS] * opu.iMx_b * opu.Hx_b .+ opu.HyT[iLS] * opu.iMy_b * opu.Hy_b
                )
                A[interfacial_nb_1_v_velocity,border_v_velocity] = bv * (
                    opv.HxT[iLS] * opv.iMx_b * opv.Hx_b .+ opv.HyT[iLS] * opv.iMy_b * opv.Hy_b
                )
                # Boundary conditions for outer boundaries
                A[border_u_velocity,interfacial_nb_1_u_velocity] = b_bu * (
                    opu.HxT_b * opu.iMx_b' * opu.Hx[iLS] .+ opu.HyT_b * opu.iMy_b' * opu.Hy[iLS]
                )
                A[border_v_velocity,interfacial_nb_1_v_velocity] = b_bv * (
                    opv.HxT_b * opv.iMx_b' * opv.Hx[iLS] .+ opv.HyT_b * opv.iMy_b' * opv.Hy[iLS]
                )

                _iLS += 1
                #endregion not Navier

            else 
                #region Navier
                # Tangential component of velocity if Navier BC #if !is_navier_cl(bc_interface[iLS]) && !is_navier(bc_interface[iLS])
                sinα_p = Diagonal(vec(sin.(gp.LS[iLS].α)))
                cosα_p = Diagonal(vec(cos.(gp.LS[iLS].α)))
                sinα_u = Diagonal(vec(sin.(gu.LS[iLS].α)))
                cosα_v = Diagonal(vec(cos.(gv.LS[iLS].α)))

                if any(isnan, sinα_p.diag) || any(isnan, cosα_p.diag) || any(isnan, sinα_u.diag) || any(isnan, cosα_v.diag)
                    @error("NaN FE_set_momentum_coupled")
                    replace!(sinα.diag, NaN=>0.0)
                    replace!(cosα.diag, NaN=>0.0)
                    replace!(sinα_u.diag, NaN=>0.0)
                    replace!(cosα_v.diag, NaN=>0.0)
                end

                # Contribution to implicit part of viscous term from inner boundaries
           
                # Navier BC, averaging coefficient computed differently
                interpolate_x = interpolating_coefficient_Navier(gu,gp,bc_interface[iLS])
                interpolate_y = interpolating_coefficient_Navier(gv,gp,bc_interface[iLS])
                # implicit part of viscous stress, from tangential velocity at wall with Navier, 
                # interpolated from scalar to u, v grids
                A[bulk_u_velocity,range_nb_iLS_Navier] = - iRe * timestep_n .* (
                    opu.BxT * opu.iMx * opu.Hx[iLS] .+
                    opu.ByT * opu.iMy * opu.Hy[iLS]
                ) * sinα_u * interpolate_x
                A[bulk_v_velocity,range_nb_iLS_Navier] = - iRe * timestep_n .* (
                    opv.BxT * opv.iMx * opv.Hx[iLS] .+
                    opv.ByT * opv.iMy * opv.Hy[iLS]
                ) * (-cosα_v) * interpolate_y

                # Boundary conditions for inner boundaries
                
                interpolate_u_to_p, interpolate_v_to_p = interpolating_coefficient_Navier_uv_grids_to_p_grid_volume(num,gp,gu,gv,iLS)

                A[range_nb_iLS_Navier,bulk_u_velocity] = bp * (
                    opp.HxT[iLS] * opp.iMx * opp.Bx .+
                    opp.HyT[iLS] * opp.iMy * opp.By
                ) * sinα_p * interpolate_u_to_p
                A[range_nb_iLS_Navier,bulk_v_velocity] = bp * (
                    opp.HxT[iLS] * opp.iMx * opp.Bx .+
                    opp.HyT[iLS] * opp.iMy * opp.By
                ) * (-cosα_p) * interpolate_v_to_p

                for i in 1:num.nLS
                    if i != iLS && (!is_navier_cl(bc_interface[i]) && !is_navier(bc_interface[i]))
                        #TODO why isn't it with volume like in interpolating_coefficient_Navier_uv_grids_to_p_grid_volume ?
                        interpolate_u_to_p, interpolate_v_to_p = interpolating_coefficient_Navier_uv_grids_to_p_grid_height(num,gp,gu,gv,i)

                        A[range_nb_iLS_Navier,i*niu+1:(i+1)*niu] = bp * (
                            opp.HxT[iLS] * opp.iMx * opp.Hx[i] .+
                            opp.HyT[iLS] * opp.iMy * opp.Hy[i]
                        ) * sinα_p * interpolate_u_to_p
                        A[range_nb_iLS_Navier,ntu+i*niv+1:ntu+(i+1)*niv] = bp * (
                            opp.HxT[iLS] * opp.iMx * opp.Hx[i] .+
                            opp.HyT[iLS] * opp.iMy * opp.Hy[i]
                        ) * (-cosα_p) * interpolate_v_to_p
                    end
                end
                A[range_nb_iLS_Navier,range_nb_iLS_Navier] = pad(bp * (
                    opp.HxT[iLS] * opp.iMx * opp.Hx[iLS] .+
                    opp.HyT[iLS] * opp.iMy * opp.Hy[iLS]
                ) .+ opp.χ[iLS])

                # coefficients (u,v) frame to (normal, tangent)
                sin_alpha_border = Diagonal(zeros(nbp))
                sin_alpha_border.diag[1:gp.ny] .= sin.(gp.LS[iLS].α[:,1])
                sin_alpha_border.diag[gp.ny+1:gp.ny+gp.nx] .= sin.(gp.LS[iLS].α[1,:])
                sin_alpha_border.diag[gp.ny+gp.nx+1:2gp.ny+gp.nx] .= sin.(gp.LS[iLS].α[:,end])
                sin_alpha_border.diag[2gp.ny+gp.nx+1:end] .= sin.(gp.LS[iLS].α[end,:])
                cos_alpha_border = Diagonal(zeros(nbp))
                cos_alpha_border.diag[1:gp.ny] .= cos.(gp.LS[iLS].α[:,1])
                cos_alpha_border.diag[gp.ny+1:gp.ny+gp.nx] .= cos.(gp.LS[iLS].α[1,:])
                cos_alpha_border.diag[gp.ny+gp.nx+1:2gp.ny+gp.nx] .= cos.(gp.LS[iLS].α[:,end])
                cos_alpha_border.diag[2gp.ny+gp.nx+1:end] .= cos.(gp.LS[iLS].α[end,:])

                #region interpolation coefficients
                interpolate_u_to_p = spdiagm(nbp, nbu, 0 => zeros(nbp), 1 => zeros(nbp-1))
                for ii in 1:gp.ny
                    interpolate_u_to_p[ii,ii] = 1.0
                    interpolate_u_to_p[gp.ny+gp.nx+ii,gu.ny+gu.nx+ii] = 1.0
                end
                for ii in 1:gp.nx
                    interpolate_u_to_p[ii+gp.ny,ii+gu.ny] = 0.5
                    interpolate_u_to_p[ii+gp.ny,ii+gu.ny+1] = 0.5
                    interpolate_u_to_p[ii+2gp.ny+gp.nx,ii+2gu.ny+gu.nx] = 0.5
                    interpolate_u_to_p[ii+2gp.ny+gp.nx,ii+2gu.ny+gu.nx+1] = 0.5
                end
                interpolate_v_to_p = spdiagm(nbp, nbv, 0 => zeros(nbp), 1 => zeros(nbp-1))
                for ii in 1:gp.ny
                    interpolate_v_to_p[ii,ii] = 0.5
                    interpolate_v_to_p[ii,ii+1] = 0.5
                    interpolate_v_to_p[gp.ny+gp.nx+ii,gv.ny+gv.nx+ii] = 0.5
                    interpolate_v_to_p[gp.ny+gp.nx+ii,gv.ny+gv.nx+ii+1] = 0.5
                end
                for ii in 1:gp.nx
                    interpolate_v_to_p[ii+gp.ny,ii+gv.ny] = 1.0
                    interpolate_v_to_p[ii+2gp.ny+gp.nx,ii+2gv.ny+gv.nx] = 1.0
                end
                #endregion interpolation coefficients

                #region Navier slip 
                A[range_nb_iLS_Navier,border_u_velocity] = bp * (
                    opp.HxT[iLS] * opp.iMx_b * opp.Hx_b .+ 
                    opp.HyT[iLS] * opp.iMy_b * opp.Hy_b
                ) * sin_alpha_border * interpolate_u_to_p
                A[range_nb_iLS_Navier,border_v_velocity] = bp * (
                    opp.HxT[iLS] * opp.iMx_b * opp.Hx_b .+ 
                    opp.HyT[iLS] * opp.iMy_b * opp.Hy_b
                ) * (-cos_alpha_border) * interpolate_v_to_p
                #endregion Navier slip
                
                # Boundary conditions for outer boundaries
                #region interpolation coefficients
                interpolate_u_to_p = spdiagm(nbu, nbp, 0 => zeros(nbp), 1 => zeros(nbp-1))
                for ii in 1:gu.ny
                    interpolate_u_to_p[ii,ii] = 1.0
                    interpolate_u_to_p[gu.ny+gu.nx+ii,gp.ny+gp.nx+ii] = 1.0
                end
                interpolate_u_to_p[gu.ny+1,gp.ny+1]  = 1.0
                interpolate_u_to_p[gu.ny+gu.nx,gp.ny+gp.nx]  = 1.0
                interpolate_u_to_p[2gu.ny+gu.nx+1,2gp.ny+gp.nx+1]  = 1.0
                interpolate_u_to_p[end,end]  = 1.0
                for ii in 2:(gu.nx-1)
                    interpolate_u_to_p[ii+gu.ny,ii+gp.ny-1] = 0.5
                    interpolate_u_to_p[ii+gu.ny,ii+gp.ny] = 0.5
                    interpolate_u_to_p[ii+2gu.ny+gu.nx,ii+2gp.ny+gp.nx-1] = 0.5
                    interpolate_u_to_p[ii+2gu.ny+gu.nx,ii+2gp.ny+gp.nx] = 0.5
                end
                interpolate_v_to_p = spdiagm(nbv, nbp, 0 => zeros(nbp), 1 => zeros(nbp-1))
                interpolate_v_to_p[1,1]  = 1.0
                interpolate_v_to_p[gv.ny,gp.ny]  = 1.0
                interpolate_v_to_p[gv.ny+gv.nx+1,gp.ny+gp.nx+1]  = 1.0
                interpolate_v_to_p[2gv.ny+gv.nx,2gp.ny+gp.nx]  = 1.0
                for ii in 2:(gv.ny-1)
                    interpolate_v_to_p[ii,ii-1] = 0.5
                    interpolate_v_to_p[ii,ii] = 0.5
                    interpolate_v_to_p[gv.ny+gv.nx+ii,gp.ny+gp.nx+ii] = 0.5
                    interpolate_v_to_p[gv.ny+gv.nx+ii,gp.ny+gp.nx+ii] = 0.5
                end
                for ii in 1:gv.nx
                    interpolate_v_to_p[ii+gv.ny,ii+gp.ny] = 1.0
                    interpolate_v_to_p[ii+2gv.ny+gv.nx,ii+2gp.ny+gp.nx] = 1.0
                end
                #endregion interpolation coefficients

                A[border_u_velocity,range_nb_iLS_Navier] = b_bu * (
                    opu.HxT_b * opu.iMx_b' * opu.Hx[iLS] .+ opu.HyT_b * opu.iMy_b' * opu.Hy[iLS]
                ) * interpolate_x * sinα_p
                A[border_v_velocity,range_nb_iLS_Navier] = b_bv * (
                    opv.HxT_b * opv.iMx_b' * opv.Hx[iLS] .+ opv.HyT_b * opv.iMy_b' * opv.Hy[iLS]
                ) * interpolate_y * (-cosα_p)
                
                #endregion Navier

            end # if !is_navier_cl(bc_interface[iLS]) && !is_navier(bc_interface[iLS])
        end #advection

        if !is_navier_cl(bc_interface[iLS]) && !is_navier(bc_interface[iLS])
            @inbounds rhs[interfacial_nb_1_u_velocity] .= opu.χ[iLS] * vec(a0u)
            @inbounds rhs[interfacial_nb_1_v_velocity] .= opv.χ[iLS] * vec(a0v)
        else
            @inbounds rhs[range_nb_iLS_Navier] .= opp.χ[iLS] * vec(a0p)
            index_Navier += 1
        end #!is_navier_cl(bc_interface[iLS]) && !is_navier(bc_interface[iLS])

    end #for iLS in 1:num.nLS

    
    @inbounds rhs[border_u_velocity] .= opu.χ_b * vec(a0_bu) #BC for u component on borders
    @inbounds rhs[border_v_velocity] .= opv.χ_b * vec(a0_bv) #BC for v component on borders

    if num.pressure_velocity_coupling == 2
        @inbounds rhs[border_pressure] .= opp.χ_b * vec(a0_bp) #BC for u component on borders
    end


    #region set first cells to boundary velocity
    if num.pressure_velocity_coupling == 3
        set_first_cells!(A,rhs,gu,ntu-nbu,0,true,false,true,false)
        set_first_cells!(A,rhs,gv,ntu+ntv-nbv,ntu,false,true,false,true)
    end
    #endregion set first cells to boundary velocity



    #rhs div u : 0 for Incompressible

    return rhs
end #end FE_set_momentum_coupled2


"""
Set the system matrix for Forward-Euler scheme

* bc_interface: interface boundary condition
* B contains Mum1 and Mvm1 (cell volumes)

## Modifications to `rhs`

1. **For each interface `iLS`:**
   - If the boundary condition is not Navier and not Navier-CL:
     ```julia
     @inbounds rhs[sbu] .= opu.χ[iLS] * vec(a0u)
     @inbounds rhs[sbv] .= opv.χ[iLS] * vec(a0v)
     ```
   - Otherwise (if the boundary condition is Navier or Navier-CL):
     ```julia
     @inbounds rhs[ntu+ntv+1+nNav1*nip:ntu+ntv+(nNav1+1)*nip] .= opp.χ[iLS] * vec(a0p)
     nNav1 += 1
     ```

2. **For the outer boundaries/borders:**
    ```julia
    @inbounds rhs[ntu-nbu+1:ntu] .= opu.χ_b * vec(a0_bu) #BC for u component on borders
    @inbounds rhs[ntu+ntv-nbv+1:ntu+ntv] .= opv.χ_b * vec(a0_bv) #BC for v component on borders
    ```
robin BC : source term a0
                
At the moment, the Levelset is not computed at borders/interfaces ? Only bulk
"""
function FE_set_momentum_coupled_two_phases(
    bc_interface, num, gp, gu, gv,
    op,
    A, B,
    rhs,
    Lu, bc_Lu, bc_Lu_b, Mum1, BCu,
    Lv, bc_Lv, bc_Lv_b, Mvm1, BCv,
    Lu_S, bc_Lu_S, bc_Lu_b_S, Mum1_S, BCu_S,
    Lv_S, bc_Lv_S, bc_Lv_b_S, Mvm1_S, BCv_S,
    ls_advection,BCp,ph=nothing,phS=nothing,
    )
    @unpack timestep_n, Re, nLS, nNavier = num


    printstyled(color=:red, @sprintf "\n coupled pressure-velocity FE_set_momentum_coupled_two_phases\n")


    #liquid phase 
    opp, opu, opv = op.opC_pL, op.opC_uL, op.opC_vL


    #region pressure eq for interface: jump

    if jump_mass_transfer_rate
        for iLS in 1:nLS
            if is_fs(bc_int[iLS])
                Smat = strain_rate(iLS, opC_u, opC_v, opC_p)
                S = Smat[1,1] * vec1(u_predictionD,grid_u) .+ Smat[1,2] * veci(u_predictionD,grid_u,iLS+1) .+
                    Smat[2,1] * vec1(v_predictionD,grid_v) .+ Smat[2,2] * veci(v_predictionD,grid_v,iLS+1)
    
                fs_mat = opC_p.HxT[iLS] * opC_p.Hx[iLS] .+ opC_p.HyT[iLS] * opC_p.Hy[iLS]
                veci(rhs_ϕ,grid,iLS+1) .= -2.0 .* mu1_over_rho1 .* S .+ Diagonal(diag(fs_mat)) * ( σ .* vec(grid.LS[iLS].κ) .- pres_free_suface .- diff_inv_rho * mass_transfer_rate ^ 2)
            end
        end
    else
        for iLS in 1:nLS
            if is_fs(bc_int[iLS])
                Smat = strain_rate(iLS, opC_u, opC_v, opC_p)
                S = Smat[1,1] * vec1(u_predictionD,grid_u) .+ Smat[1,2] * veci(u_predictionD,grid_u,iLS+1) .+
                    Smat[2,1] * vec1(v_predictionD,grid_v) .+ Smat[2,2] * veci(v_predictionD,grid_v,iLS+1)

                fs_mat = opC_p.HxT[iLS] * opC_p.Hx[iLS] .+ opC_p.HyT[iLS] * opC_p.Hy[iLS]
                veci(rhs_ϕ,grid,iLS+1) .= -2.0 .* mu1_over_rho1 .* S .+ Diagonal(diag(fs_mat)) * ( σ .* vec(grid.LS[iLS].κ) .- pres_free_suface )
            end
        end
    end

    #endregion pressure eq for interface: jump

    #region init
    iRe = num.visc_coeff

    nip = gp.nx * gp.ny
    nbp = 2 * gp.nx + 2 * gp.ny

    niu = gu.nx * gu.ny
    nbu = 2 * gu.nx + 2 * gu.ny
    ntu = (nLS - nNavier + 1) * niu + nbu

    niv = gv.nx * gv.ny
    nbv = 2 * gv.nx + 2 * gv.ny
    ntv = (nLS - nNavier + 1) * niv + nbv

    #Reset to zero
    rhs .= 0.0 

    #region BC borders u
    a0_bu = zeros(nbu)
    _a1_bu = zeros(nbu)
    _b_bu = zeros(nbu)
    for iLS in 1:num.nLS
        set_borders!(gu, gu.LS[iLS].cl, gu.LS[iLS].u, a0_bu, _a1_bu, _b_bu, BCu, num.n_ext_cl)
    end
    a1_bu = Diagonal(vec(_a1_bu))
    b_bu = Diagonal(vec(_b_bu))
    #endregion BC borders u

    #region BC borders v
    a0_bv = zeros(nbv)
    _a1_bv = zeros(nbv)
    _b_bv = zeros(nbv)
    for iLS in 1:num.nLS
        set_borders!(gv, gv.LS[iLS].cl, gv.LS[iLS].u, a0_bv, _a1_bv, _b_bv, BCv, num.n_ext_cl)
    end
    a1_bv = Diagonal(vec(_a1_bv))
    b_bv = Diagonal(vec(_b_bv))
    #endregion BC borders v


    #region BC borders p
    if num.pressure_velocity_coupling == 2
        # rhs = fnzeros(gp, num)

        a0_bp = zeros(nbp)
        _a1_bp = zeros(nbp)
        _b_bp = zeros(nbp)
        for iLS in 1:num.nLS
            set_borders_poisson!(gp, gp.LS[iLS].cl, gp.LS[iLS].u, a0_bp, _a1_bp, _b_bp, BCp, num.n_ext_cl)
        end
        a1_bp = Diagonal(vec(_a1_bp))
        b_bp = Diagonal(vec(_b_bp))
    end

    # if ls_advection
    #     # Poisson equation
    #     # A[1:ni,1:ni] = pad(L, -4.0)
    #     # A[1:ni,end-nb+1:end] = bc_L_b

    #     # Boundary conditions for outer boundaries
    #     A[end-nb+1:end,1:ni] = b_bp * (HxT_b * iMx_b' * Bx .+ HyT_b * iMy_b' * By) 
    #     A[end-nb+1:end,end-nb+1:end] = pad(b_bp * (HxT_b * iMx_bd * Hx_b .+ HyT_b * iMy_bd * Hy_b) .+ χ_b * a1_bp,- 4.0)
    # end

    # for iLS in 1:num.nLS
    #     if ls_advection
    #         if is_dirichlet(bc_type[iLS])
    #             __a1 = 1.0
    #             __a2 = 0.0
    #             __b = 0.0
    #         elseif is_neumann(bc_type[iLS])
    #             __a1 = 0.0
    #             __a2 = 0.0
    #             __b = 1.0
    #         elseif is_robin(bc_type[iLS])
    #             __a1 = 1.0
    #             __a2 = 0.0
    #             __b = 1.0
    #         elseif is_fs(bc_type[iLS])
    #             __a1 = 0.0
    #             __a2 = 1.0
    #             __b = 0.0
    #         elseif is_wall_no_slip(bc_type[iLS])
    #             __a1 = 0.0
    #             __a2 = 0.0
    #             __b = 1.0
    #         elseif is_navier(bc_type[iLS])
    #             __a1 = 0.0
    #             __a2 = 0.0
    #             __b = 1.0
    #         elseif is_navier_cl(bc_type[iLS])
    #             __a1 = 0.0
    #             __a2 = 0.0
    #             __b = 1.0
    #         else
    #             __a1 = 0.0
    #             __a2 = 0.0
    #             __b = 1.0
    #         end
    
    #         _a1 = ones(gp) .* __a1
    #         a1 = Diagonal(vec(_a1))
    #         _a2 = ones(gp) .* __a2
    #         a2 = Diagonal(vec(_a2))
    #         _b = ones(gp) .* __b
    #         b = Diagonal(vec(_b))

    #         fs_mat = HxT[iLS] * Hx[iLS] .+ HyT[iLS] * Hy[iLS]

    #         sb = iLS*ni+1:(iLS+1)*ni
            
    #         # Poisson equation
    #         A[1:ni,sb] = bc_L[iLS]
    #         # Boundary conditions for inner boundaries
    #         A[sb,1:ni] = b * (HxT[iLS] * iMx * Bx .+ HyT[iLS] * iMy * By)
    #         # Contribution to Neumann BC from other boundaries
    #         for i in 1:num.nLS
    #             if i != iLS
    #                 A[sb,i*ni+1:(i+1)*ni] = b * (HxT[iLS] * iMx * Hx[i] .+ HyT[iLS] * iMy * Hy[i])
    #             end
    #         end
    #         A[sb,sb] = pad(
    #             b *(HxT[iLS] * iMx * Hx[iLS] .+ HyT[iLS] * iMy * Hy[iLS]) .+ χ[iLS] * a1 .-
    #             a2 * Diagonal(diag(fs_mat)), -4.0
    #         )
    #         #TODO was +b... here in set_poisson 
    #         A[sb,end-nb+1:end] = b * (HxT[iLS] * iMx_b * Hx_b .+ HyT[iLS] * iMy_b * Hy_b) #why was it different in set_poisson ?
    #         # Boundary conditions for outer boundaries
    #         A[end-nb+1:end,sb] = b_b * (HxT_b * iMx_b' * Hx[iLS] .+ HyT_b * iMy_b' * Hy[iLS])
    #     end

    #     veci(rhs,gp,iLS+1) .= +χ[iLS] * vec(a0[iLS]) #was - in set_poisson
    # end

    # vecb(rhs,gp) .= +χ_b * vec(a0_b) #was - in set_poisson
    #endregion BC borders p


    #endregion init

    # Array indices 
    bulk_u_velocity = 1:niu
    bulk_v_velocity = ntu+1:ntu+niv
    # bulk_tangential_velocity = 

    border_u_velocity = ntu-nbu+1:ntu
    border_v_velocity = ntu+ntv-nbv+1:ntu+ntv

    # nt = (num.nLS - num.nNavier + 1) * ni_uv + num.nNavier * nip + nb_uv + (num.nLS + 1) * nip + nbp
    ntNavier = num.nNavier * nip
    bulk_pressure = ntu+ntv+ntNavier+1:ntu+ntv+ntNavier+nip
    #ntu+ntv+1:ntu+ntv+nip
    border_pressure = ntu+ntv+ntNavier+(num.nLS + 1)*nip+1:ntu+ntv+ntNavier+(num.nLS + 1)*nip+nbp

    print("\n len rhs_uv ",size(rhs))

    print("bulk_pressure ",bulk_pressure)
    print("border_pressure ",border_pressure)


    if ls_advection
        A.nzval .= 0.0
        # Implicit part of viscous term
        A[bulk_u_velocity,bulk_u_velocity] = pad_crank_nicolson(opu.M .- timestep_n .* Lu, gu, timestep_n)
        # Contribution to implicit part of viscous term from outer boundaries
        A[bulk_u_velocity,border_u_velocity] = - timestep_n .* bc_Lu_b

        # Boundary conditions for outer boundaries
        A[border_u_velocity,bulk_u_velocity] = b_bu * (opu.HxT_b * opu.iMx_b' * opu.Bx .+ opu.HyT_b * opu.iMy_b' * opu.By)
        A[border_u_velocity,border_u_velocity] = pad(b_bu * (
            opu.HxT_b * opu.iMx_bd * opu.Hx_b .+ 
            opu.HyT_b * opu.iMy_bd * opu.Hy_b
        ) .- opu.χ_b * a1_bu)

        # Implicit part of viscous term
        A[bulk_v_velocity,bulk_v_velocity] = pad_crank_nicolson(opv.M .- timestep_n .* Lv, gv, timestep_n)
        # Contribution to implicit part of viscous term from outer boundaries
        A[bulk_v_velocity,border_v_velocity] = - timestep_n .* bc_Lv_b

        
        # Boundary conditions for outer boundaries
        A[border_v_velocity,bulk_v_velocity] = b_bv * (opv.HxT_b * opv.iMx_b' * opv.Bx .+ opv.HyT_b * opv.iMy_b' * opv.By)
        A[border_v_velocity,border_v_velocity] = pad(b_bv * (
            opv.HxT_b * opv.iMx_bd * opv.Hx_b .+ 
            opv.HyT_b * opv.iMy_bd * opv.Hy_b
        ) .- opv.χ_b * a1_bv)

        # TODO pad 1 or -4
        #TODO sign divergence not same u v and p

        #region coupled pression
        
        if num.pressure_velocity_coupling == 2

            # Poisson equation
            # A[1:ni,1:ni] = pad(L, -4.0)
            # A[1:ni,end-nb+1:end] = bc_L_b
    
            # Boundary conditions for outer boundaries
            A[border_pressure,bulk_pressure] = b_bp * (opp.HxT_b * opp.iMx_b' * opp.Bx .+ opp.HyT_b * opp.iMy_b' * opp.By) 
            A[border_pressure,border_pressure] = pad(b_bp * (
                opp.HxT_b * opp.iMx_bd * opp.Hx_b .+ 
                opp.HyT_b * opp.iMy_bd * opp.Hy_b) .+ opp.χ_b * a1_bp) #-4.0
            
        end
        
        # Implicit gradient of pressure (volume integrated)

        # timestep_n .* irho1 .*
        irho1 = 1.0/num.rho1
        factor = num.timestep_n * irho1
        A[bulk_u_velocity,bulk_pressure] = factor * opu.AxT * opu.Rx #TODO + or multiply by cell volume ? + not required, only component of matrix
        A[bulk_v_velocity,bulk_pressure] = factor * opv.AyT * opv.Ry 
       
        if num.pressure_velocity_coupling !=3
            #Outer boundaries
            A[bulk_u_velocity,border_pressure] = factor * opu.Gx_b
            A[bulk_v_velocity,border_pressure] = factor * opv.Gy_b
        

            for iLS in 1:nLS
                interfacial_nb_iLS_pressure = ntu+ntv+nNavier*nip+nip*iLS+1:ntu+ntv+nNavier*nip+(iLS+1)*nip
                A[bulk_u_velocity,interfacial_nb_iLS_pressure] .+= factor * opu.Gx[iLS] 
                A[bulk_v_velocity,interfacial_nb_iLS_pressure] .+= factor * opv.Gy[iLS] 
            end
        
        end #num.pressure_velocity_coupling !=3

        #region check coupled matrix
        #Test pressure part
        # Adummy = copy(A)
        # Adummy.nzval .= 0.0

        # nt = ntu+ntv + nNavier * nip + (num.nLS + 1) * nip + nbp
        # ncol_A = ntu+ntv + nNavier * nip + nip

        ni_p = gp.nx * gp.ny
        nb_p = 2 * gp.nx + 2 * gp.ny
    
        ni_u = gu.nx * gu.ny
        nb_u = 2 * gu.nx + 2 * gu.ny
    
        ni_v = gv.nx * gv.ny
        nb_v = 2 * gv.nx + 2 * gv.ny

        ni_uv = ni_u + ni_v
        nb_uv = nb_u + nb_v

        nt = (num.nLS - num.nNavier + 1) * ni_uv + num.nNavier * ni_p + nb_uv + (num.nLS + 1) * ni_p + nb_p
        ncol_A = (num.nLS - num.nNavier + 1) * ni_uv + num.nNavier * ni_p + nb_uv + ni_p

        Adummy = spzeros(ncol_A, nt)


        # length_colptr = length(Adummy.colptr)
        # length_rowval = length(Adummy.rowval)
        # length_nzval = length(Adummy.nzval)

        # PDI_status = @ccall "libpdi".PDI_multi_expose("print_matrix"::Cstring,
        # "Auv_n"::Cstring, Adummy.n::Ref{Clonglong}, PDI_OUT::Cint,
        # "Auv_m"::Cstring, Adummy.m::Ref{Clonglong}, PDI_OUT::Cint,
        # "Auv_colptr_len"::Cstring, length_colptr::Ref{Clonglong}, PDI_OUT::Cint,
        # "Auv_rowval_len"::Cstring, length_rowval::Ref{Clonglong}, PDI_OUT::Cint,
        # "Auv_nzval_len"::Cstring, length_nzval::Ref{Clonglong}, PDI_OUT::Cint,
        # # "Auv_colptr_len"::Cstring, length(Adummy.colptr)::Ref{Clonglong}, PDI_OUT::Cint,
        # # "Auv_rowval_len"::Cstring, length(Adummy.rowval)::Ref{Clonglong}, PDI_OUT::Cint,
        # # "Auv_nzval_len"::Cstring, length(Adummy.nzval)::Ref{Clonglong}, PDI_OUT::Cint,
        # "Auv_colptr_1D"::Cstring, Adummy.colptr::Ptr{Clonglong}, PDI_OUT::Cint,
        # "Auv_rowval_1D"::Cstring, Adummy.rowval::Ptr{Clonglong}, PDI_OUT::Cint,
        # "Auv_nzval_1D"::Cstring, Adummy.nzval::Ptr{Cdouble}, PDI_OUT::Cint,
        # C_NULL::Ptr{Cvoid})::Cint

        irho1 = 1.0/num.rho1
        factor = num.timestep_n * irho1
        Adummy[bulk_u_velocity,bulk_pressure] = factor * opu.AxT * opu.Rx #TODO + or multiply by cell volume ? + not required, only component of matrix
        Adummy[bulk_v_velocity,bulk_pressure] = factor * opv.AyT * opv.Ry 
        #Outer boundaries
        Adummy[bulk_u_velocity,border_pressure] = factor * opu.Gx_b
        Adummy[bulk_v_velocity,border_pressure] = factor * opv.Gy_b


        for iLS in 1:nLS
            interfacial_nb_iLS_pressure = ntu+ntv+nNavier*nip+nip*iLS+1:ntu+ntv+nNavier*nip+(iLS+1)*nip
            Adummy[bulk_u_velocity,interfacial_nb_iLS_pressure] .+= factor * opu.Gx[iLS] 
            Adummy[bulk_v_velocity,interfacial_nb_iLS_pressure] .+= factor * opv.Gy[iLS] 
        end

        # print("\n A before ")

        # # print("\n  Auv.n ", A.n)
        # # print("\n  Auv.m ", A.m)

        # # print("\n  length(A.colptr) ", length(A.colptr))
        # # print("\n  length(A.rowval) ", length(A.rowval))
        # # print("\n  length(A.nzval) ", length(A.nzval))

        # # print("\n  Auv.colptr ", A.colptr)
        # # print("\n  Auv.rowval ", A.rowval)
        # # print("\n  Auv.nzval ", A.nzval)

        # length_colptr = length(A.colptr)
        # length_rowval = length(A.rowval)
        # length_nzval = length(A.nzval)


        # # nx: $nx
        # # ny: $ny
        # # nNavier: $nb_Navier_slip_BC
        # # nb_levelsets: $nb_levelsets

        # PDI_status = @ccall "libpdi".PDI_multi_expose("print_matrix_test"::Cstring,
        # "Auv_n"::Cstring, A.n::Ref{Clonglong}, PDI_OUT::Cint,
        # "Auv_m"::Cstring, A.m::Ref{Clonglong}, PDI_OUT::Cint,
        # "nx"::Cstring, gp.nx::Ref{Clonglong}, PDI_OUT::Cint,
        # "ny"::Cstring, gp.ny::Ref{Clonglong}, PDI_OUT::Cint,
        # "nb_Navier_slip_BC"::Cstring, num.nNavier::Ref{Clonglong}, PDI_OUT::Cint,
        # # "Auv_colptr_len"::Cstring, length_colptr::Ref{Clonglong}, PDI_OUT::Cint,
        # # "Auv_rowval_len"::Cstring, length_rowval::Ref{Clonglong}, PDI_OUT::Cint,
        # # "Auv_nzval_len"::Cstring, length_nzval::Ref{Clonglong}, PDI_OUT::Cint,
        # # # "Auv_colptr_len"::Cstring, length(A.colptr)::Ref{Clonglong}, PDI_OUT::Cint,
        # # # "Auv_rowval_len"::Cstring, length(A.rowval)::Ref{Clonglong}, PDI_OUT::Cint,
        # # # "Auv_nzval_len"::Cstring, length(A.nzval)::Ref{Clonglong}, PDI_OUT::Cint,
        # # "Auv_colptr_1D"::Cstring, A.colptr::Ptr{Clonglong}, PDI_OUT::Cint,
        # # "Auv_rowval_1D"::Cstring, A.rowval::Ptr{Clonglong}, PDI_OUT::Cint,
        # # "Auv_nzval_1D"::Cstring, A.nzval::Ptr{Cdouble}, PDI_OUT::Cint,
        # C_NULL::Ptr{Cvoid})::Cint

        # PDI_status = @ccall "libpdi".PDI_multi_expose("print_matrix"::Cstring,
        # "Auv_n"::Cstring, A.n::Ref{Clonglong}, PDI_OUT::Cint,
        # "Auv_m"::Cstring, A.m::Ref{Clonglong}, PDI_OUT::Cint,
        # "Auv_colptr_len"::Cstring, length_colptr::Ref{Clonglong}, PDI_OUT::Cint,
        # "Auv_rowval_len"::Cstring, length_rowval::Ref{Clonglong}, PDI_OUT::Cint,
        # "Auv_nzval_len"::Cstring, length_nzval::Ref{Clonglong}, PDI_OUT::Cint,
        # # "Auv_colptr_len"::Cstring, length(A.colptr)::Ref{Clonglong}, PDI_OUT::Cint,
        # # "Auv_rowval_len"::Cstring, length(A.rowval)::Ref{Clonglong}, PDI_OUT::Cint,
        # # "Auv_nzval_len"::Cstring, length(A.nzval)::Ref{Clonglong}, PDI_OUT::Cint,
        # "Auv_colptr_1D"::Cstring, A.colptr::Ptr{Clonglong}, PDI_OUT::Cint,
        # "Auv_rowval_1D"::Cstring, A.rowval::Ptr{Clonglong}, PDI_OUT::Cint,
        # "Auv_nzval_1D"::Cstring, A.nzval::Ptr{Cdouble}, PDI_OUT::Cint,
        # C_NULL::Ptr{Cvoid})::Cint

        # print("\n Adummy ",Adummy)
        print("\n Adummy ")

        # print("\n  Auv.n ", Adummy.n)
        # print("\n  Auv.m ", Adummy.m)

        # print("\n  length(Adummy.colptr) ", length(Adummy.colptr))
        # print("\n  length(Adummy.rowval) ", length(Adummy.rowval))
        # print("\n  length(Adummy.nzval) ", length(Adummy.nzval))

        # print("\n  Auv.colptr ", Adummy.colptr)
        # print("\n  Auv.rowval ", Adummy.rowval)
        # print("\n  Auv.nzval ", Adummy.nzval)

      
        #region print grad p part of matrix
        PDI_status = @ccall "libpdi".PDI_multi_expose("print_matrix"::Cstring,
        "Auv_n"::Cstring, Adummy.n::Ref{Clonglong}, PDI_OUT::Cint,
        "Auv_m"::Cstring, Adummy.m::Ref{Clonglong}, PDI_OUT::Cint,
        # "Auv_colptr_len"::Cstring, length_colptr::Ref{Clonglong}, PDI_OUT::Cint,
        # "Auv_rowval_len"::Cstring, length_rowval::Ref{Clonglong}, PDI_OUT::Cint,
        # "Auv_nzval_len"::Cstring, length_nzval::Ref{Clonglong}, PDI_OUT::Cint,
        "Auv_colptr_len"::Cstring, length(Adummy.colptr)::Ref{Clonglong}, PDI_OUT::Cint,
        "Auv_rowval_len"::Cstring, length(Adummy.rowval)::Ref{Clonglong}, PDI_OUT::Cint,
        "Auv_nzval_len"::Cstring, length(Adummy.nzval)::Ref{Clonglong}, PDI_OUT::Cint,
        "Auv_colptr_1D"::Cstring, Adummy.colptr::Ptr{Clonglong}, PDI_OUT::Cint,
        "Auv_rowval_1D"::Cstring, Adummy.rowval::Ptr{Clonglong}, PDI_OUT::Cint,
        "Auv_nzval_1D"::Cstring, Adummy.nzval::Ptr{Cdouble}, PDI_OUT::Cint,
        C_NULL::Ptr{Cvoid})::Cint
        #endregion print grad p part of matrix
   


        uvD = zeros(ntu + ntv + nNavier * nip + (num.nLS + 1) * nip + nbp)

        uvD[1:ntu] .= ph.uD
        uvD[ntu+1:ntu+ntv] .= ph.vD 
        uvD[ntu+ntv+ntNavier+1:ntu+ntv+ntNavier+(num.nLS+1)*nip+nbp] .= ph.pD

        PDI_status = @ccall "libpdi".PDI_multi_expose("rhs_uv"::Cstring,
        "rhs_uv_len"::Cstring, length(uvD)::Ref{Clonglong}, PDI_OUT::Cint,
        "rhs_uv_1D"::Cstring, uvD::Ptr{Cdouble}, PDI_OUT::Cint,
        C_NULL::Ptr{Cvoid})::Cint
    
        # print("\n test Adummy\n",Adummy*uvD)

        # print("\n test Adummy\n",Adummy*uvD/gp.dx[1,1]^2)

        # print("\n test grid ",gp.dx[1,1])
        # print("\n factor ", factor )

        control_volumes = ones(ncol_A)

        # control_volumes[1:ntu] .= vec(gu.LS[1].geoL.dcap[:,:,5])
        # control_volumes[[ntu+1:ntu+ntv]] .= vec(gv.LS[1].geoL.dcap[:,:,5])
        # control_volumes[ntu+ntv+ntNavier+1:ntu+ntv+ntNavier+(num.nLS+1)*nip+nbp] .= vec(gp.LS[1].geoL.dcap[:,:,5])
        
        control_volumes[1:ni_u] .= vec(gu.LS[1].geoL.dcap[:,:,5]) #non-dimensionalise bulk rhs

        control_volumes[ntu+1:ntu+ni_v] .= vec(gv.LS[1].geoL.dcap[:,:,5]) #non-dimensionalise bulk rhs
        # control_volumes[ntu+ntv+ntNavier+1:ntu+ntv+ntNavier+(num.nLS+1)*nip+nbp] .= vec(gp.LS[1].geoL.dcap[:,:,5])

        #niu 

        # print("\n test Adummy\n",Adummy*uvD/factor/gp.dx[1,1]^2)

        # print("\n test volume\n",vec(gu.LS[1].geoL.dcap[:,:,5]))

        # print("\n test Adummy\n",Adummy*uvD/factor ./ control_volumes)

        

        test_matrix = zeros(ntu + ntv + nNavier * nip + (num.nLS + 1) * nip + nbp)

        test_matrix = Adummy*uvD/factor ./ control_volumes

        PDI_status = @ccall "libpdi".PDI_multi_expose("check_coupled_matrix"::Cstring,
        "vec_1D_len"::Cstring, length(test_matrix)::Ref{Clonglong}, PDI_OUT::Cint,
        "vec_1D"::Cstring, test_matrix::Ptr{Cdouble}, PDI_OUT::Cint,
        C_NULL::Ptr{Cvoid})::Cint

        #endregion check coupled matrix



        # Explicit gradient of pressure
        # ∇ϕ_x = opu.AxT * opu.Rx * vec1(pD,grid) .+ opu.Gx_b * vecb(pD,grid)
        # ∇ϕ_y = opv.AyT * opv.Ry * vec1(pD,grid) .+ opv.Gy_b * vecb(pD,grid)
        # for iLS in 1:nLS
        #     ∇ϕ_x .+= opu.Gx[iLS] * veci(pD,grid,iLS+1)
        #     ∇ϕ_y .+= opv.Gy[iLS] * veci(pD,grid,iLS+1)
        # end
        #endregion coupled pression

        #TODO sbu sbv interfacial_nb_iLS_pressure

        #region divergence of velocity: -div U for symmetry
        #equation written on lines ntu+ntv+1:ntu+ntv+nip 

        range_divergence = ntu+ntv+ntNavier+1:ntu+ntv+ntNavier+nip
        # u bulk
        A[range_divergence,bulk_u_velocity] = -opp.AxT 
        # u border
        A[range_divergence,border_u_velocity] = -opp.Gx_b
        # v bulk
        A[range_divergence,bulk_v_velocity] = -opp.AyT 
        # v border
        A[range_divergence,border_v_velocity] = -opp.Gy_b

        for iLS in 1:nLS
            sbu = iLS*niu+1:(iLS+1)*niu
            sbv = ntu+iLS*niv+1:ntu+(iLS+1)*niv
            # hypothesis no blowing so normal velocity null if Navier
            if !is_navier(bc_interface[iLS]) && !is_navier_cl(bc_interface[iLS]) 
                A[range_divergence,sbu] = -opp.Gx[iLS] 
                A[range_divergence,sbv] = -opp.Gy[iLS] 
            end
        end

        # divergence of velocity explicit
        # Duv = opp.AxT * vec1(u_predictionD,grid_u) .+ opp.Gx_b * vecb(u_predictionD,grid_u) .+
        #       opp.AyT * vec1(v_predictionD,grid_v) .+ opp.Gy_b * vecb(v_predictionD,grid_v)
        # for iLS in 1:nLS
        #     if !is_navier(bc_int[iLS]) && !is_navier_cl(bc_int[iLS])
        #         Duv .+= opp.Gx[iLS] * veci(u_predictionD,grid_u,iLS+1) .+ 
        #                 opp.Gy[iLS] * veci(v_predictionD,grid_v,iLS+1)
        #     end
        # end
        #endregion divergence of velocity: -div U for symmetry


        B[bulk_u_velocity,bulk_u_velocity] = Mum1
        B[bulk_v_velocity,bulk_v_velocity] = Mvm1
    end

    index_Navier = 0
    _iLS = 1
    for iLS in 1:num.nLS
        #TODO can be improved for readability / compilation:
        # allocates ones...
        # better tot do mapping function if cte... else if matrix ...

        range_nb_iLS_Navier = ntu+ntv+1+index_Navier*nip:ntu+ntv+(index_Navier+1)*nip

        #region BC iLS
        #for Navier: tangential velocity = \lambda grad tangential vel 
        a0u, a1u, bu, a0v, a1v, bv, a0p, bp = set_velocity_boundary_conditions(bc_interface, iLS, gu, gv, gp, num)
        #endregion BC iLS


        interfacial_nb_1_u_velocity = _iLS*niu+1:(_iLS+1)*niu
        interfacial_nb_1_v_velocity = ntu+_iLS*niv+1:ntu+(_iLS+1)*niv

        if ls_advection
            if !is_navier_cl(bc_interface[iLS]) && !is_navier(bc_interface[iLS])
                #region not Navier
                # Contribution to implicit part of viscous term from inner boundaries
                A[bulk_u_velocity,interfacial_nb_1_u_velocity] = - timestep_n .* bc_Lu[iLS]
                A[bulk_v_velocity,interfacial_nb_1_v_velocity] = - timestep_n .* bc_Lv[iLS] 
                # Boundary conditions for inner boundaries
                A[interfacial_nb_1_u_velocity,bulk_u_velocity] = bu * (opu.HxT[iLS] * opu.iMx * opu.Bx .+ opu.HyT[iLS] * opu.iMy * opu.By)
                A[interfacial_nb_1_v_velocity,bulk_v_velocity] = bv * (opv.HxT[iLS] * opv.iMx * opv.Bx .+ opv.HyT[iLS] * opv.iMy * opv.By)
                # Contribution to Neumann BC from other boundaries
                nNav2 = 0
                for i in 1:num.nLS
                    if i != iLS && (!is_navier_cl(bc_interface[i]) && !is_navier(bc_interface[i]))
                        A[interfacial_nb_1_u_velocity,i*niu+1:(i+1)*niu] = bu * (
                            opu.HxT[iLS] * opu.iMx * opu.Hx[i] .+
                            opu.HyT[iLS] * opu.iMy * opu.Hy[i]
                        )
                        A[interfacial_nb_1_v_velocity,ntu+i*niv+1:ntu+(i+1)*niv] = bv * (
                            opv.HxT[iLS] * opv.iMx * opv.Hx[i] .+
                            opv.HyT[iLS] * opv.iMy * opv.Hy[i]
                        )
                    elseif i != iLS
                        sinα = Diagonal(vec(sin.(gp.LS[i].α)))
                        # replace!(sinα.diag, NaN=>0.0)
                        cosα = Diagonal(vec(cos.(gp.LS[i].α)))
                        # replace!(cosα.diag, NaN=>0.0)

                        if any(isnan, sinα.diag) || any(isnan, cosα.diag)
                            @error("NaN FE_set_momentum_coupled")
                            replace!(sinα.diag, NaN=>0.0)
                            replace!(cosα.diag, NaN=>0.0)
                        end

                        # not Navier, averaging coefficients differently computed 
                        interpolate_x = interpolating_coefficient_Navier(gu,gp,i,bc_interface[iLS])
                        interpolate_y = interpolating_coefficient_Navier(gv,gp,i,bc_interface[iLS])

                        A[interfacial_nb_1_u_velocity,ntu+ntv+1+nNav2*nip:ntu+ntv+(nNav2+1)*nip] = bu * (
                            opu.HxT[iLS] * opu.iMx * opu.Hx[i] .+
                            opu.HyT[iLS] * opu.iMy * opu.Hy[i]
                        ) * interpolate_x * sinα
                        A[interfacial_nb_1_v_velocity,ntu+ntv+1+nNav2*nip:ntu+ntv+(nNav2+1)*nip] = bv * (
                            opv.HxT[iLS] * opv.iMx * opv.Hx[i] .+
                            opv.HyT[iLS] * opv.iMy * opv.Hy[i]
                        ) * interpolate_y * (-cosα)

                        nNav2 += 1
                    end
                end
                A[interfacial_nb_1_u_velocity,interfacial_nb_1_u_velocity] = pad(bu * (
                    opu.HxT[iLS] * opu.iMx * opu.Hx[iLS] .+
                    opu.HyT[iLS] * opu.iMy * opu.Hy[iLS]
                ) .- opu.χ[iLS] * a1u)
                A[interfacial_nb_1_v_velocity,interfacial_nb_1_v_velocity] = pad(bv * (
                    opv.HxT[iLS] * opv.iMx * opv.Hx[iLS] .+
                    opv.HyT[iLS] * opv.iMy * opv.Hy[iLS]
                ) .- opv.χ[iLS] * a1v)

                A[interfacial_nb_1_u_velocity,border_u_velocity] = bu * (
                    opu.HxT[iLS] * opu.iMx_b * opu.Hx_b .+ opu.HyT[iLS] * opu.iMy_b * opu.Hy_b
                )
                A[interfacial_nb_1_v_velocity,border_v_velocity] = bv * (
                    opv.HxT[iLS] * opv.iMx_b * opv.Hx_b .+ opv.HyT[iLS] * opv.iMy_b * opv.Hy_b
                )
                # Boundary conditions for outer boundaries
                A[border_u_velocity,interfacial_nb_1_u_velocity] = b_bu * (
                    opu.HxT_b * opu.iMx_b' * opu.Hx[iLS] .+ opu.HyT_b * opu.iMy_b' * opu.Hy[iLS]
                )
                A[border_v_velocity,interfacial_nb_1_v_velocity] = b_bv * (
                    opv.HxT_b * opv.iMx_b' * opv.Hx[iLS] .+ opv.HyT_b * opv.iMy_b' * opv.Hy[iLS]
                )

                _iLS += 1
                #endregion not Navier

            else 
                #region Navier
                # Tangential component of velocity if Navier BC #if !is_navier_cl(bc_interface[iLS]) && !is_navier(bc_interface[iLS])
                sinα_p = Diagonal(vec(sin.(gp.LS[iLS].α)))
                cosα_p = Diagonal(vec(cos.(gp.LS[iLS].α)))
                sinα_u = Diagonal(vec(sin.(gu.LS[iLS].α)))
                cosα_v = Diagonal(vec(cos.(gv.LS[iLS].α)))

                if any(isnan, sinα_p.diag) || any(isnan, cosα_p.diag) || any(isnan, sinα_u.diag) || any(isnan, cosα_v.diag)
                    @error("NaN FE_set_momentum_coupled")
                    replace!(sinα.diag, NaN=>0.0)
                    replace!(cosα.diag, NaN=>0.0)
                    replace!(sinα_u.diag, NaN=>0.0)
                    replace!(cosα_v.diag, NaN=>0.0)
                end

                # Contribution to implicit part of viscous term from inner boundaries
           
                # Navier BC, averaging coefficient computed differently
                interpolate_x = interpolating_coefficient_Navier(gu,gp,bc_interface[iLS])
                interpolate_y = interpolating_coefficient_Navier(gv,gp,bc_interface[iLS])
                # implicit part of viscous stress, from tangential velocity at wall with Navier, 
                # interpolated from scalar to u, v grids
                A[bulk_u_velocity,range_nb_iLS_Navier] = - iRe * timestep_n .* (
                    opu.BxT * opu.iMx * opu.Hx[iLS] .+
                    opu.ByT * opu.iMy * opu.Hy[iLS]
                ) * sinα_u * interpolate_x
                A[bulk_v_velocity,range_nb_iLS_Navier] = - iRe * timestep_n .* (
                    opv.BxT * opv.iMx * opv.Hx[iLS] .+
                    opv.ByT * opv.iMy * opv.Hy[iLS]
                ) * (-cosα_v) * interpolate_y

                # Boundary conditions for inner boundaries
                
                interpolate_u_to_p, interpolate_v_to_p = interpolating_coefficient_Navier_uv_grids_to_p_grid_volume(num,gp,gu,gv,iLS)

                A[range_nb_iLS_Navier,bulk_u_velocity] = bp * (
                    opp.HxT[iLS] * opp.iMx * opp.Bx .+
                    opp.HyT[iLS] * opp.iMy * opp.By
                ) * sinα_p * interpolate_u_to_p
                A[range_nb_iLS_Navier,bulk_v_velocity] = bp * (
                    opp.HxT[iLS] * opp.iMx * opp.Bx .+
                    opp.HyT[iLS] * opp.iMy * opp.By
                ) * (-cosα_p) * interpolate_v_to_p

                for i in 1:num.nLS
                    if i != iLS && (!is_navier_cl(bc_interface[i]) && !is_navier(bc_interface[i]))
                        #TODO why isn't it with volume like in interpolating_coefficient_Navier_uv_grids_to_p_grid_volume ?
                        interpolate_u_to_p, interpolate_v_to_p = interpolating_coefficient_Navier_uv_grids_to_p_grid_height(num,gp,gu,gv,i)

                        A[range_nb_iLS_Navier,i*niu+1:(i+1)*niu] = bp * (
                            opp.HxT[iLS] * opp.iMx * opp.Hx[i] .+
                            opp.HyT[iLS] * opp.iMy * opp.Hy[i]
                        ) * sinα_p * interpolate_u_to_p
                        A[range_nb_iLS_Navier,ntu+i*niv+1:ntu+(i+1)*niv] = bp * (
                            opp.HxT[iLS] * opp.iMx * opp.Hx[i] .+
                            opp.HyT[iLS] * opp.iMy * opp.Hy[i]
                        ) * (-cosα_p) * interpolate_v_to_p
                    end
                end
                A[range_nb_iLS_Navier,range_nb_iLS_Navier] = pad(bp * (
                    opp.HxT[iLS] * opp.iMx * opp.Hx[iLS] .+
                    opp.HyT[iLS] * opp.iMy * opp.Hy[iLS]
                ) .+ opp.χ[iLS])

                # coefficients (u,v) frame to (normal, tangent)
                sin_alpha_border = Diagonal(zeros(nbp))
                sin_alpha_border.diag[1:gp.ny] .= sin.(gp.LS[iLS].α[:,1])
                sin_alpha_border.diag[gp.ny+1:gp.ny+gp.nx] .= sin.(gp.LS[iLS].α[1,:])
                sin_alpha_border.diag[gp.ny+gp.nx+1:2gp.ny+gp.nx] .= sin.(gp.LS[iLS].α[:,end])
                sin_alpha_border.diag[2gp.ny+gp.nx+1:end] .= sin.(gp.LS[iLS].α[end,:])
                cos_alpha_border = Diagonal(zeros(nbp))
                cos_alpha_border.diag[1:gp.ny] .= cos.(gp.LS[iLS].α[:,1])
                cos_alpha_border.diag[gp.ny+1:gp.ny+gp.nx] .= cos.(gp.LS[iLS].α[1,:])
                cos_alpha_border.diag[gp.ny+gp.nx+1:2gp.ny+gp.nx] .= cos.(gp.LS[iLS].α[:,end])
                cos_alpha_border.diag[2gp.ny+gp.nx+1:end] .= cos.(gp.LS[iLS].α[end,:])

                #region interpolation coefficients
                interpolate_u_to_p = spdiagm(nbp, nbu, 0 => zeros(nbp), 1 => zeros(nbp-1))
                for ii in 1:gp.ny
                    interpolate_u_to_p[ii,ii] = 1.0
                    interpolate_u_to_p[gp.ny+gp.nx+ii,gu.ny+gu.nx+ii] = 1.0
                end
                for ii in 1:gp.nx
                    interpolate_u_to_p[ii+gp.ny,ii+gu.ny] = 0.5
                    interpolate_u_to_p[ii+gp.ny,ii+gu.ny+1] = 0.5
                    interpolate_u_to_p[ii+2gp.ny+gp.nx,ii+2gu.ny+gu.nx] = 0.5
                    interpolate_u_to_p[ii+2gp.ny+gp.nx,ii+2gu.ny+gu.nx+1] = 0.5
                end
                interpolate_v_to_p = spdiagm(nbp, nbv, 0 => zeros(nbp), 1 => zeros(nbp-1))
                for ii in 1:gp.ny
                    interpolate_v_to_p[ii,ii] = 0.5
                    interpolate_v_to_p[ii,ii+1] = 0.5
                    interpolate_v_to_p[gp.ny+gp.nx+ii,gv.ny+gv.nx+ii] = 0.5
                    interpolate_v_to_p[gp.ny+gp.nx+ii,gv.ny+gv.nx+ii+1] = 0.5
                end
                for ii in 1:gp.nx
                    interpolate_v_to_p[ii+gp.ny,ii+gv.ny] = 1.0
                    interpolate_v_to_p[ii+2gp.ny+gp.nx,ii+2gv.ny+gv.nx] = 1.0
                end
                #endregion interpolation coefficients

                #region Navier slip 
                A[range_nb_iLS_Navier,border_u_velocity] = bp * (
                    opp.HxT[iLS] * opp.iMx_b * opp.Hx_b .+ 
                    opp.HyT[iLS] * opp.iMy_b * opp.Hy_b
                ) * sin_alpha_border * interpolate_u_to_p
                A[range_nb_iLS_Navier,border_v_velocity] = bp * (
                    opp.HxT[iLS] * opp.iMx_b * opp.Hx_b .+ 
                    opp.HyT[iLS] * opp.iMy_b * opp.Hy_b
                ) * (-cos_alpha_border) * interpolate_v_to_p
                #endregion Navier slip
                
                # Boundary conditions for outer boundaries
                #region interpolation coefficients
                interpolate_u_to_p = spdiagm(nbu, nbp, 0 => zeros(nbp), 1 => zeros(nbp-1))
                for ii in 1:gu.ny
                    interpolate_u_to_p[ii,ii] = 1.0
                    interpolate_u_to_p[gu.ny+gu.nx+ii,gp.ny+gp.nx+ii] = 1.0
                end
                interpolate_u_to_p[gu.ny+1,gp.ny+1]  = 1.0
                interpolate_u_to_p[gu.ny+gu.nx,gp.ny+gp.nx]  = 1.0
                interpolate_u_to_p[2gu.ny+gu.nx+1,2gp.ny+gp.nx+1]  = 1.0
                interpolate_u_to_p[end,end]  = 1.0
                for ii in 2:(gu.nx-1)
                    interpolate_u_to_p[ii+gu.ny,ii+gp.ny-1] = 0.5
                    interpolate_u_to_p[ii+gu.ny,ii+gp.ny] = 0.5
                    interpolate_u_to_p[ii+2gu.ny+gu.nx,ii+2gp.ny+gp.nx-1] = 0.5
                    interpolate_u_to_p[ii+2gu.ny+gu.nx,ii+2gp.ny+gp.nx] = 0.5
                end
                interpolate_v_to_p = spdiagm(nbv, nbp, 0 => zeros(nbp), 1 => zeros(nbp-1))
                interpolate_v_to_p[1,1]  = 1.0
                interpolate_v_to_p[gv.ny,gp.ny]  = 1.0
                interpolate_v_to_p[gv.ny+gv.nx+1,gp.ny+gp.nx+1]  = 1.0
                interpolate_v_to_p[2gv.ny+gv.nx,2gp.ny+gp.nx]  = 1.0
                for ii in 2:(gv.ny-1)
                    interpolate_v_to_p[ii,ii-1] = 0.5
                    interpolate_v_to_p[ii,ii] = 0.5
                    interpolate_v_to_p[gv.ny+gv.nx+ii,gp.ny+gp.nx+ii] = 0.5
                    interpolate_v_to_p[gv.ny+gv.nx+ii,gp.ny+gp.nx+ii] = 0.5
                end
                for ii in 1:gv.nx
                    interpolate_v_to_p[ii+gv.ny,ii+gp.ny] = 1.0
                    interpolate_v_to_p[ii+2gv.ny+gv.nx,ii+2gp.ny+gp.nx] = 1.0
                end
                #endregion interpolation coefficients

                A[border_u_velocity,range_nb_iLS_Navier] = b_bu * (
                    opu.HxT_b * opu.iMx_b' * opu.Hx[iLS] .+ opu.HyT_b * opu.iMy_b' * opu.Hy[iLS]
                ) * interpolate_x * sinα_p
                A[border_v_velocity,range_nb_iLS_Navier] = b_bv * (
                    opv.HxT_b * opv.iMx_b' * opv.Hx[iLS] .+ opv.HyT_b * opv.iMy_b' * opv.Hy[iLS]
                ) * interpolate_y * (-cosα_p)
                
                #endregion Navier

            end # if !is_navier_cl(bc_interface[iLS]) && !is_navier(bc_interface[iLS])
        end #advection

        if !is_navier_cl(bc_interface[iLS]) && !is_navier(bc_interface[iLS])
            @inbounds rhs[interfacial_nb_1_u_velocity] .= opu.χ[iLS] * vec(a0u)
            @inbounds rhs[interfacial_nb_1_v_velocity] .= opv.χ[iLS] * vec(a0v)
        else
            @inbounds rhs[range_nb_iLS_Navier] .= opp.χ[iLS] * vec(a0p)
            index_Navier += 1
        end #!is_navier_cl(bc_interface[iLS]) && !is_navier(bc_interface[iLS])

    end #for iLS in 1:num.nLS

    
    @inbounds rhs[border_u_velocity] .= opu.χ_b * vec(a0_bu) #BC for u component on borders
    @inbounds rhs[border_v_velocity] .= opv.χ_b * vec(a0_bv) #BC for v component on borders

    if num.pressure_velocity_coupling == 2
        @inbounds rhs[border_pressure] .= opp.χ_b * vec(a0_bp) #BC for u component on borders
    end


    #region set first cells to boundary velocity
    if num.pressure_velocity_coupling == 3
        set_first_cells!(A,rhs,gu,ntu-nbu,0,true,false,true,false)
        set_first_cells!(A,rhs,gv,ntu+ntv-nbv,ntu,false,true,false,true)
    end
    #endregion set first cells to boundary velocity



    #rhs div u : 0 for Incompressible

    return rhs
end #end FE_set_momentum_coupled2


"""
set first cells to boundary velocity: the first cells have no longer a velocity and are set to a BC velocity
"""
function set_first_cells!(A,rhs,grid,start_index,start_index2,left,bottom,right,top)
    # print("\n set_first_cells! ",left,bottom,right,top)
    if left
        # left
        for j in 1:grid.ny
            II = CartesianIndex(j,1)
            pII = lexicographic(II, grid.ny)
            # scal = vecb_L(scalD, grid)[j]
            u_index = start_index2 + pII
            u_index_border = start_index + j
            # print("\n u_index ",u_index," ",u_index_border)
            A[u_index,:] .= 0.0
            A[u_index,u_index] = 1.0
            A[u_index,u_index_border] = -1.0
            rhs[u_index] = 0.0
            # print("\n u_index ",u_index," ",u_index_border,A[u_index,:])

        end
    end

    if bottom
        # bottom
        for i in 1:grid.nx
            II = CartesianIndex(1,i)
            # scal = vecb_B(scalD, grid)[i]
            pII = lexicographic(II, grid.ny)
            u_index = start_index2 + pII
            u_index_border = start_index + grid.ny + i
            # print("\n u_index ",u_index," ",u_index_border)
            A[u_index,:] .= 0.0
            A[u_index,u_index] =1.0
            A[u_index,u_index_border] = -1.0
            rhs[u_index] = 0.0
            # print("\n u_index ",u_index," ",u_index_border,A[u_index,:])

        end
    end

    if right
        # right
        for j in 1:grid.ny
            II = CartesianIndex(j,grid.nx)
            # scal = vecb_R(scalD, grid)[j]
            pII = lexicographic(II, grid.ny)

            u_index = start_index2 + pII
            u_index_border = start_index + grid.ny + grid.nx+ j
            A[u_index,:] .= 0.0
            A[u_index,u_index] =1.0
            A[u_index,u_index_border] = -1.0
            rhs[u_index] = 0.0
            # print("\n u_index ",u_index," ",u_index_border,A[u_index,:])

        end
    end

    if top
        # top
        for i in 1:grid.nx
            II = CartesianIndex(grid.ny,i)
            # scal = vecb_T(scalD, grid)[i]
            pII = lexicographic(II, grid.ny)

            u_index = start_index2 + pII
            u_index_border = start_index + 2*grid.ny + grid.nx + i
            A[u_index,:] .= 0.0
            A[u_index,u_index] =1.0
            A[u_index,u_index_border] = -1.0
            rhs[u_index] = 0.0
            # print("\n u_index ",u_index," ",u_index_border,A[u_index,:])

        end
    end

end


"""
set first cells to boundary velocity: the first cells have no longer a velocity and are set to a BC velocity
"""
function set_first_cells_Neumann!(A,rhs,grid,start_index,start_index2,left,bottom,right,top)
    # print("\n set_first_cells! ",left,bottom,right,top)

    # start_index = 1
    # end_index = 1
    # start_index = 0
    # end_index = 0

    # y_range = 1:grid.ny
    y_range = 2:grid.ny-1
    # x_range = 1:grid.nx
    x_range = 2:grid.nx-1

    if left
        # left
        for j in y_range
            II = CartesianIndex(j,1)
            pII = lexicographic(II, grid.ny)
            # scal = vecb_L(scalD, grid)[j]
            u_index = start_index2 + pII
            # u_index_border = start_index + j
            u_index_ngb = u_index + grid.ny
            # print("\n u_index ",u_index," ",u_index_border)
            A[u_index,:] .= 0.0
            A[u_index,u_index] = 1.0
            A[u_index,u_index_ngb] = -1.0
            rhs[u_index] = 0.0
            # print("\n u_index ",u_index," ",u_index_border,A[u_index,:])

        end
    end

    if bottom
        # bottom
        for i in x_range
            II = CartesianIndex(1,i)
            # scal = vecb_B(scalD, grid)[i]
            pII = lexicographic(II, grid.ny)
            u_index = start_index2 + pII
            # u_index_border = start_index + grid.ny + i
            # print("\n u_index ",u_index," ",u_index_border)
            u_index_ngb = u_index + 1
            A[u_index,:] .= 0.0
            A[u_index,u_index] =1.0
            A[u_index,u_index_ngb] = -1.0
            rhs[u_index] = 0.0
            # print("\n u_index ",u_index," ",u_index_border,A[u_index,:])

        end
    end

    if right
        # right
        for j in y_range
            II = CartesianIndex(j,grid.nx)
            # scal = vecb_R(scalD, grid)[j]
            pII = lexicographic(II, grid.ny)

            u_index = start_index2 + pII
            # u_index_border = start_index + grid.ny + grid.nx+ j
            u_index_ngb = u_index - grid.ny
            A[u_index,:] .= 0.0
            A[u_index,u_index] =1.0
            A[u_index,u_index_ngb] = -1.0
            rhs[u_index] = 0.0
            # print("\n u_index ",u_index," ",u_index_border,A[u_index,:])

        end
    end

    if top
        # top
        for i in x_range
            II = CartesianIndex(grid.ny,i)
            # scal = vecb_T(scalD, grid)[i]
            pII = lexicographic(II, grid.ny)

            u_index = start_index2 + pII
            # u_index_border = start_index + 2*grid.ny + grid.nx + i
            u_index_ngb = u_index - 1

            A[u_index,:] .= 0.0
            A[u_index,u_index] =1.0
            A[u_index,u_index_ngb] = -1.0
            rhs[u_index] = 0.0
            # print("\n u_index ",u_index," ",u_index_border,A[u_index,:])

        end
    end

end


"""

"""
function deactivate_merge_first_cells_capacities!(num,grid::Mesh{GridFCx, T, N}) where {T,N}

    for iLS = 1:num.nLS #or nLS+1
        for capa_index in [2,4,5,7] #bottom 2 and top 4 # 5 cell-centered volume
            grid.LS[iLS].geoL.dcap[:,2,capa_index] .+= grid.LS[iLS].geoL.dcap[:,1,capa_index]
            grid.LS[iLS].geoL.dcap[:,1,capa_index] .= 0.0

            grid.LS[iLS].geoL.dcap[:,end-1,capa_index] .+= grid.LS[iLS].geoL.dcap[:,end,capa_index]
            grid.LS[iLS].geoL.dcap[:,end,capa_index] .= 0.0

        end
    end

end


"""

"""
function deactivate_merge_first_cells_capacities!(num,grid::Mesh{GridFCy, T, N}) where {T,N}

    for iLS = 1:num.nLS #or nLS+1
        #
        for capa_index in [1,3,5,6] #left 1 and right 3 # 5 cell-centered volume
            grid.LS[iLS].geoL.dcap[2,:,capa_index] .+= grid.LS[iLS].geoL.dcap[1,:,capa_index]
            grid.LS[iLS].geoL.dcap[1,:,capa_index] .= 0.0

            grid.LS[iLS].geoL.dcap[end-1,:,capa_index] .+= grid.LS[iLS].geoL.dcap[end,:,capa_index]
            grid.LS[iLS].geoL.dcap[end,:,capa_index] .= 0.0
        end
    end

end

