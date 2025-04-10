function run_forward(
    num, grid, grid_u, grid_v, op, phS, phL, fwd, fwdS, fwdL;
    periodic_x = false,
    periodic_y = false,
    BC_TS = Boundaries(),
    BC_TL = Boundaries(),
    BC_pS = Boundaries(),
    BC_pL = Boundaries(),
    BC_uS = Boundaries(),
    BC_uL = Boundaries(),
    BC_vS = Boundaries(),
    BC_vL = Boundaries(),
    BC_u = Boundaries(),
    BC_int = [WallNoSlip()],
    time_scheme = CN,
    ls_scheme = weno5,
    auto_reinit = false,
    heat = false,
    heat_convection = false,
    heat_liquid_phase = false,
    heat_solid_phase = false,
    navier_stokes = false,
    ns_advection = false,
    ns_liquid_phase = false,
    ns_solid_phase = false,
    hill = false,
    Vmean = false,
    levelset = true,
    speed = 0,
    analytical = false,
    verbose = false,
    show_every = 100,
    save_length = false,
    save_radius = false,
    adaptative_t = false,
    breakup = false,
    Ra = 0.0,
    St = 1,
    cutoff_length = 0.8
    )
    @unpack L0, A, N, θd, ϵ_κ, ϵ_V, σ, T_inf, L0, NB, Δ, CFL, Re, max_iterations,
            current_i, save_every, reinit_every, nb_reinit, δreinit, ϵ, m, θ₀, aniso, nLS, _nLS, nNavier = num
    @unpack opS, opL, opC_TS, opC_TL, opC_pS, opC_pL, opC_uS, opC_uL, opC_vS, opC_vL = op
    @unpack x, y, nx, ny, dx, dy, ind, LS, V = grid

    if length(BC_int) != nLS
        @error ("You have to specify $(nLS) boundary conditions.")
        return nothing
    end

    free_surface = false
    stefan = false
    flapping = false
    navier = false
    if any(is_fs, BC_int)
        free_surface = true
    end
    if any(is_stefan, BC_int)
        stefan = true
    end
    if any(is_flapping, BC_int)
        flapping = true
        D = 0.0
        xc = zeros(max_iterations+1)
        yc = zeros(max_iterations+1)
    end
    if any(is_navier_cl, BC_int) || any(is_navier, BC_int)
        navier = true
    end

    if nNavier > 1
        @warn ("When using more than 1 Navier BC, the interfaces shouldn't cross")
    end

    if free_surface && stefan && flapping
        @error ("Cannot advect the levelset using both free-surface and stefan condition.")
        return nothing
    elseif free_surface || stefan || flapping
        advection = true
    else
        advection = false
    end

    # advection = false

    iRe = 1.0 / Re
    CFL_sc = num.τ / Δ^2

    local NB_indices;

    local Cum1S = fzeros(grid_u)
    local Cum1L = fzeros(grid_u)
    local Cvm1S = fzeros(grid_v)
    local Cvm1L = fzeros(grid_v)

    local Mm1_S
    local Mum1_S
    local Mvm1_S

    local Mm1_L
    local Mum1_L
    local Mvm1_L

    θ_out = zeros(grid, 4)
    utmp = copy(LS[1].u)
    rhs_LS = fzeros(grid)

    if levelset
        NB_indices = update_all_ls_data(num, grid, grid_u, grid_v, BC_int, periodic_x, periodic_y)

        if save_radius
            n_snaps = iszero(max_iterations%save_every) ? max_iterations÷save_every+1 : max_iterations÷save_every+2
            local radius = zeros(n_snaps)
            radius[1] = find_radius(grid, LS[1])
        end
        if hill
            local radius = zeros(max_iterations+1)
            a = zeros(length(LS[1].MIXED))
            for i in eachindex(LS[1].MIXED)
                a[i] = LS[1].geoL.projection[LS[1].MIXED[i]].pos.y
            end
            radius[1] = mean(a)
        end
    elseif !levelset
        LS[1].MIXED = [CartesianIndex(-1,-1)]
        grid_u.LS[1].MIXED = [CartesianIndex(-1,-1)]
        grid_v.LS[1].MIXED = [CartesianIndex(-1,-1)]
    end

    if save_length
        fwd.length[1] = arc_length2(LS[1].geoS.projection, LS[1].MIXED)
    end

    kill_dead_cells!(phS.T, grid, LS[end].geoS)
    kill_dead_cells!(phL.T, grid, LS[end].geoL)
    
    for iLS in 1:_nLS
        @views fwd.u[iLS,1,:,:] .= LS[iLS].u
        @views fwd.ux[iLS,1,:,:] .= grid_u.LS[iLS].u
        @views fwd.uy[iLS,1,:,:] .= grid_v.LS[iLS].u
        @views fwd.κ[iLS,1,:,:] .= LS[iLS].κ
    end
    @views fwd.T[1,:,:] .= phL.T.*LS[end].geoL.cap[:,:,5] .+ phS.T[:,:].*LS[end].geoS.cap[:,:,5]
    @views fwdL.T[1,:,:] .= phL.T
    @views fwdS.T[1,:,:] .= phS.T
    @views fwdS.p[1,:,:] .= phS.p
    @views fwdL.p[1,:,:] .= phL.p
    @views fwdS.u[1,:,:] .= phS.u
    @views fwdS.v[1,:,:] .= phS.v
    @views fwdL.u[1,:,:] .= phL.u
    @views fwdL.v[1,:,:] .= phL.v

    vec1(phS.TD,grid) .= vec(phS.T)
    vec2(phS.TD,grid) .= θd
    init_borders!(phS.TD, grid, BC_TS, θd)
    @views fwdS.TD[1,:] .= phS.TD

    vec1(phL.TD,grid) .= vec(phL.T)
    vec2(phL.TD,grid) .= θd
    init_borders!(phL.TD, grid, BC_TL, θd)
    @views fwdL.TD[1,:] .= phL.TD

    vec1(phS.uD,grid_u) .= vec(phS.u)
    vec2(phS.uD,grid_u) .= num.uD
    vecb(phS.uD,grid_u) .= num.u_inf
    vec1(phL.uD,grid_u) .= vec(phL.u)
    vec2(phL.uD,grid_u) .= num.uD
    vecb(phL.uD,grid_u) .= num.u_inf
    vec1(phS.ucorrD,grid_u) .= vec(phS.u)
    vec2(phS.ucorrD,grid_u) .= num.uD
    vecb(phS.ucorrD,grid_u) .= num.u_inf
    vec1(phL.ucorrD,grid_u) .= vec(phL.u)
    vec2(phL.ucorrD,grid_u) .= num.uD
    vecb(phL.ucorrD,grid_u) .= num.u_inf
    @views fwdS.ucorrD[1,:,:] .= phS.ucorrD
    @views fwdL.ucorrD[1,:,:] .= phL.ucorrD

    vec1(phS.vD,grid_v) .= vec(phS.v)
    vec2(phS.vD,grid_v) .= num.vD
    vecb(phS.vD,grid_v) .= num.v_inf
    vec1(phL.vD,grid_v) .= vec(phL.v)
    vec2(phL.vD,grid_v) .= num.vD
    vecb(phL.vD,grid_v) .= num.v_inf
    vec1(phS.vcorrD,grid_v) .= vec(phS.v)
    vec2(phS.vcorrD,grid_v) .= num.vD
    vecb(phS.vcorrD,grid_v) .= num.v_inf
    vec1(phL.vcorrD,grid_v) .= vec(phL.v)
    vec2(phL.vcorrD,grid_v) .= num.vD
    vecb(phL.vcorrD,grid_v) .= num.v_inf
    @views fwdS.vcorrD[1,:,:] .= phS.vcorrD
    @views fwdL.vcorrD[1,:,:] .= phL.vcorrD

    @views fwdL.pD[1,:] .= phL.pD
    @views fwdS.pD[1,:] .= phS.pD

    if is_FE(time_scheme) || is_CN(time_scheme)
        NB_indices = update_all_ls_data(num, grid, grid_u, grid_v, BC_int, periodic_x, periodic_y, false)

        if navier_stokes || heat
            geoS = [LS[iLS].geoS for iLS in 1:_nLS]
            geo_uS = [grid_u.LS[iLS].geoS for iLS in 1:_nLS]
            geo_vS = [grid_v.LS[iLS].geoS for iLS in 1:_nLS]
            Lpm1_S, bc_Lpm1_S, bc_Lpm1_b_S, Lum1_S, bc_Lum1_S, bc_Lum1_b_S, Lvm1_S, bc_Lvm1_S, bc_Lvm1_b_S = set_matrices!(
                num, grid, geoS, grid_u, geo_uS, grid_v, geo_vS,
                opC_pS, opC_uS, opC_vS, periodic_x, periodic_y
            )

            geoL = [LS[iLS].geoL for iLS in 1:_nLS]
            geo_uL = [grid_u.LS[iLS].geoL for iLS in 1:_nLS]
            geo_vL = [grid_v.LS[iLS].geoL for iLS in 1:_nLS]
            Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, Lum1_L, bc_Lum1_L, bc_Lum1_b_L, Lvm1_L, bc_Lvm1_L, bc_Lvm1_b_L = set_matrices!(
                num, grid, geoL, grid_u, geo_uL, grid_v, geo_vL,
                opC_pL, opC_uL, opC_vL, periodic_x, periodic_y
            )
        end

        Mm1_L = copy(opC_pL.M)
        Mm1_S = copy(opC_pS.M)
        Mum1_L = copy(opC_uL.M)
        Mum1_S = copy(opC_uS.M)
        Mvm1_L = copy(opC_vL.M)
        Mvm1_S = copy(opC_vS.M)

        if navier_stokes || heat
            ni = grid.nx * grid.ny
            nb = 2 * grid.nx + 2 * grid.ny
            nt = (num.nLS + 1) * ni + nb
            AϕS = spzeros(nt, nt)
            AϕL = spzeros(nt, nt)

            ATS = spzeros(nt, nt)
            ATL = spzeros(nt, nt)
            BTS = spzeros(nt, nt)
            BTL = spzeros(nt, nt)

            ni = grid_u.nx * grid_u.ny
            nb = 2 * grid_u.nx + 2 * grid_u.ny
            nt = (num.nLS + 1) * ni + nb
            AuS = spzeros(nt, nt)
            AuL = spzeros(nt, nt)
            BuS = spzeros(nt, nt)
            BuL = spzeros(nt, nt)

            ni = grid_v.nx * grid_v.ny
            nb = 2 * grid_v.nx + 2 * grid_v.ny
            nt = (num.nLS + 1) * ni + nb
            AvS = spzeros(nt, nt)
            AvL = spzeros(nt, nt)
            BvS = spzeros(nt, nt)
            BvL = spzeros(nt, nt)

            ni = grid_u.nx * grid_u.ny + grid_v.nx * grid_v.ny
            nb = 2 * grid_u.nx + 2 * grid_u.ny + 2 * grid_v.nx + 2 * grid_v.ny
            nt = (nLS - nNavier + 1) * ni + nNavier * grid.nx * grid.ny + nb
            AuvS = spzeros(nt, nt)
            AuvL = spzeros(nt, nt)
            BuvS = spzeros(nt, nt)
            BuvL = spzeros(nt, nt)

            if !navier
                _ = FE_set_momentum(
                    BC_int, num, grid_u, opC_uS,
                    AuS, BuS,
                    iRe.*Lum1_S, iRe.*bc_Lum1_S, iRe.*bc_Lum1_b_S, Mum1_S, BC_uS,
                    true
                )
                _ = FE_set_momentum(
                    BC_int, num, grid_v, opC_vS,
                    AvS, BvS,
                    iRe.*Lvm1_S, iRe.*bc_Lvm1_S, iRe.*bc_Lvm1_b_S, Mvm1_S, BC_vS,
                    true
                )
            else
                _ = FE_set_momentum_coupled(
                    BC_int, num, grid, grid_u, grid_v,
                    opC_pS, opC_uS, opC_vS,
                    AuvS, BuvS,
                    iRe.*Lum1_S, iRe.*bc_Lum1_S, iRe.*bc_Lum1_b_S, Mum1_S, BC_uS,
                    iRe.*Lvm1_S, iRe.*bc_Lvm1_S, iRe.*bc_Lvm1_b_S, Mvm1_S, BC_vS,
                    true
                )
            end
            a0_p = []
            for i in 1:num.nLS
                push!(a0_p, zeros(grid))
            end
            _ = set_poisson(
                BC_int, num, grid, a0_p, opC_pS, opC_uS, opC_vS,
                AϕS, Lpm1_S, bc_Lpm1_S, bc_Lpm1_b_S, BC_pS,
                true
            )

            set_heat!(
                BC_int[1], num, grid, opC_TS, LS[1].geoS, phS, θd, BC_TS, LS[1].MIXED, LS[1].geoS.projection,
                ATS, BTS,
                opS, grid_u, grid_u.LS[1].geoS, grid_v, grid_v.LS[1].geoS,
                periodic_x, periodic_y, heat_convection, true, BC_int
            )

            if !navier
                _ = FE_set_momentum(
                    BC_int, num, grid_u, opC_uL,
                    AuL, BuL,
                    iRe.*Lum1_L, iRe.*bc_Lum1_L, iRe.*bc_Lum1_b_L, Mum1_L, BC_uL,
                    true
                )
                _ = FE_set_momentum(
                    BC_int, num, grid_v, opC_vL,
                    AvL, BvL,
                    iRe.*Lvm1_L, iRe.*bc_Lvm1_L, iRe.*bc_Lvm1_b_L, Mvm1_L, BC_vL,
                    true
                )
            else
                _ = FE_set_momentum_coupled(
                    BC_int, num, grid, grid_u, grid_v,
                    opC_pL, opC_uL, opC_vL,
                    AuvL, BuvL,
                    iRe.*Lum1_L, iRe.*bc_Lum1_L, iRe.*bc_Lum1_b_L, Mum1_L, BC_uL,
                    iRe.*Lvm1_L, iRe.*bc_Lvm1_L, iRe.*bc_Lvm1_b_L, Mvm1_L, BC_vL,
                    true
                )
            end
            a0_p = []
            for i in 1:num.nLS
                push!(a0_p, zeros(grid))
            end
            _ = set_poisson(
                BC_int, num, grid, a0_p, opC_pL, opC_uL, opC_vL,
                AϕL, Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, BC_pL,
                true
            )

            set_heat!(
                BC_int[1], num, grid, opC_TL, LS[1].geoL, phL, θd, BC_TL, LS[1].MIXED, LS[1].geoL.projection,
                ATL, BTL,
                opL, grid_u, grid_u.LS[1].geoL, grid_v, grid_v.LS[1].geoL,
                periodic_x, periodic_y, heat_convection, true, BC_int
            )
        end
    else
        error("Unknown time scheme. Available options are ForwardEuler and CrankNicolson")
    end

    if heat_convection
        NB_indices = update_all_ls_data(num, grid, grid_u, grid_v, BC_int, periodic_x, periodic_y)
    end

    V0S = volume(LS[end].geoS)
    V0L = volume(LS[end].geoL)

    # Force in x
    if adaptative_t && flapping
        for iLS in 1:nLS
            if is_flapping(BC_int[iLS])
                num.τ = min(CFL*Δ^2*Re, CFL*Δ/max(
                    abs.(V)..., sqrt(max(abs.(grid_u.V)...)^2 + max(abs.(grid_v.V)...)^2),
                    abs.(phL.u)..., abs.(phL.v)..., abs.(phS.u)..., abs.(phS.v)..., 
                    BC_int[iLS].β * BC_int[iLS].KC),
                    1.0 / BC_int[iLS].β / BC_int[iLS].βmult
                )
                CFL_sc = num.τ / Δ^2
            end
        end
    end

    current_t = 0.

    while current_i < max_iterations + 1

        if !stefan
            V .= speed*ones(ny, nx)
        end

        if heat
            if heat_solid_phase
                kill_dead_cells!(phS.T, grid, LS[1].geoS)
                veci(phS.TD,grid,1) .= vec(phS.T)
                rhs = set_heat!(
                    BC_int[1], num, grid, opC_TS, LS[1].geoS, phS, θd, BC_TS, LS[1].MIXED, LS[1].geoS.projection,
                    ATS, BTS,
                    opS, grid_u, grid_u.LS[1].geoS, grid_v, grid_v.LS[1].geoS,
                    periodic_x, periodic_y, heat_convection, advection, BC_int
                )
                mul!(rhs, BTS, phS.TD, 1.0, 1.0)

                phS.TD .= ATS \ rhs
                phS.T .= reshape(veci(phS.TD,grid,1), grid)
            end
<<<<<<< HEAD
            try
                if solid_phase
                    TS .= reshape(gmres(AS,(BS*vec(TS) + 2.0*SCUT)), (n,n))
                    #TS .= reshape(AS\(BS*vec(TS) + 2.0*SCUT), (n,n))
                    bcs!(TS, BC_TS.left, Δ)
                    bcs!(TS, BC_TS.right, Δ)
                    bcs!(TS, BC_TS.bottom, Δ)
                    bcs!(TS, BC_TS.top, Δ)
                end
                if liquid_phase
                    TL .= reshape(gmres(AL,(BL*vec(TL) + 2.0*LCUT)), (n,n))
                    #TL .= reshape(AL\(BL*vec(TL) + 2.0*LCUT), (n,n))
                    bcs!(TL, BC_TL.left, Δ)
                    bcs!(TL, BC_TL.right, Δ)
                    bcs!(TL, BC_TL.bottom, Δ)
                    bcs!(TL, BC_TL.top, Δ)
                end
            catch
                @error ("Unphysical temperature field, iteration $current_i")
                break
=======
            if heat_liquid_phase
                kill_dead_cells!(phL.T, grid, LS[1].geoL)
                veci(phL.TD,grid,1) .= vec(phL.T)
                rhs = set_heat!(
                    BC_int[1], num, grid, opC_TL, LS[1].geoL, phL, θd, BC_TL, LS[1].MIXED, LS[1].geoL.projection,
                    ATL, BTL,
                    opL, grid_u, grid_u.LS[1].geoL, grid_v, grid_v.LS[1].geoL,
                    periodic_x, periodic_y, heat_convection, advection, BC_int
                )
                mul!(rhs, BTL, phL.TD, 1.0, 1.0)

                phL.TD .= ATL \ rhs
                phL.T .= reshape(veci(phL.TD,grid,1), grid)
>>>>>>> rayleigh_benard
            end
        end

        for iLS in 1:nLS
            if is_stefan(BC_int[iLS])
                update_stefan_velocity(num, grid, iLS, LS[iLS].u, phS.T, phL.T, periodic_x, periodic_y, St, Vmean)
            elseif is_fs(BC_int[iLS])
                update_free_surface_velocity(num, grid_u, grid_v, iLS, phL.uD, phL.vD, periodic_x, periodic_y)
            elseif is_flapping(BC_int[iLS])
                if BC_int[iLS].free
                    grid_u.V .+= num.τ .* D ./ (BC_int[iLS].ρs .* volume(LS[end].geoS))
                end
                grid_v.V .= BC_int[iLS].KC .* BC_int[iLS].β .* sin(2π .* BC_int[iLS].β .* current_t)
                xc[current_i+1] = xc[current_i] + num.τ * grid_u.V[1,1]
                yc[current_i+1] = yc[current_i] + num.τ * grid_v.V[1,1]
                grid.LS[iLS].u .= sqrt.((grid.x .- xc[current_i+1]) .^ 2 .+ ((grid.y .- yc[current_i+1]) ./ num.A) .^ 2) .- num.R .* ones(grid);
                println("xc = $(xc[current_i+1]) | vx = $(grid_u.V[1,1])")
                println("yc = $(yc[current_i+1]) | vy = $(grid_v.V[1,1])")
            end
        end

        if verbose && adaptative_t && (current_i-1)%show_every == 0 
            println("τ = $(num.τ)")
        end

        if advection
<<<<<<< HEAD
            IIOE(LSA, LSB, u, V, inside, CFL, Δ, n)
            try
                u .= reshape(gmres(LSA,(LSB*vec(u))), (n,n))
                #u .= reshape(LSA\(LSB*vec(u)), (n,n))
            catch
                @error ("Inadequate level set function, iteration $current_i")
                break
=======
            # print("advection")
            for (iLS, bc) in enumerate(BC_int)
                if is_stefan(bc)
                    IIOE_normal!(grid, LS[iLS].A, LS[iLS].B, LS[iLS].u, V, CFL_sc, periodic_x, periodic_y)
                    LS[iLS].u .= reshape(gmres(LS[iLS].A, LS[iLS].B * vec(LS[iLS].u)), grid)
                    # u .= sqrt.((x .- current_i*Δ/1).^ 2 + y .^ 2) - (0.5) * ones(nx, ny);
                elseif is_fs(bc)
                    # rhs_LS .= 0.0
                    # LS[iLS].A.nzval .= 0.0
                    # LS[iLS].B.nzval .= 0.0
                    # IIOE!(grid, grid_u, grid_v, LS[iLS].A, LS[iLS].B, θ_out, num.τ, periodic_x, periodic_y)
                    # BC_LS!(grid, LS[iLS].u, LS[iLS].A, LS[iLS].B, rhs_LS, BC_u)
                    # utmp .= reshape(gmres(LS[iLS].A, LS[iLS].B * vec(LS[iLS].u) .+ rhs_LS), grid)

                    # rhs_LS .= 0.0
                    # S2IIOE!(grid, grid_u, grid_v, LS[iLS].A, LS[iLS].B, utmp, LS[iLS].u, θ_out, num.τ, periodic_x, periodic_y)
                    # BC_LS!(grid, LS[iLS].u, LS[iLS].A, LS[iLS].B, rhs_LS, BC_u)
                    # LS[iLS].u .= reshape(gmres(LS[iLS].A, LS[iLS].B * vec(LS[iLS].u) .+ rhs_LS), grid)

                    # Project velocities to the normal and use advecion scheme for advection just
                    # in the normal direction
                    tmpVx = zeros(grid)
                    tmpVy = zeros(grid)
                    V .= 0.0
                    @inbounds @threads for II in grid.LS[iLS].MIXED
                        cap1 = grid_u.LS[iLS].geoL.cap[II,5]
                        cap3 = grid_u.LS[iLS].geoL.cap[δx⁺(II),5]
                        tmpVx[II] = (grid_u.V[II] * cap1 + grid_u.V[δx⁺(II)] * cap3) / (cap1 + cap3 + eps(0.01))

                        cap2 = grid_v.LS[iLS].geoL.cap[II,5]
                        cap4 = grid_v.LS[iLS].geoL.cap[δy⁺(II),5]
                        tmpVy[II] = (grid_v.V[II] * cap2 + grid_v.V[δy⁺(II)] * cap4) / (cap2 + cap4 + eps(0.01))

                        tmpV = sqrt(tmpVx[II]^2 + tmpVy[II]^2)
                        β = atan(tmpVy[II], tmpVx[II])
                        if grid.LS[iLS].α[II] > 0.0 && β < 0.0
                            β += 2π
                        end
                        if grid.LS[iLS].α[II] < 0.0 && β > 0.0
                            β -= 2π
                        end

                        V[II] = tmpV * cos(β - grid.LS[iLS].α[II])
                    end

                    i_ext, l_ext, b_ext, r_ext, t_ext = indices_extension(grid, grid.LS[iLS], grid.ind.inside, periodic_x, periodic_y)
                    field_extension!(grid, grid.LS[iLS].u, V, i_ext, l_ext, b_ext, r_ext, t_ext, num.NB, periodic_x, periodic_y)

                    rhs_LS .= 0.0
                    IIOE_normal!(grid, LS[iLS].A, LS[iLS].B, LS[iLS].u, V, CFL_sc, periodic_x, periodic_y)
                    BC_LS!(grid, LS[iLS].u, LS[iLS].A, LS[iLS].B, rhs_LS, BC_u)
                    BC_LS_interior!(num, grid, iLS, LS[iLS].A, LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
                    LS[iLS].u .= reshape(gmres(LS[iLS].A, LS[iLS].B * vec(LS[iLS].u) .+ rhs_LS), grid)
                end
>>>>>>> rayleigh_benard
            end
            if analytical
                u[ind.b_top[1]] .= sqrt.(x[ind.b_top[1]] .^ 2 + y[ind.b_top[1]] .^ 2) .- (num.R + speed*current_i*num.τ);
                u[ind.b_bottom[1]] .= sqrt.(x[ind.b_bottom[1]] .^ 2 + y[ind.b_bottom[1]] .^ 2) .- (num.R + speed*current_i*num.τ);
                u[ind.b_left[1]] .= sqrt.(x[ind.b_left[1]] .^ 2 + y[ind.b_left[1]] .^ 2) .- (num.R + speed*current_i*num.τ);
                u[ind.b_right[1]] .= sqrt.(x[ind.b_right[1]] .^ 2 + y[ind.b_right[1]] .^ 2) .- (num.R + speed*current_i*num.τ);
            elseif nb_reinit > 0
                if auto_reinit && (current_i-1)%num.reinit_every == 0
                    for iLS in 1:nLS
                        if !is_wall(BC_int[iLS])
                            ls_rg = rg(num, grid, LS[iLS].u, periodic_x, periodic_y, BC_int)[1]
                            println(ls_rg)
                            if ls_rg >= δreinit || current_i == 1
                                println("yes")
                                RK2_reinit!(ls_scheme, grid, ind, iLS, LS[iLS].u, nb_reinit, periodic_x, periodic_y, BC_u, BC_int)
                            end
                        end
                    end
                elseif (current_i-1)%num.reinit_every == 0
                    for iLS in 1:nLS
                        if !is_wall(BC_int[iLS])
                            RK2_reinit!(ls_scheme, grid, ind, iLS, LS[iLS].u, nb_reinit, periodic_x, periodic_y, BC_u, BC_int)
                        end
                    end
                else
                    for iLS in 1:nLS
                        if !is_wall(BC_int[iLS])
                            RK2_reinit!(ls_scheme, grid, ind, iLS, LS[iLS].u, nb_reinit, periodic_x, periodic_y, BC_u, BC_int, true)
                        end
                    end
                end
            end
            # numerical breakup
            if free_surface && breakup
                count = breakup_f(LS[1].u, nx, ny, dx, dy, periodic_x, periodic_y, NB_indices, 1e-5)
                if count > 0
                    RK2_reinit!(ls_scheme, grid, ind, 1, LS[1].u, nb_reinit, periodic_x, periodic_y, BC_u, BC_int)
                end
            end
        end

        if verbose
            if (current_i-1)%show_every == 0
                printstyled(color=:green, @sprintf "\n Current iteration : %d (%d%%) | t = %.2e \n" (current_i-1) 100*(current_i-1)/max_iterations current_t)
                if heat && length(LS[end].MIXED) != 0
                    print(@sprintf "V_mean = %.2e  V_max = %.2e  V_min = %.2e\n" mean(V[LS[1].MIXED]) findmax(V[LS[1].MIXED])[1] findmin(V[LS[1].MIXED])[1])
                    print(@sprintf "κ_mean = %.2e  κ_max = %.2e  κ_min = %.2e\n" mean(LS[1].κ[LS[1].MIXED]) findmax(LS[1].κ[LS[1].MIXED])[1] findmin(LS[1].κ[LS[1].MIXED])[1])
                elseif advection && length(LS[end].MIXED) != 0
                    V_mean = mean([mean(grid_u.V[LS[1].MIXED]), mean(grid_v.V[LS[1].MIXED])])
                    V_max = max(findmax(grid_u.V[LS[1].MIXED])[1], findmax(grid_v.V[LS[1].MIXED])[1])
                    V_min = min(findmin(grid_u.V[LS[1].MIXED])[1], findmin(grid_v.V[LS[1].MIXED])[1])
                    print(@sprintf "Vol_ratio = %.3f%%\n" (volume(LS[end].geoL) / V0L * 100))
                    print(@sprintf "V_mean = %.2e  V_max = %.2e  V_min = %.2e\n" V_mean V_max V_min)
                    print(@sprintf "κ_mean = %.2e  κ_max = %.2e  κ_min = %.2e\n" mean(LS[1].κ[LS[1].MIXED]) findmax(LS[1].κ[LS[1].MIXED])[1] findmin(LS[1].κ[LS[1].MIXED])[1])
                end
                if navier_stokes
                    if ns_solid_phase
                        normuS = norm(phS.u)
                        normvS = norm(phS.v)
                        normpS = norm(phS.p.*num.τ)
                        print("$(@sprintf("norm(uS) %.6e", normuS))\t$(@sprintf("norm(vS) %.6e", normvS))\t$(@sprintf("norm(pS) %.6e", normpS))\n")
                    end
                    if ns_liquid_phase
                        normuL = norm(phL.u)
                        normvL = norm(phL.v)
                        normpL = norm(phL.p.*num.τ)
                        print("$(@sprintf("norm(uL) %.6e", normuL))\t$(@sprintf("norm(vL) %.6e", normvL))\t$(@sprintf("norm(pL) %.6e", normpL))\n")
                    end
                end
            end
        end


        if levelset && (advection || current_i<2)
            # print("level set and advection")
            NB_indices = update_all_ls_data(num, grid, grid_u, grid_v, BC_int, periodic_x, periodic_y)

            if heat && heat_convection && heat_liquid_phase && navier_stokes && ns_liquid_phase
                fwd.RB[1,current_i] = current_t
                # h = 0.
                # for jj in grid.LS[1].MIXED
                #     h += grid.y[jj] + num.Δ * grid.LS[1].mid_point[jj].y + 0.5
                # end
                # tmp = h / length(grid.LS[1].MIXED) 
                sum = 0
                count = 0
                for i in 1:grid.nx
                    for j in 1:(grid.ny-1)
                        diff = grid.LS[1].u[j,i]*grid.LS[1].u[j+1,i]
                        if diff < 0
                            pos = grid.y[j+1,1] + grid.LS[1].u[j+1,i]
                            sum += pos 
                            count += 1
                        end
                    end
                end
                tmp = (sum / count + 0.5)
                fwd.RB[2,current_i] = tmp
                fwd.RB[3,current_i] = Ra * (1. - num.θd) * tmp^3     
                if verbose
                    if (current_i-1)%show_every == 0 
                        print("$(@sprintf("height %g", fwd.RB[2,current_i]))\t$(@sprintf("Ra_eff %g", fwd.RB[3,current_i]))\n")
                    end
                end
                if tmp > cutoff_length
                    return nothing
                end
            end

            LS[end].geoL.fresh .= false
            LS[end].geoS.fresh .= false
            grid_u.LS[end].geoL.fresh .= false
            grid_u.LS[end].geoS.fresh .= false
            grid_v.LS[end].geoL.fresh .= false
            grid_v.LS[end].geoS.fresh .= false

            get_fresh_cells!(grid, grid.LS[end].geoS, Mm1_S, grid.ind.all_indices)
            get_fresh_cells!(grid, grid.LS[end].geoL, Mm1_L, grid.ind.all_indices)
            get_fresh_cells!(grid_u, grid_u.LS[end].geoS, Mum1_S, grid_u.ind.all_indices)
            get_fresh_cells!(grid_u, grid_u.LS[end].geoL, Mum1_L, grid_u.ind.all_indices)
            get_fresh_cells!(grid_v, grid_v.LS[end].geoS, Mvm1_S, grid_v.ind.all_indices)
            get_fresh_cells!(grid_v, grid_v.LS[end].geoL, Mvm1_L, grid_v.ind.all_indices)

            FRESH_L_u = findall(grid_u.LS[end].geoL.fresh)
            FRESH_S_u = findall(grid_u.LS[end].geoS.fresh)
            FRESH_L_v = findall(grid_v.LS[end].geoL.fresh)
            FRESH_S_v = findall(grid_v.LS[end].geoS.fresh)

            if navier_stokes
                init_fresh_cells!(grid_u, veci(phS.uD,grid_u,1), veci(phS.uD,grid_u,1),
                    grid_u.LS[end].geoS.projection, FRESH_S_u, periodic_x, periodic_y)
                init_fresh_cells!(grid_v, veci(phS.vD,grid_v,1), veci(phS.vD,grid_v,1),
                    grid_v.LS[end].geoS.projection, FRESH_S_v, periodic_x, periodic_y)
                init_fresh_cells!(grid_u, veci(phS.uD,grid_u,2), veci(phS.uD,grid_u,1),
                    grid_u.LS[end].geoS.projection, FRESH_S_u, periodic_x, periodic_y)
                init_fresh_cells!(grid_v, veci(phS.vD,grid_v,2), veci(phS.vD,grid_v,1),
                    grid_v.LS[end].geoS.projection, FRESH_S_v, periodic_x, periodic_y)

                init_fresh_cells!(grid_u, veci(phL.uD,grid_u,1), veci(phL.uD,grid_u,1),
                    grid_u.LS[end].geoL.projection, FRESH_L_u, periodic_x, periodic_y)
                init_fresh_cells!(grid_v, veci(phL.vD,grid_v,1), veci(phL.vD,grid_v,1),
                    grid_v.LS[end].geoL.projection, FRESH_L_v, periodic_x, periodic_y)
                init_fresh_cells!(grid_u, veci(phL.uD,grid_u,2), veci(phL.uD,grid_u,1),
                    grid_u.LS[end].geoL.projection, FRESH_L_u, periodic_x, periodic_y)
                init_fresh_cells!(grid_v, veci(phL.vD,grid_v,2), veci(phL.vD,grid_v,1),
                    grid_v.LS[end].geoL.projection, FRESH_L_v, periodic_x, periodic_y)
            end

            if iszero(current_i%save_every) || current_i==max_iterations
                snap = current_i÷save_every+1
                if save_radius
                    radius[snap] = find_radius(grid, LS[1])
                end
                if hill
                    a = zeros(length(LS[1].MIXED))
                    for i in eachindex(LS[1].MIXED)
                        a[i] = LS[1].geoL.projection[LS[1].MIXED[i]].pos.y
                    end
                    radius[snap] = mean(a)
                end
                if save_length
                    fwd.length[snap] = arc_length2(LS[1].geoS.projection, LS[1].MIXED)
                end
            end
        end

        if navier_stokes
            if !advection
                no_slip_condition!(num, grid, grid_u, grid_u.LS[1], grid_v, grid_v.LS[1], periodic_x, periodic_y)
                # grid_u.V .= Δ / (1 * num.τ)
                # grid_v.V .= 0.0
            end

            if ns_solid_phase
                geoS = [LS[iLS].geoS for iLS in 1:_nLS]
                geo_uS = [grid_u.LS[iLS].geoS for iLS in 1:_nLS]
                geo_vS = [grid_v.LS[iLS].geoS for iLS in 1:_nLS]
                Lpm1_S, bc_Lpm1_S, bc_Lpm1_b_S, Lum1_S, bc_Lum1_S, bc_Lum1_b_S, Lvm1_S, bc_Lvm1_S, bc_Lvm1_b_S,Mm1_S, Mum1_S, Mvm1_S, Cum1S, Cvm1S = pressure_projection!(
                    time_scheme, BC_int,
                    num, grid, geoS, grid_u, geo_uS, grid_v, geo_vS, phS,
                    BC_uS, BC_vS, BC_pS,
                    opC_pS, opC_uS, opC_vS, opS,
                    AuS, BuS, AvS, BvS, AϕS, AuvS, BuvS,
                    Lpm1_S, bc_Lpm1_S, bc_Lpm1_b_S, Lum1_S, bc_Lum1_S, bc_Lum1_b_S, Lvm1_S, bc_Lvm1_S, bc_Lvm1_b_S,
                    Cum1S, Cvm1S, Mum1_S, Mvm1_S,
                    periodic_x, periodic_y, ns_advection, advection, current_i, Ra, navier
                )
            end
            if ns_liquid_phase
                geoL = [LS[iLS].geoL for iLS in 1:_nLS]
                geo_uL = [grid_u.LS[iLS].geoL for iLS in 1:_nLS]
                geo_vL = [grid_v.LS[iLS].geoL for iLS in 1:_nLS]
                Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, Lum1_L, bc_Lum1_L, bc_Lum1_b_L, Lvm1_L, bc_Lvm1_L, bc_Lvm1_b_L, Mm1_L, Mum1_L, Mvm1_L, Cum1L, Cvm1L = pressure_projection!(
                    time_scheme, BC_int,
                    num, grid, geoL, grid_u, geo_uL, grid_v, geo_vL, phL,
                    BC_uL, BC_vL, BC_pL,
                    opC_pL, opC_uL, opC_vL, opL,
                    AuL, BuL, AvL, BvL, AϕL, AuvL, BuvL,
                    Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, Lum1_L, bc_Lum1_L, bc_Lum1_b_L, Lvm1_L, bc_Lvm1_L, bc_Lvm1_b_L,
                    Cum1L, Cvm1L, Mum1_L, Mvm1_L,
                    periodic_x, periodic_y, ns_advection, advection, current_i, Ra, navier
                )
                # if current_i == 1
                #     phL.u .= -0.5 .* grid_u.y .+ getproperty.(grid_u.LS[1].geoL.centroid, :y) .* grid_u.dy
                #     phL.v .= 0.5 .* grid_v.x .+ getproperty.(grid_v.LS[1].geoL.centroid, :x) .* grid_v.dx
                #     phL.u[grid_u.LS[1].SOLID] .= 0.0
                #     phL.v[grid_v.LS[1].SOLID] .= 0.0
                # end
                # linear_advection!(
                #     num, grid, LS[1].geoL, grid_u, grid_u.LS[1].geoL, grid_v, grid_v.LS[1].geoL, phL,
                #     BC_uL, BC_vL, opL
                # )
            end
        end

        cD, cL, D, L = force_coefficients!(num, grid, grid_u, grid_v, opL, fwd, phL; step = current_i+1, saveCoeffs = false)

        current_t += num.τ
        if iszero(current_i%save_every) || current_i==max_iterations
            snap = current_i÷save_every+1
            if current_i==max_iterations
                snap = size(fwd.T,1)
            end
            fwd.t[snap] = current_t
            @views fwd.V[snap,:,:] .= V
            for iLS in 1:_nLS
                @views fwd.u[iLS,snap,:,:] .= LS[iLS].u
                @views fwd.ux[iLS,snap,:,:] .= grid_u.LS[iLS].u
                @views fwd.uy[iLS,snap,:,:] .= grid_v.LS[iLS].u
                @views fwd.κ[iLS,snap,:,:] .= LS[iLS].κ
            end

            if heat_solid_phase && heat_liquid_phase
                @views fwd.T[snap,:,:] .= phL.T.*LS[end].geoL.cap[:,:,5] .+ phS.T.*LS[end].geoS.cap[:,:,5]
            end
            if heat_solid_phase
                @views fwd.T[snap,:,:] .= phL.T + phS.T
                @views fwdS.T[snap,:,:] .= phS.T
                @views fwdS.TD[snap,:] .= phS.TD
            end
            if heat_liquid_phase
                @views fwd.T[snap,:,:] .= phL.T + phS.T
                @views fwdL.T[snap,:,:] .= phL.T
                @views fwdL.TD[snap,:] .= phL.TD
            end

            if ns_solid_phase
                @views fwdS.p[snap,:,:] .= phS.p
                @views fwdS.pD[snap,:] .= phS.pD
                @views fwdS.ϕ[snap,:,:] .= phS.ϕ
                @views fwdS.u[snap,:,:] .= phS.u
                @views fwdS.v[snap,:,:] .= phS.v
                @views fwdS.ucorrD[snap,:,:] .= phS.ucorrD
                @views fwdS.vcorrD[snap,:,:] .= phS.vcorrD
            end
            if ns_liquid_phase
                @views fwdL.p[snap,:,:] .= phL.p
                @views fwdL.pD[snap,:] .= phL.pD
                @views fwdL.ϕ[snap,:,:] .= phL.ϕ
                @views fwdL.u[snap,:,:] .= phL.u
                @views fwdL.v[snap,:,:] .= phL.v
                @views fwdL.ucorrD[snap,:,:] .= phL.ucorrD
                @views fwdL.vcorrD[snap,:,:] .= phL.vcorrD
                # @views fwd.Cd[snap] = cD
                # @views fwd.Cl[snap] = cL
            end
            if advection
                fwdS.Vratio[snap] = volume(LS[end].geoS) / V0S
                fwdL.Vratio[snap] = volume(LS[end].geoL) / V0L
            end
        end
        @views fwd.Cd[current_i+1] = cD
        @views fwd.Cl[current_i+1] = cL

        if (any(isnan, phL.uD) || any(isnan, phL.vD) || any(isnan, phL.TD) || any(isnan, phS.uD) || any(isnan, phS.vD) || any(isnan, phS.TD) ||
            norm(phL.u) > 1e8 || norm(phS.u) > 1e8 || norm(phL.T) > 1e8 || norm(phS.T) > 1e8)
            println(@sprintf "\n CRASHED after %d iterations \n" current_i)
            if flapping
                return xc, yc
            else
                return nothing
            end
        end

        current_i += 1

        if adaptative_t && flapping
            for iLS in 1:nLS
                if is_flapping(BC_int[iLS])
                    num.τ = min(CFL*Δ^2*Re, CFL*Δ/max(
                        abs.(V)..., sqrt(max(abs.(grid_u.V)...)^2 + max(abs.(grid_v.V)...)^2),
                        abs.(phL.u)..., abs.(phL.v)..., abs.(phS.u)..., abs.(phS.v)..., 
                        BC_int[iLS].β * BC_int[iLS].KC),
                        1.0 / BC_int[iLS].β / BC_int[iLS].βmult
                    )
                    CFL_sc = num.τ / Δ^2
                end
            end
        elseif adaptative_t
            num.τ = min(CFL*Δ^2*Re, CFL*Δ/max(
                abs.(V)..., abs.(grid_u.V)..., abs.(grid_v.V)..., 
                abs.(phL.u)..., abs.(phL.v)..., abs.(phS.u)..., abs.(phS.v)...)
            )
            CFL_sc = num.τ / Δ^2
        end
    end

    if verbose
        try
            printstyled(color=:blue, @sprintf "\n Final iteration : %d (%d%%) | t = %.2e \n" (current_i-1) 100*(current_i-1)/max_iterations current_t)
            if stefan && advection
                print(@sprintf "V_mean = %.2e  V_max = %.2e  V_min = %.2e  V_stdev = %.5f\n" mean(V[LS[1].MIXED]) findmax(V[LS[1].MIXED])[1] findmin(V[LS[1].MIXED])[1] std(V[LS[1].MIXED]))
                print(@sprintf "κ_mean = %.2e  κ_max = %.2e  κ_min = %.2e  κ_stdev = %.5f\n" mean(LS[1].κ[LS[1].MIXED]) findmax(LS[1].κ[LS[1].MIXED])[1] findmin(LS[1].κ[LS[1].MIXED])[1] std(LS[1].κ[LS[1].MIXED]))
            end
            if free_surface && advection
                print(@sprintf "Vol_ratio = %.3f%%\n" (volume(LS[end].geoL) / V0L * 100))
                print(@sprintf "V_mean = %.2e  V_max = %.2e  V_min = %.2e  V_stdev = %.5f\n" mean(V[LS[1].MIXED]) findmax(V[LS[1].MIXED])[1] findmin(V[LS[1].MIXED])[1] std(V[LS[1].MIXED]))
                print(@sprintf "κ_mean = %.2e  κ_max = %.2e  κ_min = %.2e  κ_stdev = %.5f\n" mean(LS[1].κ[LS[1].MIXED]) findmax(LS[1].κ[LS[1].MIXED])[1] findmin(LS[1].κ[LS[1].MIXED])[1] std(LS[1].κ[LS[1].MIXED]))
            end
            if navier_stokes
                if ns_solid_phase
                    normuS = norm(phS.u)
                    normvS = norm(phS.v)
                    normpS = norm(phS.p.*num.τ)
                    print("$(@sprintf("norm(uS) %.6e", normuS))\t$(@sprintf("norm(vS) %.6e", normvS))\t$(@sprintf("norm(pS) %.6e", normpS))\n")
                end
                if ns_liquid_phase
                    normuL = norm(phL.u)
                    normvL = norm(phL.v)
                    normpL = norm(phL.p.*num.τ)
                    print("$(@sprintf("norm(uL) %.6e", normuL))\t$(@sprintf("norm(vL) %.6e", normvL))\t$(@sprintf("norm(pL) %.6e", normpL))\n")
                end
            end
            print("\n\n")
        catch
            @show (length(LS[end].MIXED))
        end
    end

    if levelset && (save_radius || hill)
        return radius
    elseif flapping
        return xc, yc
    else
        return nothing
    end
end












function run_backward2(
    num, grid, grid_u, grid_v, op, phS, phL, fwd, fwdS, fwdL;
    periodic_x = false,
    periodic_y = false,
    BC_TS = Boundaries(),
    BC_TL = Boundaries(),
    BC_pS = Boundaries(),
    BC_pL = Boundaries(),
    BC_uS = Boundaries(),
    BC_uL = Boundaries(),
    BC_vS = Boundaries(),
    BC_vL = Boundaries(),
    BC_u = Boundaries(),
    BC_int = [WallNoSlip()],
    time_scheme = CN,
    ls_scheme = weno5,
    auto_reinit = false,
    heat = false,
    heat_convection = false,
    heat_liquid_phase = false,
    heat_solid_phase = false,
    navier_stokes = false,
    ns_advection = false,
    ns_liquid_phase = false,
    ns_solid_phase = false,
    hill = false,
    Vmean = false,
    levelset = true,
    speed = 0,
    analytical = false,
    verbose = false,
    show_every = 100,
    save_length = false,
    save_radius = false,
    adaptative_t = false,
    breakup = false,
    Ra = 0.0,
    St = 1,
    cutoff_length = 0.8
    )
    @unpack L0, A, N, θd, ϵ_κ, ϵ_V, σ, T_inf, L0, NB, Δ, CFL, Re, max_iterations,
            current_i, save_every, reinit_every, nb_reinit, δreinit, ϵ, m, θ₀, aniso, nLS, _nLS, nNavier = num
    @unpack opS, opL, opC_TS, opC_TL, opC_pS, opC_pL, opC_uS, opC_uL, opC_vS, opC_vL = op
    @unpack x, y, nx, ny, dx, dy, ind, LS, V = grid

    if length(BC_int) != nLS
        @error ("You have to specify $(nLS) boundary conditions.")
        return nothing
    end

    free_surface = false
    stefan = false
    flapping = false
    navier = false
    if any(is_fs, BC_int)
        free_surface = true
    end
    if any(is_stefan, BC_int)
        stefan = true
    end
    if any(is_flapping, BC_int)
        flapping = true
        D = 0.0
        xc = zeros(max_iterations+1)
        yc = zeros(max_iterations+1)
    end
    if any(is_navier_cl, BC_int) || any(is_navier, BC_int)
        navier = true
    end

    if nNavier > 1
        @warn ("When using more than 1 Navier BC, the interfaces shouldn't cross")
    end

    if free_surface && stefan && flapping
        @error ("Cannot advect the levelset using both free-surface and stefan condition.")
        return nothing
    elseif free_surface || stefan || flapping
        advection = true
    else
        advection = false
    end

    # advection = false

    iRe = 1.0 / Re
    CFL_sc = num.τ / Δ^2

    local NB_indices;

    local Cum1S = fzeros(grid_u)
    local Cum1L = fzeros(grid_u)
    local Cvm1S = fzeros(grid_v)
    local Cvm1L = fzeros(grid_v)

    local Mm1_S
    local Mum1_S
    local Mvm1_S

    local Mm1_L
    local Mum1_L
    local Mvm1_L

    θ_out = zeros(grid, 4)
    utmp = copy(LS[1].u)
    rhs_LS = fzeros(grid)

    if levelset
        NB_indices = update_all_ls_data(num, grid, grid_u, grid_v, BC_int, periodic_x, periodic_y)

        if save_radius
            n_snaps = iszero(max_iterations%save_every) ? max_iterations÷save_every+1 : max_iterations÷save_every+2
            local radius = zeros(n_snaps)
            radius[1] = find_radius(grid, LS[1])
        end
        if hill
            local radius = zeros(max_iterations+1)
            a = zeros(length(LS[1].MIXED))
            for i in eachindex(LS[1].MIXED)
                a[i] = LS[1].geoL.projection[LS[1].MIXED[i]].pos.y
            end
            radius[1] = mean(a)
        end
    elseif !levelset
        LS[1].MIXED = [CartesianIndex(-1,-1)]
        grid_u.LS[1].MIXED = [CartesianIndex(-1,-1)]
        grid_v.LS[1].MIXED = [CartesianIndex(-1,-1)]
    end

    if save_length
        fwd.length[1] = arc_length2(LS[1].geoS.projection, LS[1].MIXED)
    end

    kill_dead_cells!(phS.T, grid, LS[end].geoS)
    kill_dead_cells!(phL.T, grid, LS[end].geoL)
    
    # for iLS in 1:_nLS
    #     @views fwd.u[iLS,1,:,:] .= LS[iLS].u
    #     @views fwd.ux[iLS,1,:,:] .= grid_u.LS[iLS].u
    #     @views fwd.uy[iLS,1,:,:] .= grid_v.LS[iLS].u
    #     @views fwd.κ[iLS,1,:,:] .= LS[iLS].κ
    # end
    # @views fwd.T[1,:,:] .= phL.T.*LS[end].geoL.cap[:,:,5] .+ phS.T[:,:].*LS[end].geoS.cap[:,:,5]
    # @views fwdL.T[1,:,:] .= phL.T
    # @views fwdS.T[1,:,:] .= phS.T
    # @views fwdS.p[1,:,:] .= phS.p
    # @views fwdL.p[1,:,:] .= phL.p
    # @views fwdS.u[1,:,:] .= phS.u
    # @views fwdS.v[1,:,:] .= phS.v
    # @views fwdL.u[1,:,:] .= phL.u
    # @views fwdL.v[1,:,:] .= phL.v

    vec1(phS.TD,grid) .= vec(phS.T)
    vec2(phS.TD,grid) .= θd
    init_borders!(phS.TD, grid, BC_TS, θd)
    @views fwdS.TD[1,:] .= phS.TD

    vec1(phL.TD,grid) .= vec(phL.T)
    vec2(phL.TD,grid) .= θd
    init_borders!(phL.TD, grid, BC_TL, θd)
    @views fwdL.TD[1,:] .= phL.TD

    vec1(phS.uD,grid_u) .= vec(phS.u)
    vec2(phS.uD,grid_u) .= num.uD
    vecb(phS.uD,grid_u) .= num.u_inf
    vec1(phL.uD,grid_u) .= vec(phL.u)
    vec2(phL.uD,grid_u) .= num.uD
    vecb(phL.uD,grid_u) .= num.u_inf
    vec1(phS.ucorrD,grid_u) .= vec(phS.u)
    vec2(phS.ucorrD,grid_u) .= num.uD
    vecb(phS.ucorrD,grid_u) .= num.u_inf
    vec1(phL.ucorrD,grid_u) .= vec(phL.u)
    vec2(phL.ucorrD,grid_u) .= num.uD
    vecb(phL.ucorrD,grid_u) .= num.u_inf
    @views fwdS.ucorrD[1,:,:] .= phS.ucorrD
    @views fwdL.ucorrD[1,:,:] .= phL.ucorrD

    vec1(phS.vD,grid_v) .= vec(phS.v)
    vec2(phS.vD,grid_v) .= num.vD
    vecb(phS.vD,grid_v) .= num.v_inf
    vec1(phL.vD,grid_v) .= vec(phL.v)
    vec2(phL.vD,grid_v) .= num.vD
    vecb(phL.vD,grid_v) .= num.v_inf
    vec1(phS.vcorrD,grid_v) .= vec(phS.v)
    vec2(phS.vcorrD,grid_v) .= num.vD
    vecb(phS.vcorrD,grid_v) .= num.v_inf
    vec1(phL.vcorrD,grid_v) .= vec(phL.v)
    vec2(phL.vcorrD,grid_v) .= num.vD
    vecb(phL.vcorrD,grid_v) .= num.v_inf
    @views fwdS.vcorrD[1,:,:] .= phS.vcorrD
    @views fwdL.vcorrD[1,:,:] .= phL.vcorrD

    @views fwdL.pD[1,:] .= phL.pD
    @views fwdS.pD[1,:] .= phS.pD

    if is_FE(time_scheme) || is_CN(time_scheme)
        NB_indices = update_all_ls_data(num, grid, grid_u, grid_v, BC_int, periodic_x, periodic_y, false)

        if navier_stokes || heat
            geoS = [LS[iLS].geoS for iLS in 1:_nLS]
            geo_uS = [grid_u.LS[iLS].geoS for iLS in 1:_nLS]
            geo_vS = [grid_v.LS[iLS].geoS for iLS in 1:_nLS]
            Lpm1_S, bc_Lpm1_S, bc_Lpm1_b_S, Lum1_S, bc_Lum1_S, bc_Lum1_b_S, Lvm1_S, bc_Lvm1_S, bc_Lvm1_b_S = set_matrices!(
                num, grid, geoS, grid_u, geo_uS, grid_v, geo_vS,
                opC_pS, opC_uS, opC_vS, periodic_x, periodic_y
            )

            geoL = [LS[iLS].geoL for iLS in 1:_nLS]
            geo_uL = [grid_u.LS[iLS].geoL for iLS in 1:_nLS]
            geo_vL = [grid_v.LS[iLS].geoL for iLS in 1:_nLS]
            Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, Lum1_L, bc_Lum1_L, bc_Lum1_b_L, Lvm1_L, bc_Lvm1_L, bc_Lvm1_b_L = set_matrices!(
                num, grid, geoL, grid_u, geo_uL, grid_v, geo_vL,
                opC_pL, opC_uL, opC_vL, periodic_x, periodic_y
            )
        end

        Mm1_L = copy(opC_pL.M)
        Mm1_S = copy(opC_pS.M)
        Mum1_L = copy(opC_uL.M)
        Mum1_S = copy(opC_uS.M)
        Mvm1_L = copy(opC_vL.M)
        Mvm1_S = copy(opC_vS.M)

        if navier_stokes || heat
            ni = grid.nx * grid.ny
            nb = 2 * grid.nx + 2 * grid.ny
            nt = (num.nLS + 1) * ni + nb
            AϕS = spzeros(nt, nt)
            AϕL = spzeros(nt, nt)

            ATS = spzeros(nt, nt)
            ATL = spzeros(nt, nt)
            BTS = spzeros(nt, nt)
            BTL = spzeros(nt, nt)

            ni = grid_u.nx * grid_u.ny
            nb = 2 * grid_u.nx + 2 * grid_u.ny
            nt = (num.nLS + 1) * ni + nb
            AuS = spzeros(nt, nt)
            AuL = spzeros(nt, nt)
            BuS = spzeros(nt, nt)
            BuL = spzeros(nt, nt)

            ni = grid_v.nx * grid_v.ny
            nb = 2 * grid_v.nx + 2 * grid_v.ny
            nt = (num.nLS + 1) * ni + nb
            AvS = spzeros(nt, nt)
            AvL = spzeros(nt, nt)
            BvS = spzeros(nt, nt)
            BvL = spzeros(nt, nt)

            ni = grid_u.nx * grid_u.ny + grid_v.nx * grid_v.ny
            nb = 2 * grid_u.nx + 2 * grid_u.ny + 2 * grid_v.nx + 2 * grid_v.ny
            nt = (nLS - nNavier + 1) * ni + nNavier * grid.nx * grid.ny + nb
            AuvS = spzeros(nt, nt)
            AuvL = spzeros(nt, nt)
            BuvS = spzeros(nt, nt)
            BuvL = spzeros(nt, nt)

            if !navier
                _ = FE_set_momentum(
                    BC_int, num, grid_u, opC_uS,
                    AuS, BuS,
                    iRe.*Lum1_S, iRe.*bc_Lum1_S, iRe.*bc_Lum1_b_S, Mum1_S, BC_uS,
                    true
                )
                _ = FE_set_momentum(
                    BC_int, num, grid_v, opC_vS,
                    AvS, BvS,
                    iRe.*Lvm1_S, iRe.*bc_Lvm1_S, iRe.*bc_Lvm1_b_S, Mvm1_S, BC_vS,
                    true
                )
            else
                _ = FE_set_momentum_coupled(
                    BC_int, num, grid, grid_u, grid_v,
                    opC_pS, opC_uS, opC_vS,
                    AuvS, BuvS,
                    iRe.*Lum1_S, iRe.*bc_Lum1_S, iRe.*bc_Lum1_b_S, Mum1_S, BC_uS,
                    iRe.*Lvm1_S, iRe.*bc_Lvm1_S, iRe.*bc_Lvm1_b_S, Mvm1_S, BC_vS,
                    true
                )
            end
            a0_p = []
            for i in 1:num.nLS
                push!(a0_p, zeros(grid))
            end
            _ = set_poisson(
                BC_int, num, grid, a0_p, opC_pS, opC_uS, opC_vS,
                AϕS, Lpm1_S, bc_Lpm1_S, bc_Lpm1_b_S, BC_pS,
                true
            )

            set_heat!(
                BC_int[1], num, grid, opC_TS, LS[1].geoS, phS, θd, BC_TS, LS[1].MIXED, LS[1].geoS.projection,
                ATS, BTS,
                opS, grid_u, grid_u.LS[1].geoS, grid_v, grid_v.LS[1].geoS,
                periodic_x, periodic_y, heat_convection, true, BC_int
            )

            if !navier
                _ = FE_set_momentum(
                    BC_int, num, grid_u, opC_uL,
                    AuL, BuL,
                    iRe.*Lum1_L, iRe.*bc_Lum1_L, iRe.*bc_Lum1_b_L, Mum1_L, BC_uL,
                    true
                )
                _ = FE_set_momentum(
                    BC_int, num, grid_v, opC_vL,
                    AvL, BvL,
                    iRe.*Lvm1_L, iRe.*bc_Lvm1_L, iRe.*bc_Lvm1_b_L, Mvm1_L, BC_vL,
                    true
                )
            else
                _ = FE_set_momentum_coupled(
                    BC_int, num, grid, grid_u, grid_v,
                    opC_pL, opC_uL, opC_vL,
                    AuvL, BuvL,
                    iRe.*Lum1_L, iRe.*bc_Lum1_L, iRe.*bc_Lum1_b_L, Mum1_L, BC_uL,
                    iRe.*Lvm1_L, iRe.*bc_Lvm1_L, iRe.*bc_Lvm1_b_L, Mvm1_L, BC_vL,
                    true
                )
            end
            a0_p = []
            for i in 1:num.nLS
                push!(a0_p, zeros(grid))
            end
            _ = set_poisson(
                BC_int, num, grid, a0_p, opC_pL, opC_uL, opC_vL,
                AϕL, Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, BC_pL,
                true
            )

            set_heat!(
                BC_int[1], num, grid, opC_TL, LS[1].geoL, phL, θd, BC_TL, LS[1].MIXED, LS[1].geoL.projection,
                ATL, BTL,
                opL, grid_u, grid_u.LS[1].geoL, grid_v, grid_v.LS[1].geoL,
                periodic_x, periodic_y, heat_convection, true, BC_int
            )
        end
    else
        error("Unknown time scheme. Available options are ForwardEuler and CrankNicolson")
    end

    if heat_convection
        NB_indices = update_all_ls_data(num, grid, grid_u, grid_v, BC_int, periodic_x, periodic_y)
    end

    V0S = volume(LS[end].geoS)
    V0L = volume(LS[end].geoL)

    # Force in x
    if adaptative_t && flapping
        for iLS in 1:nLS
            if is_flapping(BC_int[iLS])
                num.τ = min(CFL*Δ^2*Re, CFL*Δ/max(
                    abs.(V)..., sqrt(max(abs.(grid_u.V)...)^2 + max(abs.(grid_v.V)...)^2),
                    abs.(phL.u)..., abs.(phL.v)..., abs.(phS.u)..., abs.(phS.v)..., 
                    BC_int[iLS].β * BC_int[iLS].KC),
                    1.0 / BC_int[iLS].β / BC_int[iLS].βmult
                )
                CFL_sc = num.τ / Δ^2
            end
        end
    end

    current_t = 0.

    while current_i < max_iterations + 1

        if !stefan
            V .= speed*ones(ny, nx)
        end

        if heat
            if heat_solid_phase
                kill_dead_cells!(phS.T, grid, LS[1].geoS)
                veci(phS.TD,grid,1) .= vec(phS.T)
                rhs = set_heat!(
                    BC_int[1], num, grid, opC_TS, LS[1].geoS, phS, θd, BC_TS, LS[1].MIXED, LS[1].geoS.projection,
                    ATS, BTS,
                    opS, grid_u, grid_u.LS[1].geoS, grid_v, grid_v.LS[1].geoS,
                    periodic_x, periodic_y, heat_convection, advection, BC_int
                )
                mul!(rhs, BTS, phS.TD, 1.0, 1.0)

                phS.TD .= ATS \ rhs
                phS.T .= reshape(veci(phS.TD,grid,1), grid)
            end
            if heat_liquid_phase
                kill_dead_cells!(phL.T, grid, LS[1].geoL)
                veci(phL.TD,grid,1) .= vec(phL.T)
                rhs = set_heat!(
                    BC_int[1], num, grid, opC_TL, LS[1].geoL, phL, θd, BC_TL, LS[1].MIXED, LS[1].geoL.projection,
                    ATL, BTL,
                    opL, grid_u, grid_u.LS[1].geoL, grid_v, grid_v.LS[1].geoL,
                    periodic_x, periodic_y, heat_convection, advection, BC_int
                )
                mul!(rhs, BTL, phL.TD, 1.0, 1.0)

                phL.TD .= ATL \ rhs
                phL.T .= reshape(veci(phL.TD,grid,1), grid)
            end
        end

        for iLS in 1:nLS
            if is_stefan(BC_int[iLS])
                update_stefan_velocity(num, grid, iLS, LS[iLS].u, phS.T, phL.T, periodic_x, periodic_y, St, Vmean)
            elseif is_fs(BC_int[iLS])
                update_free_surface_velocity(num, grid_u, grid_v, iLS, phL.uD, phL.vD, periodic_x, periodic_y)
            elseif is_flapping(BC_int[iLS])
                if BC_int[iLS].free
                    grid_u.V .+= num.τ .* D ./ (BC_int[iLS].ρs .* volume(LS[end].geoS))
                end
                grid_v.V .= BC_int[iLS].KC .* BC_int[iLS].β .* sin(2π .* BC_int[iLS].β .* current_t)
                xc[current_i+1] = xc[current_i] + num.τ * grid_u.V[1,1]
                yc[current_i+1] = yc[current_i] + num.τ * grid_v.V[1,1]
                grid.LS[iLS].u .= sqrt.((grid.x .- xc[current_i+1]) .^ 2 .+ ((grid.y .- yc[current_i+1]) ./ num.A) .^ 2) .- num.R .* ones(grid);
                println("xc = $(xc[current_i+1]) | vx = $(grid_u.V[1,1])")
                println("yc = $(yc[current_i+1]) | vy = $(grid_v.V[1,1])")
            end
        end

        if verbose && adaptative_t && (current_i-1)%show_every == 0 
            println("τ = $(num.τ)")
        end

        if advection
            # print("advection")
            for (iLS, bc) in enumerate(BC_int)
                if is_stefan(bc)
                    # IIOE_normal!(grid, LS[iLS].A, LS[iLS].B, LS[iLS].u, V, CFL_sc, periodic_x, periodic_y)
                    # LS[iLS].u .= reshape(gmres(LS[iLS].A, LS[iLS].B * vec(LS[iLS].u)), grid)
                    LS[iLS].u .= fwd.u[iLS,end-current_i,:,:]
                    # u .= sqrt.((x .- current_i*Δ/1).^ 2 + y .^ 2) - (0.5) * ones(nx, ny);
                elseif is_fs(bc)
                    # rhs_LS .= 0.0
                    # LS[iLS].A.nzval .= 0.0
                    # LS[iLS].B.nzval .= 0.0
                    # IIOE!(grid, grid_u, grid_v, LS[iLS].A, LS[iLS].B, θ_out, num.τ, periodic_x, periodic_y)
                    # BC_LS!(grid, LS[iLS].u, LS[iLS].A, LS[iLS].B, rhs_LS, BC_u)
                    # utmp .= reshape(gmres(LS[iLS].A, LS[iLS].B * vec(LS[iLS].u) .+ rhs_LS), grid)

                    # rhs_LS .= 0.0
                    # S2IIOE!(grid, grid_u, grid_v, LS[iLS].A, LS[iLS].B, utmp, LS[iLS].u, θ_out, num.τ, periodic_x, periodic_y)
                    # BC_LS!(grid, LS[iLS].u, LS[iLS].A, LS[iLS].B, rhs_LS, BC_u)
                    # LS[iLS].u .= reshape(gmres(LS[iLS].A, LS[iLS].B * vec(LS[iLS].u) .+ rhs_LS), grid)

                    # Project velocities to the normal and use advecion scheme for advection just
                    # in the normal direction
                    tmpVx = zeros(grid)
                    tmpVy = zeros(grid)
                    V .= 0.0
                    @inbounds @threads for II in grid.LS[iLS].MIXED
                        cap1 = grid_u.LS[iLS].geoL.cap[II,5]
                        cap3 = grid_u.LS[iLS].geoL.cap[δx⁺(II),5]
                        tmpVx[II] = (grid_u.V[II] * cap1 + grid_u.V[δx⁺(II)] * cap3) / (cap1 + cap3 + eps(0.01))

                        cap2 = grid_v.LS[iLS].geoL.cap[II,5]
                        cap4 = grid_v.LS[iLS].geoL.cap[δy⁺(II),5]
                        tmpVy[II] = (grid_v.V[II] * cap2 + grid_v.V[δy⁺(II)] * cap4) / (cap2 + cap4 + eps(0.01))

                        tmpV = sqrt(tmpVx[II]^2 + tmpVy[II]^2)
                        β = atan(tmpVy[II], tmpVx[II])
                        if grid.LS[iLS].α[II] > 0.0 && β < 0.0
                            β += 2π
                        end
                        if grid.LS[iLS].α[II] < 0.0 && β > 0.0
                            β -= 2π
                        end

                        V[II] = tmpV * cos(β - grid.LS[iLS].α[II])
                    end

                    i_ext, l_ext, b_ext, r_ext, t_ext = indices_extension(grid, grid.LS[iLS], grid.ind.inside, periodic_x, periodic_y)
                    field_extension!(grid, grid.LS[iLS].u, V, i_ext, l_ext, b_ext, r_ext, t_ext, num.NB, periodic_x, periodic_y)

                    rhs_LS .= 0.0
                    IIOE_normal!(grid, LS[iLS].A, LS[iLS].B, LS[iLS].u, V, CFL_sc, periodic_x, periodic_y)
                    BC_LS!(grid, LS[iLS].u, LS[iLS].A, LS[iLS].B, rhs_LS, BC_u)
                    BC_LS_interior!(num, grid, iLS, LS[iLS].A, LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
                    LS[iLS].u .= reshape(gmres(LS[iLS].A, LS[iLS].B * vec(LS[iLS].u) .+ rhs_LS), grid)
                end
            end
            if analytical
                u[ind.b_top[1]] .= sqrt.(x[ind.b_top[1]] .^ 2 + y[ind.b_top[1]] .^ 2) .- (num.R + speed*current_i*num.τ);
                u[ind.b_bottom[1]] .= sqrt.(x[ind.b_bottom[1]] .^ 2 + y[ind.b_bottom[1]] .^ 2) .- (num.R + speed*current_i*num.τ);
                u[ind.b_left[1]] .= sqrt.(x[ind.b_left[1]] .^ 2 + y[ind.b_left[1]] .^ 2) .- (num.R + speed*current_i*num.τ);
                u[ind.b_right[1]] .= sqrt.(x[ind.b_right[1]] .^ 2 + y[ind.b_right[1]] .^ 2) .- (num.R + speed*current_i*num.τ);
            elseif nb_reinit > 0
                if auto_reinit && (current_i-1)%num.reinit_every == 0
                    for iLS in 1:nLS
                        if !is_wall(BC_int[iLS])
                            ls_rg = rg(num, grid, LS[iLS].u, periodic_x, periodic_y, BC_int)[1]
                            println(ls_rg)
                            if ls_rg >= δreinit || current_i == 1
                                println("yes")
                                RK2_reinit!(ls_scheme, grid, ind, iLS, LS[iLS].u, nb_reinit, periodic_x, periodic_y, BC_u, BC_int)
                            end
                        end
                    end
                elseif (current_i-1)%num.reinit_every == 0
                    for iLS in 1:nLS
                        if !is_wall(BC_int[iLS])
                            RK2_reinit!(ls_scheme, grid, ind, iLS, LS[iLS].u, nb_reinit, periodic_x, periodic_y, BC_u, BC_int)
                        end
                    end
                else
                    for iLS in 1:nLS
                        if !is_wall(BC_int[iLS])
                            RK2_reinit!(ls_scheme, grid, ind, iLS, LS[iLS].u, nb_reinit, periodic_x, periodic_y, BC_u, BC_int, true)
                        end
                    end
                end
            end
            # numerical breakup
            if free_surface && breakup
                count = breakup_f(LS[1].u, nx, ny, dx, dy, periodic_x, periodic_y, NB_indices, 1e-5)
                if count > 0
                    RK2_reinit!(ls_scheme, grid, ind, 1, LS[1].u, nb_reinit, periodic_x, periodic_y, BC_u, BC_int)
                end
            end
        end

        if verbose
            if (current_i-1)%show_every == 0
                printstyled(color=:green, @sprintf "\n Current iteration : %d (%d%%) | t = %.2e \n" (current_i-1) 100*(current_i-1)/max_iterations current_t)
                if heat && length(LS[end].MIXED) != 0
                    print(@sprintf "V_mean = %.2e  V_max = %.2e  V_min = %.2e\n" mean(V[LS[1].MIXED]) findmax(V[LS[1].MIXED])[1] findmin(V[LS[1].MIXED])[1])
                    print(@sprintf "κ_mean = %.2e  κ_max = %.2e  κ_min = %.2e\n" mean(LS[1].κ[LS[1].MIXED]) findmax(LS[1].κ[LS[1].MIXED])[1] findmin(LS[1].κ[LS[1].MIXED])[1])
                elseif advection && length(LS[end].MIXED) != 0
                    V_mean = mean([mean(grid_u.V[LS[1].MIXED]), mean(grid_v.V[LS[1].MIXED])])
                    V_max = max(findmax(grid_u.V[LS[1].MIXED])[1], findmax(grid_v.V[LS[1].MIXED])[1])
                    V_min = min(findmin(grid_u.V[LS[1].MIXED])[1], findmin(grid_v.V[LS[1].MIXED])[1])
                    print(@sprintf "Vol_ratio = %.3f%%\n" (volume(LS[end].geoL) / V0L * 100))
                    print(@sprintf "V_mean = %.2e  V_max = %.2e  V_min = %.2e\n" V_mean V_max V_min)
                    print(@sprintf "κ_mean = %.2e  κ_max = %.2e  κ_min = %.2e\n" mean(LS[1].κ[LS[1].MIXED]) findmax(LS[1].κ[LS[1].MIXED])[1] findmin(LS[1].κ[LS[1].MIXED])[1])
                end
                if navier_stokes
                    if ns_solid_phase
                        normuS = norm(phS.u)
                        normvS = norm(phS.v)
                        normpS = norm(phS.p.*num.τ)
                        print("$(@sprintf("norm(uS) %.6e", normuS))\t$(@sprintf("norm(vS) %.6e", normvS))\t$(@sprintf("norm(pS) %.6e", normpS))\n")
                    end
                    if ns_liquid_phase
                        normuL = norm(phL.u)
                        normvL = norm(phL.v)
                        normpL = norm(phL.p.*num.τ)
                        print("$(@sprintf("norm(uL) %.6e", normuL))\t$(@sprintf("norm(vL) %.6e", normvL))\t$(@sprintf("norm(pL) %.6e", normpL))\n")
                    end
                end
            end
        end


        if levelset && (advection || current_i<2)
            # print("level set and advection")
            NB_indices = update_all_ls_data(num, grid, grid_u, grid_v, BC_int, periodic_x, periodic_y)

            if heat && heat_convection && heat_liquid_phase && navier_stokes && ns_liquid_phase
                fwd.RB[1,current_i] = current_t
                # h = 0.
                # for jj in grid.LS[1].MIXED
                #     h += grid.y[jj] + num.Δ * grid.LS[1].mid_point[jj].y + 0.5
                # end
                # tmp = h / length(grid.LS[1].MIXED) 
                sum = 0
                count = 0
                for i in 1:grid.nx
                    for j in 1:(grid.ny-1)
                        diff = grid.LS[1].u[j,i]*grid.LS[1].u[j+1,i]
                        if diff < 0
                            pos = grid.y[j+1,1] + grid.LS[1].u[j+1,i]
                            sum += pos 
                            count += 1
                        end
                    end
                end
                tmp = (sum / count + 0.5)
                fwd.RB[2,current_i] = tmp
                fwd.RB[3,current_i] = Ra * (1. - num.θd) * tmp^3     
                if verbose
                    if (current_i-1)%show_every == 0 
                        print("$(@sprintf("height %g", fwd.RB[2,current_i]))\t$(@sprintf("Ra_eff %g", fwd.RB[3,current_i]))\n")
                    end
                end
                if tmp > cutoff_length
                    return nothing
                end
            end

            LS[end].geoL.fresh .= false
            LS[end].geoS.fresh .= false
            grid_u.LS[end].geoL.fresh .= false
            grid_u.LS[end].geoS.fresh .= false
            grid_v.LS[end].geoL.fresh .= false
            grid_v.LS[end].geoS.fresh .= false

            get_fresh_cells!(grid, grid.LS[end].geoS, Mm1_S, grid.ind.all_indices)
            get_fresh_cells!(grid, grid.LS[end].geoL, Mm1_L, grid.ind.all_indices)
            get_fresh_cells!(grid_u, grid_u.LS[end].geoS, Mum1_S, grid_u.ind.all_indices)
            get_fresh_cells!(grid_u, grid_u.LS[end].geoL, Mum1_L, grid_u.ind.all_indices)
            get_fresh_cells!(grid_v, grid_v.LS[end].geoS, Mvm1_S, grid_v.ind.all_indices)
            get_fresh_cells!(grid_v, grid_v.LS[end].geoL, Mvm1_L, grid_v.ind.all_indices)

            FRESH_L_u = findall(grid_u.LS[end].geoL.fresh)
            FRESH_S_u = findall(grid_u.LS[end].geoS.fresh)
            FRESH_L_v = findall(grid_v.LS[end].geoL.fresh)
            FRESH_S_v = findall(grid_v.LS[end].geoS.fresh)

            if navier_stokes
                init_fresh_cells!(grid_u, veci(phS.uD,grid_u,1), veci(phS.uD,grid_u,1),
                    grid_u.LS[end].geoS.projection, FRESH_S_u, periodic_x, periodic_y)
                init_fresh_cells!(grid_v, veci(phS.vD,grid_v,1), veci(phS.vD,grid_v,1),
                    grid_v.LS[end].geoS.projection, FRESH_S_v, periodic_x, periodic_y)
                init_fresh_cells!(grid_u, veci(phS.uD,grid_u,2), veci(phS.uD,grid_u,1),
                    grid_u.LS[end].geoS.projection, FRESH_S_u, periodic_x, periodic_y)
                init_fresh_cells!(grid_v, veci(phS.vD,grid_v,2), veci(phS.vD,grid_v,1),
                    grid_v.LS[end].geoS.projection, FRESH_S_v, periodic_x, periodic_y)

                init_fresh_cells!(grid_u, veci(phL.uD,grid_u,1), veci(phL.uD,grid_u,1),
                    grid_u.LS[end].geoL.projection, FRESH_L_u, periodic_x, periodic_y)
                init_fresh_cells!(grid_v, veci(phL.vD,grid_v,1), veci(phL.vD,grid_v,1),
                    grid_v.LS[end].geoL.projection, FRESH_L_v, periodic_x, periodic_y)
                init_fresh_cells!(grid_u, veci(phL.uD,grid_u,2), veci(phL.uD,grid_u,1),
                    grid_u.LS[end].geoL.projection, FRESH_L_u, periodic_x, periodic_y)
                init_fresh_cells!(grid_v, veci(phL.vD,grid_v,2), veci(phL.vD,grid_v,1),
                    grid_v.LS[end].geoL.projection, FRESH_L_v, periodic_x, periodic_y)
            end

            if iszero(current_i%save_every) || current_i==max_iterations
                snap = current_i÷save_every+1
                if save_radius
                    radius[snap] = find_radius(grid, LS[1])
                end
                if hill
                    a = zeros(length(LS[1].MIXED))
                    for i in eachindex(LS[1].MIXED)
                        a[i] = LS[1].geoL.projection[LS[1].MIXED[i]].pos.y
                    end
                    radius[snap] = mean(a)
                end
                if save_length
                    fwd.length[snap] = arc_length2(LS[1].geoS.projection, LS[1].MIXED)
                end
            end
        end

        if navier_stokes
            if !advection
                no_slip_condition!(num, grid, grid_u, grid_u.LS[1], grid_v, grid_v.LS[1], periodic_x, periodic_y)
                # grid_u.V .= Δ / (1 * num.τ)
                # grid_v.V .= 0.0
            end

            if ns_solid_phase
                geoS = [LS[iLS].geoS for iLS in 1:_nLS]
                geo_uS = [grid_u.LS[iLS].geoS for iLS in 1:_nLS]
                geo_vS = [grid_v.LS[iLS].geoS for iLS in 1:_nLS]
                Lpm1_S, bc_Lpm1_S, bc_Lpm1_b_S, Lum1_S, bc_Lum1_S, bc_Lum1_b_S, Lvm1_S, bc_Lvm1_S, bc_Lvm1_b_S,Mm1_S, Mum1_S, Mvm1_S, Cum1S, Cvm1S = pressure_projection!(
                    time_scheme, BC_int,
                    num, grid, geoS, grid_u, geo_uS, grid_v, geo_vS, phS,
                    BC_uS, BC_vS, BC_pS,
                    opC_pS, opC_uS, opC_vS, opS,
                    AuS, BuS, AvS, BvS, AϕS, AuvS, BuvS,
                    Lpm1_S, bc_Lpm1_S, bc_Lpm1_b_S, Lum1_S, bc_Lum1_S, bc_Lum1_b_S, Lvm1_S, bc_Lvm1_S, bc_Lvm1_b_S,
                    Cum1S, Cvm1S, Mum1_S, Mvm1_S,
                    periodic_x, periodic_y, ns_advection, advection, current_i, Ra, navier
                )
            end
            if ns_liquid_phase
                geoL = [LS[iLS].geoL for iLS in 1:_nLS]
                geo_uL = [grid_u.LS[iLS].geoL for iLS in 1:_nLS]
                geo_vL = [grid_v.LS[iLS].geoL for iLS in 1:_nLS]
                Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, Lum1_L, bc_Lum1_L, bc_Lum1_b_L, Lvm1_L, bc_Lvm1_L, bc_Lvm1_b_L, Mm1_L, Mum1_L, Mvm1_L, Cum1L, Cvm1L = pressure_projection!(
                    time_scheme, BC_int,
                    num, grid, geoL, grid_u, geo_uL, grid_v, geo_vL, phL,
                    BC_uL, BC_vL, BC_pL,
                    opC_pL, opC_uL, opC_vL, opL,
                    AuL, BuL, AvL, BvL, AϕL, AuvL, BuvL,
                    Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, Lum1_L, bc_Lum1_L, bc_Lum1_b_L, Lvm1_L, bc_Lvm1_L, bc_Lvm1_b_L,
                    Cum1L, Cvm1L, Mum1_L, Mvm1_L,
                    periodic_x, periodic_y, ns_advection, advection, current_i, Ra, navier
                )
                # if current_i == 1
                #     phL.u .= -0.5 .* grid_u.y .+ getproperty.(grid_u.LS[1].geoL.centroid, :y) .* grid_u.dy
                #     phL.v .= 0.5 .* grid_v.x .+ getproperty.(grid_v.LS[1].geoL.centroid, :x) .* grid_v.dx
                #     phL.u[grid_u.LS[1].SOLID] .= 0.0
                #     phL.v[grid_v.LS[1].SOLID] .= 0.0
                # end
                # linear_advection!(
                #     num, grid, LS[1].geoL, grid_u, grid_u.LS[1].geoL, grid_v, grid_v.LS[1].geoL, phL,
                #     BC_uL, BC_vL, opL
                # )
            end
        end

        cD, cL, D, L = force_coefficients!(num, grid, grid_u, grid_v, opL, fwd, phL; step = current_i+1, saveCoeffs = false)

        current_t += num.τ
        # if iszero(current_i%save_every) || current_i==max_iterations
        #     snap = current_i÷save_every+1
        #     if current_i==max_iterations
        #         snap = size(fwd.T,1)
        #     end
        #     fwd.t[snap] = current_t
        #     @views fwd.V[snap,:,:] .= V
        #     for iLS in 1:_nLS
        #         @views fwd.u[iLS,snap,:,:] .= LS[iLS].u
        #         @views fwd.ux[iLS,snap,:,:] .= grid_u.LS[iLS].u
        #         @views fwd.uy[iLS,snap,:,:] .= grid_v.LS[iLS].u
        #         @views fwd.κ[iLS,snap,:,:] .= LS[iLS].κ
        #     end

        #     if heat_solid_phase && heat_liquid_phase
        #         @views fwd.T[snap,:,:] .= phL.T.*LS[end].geoL.cap[:,:,5] .+ phS.T.*LS[end].geoS.cap[:,:,5]
        #     end
        #     if heat_solid_phase
        #         @views fwd.T[snap,:,:] .= phL.T + phS.T
        #         @views fwdS.T[snap,:,:] .= phS.T
        #         @views fwdS.TD[snap,:] .= phS.TD
        #     end
        #     if heat_liquid_phase
        #         @views fwd.T[snap,:,:] .= phL.T + phS.T
        #         @views fwdL.T[snap,:,:] .= phL.T
        #         @views fwdL.TD[snap,:] .= phL.TD
        #     end

        #     if ns_solid_phase
        #         @views fwdS.p[snap,:,:] .= phS.p
        #         @views fwdS.pD[snap,:] .= phS.pD
        #         @views fwdS.ϕ[snap,:,:] .= phS.ϕ
        #         @views fwdS.u[snap,:,:] .= phS.u
        #         @views fwdS.v[snap,:,:] .= phS.v
        #         @views fwdS.ucorrD[snap,:,:] .= phS.ucorrD
        #         @views fwdS.vcorrD[snap,:,:] .= phS.vcorrD
        #     end
        #     if ns_liquid_phase
        #         @views fwdL.p[snap,:,:] .= phL.p
        #         @views fwdL.pD[snap,:] .= phL.pD
        #         @views fwdL.ϕ[snap,:,:] .= phL.ϕ
        #         @views fwdL.u[snap,:,:] .= phL.u
        #         @views fwdL.v[snap,:,:] .= phL.v
        #         @views fwdL.ucorrD[snap,:,:] .= phL.ucorrD
        #         @views fwdL.vcorrD[snap,:,:] .= phL.vcorrD
        #         # @views fwd.Cd[snap] = cD
        #         # @views fwd.Cl[snap] = cL
        #     end
        #     if advection
        #         fwdS.Vratio[snap] = volume(LS[end].geoS) / V0S
        #         fwdL.Vratio[snap] = volume(LS[end].geoL) / V0L
        #     end
        # end
        # @views fwd.Cd[current_i+1] = cD
        # @views fwd.Cl[current_i+1] = cL

        if (any(isnan, phL.uD) || any(isnan, phL.vD) || any(isnan, phL.TD) || any(isnan, phS.uD) || any(isnan, phS.vD) || any(isnan, phS.TD) ||
            norm(phL.u) > 1e8 || norm(phS.u) > 1e8 || norm(phL.T) > 1e8 || norm(phS.T) > 1e8)
            println(@sprintf "\n CRASHED after %d iterations \n" current_i)
            if flapping
                return xc, yc
            else
                return nothing
            end
        end

        current_i += 1

        if adaptative_t && flapping
            for iLS in 1:nLS
                if is_flapping(BC_int[iLS])
                    num.τ = min(CFL*Δ^2*Re, CFL*Δ/max(
                        abs.(V)..., sqrt(max(abs.(grid_u.V)...)^2 + max(abs.(grid_v.V)...)^2),
                        abs.(phL.u)..., abs.(phL.v)..., abs.(phS.u)..., abs.(phS.v)..., 
                        BC_int[iLS].β * BC_int[iLS].KC),
                        1.0 / BC_int[iLS].β / BC_int[iLS].βmult
                    )
                    CFL_sc = num.τ / Δ^2
                end
            end
        elseif adaptative_t
            num.τ = min(CFL*Δ^2*Re, CFL*Δ/max(
                abs.(V)..., abs.(grid_u.V)..., abs.(grid_v.V)..., 
                abs.(phL.u)..., abs.(phL.v)..., abs.(phS.u)..., abs.(phS.v)...)
            )
            CFL_sc = num.τ / Δ^2
        end
    end

    if verbose
        try
            printstyled(color=:blue, @sprintf "\n Final iteration : %d (%d%%) | t = %.2e \n" (current_i-1) 100*(current_i-1)/max_iterations current_t)
            if stefan && advection
                print(@sprintf "V_mean = %.2e  V_max = %.2e  V_min = %.2e  V_stdev = %.5f\n" mean(V[LS[1].MIXED]) findmax(V[LS[1].MIXED])[1] findmin(V[LS[1].MIXED])[1] std(V[LS[1].MIXED]))
                print(@sprintf "κ_mean = %.2e  κ_max = %.2e  κ_min = %.2e  κ_stdev = %.5f\n" mean(LS[1].κ[LS[1].MIXED]) findmax(LS[1].κ[LS[1].MIXED])[1] findmin(LS[1].κ[LS[1].MIXED])[1] std(LS[1].κ[LS[1].MIXED]))
            end
            if free_surface && advection
                print(@sprintf "Vol_ratio = %.3f%%\n" (volume(LS[end].geoL) / V0L * 100))
                print(@sprintf "V_mean = %.2e  V_max = %.2e  V_min = %.2e  V_stdev = %.5f\n" mean(V[LS[1].MIXED]) findmax(V[LS[1].MIXED])[1] findmin(V[LS[1].MIXED])[1] std(V[LS[1].MIXED]))
                print(@sprintf "κ_mean = %.2e  κ_max = %.2e  κ_min = %.2e  κ_stdev = %.5f\n" mean(LS[1].κ[LS[1].MIXED]) findmax(LS[1].κ[LS[1].MIXED])[1] findmin(LS[1].κ[LS[1].MIXED])[1] std(LS[1].κ[LS[1].MIXED]))
            end
            if navier_stokes
                if ns_solid_phase
                    normuS = norm(phS.u)
                    normvS = norm(phS.v)
                    normpS = norm(phS.p.*num.τ)
                    print("$(@sprintf("norm(uS) %.6e", normuS))\t$(@sprintf("norm(vS) %.6e", normvS))\t$(@sprintf("norm(pS) %.6e", normpS))\n")
                end
                if ns_liquid_phase
                    normuL = norm(phL.u)
                    normvL = norm(phL.v)
                    normpL = norm(phL.p.*num.τ)
                    print("$(@sprintf("norm(uL) %.6e", normuL))\t$(@sprintf("norm(vL) %.6e", normvL))\t$(@sprintf("norm(pL) %.6e", normpL))\n")
                end
            end
            print("\n\n")
        catch
            @show (length(LS[end].MIXED))
        end
    end

    if levelset && (save_radius || hill)
        return radius
    elseif flapping
        return xc, yc
    else
        return nothing
    end
end
























function run_backward(num, grid, opS, opL, fwd, adj;
    periodic_x = false,
    periodic_y = false,
    BC_TS = Boundaries(
        left = Boundary(),
        right = Boundary(),
        bottom = Boundary(),
        top = Boundary()),
    BC_TL = Boundaries(
        left = Boundary(),
        right = Boundary(),
        bottom = Boundary(),
        top = Boundary()),
    BC_u = Boundaries(
        left = Boundary(),
        right = Boundary(),
        bottom = Boundary(),
        top = Boundary()),
    stefan = false,
    advection = false,
    heat = false,
    heat_convection = false,
    liquid_phase = false,
    solid_phase = false,
    hill = false,
    Vmean = false,
    levelset = true,
    speed = 0,
    analytical = false,
    verbose = false,
    show_every = 100,
    save_length = false,
    save_radius = false
    )

    @unpack L0, A, N, θd, ϵ_κ, ϵ_V, T_inf, τ, L0, NB, max_iterations, current_i, reinit_every, nb_reinit, ϵ, m, θ₀, aniso = num
    @unpack nx, ny, ind, faces, geoS, geoL, mid_point = grid
    @unpack all_indices, inside, b_left, b_bottom, b_right, b_top = ind
    @unpack usave, TSsave, TLsave, Tsave, Vsave, κsave = fwd
    @unpack iso, u, TS, TL, DTS, DTL, κ, V = adj

    local MIXED; local SOLID; local LIQUID;
    local WAS_SOLID; local WAS_LIQUID;
    local NB_indices_base; local NB_indices;
    local FRESH_L; local FRESH_S;

    if periodic_x
        BC_u.left.ind = idx.periodicL;
        BC_u.right.ind = idx.periodicR;
        BC_u.left.f = BC_u.right.f = periodic
    else
        BC_u.left.ind = b_left;
        BC_u.right.ind = b_right;
    end

    if periodic_y
        BC_u.bottom.ind = idx.periodicB;
        BC_u.top.ind = idx.periodicT;
        BC_u.bottom.f = BC_u.top.f = periodic
    else
        BC_u.bottom.ind = b_bottom;
        BC_u.top.ind = b_top;
    end

    current_i = max_iterations + 1

    if levelset
        marching_squares!(num, grid, u, periodic_x, periodic_y)

        NB_indices_base = get_NB_width_indices_base(NB)

        MIXED, SOLID, LIQUID = get_cells_indices(iso, inside)
        NB_indices = get_NB_width(MIXED, NB_indices_base)

        get_iterface_location!(grid, MIXED)
        get_curvature(num, grid, u, MIXED, periodic_x, periodic_y)
    elseif !levelset
        MIXED = [CartesianIndex(-1,-1)]
    end

    HS = zeros(ny, nx)
    for II in vcat(b_left[1], b_bottom[1], b_right[1], b_top[1])
        HS[II] = distance(mid_point[II], geoS.centroid[II], dx[II], dy[II])
    end
    
    HL = zeros(ny, nx)
    for II in vcat(b_left[1], b_bottom[1], b_right[1], b_top[1])
        HL[II] = distance(mid_point[II], geoL.centroid[II], dx[II], dy[II])
    end

    DTS .= θd
    DTL .= θd
    bcS = similar(DTS)
    bcL = similar(DTL)
    apply_curvature(bcT, DTS, num, grid, all_indices)
    apply_curvature(bcT, DTL, num, grid, all_indices)
    if aniso
        apply_anisotropy(bcS, DTS, MIXED, num, grid, geoS.projection)
        apply_anisotropy(bcL, DTL, MIXED, num, grid, geoS.projection)
    end
    bcSx, bcSy = set_bc_bnds(dir, bcS, HS, BC_TS)
    bcLx, bcLy = set_bc_bnds(dir, bcL, HL, BC_TL)

    laplacian!(dir, opS.LT, opS.CUTT, bcSx, bcSy, geoS.dcap, ny, BC_TS, inside, LIQUID,
                MIXED, b_left[1], b_bottom[1], b_right[1], b_top[1])
    laplacian!(dir, opL.LT, opL.CUTT, bcLx, bcLy, geoL.dcap, ny, BC_TL, inside, SOLID,
                MIXED, b_left[1], b_bottom[1], b_right[1], b_top[1])

    while current_i > 1

        if heat
            opL.CUTT .= zeros(nx*ny)
            opS.CUTT .= zeros(nx*ny)

            try
                if solid_phase
<<<<<<< HEAD
                    TS .= reshape(gmres(AS,(BS*vec(TS) + 2.0*SCUT)), (n,n))
                    #TS .= reshape(AS\(BS*vec(TS) + 2.0*SCUT), (n,n))
                    bcs!(TS, BC_TS.left, Δ)
                    bcs!(TS, BC_TS.right, Δ)
                    bcs!(TS, BC_TS.bottom, Δ)
                    bcs!(TS, BC_TS.top, Δ)
                end
                if liquid_phase
                    #@show (TL[b_top[1]])
                    TL .= reshape(gmres(AL,(BL*vec(TL) + 2.0*LCUT)), (n,n))
                    # TL .= reshape(AL\(BL*vec(TL) + 2.0*LCUT), (n,n))
                    bcs!(TL, BC_TL.left, Δ)
                    bcs!(TL, BC_TL.right, Δ)
                    bcs!(TL, BC_TL.bottom, Δ)
                    bcs!(TL, BC_TL.top, Δ)
=======
                    HS .= 0.
                    for II in vcat(ind.b_left[1], ind.b_bottom[1], ind.b_right[1], ind.b_top[1])
                        HS[II] = distance(mid_point[II], geoS.centroid[II], dx[II], dy[II])
                    end

                    DTS .= θd
                    apply_curvature(bcT, DTS, num, grid, all_indices)
                    if aniso
                        apply_anisotropy(bcS, DTS, MIXED, num, grid, geoS.projection)
                    end
                    bcSx, bcSy = set_bc_bnds(dir, bcS, HS, BC_TS)

                    laplacian!(dir, opS.LT, opS.CUTT, bcSx, bcSy, geoS.dcap, ny, BC_TS, inside, LIQUID,
                                MIXED, b_left[1], b_bottom[1], b_right[1], b_top[1])
                    crank_nicolson!(num, grid, geoS, opS)
                    TS .= reshape(gmres(opS.A,(opS.B*vec(TS) + 2.0*τ*opS.CUTT)), (ny,nx))
                end
                if liquid_phase
                    HL .= 0.
                    for II in vcat(b_left[1], b_bottom[1], b_right[1], b_top[1])
                        HL[II] = distance(mid_point[II], geoL.centroid[II], dx[II], dy[II])
                    end

                    DTL .= θd
                    apply_curvature(bcT, DTL, num, grid, all_indices)
                    if aniso
                        apply_anisotropy(bcL, DTL, MIXED, num, grid, geoS.projection)
                    end
                    bcLx, bcLy = set_bc_bnds(dir, bcL, HL, BC_TL)

                    laplacian!(dir, opL.LT, opL.CUTT, bcLx, bcLy, geoL.dcap, ny, BC_TL, inside, SOLID,
                                MIXED, b_left[1], b_bottom[1], b_right[1], b_top[1])
                    crank_nicolson!(num, grid, geoL, opL)
                    TL .= reshape(gmres(opL.A,(opL.B*vec(TL) + 2.0*τ*opL.CUTT)), (ny,nx))
>>>>>>> rayleigh_benard
                end
            catch
                @error ("Unphysical temperature field, iteration $current_i")
                break
            end
        end

        if verbose
            if current_i%show_every == 0
                try
                    printstyled(color=:green, @sprintf "\n Current iteration : %d (%d%%) \n" (current_i-1) 100*(current_i-1)/max_iterations)
                    print(@sprintf "V_mean = %.2f  V_max = %.2f  V_min = %.2f\n" mean(V[MIXED]) findmax(V[MIXED])[1] findmin(V[MIXED])[1])
                    print(@sprintf "κ_mean = %.2f  κ_max = %.2f  κ_min = %.2f\n" mean(κ[MIXED]) findmax(κ[MIXED])[1] findmin(κ[MIXED])[1])
                catch
                    @show (MIXED)
                end
            end
        end

        if levelset
            marching_squares!(num, grid, u, periodic_x, periodic_y)

            WAS_LIQUID = copy(LIQUID)
            WAS_SOLID = copy(SOLID)

            MIXED, SOLID, LIQUID = get_cells_indices(iso, inside)
            NB_indices = Flower.get_NB_width(MIXED, NB_indices_base)

            get_iterface_location!(grid, MIXED)

            FRESH_L = intersect(MIXED, WAS_SOLID)
            FRESH_S = intersect(MIXED, WAS_LIQUID)

            init_fresh_cells!(grid, TS, geoS.projection, FRESH_S, periodic_x, periodic_y)
            init_fresh_cells!(grid, TL, geoL.projection, FRESH_L, periodic_x, periodic_y)

            get_curvature(num, grid, u, MIXED, periodic_x, periodic_y)
        end

        current_i -= 1
        κ .= κsave[current_i,:,:]
        u .= usave[current_i,:,:]
    end

    if verbose
        try
            printstyled(color=:blue, @sprintf "\n Final iteration : %d (%d%%) \n" (current_i-1) 100*(current_i-1)/max_iterations)
            print(@sprintf "V_mean = %.2f  V_max = %.2f  V_min = %.2f  V_stdev = %.5f\n" mean(V[MIXED]) findmax(V[MIXED])[1] findmin(V[MIXED])[1] std(V[MIXED]))
            print(@sprintf "κ_mean = %.2f  κ_max = %.2f  κ_min = %.2f  κ_stdev = %.5f\n" mean(κ[MIXED]) findmax(κ[MIXED])[1] findmin(κ[MIXED])[1] std(κ[MIXED]))
            print("\n \n")
        catch
            @show (length(MIXED))
        end
    end
end
