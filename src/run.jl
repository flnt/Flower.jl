using MPI

function run_forward(num, grid, grid_u, grid_v,
    opS, opL, phS, phL, fwd;
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
    BC_pS = Boundaries(
        left = Boundary(),
        right = Boundary(),
        bottom = Boundary(),
        top = Boundary()),
    BC_pL = Boundaries(
        left = Boundary(),
        right = Boundary(),
        bottom = Boundary(),
        top = Boundary()),
    BC_uS = Boundaries(
        left = Boundary(),
        right = Boundary(),
        bottom = Boundary(),
        top = Boundary()),
    BC_uL = Boundaries(
        left = Boundary(),
        right = Boundary(),
        bottom = Boundary(),
        top = Boundary()),
    BC_vS = Boundaries(
        left = Boundary(),
        right = Boundary(),
        bottom = Boundary(),
        top = Boundary()),
    BC_vL = Boundaries(
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
    ns_advection = false,
    navier_stokes = false,
    heat_liquid_phase = false,
    heat_solid_phase = false,
    ns_liquid_phase = false,
    ns_solid_phase = false,
    free_surface = false,
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
    Ra = 0,
    λ = 1,
    )

    if free_surface && stefan
        @error ("Cannot advect the levelset using both free-surface and stefan condition.")
    end

    @unpack L0, A, N, θd, ϵ_κ, ϵ_V, σ, T_inf, τ, L0, NB, Δ, CFL, Re, max_iterations, current_i, save_every, reinit_every, nb_reinit, ϵ, m, θ₀, aniso = num
    @unpack x, y, nx, ny, dx, dy, ind, u, iso, faces, geoS, geoL, V, κ, LSA, LSB = grid
    @unpack Tall, usave, uusave, uvsave, TSsave, TLsave, Tsave, psave, ϕsave, Uxsave, Uysave, Uxcorrsave, Uycorrsave, Vsave, κsave, lengthsave, time, Cd, Cl = fwd

    MPI.Initialized() || MPI.Init()
    PETSc.initialize()

    local ksppS
    local nsS
    local ksppL
    local nsL

    local MIXED; local SOLID; local LIQUID;
    local MIXED_vel_ext; local SOLID_vel_ext; local LIQUID_vel_ext;
    local MIXED_u; local SOLID_u; local LIQUID_u;
    local MIXED_v; local SOLID_v; local LIQUID_v;
    local WAS_MIXED; local WAS_SOLID; local WAS_LIQUID;
    local WAS_MIXED_u; local WAS_SOLID_u; local WAS_LIQUID_u;
    local WAS_MIXED_v; local WAS_SOLID_v; local WAS_LIQUID_v;
    local NB_indices_base; local NB_indices;
    local FRESH_L = []; local FRESH_S = [];
    local FRESH_L_u = []; local FRESH_S_u = [];
    local FRESH_L_v = []; local FRESH_S_v = [];

    local Cum1S = zeros(grid_u.nx*grid_u.ny)
    local Cum1L = zeros(grid_u.nx*grid_u.ny)
    local Cvm1S = zeros(grid_v.nx*grid_v.ny)
    local Cvm1L = zeros(grid_v.nx*grid_v.ny)

    GxpSm1 = copy(opS.Gxp)
    GypSm1 = copy(opS.Gyp)
    GxpLm1 = copy(opL.Gxp)
    GypLm1 = copy(opL.Gyp)
    LpSm1 = copy(opS.Lp)
    LpLm1 = copy(opL.Lp)
    LuSm1 = copy(opS.Lu)
    LvSm1 = copy(opS.Lv)
    LuLm1 = copy(opL.Lu)
    LvLm1 = copy(opL.Lv)
    SCUTpm1 = copy(opS.CUTp)
    LCUTpm1 = copy(opL.CUTp)
    SCUTum1 = copy(opS.CUTu)
    SCUTvm1 = copy(opS.CUTv)
    LCUTum1 = copy(opL.CUTu)
    LCUTvm1 = copy(opL.CUTv)

    HTS = zeros(ny,nx)
    HTL = zeros(ny,nx)
    phS.DT .= phS.T
    phL.DT .= phL.T
    bcTS = copy(phS.DT)
    bcTL = copy(phL.DT)

    bcϕSx = zeros(ny,nx)
    bcϕSy = zeros(ny,nx)
    bcϕLx = zeros(ny,nx)
    bcϕLy = zeros(ny,nx)

    HϕS = zeros(ny,nx)
    HϕL = zeros(ny,nx)
    bcgpSx = zeros(ny,nx)
    bcgpLx = zeros(ny,nx)
    bcgϕSx = zeros(ny,nx)
    bcgϕLx = zeros(ny,nx)
    bcgpSy = zeros(ny,nx)
    bcgpLy = zeros(ny,nx)
    bcgϕSy = zeros(ny,nx)
    bcgϕLy = zeros(ny,nx)

    HuS = zeros(grid_u.ny,grid_u.nx)
    HuL = zeros(grid_u.ny,grid_u.nx)
    bcuSx = zeros(grid_u.ny,grid_u.nx)
    bcuLx = zeros(grid_u.ny,grid_u.nx)
    bcuSy = zeros(grid_u.ny,grid_u.nx)
    bcuLy = zeros(grid_u.ny,grid_u.nx)

    HvS = zeros(grid_v.ny,grid_v.nx)
    HvL = zeros(grid_v.ny,grid_v.nx)
    bcvSx = zeros(grid_v.ny,grid_v.nx)
    bcvLx = zeros(grid_v.ny,grid_v.nx)
    bcvSy = zeros(grid_v.ny,grid_v.nx)
    bcvLy = zeros(grid_v.ny,grid_v.nx)

    Ip = Diagonal(ones(ny*nx))
    Iu = Diagonal(ones(grid_u.ny*grid_u.nx))
    Iv = Diagonal(ones(grid_v.ny*grid_v.nx))
    MpS = copy(Ip)
    MpL = copy(Ip)
    iMpS = copy(Ip)
    iMpL = copy(Ip)
    iMDxS = copy(Ip)
    iMDyS = copy(Ip)
    iMDxL = copy(Ip)
    iMDyL = copy(Ip)
    MuS = copy(Iu)
    MuL = copy(Iu)
    MvS = copy(Iv)
    MvL = copy(Iv)
    iMuS = copy(Iu)
    iMuL = copy(Iu)
    iMvS = copy(Iv)
    iMvL = copy(Iv)
    iMGxS = copy(Iu)
    iMGyS = copy(Iv)
    iMGxL = copy(Iu)
    iMGyL = copy(Iv)

    iMpSm1 = copy(Ip)
    iMpLm1 = copy(Ip)
    MuSm1 = copy(Iu)
    MuLm1 = copy(Iu)
    MvSm1 = copy(Iv)
    MvLm1 = copy(Iv)
    iMuSm1 = copy(Iu)
    iMuLm1 = copy(Iu)
    iMvSm1 = copy(Iv)
    iMvLm1 = copy(Iv)
    iMGxSm1 = copy(Iu)
    iMGySm1 = copy(Iv)
    iMGxLm1 = copy(Iu)
    iMGyLm1 = copy(Iv)

    nullspaceS = true
    if BC_pS.left.t == dir || BC_pS.bottom.t == dir || BC_pS.right.t == dir || BC_pS.top.t == dir || free_surface
        nullspaceS = false
    end
    nullspaceL = true
    if BC_pL.left.t == dir || BC_pL.bottom.t == dir || BC_pL.right.t == dir || BC_pL.top.t == dir || free_surface
        nullspaceL = false
    end

    ns_vecS = ones(ny*nx)
    ns_vecL = ones(ny*nx)

    if periodic_x
        # BC_u.left.ind = ind.periodicL;
        # BC_u.right.ind = ind.periodicR;
        BC_u.left.ind = ind.b_left;
        BC_u.right.ind = ind.b_right;
        BC_u.left.f = BC_u.right.f = periodic
    else
        BC_u.left.ind = ind.b_left;
        BC_u.right.ind = ind.b_right;
    end

    if periodic_y
        # BC_u.bottom.ind = ind.periodicB;
        # BC_u.top.ind = ind.periodicT;
        BC_u.bottom.ind = ind.b_bottom;
        BC_u.top.ind = ind.b_top;
        BC_u.bottom.f = BC_u.top.f = periodic
    else
        BC_u.bottom.ind = ind.b_bottom;
        BC_u.top.ind = ind.b_top;
    end

    usave[1,:,:] .= u
    Tsave[1,:,:] .= phL.T .+ phS.T
    psave[1,:,:] .= phL.p .+ phS.p
    Uxsave[1,:,:] .= phL.u .+ phS.u
    Uysave[1,:,:] .= phL.v .+ phS.v

    if levelset
        marching_squares!(num, grid)
        interpolate_scalar!(grid, grid_u, grid_v, u, grid_u.u, grid_v.u)
        uusave[1,:,:] .= grid_u.u
        uvsave[1,:,:] .= grid_v.u
        marching_squares!(num, grid_u)        
        marching_squares!(num, grid_v)

        NB_indices_base = get_NB_width_indices_base(NB)

        MIXED_vel_ext, SOLID_vel_ext, LIQUID_vel_ext = get_cells_indices(iso, ind.inside)
        MIXED, SOLID, LIQUID = get_cells_indices(iso, ind.all_indices)
        MIXED_u, SOLID_u, LIQUID_u = get_cells_indices(grid_u.iso, grid_u.ind.all_indices)
        MIXED_v, SOLID_v, LIQUID_v = get_cells_indices(grid_v.iso, grid_v.ind.all_indices)

        kill_dead_cells!(phS.T, opS.LT, LIQUID, MIXED, ny)
        kill_dead_cells!(phL.T, opL.LT, SOLID, MIXED, ny)
        kill_dead_cells!(phS.p, opS.Lp, LIQUID, MIXED, ny)
        kill_dead_cells!(phL.p, opL.Lp, SOLID, MIXED, ny)
        kill_dead_cells!(phS.u, opS.Lu, LIQUID_u, MIXED_u, grid_u.ny)
        kill_dead_cells!(phL.u, opL.Lu, SOLID_u, MIXED_u, grid_u.ny)
        kill_dead_cells!(phS.v, opS.Lv, LIQUID_v, MIXED_v, grid_v.ny)
        kill_dead_cells!(phL.v, opL.Lv, SOLID_v, MIXED_v, grid_v.ny)

        # @inbounds @threads for II in MIXED
        #     TS[II] = θd
        #     TL[II] = θd
        # end

        NB_indices = get_NB_width(MIXED, NB_indices_base)
        get_iterface_location!(grid, MIXED)
        get_iterface_location!(grid_u, MIXED_u)
        get_iterface_location!(grid_v, MIXED_v)
        get_interface_location_borders!(grid_u, periodic_x, periodic_y)
        get_interface_location_borders!(grid_v, periodic_x, periodic_y)

        postprocess_grids!(grid, grid_u, grid_v, MIXED, periodic_x, periodic_y, advection, ϵ)

        @inbounds @threads for II in ind.all_indices
            pII = lexicographic(II, ny)
            pJJ = lexicographic(II, grid_v.ny)

            @inbounds iMpSm1[pII,pII] = iMpS[pII,pII]
            @inbounds iMpLm1[pII,pII] = iMpL[pII,pII]
            @inbounds iMGxSm1[pII,pII] = iMGxS[pII,pII]
            @inbounds iMGySm1[pJJ,pJJ] = iMGyS[pJJ,pJJ]
            @inbounds iMGxLm1[pII,pII] = iMGxL[pII,pII]
            @inbounds iMGyLm1[pJJ,pJJ] = iMGyL[pJJ,pJJ]
                
            @inbounds MpS[pII,pII] = geoS.dcap[II,5]
            @inbounds MpL[pII,pII] = geoL.dcap[II,5]
            @inbounds iMpS[pII,pII] = 1. / (geoS.dcap[II,5] + eps(0.01))
            @inbounds iMpL[pII,pII] = 1. / (geoL.dcap[II,5] + eps(0.01))
            @inbounds iMDxS[pII,pII] = 1. / (geoS.dcap[II,5] + eps(0.01))
            @inbounds iMDyS[pII,pII] = 1. / (geoS.dcap[II,5] + eps(0.01))
            @inbounds iMDxL[pII,pII] = 1. / (geoL.dcap[II,5] + eps(0.01))
            @inbounds iMDyL[pII,pII] = 1. / (geoL.dcap[II,5] + eps(0.01))
            @inbounds iMGxS[pII,pII] = 1. / (geoS.dcap[II,8] + eps(0.01))
            @inbounds iMGyS[pJJ,pJJ] = 1. / (geoS.dcap[II,9] + eps(0.01))
            @inbounds iMGxL[pII,pII] = 1. / (geoL.dcap[II,8] + eps(0.01))
            @inbounds iMGyL[pJJ,pJJ] = 1. / (geoL.dcap[II,9] + eps(0.01))
            # @inbounds iMGxS[pII,pII] = 1. / (grid_u.geoS.dcap[II,5] + eps(0.01))
            # @inbounds iMGyS[pJJ,pJJ] = 1. / (grid_v.geoS.dcap[II,5] + eps(0.01))
            # @inbounds iMGxL[pII,pII] = 1. / (grid_u.geoL.dcap[II,5] + eps(0.01))
            # @inbounds iMGyL[pJJ,pJJ] = 1. / (grid_v.geoL.dcap[II,5] + eps(0.01))
        end
        @inbounds @threads for II in grid_u.ind.b_right[1]
            pII = lexicographic(II, grid_u.ny)

            @inbounds iMGxSm1[pII,pII] = iMGxS[pII,pII]
            @inbounds iMGxLm1[pII,pII] = iMGxL[pII,pII]

            @inbounds iMGxS[pII,pII] = 1. / (geoS.dcap[δx⁻(II),10] + eps(0.01))
            @inbounds iMGxL[pII,pII] = 1. / (geoL.dcap[δx⁻(II),10] + eps(0.01))
            # @inbounds iMGxS[pII,pII] = 1. / (grid_u.geoS.dcap[II,5] + eps(0.01))
            # @inbounds iMGxL[pII,pII] = 1. / (grid_u.geoL.dcap[II,5] + eps(0.01))
        end
        @inbounds @threads for II in grid_v.ind.b_top[1]
            pII = lexicographic(II, grid_v.ny)

            @inbounds iMGySm1[pII,pII] = iMGyS[pII,pII]
            @inbounds iMGyLm1[pII,pII] = iMGyL[pII,pII]

            @inbounds iMGyS[pII,pII] = 1. / (geoS.dcap[δy⁻(II),11] + eps(0.01))
            @inbounds iMGyL[pII,pII] = 1. / (geoL.dcap[δy⁻(II),11] + eps(0.01))
            # @inbounds iMGyS[pII,pII] = 1. / (grid_v.geoS.dcap[II,5] + eps(0.01))
            # @inbounds iMGyL[pII,pII] = 1. / (grid_v.geoL.dcap[II,5] + eps(0.01))
        end
        @inbounds @threads for II in grid_u.ind.all_indices
            pII = lexicographic(II, grid_u.ny)

            @inbounds MuSm1[pII,pII] = MuS[pII,pII]
            @inbounds MuLm1[pII,pII] = MuL[pII,pII]
            @inbounds iMuSm1[pII,pII] = 1. / (MuS[pII,pII] + eps(0.01))
            @inbounds iMuLm1[pII,pII] = 1. / (MuL[pII,pII] + eps(0.01))

            @inbounds MuS[pII,pII] = grid_u.geoS.dcap[II,5]
            @inbounds MuL[pII,pII] = grid_u.geoL.dcap[II,5]
            @inbounds iMuS[pII,pII] = 1. / (MuS[pII,pII] + eps(0.01))
            @inbounds iMuL[pII,pII] = 1. / (MuL[pII,pII] + eps(0.01))
        end
        @inbounds @threads for II in grid_v.ind.all_indices
            pII = lexicographic(II, grid_v.ny)

            @inbounds MvSm1[pII,pII] = MvS[pII,pII]
            @inbounds MvLm1[pII,pII] = MvL[pII,pII]
            @inbounds iMvSm1[pII,pII] = 1. / (MvS[pII,pII] + eps(0.01))
            @inbounds iMvLm1[pII,pII] = 1. / (MvL[pII,pII] + eps(0.01))

            @inbounds MvS[pII,pII] = grid_v.geoS.dcap[II,5]
            @inbounds MvL[pII,pII] = grid_v.geoL.dcap[II,5]
            @inbounds iMvS[pII,pII] = 1. / (MvS[pII,pII] + eps(0.01))
            @inbounds iMvL[pII,pII] = 1. / (MvL[pII,pII] + eps(0.01))
        end

        get_curvature(num, grid, MIXED)
        if save_radius
            n_snaps = iszero(max_iterations%save_every) ? max_iterations÷save_every+1 : max_iterations÷save_every+2
            local radius = zeros(n_snaps)
            radius[1] = find_radius(grid, MIXED)
        end
        if hill
            local radius = zeros(max_iterations+1)
            a = zeros(length(MIXED))
            for i in 1:length(MIXED)
                a[i] = geoL.projection[MIXED[i]].pos.y
            end
            radius[1] = mean(a)
        end
    elseif !levelset
        MIXED = [CartesianIndex(-1,-1)]
        MIXED_u = [CartesianIndex(-1,-1)]
        MIXED_v = [CartesianIndex(-1,-1)]
    end

    # Heat
    set_heat!(num, grid, geoS, geoS.projection,
            opS, phS, HTS, bcTS, HuS, HvS,
            BC_TS, BC_uS, BC_vS,
            MIXED, LIQUID, heat_convection)
    
    set_heat!(num, grid, geoL, geoS.projection,
            opL, phL, HTL, bcTL, HuL, HvL,
            BC_TL, BC_uL, BC_vL,
            MIXED, SOLID, heat_convection)

    # Navier-Stokes
    if free_surface
        set_free_surface!(grid, geoS, grid_u, grid_u.geoS, grid_v, grid_v.geoS, opS, phS,
                    HϕS, BC_pS,
                    LpSm1, SCUTpm1, GxpSm1, GypSm1, LIQUID,
                    HuS, bcuSx, bcuSy, BC_uS, LuSm1, SCUTum1, LIQUID_u,
                    HvS, bcvSx, bcvSy, BC_vS, LvSm1, SCUTvm1, LIQUID_v,
                    ns_vecS, MIXED, MIXED_u, MIXED_v,
                    iMuS, iMvS, FRESH_S_u, FRESH_S_v,
                    num.Re, σ, true, ns_advection
        )
    
        set_free_surface!(grid, geoL, grid_u, grid_u.geoL, grid_v, grid_v.geoL, opL, phL,
                    HϕL, BC_pL,
                    LpLm1, LCUTpm1, GxpLm1, GypLm1, SOLID,
                    HuL, bcuLx, bcuLy, BC_uL, LuLm1, LCUTum1, SOLID_u,
                    HvL, bcvLx, bcvLy, BC_vL, LvLm1, LCUTvm1, SOLID_v,
                    ns_vecL, MIXED, MIXED_u, MIXED_v,
                    iMuL, iMvL, FRESH_L_u, FRESH_L_v,
                    num.Re, σ, true, ns_advection
        )
    else
        set_stokes!(grid, geoS, grid_u, grid_u.geoS, grid_v, grid_v.geoS, opS, phS,
                    HϕS, bcgpSx, bcgpSy, bcϕSx, bcϕSy, bcgϕSx, bcgϕSy, BC_pS,
                    LpSm1, SCUTpm1, GxpSm1, GypSm1, LIQUID,
                    HuS, BC_uS, LuSm1, SCUTum1, LIQUID_u,
                    HvS, BC_vS, LvSm1, SCUTvm1, LIQUID_v,
                    ns_vecS, MIXED, MIXED_u, MIXED_v,
                    iMuS, iMvS, FRESH_S_u, FRESH_S_v,
                    num.Re, true, ns_advection, periodic_x, periodic_y
        )
        
        set_stokes!(grid, geoL, grid_u, grid_u.geoL, grid_v, grid_v.geoL, opL, phL,
                    HϕL, bcgpLx, bcgpLy, bcϕLx, bcϕLy, bcgϕLx, bcgϕLy, BC_pL,
                    LpLm1, LCUTpm1, GxpLm1, GypLm1, SOLID,
                    HuL, BC_uL, LuLm1, LCUTum1, SOLID_u,
                    HvL, BC_vL, LvLm1, LCUTvm1, SOLID_v,
                    ns_vecL, MIXED, MIXED_u, MIXED_v,
                    iMuL, iMvL, FRESH_L_u, FRESH_L_v,
                    num.Re, true, ns_advection, periodic_x, periodic_y
        )
    end

    if ns_advection
        Cum1S .= opS.Cu * vec(phS.u) .+ opS.CUTCu
        Cvm1S .= opS.Cv * vec(phS.v) .+ opS.CUTCv
        Cum1L .= opL.Cu * vec(phL.u) .+ opL.CUTCu
        Cvm1L .= opL.Cv * vec(phL.v) .+ opL.CUTCv
    end

    ksppS, nsS = init_ksp_solver(opS.Lp, nx, nullspaceS, ns_vecS)
    ksppL, nsL = init_ksp_solver(opL.Lp, nx, nullspaceL, ns_vecL)

    CFL_sc = τ / Δ^2
    IIOE(LSA, LSB, u, V, ind.inside, CFL_sc, ny)

    if save_length
        lengthsave[1] = arc_length2(geoS.projection, MIXED)
        κsave[1,:,:] .= κ
    end

    current_t = 0.
    while current_i < max_iterations + 1

        if !stefan
            V .= speed*ones(ny, nx)
        end

        if heat
            try
                if heat_solid_phase
                    set_heat!(num, grid, geoS, geoS.projection,
                            opS, phS, HTS, bcTS, HuS, HvS,
                            BC_TS, BC_uS, BC_vS,
                            MIXED, LIQUID, heat_convection)
                    phS.T .= reshape(gmres(opS.A,(opS.B*vec(phS.T) .+ 2.0.*τ.*opS.CUTT .- τ.*opS.CUTCT)), (ny,nx))
                    # TS .= reshape(gmres(AS,(rhs .+ τ*SCUTT)), (n,n))
                end
                if heat_liquid_phase
                    set_heat!(num, grid, geoL, geoS.projection,
                            opL, phL, HTL, bcTL, HuL, HvL,
                            BC_TL, BC_uL, BC_vL,
                            MIXED, SOLID, heat_convection)
                    phL.T .= reshape(gmres(opL.A,(opL.B*vec(phL.T) .+ 2.0.*τ.*opL.CUTT .- τ.*opL.CUTCT)), (ny,nx))
                    # TL .= reshape(gmres(AL,(rhs .+ τ*LCUTT)), (n,n))
                end
            catch
                @error ("Unphysical temperature field, iteration $current_i")
                break
            end
        end

        if stefan
            Stefan_velocity!(num, grid, phS.T, phL.T, MIXED, periodic_x, periodic_y)
            V[MIXED] .*= 1. ./ λ
            if Vmean
                a = mean(V[MIXED])
                V[MIXED] .= a
            end
            velocity_extension!(grid, vcat(SOLID_vel_ext,LIQUID_vel_ext), NB, periodic_x, periodic_y)
        end

        if free_surface
            free_surface_velocity!(grid, grid_u, grid_u.geoL.projection, grid_v, grid_v.geoL.projection,
                                phL.u, phL.v, MIXED, MIXED_u, MIXED_v, periodic_x, periodic_y)

            velocity_extension!(grid, vcat(SOLID_vel_ext,LIQUID_vel_ext), NB, periodic_x, periodic_y)
        end

        if verbose && adaptative_t
            println("τ = $τ")
        end

        if advection
            CFL_sc = τ / Δ^2
            IIOE(grid, LSA, LSB, u, V, CFL_sc, periodic_x, periodic_y)
            try
                u .= reshape(gmres(LSA,(LSB*vec(u))), (ny,nx))
                # u .= sqrt.((x .- current_i*Δ/1).^ 2 + y .^ 2) - (0.5) * ones(nx, ny);
            catch
                @error ("Inadequate level set function, iteration $current_i")
                break
            end
            if analytical
                u[ind.b_top[1]] .= sqrt.(x[ind.b_top[1]] .^ 2 + y[ind.b_top[1]] .^ 2) .- (num.R + speed*current_i*τ);
                u[ind.b_bottom[1]] .= sqrt.(x[ind.b_bottom[1]] .^ 2 + y[ind.b_bottom[1]] .^ 2) .- (num.R + speed*current_i*τ);
                u[ind.b_left[1]] .= sqrt.(x[ind.b_left[1]] .^ 2 + y[ind.b_left[1]] .^ 2) .- (num.R + speed*current_i*τ);
                u[ind.b_right[1]] .= sqrt.(x[ind.b_right[1]] .^ 2 + y[ind.b_right[1]] .^ 2) .- (num.R + speed*current_i*τ);
            elseif nb_reinit > 0
                if current_i%num.reinit_every == 0
                    FE_reinit(grid, ind, u, nb_reinit, BC_u, periodic_x, periodic_y)
                end
            end
        end

        if verbose
            if current_i%show_every == 0
                try
                    printstyled(color=:green, @sprintf "\n Current iteration : %d (%d%%) \n" (current_i-1) 100*(current_i-1)/max_iterations)
                    if heat || free_surface
                        print(@sprintf "V_mean = %.2f  V_max = %.2f  V_min = %.2f\n" mean(V[MIXED]) findmax(V[MIXED])[1] findmin(V[MIXED])[1])
                        print(@sprintf "κ_mean = %.2f  κ_max = %.2f  κ_min = %.2f\n" mean(κ[MIXED]) findmax(κ[MIXED])[1] findmin(κ[MIXED])[1])
                    end
                    if navier_stokes
                        if ns_solid_phase
                            normuS = norm(phS.u)
                            normvS = norm(phS.v)
                            normpS = norm(phS.p.*τ)
                            print("$(@sprintf("norm(uS) %.6e", normuS))\t$(@sprintf("norm(vS) %.6e", normvS))\t$(@sprintf("norm(pS) %.6e", normpS))\n")
                        end
                        if ns_liquid_phase
                            normuL = norm(phL.u)
                            normvL = norm(phL.v)
                            normpL = norm(phL.p.*τ)
                            print("$(@sprintf("norm(uL) %.6e", normuL))\t$(@sprintf("norm(vL) %.6e", normvL))\t$(@sprintf("norm(pL) %.6e", normpL))\n")
                        end
                    end
                catch
                    @show (MIXED)
                end
            end
        end


        if levelset && (advection || current_i<2)
            grid.α .= NaN
            grid_u.α .= NaN
            grid_v.α .= NaN
            faces .= 0.
            grid_u.faces .= 0.
            grid_v.faces .= 0.
            marching_squares!(num, grid)
            interpolate_scalar!(grid, grid_u, grid_v, u, grid_u.u, grid_v.u)
            marching_squares!(num, grid_u)
            marching_squares!(num, grid_v)

            WAS_MIXED = copy(MIXED)
            WAS_LIQUID = copy(LIQUID)
            WAS_SOLID = copy(SOLID)
            WAS_MIXED_u = copy(MIXED_u)
            WAS_LIQUID_u = copy(LIQUID_u)
            WAS_SOLID_u = copy(SOLID_u)
            WAS_MIXED_v = copy(MIXED_v)
            WAS_LIQUID_v = copy(LIQUID_v)
            WAS_SOLID_v = copy(SOLID_v)

            MIXED_vel_ext, SOLID_vel_ext, LIQUID_vel_ext = get_cells_indices(iso, ind.inside)
            MIXED, SOLID, LIQUID = get_cells_indices(iso, ind.all_indices)
            MIXED_u, SOLID_u, LIQUID_u = get_cells_indices(grid_u.iso, grid_u.ind.all_indices)
            MIXED_v, SOLID_v, LIQUID_v = get_cells_indices(grid_v.iso, grid_v.ind.all_indices)

            kill_dead_cells!(phS.T, opS.LT, LIQUID, MIXED, ny)
            kill_dead_cells!(phL.T, opL.LT, SOLID, MIXED, ny)

            # @inbounds @threads for II in MIXED
            #     TS[II] = θd
            #     TL[II] = θd
            # end

            NB_indices = get_NB_width(MIXED, NB_indices_base)

            get_iterface_location!(grid, MIXED)
            get_iterface_location!(grid_u, MIXED_u)
            get_iterface_location!(grid_v, MIXED_v)
            get_interface_location_borders!(grid_u, periodic_x, periodic_y)
            get_interface_location_borders!(grid_v, periodic_x, periodic_y)

            tmpFRESH_L = intersect(findall(geoL.emptied), WAS_MIXED)
            tmpFRESH_S = intersect(findall(geoS.emptied), WAS_MIXED)
            tmpFRESH_L_u = intersect(findall(grid_u.geoL.emptied), WAS_MIXED_u)
            tmpFRESH_S_u = intersect(findall(grid_u.geoS.emptied), WAS_MIXED_u)
            tmpFRESH_L_v = intersect(findall(grid_v.geoL.emptied), WAS_MIXED_v)
            tmpFRESH_S_v = intersect(findall(grid_v.geoS.emptied), WAS_MIXED_v)

            geoL.emptied .= false
            geoS.emptied .= false
            grid_u.geoL.emptied .= false
            grid_u.geoS.emptied .= false
            grid_v.geoL.emptied .= false
            grid_v.geoS.emptied .= false

            postprocess_grids!(grid, grid_u, grid_v, MIXED, periodic_x, periodic_y, advection, ϵ)

            @inbounds @threads for II in ind.all_indices
                pII = lexicographic(II, ny)
                pJJ = lexicographic(II, grid_v.ny)
    
                @inbounds iMpSm1[pII,pII] = iMpS[pII,pII]
                @inbounds iMpLm1[pII,pII] = iMpL[pII,pII]
                @inbounds iMGxSm1[pII,pII] = iMGxS[pII,pII]
                @inbounds iMGySm1[pJJ,pJJ] = iMGyS[pJJ,pJJ]
                @inbounds iMGxLm1[pII,pII] = iMGxL[pII,pII]
                @inbounds iMGyLm1[pJJ,pJJ] = iMGyL[pJJ,pJJ]
                
                @inbounds MpS[pII,pII] = geoS.dcap[II,5]
                @inbounds MpL[pII,pII] = geoL.dcap[II,5]
                @inbounds iMpS[pII,pII] = 1. / (geoS.dcap[II,5] + eps(0.01))
                @inbounds iMpL[pII,pII] = 1. / (geoL.dcap[II,5] + eps(0.01))
                @inbounds iMDxS[pII,pII] = 1. / (geoS.dcap[II,5] + eps(0.01))
                @inbounds iMDyS[pII,pII] = 1. / (geoS.dcap[II,5] + eps(0.01))
                @inbounds iMDxL[pII,pII] = 1. / (geoL.dcap[II,5] + eps(0.01))
                @inbounds iMDyL[pII,pII] = 1. / (geoL.dcap[II,5] + eps(0.01))
                @inbounds iMGxS[pII,pII] = 1. / (geoS.dcap[II,8] + eps(0.01))
                @inbounds iMGyS[pJJ,pJJ] = 1. / (geoS.dcap[II,9] + eps(0.01))
                @inbounds iMGxL[pII,pII] = 1. / (geoL.dcap[II,8] + eps(0.01))
                @inbounds iMGyL[pJJ,pJJ] = 1. / (geoL.dcap[II,9] + eps(0.01))
                # @inbounds iMGxS[pII,pII] = 1. / (grid_u.geoS.dcap[II,5] + eps(0.01))
                # @inbounds iMGyS[pJJ,pJJ] = 1. / (grid_v.geoS.dcap[II,5] + eps(0.01))
                # @inbounds iMGxL[pII,pII] = 1. / (grid_u.geoL.dcap[II,5] + eps(0.01))
                # @inbounds iMGyL[pJJ,pJJ] = 1. / (grid_v.geoL.dcap[II,5] + eps(0.01))
            end
            @inbounds @threads for II in grid_u.ind.b_right[1]
                pII = lexicographic(II, grid_u.ny)

                @inbounds iMGxSm1[pII,pII] = iMGxS[pII,pII]
                @inbounds iMGxLm1[pII,pII] = iMGxL[pII,pII]

                @inbounds iMGxS[pII,pII] = 1. / (geoS.dcap[δx⁻(II),10] + eps(0.01))
                @inbounds iMGxL[pII,pII] = 1. / (geoL.dcap[δx⁻(II),10] + eps(0.01))
                # @inbounds iMGxS[pII,pII] = 1. / (grid_u.geoS.dcap[II,5] + eps(0.01))
                # @inbounds iMGxL[pII,pII] = 1. / (grid_u.geoL.dcap[II,5] + eps(0.01))
            end
            @inbounds @threads for II in grid_v.ind.b_top[1]
                pII = lexicographic(II, grid_v.ny)

                @inbounds iMGySm1[pII,pII] = iMGyS[pII,pII]
                @inbounds iMGyLm1[pII,pII] = iMGyL[pII,pII]

                @inbounds iMGyS[pII,pII] = 1. / (geoS.dcap[δy⁻(II),11] + eps(0.01))
                @inbounds iMGyL[pII,pII] = 1. / (geoL.dcap[δy⁻(II),11] + eps(0.01))
                # @inbounds iMGyS[pII,pII] = 1. / (grid_v.geoS.dcap[II,5] + eps(0.01))
                # @inbounds iMGyL[pII,pII] = 1. / (grid_v.geoL.dcap[II,5] + eps(0.01))
            end
            @inbounds @threads for II in grid_u.ind.all_indices
                pII = lexicographic(II, grid_u.ny)
    
                @inbounds MuSm1[pII,pII] = MuS[pII,pII]
                @inbounds MuLm1[pII,pII] = MuL[pII,pII]
                @inbounds iMuSm1[pII,pII] = 1. / (MuS[pII,pII] + eps(0.01))
                @inbounds iMuLm1[pII,pII] = 1. / (MuL[pII,pII] + eps(0.01))
    
                @inbounds MuS[pII,pII] = grid_u.geoS.dcap[II,5]
                @inbounds MuL[pII,pII] = grid_u.geoL.dcap[II,5]
                @inbounds iMuS[pII,pII] = 1. / (MuS[pII,pII] + eps(0.01))
                @inbounds iMuL[pII,pII] = 1. / (MuL[pII,pII] + eps(0.01))
            end
            @inbounds @threads for II in grid_v.ind.all_indices
                pII = lexicographic(II, grid_v.ny)
    
                @inbounds MvSm1[pII,pII] = MvS[pII,pII]
                @inbounds MvLm1[pII,pII] = MvL[pII,pII]
                @inbounds iMvSm1[pII,pII] = 1. / (MvS[pII,pII] + eps(0.01))
                @inbounds iMvLm1[pII,pII] = 1. / (MvL[pII,pII] + eps(0.01))
    
                @inbounds MvS[pII,pII] = grid_v.geoS.dcap[II,5]
                @inbounds MvL[pII,pII] = grid_v.geoL.dcap[II,5]
                @inbounds iMvS[pII,pII] = 1. / (MvS[pII,pII] + eps(0.01))
                @inbounds iMvL[pII,pII] = 1. / (MvL[pII,pII] + eps(0.01))
            end

            get_curvature(num, grid, MIXED)

            FRESH_L = union(tmpFRESH_L, intersect(MIXED, WAS_SOLID))
            FRESH_S = union(tmpFRESH_S, intersect(MIXED, WAS_LIQUID))
            FRESH_L_u = union(tmpFRESH_L_u, intersect(MIXED_u, WAS_SOLID_u))
            FRESH_S_u = union(tmpFRESH_S_u, intersect(MIXED_u, WAS_LIQUID_u))
            FRESH_L_v = union(tmpFRESH_L_v, intersect(MIXED_v, WAS_SOLID_v))
            FRESH_S_v = union(tmpFRESH_S_v, intersect(MIXED_v, WAS_LIQUID_v))

            init_fresh_cells!(grid, phS.T, geoS.projection, FRESH_S, periodic_x, periodic_y)
            init_fresh_cells!(grid, phL.T, geoL.projection, FRESH_L, periodic_x, periodic_y)

            if iszero(current_i%save_every) || current_i==max_iterations
                snap = current_i÷save_every+1
                if save_radius
                    radius[snap] = find_radius(grid, MIXED)
                end
                if hill
                    a = zeros(length(MIXED))
                    for i in 1:length(MIXED)
                        a[i] = geoL.projection[MIXED[i]].pos.y
                    end
                    radius[snap] = mean(a)
                end
                if save_length
                    lengthsave[snap] = arc_length2(geoS.projection, MIXED)
                    κsave[snap,:,:] .= κ
                end
            end
        end

        if navier_stokes
            if !free_surface
                no_slip_condition!(grid, grid_u, grid_v)
                grid_u.V .= imfilter(grid_u.V, Kernel.gaussian(2))
                grid_v.V .= imfilter(grid_v.V, Kernel.gaussian(2))
                # grid_u.V .= Δ / (1 * τ)
                # grid_v.V .= 0.0
            end

            if ns_solid_phase
                if free_surface
                    set_free_surface!(grid, geoS, grid_u, grid_u.geoS, grid_v, grid_v.geoS, opS, phS,
                                HϕS, BC_pS,
                                LpSm1, SCUTpm1, GxpSm1, GypSm1, LIQUID,
                                HuS, bcuSx, bcuSy, BC_uS, LuSm1, SCUTum1, LIQUID_u,
                                HvS, bcvSx, bcvSy, BC_vS, LvSm1, SCUTvm1, LIQUID_v,
                                ns_vecS, MIXED, MIXED_u, MIXED_v,
                                iMuS, iMvS, FRESH_S_u, FRESH_S_v,
                                num.Re, σ, true, ns_advection
                    )
                else
                    set_stokes!(grid, geoS, grid_u, grid_u.geoS, grid_v, grid_v.geoS, opS, phS,
                                HϕS, bcgpSx, bcgpSy, bcϕSx, bcϕSy, bcgϕSx, bcgϕSy, BC_pS,
                                LpSm1, SCUTpm1, GxpSm1, GypSm1, LIQUID,
                                HuS, BC_uS, LuSm1, SCUTum1, LIQUID_u,
                                HvS, BC_vS, LvSm1, SCUTvm1, LIQUID_v,
                                ns_vecS, MIXED, MIXED_u, MIXED_v,
                                iMuS, iMvS, FRESH_S_u, FRESH_S_v,
                                num.Re, advection, ns_advection, periodic_x, periodic_y)
                end
                pressure_projection!(num, grid, geoS, grid_u, grid_u.geoS, grid_v, grid_v.geoS, opS, phS,
                                LuSm1, LvSm1, SCUTum1, SCUTvm1, Cum1S, Cvm1S,
                                ksppS, nsS, ns_vecS, GxpSm1, GypSm1,
                                MpS, iMpS, MuS, MvS, iMGxS, iMGyS, iMDxS, iMDyS,
                                iMuSm1, iMvSm1, iMGxSm1, iMGySm1, MuSm1, MvSm1,
                                MIXED, MIXED_u, MIXED_v, SOLID, LIQUID, LIQUID_u, LIQUID_v,
                                FRESH_S, FRESH_S_u, FRESH_S_v, nullspaceS, ns_advection, Ra,
                                periodic_x, periodic_y)
            end
            if ns_liquid_phase
                if free_surface
                    set_free_surface!(grid, geoL, grid_u, grid_u.geoL, grid_v, grid_v.geoL, opL, phL,
                                HϕL, BC_pL,
                                LpLm1, LCUTpm1, GxpLm1, GypLm1, SOLID,
                                HuL, bcuLx, bcuLy, BC_uL, LuLm1, LCUTum1, SOLID_u,
                                HvL, bcvLx, bcvLy, BC_vL, LvLm1, LCUTvm1, SOLID_v,
                                ns_vecL, MIXED, MIXED_u, MIXED_v,
                                iMuL, iMvL, FRESH_L_u, FRESH_L_v,
                                num.Re, σ, true, ns_advection)
                else
                    set_stokes!(grid, geoL, grid_u, grid_u.geoL, grid_v, grid_v.geoL, opL, phL,
                                HϕL, bcgpLx, bcgpLy, bcϕLx, bcϕLy, bcgϕLx, bcgϕLy, BC_pL,
                                LpLm1, LCUTpm1, GxpLm1, GypLm1, SOLID,
                                HuL, BC_uL, LuLm1, LCUTum1, SOLID_u,
                                HvL, BC_vL, LvLm1, LCUTvm1, SOLID_v,
                                ns_vecL, MIXED, MIXED_u, MIXED_v,
                                iMuL, iMvL, FRESH_L_u, FRESH_L_v,
                                num.Re, advection, ns_advection, periodic_x, periodic_y)
                end

                pressure_projection!(num, grid, geoL, grid_u, grid_u.geoL, grid_v, grid_v.geoL, opL, phL,
                            LuLm1, LvLm1, LCUTum1, LCUTvm1, Cum1L, Cvm1L,
                            ksppL, nsL, ns_vecL, GxpLm1, GypLm1,
                            MpL, iMpL, MuL, MvL, iMGxL, iMGyL, iMDxL, iMDyL,
                            iMuLm1, iMvLm1, iMGxLm1, iMGyLm1, MuLm1, MvLm1,
                            MIXED, MIXED_u, MIXED_v, LIQUID, SOLID, SOLID_u, SOLID_v,
                            FRESH_L, FRESH_L_u, FRESH_L_v, nullspaceL, ns_advection, Ra,
                            periodic_x, periodic_y)
            end
        end

        current_t += τ
        if iszero(current_i%save_every) || current_i==max_iterations
            snap = current_i÷save_every+1
            if current_i==max_iterations
                snap = size(Tsave,1)
            end
            time[snap] = current_t
            Vsave[snap,:,:] .= V
            usave[snap,:,:] .= u
            uusave[snap,:,:] .= grid_u.u
            uvsave[snap,:,:] .= grid_v.u
            Tsave[snap,:,:] .= phL.T .+ phS.T
            if ns_solid_phase && ns_liquid_phase
                psave[snap,:,:] .= phL.p[:,:].*geoL.cap[:,:,5] .+ phS.p[:,:].*geoS.cap[:,:,5]
                ϕsave[snap,:,:] .= phL.ϕ[:,:].*geoL.cap[:,:,5] .+ phS.ϕ[:,:].*geoS.cap[:,:,5]
                Uxsave[snap,:,:] .= phL.u[:,:].*grid_u.geoL.cap[:,:,5] .+ phS.u[:,:].*geoS.capu[:,:,5]
                Uysave[snap,:,:] .= phL.v[:,:].*grid_v.geoL.cap[:,:,5] .+ phS.v[:,:].*geoS.capv[:,:,5]
                Uxcorrsave[snap,:,:] .= phL.ucorr[:,:].*grid_u.geoL.cap[:,:,5] .+ phS.ucorr[:,:].*geoS.capu[:,:,5]
                Uycorrsave[snap,:,:] .= phL.vcorr[:,:].*grid_v.geoL.cap[:,:,5] .+ phS.vcorr[:,:].*geoS.capv[:,:,5]
                force_coefficients!(num, grid, grid_u, grid_v, opL, fwd; step=snap)
            elseif ns_solid_phase
                psave[snap,:,:] .= phS.p[:,:]
                ϕsave[snap,:,:] .= phS.ϕ[:,:]
                Uxsave[snap,:,:] .= phS.u[:,:]
                Uysave[snap,:,:] .= phS.v[:,:]
                Uxcorrsave[snap,:,:] .= phS.ucorr[:,:]
                Uycorrsave[snap,:,:] .= phS.vcorr[:,:]
            elseif ns_liquid_phase
                psave[snap,:,:] .= phL.p[:,:]
                ϕsave[snap,:,:] .= phL.ϕ[:,:]
                Uxsave[snap,:,:] .= phL.u[:,:]
                Uysave[snap,:,:] .= phL.v[:,:]
                Uxcorrsave[snap,:,:] .= phL.ucorr[:,:]
                Uycorrsave[snap,:,:] .= phL.vcorr[:,:]
                force_coefficients!(num, grid, grid_u, grid_v, opL, fwd; step=snap)
            end
        end

        current_i += 1

        if adaptative_t
            τ = min(CFL*Δ^2*Re, CFL*Δ/max(abs.(V)..., abs.(phL.u)..., abs.(phL.v)..., abs.(phS.u)..., abs.(phS.v)...))
        end
    end

    if verbose
        try
            printstyled(color=:blue, @sprintf "\n Final iteration : %d (%d%%) \n" (current_i-1) 100*(current_i-1)/max_iterations)
            if heat || free_surface
                print(@sprintf "V_mean = %.2f  V_max = %.2f  V_min = %.2f  V_stdev = %.5f\n" mean(V[MIXED]) findmax(V[MIXED])[1] findmin(V[MIXED])[1] std(V[MIXED]))
                print(@sprintf "κ_mean = %.2f  κ_max = %.2f  κ_min = %.2f  κ_stdev = %.5f\n" mean(κ[MIXED]) findmax(κ[MIXED])[1] findmin(κ[MIXED])[1] std(κ[MIXED]))
            end
            if navier_stokes
                if ns_solid_phase
                    normuS = norm(phS.u)
                    normvS = norm(phS.v)
                    normpS = norm(phS.p.*τ)
                    print("$(@sprintf("norm(uS) %.6e", normuS))\t$(@sprintf("norm(vS) %.6e", normvS))\t$(@sprintf("norm(pS) %.6e", normpS))\n")
                end
                if ns_liquid_phase
                    normuL = norm(phL.u)
                    normvL = norm(phL.v)
                    normpL = norm(phL.p.*τ)
                    print("$(@sprintf("norm(uL) %.6e", normuL))\t$(@sprintf("norm(vL) %.6e", normvL))\t$(@sprintf("norm(pL) %.6e", normpL))\n")
                end
            end
            print("\n\n")
        catch
            @show (length(MIXED))
        end
    end

    if nullspaceS
        PETSc.destroy(nsS)
    end
    PETSc.destroy(ksppS)

    if nullspaceL
        PETSc.destroy(nsL)
    end
    PETSc.destroy(ksppL)

    if levelset
        if save_radius || hill
            return MIXED, SOLID, LIQUID, radius
        end
        return MIXED, SOLID, LIQUID
    else
        return MIXED
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
        marching_squares!(num, grid)

        bcs!(faces, BC_u.left, dx[1,1])
        bcs!(faces, BC_u.right, dx[1,end])
        bcs!(faces, BC_u.bottom, dy[1,1])
        bcs!(faces, BC_u.top, dy[end,1])

        NB_indices_base = get_NB_width_indices_base(NB)

        MIXED, SOLID, LIQUID = get_cells_indices(iso, inside)
        NB_indices = get_NB_width(MIXED, NB_indices_base)

        get_iterface_location!(grid, MIXED)
        get_curvature(num, grid, MIXED)
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
            marching_squares!(num, grid)

            bcs!(faces, BC_u.left, dx[1,1])
            bcs!(faces, BC_u.right, dx[1,end])
            bcs!(faces, BC_u.bottom, dy[1,1])
            bcs!(faces, BC_u.top, dy[end,1])

            WAS_LIQUID = copy(LIQUID)
            WAS_SOLID = copy(SOLID)

            MIXED, SOLID, LIQUID = get_cells_indices(iso, inside)
            NB_indices = Flower.get_NB_width(MIXED, NB_indices_base)

            get_iterface_location!(grid, MIXED)

            FRESH_L = intersect(MIXED, WAS_SOLID)
            FRESH_S = intersect(MIXED, WAS_LIQUID)

            init_fresh_cells!(grid, TS, geoS.projection, FRESH_S)
            init_fresh_cells!(grid, TL, geoL.projection, FRESH_L)

            get_curvature(num, grid, MIXED)
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
