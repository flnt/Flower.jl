using MPI

function run_forward(num, idx, idxu, idxv, tmp, fwd;
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
        left = Boundary(t = dir, f = dirichlet),
        right = Boundary(t = dir, f = dirichlet),
        bottom = Boundary(t = dir, f = dirichlet),
        top = Boundary(t = dir, f = dirichlet)),
    BC_uL = Boundaries(
        left = Boundary(t = dir, f = dirichlet),
        right = Boundary(t = dir, f = dirichlet),
        bottom = Boundary(t = dir, f = dirichlet),
        top = Boundary(t = dir, f = dirichlet)),
    BC_vS = Boundaries(
        left = Boundary(t = dir, f = dirichlet),
        right = Boundary(t = dir, f = dirichlet),
        bottom = Boundary(t = dir, f = dirichlet),
        top = Boundary(t = dir, f = dirichlet)),
    BC_vL = Boundaries(
        left = Boundary(t = dir, f = dirichlet),
        right = Boundary(t = dir, f = dirichlet),
        bottom = Boundary(t = dir, f = dirichlet),
        top = Boundary(t = dir, f = dirichlet)),
    BC_u = Boundaries(
        left = Boundary(),
        right = Boundary(),
        bottom = Boundary(),
        top = Boundary()),
    stefan = false,
    advection = false,
    heat = false,
    navier_stokes = false,
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

    @unpack L0, A, N, θd, ϵ_κ, ϵ_V, T_inf, τ, L0, NB, n, Δ, CFL, Re, max_iterations, current_i, reinit_every, nb_reinit, ϵ, X, Y, Xu, Yu, Xv, Yv, H, B, BT, m, θ₀, aniso = num
    @unpack all_indices, inside, b_left, b_bottom, b_right, b_top = idxu
    all_indices_u, inside_u, b_left_u, b_bottom_u, b_right_u, b_top_u = all_indices, inside, b_left, b_bottom, b_right, b_top
    @unpack all_indices, inside, b_left, b_bottom, b_right, b_top = idxv
    all_indices_v, inside_v, b_left_v, b_bottom_v, b_right_v, b_top_v = all_indices, inside, b_left, b_bottom, b_right, b_top
    @unpack all_indices, inside, b_left, b_bottom, b_right, b_top = idx
    @unpack SCUTT, LCUTT, SCUTp, LCUTp, SCUTu, LCUTu, SCUTv, LCUTv, SCUTDx, SCUTDy, LCUTDx, LCUTDy, LTS, LTL, LpS, LpL, LuS, LuL, LvS, LvL, AS, AL, BS, BL, LSA, LSB, GxpS, GxpL, GypS, GypL, DxuS, DxuL, DyvS, DyvL, ApS, ApL, AuS, AuL, AvS, AvL, SOL, LIQ, SOLu, LIQu, SOLv, LIQv, sol_projection, liq_projection, sol_projectionu, liq_projectionu, sol_projectionv, liq_projectionv, sol_centroid, liq_centroid, mid_point, sol_centroidu, liq_centroidu, mid_pointu, sol_centroidv, liq_centroidv, mid_pointv = tmp
    @unpack iso, isou, isov, u, uu, uv, TS, TL, pS, pL, uS, uL, vS, vL, Tall, DTS, DTL, V, κ, κu, κv, usave, TSsave, TLsave, Tsave, psave, Uxsave, Uysave, Vsave, κsave, lengthsave = fwd

    MPI.Initialized() || MPI.Init()
    PETSc.initialize()

    local ksppS
    local kspuS
    local kspvS
    local nsS
    local ksppL
    local kspuL
    local kspvL
    local nsL

    local MIXED; local SOLID; local LIQUID;
    local MIXED_u; local SOLID_u; local LIQUID_u;
    local MIXED_v; local SOLID_v; local LIQUID_v;
    local WAS_SOLID; local WAS_LIQUID;
    local WAS_SOLID_u; local WAS_LIQUID_u;
    local WAS_SOLID_v; local WAS_LIQUID_v;
    local NB_indices_base; local NB_indices;
    local FRESH_L; local FRESH_S;
    local FRESH_L_u; local FRESH_S_u;
    local FRESH_L_v; local FRESH_S_v;

    local faces = zeros(n,n,4);
    local facesu = zeros(n,n+1,4);
    local facesv = zeros(n+1,n,4);

    HTS = zeros(n,n)
    HTL = zeros(n,n)
    DTS .= TS
    DTL .= TL
    bcTS = copy(DTS)
    bcTL = copy(DTL)

    bcpSx = zeros(n,n)
    bcpSy = zeros(n,n)
    bcpLx = zeros(n,n)
    bcpLy = zeros(n,n)

    HuS = zeros(n,n+1)
    HuL = zeros(n,n+1)
    DuS = zeros(n,n+1)
    DuL = zeros(n,n+1)
    bcuS = zeros(n,n+1)
    bcuL = zeros(n,n+1)

    HvS = zeros(n+1,n)
    HvL = zeros(n+1,n)
    DvS = zeros(n+1,n)
    DvL = zeros(n+1,n)
    bcvS = zeros(n+1,n)
    bcvL = zeros(n+1,n)

    Ip = spdiagm(0 => ones(n^2))
    Iuv = spdiagm(0 => ones(n*(n+1)))
    MpS = copy(Ip)
    MpL = copy(Ip)
    iMpS = copy(Ip)
    iMpL = copy(Ip)
    iMDxS = copy(Ip)
    iMDyS = copy(Ip)
    iMDxL = copy(Ip)
    iMDyL = copy(Ip)
    MuS = copy(Iuv)
    MuL = copy(Iuv)
    MvS = copy(Iuv)
    MvL = copy(Iuv)
    iMGxS = copy(Iuv)
    iMGyS = copy(Iuv)
    iMGxL = copy(Iuv)
    iMGyL = copy(Iuv)

    ns_vecS = ones(n^2)
    ns_vecL = ones(n^2)

    iRe = 1 / Re

    if periodic_x
        BC_u.left.ind = idx.periodicL;
        BC_u.right.ind = idx.periodicR;
        BC_u.left.f = BC_u.right.f = periodic
    else
        BC_u.left.ind =b_left;
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

    usave[1, :, :] .= u[:,:]
    Tsave[1, :, :] .= TL[:,:] .+ TS[:,:]

    if levelset
        marching_squares!(H, iso, u, TS, TL, κ, SOL, LIQ, sol_projection, liq_projection, Δ, L0, B, BT, inside, ϵ, n, n, faces)
        interpolate_LS!(u, uu, uv, B, BT, n, inside,b_left[1], b_bottom[1], b_right[1], b_top[1])
        marching_squares!(H, isou, uu, uS, uL, κu, SOLu, LIQu, sol_projectionu, liq_projectionu, Δ, L0, B, BT, inside_u, ϵ, n+1, n, facesu)
        marching_squares!(H, isov, uv, vS, vL, κv, SOLv, LIQv, sol_projectionv, liq_projectionv, Δ, L0, B, BT, inside_v, ϵ, n, n+1, facesv)

        bcs!(faces, BC_u.left, Δ)
        bcs!(faces, BC_u.right, Δ)
        bcs!(faces, BC_u.bottom, Δ)
        bcs!(faces, BC_u.top, Δ)

        bcs!(facesu, BC_u.left, Δ)
        bcs!(facesu, BC_u.right, Δ)
        bcs!(facesu, BC_u.bottom, Δ)
        bcs!(facesu, BC_u.top, Δ)

        bcs!(facesv, BC_u.left, Δ)
        bcs!(facesv, BC_u.right, Δ)
        bcs!(facesv, BC_u.bottom, Δ)
        bcs!(facesv, BC_u.top, Δ)

        NB_indices_base = get_NB_width_indices_base(NB)

        MIXED, SOLID, LIQUID = get_cells_indices(iso, inside)
        MIXED_u, SOLID_u, LIQUID_u = get_cells_indices(isou, inside_u)
        MIXED_v, SOLID_v, LIQUID_v = get_cells_indices(isov, inside_v)
        
        kill_dead_cells!(TS, LIQUID)
        kill_dead_cells!(TL, SOLID)
        kill_dead_cells!(pS, LIQUID)
        kill_dead_cells!(pL, SOLID)
        kill_dead_cells!(uS, LIQUID_u)
        kill_dead_cells!(uL, SOLID_u)
        kill_dead_cells!(vS, LIQUID_v)
        kill_dead_cells!(vL, SOLID_v)

        @inbounds @threads for II in MIXED
            TS[II] = θd
            TL[II] = θd
        end
        @inbounds @threads for II in LIQUID
            pII = lexicographic(II, n)
            @inbounds ns_vecS[pII] = 0.
        end
        @inbounds @threads for II in SOLID
            pII = lexicographic(II, n)
            @inbounds ns_vecL[pII] = 0.
        end

        NB_indices = get_NB_width(MIXED, NB_indices_base)

        get_iterface_location!(gcc, X, Y, iso, u, TS, TL, κ, SOL, LIQ, sol_projection, liq_projection, sol_centroid, liq_centroid, mid_point, Δ, L0, B, BT, idx, MIXED, ϵ, n, faces, periodic_x, periodic_y)
        get_iterface_location!(gfcx, Xu, Yu, isou, uu, uS, uL, κu, SOLu, LIQu, sol_projectionu, liq_projectionu, sol_centroidu, liq_centroidu, mid_pointu, Δ, L0, B, BT, idxu, MIXED_u, ϵ, n, facesu, periodic_x, periodic_y)
        get_iterface_location!(gfcy, Xv, Yv, isov, uv, vS, vL, κv, SOLv, LIQv, sol_projectionv, liq_projectionv, sol_centroidv, liq_centroidv, mid_pointv, Δ, L0, B, BT, idxv, MIXED_v, ϵ, n, facesv, periodic_x, periodic_y)
        
        @inbounds @threads for II in all_indices
            pII = lexicographic(II, n)
            pJJ = lexicographic(II, n+1)
            @inbounds MpS[pII,pII] = SOL[II,5] * Δ^2
            @inbounds MpL[pII,pII] = LIQ[II,5] * Δ^2
            @inbounds iMpS[pII,pII] = 1. / (SOL[II,5] * Δ^2 + eps(0.01))
            @inbounds iMpL[pII,pII] = 1. / (LIQ[II,5] * Δ^2 + eps(0.01))
            @inbounds iMDxS[pII,pII] = 1. / (SOLu[II,10] * Δ^2 + eps(0.01))
            @inbounds iMDyS[pII,pII] = 1. / (SOLv[II,11] * Δ^2 + eps(0.01))
            @inbounds iMDxL[pII,pII] = 1. / (LIQu[II,10] * Δ^2 + eps(0.01))
            @inbounds iMDyL[pII,pII] = 1. / (LIQv[II,11] * Δ^2 + eps(0.01))
            @inbounds iMGxS[pII,pII] = 1. / (SOL[II,8] * Δ^2 + eps(0.01))
            @inbounds iMGyS[pJJ,pJJ] = 1. / (SOL[II,9] * Δ^2 + eps(0.01))
            @inbounds iMGxL[pII,pII] = 1. / (LIQ[II,8] * Δ^2 + eps(0.01))
            @inbounds iMGyL[pJJ,pJJ] = 1. / (LIQ[II,9] * Δ^2 + eps(0.01))
        end
        @inbounds @threads for II in b_right_u[1]
            pII = lexicographic(II, n)
            @inbounds iMGxS[pII,pII] = 1. / (SOL[δx⁻(II),10] * Δ^2 + eps(0.01))
            @inbounds iMGxL[pII,pII] = 1. / (LIQ[δx⁻(II),10] * Δ^2 + eps(0.01))
        end
        @inbounds @threads for II in b_top_v[1]
            pII = lexicographic(II, n+1)
            @inbounds iMGyS[pII,pII] = 1. / (SOL[δy⁻(II),11] * Δ^2 + eps(0.01))
            @inbounds iMGyL[pII,pII] = 1. / (LIQ[δy⁻(II),11] * Δ^2 + eps(0.01))
        end
        @inbounds @threads for II in all_indices_u
            pII = lexicographic(II, n)
            @inbounds MuS[pII,pII] = SOLu[II,5] * Δ^2
            @inbounds MuL[pII,pII] = LIQu[II,5] * Δ^2
        end
        @inbounds @threads for II in all_indices_v
            pII = lexicographic(II, n+1)
            @inbounds MvS[pII,pII] = SOLv[II,5] * Δ^2
            @inbounds MvL[pII,pII] = LIQv[II,5] * Δ^2
        end

        get_curvature(u, liq_projection, κ, B, BT, Δ, MIXED)
        if save_radius
            local radius = zeros(max_iterations+1)
            radius[1] = find_radius(u, MIXED, iso, B, BT, H, inside, Δ)
        end
        if hill
            local radius = zeros(max_iterations+1)
            a = zeros(length(MIXED))
            for i in 1:length(MIXED)
                a[i] = liq_projection[MIXED[i]].pos.y
            end
            radius[1] = mean(a)
        end
    elseif !levelset
        MIXED = [CartesianIndex(-1,-1)]
        MIXED_u = [CartesianIndex(-1,-1)]
        MIXED_v = [CartesianIndex(-1,-1)]
    end

    # Heat
    set_heat!(HTS, DTS, TS, bcTS, BC_TS, κ, ϵ_κ, ϵ_V, V, m, θ₀, τ, aniso,
               LTS, SCUTT, SOL, AS, BS,
               sol_centroid, mid_point, sol_projection,
               all_indices, inside, MIXED, LIQUID, b_left, b_bottom, b_right, b_top, n, Δ)
    set_heat!(HTL, DTL, TL, bcTL, BC_TL, κ, ϵ_κ, ϵ_V, V, m, θ₀, τ, aniso,
               LTL, LCUTT, LIQ, AL, BL,
               liq_centroid, mid_point, sol_projection,
               all_indices, inside, MIXED, SOLID, b_left, b_bottom, b_right, b_top, n, Δ)

    # Navier-Stokes
    set_stokes!(bcpSx, bcpSy, BC_pS, LpS, SCUTp, SOL, n, Δ,
                inside, LIQUID, b_left, b_bottom, b_right, b_top,
                HuS, sol_centroidu, mid_pointu, DuS, bcuS, BC_uS, LuS, SCUTu, SOLu, 
                inside_u, LIQUID_u, b_left_u, b_bottom_u, b_right_u, b_top_u,
                HvS, sol_centroidv, mid_pointv, DvS, bcvS, BC_vS, LvS, SCUTv, SOLv, 
                inside_v, LIQUID_v, b_left_v, b_bottom_v, b_right_v, b_top_v,
                DxuS, DyvS, SCUTDx, SCUTDy, all_indices,
                GxpS, GypS)
                            
    set_stokes!(bcpLx, bcpLy, BC_pL, LpL, LCUTp, LIQ, n, Δ,
                inside, SOLID, b_left, b_bottom, b_right, b_top,
                HuL, liq_centroidu, mid_pointu, DuL, bcuL, BC_uL, LuL, LCUTu, LIQu, 
                inside_u, SOLID_u, b_left_u, b_bottom_u, b_right_u, b_top_u,
                HvL, liq_centroidv, mid_pointv, DvL, bcvL, BC_vL, LvL, LCUTv, LIQv, 
                inside_v, SOLID_v, b_left_v, b_bottom_v, b_right_v, b_top_v,
                DxuL, DyvL, LCUTDx, LCUTDy, all_indices,
                GxpL, GypL)

    ksppS, nsS = init_ksp_solver(LpS, true, ns_vecS)
    kspuS = init_ksp_solver(LuS)
    kspvS = init_ksp_solver(LvS)

    ksppL, nsL = init_ksp_solver(LpL, true, ns_vecL)
    kspuL = init_ksp_solver(LuL)
    kspvL = init_ksp_solver(LvL)

    IIOE(LSA, LSB, u, V, inside, CFL, Δ, n)

    if save_length
        lengthsave[1] = arc_length2(sol_projection, MIXED, Δ)
        κsave[1, :, :] .= κ
    end

    while current_i < max_iterations + 1

        if !stefan
            V .= speed*ones(n,n)
        end

        if heat
            try
                if solid_phase
                    set_heat!(HTS, DTS, θd, bcTS, BC_TS, κ, ϵ_κ, ϵ_V, V, m, θ₀, τ, aniso,
                            LTS, SCUTT, SOL, AS, BS,
                            sol_centroid, mid_point, sol_projection,
                            all_indices, inside, MIXED, LIQUID, b_left, b_bottom, b_right, b_top, n, Δ)
                    TS .= reshape(gmres(AS,(BS*vec(TS) .+ 2.0*τ*SCUTT)), (n,n))
                    # TS .= reshape(gmres(AS,(rhs .+ τ*SCUTT)), (n,n))
                end
                if liquid_phase
                    set_heat!(HTL, DTL, θd, bcTL, BC_TL, κ, ϵ_κ, ϵ_V, V, m, θ₀, τ, aniso,
                            LTL, LCUTT, LIQ, AL, BL,
                            liq_centroid, mid_point, sol_projection,
                            all_indices, inside, MIXED, SOLID, b_left, b_bottom, b_right, b_top, n, Δ)
                    TL .= reshape(gmres(AL,(BL*vec(TL) .+ 2.0*τ*LCUTT)), (n,n))
                    # TL .= reshape(gmres(AL,(rhs .+ τ*LCUTT)), (n,n))
                end
            catch
                @error ("Unphysical temperature field, iteration $current_i")
                break
            end
        end

        if navier_stokes
            if solid_phase
                set_stokes!(bcpSx, bcpSy, BC_pS, LpS, SCUTp, SOL, n, Δ,
                            inside, LIQUID, b_left, b_bottom, b_right, b_top,
                            HuS, sol_centroidu, mid_pointu, DuS, bcuS, BC_uS, LuS, SCUTu, SOLu, 
                            inside_u, LIQUID_u, b_left_u, b_bottom_u, b_right_u, b_top_u,
                            HvS, sol_centroidv, mid_pointv, DvS, bcvS, BC_vS, LvS, SCUTv, SOLv, 
                            inside_v, LIQUID_v, b_left_v, b_bottom_v, b_right_v, b_top_v,
                            DxuS, DyvS, SCUTDx, SCUTDy, all_indices,
                            GxpS, GypS)

                pressure_projection!(pS, uS, vS, ApS, AuS, AvS, LpS, LuS, LvS,
                            SCUTp, SCUTu, SCUTv,
                            ksppS, kspuS, kspvS, nsS, ns_vecS,
                            GxpS, GypS, DxuS, DyvS, SCUTDx, SCUTDy,
                            MpS, iMpS, MuS, MvS, iMGxS, iMGyS, iMDxS, iMDyS,
                            τ, iRe, n, b_right)
                
            end
            if liquid_phase
                set_stokes!(bcpLx, bcpLy, BC_pL, LpL, LCUTp, LIQ, n, Δ,
                           inside, SOLID, b_left, b_bottom, b_right, b_top,
                           HuL, liq_centroidu, mid_pointu, DuL, bcuL, BC_uL, LuL, LCUTu, LIQu, 
                           inside_u, SOLID_u, b_left_u, b_bottom_u, b_right_u, b_top_u,
                           HvL, liq_centroidv, mid_pointv, DvL, bcvL, BC_vL, LvL, LCUTv, LIQv, 
                           inside_v, SOLID_v, b_left_v, b_bottom_v, b_right_v, b_top_v,
                           DxuL, DyvL, LCUTDx, LCUTDy, all_indices,
                           GxpL, GypL)

                pressure_projection!(pL, uL, vL, ApL, AuL, AvL, LpL, LuL, LvL,
                            LCUTp, LCUTu, LCUTv,
                            ksppL, kspuL, kspvL, nsL, ns_vecL,
                            GxpL, GypL, DxuL, DyvL, LCUTDx, LCUTDy,
                            MpL, iMpL, MuL, MvL, iMGxL, iMGyL, iMDxL, iMDyL,
                            τ, iRe, n, b_right)
            end
        end

        if stefan
            Stefan_velocity!(TS, TL, sol_projection, liq_projection, V, MIXED, κ, ϵ_κ, ϵ_V, θd, Δ, m, θ₀, aniso)
            if Vmean
                a = mean(V[MIXED])
                V[MIXED] .= a
            end
            velocity_extension!(V, u, vcat(SOLID,LIQUID), n, Δ, NB, BC_u)
        end

        if advection
            IIOE(LSA, LSB, u, V, inside, CFL, Δ, n)
            try
                u .= reshape(gmres(LSA,(LSB*vec(u))), (n,n))
            catch
                @error ("Inadequate level set function, iteration $current_i")
                break
            end
            if analytical
                u[b_top[1]] .= sqrt.(num.X[b_top[1]] .^ 2 + num.Y[b_top[1]] .^ 2) .- (num.R + speed*current_i*num.τ);
                u[b_bottom[1]] .= sqrt.(num.X[b_bottom[1]] .^ 2 + num.Y[b_bottom[1]] .^ 2) .- (num.R + speed*current_i*num.τ);
                u[b_left[1]] .= sqrt.(num.X[b_left[1]] .^ 2 + num.Y[b_left[1]] .^ 2) .- (num.R + speed*current_i*num.τ);
                u[b_right[1]] .= sqrt.(num.X[b_right[1]] .^ 2 + num.Y[b_right[1]] .^ 2) .- (num.R + speed*current_i*num.τ);
            elseif nb_reinit > 0
                FE_reinit(u, Δ, n, nb_reinit, BC_u, idx)
            end
        end

        if verbose
            if current_i%show_every == 0
                try
                    printstyled(color=:green, @sprintf "\n Current iteration : %d (%d%%) \n" (current_i-1) 100*(current_i-1)/max_iterations)
                    if heat
                        print(@sprintf "V_mean = %.2f  V_max = %.2f  V_min = %.2f\n" mean(V[MIXED]) findmax(V[MIXED])[1] findmin(V[MIXED])[1])
                        print(@sprintf "κ_mean = %.2f  κ_max = %.2f  κ_min = %.2f\n" mean(κ[MIXED]) findmax(κ[MIXED])[1] findmin(κ[MIXED])[1])
                    end
                    if navier_stokes
                        if solid_phase
                            normuS = norm(uS)
                            normvS = norm(vS)
                            normpS = norm(pS.*τ)
                            print("$(@sprintf("norm(uS) %.6e", normuS))\t$(@sprintf("norm(vS) %.6e", normvS))\t$(@sprintf("norm(pS) %.6e", normpS))\n")
                        end
                        if liquid_phase
                            normuL = norm(uL)
                            normvL = norm(vL)
                            normpL = norm(pL.*τ)
                            print("$(@sprintf("norm(uL) %.6e", normuL))\t$(@sprintf("norm(vL) %.6e", normvL))\t$(@sprintf("norm(pL) %.6e", normpL))\n")
                        end
                    end
                catch
                    @show (MIXED)
                end
            end
        end


        if levelset
            marching_squares!(H, iso, u, TS, TL, κ, SOL, LIQ, sol_projection, liq_projection, Δ, L0, B, BT, inside, ϵ, n, n, faces)
            interpolate_LS!(u, uu, uv, B, BT, n, inside,b_left[1], b_bottom[1], b_right[1], b_top[1])
            marching_squares!(H, isou, uu, uS, uL, κu, SOLu, LIQu, sol_projectionu, liq_projectionu, Δ, L0, B, BT, inside_u, ϵ, n+1, n, facesu)
            marching_squares!(H, isov, uv, vS, vL, κv, SOLv, LIQv, sol_projectionv, liq_projectionv, Δ, L0, B, BT, inside_v, ϵ, n, n+1, facesv)

            bcs!(faces, BC_u.left, Δ)
            bcs!(faces, BC_u.right, Δ)
            bcs!(faces, BC_u.bottom, Δ)
            bcs!(faces, BC_u.top, Δ)

            bcs!(facesu, BC_u.left, Δ)
            bcs!(facesu, BC_u.right, Δ)
            bcs!(facesu, BC_u.bottom, Δ)
            bcs!(facesu, BC_u.top, Δ)

            bcs!(facesv, BC_u.left, Δ)
            bcs!(facesv, BC_u.right, Δ)
            bcs!(facesv, BC_u.bottom, Δ)
            bcs!(facesv, BC_u.top, Δ)

            WAS_LIQUID = copy(LIQUID)
            WAS_SOLID = copy(SOLID)
            WAS_LIQUID_u = copy(LIQUID_u)
            WAS_SOLID_u = copy(SOLID_u)
            WAS_LIQUID_v = copy(LIQUID_v)
            WAS_SOLID_v = copy(SOLID_v)

            MIXED, SOLID, LIQUID = get_cells_indices(iso, inside)
            MIXED_u, SOLID_u, LIQUID_u = get_cells_indices(isou, inside_u)
            MIXED_v, SOLID_v, LIQUID_v = get_cells_indices(isov, inside_v)
            
            kill_dead_cells!(TS, LIQUID)
            kill_dead_cells!(TL, SOLID)
            kill_dead_cells!(pS, LIQUID)
            kill_dead_cells!(pL, SOLID)
            kill_dead_cells!(uS, LIQUID_u)
            kill_dead_cells!(uL, SOLID_u)
            kill_dead_cells!(vS, LIQUID_v)
            kill_dead_cells!(vL, SOLID_v)

            @inbounds @threads for II in MIXED
                TS[II] = θd
                TL[II] = θd
            end
            ns_vecS .= 1.
            ns_vecL .= 1.
            @inbounds @threads for II in LIQUID
                pII = lexicographic(II, n)
                @inbounds ns_vecS[pII] = 0.
            end
            @inbounds @threads for II in SOLID
                pII = lexicographic(II, n)
                @inbounds ns_vecL[pII] = 0.
            end

            NB_indices = get_NB_width(MIXED, NB_indices_base)

            get_iterface_location!(gcc, X, Y, iso, u, TS, TL, κ, SOL, LIQ, sol_projection, liq_projection, sol_centroid, liq_centroid, mid_point, Δ, L0, B, BT, idx, MIXED, ϵ, n, faces, periodic_x, periodic_y)
            get_iterface_location!(gfcx, Xu, Yu, isou, uu, uS, uL, κu, SOLu, LIQu, sol_projectionu, liq_projectionu, sol_centroidu, liq_centroidu, mid_pointu, Δ, L0, B, BT, idxu, MIXED_u, ϵ, n, facesu, periodic_x, periodic_y)
            get_iterface_location!(gfcy, Xv, Yv, isov, uv, vS, vL, κv, SOLv, LIQv, sol_projectionv, liq_projectionv, sol_centroidv, liq_centroidv, mid_pointv, Δ, L0, B, BT, idxv, MIXED_v, ϵ, n, facesv, periodic_x, periodic_y)
            
            @inbounds @threads for II in all_indices
                pII = lexicographic(II, n)
                pJJ = lexicographic(II, n+1)
                @inbounds MpS[pII,pII] = SOL[II,5] * Δ^2
                @inbounds MpL[pII,pII] = LIQ[II,5] * Δ^2
                @inbounds iMpS[pII,pII] = 1. / (SOL[II,5] * Δ^2 + eps(0.01))
                @inbounds iMpL[pII,pII] = 1. / (LIQ[II,5] * Δ^2 + eps(0.01))
                @inbounds iMDxS[pII,pII] = 1. / (SOLu[II,10] * Δ^2 + eps(0.01))
                @inbounds iMDyS[pII,pII] = 1. / (SOLv[II,11] * Δ^2 + eps(0.01))
                @inbounds iMDxL[pII,pII] = 1. / (LIQu[II,10] * Δ^2 + eps(0.01))
                @inbounds iMDyL[pII,pII] = 1. / (LIQv[II,11] * Δ^2 + eps(0.01))
                @inbounds iMGxS[pII,pII] = 1. / (SOL[II,8] * Δ^2 + eps(0.01))
                @inbounds iMGyS[pJJ,pJJ] = 1. / (SOL[II,9] * Δ^2 + eps(0.01))
                @inbounds iMGxL[pII,pII] = 1. / (LIQ[II,8] * Δ^2 + eps(0.01))
                @inbounds iMGyL[pJJ,pJJ] = 1. / (LIQ[II,9] * Δ^2 + eps(0.01))
            end
            @inbounds @threads for II in b_right_u[1]
                pII = lexicographic(II, n)
                @inbounds iMGxS[pII,pII] = 1. / (SOL[δx⁻(II),10] * Δ^2 + eps(0.01))
                @inbounds iMGxL[pII,pII] = 1. / (LIQ[δx⁻(II),10] * Δ^2 + eps(0.01))
            end
            @inbounds @threads for II in b_top_v[1]
                pII = lexicographic(II, n+1)
                @inbounds iMGyS[pII,pII] = 1. / (SOL[δy⁻(II),11] * Δ^2 + eps(0.01))
                @inbounds iMGyL[pII,pII] = 1. / (LIQ[δy⁻(II),11] * Δ^2 + eps(0.01))
            end
            @inbounds @threads for II in all_indices_u
                pII = lexicographic(II, n)
                @inbounds MuS[pII,pII] = SOLu[II,5] * Δ^2
                @inbounds MuL[pII,pII] = LIQu[II,5] * Δ^2
            end
            @inbounds @threads for II in all_indices_v
                pII = lexicographic(II, n+1)
                @inbounds MvS[pII,pII] = SOLv[II,5] * Δ^2
                @inbounds MvL[pII,pII] = LIQv[II,5] * Δ^2
            end

            get_curvature(u, liq_projection, κ, B, BT, Δ, MIXED)

            FRESH_L = intersect(MIXED, WAS_SOLID)
            FRESH_S = intersect(MIXED, WAS_LIQUID)
            FRESH_L_u = intersect(MIXED_u, WAS_SOLID_u)
            FRESH_S_u = intersect(MIXED_u, WAS_LIQUID_u)
            FRESH_L_v = intersect(MIXED_v, WAS_SOLID_v)
            FRESH_S_v = intersect(MIXED_v, WAS_LIQUID_v)

            init_fresh_cells!(TS, sol_projection, FRESH_S)
            init_fresh_cells!(TL, liq_projection, FRESH_L)

            if save_radius
                radius[current_i+1] = find_radius(u, MIXED, iso, B, BT, H, inside, Δ)
            end
            if hill
                a = zeros(length(MIXED))
                for i in 1:length(MIXED)
                    a[i] = liq_projection[MIXED[i]].pos.y
                end
                radius[current_i+1] = mean(a)
            end
            if save_length
                lengthsave[current_i+1] = arc_length2(sol_projection, MIXED, Δ)
                κsave[current_i+1, :, :] .= κ
            end
        end

        Vsave[current_i+1, :, :] .= V[:,:]
        usave[current_i+1, :, :] .= u[:,:]
        Tsave[current_i+1, :, :] .= TL[:,:] .+ TS[:,:]
        psave[current_i+1, :, :] .= pL[:,:] .+ pS[:,:]
        Uxsave[current_i+1, :, :] .= uL[:,:] .+ uS[:,:]
        Uysave[current_i+1, :, :] .= vL[:,:] .+ vS[:,:]

        current_i += 1
    end

    if verbose
        try
            printstyled(color=:blue, @sprintf "\n Final iteration : %d (%d%%) \n" (current_i-1) 100*(current_i-1)/max_iterations)
            if heat
                print(@sprintf "V_mean = %.2f  V_max = %.2f  V_min = %.2f  V_stdev = %.5f\n" mean(V[MIXED]) findmax(V[MIXED])[1] findmin(V[MIXED])[1] std(V[MIXED]))
                print(@sprintf "κ_mean = %.2f  κ_max = %.2f  κ_min = %.2f  κ_stdev = %.5f\n" mean(κ[MIXED]) findmax(κ[MIXED])[1] findmin(κ[MIXED])[1] std(κ[MIXED]))
            end
            if navier_stokes
                if solid_phase
                    normuS = norm(uS)
                    normvS = norm(vS)
                    normpS = norm(pS.*τ)
                    print("$(@sprintf("norm(uS) %.6e", normuS))\t$(@sprintf("norm(vS) %.6e", normvS))\t$(@sprintf("norm(pS) %.6e", normpS))\n")
                end
                if liquid_phase
                    normuL = norm(uL)
                    normvL = norm(vL)
                    normpL = norm(pL.*τ)
                    print("$(@sprintf("norm(uL) %.6e", normuL))\t$(@sprintf("norm(vL) %.6e", normvL))\t$(@sprintf("norm(pL) %.6e", normpL))\n")
                end
            end
            print("\n \n")
        catch
            @show (length(MIXED))
        end
    end

    if levelset
        if save_radius || hill

            PETSc.destroy(nsS)
            PETSc.destroy(ksppS)
            PETSc.destroy(kspuS)
            PETSc.destroy(kspvS)

            PETSc.destroy(nsL)
            PETSc.destroy(ksppL)
            PETSc.destroy(kspuL)
            PETSc.destroy(kspvL)

            return MIXED, SOLID, LIQUID, radius
        end

        PETSc.destroy(nsS)
        PETSc.destroy(ksppS)
        PETSc.destroy(kspuS)
        PETSc.destroy(kspvS)

        PETSc.destroy(nsL)
        PETSc.destroy(ksppL)
        PETSc.destroy(kspuL)
        PETSc.destroy(kspvL)
        
        return MIXED, SOLID, LIQUID
    else
        PETSc.destroy(nsS)
        PETSc.destroy(ksppS)
        PETSc.destroy(kspuS)
        PETSc.destroy(kspvS)

        PETSc.destroy(nsL)
        PETSc.destroy(ksppL)
        PETSc.destroy(kspuL)
        PETSc.destroy(kspvL)
        return MIXED
    end
end

function run_backward(num, idx, tmp, fwd, adj;
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

    @unpack L0, A, N, θd, ϵ_κ, ϵ_V, T_inf, τ, L0, NB, n, Δ, CFL, max_iterations, current_i, reinit_every, nb_reinit, ϵ, X, Y, H, B, BT, m, θ₀, aniso = num
    @unpack all_indices, inside, b_left, b_bottom, b_right, b_top = idx
    @unpack SCUTT, LCUTT, LTS, LTL, AS, AL, BS, BL, LSA, LSB, SOL, LIQ, sol_projection, liq_projection, sol_centroid, liq_centroid, mid_point = tmp
    @unpack usave, TSsave, TLsave, Tsave, Vsave, κsave = fwd
    @unpack iso, u, TS, TL, DTS, DTL, κ, V = adj

    local MIXED; local SOLID; local LIQUID;
    local WAS_SOLID; local WAS_LIQUID;
    local NB_indices_base; local NB_indices;
    local FRESH_L; local FRESH_S;

    local faces = zeros(n,n,4);

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
        marching_squares!(H, iso, u, TS, TL, κ, SOL, LIQ, sol_projection, liq_projection, Δ, L0, B, BT, inside, ϵ, n, n, faces)

        bcs!(faces, BC_u.left, Δ)
        bcs!(faces, BC_u.right, Δ)
        bcs!(faces, BC_u.bottom, Δ)
        bcs!(faces, BC_u.top, Δ)

        NB_indices_base = get_NB_width_indices_base(NB)

        MIXED, SOLID, LIQUID = get_cells_indices(iso, inside)
        NB_indices = get_NB_width(MIXED, NB_indices_base)

        get_iterface_location!(gcc, X, Y, iso, u, TS, TL, κ, SOL, LIQ, sol_projection, liq_projection, sol_centroid, liq_centroid, mid_point, Δ, L0, B, BT, idx, MIXED, ϵ, n, faces, periodic_x, periodic_y)
        get_curvature(u, liq_projection, κ, B, BT, Δ, MIXED)
    elseif !levelset
        MIXED = [CartesianIndex(-1,-1)]
    end

    HS = zeros(n,n)
    for II in vcat(b_left[1], b_bottom[1], b_right[1], b_top[1])
        HS[II] = distance(mid_point[II], sol_centroid[II]) * Δ
    end
    
    HL = zeros(n,n)
    for II in vcat(b_left[1], b_bottom[1], b_right[1], b_top[1])
        HL[II] = distance(mid_point[II], liq_centroid[II]) * Δ
    end

    DTS .= θd
    DTL .= θd
    bcS = similar(DTS)
    bcL = similar(DTL)
    apply_curvature(bcS, DTS, κ, ϵ_κ, ϵ_V, V, all_indices)
    apply_curvature(bcL, DTL, κ, ϵ_κ, ϵ_V, V, all_indices)
    if aniso
        apply_anisotropy(bcS, DTS, MIXED, κ, ϵ_κ, ϵ_V, V, m, θ₀, sol_projection)
        apply_anisotropy(bcL, DTL, MIXED, κ, ϵ_κ, ϵ_V, V, m, θ₀, sol_projection)
    end
    bcSx, bcSy = set_bc_bnds(dir, bcS, HS, BC_TS)
    bcLx, bcLy = set_bc_bnds(dir, bcL, HL, BC_TL)

    laplacian!(dir, LTS, SCUTT, bcSx, bcSy, SOL, n, num.Δ, BC_TS, inside, LIQUID,
                b_left[1], b_bottom[1], b_right[1], b_top[1])
    laplacian!(dir, LTL, LCUTT, bcLx, bcLy, LIQ, n, num.Δ, BC_TL, inside, SOLID,
                b_left[1], b_bottom[1], b_right[1], b_top[1])

    while current_i > 1

        if heat
            LCUTT .= zeros(n^2)
            SCUTT .= zeros(n^2)

            try
                if solid_phase
                    HS .= 0.
                    for II in vcat(b_left[1], b_bottom[1], b_right[1], b_top[1])
                        HS[II] = distance(mid_point[II], sol_centroid[II]) * Δ
                    end

                    DTS .= θd
                    apply_curvature(bcS, DTS, κ, ϵ_κ, ϵ_V, V, all_indices)
                    if aniso
                        apply_anisotropy(bcS, DTS, MIXED, κ, ϵ_κ, ϵ_V, V, m, θ₀, sol_projection)
                    end
                    bcSx, bcSy = set_bc_bnds(dir, bcS, HS, BC_TS)

                    laplacian!(dir, LTS, SCUTT, bcSx, bcSy, SOL, n, num.Δ, BC_TS, inside, LIQUID,
                                b_left[1], b_bottom[1], b_right[1], b_top[1])
                    crank_nicolson!(LTS, AS, BS, SOL, τ, n, Δ, all_indices)
                    TS .= reshape(gmres(AS,(BS*vec(TS) + 2.0*τ*SCUTT)), (n,n))
                end
                if liquid_phase
                    HL .= 0.
                    for II in vcat(b_left[1], b_bottom[1], b_right[1], b_top[1])
                        HL[II] = distance(mid_point[II], liq_centroid[II]) * Δ
                    end

                    DTL .= θd
                    apply_curvature(bcL, DTL, κ, ϵ_κ, ϵ_V, V, all_indices)
                    if aniso
                        apply_anisotropy(bcL, DTL, MIXED, κ, ϵ_κ, ϵ_V, V, m, θ₀, sol_projection)
                    end
                    bcLx, bcLy = set_bc_bnds(dir, bcL, HL, BC_TL)

                    laplacian!(dir, LTL, LCUTT, bcLx, bcLy, LIQ, n, num.Δ, BC_TL, inside, SOLID,
                                b_left[1], b_bottom[1], b_right[1], b_top[1])
                    crank_nicolson!(LTL, AL, BL, LIQ, τ, n, Δ, all_indices)
                    TL .= reshape(gmres(AL,(BL*vec(TL) + 2.0*τ*LCUTT)), (n,n))
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
            marching_squares!(H, iso, u, TS, TL, κ, SOL, LIQ, sol_projection, liq_projection, Δ, L0, B, BT, inside, ϵ, n, n, faces)

            bcs!(faces, BC_u.left, Δ)
            bcs!(faces, BC_u.right, Δ)
            bcs!(faces, BC_u.bottom, Δ)
            bcs!(faces, BC_u.top, Δ)

            WAS_LIQUID = copy(LIQUID)
            WAS_SOLID = copy(SOLID)

            MIXED, SOLID, LIQUID = get_cells_indices(iso, inside)
            NB_indices = Flower.get_NB_width(MIXED, NB_indices_base)

            get_iterface_location!(gcc, X, Y, iso, u, TS, TL, κ, SOL, LIQ, sol_projection, liq_projection, sol_centroid, liq_centroid, mid_point, Δ, L0, B, BT, idx, MIXED, ϵ, n, faces, periodic_x, periodic_y)

            FRESH_L = intersect(MIXED, WAS_SOLID)
            FRESH_S = intersect(MIXED, WAS_LIQUID)

            init_fresh_cells!(TS, sol_projection, FRESH_S)
            init_fresh_cells!(TL, liq_projection, FRESH_L)

            get_curvature(u, liq_projection, κ, B, BT, Δ, MIXED)
        end

        current_i -= 1
        κ .= κsave[current_i, :, :]
        u .= usave[current_i, :, :]
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
