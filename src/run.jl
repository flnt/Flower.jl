function run_forward(num, idx, tmp, fwd;
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

    @unpack L0, A, N, θd, ϵ_κ, ϵ_V, T_inf, τ, L0, NB, n, Δ, CFL, max_iterations, current_i, reinit_every, nb_reinit, ϵ, H, B, BT, m, θ₀, aniso = num
    @unpack inside, b_left, b_bottom, b_right, b_top = idx
    @unpack SCUT, LCUT, AS, AL, BS, BL, LSA, LSB, SOL, LIQ, sol_projection, liq_projection = tmp
    @unpack iso, u, TS, TL, Tall, V, κ, usave, TSsave, TLsave, Tsave, Vsave, κsave, lengthsave = fwd

    local MIXED; local SOLID; local LIQUID;
    local WAS_SOLID; local WAS_LIQUID;
    local NB_indices_base; local NB_indices;
    local FRESH_L; local FRESH_S;

    local faces = zeros(n,n,4);

    if periodic_x
        BC_TS.left.ind = BC_TL.left.ind = BC_u.left.ind = idx.periodicL;
        BC_TS.right.ind = BC_TL.right.ind = BC_u.right.ind = idx.periodicR;
        BC_TS.left.f = BC_TL.left.f = BC_u.left.f = BC_TS.right.f = BC_TL.right.f = BC_u.right.f = periodic
    else
        BC_TS.left.ind = BC_TL.left.ind = BC_u.left.ind = b_left;
        BC_TS.right.ind = BC_TL.right.ind = BC_u.right.ind = b_right;
    end

    if periodic_y
        BC_TS.bottom.ind = BC_TL.bottom.ind = BC_u.bottom.ind = idx.periodicB;
        BC_TS.top.ind = BC_TL.top.ind = BC_u.top.ind = idx.periodicT;
        BC_TS.bottom.f = BC_TL.bottom.f = BC_u.bottom.f = BC_TS.top.f = BC_TL.top.f = BC_u.top.f = periodic
    else
        BC_TS.bottom.ind = BC_TL.bottom.ind = BC_u.bottom.ind = b_bottom;
        BC_TS.top.ind = BC_TL.top.ind = BC_u.top.ind = b_top;
    end

    usave[1, :, :] .= u[:,:]
    Tsave[1, :, :] .= TL[:,:] .+ TS[:,:]

    if levelset
        marching_squares!(H, iso, u, TS, TL, κ, SOL, LIQ, sol_projection, liq_projection, Δ, L0, B, BT, inside, ϵ, n, faces)

        bcs!(faces, BC_u.left, Δ)
        bcs!(faces, BC_u.right, Δ)
        bcs!(faces, BC_u.bottom, Δ)
        bcs!(faces, BC_u.top, Δ)

        NB_indices_base = get_NB_width_indices_base(NB)

        MIXED, SOLID, LIQUID = get_cells_indices(iso, inside)
        NB_indices = get_NB_width(MIXED, NB_indices_base)

        get_iterface_location!(H, iso, u, TS, TL, κ, SOL, LIQ, sol_projection, liq_projection, Δ, L0, B, BT, MIXED, ϵ, n, faces)
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
    end

    @sync begin
        @spawn bcs!(SOL, BC_u.left, Δ)
        @spawn bcs!(SOL, BC_u.right, Δ)
        @spawn bcs!(SOL, BC_u.bottom, Δ)
        @spawn bcs!(SOL, BC_u.top, Δ)

        @spawn bcs!(LIQ, BC_u.left, Δ)
        @spawn bcs!(LIQ, BC_u.right, Δ)
        @spawn bcs!(LIQ, BC_u.bottom, Δ)
        @spawn bcs!(LIQ, BC_u.top, Δ)
    end

    crank_nicolson(SCUT, LCUT, SOL, LIQ, AS, AL, BS, BL, n, inside, CFL, θd, κ, ϵ_κ, ϵ_V, V)

    if aniso
        crank_nicolson(SCUT, LCUT, SOL, LIQ, AS, AL, BS, BL, n, MIXED, CFL, θd, κ, ϵ_κ, ϵ_V, V, m, θ₀, aniso, sol_projection)
    end

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
            LCUT .= zeros(n^2)
            SCUT .= zeros(n^2)
            crank_nicolson(SCUT, LCUT, SOL, LIQ, AS, AL, BS, BL, n, inside, CFL, θd, κ, ϵ_κ, ϵ_V, V)
            if aniso
                crank_nicolson(SCUT, LCUT, SOL, LIQ, AS, AL, BS, BL, n, MIXED, CFL, θd, κ, ϵ_κ, ϵ_V, V, m, θ₀, aniso, sol_projection)
            end
            try
                if solid_phase
                    TS .= reshape(cg(AS,(BS*vec(TS) + 2.0*SCUT)), (n,n))
                    #TS .= reshape(AS\(BS*vec(TS) + 2.0*SCUT), (n,n))
                    bcs!(TS, BC_TS.left, Δ)
                    bcs!(TS, BC_TS.right, Δ)
                    bcs!(TS, BC_TS.bottom, Δ)
                    bcs!(TS, BC_TS.top, Δ)
                end
                if liquid_phase
                    TL .= reshape(cg(AL,(BL*vec(TL) + 2.0*LCUT)), (n,n))
                    #TL .= reshape(AL\(BL*vec(TL) + 2.0*LCUT), (n,n))
                    bcs!(TL, BC_TL.left, Δ)
                    bcs!(TL, BC_TL.right, Δ)
                    bcs!(TL, BC_TL.bottom, Δ)
                    bcs!(TL, BC_TL.top, Δ)
                end
            catch
                @error ("Unphysical temperature field, iteration $current_i")
                break
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
                u .= reshape(cg(LSA,(LSB*vec(u))), (n,n))
                #u .= reshape(LSA\(LSB*vec(u)), (n,n))
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
                    print(@sprintf "V_mean = %.2f  V_max = %.2f  V_min = %.2f\n" mean(V[MIXED]) findmax(V[MIXED])[1] findmin(V[MIXED])[1])
                    print(@sprintf "κ_mean = %.2f  κ_max = %.2f  κ_min = %.2f\n" mean(κ[MIXED]) findmax(κ[MIXED])[1] findmin(κ[MIXED])[1])
                catch
                    @show (MIXED)
                end
            end
        end


        if levelset
            marching_squares!(H, iso, u, TS, TL, κ, SOL, LIQ, sol_projection, liq_projection, Δ, L0, B, BT, inside, ϵ, n, faces)

            bcs!(faces, BC_u.left, Δ)
            bcs!(faces, BC_u.right, Δ)
            bcs!(faces, BC_u.bottom, Δ)
            bcs!(faces, BC_u.top, Δ)

            WAS_LIQUID = copy(LIQUID)
            WAS_SOLID = copy(SOLID)

            MIXED, SOLID, LIQUID = get_cells_indices(iso, inside)
            NB_indices = get_NB_width(MIXED, NB_indices_base)

            get_iterface_location!(H, iso, u, TS, TL, κ, SOL, LIQ, sol_projection, liq_projection, Δ, L0, B, BT, MIXED, ϵ, n, faces)
            get_curvature(u, liq_projection, κ, B, BT, Δ, MIXED)

            FRESH_L = intersect(MIXED, WAS_SOLID)
            FRESH_S = intersect(MIXED, WAS_LIQUID)

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

        @sync begin
            @spawn bcs!(SOL, BC_u.left, Δ)
            @spawn bcs!(SOL, BC_u.right, Δ)
            @spawn bcs!(SOL, BC_u.bottom, Δ)
            @spawn bcs!(SOL, BC_u.top, Δ)

            @spawn bcs!(LIQ, BC_u.left, Δ)
            @spawn bcs!(LIQ, BC_u.right, Δ)
            @spawn bcs!(LIQ, BC_u.bottom, Δ)
            @spawn bcs!(LIQ, BC_u.top, Δ)
        end

        Vsave[current_i+1, :, :] .= V[:,:]
        usave[current_i+1, :, :] .= u[:,:]
        Tsave[current_i+1, :, :] .= TL[:,:] .+ TS[:,:]

        current_i += 1
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

    if levelset
        if save_radius || hill
            return MIXED, SOLID, LIQUID, radius
        end
        return MIXED, SOLID, LIQUID
    else
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

    @unpack L0, A, N, θd, ϵ_κ, ϵ_V, T_inf, τ, L0, NB, n, Δ, CFL, max_iterations, current_i, reinit_every, nb_reinit, ϵ, H, B, BT, m, θ₀, aniso = num
    @unpack inside, b_left, b_bottom, b_right, b_top = idx
    @unpack SCUT, LCUT, AS, AL, BS, BL, LSA, LSB, SOL, LIQ, sol_projection, liq_projection = tmp
    @unpack usave, TSsave, TLsave, Tsave, Vsave, κsave = fwd
    @unpack iso, u, TS, TL, κ, V = adj

    local MIXED; local SOLID; local LIQUID;
    local WAS_SOLID; local WAS_LIQUID;
    local NB_indices_base; local NB_indices;
    local FRESH_L; local FRESH_S;

    local faces = zeros(n,n,4);

    if periodic_x
        BC_TS.left.ind = BC_TL.left.ind = BC_u.left.ind = idx.periodicL;
        BC_TS.right.ind = BC_TL.right.ind = BC_u.right.ind = idx.periodicR;
        BC_TS.left.f = BC_TL.left.f = BC_u.left.f = BC_TS.right.f = BC_TL.right.f = BC_u.right.f = periodic
    else
        BC_TS.left.ind = BC_TL.left.ind = BC_u.left.ind = b_left;
        BC_TS.right.ind = BC_TL.right.ind = BC_u.right.ind = b_right;
    end

    if periodic_y
        BC_TS.bottom.ind = BC_TL.bottom.ind = BC_u.bottom.ind = idx.periodicB;
        BC_TS.top.ind = BC_TL.top.ind = BC_u.top.ind = idx.periodicT;
        BC_TS.bottom.f = BC_TL.bottom.f = BC_u.bottom.f = BC_TS.top.f = BC_TL.top.f = BC_u.top.f = periodic
    else
        BC_TS.bottom.ind = BC_TL.bottom.ind = BC_u.bottom.ind = b_bottom;
        BC_TS.top.ind = BC_TL.top.ind = BC_u.top.ind = b_top;
    end

    current_i = max_iterations + 1

    if levelset
        marching_squares!(H, iso, u, TS, TL, κ, SOL, LIQ, sol_projection, liq_projection, Δ, L0, B, BT, inside, ϵ, n, faces)

        bcs!(faces, BC_u.left, Δ)
        bcs!(faces, BC_u.right, Δ)
        bcs!(faces, BC_u.bottom, Δ)
        bcs!(faces, BC_u.top, Δ)

        NB_indices_base = get_NB_width_indices_base(NB)

        MIXED, SOLID, LIQUID = get_cells_indices(iso, inside)
        NB_indices = get_NB_width(MIXED, NB_indices_base)

        get_iterface_location!(H, iso, u, TS, TL, κ, SOL, LIQ, sol_projection, liq_projection, Δ, L0, B, BT, MIXED, ϵ, n, faces)
        get_curvature(u, liq_projection, κ, B, BT, Δ, MIXED)
    elseif !levelset
        MIXED = [CartesianIndex(-1,-1)]
    end

    @sync begin
        @spawn bcs!(SOL, BC_u.left, Δ)
        @spawn bcs!(SOL, BC_u.right, Δ)
        @spawn bcs!(SOL, BC_u.bottom, Δ)
        @spawn bcs!(SOL, BC_u.top, Δ)

        @spawn bcs!(LIQ, BC_u.left, Δ)
        @spawn bcs!(LIQ, BC_u.right, Δ)
        @spawn bcs!(LIQ, BC_u.bottom, Δ)
        @spawn bcs!(LIQ, BC_u.top, Δ)
    end

    crank_nicolson(SCUT, LCUT, SOL, LIQ, AS, AL, BS, BL, n, inside, CFL, θd, κ, ϵ_κ, ϵ_V, V)

    if aniso
        crank_nicolson(SCUT, LCUT, SOL, LIQ, AS, AL, BS, BL, n, MIXED, CFL, θd, κ, ϵ_κ, ϵ_V, V, m, θ₀, aniso, sol_projection)
    end

    while current_i > 1

        if heat
            LCUT .= zeros(n^2)
            SCUT .= zeros(n^2)
            crank_nicolson(SCUT, LCUT, SOL, LIQ, AS, AL, BS, BL, n, inside, CFL, θd, κ, ϵ_κ, ϵ_V, V)
            if aniso
                crank_nicolson(SCUT, LCUT, SOL, LIQ, AS, AL, BS, BL, n, MIXED, CFL, θd, κ, ϵ_κ, ϵ_V, V, m, θ₀, aniso, sol_projection)
            end
            try
                if solid_phase
                    TS .= reshape(cg(AS,(BS*vec(TS) + 2.0*SCUT)), (n,n))
                    #TS .= reshape(AS\(BS*vec(TS) + 2.0*SCUT), (n,n))
                    bcs!(TS, BC_TS.left, Δ)
                    bcs!(TS, BC_TS.right, Δ)
                    bcs!(TS, BC_TS.bottom, Δ)
                    bcs!(TS, BC_TS.top, Δ)
                end
                if liquid_phase
                    #@show (TL[b_top[1]])
                    #TL .= reshape(cg(AL,(BL*vec(TL) + 2.0*LCUT)), (n,n))
                    TL .= reshape(AL\(BL*vec(TL) + 2.0*LCUT), (n,n))
                    bcs!(TL, BC_TL.left, Δ)
                    bcs!(TL, BC_TL.right, Δ)
                    bcs!(TL, BC_TL.bottom, Δ)
                    bcs!(TL, BC_TL.top, Δ)
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
            marching_squares!(H, iso, u, TS, TL, κ, SOL, LIQ, sol_projection, liq_projection, Δ, L0, B, BT, inside, ϵ, n, faces)

            bcs!(faces, BC_u.left, Δ)
            bcs!(faces, BC_u.right, Δ)
            bcs!(faces, BC_u.bottom, Δ)
            bcs!(faces, BC_u.top, Δ)

            WAS_LIQUID = copy(LIQUID)
            WAS_SOLID = copy(SOLID)

            MIXED, SOLID, LIQUID = get_cells_indices(iso, inside)
            NB_indices = Flower.get_NB_width(MIXED, NB_indices_base)

            get_iterface_location!(H, iso, u, TS, TL, κ, SOL, LIQ, sol_projection, liq_projection, Δ, L0, B, BT, MIXED, ϵ, n, faces)

            FRESH_L = intersect(MIXED, WAS_SOLID)
            FRESH_S = intersect(MIXED, WAS_LIQUID)

            init_fresh_cells!(TS, sol_projection, FRESH_S)
            init_fresh_cells!(TL, liq_projection, FRESH_L)

            get_curvature(u, liq_projection, κ, B, BT, Δ, MIXED)
        end

        @sync begin
            @spawn bcs!(SOL, BC_u.left, Δ)
            @spawn bcs!(SOL, BC_u.right, Δ)
            @spawn bcs!(SOL, BC_u.bottom, Δ)
            @spawn bcs!(SOL, BC_u.top, Δ)

            @spawn bcs!(LIQ, BC_u.left, Δ)
            @spawn bcs!(LIQ, BC_u.right, Δ)
            @spawn bcs!(LIQ, BC_u.bottom, Δ)
            @spawn bcs!(LIQ, BC_u.top, Δ)
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
