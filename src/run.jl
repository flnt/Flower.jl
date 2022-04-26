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
    @unpack all_indices, inside, b_left, b_bottom, b_right, b_top = idx
    @unpack SCUT, LCUT, LTS, LTL, AS, AL, BS, BL, LSA, LSB, SOL, LIQ, sol_projection, liq_projection, sol_centroid, liq_centroid, mid_point = tmp
    @unpack iso, u, TS, TL, Tall, DTS, DTL, V, κ, usave, TSsave, TLsave, Tsave, Vsave, κsave, lengthsave = fwd

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

        get_iterface_location!(H, iso, u, TS, TL, κ, SOL, LIQ, sol_projection, liq_projection, sol_centroid, liq_centroid, mid_point, Δ, L0, B, BT, idx, MIXED, ϵ, n, faces, periodic_x, periodic_y)
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

    HS = similar(DTS)
    HS .= 0.
    for II in vcat(b_left[1], b_bottom[1], b_right[1], b_top[1])
        HS[II] = distance(mid_point[II], sol_centroid[II]) * Δ
    end
    
    HL = similar(DTL)
    HL .= 0.
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
    bcSx, bcSy = set_bc_bnds(bcS, HS, BC_TS)
    bcLx, bcLy = set_bc_bnds(bcL, HL, BC_TL)

    laplacian!(dir, LTS, SCUT, bcSx, bcSy, SOL, n, num.Δ, BC_TS, all_indices, LIQUID,
                b_left[1], b_bottom[1], b_right[1], b_top[1])
    laplacian!(dir, LTL, LCUT, bcLx, bcLy, LIQ, n, num.Δ, BC_TL, all_indices, SOLID,
                b_left[1], b_bottom[1], b_right[1], b_top[1])
    crank_nicolson!(LTS, AS, BS, SOL, τ, n, Δ, all_indices)
    crank_nicolson!(LTL, AL, BL, LIQ, τ, n, Δ, all_indices)

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
                    bcSx, bcSy = set_bc_bnds(bcS, HS, BC_TS)

                    laplacian!(dir, LTS, SCUT, bcSx, bcSy, SOL, n, num.Δ, BC_TS, all_indices, LIQUID,
                                b_left[1], b_bottom[1], b_right[1], b_top[1])
                    crank_nicolson!(LTS, AS, BS, SOL, τ, n, Δ, all_indices)
                    TS .= reshape(gmres(AS,(BS*vec(TS) .+ 2.0*τ*SCUT)), (n,n))
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
                    bcLx, bcLy = set_bc_bnds(bcL, HL, BC_TL)

                    laplacian!(dir, LTL, LCUT, bcLx, bcLy, LIQ, n, num.Δ, BC_TL, all_indices, SOLID,
                                b_left[1], b_bottom[1], b_right[1], b_top[1])
                    crank_nicolson!(LTL, AL, BL, LIQ, τ, n, Δ, all_indices)
                    TL .= reshape(gmres(AL,(BL*vec(TL) .+ 2.0*τ*LCUT)), (n,n))
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

            get_iterface_location!(H, iso, u, TS, TL, κ, SOL, LIQ, sol_projection, liq_projection, sol_centroid, liq_centroid, mid_point, Δ, L0, B, BT, idx, MIXED, ϵ, n, faces, periodic_x, periodic_y)
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
    @unpack all_indices, inside, b_left, b_bottom, b_right, b_top = idx
    @unpack SCUT, LCUT, LTS, LTL, AS, AL, BS, BL, LSA, LSB, SOL, LIQ, sol_projection, liq_projection, sol_centroid, liq_centroid, mid_point = tmp
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
        marching_squares!(H, iso, u, TS, TL, κ, SOL, LIQ, sol_projection, liq_projection, Δ, L0, B, BT, inside, ϵ, n, faces)

        bcs!(faces, BC_u.left, Δ)
        bcs!(faces, BC_u.right, Δ)
        bcs!(faces, BC_u.bottom, Δ)
        bcs!(faces, BC_u.top, Δ)

        NB_indices_base = get_NB_width_indices_base(NB)

        MIXED, SOLID, LIQUID = get_cells_indices(iso, inside)
        NB_indices = get_NB_width(MIXED, NB_indices_base)

        get_iterface_location!(H, iso, u, TS, TL, κ, SOL, LIQ, sol_projection, liq_projection, sol_centroid, liq_centroid, mid_point, Δ, L0, B, BT, idx, MIXED, ϵ, n, faces, periodic_x, periodic_y)
        get_curvature(u, liq_projection, κ, B, BT, Δ, MIXED)
    elseif !levelset
        MIXED = [CartesianIndex(-1,-1)]
    end

    HS = similar(DTS)
    HS .= 0.
    for II in vcat(b_left[1], b_bottom[1], b_right[1], b_top[1])
        HS[II] = distance(mid_point[II], sol_centroid[II]) * Δ
    end
    
    HL = similar(DTL)
    HL .= 0.
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
    bcSx, bcSy = set_bc_bnds(bcS, HS, BC_TS)
    bcLx, bcLy = set_bc_bnds(bcL, HL, BC_TL)

    laplacian!(dir, LTS, SCUT, bcSx, bcSy, SOL, n, num.Δ, BC_TS, all_indices, LIQUID,
                b_left[1], b_bottom[1], b_right[1], b_top[1])
    laplacian!(dir, LTL, LCUT, bcLx, bcLy, LIQ, n, num.Δ, BC_TL, all_indices, SOLID,
                b_left[1], b_bottom[1], b_right[1], b_top[1])
    crank_nicolson!(LTS, AS, BS, SOL, τ, n, Δ, all_indices)
    crank_nicolson!(LTL, AL, BL, LIQ, τ, n, Δ, all_indices)

    while current_i > 1

        if heat
            LCUT .= zeros(n^2)
            SCUT .= zeros(n^2)

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
                    bcSx, bcSy = set_bc_bnds(bcS, HS, BC_TS)

                    laplacian!(dir, LTS, SCUT, bcSx, bcSy, SOL, n, num.Δ, BC_TS, all_indices, LIQUID,
                                b_left[1], b_bottom[1], b_right[1], b_top[1])
                    crank_nicolson!(LTS, AS, BS, SOL, τ, n, Δ, all_indices)
                    TS .= reshape(gmres(AS,(BS*vec(TS) + 2.0*τ*SCUT)), (n,n))
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
                    bcLx, bcLy = set_bc_bnds(bcL, HL, BC_TL)

                    laplacian!(dir, LTL, LCUT, bcLx, bcLy, LIQ, n, num.Δ, BC_TL, all_indices, SOLID,
                                b_left[1], b_bottom[1], b_right[1], b_top[1])
                    crank_nicolson!(LTL, AL, BL, LIQ, τ, n, Δ, all_indices)
                    TL .= reshape(gmres(AL,(BL*vec(TL) + 2.0*τ*LCUT)), (n,n))
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

        pII = lexicographic(II, n)
            get_iterface_location!(H, iso, u, TS, TL, κ, SOL, LIQ, sol_projection, liq_projection, sol_centroid, liq_centroid, mid_point, Δ, L0, B, BT, idx, MIXED, ϵ, n, faces, periodic_x, periodic_y)

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
