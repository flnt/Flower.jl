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
    )

    @unpack L0, A, N, θd, ϵ_κ, ϵ_V, T_inf, num.τ, L0, NB, max_iterations, num.current_i, reinit_every, nb_reinit, ϵ, m, θ₀, aniso = num
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

    num.current_i = max_iterations + 1

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

    laplacian!(dir, num, opS.LT, opS.CUTT, bcSx, bcSy, geoS.dcap, ny, BC_TS, inside, LIQUID,
                MIXED, b_left[1], b_bottom[1], b_right[1], b_top[1])
    laplacian!(dir, num, opL.LT, opL.CUTT, bcLx, bcLy, geoL.dcap, ny, BC_TL, inside, SOLID,
                MIXED, b_left[1], b_bottom[1], b_right[1], b_top[1])

    while num.current_i > 1

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

                    laplacian!(dir, num, opS.LT, opS.CUTT, bcSx, bcSy, geoS.dcap, ny, BC_TS, inside, LIQUID,
                                MIXED, b_left[1], b_bottom[1], b_right[1], b_top[1])
                    crank_nicolson!(num, grid, geoS, opS)
                    TS .= reshape(gmres(opS.A,(opS.B*vec(TS) + 2.0*num.τ*opS.CUTT)), (ny,nx))
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

                    laplacian!(dir, num, opL.LT, opL.CUTT, bcLx, bcLy, geoL.dcap, ny, BC_TL, inside, SOLID,
                                MIXED, b_left[1], b_bottom[1], b_right[1], b_top[1])
                    crank_nicolson!(num, grid, geoL, opL)
                    TL .= reshape(gmres(opL.A,(opL.B*vec(TL) + 2.0*num.τ*opL.CUTT)), (ny,nx))
                end
            catch
                @error ("Unphysical temperature field, iteration $num.current_i")
                break
            end
        end

        if verbose
            if num.current_i%show_every == 0
                try
                    printstyled(color=:green, @sprintf "\n Current iteration : %d (%d%%) \n" (num.current_i-1) 100*(num.current_i-1)/max_iterations)
                    printstyled(color=:green, @sprintf "\n CFL : %.2e CFL : %.2e num.τ : %.2e\n" CFL max(abs.(V)..., abs.(phL.u)..., abs.(phL.v)..., abs.(phS.u)..., abs.(phS.v)...)*num.τ/Δ num.τ)

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

        num.current_i -= 1
        κ .= κsave[num.current_i,:,:]
        u .= usave[num.current_i,:,:]
    end

    if verbose
        try
            printstyled(color=:blue, @sprintf "\n Final iteration : %d (%d%%) \n" (num.current_i-1) 100*(num.current_i-1)/max_iterations)
            print(@sprintf "V_mean = %.2f  V_max = %.2f  V_min = %.2f  V_stdev = %.5f\n" mean(V[MIXED]) findmax(V[MIXED])[1] findmin(V[MIXED])[1] std(V[MIXED]))
            print(@sprintf "κ_mean = %.2f  κ_max = %.2f  κ_min = %.2f  κ_stdev = %.5f\n" mean(κ[MIXED]) findmax(κ[MIXED])[1] findmin(κ[MIXED])[1] std(κ[MIXED]))
            print("\n \n")
        catch
            @show (length(MIXED))
        end
    end
end
