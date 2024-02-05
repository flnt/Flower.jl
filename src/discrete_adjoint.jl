function coloring(grid, stencil, periodic_x, periodic_y)
    @unpack nx, ny, ind = grid

    st = (stencil - 1) ÷ 2

    ix = last.(getproperty.(ind.all_indices, :I))
    iy = first.(getproperty.(ind.all_indices, :I))
    colors = (ix .+ 2stencil*iy) .% (stencil^2)
    graph = zeros(Int64,2*nx*ny)
    px0 = ones(Int64,4)
    py0 = ones(Int64,4)
    pxf = nx.*ones(Int64,4)
    pyf = ny.*ones(Int64,4)
    if !periodic_x && !periodic_y
        nloops = 1
    elseif periodic_x && periodic_y
        nloops = 4
        pxf[1] = nx-st
        pyf[1] = ny-st
        px0[2] = nx-st+1
        pyf[2] = ny-st
        pxf[3] = nx-st
        py0[3] = ny-st+1
        px0[4] = nx-st+1
        py0[4] = ny-st+1
    elseif periodic_x
        nloops = 2
        pxf[1] = nx-st
        px0[2] = nx-st+1
    else
        nloops = 2
        pyf[1] = ny-st
        py0[2] = ny-st+1
    end

    return colors, graph, nloops, px0, pxf, py0, pyf
end

function assign_color!(grid, stencil, graph, ind_pert, per_x, per_y)
    @unpack nx, ny, ind = grid
    @unpack all_indices = ind

    st = (stencil - 1) ÷ 2

    graph .= 0
    if !per_x && !per_y
        @inbounds for II in all_indices
            i = lexicographic(II, ny)
            @inbounds for JJ in ind_pert
                j = lexicographic(JJ, ny)
                if abs(II[1] - JJ[1]) <= st && abs(II[2] - JJ[2]) <= st
                    graph[i] = j
                    graph[i+nx*ny] = j
                    break
                end
            end
        end
    elseif !per_x && per_y
        @inbounds for II in all_indices
            i = lexicographic(II, ny)
            @inbounds for JJ in ind_pert
                j = lexicographic(JJ, ny)
                if (abs(II[1] - JJ[1]) <= st || abs(II[1] - JJ[1]) >= ny-st) && abs(II[2] - JJ[2]) <= st
                    graph[i] = j
                    graph[i+nx*ny] = j
                    break
                end
            end
        end
    elseif per_x && !per_y
        @inbounds for II in all_indices
            i = lexicographic(II, ny)
            @inbounds for JJ in ind_pert
                j = lexicographic(JJ, ny)
                if abs(II[1] - JJ[1]) <= st && (abs(II[2] - JJ[2]) <= st || abs(II[2] - JJ[2]) >= nx-st)
                    graph[i] = j
                    graph[i+nx*ny] = j
                    break
                end
            end
        end
    else @inbounds for II in all_indices
            i = lexicographic(II, ny)
            @inbounds for JJ in ind_pert
                j = lexicographic(JJ, ny)
                if (abs(II[1] - JJ[1]) <= st || abs(II[1] - JJ[1]) >= ny-st) && (abs(II[2] - JJ[2]) <= st || abs(II[2] - JJ[2]) >= nx-st)
                    graph[i] = j
                    graph[i+nx*ny] = j
                    break
                end
            end
        end
    end

    return nothing
end

function assign_color!(grid, grid_pert, stencil, graph, ind_pert, per_x, per_y)
    @unpack nx, ny, ind = grid
    @unpack all_indices = ind

    st = (stencil - 1) ÷ 2

    graph .= 0
    if !per_x && !per_y
        @inbounds for II in all_indices
            i = lexicographic(II, ny)
            @inbounds for JJ in ind_pert
                j = lexicographic(JJ, grid_pert.ny)
                if abs(II[1] - JJ[1]) <= st && abs(II[2] - JJ[2]) <= st
                    graph[i] = j
                    graph[i+nx*ny] = j
                    break
                end
            end
        end
    elseif !per_x && per_y
        @inbounds for II in all_indices
            i = lexicographic(II, ny)
            @inbounds for JJ in ind_pert
                j = lexicographic(JJ, grid_pert.ny)
                if (abs(II[1] - JJ[1]) <= st || abs(II[1] - JJ[1]) >= ny-st) && abs(II[2] - JJ[2]) <= st
                    graph[i] = j
                    graph[i+nx*ny] = j
                    break
                end
            end
        end
    elseif per_x && !per_y
        @inbounds for II in all_indices
            i = lexicographic(II, ny)
            @inbounds for JJ in ind_pert
                j = lexicographic(JJ, grid_pert.ny)
                if abs(II[1] - JJ[1]) <= st && (abs(II[2] - JJ[2]) <= st || abs(II[2] - JJ[2]) >= nx-st)
                    graph[i] = j
                    graph[i+nx*ny] = j
                    break
                end
            end
        end
    else @inbounds for II in all_indices
            i = lexicographic(II, ny)
            @inbounds for JJ in ind_pert
                j = lexicographic(JJ, grid_pert.ny)
                if (abs(II[1] - JJ[1]) <= st || abs(II[1] - JJ[1]) >= ny-st) && (abs(II[2] - JJ[2]) <= st || abs(II[2] - JJ[2]) >= nx-st)
                    graph[i] = j
                    graph[i+nx*ny] = j
                    break
                end
            end
        end
    end

    return nothing
end

function relocate!(indices, px0, py0)
    @inbounds @threads for i in eachindex(indices)
        indices[i] += CartesianIndex(py0 - 1, 0)
        indices[i] += CartesianIndex(0, px0 - 1)
    end
end

function Rheat_q0(num, grid, grid_u, grid_v, adj_der,
    TD0_S, TD1_S, A_S, B_S, opC_TS, BC_TS,
    TD0_L, TD1_L, A_L, B_L, opC_TL, BC_TL,
    u0, u1, LSA, LSB, tmpχ_S, tmpχ_L,
    CFL_sc, periodic_x, periodic_y, ϵ_adj, λ, Vmean,
    opS, phS, opL, phL,
    heat_solid_phase, heat_liquid_phase, heat_convection)

    @unpack ϵ, NB = num
    @unpack nx, ny, ind, V, iso, faces, geoS, geoL = grid
    @unpack RheatS_ls, RheatL_ls, RlsS_ls = adj_der

    κ = copy(grid.κ)
    
    RheatS_ls.nzval .= 0.
    RheatL_ls.nzval .= 0.
    RlsS_ls.nzval .= 0.
    uj = copy(u0)

    TS = reshape(veci(TD1_S, grid, 1), (ny, nx))
    TL = reshape(veci(TD1_L, grid, 1), (ny, nx))

    derA = copy(A_S)
    derA.nzval .= 0.
    derB = copy(B_S)
    derB.nzval .= 0.

    derLSA = copy(LSA)
    derLSA.nzval .= 0.
    derLSB = copy(LSB)
    derLSB.nzval .= 0.
    LSAj = copy(LSA)
    LSBj = copy(LSB)

    derχ = copy(opC_TS.χ)
    bc = zeros(2*ny*nx)
    a0 = num.θd .* ones(nx*ny)

    # graph coloring
    stencil = 7
    colors, graph, nloops, px0, pxf, py0, pyf = coloring(grid, stencil, periodic_x, periodic_y)
    @inbounds for color = 0:stencil^2-1
        @inbounds for p = 1:nloops
            indices = findall(colors[py0[p]:pyf[p], px0[p]:pxf[p]] .== color)
            relocate!(indices, px0[p], py0[p])
            assign_color!(grid, stencil, graph, indices, periodic_x, periodic_y)

            uj[indices] .+= ϵ_adj

            # Compute capacities
            update_ls_data(num, grid, grid_u, grid_v, uj, κ, periodic_x, periodic_y)
            update_stefan_velocity(num, grid, uj, TS, TL, periodic_x, periodic_y, λ, Vmean)

            # get perturbed matrices
            if heat_solid_phase
                Aj, Bj, _ = set_heat!(dir, num, grid, opC_TS, geoS, phS, num.θd, BC_TS, grid.ind.MIXED, geoS.projection,
                                    opS, grid_u, grid_u.geoS, grid_v, grid_v.geoS,
                                    periodic_x, periodic_y, heat_convection)
                derA .= (Aj .- A_S) ./ ϵ_adj
                derB .= (Bj .- B_S) ./ ϵ_adj
                derχ .= (opC_TS.χ .-  tmpχ_S) ./ ϵ_adj
                bc[ny*nx+1:end] .= derχ * a0
                Rj = sparse(derA * TD1_S .- derB * TD0_S .- bc)

                # Do NOT parallelize. Sparsity pattern changes due to periodic BCs
                # not being preallocated
                rows = rowvals(Rj)
                for i in nzrange(Rj, 1)
                    @inbounds row = rows[i]
                    j = graph[row]
                    @inbounds RheatS_ls[row,j] = Rj[row]
                end
            end

            if heat_liquid_phase
                Aj, Bj, _ = set_heat!(dir, num, grid, opC_TL, geoL, phL, num.θd, BC_TL, grid.ind.MIXED, geoL.projection,
                                    opL, grid_u, grid_u.geoL, grid_v, grid_v.geoL,    
                                    periodic_x, periodic_y, heat_convection)
                derA .= (Aj .- A_L) ./ ϵ_adj
                derB .= (Bj .- B_L) ./ ϵ_adj
                derχ .= (opC_TL.χ .-  tmpχ_L) ./ ϵ_adj
                bc[ny*nx+1:end] .= derχ * a0
                Rj = sparse(derA * TD1_L .- derB * TD0_L .- bc)

                # Do NOT parallelize. Sparsity pattern changes due to periodic BCs
                # not being preallocated
                rows = rowvals(Rj)
                for i in nzrange(Rj, 1)
                    @inbounds row = rows[i]
                    j = graph[row]
                    @inbounds RheatL_ls[row,j] = Rj[row]
                end
            end

            IIOE_normal!(grid, LSAj, LSBj, uj, V, CFL_sc, periodic_x, periodic_y)
            derLSA .= (LSAj .- LSA) ./ ϵ_adj
            derLSB .= (LSBj .- LSB) ./ ϵ_adj
            utmp = zeros(ny, nx)
            utmp[indices] .= 1.0
            Rj = sparse(derLSA * vec(u1) .- derLSB * vec(u0) .- LSB * vec(utmp))
            
            # Do NOT parallelize. Sparsity pattern changes due to periodic BCs
            # not being preallocated
            rows = rowvals(Rj)
            for i in nzrange(Rj, 1)
                @inbounds row = rows[i]
                j = graph[row]
                @inbounds RlsS_ls[row,j] = Rj[row]
            end

            uj .= u0
        end
    end

    return nothing
end

function Rheat_q1(num, grid, grid_u, grid_v, adj_der, 
    TD_S, TD_L,
    u0, u1, LSA, LSB,
    CFL_sc, periodic_x, periodic_y, ϵ_adj, λ, Vmean,
    heat_solid_phase, heat_liquid_phase)

    @unpack ϵ, NB = num
    @unpack nx, ny, ind, V, faces, iso, geoS, geoL = grid
    @unpack RlsS_TS, RlsS_TL = adj_der

    κ = copy(grid.κ)
    
    RlsS_TS.nzval .= 0.
    RlsS_TL.nzval .= 0.

    TS = reshape(veci(TD_S, grid, 1), (ny, nx))
    TL = reshape(veci(TD_L, grid, 1), (ny, nx))
    TSj = copy(TS)
    TLj = copy(TL)
    
    derLSA = copy(LSA)
    derLSA.nzval .= 0.
    derLSB = copy(LSB)
    derLSB.nzval .= 0.

    LSAj = copy(LSA)
    LSBj = copy(LSB)

    # graph coloring
    stencil = 7
    colors, graph, nloops, px0, pxf, py0, pyf = coloring(grid, stencil, periodic_x, periodic_y)
    @inbounds for color = 0:stencil^2-1
        @inbounds for p = 1:nloops
            indices = findall(colors[py0[p]:pyf[p], px0[p]:pxf[p]] .== color)
            relocate!(indices, px0[p], py0[p])
            js = lexicographic.(indices, ny)
            assign_color!(grid, stencil, graph, indices, periodic_x, periodic_y)

            TSj[js] .+= ϵ_adj
            TLj[js] .+= ϵ_adj

            # Compute capacities
            update_ls_data(num, grid, grid_u, grid_v, u0, κ, periodic_x, periodic_y)

            # get perturbed matrices
            if heat_solid_phase
                update_stefan_velocity(num, grid, u0, TSj, TL, periodic_x, periodic_y, λ, Vmean)
                IIOE_normal!(grid, LSAj, LSBj, u0, V, CFL_sc, periodic_x, periodic_y)
                derLSA .= (LSAj .- LSA) ./ ϵ_adj
                derLSB .= (LSBj .- LSB) ./ ϵ_adj
                Rj = sparse(derLSA * vec(u1) .- derLSB * vec(u0))

                # Do NOT parallelize. Sparsity pattern changes due to periodic BCs
                # not being preallocated
                rows = rowvals(Rj)
                for i in nzrange(Rj, 1)
                    @inbounds row = rows[i]
                    j = graph[row]
                    @inbounds RlsS_TS[row,j] = Rj[row]
                end
            end

            if heat_liquid_phase
                update_stefan_velocity(num, grid, u0, TS, TLj, periodic_x, periodic_y, λ, Vmean)
                IIOE_normal!(grid, LSAj, LSBj, u0, V, CFL_sc, periodic_x, periodic_y)
                derLSA .= (LSAj .- LSA) ./ ϵ_adj
                derLSB .= (LSBj .- LSB) ./ ϵ_adj
                Rj = sparse(derLSA * vec(u1) .- derLSB * vec(u0))

                # Do NOT parallelize. Sparsity pattern changes due to periodic BCs
                # not being preallocated
                rows = rowvals(Rj)
                for i in nzrange(Rj, 1)
                    @inbounds row = rows[i]
                    j = graph[row]
                    @inbounds RlsS_TL[row,j] = Rj[row]
                end
            end

            TSj .= TS
            TLj .= TL
        end
    end

    return nothing
end

function Rproj_q1(num, grid, grid_u, grid_v, adj_der,
    ϕD1_S, opC_pS, BC_pS, ϕD1_L, opC_pL, BC_pL,
    uD0_S, opC_uS, BC_uS, uD0_L, opC_uL, BC_uL,
    vD0_S, opC_vS, BC_vS, vD0_L, opC_vL, BC_vL,
    ucorrD1_S, ucorrD1_L, vcorrD1_S, vcorrD1_L,
    Lum1_S, bc_Lum1_S, Lvm1_S, bc_Lvm1_S,
    Lum1_L, bc_Lum1_L, Lvm1_L, bc_Lvm1_L,
    Mum1_S, Mum1_L, Mvm1_S, Mvm1_L,
    u1, periodic_x, periodic_y, ϵ_adj,
    opS, phS, opL, phL,
    ns_solid_phase, ns_liquid_phase,
    free_surface, ns_advection)

    @unpack Re, τ, σ, g, β, ϵ, NB = num
    @unpack nx, ny, ind, V, iso, faces, geoS, geoL = grid
    @unpack RuS_ls, RuL_ls, RvS_ls, RvL_ls, RpS_ls, RpL_ls,
            RucorrS_ls1, RvcorrS_ls1,
            RucorrL_ls1, RvcorrL_ls1 = adj_der

    iRe = 1.0 / Re
    iτ = 1.0 / τ

    κj = copy(grid.κ)
    uj = copy(u1)

    opC_pj = copy(opC_pL)
    opC_uj = copy(opC_uL)
    opC_vj = copy(opC_vL)

    if ns_solid_phase
        RuS_ls.nzval .= 0.
        RvS_ls.nzval .= 0.
        
        ∇ϕ_x = opC_uS.AxT * opC_uS.Rx * veci(ϕD1_S,grid,1) .+ opC_uS.Gx * veci(ϕD1_S,grid,2)
        ∇ϕ_y = opC_vS.AyT * opC_vS.Ry * veci(ϕD1_S,grid,1) .+ opC_vS.Gy * veci(ϕD1_S,grid,2)
        iMu = Diagonal(1 ./ (opC_uS.M.diag .+ eps(0.01)))
        iMv = Diagonal(1 ./ (opC_vS.M.diag .+ eps(0.01)))
        GxS = τ .* iMu * ∇ϕ_x
        GyS = τ .* iMv * ∇ϕ_y
        
        RpS_ls.nzval .= 0.

        if free_surface
            AuS, _, _, AvS, _, _, AϕS = set_navier_stokes(neu, num, grid, grid.geoS, grid_u, grid_u.geoS, grid_v, grid_v.geoS,
                                                        opC_pS, opC_uS, opC_vS, BC_pS, BC_uS, BC_vS,
                                                        Mum1_S, Mvm1_S, iRe,
                                                        opS, phS,
                                                        periodic_x, periodic_y, ns_advection)[1:7]
        else
            AuS, _, _, AvS, _, _, AϕS = set_navier_stokes(dir, num, grid, grid.geoS, grid_u, grid_u.geoS, grid_v, grid_v.geoS,
                                                        opC_pS, opC_uS, opC_vS, BC_pS, BC_uS, BC_vS,
                                                        Lum1_S, bc_Lum1_S, Lvm1_S, bc_Lvm1_S, Mum1_S, Mvm1_S, iRe,
                                                        opS, phS,
                                                        periodic_x, periodic_y, ns_advection)[1:7]
        end
    
        SmatS = iRe .* strain_rate(opC_uS, opC_vS)
    
        GxT = opC_uS.Gx'
        GyT = opC_vS.Gy'
        BϕκS = [spdiagm(fzeros(grid)); σ .* (GxT * opC_uS.Gx .+ GyT * opC_vS.Gy)]

        RucorrS_ls1.nzval .= 0.
        RvcorrS_ls1.nzval .= 0.
    end

    if ns_liquid_phase
        RuL_ls.nzval .= 0.
        RvL_ls.nzval .= 0.
        
        ∇ϕ_x = opC_uL.AxT * opC_uL.Rx * veci(ϕD1_L,grid,1) .+ opC_uL.Gx * veci(ϕD1_L,grid,2)
        ∇ϕ_y = opC_vL.AyT * opC_vL.Ry * veci(ϕD1_L,grid,1) .+ opC_vL.Gy * veci(ϕD1_L,grid,2)
        iMu = Diagonal(1 ./ (opC_uL.M.diag .+ eps(0.01)))
        iMv = Diagonal(1 ./ (opC_vL.M.diag .+ eps(0.01)))
        GxL = τ .* iMu * ∇ϕ_x
        GyL = τ .* iMv * ∇ϕ_y

        RpL_ls.nzval .= 0.

        if free_surface
            AuL, _, _, AvL, _, _, AϕL = set_navier_stokes(neu, num, grid, grid.geoL, grid_u, grid_u.geoL,      grid_v, grid_v.geoL,
                                                        opC_pL, opC_uL, opC_vL, BC_pL, BC_uL, BC_vL,
                                                        Mum1_L, Mvm1_L, iRe,
                                                        opL, phL,
                                                        periodic_x, periodic_y, ns_advection)[1:7]
        else
            AuL, _, _, AvL, _, _, AϕL = set_navier_stokes(dir, num, grid, grid.geoL, grid_u, grid_u.geoL, grid_v, grid_v.geoL,
                                                        opC_pL, opC_uL, opC_vL, BC_pL, BC_uL, BC_vL,
                                                        Lum1_L, bc_Lum1_L, Lvm1_L, bc_Lvm1_L, Mum1_L, Mvm1_L, iRe,
                                                        opL, phL,
                                                        periodic_x, periodic_y, ns_advection)[1:7]
        end

        SmatL = iRe .* strain_rate(opC_uL, opC_vL)

        GxT = opC_uL.Gx'
        GyT = opC_vL.Gy'
        BϕκL = [spdiagm(fzeros(grid)); σ .* (GxT * opC_uL.Gx .+ GyT * opC_vL.Gy)]

        RucorrL_ls1.nzval .= 0.
        RvcorrL_ls1.nzval .= 0.
    end

    a0_p = zeros(grid)
    _a1_p = zeros(grid)
    _b0_p = ones(grid)
    _b1_p = zeros(grid)
    set_borders!(grid, a0_p, _a1_p, _b0_p, _b1_p, BC_pL, periodic_x, periodic_y)
    b0_p = Diagonal(vec(_b0_p))

    derAϕ = copy(AϕL)
    derAϕ.nzval .= 0.
    derBϕuD1 = copy(opC_pL.AxT)
    derBϕuD1.nzval .= 0.
    derBϕvD1 = copy(opC_pL.AyT)
    derBϕvD1.nzval .= 0.
    derBϕuD2 = copy(opC_pL.Gx)
    derBϕuD2.nzval .= 0.
    derBϕvD2 = copy(opC_pL.Gy)
    derBϕvD2.nzval .= 0.
    derSmat = copy(SmatL)
    derBϕκ = copy(BϕκL)
    derBϕκ.nzval .= 0.

    derAu = copy(AuL)
    derAu.nzval .= 0.
    derAv = copy(AvL)
    derAv.nzval .= 0.
    derMu = copy(opC_uL.M)
    derMv = copy(opC_vL.M)

    # graph coloring
    stencil = 15
    colors, graph, nloops, px0, pxf, py0, pyf = coloring(grid, stencil, periodic_x, periodic_y)
    graph_u = coloring(grid_u, stencil, periodic_x, periodic_y)[2]
    graph_v = coloring(grid_v, stencil, periodic_x, periodic_y)[2]
    @inbounds for color = 0:stencil^2-1
        @inbounds for p = 1:nloops
            indices = findall(colors[py0[p]:pyf[p], px0[p]:pxf[p]] .== color)
            relocate!(indices, px0[p], py0[p])
            assign_color!(grid, stencil, graph, indices, periodic_x, periodic_y)
            assign_color!(grid_u, grid, stencil, graph_u, indices, periodic_x, periodic_y)
            assign_color!(grid_v, grid, stencil, graph_v, indices, periodic_x, periodic_y)

            uj[indices] .+= ϵ_adj

            if ns_solid_phase
                # Compute capacities
                update_ls_data(num, grid, grid_u, grid_v, uj, κj, periodic_x, periodic_y)
                update_free_surface_velocity(num, grid_u, grid_v, uD0_S, vD0_S, periodic_x, periodic_y)


                if free_surface
                    Auj, _, _, Avj, _, _, Aϕj = set_navier_stokes(neu, num, grid, grid.geoS, grid_u, grid_u.geoS, grid_v, grid_v.geoS,
                                                                opC_pj, opC_uj, opC_vj, BC_pS, BC_uS, BC_vS,
                                                                Mum1_S, Mvm1_S, iRe,
                                                                opS, phS,
                                                                periodic_x, periodic_y, ns_advection)[1:7]
                else
                    Auj, _, _, Avj, _, _, Aϕj = set_navier_stokes(dir, num, grid, grid.geoS, grid_u, grid_u.geoS, grid_v, grid_v.geoS,
                                                                opC_pj, opC_uj, opC_vj, BC_pj, BC_uS, BC_vS,
                                                                Lum1_S, bc_Lum1_S, Lvm1_S, bc_Lvm1_S, Mum1_S, Mvm1_S, iRe,
                                                                opS, phS,
                                                                periodic_x, periodic_y, ns_advection)[1:7]
                end

                Smatj = iRe .* strain_rate(opC_uj, opC_vj)

                GxTj = opC_uj.Gx'
                GyTj = opC_vj.Gy'
                Bϕκj = [spdiagm(fzeros(grid)); σ .* (GxTj * opC_uj.Gx .+ GyTj * opC_vj.Gy)]

                # projection step derivatives
                ∇ϕ_x = opC_uj.AxT * opC_uj.Rx * veci(ϕD1_S,grid,1) .+ opC_uj.Gx * veci(ϕD1_S,grid,2)
                ∇ϕ_y = opC_vj.AyT * opC_vj.Ry * veci(ϕD1_S,grid,1) .+ opC_vj.Gy * veci(ϕD1_S,grid,2)
                iMu = Diagonal(1 ./ (opC_uj.M.diag .+ eps(0.01)))
                iMv = Diagonal(1 ./ (opC_vj.M.diag .+ eps(0.01)))
                Gxj = τ .* iMu * ∇ϕ_x
                Gyj = τ .* iMv * ∇ϕ_y
                Ruj = sparse((Gxj .- GxS) ./ ϵ_adj)
                Rvj = sparse((Gyj .- GyS) ./ ϵ_adj)

                rows = rowvals(Ruj)
                for i in nzrange(Ruj, 1)
                    @inbounds row = rows[i]
                    j = graph_u[row]
                    @inbounds RuS_ls[row,j] = Ruj[row]
                end
                rows = rowvals(Rvj)
                for i in nzrange(Rvj, 1)
                    @inbounds row = rows[i]
                    j = graph_v[row]
                    @inbounds RvS_ls[row,j] = Rvj[row]
                end

                # poisson equation derivatives
                derAϕ .= (Aϕj .- AϕS) ./ ϵ_adj
                derBϕuD1 .= (opC_pj.AxT .- opC_pS.AxT) ./ ϵ_adj
                derBϕuD2 .= (opC_pj.Gx .- opC_pS.Gx) ./ ϵ_adj
                derSmat .= (Smatj .- SmatS) ./ ϵ_adj
                derBϕκ .= (Bϕκj .- BϕκS) ./ ϵ_adj
                derBϕvD1 .= (opC_pj.AyT .- opC_pS.AyT) ./ ϵ_adj
                derBϕvD2 .= (opC_pj.Gy .- opC_pS.Gy) ./ ϵ_adj
                derBϕu = [derBϕuD1 derBϕuD2;
                        b0_p*derSmat[1,1] b0_p*derSmat[1,2]]
                derBϕv = [derBϕvD1 derBϕvD2;
                        b0_p*derSmat[2,1] b0_p*derSmat[2,2]]
                derκ = vec(κj .- grid.κ) ./ ϵ_adj

                Rj = sparse(derAϕ * ϕD1_S .-
                            iτ .* derBϕu * ucorrD1_S .-
                            iτ .* derBϕv * vcorrD1_S .+
                            derBϕκ * vec(grid.κ) .+
                            BϕκS * derκ)

                rows = rowvals(Rj)
                for i in nzrange(Rj, 1)
                    @inbounds row = rows[i]
                    j = graph[row]
                    @inbounds RpS_ls[row,j] = Rj[row]
                end

                # prediction step derivatives
                derAu .= (Auj .- AuS) ./ ϵ_adj
                derAv .= (Avj .- AvS) ./ ϵ_adj
                derMu .= (opC_uj.M .- opC_uS.M) ./ ϵ_adj
                derMv .= (opC_vj.M .- opC_vS.M) ./ ϵ_adj

                Ruj = sparse(derAu * ucorrD1_S .- τ.*g.*sin(β) .* derMu * fones(grid_u))
                Rvj = sparse(derAv * vcorrD1_S .+ τ.*g.*cos(β) .* derMv * fones(grid_v))

                rows = rowvals(Ruj)
                for i in nzrange(Ruj, 1)
                    @inbounds row = rows[i]
                    j = graph_u[row]
                    @inbounds RucorrS_ls1[row,j] = Ruj[row]
                end
                rows = rowvals(Rvj)
                for i in nzrange(Rvj, 1)
                    @inbounds row = rows[i]
                    j = graph_v[row]
                    @inbounds RvcorrS_ls1[row,j] = Rvj[row]
                end
            end

            if ns_liquid_phase
                # Compute capacities
                update_ls_data(num, grid, grid_u, grid_v, uj, κj, periodic_x, periodic_y)
                update_free_surface_velocity(num, grid_u, grid_v, uD0_L, vD0_L, periodic_x, periodic_y)

                if free_surface
                    Auj, _, _, Avj, _, _, Aϕj = set_navier_stokes(neu, num, grid, grid.geoL, grid_u, grid_u.geoL,      grid_v, grid_v.geoL,
                                                                opC_pj, opC_uj, opC_vj, BC_pL, BC_uL, BC_vL,
                                                                Mum1_L, Mvm1_L, iRe,
                                                                opL, phL,
                                                                periodic_x, periodic_y, ns_advection)[1:7]
                else
                    Auj, _, _, Avj, _, _, Aϕj = set_navier_stokes(dir, num, grid, grid.geoL, grid_u, grid_u.geoL, grid_v, grid_v.geoL,
                                                                opC_pj, opC_uj, opC_vj, BC_pL, BC_uL, BC_vL,
                                                                Lum1_L, bc_Lum1_L, Lvm1_L, bc_Lvm1_L, Mum1_L, Mvm1_L, iRe,
                                                                opL, phL,
                                                                periodic_x, periodic_y, ns_advection)[1:7]
                end

                Smatj = iRe .* strain_rate(opC_uj, opC_vj)

                GxTj = opC_uj.Gx'
                GyTj = opC_vj.Gy'
                Bϕκj = [spdiagm(fzeros(grid)); σ .* (GxTj * opC_uj.Gx .+ GyTj * opC_vj.Gy)]

                # projection step derivatives
                ∇ϕ_x = opC_uj.AxT * opC_uj.Rx * veci(ϕD1_L,grid,1) .+ opC_uj.Gx * veci(ϕD1_L,grid,2)
                ∇ϕ_y = opC_vj.AyT * opC_vj.Ry * veci(ϕD1_L,grid,1) .+ opC_vj.Gy * veci(ϕD1_L,grid,2)
                iMu = Diagonal(1 ./ (opC_uj.M.diag .+ eps(0.01)))
                iMv = Diagonal(1 ./ (opC_vj.M.diag .+ eps(0.01)))
                Gxj = τ .* iMu * ∇ϕ_x
                Gyj = τ .* iMv * ∇ϕ_y
                Ruj = sparse((Gxj .- GxL) ./ ϵ_adj)
                Rvj = sparse((Gyj .- GyL) ./ ϵ_adj)

                rows = rowvals(Ruj)
                for i in nzrange(Ruj, 1)
                    @inbounds row = rows[i]
                    j = graph_u[row]
                    @inbounds RuL_ls[row,j] = Ruj[row]
                end
                rows = rowvals(Rvj)
                for i in nzrange(Rvj, 1)
                    @inbounds row = rows[i]
                    j = graph_v[row]
                    @inbounds RvL_ls[row,j] = Rvj[row]
                end

                # poisson equation derivatives
                derAϕ .= (Aϕj .- AϕL) ./ ϵ_adj
                derBϕuD1 .= (opC_pj.AxT .- opC_pL.AxT) ./ ϵ_adj
                derBϕuD2 .= (opC_pj.Gx .- opC_pL.Gx) ./ ϵ_adj
                derSmat .= (Smatj .- SmatL) ./ ϵ_adj
                derBϕκ .= (Bϕκj .- BϕκL) ./ ϵ_adj
                derBϕvD1 .= (opC_pj.AyT .- opC_pL.AyT) ./ ϵ_adj
                derBϕvD2 .= (opC_pj.Gy .- opC_pL.Gy) ./ ϵ_adj
                derBϕu = [derBϕuD1 derBϕuD2;
                        b0_p*derSmat[1,1] b0_p*derSmat[1,2]]
                derBϕv = [derBϕvD1 derBϕvD2;
                        b0_p*derSmat[2,1] b0_p*derSmat[2,2]]
                derκ = vec(κj .- grid.κ) ./ ϵ_adj

                Rj = sparse(derAϕ * ϕD1_L .-
                            iτ .* derBϕu * ucorrD1_L .-
                            iτ .* derBϕv * vcorrD1_L .+
                            derBϕκ * vec(grid.κ) .+
                            BϕκL * derκ)

                rows = rowvals(Rj)
                for i in nzrange(Rj, 1)
                    @inbounds row = rows[i]
                    j = graph[row]
                    @inbounds RpL_ls[row,j] = Rj[row]
                end

                # prediction step derivatives
                derAu .= (Auj .- AuL) ./ ϵ_adj
                derAv .= (Avj .- AvL) ./ ϵ_adj
                derMu .= (opC_uj.M .- opC_uL.M) ./ ϵ_adj
                derMv .= (opC_vj.M .- opC_vL.M) ./ ϵ_adj

                Ruj = sparse(derAu * ucorrD1_L .-
                    vcat(τ.*g.*sin(β) .* derMu * fones(grid_u), fzeros(grid_u)))
                Rvj = sparse(derAv * vcorrD1_L .+
                    vcat(τ.*g.*cos(β) .* derMv * fones(grid_v), fzeros(grid_v)))

                rows = rowvals(Ruj)
                for i in nzrange(Ruj, 1)
                    @inbounds row = rows[i]
                    j = graph_u[row]
                    @inbounds RucorrL_ls1[row,j] = Ruj[row]
                end
                rows = rowvals(Rvj)
                for i in nzrange(Rvj, 1)
                    @inbounds row = rows[i]
                    j = graph_v[row]
                    @inbounds RvcorrL_ls1[row,j] = Rvj[row]
                end
            end

            uj .= u1
        end
    end

    return nothing
end

function Rproj_q0(num, grid, grid_u, grid_v, adj_der,
    opC_pS, BC_pS, BC_pL,
    uD0_S, opC_uS, BC_uS,
    uD0_L, BC_uL,
    vD0_S, opC_vS, BC_vS,
    vD0_L, BC_vL,
    BuSm1, BuLm1, BvSm1, BvLm1,
    Mum1_S, Mum1_L, Mvm1_S, Mvm1_L,
    u0, u1, LSA, LSB,
    periodic_x, periodic_y, ϵ_adj,
    ns_solid_phase,
    ns_liquid_phase)

    @unpack Re, τ, σ, g, β, ϵ, NB = num
    @unpack nx, ny, ind, V, iso, faces, geoS, geoL = grid
    @unpack RucorrS_ls0, RvcorrS_ls0, RucorrL_ls0, RvcorrL_ls0,
            RlsFS_ls, RlsFS_ucorrS, RlsFS_ucorrL,
            RlsFS_vcorrS, RlsFS_vcorrL = adj_der

    iRe = 1.0 / Re

    κj = copy(grid.κ)
    uj = copy(u0)
    
    uD0_Sj = copy(uD0_S)
    vD0_Sj = copy(vD0_S)
    uD0_Lj = copy(uD0_L)
    vD0_Lj = copy(vD0_L)

    opC_pj = copy(opC_pS)
    opC_uj = copy(opC_uS)
    opC_vj = copy(opC_vS)

    θ_out = zeros(grid,4)

    # prediction step derivatives
    RucorrS_ls0.nzval .= 0.
    RvcorrS_ls0.nzval .= 0.

    RucorrL_ls0.nzval .= 0.
    RvcorrL_ls0.nzval .= 0.

    derBu = copy(BuLm1)
    derBu.nzval .= 0.
    derBv = copy(BvLm1)
    derBv.nzval .= 0.

    # levelset advection derivatives
    RlsFS_ls.nzval .= 0.
    RlsFS_ucorrS .= 0.
    RlsFS_ucorrL .= 0.
    RlsFS_vcorrL .= 0.
    RlsFS_vcorrS .= 0.

    derLSA = copy(LSA)
    derLSA.nzval .= 0.
    derLSB = copy(LSB)
    derLSB.nzval .= 0.
    LSAj = copy(LSA)
    LSBj = copy(LSB)

    # graph coloring
    stencil = 15
    colors, graph, nloops, px0, pxf, py0, pyf = coloring(grid, stencil, periodic_x, periodic_y)
    colors_u, graph_u, _, px0_u, pxf_u, py0_u, pyf_u = coloring(grid_u, stencil, periodic_x, periodic_y)
    colors_v, graph_v, _, px0_v, pxf_v, py0_v, pyf_v = coloring(grid_v, stencil, periodic_x, periodic_y)
    graph_ls_u = coloring(grid, stencil, periodic_x, periodic_y)[2]
    graph_ls_v = coloring(grid, stencil, periodic_x, periodic_y)[2]
    # @inbounds for color = 0:stencil^2-1
    #     @inbounds for p = 1:nloops
        for (index_u, index_v) in zip(grid_u.ind.all_indices, grid_v.ind.all_indices)
            # Perturbations of the levelset
            # indices = findall(colors[py0[p]:pyf[p], px0[p]:pxf[p]] .== color)
            # relocate!(indices, px0[p], py0[p])
            # assign_color!(grid, stencil, graph, indices, periodic_x, periodic_y)
            # assign_color!(grid_u, grid, stencil, graph_u, indices, periodic_x, periodic_y)
            # assign_color!(grid_v, grid, stencil, graph_v, indices, periodic_x, periodic_y)

            # uj[indices] .+= ϵ_adj
            if index_u[2] <= grid.nx
                j = lexicographic(index_u, grid.ny)
                uj[index_u] += ϵ_adj
            end

            # Perturbations of velocity components
            # indices_u = findall(colors_u[py0_u[p]:pyf_u[p], px0_u[p]:pxf_u[p]] .== color)
            # relocate!(indices_u, px0_u[p], py0_u[p])
            # js_u = lexicographic.(indices_u, grid_u.ny)
            # assign_color!(grid, grid_u, stencil, graph_ls_u, indices_u, periodic_x, periodic_y)

            # veci(uD0_Sj,grid_u,2)[js_u] .+= ϵ_adj
            # veci(uD0_Lj,grid_u,2)[js_u] .+= ϵ_adj
            j_u = lexicographic(index_u, grid_u.ny)
            veci(uD0_Sj,grid_u,2)[j_u] += ϵ_adj
            veci(uD0_Lj,grid_u,2)[j_u] += ϵ_adj

            # indices_v = findall(colors_v[py0_v[p]:pyf_v[p], px0_v[p]:pxf_v[p]] .== color)
            # relocate!(indices_v, px0_v[p], py0_v[p])
            # js_v = lexicographic.(indices_v, grid_v.ny)
            # assign_color!(grid, grid_v, stencil, graph_ls_v, indices_v, periodic_x, periodic_y)

            # veci(vD0_Sj,grid_v,2)[js_v] .+= ϵ_adj
            # veci(vD0_Lj,grid_v,2)[js_v] .+= ϵ_adj
            j_v = lexicographic(index_v, grid_v.ny)
            veci(vD0_Sj,grid_v,2)[j_v] += ϵ_adj
            veci(vD0_Lj,grid_v,2)[j_v] += ϵ_adj


            if ns_solid_phase
                # Compute capacities
                update_ls_data(num, grid, grid_u, grid_v, uj, κj, periodic_x, periodic_y)
                update_free_surface_velocity(num, grid_u, grid_v, uD0_S, vD0_S, periodic_x, periodic_y)

                Lp, bc_Lp, Lu, bc_Lu, Lv, bc_Lv = set_matrices!(grid, grid.geoS, grid_u, grid_u.geoS, grid_v, grid_v.geoS,
                    opC_pj, opC_uj, opC_vj,
                    periodic_x, periodic_y)
            
                _, Buj, _, _, Bvj, _, _, _ = set_crank_nicolson_block(neu, num,
                    grid, opC_pj, Lp, bc_Lp, Lp_fs, bc_Lp_fs, BC_pS,
                    grid_u, opC_uj, iRe.*Lu, iRe.*bc_Lu, Mum1_S, BC_uS,
                    grid_v, opC_vj, iRe.*Lv, iRe.*bc_Lv, Mvm1_S, BC_vS,
                    periodic_x, periodic_y)

                # prediction step derivatives
                derBu .= (Buj .- BuSm1) ./ ϵ_adj
                derBv .= (Bvj .- BvSm1) ./ ϵ_adj

                Ruj = sparse(- derBu * uD0_S)
                Rvj = sparse(- derBv * vD0_S)

                rows = rowvals(Ruj)
                for i in nzrange(Ruj, 1)
                    @inbounds row = rows[i]
                    j = graph_u[row]
                    @inbounds RucorrS_ls0[row,j] = Ruj[row]
                end
                rows = rowvals(Rvj)
                for i in nzrange(Rvj, 1)
                    @inbounds row = rows[i]
                    j = graph_v[row]
                    @inbounds RvcorrS_ls0[row,j] = Rvj[row]
                end

                # levelset advection derivatives
                update_ls_data(num, grid, grid_u, grid_v, u0, grid.κ, periodic_x, periodic_y)
                update_free_surface_velocity(num, grid_u, grid_v, uD0_Sj, vD0_S, periodic_x, periodic_y)
                IIOE!(grid, grid_u, grid_v, LSAj, LSBj, θ_out, τ, periodic_x, periodic_y)
                utmp = reshape(gmres(LSAj,(LSBj*vec(u0))), (ny,nx))
                S2IIOE!(grid, grid_u, grid_v, LSAj, LSBj, utmp, u0, θ_out, τ, periodic_x, periodic_y)

                derLSA .= (LSAj .- LSA) ./ ϵ_adj
                derLSB .= (LSBj .- LSB) ./ ϵ_adj
                Rj = sparse(derLSA * vec(u1) .- derLSB * vec(u0))

                rows = rowvals(Rj)
                for i in nzrange(Rj, 1)
                    @inbounds row = rows[i]
                    j = graph_ls_u[row] + grid_u.ny*grid_u.nx
                    @inbounds RlsFS_ucorrS[row,j] = Rj[row]
                end

                update_ls_data(num, grid, grid_u, grid_v, u0, grid.κ, periodic_x, periodic_y)
                update_free_surface_velocity(num, grid_u, grid_v, uD0_S, vD0_Sj, periodic_x, periodic_y)
                IIOE!(grid, grid_u, grid_v, LSAj, LSBj, θ_out, τ, periodic_x, periodic_y)
                utmp = reshape(gmres(LSAj,(LSBj*vec(u0))), (ny,nx))
                S2IIOE!(grid, grid_u, grid_v, LSAj, LSBj, utmp, u0, θ_out, τ, periodic_x, periodic_y)

                derLSA .= (LSAj .- LSA) ./ ϵ_adj
                derLSB .= (LSBj .- LSB) ./ ϵ_adj
                Rj = sparse(derLSA * vec(u1) .- derLSB * vec(u0))

                rows = rowvals(Rj)
                for i in nzrange(Rj, 1)
                    @inbounds row = rows[i]
                    j = graph_ls_v[row] + grid_v.ny*grid_v.nx
                    @inbounds RlsFS_vcorrS[row,j] = Rj[row]
                end
            end

            if ns_liquid_phase
                # Compute capacities
                update_ls_data(num, grid, grid_u, grid_v, uj, κj, periodic_x, periodic_y)
                update_free_surface_velocity(num, grid_u, grid_v, uD0_L, vD0_L, periodic_x, periodic_y)

                Lp, bc_Lp, Lu, bc_Lu, Lv, bc_Lv, Lp_fs, bc_Lp_fs = set_matrices!(grid, grid.geoL, grid_u, grid_u.geoL, grid_v, grid_v.geoL,
                    opC_pj, opC_uj, opC_vj,
                    periodic_x, periodic_y)
            
                _, Buj, _, _, Bvj, _, _, _ = set_crank_nicolson_block(neu, num,
                    grid, opC_pj, Lp, bc_Lp, Lp_fs, bc_Lp_fs, BC_pL,
                    grid_u, opC_uj, iRe.*Lu, iRe.*bc_Lu, Mum1_L, BC_uL,
                    grid_v, opC_vj, iRe.*Lv, iRe.*bc_Lv, Mvm1_L, BC_vL,
                    periodic_x, periodic_y)

                # prediction step derivatives
                derBu .= (Buj .- BuLm1) ./ ϵ_adj
                derBv .= (Bvj .- BvLm1) ./ ϵ_adj

                Ruj = sparse(- derBu * uD0_L)
                Rvj = sparse(- derBv * vD0_L)

                if index_u[2] <= grid.nx
                    rows = rowvals(Ruj)
                    for i in nzrange(Ruj, 1)
                        @inbounds row = rows[i]
                        # j = graph_u[row]
                        @inbounds RucorrL_ls0[row,j] = Ruj[row]
                    end
                    rows = rowvals(Rvj)
                    for i in nzrange(Rvj, 1)
                        @inbounds row = rows[i]
                        # j = graph_v[row]
                        @inbounds RvcorrL_ls0[row,j] = Rvj[row]
                    end
                end

                # levelset advection derivatives
                IIOE!(grid, grid_u, grid_v, LSAj, LSBj, θ_out, τ, periodic_x, periodic_y)
                utmp = reshape(gmres(LSAj,(LSBj*vec(uj))), (ny,nx))
                S2IIOE!(grid, grid_u, grid_v, LSAj, LSBj, utmp, uj, θ_out, τ, periodic_x, periodic_y)

                derLSA .= (LSAj .- LSA) ./ ϵ_adj
                derLSB .= (LSBj .- LSB) ./ ϵ_adj

                utmp = zeros(grid)
                # utmp[indices] .= 1.0
                if index_u[2] <= grid.nx
                    utmp[index_u] = 1.0
                end
                Rj = sparse(derLSA * vec(u1) .- derLSB * vec(u0) .- LSB * vec(utmp))

                # We should get exactly the same with both solid and liquid phases,
                # since the jump of velocity fields at the interface should be 0
                if index_u[2] <= grid.nx
                    rows = rowvals(Rj)
                    for i in nzrange(Rj, 1)
                        @inbounds row = rows[i]
                        # j = graph[row]
                        @inbounds RlsFS_ls[row,j] = Rj[row]
                    end
                end

                update_ls_data(num, grid, grid_u, grid_v, u0, grid.κ, periodic_x, periodic_y)
                update_free_surface_velocity(num, grid_u, grid_v, uD0_Lj, vD0_L, periodic_x, periodic_y)
                IIOE!(grid, grid_u, grid_v, LSAj, LSBj, θ_out, τ, periodic_x, periodic_y)
                utmp = reshape(gmres(LSAj,(LSBj*vec(u0))), (ny,nx))
                S2IIOE!(grid, grid_u, grid_v, LSAj, LSBj, utmp, u0, θ_out, τ, periodic_x, periodic_y)

                derLSA .= (LSAj .- LSA) ./ ϵ_adj
                derLSB .= (LSBj .- LSB) ./ ϵ_adj
                Rj = sparse(derLSA * vec(u1) .- derLSB * vec(u0))

                rows = rowvals(Rj)
                for i in nzrange(Rj, 1)
                    @inbounds row = rows[i]
                    # j = graph_ls_u[row] + grid_u.ny*grid_u.nx
                    @inbounds RlsFS_ucorrL[row,j_u + grid_u.ny*grid_u.nx] = Rj[row]
                end

                update_ls_data(num, grid, grid_u, grid_v, u0, grid.κ, periodic_x, periodic_y)
                update_free_surface_velocity(num, grid_u, grid_v, uD0_L, vD0_Lj, periodic_x, periodic_y)
                IIOE!(grid, grid_u, grid_v, LSAj, LSBj, θ_out, τ, periodic_x, periodic_y)
                utmp = reshape(gmres(LSAj,(LSBj*vec(u0))), (ny,nx))
                S2IIOE!(grid, grid_u, grid_v, LSAj, LSBj, utmp, u0, θ_out, τ, periodic_x, periodic_y)

                derLSA .= (LSAj .- LSA) ./ ϵ_adj
                derLSB .= (LSBj .- LSB) ./ ϵ_adj
                Rj = sparse(derLSA * vec(u1) .- derLSB * vec(u0))

                rows = rowvals(Rj)
                for i in nzrange(Rj, 1)
                    @inbounds row = rows[i]
                    # j = graph_ls_v[row] + grid_v.ny*grid_v.nx
                    @inbounds RlsFS_vcorrL[row,j_v + grid_v.ny*grid_v.nx] = Rj[row]
                end
            end

            uj .= u0
            uD0_Sj .= uD0_S
            vD0_Sj .= vD0_S
            uD0_Lj .= uD0_L
            vD0_Lj .= vD0_L
        end
    # end

    return nothing
end