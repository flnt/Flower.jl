function coloring(grid, stencil, periodic_x, periodic_y)
    @unpack nx, ny, ind = grid
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
        pxf[1] = nx-3
        pyf[1] = ny-3
        px0[2] = nx-2
        pyf[2] = ny-3
        pxf[3] = nx-3
        py0[3] = ny-2
        px0[4] = nx-2
        py0[4] = ny-2
    elseif periodic_x
        nloops = 2
        pxf[1] = nx-3
        px0[2] = nx-2
    else
        nloops = 2
        pyf[1] = ny-3
        py0[2] = ny-2
    end

    return colors, graph, nloops, px0, pxf, py0, pyf
end

function assign_color!(grid, graph, ind_pert, per_x, per_y)
    @unpack nx, ny, ind = grid
    @unpack all_indices = ind

    graph .= 0
    if !per_x && !per_y
        @inbounds for II in all_indices
            i = lexicographic(II, ny)
            @inbounds for JJ in ind_pert
                j = lexicographic(JJ, ny)
                if abs(II[1] - JJ[1]) <= 3 && abs(II[2] - JJ[2]) <= 3
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
                if (abs(II[1] - JJ[1]) <= 3 || abs(II[1] - JJ[1]) >= ny-3) && abs(II[2] - JJ[2]) <= 3
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
                if abs(II[1] - JJ[1]) <= 3 && (abs(II[2] - JJ[2]) <= 3 || abs(II[2] - JJ[2]) >= nx-3)
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
                if (abs(II[1] - JJ[1]) <= 3 || abs(II[1] - JJ[1]) >= ny-3) && (abs(II[2] - JJ[2]) <= 3 || abs(II[2] - JJ[2]) >= nx-3)
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
    heat_solid_phase, heat_liquid_phase)

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
            assign_color!(grid, graph, indices, periodic_x, periodic_y)

            uj[indices] .+= ϵ_adj

            # Compute capacities
            update_ls_data(num, grid, grid_u, grid_v, uj, κ, periodic_x, periodic_y)
            update_stefan_velocity(num, grid, uj, TS, TL, periodic_x, periodic_y, λ, Vmean)

            # get perturbed matrices
            if heat_solid_phase
                Aj, Bj, _ = set_heat!(dir, num, grid, opC_TS, geoS, BC_TS, grid.ind.MIXED, geoS.projection,
                                        periodic_x, periodic_y)
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
                Aj, Bj, _ = set_heat!(dir, num, grid, opC_TL, geoL, BC_TL, grid.ind.MIXED, geoL.projection,
                                        periodic_x, periodic_y)
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

            IIOE(grid, LSAj, LSBj, uj, V, CFL_sc, periodic_x, periodic_y)
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
    @unpack nx, ny, ind, u, V, faces, iso, geoS, geoL = grid
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
            assign_color!(grid, graph, indices, periodic_x, periodic_y)

            TSj[js] .+= ϵ_adj
            TLj[js] .+= ϵ_adj

            # Compute capacities
            update_ls_data(num, grid, grid_u, grid_v, u, κ, periodic_x, periodic_y)

            # get perturbed matrices
            if heat_solid_phase
                update_stefan_velocity(num, grid, u, TSj, TL, periodic_x, periodic_y, λ, Vmean)
                IIOE(grid, LSAj, LSBj, u, V, CFL_sc, periodic_x, periodic_y)
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
                update_stefan_velocity(num, grid, u, TS, TLj, periodic_x, periodic_y, λ, Vmean)
                IIOE(grid, LSAj, LSBj, u, V, CFL_sc, periodic_x, periodic_y)
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
    MuSm1, MuLm1, MvSm1, MvLm1,
    u1, periodic_x, periodic_y, ϵ_adj)

    @unpack Re, τ, σ, g, β, ϵ, NB = num
    @unpack nx, ny, ind, V, iso, faces, geoS, geoL = grid
    @unpack RuS_ls, RuL_ls, RuS_ls, RuL_ls, RpS_ls, RpL_ls,
            RucorrS_ls0, RucorrS_ls1, RvcorrS_ls0, RvcorrS_ls1,
            RucorrL_ls0, RucorrL_ls1, RvcorrL_ls0, RvcorrL_ls1 = adj_der

    iRe = 1.0 / Re
    iτ = 1.0 / τ

    κj = copy(grid.κ)
    uj = copy(u1)

    opC_pj = copy(opC_pS)
    opC_uj = copy(opC_uS)
    opC_vj = copy(opC_vS)

    # projection step derivatives
    RuS_ls.nzval .= 0.
    RvS_ls.nzval .= 0.
    
    ∇ϕ_x = opC_uS.AxT * opC_uS.Rx * veci(ϕD1_S,grid,1) .+ opC_uS.Gx * veci(ϕD1_S,grid,2)
    ∇ϕ_y = opC_vS.AyT * opC_vS.Ry * veci(ϕD1_S,grid,1) .+ opC_vS.Gy * veci(ϕD1_S,grid,2)
    iMu = Diagonal(1 ./ (opC_uS.M.diag .+ eps(0.01)))
    iMv = Diagonal(1 ./ (opC_vS.M.diag .+ eps(0.01)))
    GxS = τ .* iMu * ∇ϕ_x
    GyS = τ .* iMv * ∇ϕ_y

    RuL_ls.nzval .= 0.
    RvL_ls.nzval .= 0.
    
    ∇ϕ_x = opC_uL.AxT * opC_uL.Rx * veci(ϕD1_L,grid,1) .+ opC_uL.Gx * veci(ϕD1_L,grid,2)
    ∇ϕ_y = opC_vL.AyT * opC_vL.Ry * veci(ϕD1_L,grid,1) .+ opC_vL.Gy * veci(ϕD1_L,grid,2)
    iMu = Diagonal(1 ./ (opC_uL.M.diag .+ eps(0.01)))
    iMv = Diagonal(1 ./ (opC_vL.M.diag .+ eps(0.01)))
    GxL = τ .* iMu * ∇ϕ_x
    GyL = τ .* iMv * ∇ϕ_y

    # poisson equation derivatives
    RpS_ls.nzval .= 0.

    Lp, bc_Lp, Lu, bc_Lu, Lv, bc_Lv, Lp_fs, bc_Lp_fs, Lu_fs, bc_Lu_fs, Lv_fs, bc_Lv_fs = set_laplacians!(grid, grid.geoS, grid_u, grid_u.geoS, grid_v, grid_v.geoS,
        opC_pS, opC_uS, opC_vS,
        periodic_x, periodic_y, true)

    AuS, _, _, AvS, _, _, AϕS, _ = set_crank_nicolson_block(neu, num,
        grid, opC_pS, Lp, bc_Lp, Lp_fs, bc_Lp_fs, BC_pS,
        grid_u, opC_uS, iRe.*Lu, iRe.*bc_Lu, iRe.*Lum1, iRe.*bc_Lum1, MuSm1, BC_uS,
        grid_v, opC_vS, iRe.*Lv, iRe.*bc_Lv, iRe.*Lvm1, iRe.*bc_Lvm1, MvSm1, BC_vS)

    SmatS = iRe .* strain_rate(opC_uS, opC_vS)

    GxT = opC_uS.Gx'
    GyT = opC_vS.Gy'
    BϕκS = [spdiagm(fones(grid)); σ .* (GxT * opC_uS.Gx .+ GyT * opC_vS.Gy)]

    RpL_ls.nzval .= 0.

    Lp, bc_Lp, Lu, bc_Lu, Lv, bc_Lv, Lp_fs, bc_Lp_fs, Lu_fs, bc_Lu_fs, Lv_fs, bc_Lv_fs = set_laplacians!(grid, grid.geoL, grid_u, grid_u.geoL, grid_v, grid_v.geoL,
        opC_pL, opC_uL, opC_vL,
        periodic_x, periodic_y, true)

    AuL, _, _, AvL, _, _, AϕL, _ = set_crank_nicolson_block(neu, num,
        grid, opC_pL, Lp, bc_Lp, Lp_fs, bc_Lp_fs, BC_pL,
        grid_u, opC_uL, iRe.*Lu, iRe.*bc_Lu, iRe.*Lum1, iRe.*bc_Lum1, MuLm1, BC_uL,
        grid_v, opC_vL, iRe.*Lv, iRe.*bc_Lv, iRe.*Lvm1, iRe.*bc_Lvm1, MvLm1, BC_vL)

    SmatL = iRe .* strain_rate(opC_uL, opC_vL)

    GxT = opC_uL.Gx'
    GyT = opC_vL.Gy'
    BϕκL = [spdiagm(fones(grid)); σ .* (GxT * opC_uL.Gx .+ GyT * opC_vL.Gy)]

    a0_p = zeros(grid)
    _a1_p = zeros(grid)
    _b0_p = ones(grid)
    _b1_p = zeros(grid)
    set_borders!(grid, a0_p, _a1_p, _b0_p, _b1_p, BC_p)
    b0_p = Diagonal(vec(_b0_p))

    derAϕ = copy(AϕS)
    derAϕ.nzval .= 0.
    derBϕuD1 = copy(opC_pS.AxT)
    derBϕuD1.nzval .= 0.
    derBϕvD1 = copy(opC_pS.AyT)
    derBϕvD1.nzval .= 0.
    derBϕuD2 = copy(opC_pS.Gx)
    derBϕuD2.nzval .= 0.
    derBϕvD2 = copy(opC_pS.Gy)
    derBϕvD2.nzval .= 0.
    derSmat = copy(SmatS)
    derBϕκ = copy(Bϕκ)
    derBϕκ.nzval = 0.

    # prediction step derivatives
    RucorrS_ls1.nzval .= 0.
    RvcorrS_ls1.nzval .= 0.

    RucorrL_ls1.nzval .= 0.
    RvcorrL_ls1.nzval .= 0.

    derAu = copy(AuS)
    derAu.nzval .= 0.
    derAv = copy(AvS)
    derAv.nzval .= 0.
    derMu = copy(opC_uS.M)
    derMv = copy(opC_vS.M)

    # graph coloring
    stencil = 7
    colors, graph, nloops, px0, pxf, py0, pyf = coloring(grid, stencil, periodic_x, periodic_y)
    @inbounds for color = 0:stencil^2-1
        @inbounds for p = 1:nloops
            indices = findall(colors[py0[p]:pyf[p], px0[p]:pxf[p]] .== color)
            relocate!(indices, px0[p], py0[p])
            assign_color!(grid, graph, indices, periodic_x, periodic_y)

            uj[indices] .+= ϵ_adj

            # SOLID PHASE
            # Compute capacities
            update_ls_data(num, grid, grid_u, grid_v, uj, κj, periodic_x, periodic_y)
            update_free_surface_velocity(num, grid_u, grid_v, uD0_S, vD0_S, periodic_x, periodic_y)

            Lp, bc_Lp, Lu, bc_Lu, Lv, bc_Lv, Lp_fs, bc_Lp_fs, Lu_fs, bc_Lu_fs, Lv_fs, bc_Lv_fs = set_laplacians!(grid, grid.geoS, grid_u, grid_u.geoS, grid_v, grid_v.geoS,
                opC_pj, opC_uj, opC_vj,
                periodic_x, periodic_y, true)
        
            Auj, _, _, Avj, _, _, Aϕj, _ = set_crank_nicolson_block(neu, num,
                grid, opC_pj, Lp, bc_Lp, Lp_fs, bc_Lp_fs, BC_pS,
                grid_u, opC_uj, iRe.*Lu, iRe.*bc_Lu, iRe.*Lum1, iRe.*bc_Lum1, MuSm1, BC_uS,
                grid_v, opC_vj, iRe.*Lv, iRe.*bc_Lv, iRe.*Lvm1, iRe.*bc_Lvm1, MvSm1, BC_vS)

            Smatj = iRe .* strain_rate(opC_uj, opC_vj)

            GxTj = opC_uS.Gx'
            GyTj = opC_vS.Gy'
            Bϕκj = [spdiagm(fones(grid)); σ .* (GxTj * opC_uj.Gx .+ GyTj * opC_vj.Gy)]

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
                j = graph[row]
                @inbounds RuS_ls[row,j] = Ruj[row]
            end
            rows = rowvals(Rvj)
            for i in nzrange(Rvj, 1)
                @inbounds row = rows[i]
                j = graph[row]
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
            derκ = vcat(fzeros(grid), vec(κj .- grid.κ) ./ ϵ_adj)

            Rj = sparse(derAϕ * ϕD1_S .-
                        iτ .* derBϕu * ucorrD1_S .-
                        iτ .* derBϕv * vcorrD1_S .+
                        derBϕκ * vcat(fzeros(grid), vec(κ)) .+
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
                j = graph[row]
                @inbounds RucorrS_ls1[row,j] = Ruj[row]
            end
            rows = rowvals(Rvj)
            for i in nzrange(Rvj, 1)
                @inbounds row = rows[i]
                j = graph[row]
                @inbounds RvcorrS_ls1[row,j] = Rvj[row]
            end

            # LIQUID PHASE
            # Compute capacities
            update_ls_data(num, grid, grid_u, grid_v, uj, κj, periodic_x, periodic_y)
            update_free_surface_velocity(num, grid_u, grid_v, uD0_L, vD0_L, periodic_x, periodic_y)

            Lp, bc_Lp, Lu, bc_Lu, Lv, bc_Lv, Lp_fs, bc_Lp_fs, Lu_fs, bc_Lu_fs, Lv_fs, bc_Lv_fs = set_laplacians!(grid, grid.geoL, grid_u, grid_u.geoL, grid_v, grid_v.geoL,
                opC_pj, opC_uj, opC_vj,
                periodic_x, periodic_y, true)
        
            Auj, _, _, Avj, _, _, Aϕj, _ = set_crank_nicolson_block(neu, num,
                grid, opC_pj, Lp, bc_Lp, Lp_fs, bc_Lp_fs, BC_pL,
                grid_u, opC_uj, iRe.*Lu, iRe.*bc_Lu, iRe.*Lum1, iRe.*bc_Lum1, MuLm1, BC_uL,
                grid_v, opC_vj, iRe.*Lv, iRe.*bc_Lv, iRe.*Lvm1, iRe.*bc_Lvm1, MvLm1, BC_vL)

            Smatj = iRe .* strain_rate(opC_uj, opC_vj)

            GxTj = opC_uL.Gx'
            GyTj = opC_vL.Gy'
            Bϕκj = [spdiagm(fones(grid)); σ .* (GxTj * opC_uj.Gx .+ GyTj * opC_vj.Gy)]

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
                j = graph[row]
                @inbounds RuL_ls[row,j] = Ruj[row]
            end
            rows = rowvals(Rvj)
            for i in nzrange(Rvj, 1)
                @inbounds row = rows[i]
                j = graph[row]
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
            derκ = vcat(fzeros(grid), vec(κj .- grid.κ) ./ ϵ_adj)

            Rj = sparse(derAϕ * ϕD1_L .-
                        iτ .* derBϕu * ucorrD1_L .-
                        iτ .* derBϕv * vcorrD1_L .+
                        derBϕκ * vcat(fzeros(grid), vec(κ)) .+
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

            Ruj = sparse(derAu * ucorrD1_L .- τ.*g.*sin(β) .* derMu * fones(grid_u))
            Rvj = sparse(derAv * vcorrD1_L .+ τ.*g.*cos(β) .* derMv * fones(grid_v))

            rows = rowvals(Ruj)
            for i in nzrange(Ruj, 1)
                @inbounds row = rows[i]
                j = graph[row]
                @inbounds RucorrL_ls1[row,j] = Ruj[row]
            end
            rows = rowvals(Rvj)
            for i in nzrange(Rvj, 1)
                @inbounds row = rows[i]
                j = graph[row]
                @inbounds RvcorrL_ls1[row,j] = Rvj[row]
            end

            uj .= u1
        end
    end

    return nothing
end

function Rproj_q0(num, grid, grid_u, grid_v, adj_der,
    ϕD1_S, opC_pS, BC_pS,
    ϕD0_L, ϕD1_L, p0_L, p1_L, opC_pL, BC_pL,
    uD0_S, opC_uS, BC_uS,
    uD0_L, uD1_L, opC_uL, BC_uL,
    vD0_S, opC_vS, BC_vS,
    vD0_L, vD1_L, opC_vL, BC_vL,
    MuSm1, MuLm1, MvSm1, MvLm1,
    u0, u1, LSA, LSB,
    periodic_x, periodic_y, ϵ_adj)

    @unpack Re, τ, σ, g, β, ϵ, NB = num
    @unpack nx, ny, ind, V, iso, faces, geoS, geoL = grid
    @unpack RuS_ls, RuL_ls, RuS_ls, RuL_ls, RpS_ls, RpL_ls,
            RucorrS_ls0, RucorrS_ls1, RvcorrS_ls0, RvcorrS_ls1,
            RucorrL_ls0, RucorrL_ls1, RvcorrL_ls0, RvcorrL_ls1 = adj_der

    iRe = 1.0 / Re
    iτ = 1.0 / τ

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

    Lp, bc_Lp, Lu, bc_Lu, Lv, bc_Lv, Lp_fs, bc_Lp_fs, Lu_fs, bc_Lu_fs, Lv_fs, bc_Lv_fs = set_laplacians!(grid, grid.geoS, grid_u, grid_u.geoS, grid_v, grid_v.geoS,
        opC_pS, opC_uS, opC_vS,
        periodic_x, periodic_y, true)

    _, BuS, _, _, BvS, _, _, _ = set_crank_nicolson_block(neu, num,
        grid, opC_pS, Lp, bc_Lp, Lp_fs, bc_Lp_fs, BC_pS,
        grid_u, opC_uS, iRe.*Lu, iRe.*bc_Lu, iRe.*Lum1, iRe.*bc_Lum1, MuSm1, BC_uS,
        grid_v, opC_vS, iRe.*Lv, iRe.*bc_Lv, iRe.*Lvm1, iRe.*bc_Lvm1, MvSm1, BC_vS)

    
    Lp, bc_Lp, Lu, bc_Lu, Lv, bc_Lv, Lp_fs, bc_Lp_fs, Lu_fs, bc_Lu_fs, Lv_fs, bc_Lv_fs = set_laplacians!(grid, grid.geoL, grid_u, grid_u.geoL, grid_v, grid_v.geoL,
        opC_pL, opC_uL, opC_vL,
        periodic_x, periodic_y, true)

    _, BuL, _, _, BvL, _, _, _ = set_crank_nicolson_block(neu, num,
        grid, opC_pL, Lp, bc_Lp, Lp_fs, bc_Lp_fs, BC_pL,
        grid_u, opC_uL, iRe.*Lu, iRe.*bc_Lu, iRe.*Lum1, iRe.*bc_Lum1, MuLm1, BC_uL,
        grid_v, opC_vL, iRe.*Lv, iRe.*bc_Lv, iRe.*Lvm1, iRe.*bc_Lvm1, MvLm1, BC_vL)

    derBu = copy(BuS)
    derBu.nzval .= 0.
    derBv = copy(BvS)
    derBv.nzval .= 0.

    # levelset advection derivatives
    RlsFS_ls.nzval .= 0.

    derLSA = copy(LSA)
    derLSA.nzval .= 0.
    derLSB = copy(LSB)
    derLSB.nzval .= 0.
    LSAj = copy(LSA)
    LSBj = copy(LSB)

    # graph coloring
    stencil = 7
    colors, graph, nloops, px0, pxf, py0, pyf = coloring(grid, stencil, periodic_x, periodic_y)
    colors_u, graph_u, _, px0_u, pxf_u, py0_u, pyf_u = coloring(grid_u, stencil, periodic_x, periodic_y)
    colors_v, graph_v, _, px0_v, pxf_v, py0_v, pyf_v = coloring(grid_v, stencil, periodic_x, periodic_y)
    @inbounds for color = 0:stencil^2-1
        @inbounds for p = 1:nloops
            indices = findall(colors[py0[p]:pyf[p], px0[p]:pxf[p]] .== color)
            relocate!(indices, px0[p], py0[p])
            assign_color!(grid, graph, indices, periodic_x, periodic_y)
            uj[indices] .+= ϵ_adj

            indices_u = findall(colors_u[py0_u[p]:pyf_u[p], px0_u[p]:pxf_u[p]] .== color)
            relocate!(indices_u, px0_u[p], py0_u[p])
            js_u = lexicographic.(indices_u, grid_u.ny)
            assign_color!(grid_u, graph_u, indices_u, periodic_x, periodic_y)
            veci(uD0_Sj,grid_u,2)[js_u] .+= ϵ_adj
            veci(uD0_Lj,grid_u,2)[js_u] .+= ϵ_adj

            indices_v = findall(colors_v[py0_v[p]:pyf_v[p], px0_v[p]:pxf_v[p]] .== color)
            relocate!(indices_v, px0_v[p], py0_v[p])
            js_v = lexicographic.(indices_v, grid_v.ny)
            assign_color!(grid_v, graph_v, indices_v, periodic_x, periodic_y)
            veci(vD0_Sj,grid_v,2)[js_v] .+= ϵ_adj
            veci(vD0_Lj,grid_v,2)[js_v] .+= ϵ_adj

            # SOLID PHASE
            # Compute capacities
            update_ls_data(num, grid, grid_u, grid_v, uj, κj, periodic_x, periodic_y)
            update_free_surface_velocity(num, grid_u, grid_v, uD0_S, vD0_S, periodic_x, periodic_y)

            Lp, bc_Lp, Lu, bc_Lu, Lv, bc_Lv, Lp_fs, bc_Lp_fs, Lu_fs, bc_Lu_fs, Lv_fs, bc_Lv_fs = set_laplacians!(grid, grid.geoS, grid_u, grid_u.geoS, grid_v, grid_v.geoS,
                opC_pj, opC_uj, opC_vj,
                periodic_x, periodic_y, true)
        
            _, Buj, _, _, Bvj, _, _, _ = set_crank_nicolson_block(neu, num,
                grid, opC_pj, Lp, bc_Lp, Lp_fs, bc_Lp_fs, BC_pS,
                grid_u, opC_uj, iRe.*Lu, iRe.*bc_Lu, iRe.*Lum1, iRe.*bc_Lum1, MuSm1, BC_uS,
                grid_v, opC_vj, iRe.*Lv, iRe.*bc_Lv, iRe.*Lvm1, iRe.*bc_Lvm1, MvSm1, BC_vS)

            # prediction step derivatives
            derBu .= (Buj .- BuS) ./ ϵ_adj
            derBv .= (Bvj .- BvS) ./ ϵ_adj

            Ruj = sparse(- derBu * uD0_S)
            Rvj = sparse(- derBv * vD0_S)

            rows = rowvals(Ruj)
            for i in nzrange(Ruj, 1)
                @inbounds row = rows[i]
                j = graph[row]
                @inbounds RucorrS_ls0[row,j] = Ruj[row]
            end
            rows = rowvals(Rvj)
            for i in nzrange(Rvj, 1)
                @inbounds row = rows[i]
                j = graph[row]
                @inbounds RvcorrS_ls0[row,j] = Rvj[row]
            end

            # levelset advection derivatives
            update_ls_data(num, grid, grid_u, grid_v, u0, grid.κ, periodic_x, periodic_y)
            update_free_surface_velocity(num, grid_u, grid_v, uD0_Sj, vD0_S, periodic_x, periodic_y)
            level_update_IIOE!(grid, grid_u, grid_v, LSAj, LSBj, θ_out, ind.MIXED, τ, periodic_x, periodic_y)
            utmp .= reshape(gmres(LSAj,(LSBj*vec(u0))), (ny,nx))
            S2IIOE!(grid, grid_u, grid_v, LSAj, LSBj, utmp, u0, θ_out, ind.MIXED, τ, periodic_x, periodic_y)

            derLSA .= (LSAj .- LSA) ./ ϵ_adj
            derLSB .= (LSBj .- LSB) ./ ϵ_adj
            Rj = sparse(derLSA * vec(u1) .- derLSB * vec(u0))

            rows = rowvals(Rj)
            for i in nzrange(Rj, 1)
                @inbounds row = rows[i]
                j = graph[row] + grid_u.ny*grid_u.nx
                @inbounds RlsFS_uS[row,j] = Rj[row]
            end

            update_ls_data(num, grid, grid_u, grid_v, u0, grid.κ, periodic_x, periodic_y)
            update_free_surface_velocity(num, grid_u, grid_v, uD0_S, vD0_Sj, periodic_x, periodic_y)
            level_update_IIOE!(grid, grid_u, grid_v, LSAj, LSBj, θ_out, ind.MIXED, τ, periodic_x, periodic_y)
            utmp .= reshape(gmres(LSAj,(LSBj*vec(u0))), (ny,nx))
            S2IIOE!(grid, grid_u, grid_v, LSAj, LSBj, utmp, u0, θ_out, ind.MIXED, τ, periodic_x, periodic_y)

            derLSA .= (LSAj .- LSA) ./ ϵ_adj
            derLSB .= (LSBj .- LSB) ./ ϵ_adj
            Rj = sparse(derLSA * vec(u1) .- derLSB * vec(u0))

            rows = rowvals(Rj)
            for i in nzrange(Rj, 1)
                @inbounds row = rows[i]
                j = graph[row] + grid_v.ny*grid_v.nx
                @inbounds RlsFS_vS[row,j] = Rj[row]
            end

            # LIQUID PHASE
            # Compute capacities
            update_ls_data(num, grid, grid_u, grid_v, uj, κj, periodic_x, periodic_y)
            update_free_surface_velocity(num, grid_u, grid_v, uD0_L, vD0_L, periodic_x, periodic_y)

            Lp, bc_Lp, Lu, bc_Lu, Lv, bc_Lv, Lp_fs, bc_Lp_fs, Lu_fs, bc_Lu_fs, Lv_fs, bc_Lv_fs = set_laplacians!(grid, grid.geoL, grid_u, grid_u.geoL, grid_v, grid_v.geoL,
                opC_pj, opC_uj, opC_vj,
                periodic_x, periodic_y, true)
        
            _, Buj, _, _, Bvj, _, _, _ = set_crank_nicolson_block(neu, num,
                grid, opC_pj, Lp, bc_Lp, Lp_fs, bc_Lp_fs, BC_pL,
                grid_u, opC_uj, iRe.*Lu, iRe.*bc_Lu, iRe.*Lum1, iRe.*bc_Lum1, MuLm1, BC_uL,
                grid_v, opC_vj, iRe.*Lv, iRe.*bc_Lv, iRe.*Lvm1, iRe.*bc_Lvm1, MvLm1, BC_vL)

            # prediction step derivatives
            derBu .= (Buj .- BuL) ./ ϵ_adj
            derBv .= (Bvj .- BvL) ./ ϵ_adj

            Ruj = sparse(- derBu * uD0_L)
            Rvj = sparse(- derBv * vD0_L)

            rows = rowvals(Ruj)
            for i in nzrange(Ruj, 1)
                @inbounds row = rows[i]
                j = graph[row]
                @inbounds RucorrL_ls0[row,j] = Ruj[row]
            end
            rows = rowvals(Rvj)
            for i in nzrange(Rvj, 1)
                @inbounds row = rows[i]
                j = graph[row]
                @inbounds RvcorrL_ls0[row,j] = Rvj[row]
            end

            # levelset advection derivatives
            level_update_IIOE!(grid, grid_u, grid_v, LSAj, LSBj, θ_out, ind.MIXED, τ, periodic_x, periodic_y)
            utmp .= reshape(gmres(LSAj,(LSBj*vec(uj))), (ny,nx))
            S2IIOE!(grid, grid_u, grid_v, LSAj, LSBj, utmp, uj, θ_out, ind.MIXED, τ, periodic_x, periodic_y)

            derLSA .= (LSAj .- LSA) ./ ϵ_adj
            derLSB .= (LSBj .- LSB) ./ ϵ_adj

            utmp = zeros(ny, nx)
            utmp[indices] .= 1.0
            Rj = sparse(derLSA * vec(u1) .- derLSB * vec(u0) .- LSB * vec(utmp))

            # We should get exactly the same with both solid and liquid phases,
            # since the jump of velocity fields at the interface should be 0
            rows = rowvals(Rj)
            for i in nzrange(Rj, 1)
                @inbounds row = rows[i]
                j = graph[row]
                @inbounds RlsFS_ls[row,j] = Rj[row]
            end

            update_ls_data(num, grid, grid_u, grid_v, u0, grid.κ, periodic_x, periodic_y)
            update_free_surface_velocity(num, grid_u, grid_v, uD0_Lj, vD0_L, periodic_x, periodic_y)
            level_update_IIOE!(grid, grid_u, grid_v, LSAj, LSBj, θ_out, ind.MIXED, τ, periodic_x, periodic_y)
            utmp .= reshape(gmres(LSAj,(LSBj*vec(u0))), (ny,nx))
            S2IIOE!(grid, grid_u, grid_v, LSAj, LSBj, utmp, u0, θ_out, ind.MIXED, τ, periodic_x, periodic_y)

            derLSA .= (LSAj .- LSA) ./ ϵ_adj
            derLSB .= (LSBj .- LSB) ./ ϵ_adj
            Rj = sparse(derLSA * vec(u1) .- derLSB * vec(u0))

            rows = rowvals(Rj)
            for i in nzrange(Rj, 1)
                @inbounds row = rows[i]
                j = graph[row] + grid_u.ny*grid_u.nx
                @inbounds RlsFS_uL[row,j] = Rj[row]
            end

            update_ls_data(num, grid, grid_u, grid_v, u0, grid.κ, periodic_x, periodic_y)
            update_free_surface_velocity(num, grid_u, grid_v, uD0_L, vD0_Lj, periodic_x, periodic_y)
            level_update_IIOE!(grid, grid_u, grid_v, LSAj, LSBj, θ_out, ind.MIXED, τ, periodic_x, periodic_y)
            utmp .= reshape(gmres(LSAj,(LSBj*vec(u0))), (ny,nx))
            S2IIOE!(grid, grid_u, grid_v, LSAj, LSBj, utmp, u0, θ_out, ind.MIXED, τ, periodic_x, periodic_y)

            derLSA .= (LSAj .- LSA) ./ ϵ_adj
            derLSB .= (LSBj .- LSB) ./ ϵ_adj
            Rj = sparse(derLSA * vec(u1) .- derLSB * vec(u0))

            rows = rowvals(Rj)
            for i in nzrange(Rj, 1)
                @inbounds row = rows[i]
                j = graph[row] + grid_v.ny*grid_v.nx
                @inbounds RlsFS_vL[row,j] = Rj[row]
            end

            uj .= u0
            uD0_Sj .= uD0_S
            vD0_Sj .= vD0_S
            uD0_Lj .= uD0_L
            vD0_Lj .= vD0_L
        end
    end

    return nothing
end