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

function Rheat_u(num, grid, grid_u, grid_v, adj_der, um1,
    TD0_S, TD1_S, A_S, B_S, opC_TS, BC_TS,
    TD0_L, TD1_L, A_L, B_L, opC_TL, BC_TL,
    u0, u1, LSA, LSB, tmpχ_S, tmpχ_L,
    CFL_sc, periodic_x, periodic_y, ϵ_adj, λ, Vmean)

    @unpack ϵ, NB = num
    @unpack nx, ny, ind, V, iso, faces, geoS, geoL = grid
    @unpack RheatS_u, RheatL_u, RlsS_u = adj_der
    
    RheatS_u.nzval .= 0.
    RheatL_u.nzval .= 0.
    RlsS_u.nzval .= 0.
    uj = copy(um1)

    TS = reshape(veci(TD1_S, grid, 1), (ny, nx))
    TL = reshape(veci(TD1_L, grid, 1), (ny, nx))

    derA_S = copy(A_S)
    derA_S.nzval .= 0.
    derB_S = copy(B_S)
    derB_S.nzval .= 0.

    derA_L = copy(A_L)
    derA_L.nzval .= 0.
    derB_L = copy(B_L)
    derB_L.nzval .= 0.

    derLSA = copy(LSA)
    derLSA.nzval .= 0.
    derLSB = copy(LSB)
    derLSB.nzval .= 0.
    LSAj = copy(LSA)
    LSBj = copy(LSB)

    derχ_S = copy(opC_TS.χ)
    derχ_L = copy(opC_TL.χ)
    bc_S = zeros(2*ny*nx)
    bc_L = zeros(2*ny*nx)
    a0 = num.θd .* ones(nx*ny)

    # graph coloring
    stencil = 7
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
    @inbounds for color = 0:stencil^2-1
        @inbounds for p = 1:nloops
            indices = findall(colors[py0[p]:pyf[p], px0[p]:pxf[p]] .== color)
            relocate!(indices, px0[p], py0[p])
            assign_color!(grid, graph, indices, periodic_x, periodic_y)

            uj[indices] .+= ϵ_adj

            # Compute capacities
            update_ls_data(num, grid, grid_u, grid_v, uj, periodic_x, periodic_y)
            update_stefan_velocity(num, grid, uj, TS, TL, periodic_x, periodic_y, λ, Vmean)

            # get perturbed matrices
            Aj_S, Bj_S, _ = set_heat!(dir, num, grid, opC_TS, geoS, BC_TS, grid.ind.MIXED, geoS.projection,
                                    periodic_x, periodic_y)
            derA_S .= (Aj_S .- A_S) ./ ϵ_adj
            derB_S .= (Bj_S .- B_S) ./ ϵ_adj
            derχ_S .= (opc_TS.χ .-  tmpχ_S) ./ ϵ_adj
            bc_S = derχ_S * a0
            Rj_S = sparse(derA_S * TD1_S .- derB_S * TD0_S .- bc_S)

            Aj_L, Bj_L, _ = set_heat!(dir, num, grid, opC_TL, geoL, BC_TL, grid.ind.MIXED, geoL.projection,
                                    periodic_x, periodic_y)
            derA_L .= (Aj_L .- A_L) ./ ϵ_adj
            derB_L .= (Bj_L .- B_L) ./ ϵ_adj
            derχ_L .= (opc_TL.χ .-  tmpχ_L) ./ ϵ_adj
            bc_L = derχ_L * a0
            Rj_L = sparse(derA_L * TD1_L .- derB_L * TD0_L .- bc_L)

            IIOE(grid, LSAj, LSBj, uj, V, CFL_sc, periodic_x, periodic_y)
            derLSA .= (LSAj .- LSA) ./ ϵ_adj
            derLSB .= (LSBj .- LSB) ./ ϵ_adj
            utmp = zeros(ny, nx)
            utmp[indices] .= 1.0
            Rj_u = sparse(derLSA * vec(u1) .- derLSB * vec(u0) .- LSB * vec(utmp))

            # Do NOT parallelize. Sparsity pattern changes due to periodic BCs
            # not being preallocated
            rows = rowvals(Rj_S)
            for i in nzrange(Rj_S, 1)
                @inbounds row = rows[i]
                j = graph[row]
                @inbounds RheatS_u[row,j] = Rj_S[row]
            end
            rows = rowvals(Rj_L)
            for i in nzrange(Rj_L, 1)
                @inbounds row = rows[i]
                j = graph[row]
                @inbounds RheatL_u[row,j] = Rj_L[row]
            end
            rows = rowvals(Rj_u)
            for i in nzrange(Rj_u, 1)
                @inbounds row = rows[i]
                j = graph[row]
                @inbounds RlsS_u[row,j] = Rj_u[row]
            end

            uj .= um1
        end
    end

    return nothing
end

function Rheat_T(num, grid, grid_u, grid_v, adj_der, 
    TD_S, TD_L,
    u0, u1, LSA, LSB,
    CFL_sc, periodic_x, periodic_y, ϵ_adj, λ, Vmean)

    @unpack ϵ, NB = num
    @unpack nx, ny, ind, u, V, faces, iso, geoS, geoL = grid
    @unpack RlsS_TS, RlsS_TL = adj_der
    
    RlsS_TS.nzval .= 0.
    RlsS_TL.nzval .= 0.

    TS = reshape(veci(TD_S, grid, 1), (ny, nx))
    TL = reshape(veci(TD_L, grid, 1), (ny, nx))
    TSj = copy(TS)
    TLj = copy(TL)
    
    derLSA_S = copy(LSA)
    derLSA_S.nzval .= 0.
    derLSB_S = copy(LSB)
    derLSB_S.nzval .= 0.

    derLSA_L = copy(LSA)
    derLSA_L.nzval .= 0.
    derLSB_L = copy(LSB)
    derLSB_L.nzval .= 0.

    LSAj = copy(LSA)
    LSBj = copy(LSB)

    # graph coloring
    stencil = 7
    ix = last.(getproperty.(ind.all_indices, :I))
    iy = first.(getproperty.(ind.all_indices, :I))
    colors = (ix .+ 2stencil*iy) .% (stencil^2)
    graph = zeros(Int64,2nx*ny)
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
    @inbounds for color = 0:stencil^2-1
        @inbounds for p = 1:nloops
            indices = findall(colors[py0[p]:pyf[p], px0[p]:pxf[p]] .== color)
            relocate!(indices, px0[p], py0[p])
            js = lexicographic.(indices, ny)
            assign_color!(grid, graph, indices, periodic_x, periodic_y)

            TSj[js] .+= ϵ_adj
            TLj[js] .+= ϵ_adj

            # Compute capacities
            update_ls_data(num, grid, grid_u, grid_v, u, periodic_x, periodic_y)

            # get perturbed matrices
            update_stefan_velocity(num, grid, u, TSj, TL, periodic_x, periodic_y, λ, Vmean)
            IIOE(grid, LSAj, LSBj, u, V, CFL_sc, periodic_x, periodic_y)
            derLSA_S .= (LSAj .- LSA) ./ ϵ_adj
            derLSB_S .= (LSBj .- LSB) ./ ϵ_adj
            Rj_S = sparse(derLSA_S * vec(u1) .- derLSB_S * vec(u0))

            update_stefan_velocity(num, grid, u, TS, TLj, periodic_x, periodic_y, λ, Vmean)
            IIOE(grid, LSAj, LSBj, u, V, CFL_sc, periodic_x, periodic_y)
            derLSA_L .= (LSAj .- LSA) ./ ϵ_adj
            derLSB_L .= (LSBj .- LSB) ./ ϵ_adj
            Rj_L = sparse(derLSA_L * vec(u1) .- derLSB_L * vec(u0))

            # Do NOT parallelize. Sparsity pattern changes due to periodic BCs
            # not being preallocated
            rows = rowvals(Rj_S)
            for i in nzrange(Rj_S, 1)
                @inbounds row = rows[i]
                j = graph[row]
                @inbounds RlsS_TS[row,j] = Rj_S[row]
            end
            rows = rowvals(Rj_L)
            for i in nzrange(Rj_L, 1)
                @inbounds row = rows[i]
                j = graph[row]
                @inbounds RlsS_TL[row,j] = Rj_L[row]
            end

            TSj .= TS
            TLj .= TL
        end
    end

    return nothing
end