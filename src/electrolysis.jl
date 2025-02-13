#From set_heat in heat_coupled.jl, poisson.jl

function set_convection_2!(
    num, grid, geo, grid_u, LS_u, grid_v, LS_v,
    u, v, op, ph, BC_u, BC_v
    )
    @unpack Cu, CUTCu, Cv, CUTCv = op
    @unpack uD, vD = ph

    Du_x = zeros(grid_u)
    Du_y = zeros(grid_u)
    for iLS in 1:num.nLS
        Du_x[LS_u[iLS].MIXED] .= reshape(veci(uD,grid_u,iLS+1), grid_u)[LS_u[iLS].MIXED]
        Du_y[LS_u[iLS].MIXED] .= reshape(veci(uD,grid_u,iLS+1), grid_u)[LS_u[iLS].MIXED]
    end
    Du_x[:,1] .= vecb_L(uD,grid_u)
    Du_y[1,:] .= vecb_B(uD,grid_u)
    Du_x[:,end] .= vecb_R(uD,grid_u)
    Du_y[end,:] .= vecb_T(uD,grid_u)

    # bnds_u = [grid_u.ind.b_left[1], grid_u.ind.b_bottom[1], grid_u.ind.b_right[1], grid_u.ind.b_top[1]]
    # bnds_v = [grid_v.ind.b_left[1], grid_v.ind.b_bottom[1], grid_v.ind.b_right[1], grid_v.ind.b_top[1]]
    # Δu = [grid_u.dx[1,1] * 0.25, grid_u.dy[1,1] * 0.5, grid_u.dx[end,end] * 0.25, grid_u.dy[end,end] * 0.5]
    # Δv = [grid_v.dx[1,1] * 0.5, grid_v.dy[1,1] * 0.25, grid_v.dx[end,end] * 0.5, grid_v.dy[end,end] * 0.25]

    # Hu = zeros(grid_u)
    # for i in eachindex(bnds_u)
    #     for II in bnds_u[i]
    #         Hu[II] = Δu[i]
    #     end
    # end

    # Hv = zeros(grid_v)
    # for i in eachindex(bnds_v)
    #     for II in bnds_v[i]
    #         Hv[II] = Δv[i]
    #     end
    # end


    vector_convection!(dir, GridFCx, Cu, CUTCu, u, v, Du_x, Du_y, Dv_x, Dv_y,
            geo.dcap, grid.nx, grid.ny, BC_u, grid_u.ind.inside,
            grid_u.ind.b_left[1], grid_u.ind.b_bottom[1], grid_u.ind.b_right[1], grid_u.ind.b_top[1])
    
    return nothing
end


    function set_bc_bnds_total(::Dirichlet, D, H, BC)
        Dx = copy(D)
        Dy = copy(D)

        if is_neumann(BC.left)
            @inbounds Dx[:,1] .= H[:,1] .* BC.left.val
        elseif is_dirichlet(BC.left)
            @inbounds Dx[:,1] .= BC.left.val
        end
        @inbounds Dy[:,1] .= Dx[:,1]

        if is_neumann(BC.bottom)
            @inbounds Dy[1,:] .= H[1,:] .* BC.bottom.val 
        elseif is_dirichlet(BC.bottom)
            @inbounds Dy[1,:] .= BC.bottom.val
        end
        @inbounds Dx[1,:] .= Dy[1,:]

        if is_neumann(BC.right)
            @inbounds Dx[:,end] .= H[:,end] .* BC.right.val 
        elseif is_dirichlet(BC.right)
            @inbounds Dx[:,end] .= BC.right.val
        end
        @inbounds Dy[:,end] .= Dx[:,end]

        if is_neumann(BC.top)
            @inbounds Dy[end,:] .= H[end,:] .* BC.top.val 
        elseif is_dirichlet(BC.top)
            @inbounds Dy[end,:] .= BC.top.val
        end
        @inbounds Dx[end,:] .= Dy[end,:]


        return Dx, Dy
    end


    """
    set_scalar_boundaries

    # The idea there is to have at the corners the contributions from both borders. 
    # bcTx is only used for derivatives in x and bcTy for derivatives in y, so it doesn't matter 
    # that bcTx is not filled at the top and bottom and the same for bcTy
    # Otherwise we would have to select one of the contributions

    #TODO WARNING: only first interfacial value is filled
    """
    function set_scalar_boundaries!(
        num::Numerical{Float64, Int64},
        grid::Mesh{Flower.GridCC, Float64, Int64},bc::BoundariesInt,scalD,Dx,Dy)
  
        HT = zeros(grid)

        # Dx = zeros(grid)
        # Dy = zeros(grid)
        # Dx .= 0.0
        # Dy .= 0.0

        iLS = 1 #TODO different for multiple LS?
        # @inbounds @threads for II in vcat(b_left[1], b_bottom[1], b_right[1], b_top[1])
        @inbounds for II in vcat(grid.ind.b_left[1], grid.ind.b_bottom[1], grid.ind.b_right[1], grid.ind.b_top[1])
            HT[II] = distance(grid.LS[iLS].mid_point[II], grid.LS[iLS].geoL.centroid[II], grid.dx[II], grid.dy[II])
        end  

        if num.convection_mode == 0

            #Interface boundary condition
            if is_dirichlet(bc.int)
                __a0 = bc.int.val        
            elseif is_neumann(bc.int)
                __a0 = bc.int.val
            elseif is_robin(bc.int)
                __a0 = bc_type.val  
            elseif is_stefan(bc.int)
                __a0 = num.concentration0[iscal]       
            elseif is_wall(bc.int)
                __a0 = bc.int.val         
            else
                __a0 = bc.int.val
            end

            # Flags with BCs
            Dx .= __a0
            Dy .= __a0

            # Dx, Dy = set_bc_bnds(dir, a0, HT, bc)

            if is_neumann(bc.left)
                @inbounds Dx[:,1] .= HT[:,1] .* bc.left.val
            elseif is_dirichlet(bc.left)
                @inbounds Dx[:,1] .= bc.left.val
            end
            if is_neumann(bc.bottom)
                @inbounds Dy[1,:] .= HT[1,:] .* bc.bottom.val 
            elseif is_dirichlet(bc.bottom)
                @inbounds Dy[1,:] .= bc.bottom.val
            end
            if is_neumann(bc.right)
                @inbounds Dx[:,end] .= HT[:,end] .* bc.right.val 
            elseif is_dirichlet(bc.right)
                @inbounds Dx[:,end] .= bc.right.val
            end
            if is_neumann(bc.top)
                @inbounds Dy[end,:] .= HT[end,:] .* bc.top.val 
            elseif is_dirichlet(bc.top)
                @inbounds Dy[end,:] .= bc.top.val
            end

        elseif num.convection_mode == 1

          

            # for iLS in 1:num.nLS
                
            #     χx = (grid.LS[iLS].geoL.dcap[:,:,3] .- grid.LS[iLS].geoL.dcap[:,:,1]) .^ 2
            #     χy = (grid.LS[iLS].geoL.dcap[:,:,4] .- grid.LS[iLS].geoL.dcap[:,:,2]) .^ 2
            #     op.χ[iLS].diag .= sqrt.(vec(χx .+ χy))
            #     a0 .+= __a0[iscal] * op.χ[iLS] 
            # end

        
            for iLS in 1:num.nLS
                Dx[grid.LS[iLS].MIXED] .= reshape(veci(scalD,grid,iLS+1), grid)[grid.LS[iLS].MIXED]
                Dy[grid.LS[iLS].MIXED] .= reshape(veci(scalD,grid,iLS+1), grid)[grid.LS[iLS].MIXED]
            end
        
            #either 
            # Dx[:,1] .= vecb_L(uD,grid_u) 
            # Dy[1,:] .= vecb_B(uD,grid_u)
            # Dx[:,end] .= vecb_R(uD,grid_u)
            # Dy[end,:] .= vecb_T(uD,grid_u)

            #or
            if is_neumann(bc.left)
                @inbounds Dx[:,1] .= HT[:,1] .* bc.left.val
            elseif is_dirichlet(bc.left)
                @inbounds Dx[:,1] .= bc.left.val
            end
            if is_neumann(bc.bottom)
                @inbounds Dy[1,:] .= HT[1,:] .* bc.bottom.val 
            elseif is_dirichlet(bc.bottom)
                @inbounds Dy[1,:] .= bc.bottom.val
            end
            if is_neumann(bc.right)
                @inbounds Dx[:,end] .= HT[:,end] .* bc.right.val 
            elseif is_dirichlet(bc.right)
                @inbounds Dx[:,end] .= bc.right.val
            end
            if is_neumann(bc.top)
                @inbounds Dy[end,:] .= HT[end,:] .* bc.top.val 
            elseif is_dirichlet(bc.top)
                @inbounds Dy[end,:] .= bc.top.val
            end

        end # if num.convection_mode
    end #function 


"""    
scalar transport (convection and diffusion)

interface term in rhs comes from op.opC_TL

# Arguments
- num, grid, grid_u, grid_v
- BC: BC for interface, 
- a0, 
- opC, 
- A, 
- BC: BC for wall
- ls_advection

Example
```math
\\frac{\\partial c_{H_2}}{\\partial t} +  \\nabla \\cdot \\left(c_{H_2}\\mathbf{u}\\right)=D_{H_2} \\nabla^2 c_{H_2}
```

Divergence
```julia
    BxT * iMx * Hx[iLS] .+ ByT * iMy * Hy[iLS]
```
""" 
function scalar_transport!(num::Numerical{Float64, Int64},
    grid::Mesh{Flower.GridCC, Float64, Int64}, 
    grid_u::Mesh{Flower.GridFCx, Float64, Int64},
    grid_v::Mesh{Flower.GridFCy, Float64, Int64},
    op::Operators{Float64, Int64},
    op_conv::OperatorsConvection{Float64, Int64},
    ph::Phase{Float64}, 
    bc::Vector{BoundariesInt},
    BC_int::Vector{<:BoundaryCondition}, 
    A::SparseMatrixCSC{Float64, Int64},
    B::SparseArrays.SparseMatrixCSC{Float64, Int64},
    all_CUTCT::Array{Float64, 2},
    rhs::Array{Float64, 1},
    a0::Array{Float64, 2},
    tmp_vec_u::Array{Float64, 2},
    tmp_vec_v::Array{Float64, 2},
    periodic_x::Bool, 
    periodic_y::Bool, 
    convection::Bool, 
    ls_advection::Bool)

    @unpack nx, ny, dx, dy, ind, LS  = grid
    @unpack all_indices, inside, b_left, b_bottom, b_right, b_top = ind
    @unpack Bx, By, BxT, ByT, Hx, Hy, HxT, HyT, M, iMx, iMy, χ = op
    @unpack u, v, uD, vD = ph

    # printstyled(color=:red, @sprintf "\n levelset: start scalar_transport!\n")
    # print("\n nb_transported_scalars ",num.nb_transported_scalars)
    # println("\n grid.LS[1].geoL.dcap[1,1,:]",grid.LS[1].geoL.dcap[1,1,:])

    if convection

    else
        printstyled(color=:red, @sprintf "\n no convection\n")

    end

    ni = nx * ny #Size in the bulk 
    nb = 2 * nx + 2 * ny #for cells in the borders of the domain

    # Reset arrays to zero
    #TODO use fill or nzval ? check if correct
    A .= 0.0
    B .= 0.0
    # A.nzval .= 0.0
    # B.nzval .= 0.0

    all_CUTCT .=0.0
    a0 .= 0.0

    # rhs_test = fnzeros(grid, num)

    #tmp vec scalar
    tmp_vec_p_2 = zeros(grid)
    tmp_vec_p_3 = zeros(grid)

    
    # Operators
    if convection
        
        if num.convection_mode == 0

            #old
            tmp_vec_u = zeros(grid_u)
            tmp_vec_u .= reshape(vec2(uD,grid_u), grid_u)
            tmp_vec_u[1,:] .= vecb_B(uD, grid_u)
            tmp_vec_u[end,:] .= vecb_T(uD, grid_u)
            tmp_vec_u[:,1] .= vecb_L(uD, grid_u)
            tmp_vec_u[:,end] .= vecb_R(uD, grid_u)
            
            tmp_vec_v = zeros(grid_v)
            tmp_vec_v .= reshape(vec2(vD,grid_v), grid_v)
            tmp_vec_v[:,1] .= vecb_L(vD, grid_v)
            tmp_vec_v[:,end] .= vecb_R(vD, grid_v)
            tmp_vec_v[1,:] .= vecb_B(vD, grid_v)
            tmp_vec_v[end,:] .= vecb_T(vD, grid_v)

        elseif num.convection_mode ==1

            tmp_vec_u .= 0.0 #reset to zero
            for iLS in 1:num.nLS
                tmp_vec_u[grid_u.LS[iLS].MIXED] .= reshape(veci(uD,grid_u,iLS+1), grid_u)[grid_u.LS[iLS].MIXED]
            end
            tmp_vec_u[:,1] .= vecb_L(uD,grid_u) 
            tmp_vec_u[1,:] .= vecb_B(uD,grid_u)
            tmp_vec_u[:,end] .= vecb_R(uD,grid_u)
            tmp_vec_u[end,:] .= vecb_T(uD,grid_u)
        
            tmp_vec_v .= 0.0
            for iLS in 1:num.nLS
                tmp_vec_v[grid_v.LS[iLS].MIXED] .= reshape(veci(vD,grid_v,iLS+1), grid_v)[grid_v.LS[iLS].MIXED]
            end
            tmp_vec_v[:,1] .= vecb_L(vD,grid_v)
            tmp_vec_v[1,:] .= vecb_B(vD,grid_v)
            tmp_vec_v[:,end] .= vecb_R(vD,grid_v)
            tmp_vec_v[end,:] .= vecb_T(vD,grid_v)

            # In NS
            # Du_x = zeros(grid_u)
            # Du_y = zeros(grid_u)
            # for iLS in 1:num.nLS
            #     Du_x[LS_u[iLS].MIXED] .= reshape(veci(uD,grid_u,iLS+1), grid_u)[LS_u[iLS].MIXED]
            #     Du_y[LS_u[iLS].MIXED] .= reshape(veci(uD,grid_u,iLS+1), grid_u)[LS_u[iLS].MIXED]
            # end
            # Du_x[:,1] .= vecb_L(uD,grid_u) 
            # Du_y[1,:] .= vecb_B(uD,grid_u)
            # Du_x[:,end] .= vecb_R(uD,grid_u)
            # Du_y[end,:] .= vecb_T(uD,grid_u)
        
            # Dv_x = zeros(grid_v)
            # Dv_y = zeros(grid_v)
            # for iLS in 1:num.nLS
            #     Dv_x[LS_v[iLS].MIXED] .= reshape(veci(vD,grid_v,iLS+1), grid_v)[LS_v[iLS].MIXED]
            #     Dv_y[LS_v[iLS].MIXED] .= reshape(veci(vD,grid_v,iLS+1), grid_v)[LS_v[iLS].MIXED]
            # end
            # Dv_x[:,1] .= vecb_L(vD,grid_v)
            # Dv_y[1,:] .= vecb_B(vD,grid_v)
            # Dv_x[:,end] .= vecb_R(vD,grid_v)
            # Dv_y[end,:] .= vecb_T(vD,grid_v)
        end #num.convection_mode

        
        #TODO HT multiple LS
        # HT = zeros(grid)
        # # @inbounds @threads for II in vcat(b_left[1], b_bottom[1], b_right[1], b_top[1])
        # @inbounds for II in vcat(b_left[1], b_bottom[1], b_right[1], b_top[1])
        #     HT[II] = distance(grid.LS[1].mid_point[II], grid.LS[1].geoL.centroid[II], dx[II], dy[II])
        # end  

        iscal = 1

        # Dx = zeros(grid)
        # Dy = zeros(grid)
        #values for convection
        set_scalar_boundaries!(num,grid, bc[iscal], ph.trans_scalD[:,iscal], tmp_vec_p_2, tmp_vec_p_3)

        # @views necessary
        @views scalar_convection!(dir, op_conv.CT, all_CUTCT[:,iscal], u, v, tmp_vec_p_2, tmp_vec_p_3,
        tmp_vec_u, tmp_vec_v, grid.LS[end].geoL.dcap, ny, 
        bc[iscal], inside, b_left[1], b_bottom[1], b_right[1], b_top[1]
        )
        # @views scalar_convection!(dir, op_conv.CT, all_CUTCT[:,iscal], u, v, bcTx, bcTy, bcU, bcV, grid.LS[1].geoL.dcap, ny, 
        # bc[iscal], inside, b_left[1], b_bottom[1], b_right[1], b_top[1]
        # )

        if num.nb_transported_scalars>1
            for iscal=2:num.nb_transported_scalars

                #reset vectors
                tmp_vec_p_2 .= 0.0
                tmp_vec_p_3 .= 0.0
                #values for convection
                set_scalar_boundaries!(num,grid, bc[iscal], ph.trans_scalD[:,iscal], tmp_vec_p_2, tmp_vec_p_3)

                # @views necessary
                @views scalar_convection_CUTCT!(dir, all_CUTCT[:,iscal], u, v, tmp_vec_p_2, tmp_vec_p_3, 
                tmp_vec_u, tmp_vec_v, grid.LS[end].geoL.dcap, ny, 
                bc[iscal], inside, b_left[1], b_bottom[1], b_right[1], b_top[1]
                )
                # @views scalar_convection_CUTCT!(dir, all_CUTCT[:,iscal], u, v, bcTx, bcTy, bcU, bcV, grid.LS[1].geoL.dcap, ny, 
                # bc[iscal], inside, b_left[1], b_bottom[1], b_right[1], b_top[1]
                # )
            end
        end       
    end #convection

    # printstyled(color=:red, @sprintf "\n levelset: end scalar_convection\n")
    # println(grid.LS[1].geoL.dcap[1,1,:])


    if num.nLS ==1
        #error iLS 1 
        if ls_advection
            update_all_ls_data(num, grid, grid_u, grid_v, BC_int, periodic_x, periodic_y, false)

            # Mass matrices
            M.diag .= vec(grid.LS[1].geoL.dcap[:,:,5])
            Mx = zeros(ny,nx+1)
            for II in ind.all_indices
                Mx[II] = grid.LS[1].geoL.dcap[II,8]
            end
            for II in ind.b_right[1]
                Mx[δx⁺(II)] = grid.LS[1].geoL.dcap[II,10]
            end
            My = zeros(ny+1,nx)
            for II in ind.all_indices
                My[II] = grid.LS[1].geoL.dcap[II,9]
            end
            for II in ind.b_top[1]
                My[δy⁺(II)] = grid.LS[1].geoL.dcap[II,11]
            end
        
            # iMx.diag .= 1. ./ (vec(Mx) .+ eps(0.01))
            # iMy.diag .= 1. ./ (vec(My) .+ eps(0.01))
            # iMx = Diagonal(inv_weight_eps2.(num.epsilon_mode,num.epsilon_vol,vec(Mx)))
            # iMy = Diagonal(inv_weight_eps2.(num.epsilon_mode,num.epsilon_vol,vec(My)))
        
            iMx.diag .= inv_weight_eps2.(num.epsilon_mode,num.epsilon_vol,vec(Mx))
            iMy.diag .= inv_weight_eps2.(num.epsilon_mode,num.epsilon_vol,vec(My))


            # Discrete gradient and divergence operators
            divergence_B!(BxT, ByT, grid.LS[1].geoL.dcap, ny, ind.all_indices)
            mat_assign!(Bx, sparse(-BxT'))
            mat_assign!(By, sparse(-ByT'))

            # Matrices for interior BCs
            for iLS in 1:num.nLS
                bc_matrix!(grid, Hx[iLS], Hy[iLS], grid.LS[1].geoL.dcap, grid.LS[1].geoL.dcap, ny, ind.all_indices)

                mat_assign_T!(HxT[iLS], sparse(Hx[iLS]'))
                mat_assign_T!(HyT[iLS], sparse(Hy[iLS]'))

                periodic_bcs!(grid, Bx, By, Hx[iLS], Hy[iLS], periodic_x, periodic_y)

                χx = (grid.LS[1].geoL.dcap[:,:,3] .- grid.LS[1].geoL.dcap[:,:,1]) .^ 2
                χy = (grid.LS[1].geoL.dcap[:,:,4] .- grid.LS[1].geoL.dcap[:,:,2]) .^ 2
                op.χ[iLS].diag .= sqrt.(vec(χx .+ χy))
            end
            mat_assign!(BxT, sparse(-Bx'))
            mat_assign!(ByT, sparse(-By'))

            # Matrices for borders BCs
            set_boundary_indicator!(grid, grid.LS[1].geoL, grid.LS[1].geoL, op)
            mass_matrix_borders!(num,ind, op.iMx_b, op.iMy_b, op.iMx_bd, op.iMy_bd, grid.LS[1].geoL.dcap, ny)
            bc_matrix_borders!(grid, ind, op.Hx_b, op.Hy_b, grid.LS[1].geoL.dcap)
            mat_assign_T!(op.HxT_b, sparse(op.Hx_b'))
            mat_assign_T!(op.HyT_b, sparse(op.Hy_b'))
            periodic_bcs_borders!(grid, op.Hx_b, op.Hy_b, periodic_x, periodic_y)
        end

    else
        # #TODO better alloc
        # geo = [grid.LS[iLS].geoL for iLS in 1:num._nLS]
        # geo_u = [grid_u.LS[iLS].geoL for iLS in 1:num._nLS]
        # geo_v = [grid_v.LS[iLS].geoL for iLS in 1:num._nLS]

        # laps = set_matrices!(
        #     num, grid, geo, grid_u, geo_u, grid_v, geo_v,
        #     opC_p, opC_u, opC_v,
        #     periodic_x, periodic_y)

         #error iLS 1 
         if ls_advection
            update_all_ls_data(num, grid, grid_u, grid_v, BC_int, periodic_x, periodic_y, false)

            # Mass matrices
            M.diag .= vec(grid.LS[end].geoL.dcap[:,:,5])
            Mx = zeros(ny,nx+1)
            for II in ind.all_indices
                Mx[II] = grid.LS[end].geoL.dcap[II,8]
            end
            for II in ind.b_right[1]
                Mx[δx⁺(II)] = grid.LS[end].geoL.dcap[II,10]
            end
            My = zeros(ny+1,nx)
            for II in ind.all_indices
                My[II] = grid.LS[end].geoL.dcap[II,9]
            end
            for II in ind.b_top[1]
                My[δy⁺(II)] = grid.LS[end].geoL.dcap[II,11]
            end
        
            # iMx.diag .= 1. ./ (vec(Mx) .+ eps(0.01))
            # iMy.diag .= 1. ./ (vec(My) .+ eps(0.01))
            # iMx = Diagonal(inv_weight_eps2.(num.epsilon_mode,num.epsilon_vol,vec(Mx)))
            # iMy = Diagonal(inv_weight_eps2.(num.epsilon_mode,num.epsilon_vol,vec(My)))
        
            iMx.diag .= inv_weight_eps2.(num.epsilon_mode,num.epsilon_vol,vec(Mx))
            iMy.diag .= inv_weight_eps2.(num.epsilon_mode,num.epsilon_vol,vec(My))


            # Discrete gradient and divergence operators
            divergence_B!(BxT, ByT, grid.LS[end].geoL.dcap, ny, ind.all_indices)
            mat_assign!(Bx, sparse(-BxT'))
            mat_assign!(By, sparse(-ByT'))

            # Matrices for interior BCs
            for iLS in 1:num.nLS
                bc_matrix!(grid, Hx[iLS], Hy[iLS], grid.LS[iLS].geoL.dcap, grid.LS[iLS].geoL.dcap, ny, ind.all_indices)

                mat_assign_T!(HxT[iLS], sparse(Hx[iLS]'))
                mat_assign_T!(HyT[iLS], sparse(Hy[iLS]'))

                periodic_bcs!(grid, Bx, By, Hx[iLS], Hy[iLS], periodic_x, periodic_y)

                χx = (grid.LS[iLS].geoL.dcap[:,:,3] .- grid.LS[iLS].geoL.dcap[:,:,1]) .^ 2
                χy = (grid.LS[iLS].geoL.dcap[:,:,4] .- grid.LS[iLS].geoL.dcap[:,:,2]) .^ 2
                op.χ[iLS].diag .= sqrt.(vec(χx .+ χy))

                #test χ
                # veci(rhs_test,grid,iLS+1) .+= op.χ[iLS] * vec(ones(grid))
            end
            mat_assign!(BxT, sparse(-Bx'))
            mat_assign!(ByT, sparse(-By'))

            # Matrices for borders BCs
            set_boundary_indicator!(grid, grid.LS[end].geoL, grid.LS[end].geoL, op)
            mass_matrix_borders!(num,ind, op.iMx_b, op.iMy_b, op.iMx_bd, op.iMy_bd, grid.LS[end].geoL.dcap, ny)
            bc_matrix_borders!(grid, ind, op.Hx_b, op.Hy_b, grid.LS[end].geoL.dcap)
            mat_assign_T!(op.HxT_b, sparse(op.Hx_b'))
            mat_assign_T!(op.HyT_b, sparse(op.Hy_b'))
            periodic_bcs_borders!(grid, op.Hx_b, op.Hy_b, periodic_x, periodic_y)

           
            # itest=-1
            # if num.io_pdi>0
            #     try
            #         printstyled(color=:magenta, @sprintf "\n PDI write_scalar_transport %.5i \n" num.current_i)
            #         #in YAML file: save only if iscal ==1 for example
            #         PDI_status = @ccall "libpdi".PDI_multi_expose("write_scalar_transport"::Cstring,
            #         "iscal"::Cstring, itest::Ref{Clonglong}, PDI_OUT::Cint,
            #         "rhs_1D"::Cstring, rhs_test::Ptr{Cdouble}, PDI_OUT::Cint,
            #         C_NULL::Ptr{Cvoid})::Cint
            #     catch error
            #         printstyled(color=:red, @sprintf "\n PDI error \n")
            #         print(error)
            #         printstyled(color=:red, @sprintf "\n PDI error \n")
            #     end
            # end #if io_pdi
            
        end

    end #nLS ==1

    #either update here or compute 
    # update_all_ls_data(num, grid, grid_u, grid_v, BC_int, periodic_x, periodic_y)

    for iscal=1:num.nb_transported_scalars

        # reset zero A =0 ? B=0 ?
        # check no need to reinitialize A and B
        a0 .= 0.0
        rhs .= 0.0


        # TODO define rhs from yaml with 
        # eval(Meta.parseall(macros.rhs_scalar))
        # add in macros in yml:
        # printstyled(color=:red, @sprintf "\n rhs modified \n")

        # # eval(Meta.parseall(macros.print_parameters))
        # x_centroid = grid.x .+ getproperty.(grid.LS[1].geoL.centroid, :x) .* grid.dx
        # y_centroid = grid.y .+ getproperty.(grid.LS[1].geoL.centroid, :y) .* grid.dy
        
        # b = Δf.(
        #     x_centroid,
        #     y_centroid,
        #     time
        # )
        # veci(rhs,grid_p,1) .+= op.opC_pL.M * vec(b)
        
        # print("\n rhs ",minimum(rhs)," rhs ",maximum(rhs))



        # A.nzval .= 0.0
        # B.nzval .= 0.0
        A .= 0.0
        B .= 0.0



        # printstyled(color=:red, @sprintf "\n levelset: start iscal!\n")
        # println(grid.LS[1].geoL.dcap[1,1,:])
    
        diffusion_coeff_scal = num.diffusion_coeff[iscal]




        # printstyled(color=:red, @sprintf "\n test diffusion_coeff_scal ")
        # diffusion_coeff_scal = 0.0
        # if iscal == 2
        #     # diffusion_coeff_scal = 0.0
        #     # diffusion_coeff_scal *= 2.0
        #     diffusion_coeff_scal *= 0.5
        # end
        printstyled(color=:red, @sprintf "\n test diffusion_coeff_scal %.2e \n" diffusion_coeff_scal )


        #Interface boundary condition
        if num.nLS ==1

            if is_dirichlet(bc[iscal].int)
                # printstyled(color=:green, @sprintf "\n Dirichlet : %.2e\n" bc[iscal].int.val )
                __a0 = bc[iscal].int.val
                __a1 = -1.0
                __b = 0.0
            elseif is_neumann(bc[iscal].int)
                __a0 = bc[iscal].int.val
                __a1 = 0.0
                __b = 1.0
            elseif is_robin(bc[iscal].int)
                __a0 = bc_type.val
                __a1 = -1.0
                __b = 1.0
            elseif is_stefan(bc[iscal].int)
                __a0 = num.concentration0[iscal]
                __a1 = -1.0
                __b = 0.0
            elseif is_wall(bc[iscal].int)
                __a0 = bc[iscal].int.val
                __a1 = -1.0
                __b = 0.0
            else
                __a0 = bc[iscal].int.val
                __a1 = -1.0
                __b = 0.0
            end

            # Flags with BCs
            a0 .= __a0
            # a0 = ones(grid) .* __a0

            _a1 = ones(grid) .* __a1
            a1 = Diagonal(vec(_a1))
            _b = ones(grid) .* __b
            b = Diagonal(vec(_b))

        end #if num.nLS ==1


        #Wall BC
        a0_b = zeros(nb)
        _a1_b = zeros(nb)
        _b_b = zeros(nb)

        for iLS in 1:num.nLS
            set_borders!(grid, grid.LS[iLS].cl, grid.LS[iLS].u, a0_b, _a1_b, _b_b, bc[iscal], num.n_ext_cl)
        end
        a1_b = Diagonal(vec(_a1_b))
        b_b = Diagonal(vec(_b_b))

        LT = BxT * iMx * Bx .+ ByT * iMy * By
        # LD = BxT * iMx * Hx[1] .+ ByT * iMy * Hy[1]
        LD_b = BxT * op.iMx_b * op.Hx_b .+ ByT * op.iMy_b * op.Hy_b

        if num.scalar_scheme == 0
            time_factor = 0.5 .* num.τ
        elseif num.scalar_scheme == 1
            time_factor = num.τ
        end
      

        if num.nLS ==1
            LD = BxT * iMx * Hx[1] .+ ByT * iMy * Hy[1]

            # Implicit part of heat equation
            A[1:ni,1:ni] = pad_crank_nicolson(M .- time_factor .* diffusion_coeff_scal .* LT, grid, 4 * time_factor)
            A[1:ni,ni+1:2*ni] = - time_factor .* diffusion_coeff_scal .* LD
            A[1:ni,end-nb+1:end] = - time_factor .* diffusion_coeff_scal .* LD_b

            # Interior BC
            A[ni+1:2*ni,1:ni] = b * (HxT[1] * iMx * Bx .+ HyT[1] * iMy * By)
            A[ni+1:2*ni,ni+1:2*ni] = pad(b * (HxT[1] * iMx * Hx[1] .+ HyT[1] * iMy * Hy[1]) .- op.χ[1] * a1)
            A[ni+1:2*ni,end-nb+1:end] = b * (HxT[1] * op.iMx_b * op.Hx_b .+ HyT[1] * op.iMy_b * op.Hy_b)

            # Border BCs
            A[end-nb+1:end,1:ni] = b_b * (op.HxT_b * op.iMx_b' * Bx .+ op.HyT_b * op.iMy_b' * By)
            A[end-nb+1:end,ni+1:2*ni] = b_b * (op.HxT_b * op.iMx_b' * op.Hx[1] .+ op.HyT_b * op.iMy_b' * op.Hy[1])
            A[end-nb+1:end,end-nb+1:end] = pad(b_b * (op.HxT_b * op.iMx_bd * op.Hx_b .+ op.HyT_b * op.iMy_bd * op.Hy_b) .- op.χ_b * a1_b, 4.0)
            # A[end-nb+1:end,end-nb+1:end] = pad(b_b * (op.HxT_b * op.iMx_bd * op.Hx_b .+ op.HyT_b * op.iMy_bd * op.Hy_b) .- op.χ_b * a1_b)

            # Explicit part of heat equation
            B[1:ni,1:ni] = M .+ time_factor .* diffusion_coeff_scal .* LT .- num.τ .* op_conv.CT
            B[1:ni,ni+1:2*ni] = time_factor .* diffusion_coeff_scal .* LD
            B[1:ni,end-nb+1:end] = time_factor .* diffusion_coeff_scal .* LD_b

            vec2(rhs,grid) .+= op.χ[1] * vec(a0)


        else #num.nLS != 1

            # Implicit part of heat equation
            A[1:ni,1:ni] = pad_crank_nicolson(M .- time_factor .* diffusion_coeff_scal .* LT, grid, 4 * time_factor)

            # A[1:ni,ni+1:2*ni] = - time_factor .* diffusion_coeff_scal .* LD
            A[1:ni,end-nb+1:end] = - time_factor .* diffusion_coeff_scal .* LD_b

            # Border BCs
            A[end-nb+1:end,1:ni] = b_b * (op.HxT_b * op.iMx_b' * Bx .+ op.HyT_b * op.iMy_b' * By)
            # A[end-nb+1:end,ni+1:2*ni] = b_b * (op.HxT_b * op.iMx_b' * op.Hx[1] .+ op.HyT_b * op.iMy_b' * op.Hy[1])
            A[end-nb+1:end,end-nb+1:end] = pad(b_b * (op.HxT_b * op.iMx_bd * op.Hx_b .+ op.HyT_b * op.iMy_bd * op.Hy_b) .- op.χ_b * a1_b, 4.0)
            # A[end-nb+1:end,end-nb+1:end] = pad(b_b * (op.HxT_b * op.iMx_bd * op.Hx_b .+ op.HyT_b * op.iMy_bd * op.Hy_b) .- op.χ_b * a1_b)

            # Explicit part of heat equation
            B[1:ni,1:ni] = M .+ time_factor .* diffusion_coeff_scal .* LT .- num.τ .* op_conv.CT
            # B[1:ni,ni+1:2*ni] = time_factor .* diffusion_coeff_scal .* LD
            B[1:ni,end-nb+1:end] = time_factor .* diffusion_coeff_scal .* LD_b

            for iLS in 1:num.nLS

                a0 .= 0.0

                # print("\n BC iLS ",iLS," ",bc[iscal].LS[iLS])
                # print("\n BC iLS val ",bc[iscal].LS[iLS].val)

                if is_dirichlet(bc[iscal].LS[iLS])
                    printstyled(color=:green, @sprintf "\n Dirichlet : %.2e\n" bc[iscal].LS[iLS].val )
                    __a0 = bc[iscal].LS[iLS].val
                    __a1 = -1.0
                    __b = 0.0
                elseif is_neumann(bc[iscal].LS[iLS])
                    printstyled(color=:green, @sprintf "\n Neumann : %.2e\n" bc[iscal].LS[iLS].val )
                    __a0 = bc[iscal].LS[iLS].val
                    __a1 = 0.0
                    __b = 1.0
                elseif is_robin(bc[iscal].LS[iLS])
                    __a0 = bc[iscal].LS[iLS].val
                    __a1 = -1.0
                    __b = 1.0
                elseif is_stefan(bc[iscal].LS[iLS])
                    __a0 = num.concentration0[iscal]
                    __a1 = -1.0
                    __b = 0.0
                elseif is_wall(bc[iscal].LS[iLS])
                    __a0 = bc[iscal].LS[iLS].val
                    __a1 = -1.0
                    __b = 0.0
                else
                    printstyled(color=:green, @sprintf "\n other bc : %.2e\n" bc[iscal].LS[iLS].val )
                    __a0 = bc[iscal].LS[iLS].val
                    __a1 = -1.0
                    __b = 0.0
                end
                
                #BC for LS describing the wall for electrolysis
                #Supposing phi = num.phi_ele1 in wall
                if num.scalar_bc == 1
                    iLS_elec = 2
                    if iLS == iLS_elec #only for second LS
                        # a0 = butler_volmer_no_concentration.(num.alpha_a,num.alpha_c,num.Faraday,num.i0,veci(ph.phi_eleD, grid,iLS_elec+1),
                        # num.phi_ele1,num.Ru,num.temperature0)./(2*num.Faraday*num.diffusion_coeff[iscal])

                        if iscal==1 || iscal==2 #produced: *- (i negative)
                            # a0 = -butler_volmer_no_concentration.(num.alpha_a,num.alpha_c,num.Faraday,num.i0,reshape(veci(ph.phi_eleD, grid,iLS_elec+1),grid),
                            # num.phi_ele1,num.Ru,num.temperature0)./(2*num.Faraday*num.diffusion_coeff[iscal])
                            inv_stoechiometric_coeff = -1.0/2.0 #H2 and KOH
                        
                        elseif iscal == 3
                            #H2O consummed: *+
                            # a0 = +butler_volmer_no_concentration.(num.alpha_a,num.alpha_c,num.Faraday,num.i0,reshape(veci(ph.phi_eleD, grid,iLS_elec+1),grid),
                            # num.phi_ele1,num.Ru,num.temperature0)./(2*num.Faraday*num.diffusion_coeff[iscal])
                            inv_stoechiometric_coeff = 1
                        else 
                            @error("inv_stoechiometric_coeff")
                        end

                        # if iscal ==2
                        #     @error("test inv_stoechiometric_coeff")

                        #     inv_stoechiometric_coeff *= -1
                        # end

                        # if iscal ==2
                        #     @error("test inv_stoechiometric_coeff")

                        #     inv_stoechiometric_coeff *= -1
                        # end
                        # printstyled(color=:red, @sprintf "\n test BC chem \n")

                        # inv_stoechiometric_coeff = 0

                        # printstyled(color=:red, @sprintf "\n test BC chem \n")


                        # butler_volmer_no_concentration_concentration_Neumann!.(num.alpha_a,num.alpha_c,num.Faraday,num.i0,reshape(veci(ph.phi_eleD, grid,iLS_elec+1),grid),
                        # num.phi_ele1,num.Ru,num.temperature0,num.diffusion_coeff[iscal],inv_stoechiometric_coeff,a0)

                        a0 .= butler_volmer_no_concentration_concentration_Neumann.(num.alpha_a,num.alpha_c,num.Faraday,num.i0,
                        reshape(veci(ph.phi_eleD, grid,iLS_elec+1),grid),
                        num.phi_ele1,num.Ru,num.temperature0,num.diffusion_coeff[iscal],inv_stoechiometric_coeff)
                        
                        print("\n a0 iLS_elec min max ",minimum(a0)," ",maximum(a0))


                    else #usual boundary
                        a0 .= __a0

                        print("\n a0 USUAL min max ",minimum(a0)," ",maximum(a0))


                    end #num.scalar_bc == 1

                    printstyled(color=:red, @sprintf "\n phi %.2e %.2e \n" maximum(ph.phi_eleD) minimum(ph.phi_eleD) )

                    printstyled(color=:red, @sprintf "\n BC scalar_transport %.5i %.2e %.2e %.2e \n" iscal maximum(a0) minimum(a0) bc[iscal].LS[iLS].val)

                    # if (maximum(a0) != bc[iscal].LS[iLS].val)
                    #     printstyled(color=:red, @sprintf "\n phi %.2e %.2e \n" iscal maximum(ph.phi_eleD) minimum(ph.phi_eleD) )
 
                    #     # printstyled(color=:red, @sprintf "\n error BC scalar_transport ")
                    #     # print("\n",num.alpha_a," a_c ",num.alpha_c," F ",num.Faraday," i0 ",num.i0," phi0 ",
                    #     # reshape(veci(ph.phi_eleD, grid,iLS_elec+1),grid)
                    #     # ," phi1 ",num.phi_ele1," Ru ",num.Ru," T ",num.temperature0," D ",num.diffusion_coeff[iscal]," st ",inv_stoechiometric_coeff)
                    # end

                else #if num.scalar_bc == 1
                    # Flags with BCs
                    a0 .= __a0
                    # a0 = ones(grid) .* __a0

                end #if num.scalar_bc == 1
    
                _a1 = ones(grid) .* __a1
                a1 = Diagonal(vec(_a1))
                _b = ones(grid) .* __b
                b = Diagonal(vec(_b))


                sb = iLS*ni+1:(iLS+1)*ni
            
                #Here allocate better or replace LD_i
                LD_i = BxT * iMx * Hx[iLS] .+ ByT * iMy * Hy[iLS] 

                # Implicit part of heat equation
                A[1:ni,sb] = - time_factor .* diffusion_coeff_scal .* LD_i

                # Interior BC
                A[sb,1:ni] = b * (HxT[iLS] * iMx * Bx .+ HyT[iLS] * iMy * By)
                
                # Contribution to Neumann BC from other boundaries
                for i in 1:num.nLS
                    if i != iLS
                        A[sb,i*ni+1:(i+1)*ni] = b * (HxT[iLS] * iMx * Hx[i] .+ HyT[iLS] * iMy * Hy[i]) #-b in Poisson
                    end
                end

                A[sb,sb] = pad(b * (HxT[iLS] * iMx * Hx[iLS] .+ HyT[iLS] * iMy * Hy[iLS]) .- op.χ[iLS] * a1)

                #TODO no need for fs_mat

                A[sb,end-nb+1:end] = b * (HxT[iLS] * op.iMx_b * op.Hx_b .+ HyT[iLS] * op.iMy_b * op.Hy_b)

                #Border BCs
                A[end-nb+1:end,sb] = b_b * (op.HxT_b * op.iMx_b' * op.Hx[iLS] .+ op.HyT_b * op.iMy_b' * op.Hy[iLS])

                # Explicit part of scalar equation
                B[1:ni,sb] = time_factor .* diffusion_coeff_scal .* LD_i

                veci(rhs,grid,iLS+1) .+= op.χ[iLS] * vec(a0)

                printstyled(color=:red, @sprintf "\n rhs scalar_transport %.5i %.2e %.2e %.2e \n" iscal maximum(veci(rhs,grid,iLS+1)) minimum(veci(rhs,grid,iLS+1)) bc[iscal].LS[iLS].val)
                print("\n a0 min max ",minimum(a0)," ",maximum(a0))
                print("\n _b min max ",minimum(_b)," ",maximum(_b))

                
            end #for iLS in 1:num.nLS

        end #if num.nLS ==1

        if convection
            vec1(rhs,grid) .-= num.τ .* all_CUTCT[:,iscal]
        end

        vecb(rhs,grid) .+= op.χ_b * vec(a0_b)
        
        if num.convection_Cdivu>0
            # Duv = fzeros(grid)
            # Duv = fnzeros(grid,num)

            #TODO BC missing

            Duv = op.AxT * vec1(ph.uD,grid_u) .+ op.Gx_b * vecb(ph.uD,grid_u) .+
            op.AyT * vec1(ph.vD,grid_v) .+ op.Gy_b * vecb(ph.vD,grid_v)
            for iLS in 1:num.nLS
                if !is_navier(BC_int[iLS]) && !is_navier_cl(BC_int[iLS])
                    Duv .+= op.Gx[iLS] * veci(ph.uD,grid_u,iLS+1) .+ 
                            op.Gy[iLS] * veci(ph.vD,grid_v,iLS+1)
                end
            end
        
            # rhs .+= Duv #.* ph.trans_scalD[:,iscal] multiplied just after

            # vec1(rhs_ϕ,grid) .= rho1 .* iτ .* Duv #TODO
            vec1(rhs,grid) .+= Duv 

            printstyled(color=:green, @sprintf "\n max Duv for C.div(u): %.2e\n" maximum(Duv))

        end #num.convection_Cdivu>0
        
        if num.scalar_scheme == 0
        @views mul!(rhs, B, ph.trans_scalD[:,iscal], 1.0, 1.0) #TODO @views not necessary ?
        end

        #System resolution
        @views ph.trans_scalD[:,iscal] .= A \ rhs
        
        @views ph.trans_scal[:,:,iscal] .= reshape(veci(ph.trans_scalD[:,iscal],grid,1), grid)

        if num.io_pdi>0
            try
                printstyled(color=:magenta, @sprintf "\n PDI write_scalar_transport %.5i \n" num.current_i)
                #in YAML file: save only if iscal ==1 for example
                PDI_status = @ccall "libpdi".PDI_multi_expose("write_scalar_transport"::Cstring,
                "iscal"::Cstring, iscal::Ref{Clonglong}, PDI_OUT::Cint,
                "rhs_1D"::Cstring, rhs::Ptr{Cdouble}, PDI_OUT::Cint,
                C_NULL::Ptr{Cvoid})::Cint
            catch error
                printstyled(color=:red, @sprintf "\n PDI error \n")
                print(error)
                printstyled(color=:red, @sprintf "\n PDI error \n")
            end
        end #if io_pdi



        #Checks after resolutions: is the value physical?
        iLS = 1
        nonzero = mean_intfc_non_null(ph.trans_scalD,iscal,grid,iLS) #Value at interface
        printstyled(color=:green, @sprintf "\n mean  interface : %.2e\n" nonzero)

        if num.nLS>1
            iLS = 2
            nonzero = mean_intfc_non_null(ph.trans_scalD,iscal,grid,iLS) #Value at interface
            printstyled(color=:green, @sprintf "\n mean  wall LS : %.2e\n" nonzero)

        end

        if iscal!=3 #H2O consummed at the electrode, would need to make distinction to make sure the decrease in H2O is physical or not
          
            @views kill_dead_cells_val!(ph.trans_scal[:,:,iscal], grid, LS[1].geoL,num.concentration0[iscal]) 

            # @views veci(ph.trans_scalD[:,iscal],grid,1) .= vec(ph.trans_scal[:,:,iscal])

            min_border_bulk = min(minimum(ph.trans_scal[:,:,iscal]),minimum(vecb(ph.trans_scal[:,:,iscal],grid)))

            if min_border_bulk.<num.concentration0[iscal]*(1-num.concentration_check_factor)
                print("iscal ",iscal)
                printstyled(color=:red, @sprintf "\n concentration: %.10e %.10e \n" min_border_bulk num.concentration0[iscal]*(1-num.concentration_check_factor))
                # printstyled(color=:red, @sprintf "\n concentration drop: %.2e%% \n" (minimum(ph.trans_scal[:,:,iscal])-num.concentration0[iscal])/num.concentration0[iscal]*100)
                @error("concentration too low")

               
            end #too low

            printstyled(color=:red, @sprintf "\n concentration variation vs min: %.2e%% \n" (minimum(ph.trans_scal[:,:,iscal])-num.concentration0[iscal])/num.concentration0[iscal]*100)

        else 
            @views kill_dead_cells_val!(ph.trans_scal[:,:,iscal], grid, LS[1].geoL,num.concentration0[iscal]) 

            max_border_bulk = max(maximum(ph.trans_scal[:,:,iscal]),maximum(vecb(ph.trans_scal[:,:,iscal],grid)))
            if max_border_bulk.> num.concentration0[iscal]*(1+num.concentration_check_factor)
                print("iscal ",iscal)
                printstyled(color=:red, @sprintf "\n concentration: %.10e %.10e \n" max_border_bulk num.concentration0[iscal]*(1-num.concentration_check_factor))
                # printstyled(color=:red, @sprintf "\n concentration increase: %.2e%% \n" (maximum(ph.trans_scal[:,:,iscal])-num.concentration0[iscal])/num.concentration0[iscal]*100)
                @error("concentration too high")

            end #too high
            printstyled(color=:red, @sprintf "\n concentration variation vs max: %.2e%% \n" (maximum(ph.trans_scal[:,:,iscal])-num.concentration0[iscal])/num.concentration0[iscal]*100)
        end

        @views kill_dead_cells_val!(ph.trans_scal[:,:,iscal], grid, LS[1].geoL,0.0)  #reset to zero (for plot)
        #kill_dead_cells_val! can be used for display to avoid displaying a range from 0 to c0 in python

        print("\n test concentration ", minimum(ph.trans_scal[:,:,iscal]), " ", maximum(ph.trans_scal[:,:,iscal]))

    end #end loop iscal

    return nothing
end






# """
# From Stefan_velocity!
# """
# function electrolysis_velocity!(num, grid, LS, V, TL, MIXED, periodic_x, periodic_y, concentration_scal_intfc,electrolysis_phase_change_case, mass_flux)
#     @unpack geoS, geoL, κ = LS

#     V .= 0

#     if electrolysis_phase_change_case == "levelset"

#         # @inbounds @threads for II in MIXED
#         @inbounds for II in MIXED

#             θ_d = concentration_scal_intfc

#             # dTS = 0.
#             dTL = 0.
#             if geoL.projection[II].flag
#                 T_1, T_2 = interpolated_temperature(grid, geoL.projection[II].angle, geoL.projection[II].point1, geoL.projection[II].point2, TL, II, periodic_x, periodic_y)
#                 dTL = normal_gradient(geoL.projection[II].d1, geoL.projection[II].d2, T_1, T_2, θ_d)
#             else
#                 T_1 = interpolated_temperature(grid, geoL.projection[II].angle, geoL.projection[II].point1, TL, II, periodic_x, periodic_y)
#                 dTL = normal_gradient(geoL.projection[II].d1, T_1, θ_d)
#             end
#             V[II] = dTL #+ dTS

#             print("\n grad",II,"val",dTL,"val",geoL.projection[II].flag,"val",geoL.projection[II].angle,"val", geoL.projection[II].point1,"val", geoL.projection[II].point2,"val", T_1,"val",T_2,"val",θ_d)
#         end

#     else

#         intfc_length = 0.0
#         @inbounds for II in MIXED
#             V[II] = sum(mass_flux)

#             #χx = (geo.dcap[:,:,3] .- geo.dcap[:,:,1]) .^ 2
#             #χy = (geo.dcap[:,:,4] .- geo.dcap[:,:,2]) .^ 2
#             #χ[iLS].diag .= sqrt.(vec(χx .+ χy))

#             χx = (geo.dcap[II,3] .- geo.dcap[II,1]) .^ 2
#             χy = (geo.dcap[II,4] .- geo.dcap[II,2]) .^ 2

#             intfc_length += sqrt.(vec(χx .+ χy))
#         end 

#         V./= intfc_length

#         print("\n phase-change velocity ", sum(mass_flux)/intfc_length)

#         printstyled(color=:magenta, @sprintf "\n phase-change velocity %.2e intfc_length %.2e πR %.2e\n" sum(mass_flux)/intfc_length intfc_length π*num.R)

#         i_ext, l_ext, b_ext, r_ext, t_ext = indices_extension(grid, grid.LS[iLS], grid.ind.inside, periodic_x, periodic_y)
#         field_extension!(grid, u, grid.V, i_ext, l_ext, b_ext, r_ext, t_ext, num.NB, periodic_x, periodic_y)
        

#     end

#     return nothing
# end

"""
From update_free_surface_velocity and update_stefan_velocity
"""
function update_free_surface_velocity_electrolysis!(num, grid, grid_u, grid_v, iLS, uD, vD, 
    periodic_x, periodic_y, average_velocity, concentration_scalD, concentration_scal, diffusion_coeff_scal,concentration_scal_intfc, 
    electrolysis_phase_change_case,mass_flux)

    grid.V .= 0
    num.sum_mass_flux = 0.0
    v_mean = 0.0

    factor = -(1.0/num.rho2-1.0/num.rho1).*diffusion_coeff_scal[1].*num.MWH2

    intfc_length = 0.0

    num_mixed_cells = 0

    @inbounds for II in grid.LS[iLS].MIXED
        # print("\n II update ",II, grid.LS[end].u[II], " iso end ",grid.LS[end].iso[II]," iso 1 ",grid.LS[1].iso[II])
        if grid.LS[end].iso[II] < 14.5 #15.0 -0.5 # check if inside domain defined by other LS 
        # if grid.LS[end].u[II]>0.0 # check if inside domain defined by other LS 
        # if grid.LS[2].u[II]>0.0 #second wall
            # print("\n cells for free surface", II," x ",grid.x[II]," LS[iLS] ",grid.LS[iLS].u[II]," LS[end] ",grid.LS[end].u[II]," LS[2] ",grid.LS[2].u[II])
            # grid.V[II] = mass_flux[II] * factor

            num_mixed_cells += 1

            if num.mass_flux == 0
                num.sum_mass_flux += mass_flux[II]
                #compute interface length
                χx = (grid.LS[iLS].geoL.dcap[II,3] .- grid.LS[iLS].geoL.dcap[II,1]) .^ 2
                χy = (grid.LS[iLS].geoL.dcap[II,4] .- grid.LS[iLS].geoL.dcap[II,2]) .^ 2
                intfc_length_cell = sqrt(χx + χy)
                intfc_length += intfc_length_cell

                # intfc_length_cell !=0 since mixed cell

                print("\n intfc_length_cell ", intfc_length_cell)

                grid.V[II] = mass_flux[II] * factor / intfc_length_cell
            elseif num.mass_flux == 1

                #TODO iLS or end grid.LS[iLS].geoL

                dTL = 0.0
                # print("\n II ",II," flag ",grid.LS[iLS].geoL.projection[II].flag)
                if grid.LS[iLS].geoL.projection[II].flag
                    T_1, T_2 = interpolated_temperature(grid, grid.LS[iLS].geoL.projection[II].angle, grid.LS[iLS].geoL.projection[II].point1, grid.LS[iLS].geoL.projection[II].point2, concentration_scal, II, periodic_x, periodic_y)
                    dTL = normal_gradient(grid.LS[iLS].geoL.projection[II].d1, grid.LS[iLS].geoL.projection[II].d2, T_1, T_2, concentration_scal_intfc)
                    printstyled(color=:cyan, @sprintf "\n T1 %.2e T2 %.2e \n" T_1 T_2 )
                    if isnan(T_2)
                        printstyled(color=:red, @sprintf "\n T2 NaN, resorting to other method \n")
                        print("\n P2 ",grid.LS[iLS].geoL.projection[II].point2)


                        T_1 = interpolated_temperature(grid, grid.LS[iLS].geoL.projection[II].angle, grid.LS[iLS].geoL.projection[II].point1, concentration_scal, II, periodic_x, periodic_y)
                        dTL = normal_gradient(grid.LS[iLS].geoL.projection[II].d1, T_1, concentration_scal_intfc)
                    end

                    if isnan(T_1) || isnan(T_2) #debug
                        printstyled(color=:red, @sprintf "\n T1 or T2 NaN, debug \n")

                        print("\n II ",II," flag ",grid.LS[iLS].geoL.projection[II].flag)

                        vtx_num = 2                
                        vtx_x = [grid.LS[iLS].geoL.projection[II].point1.x,grid.LS[iLS].geoL.projection[II].point2.x]
                        vtx_y = [grid.LS[iLS].geoL.projection[II].point1.y,grid.LS[iLS].geoL.projection[II].point2.y]

                        PDI_status = @ccall "libpdi".PDI_multi_expose("debug_phase_change"::Cstring,
                        "vtx_num"::Cstring, vtx_num::Ref{Clonglong}, PDI_OUT::Cint, 
                        "vtx_x"::Cstring, vtx_x::Ptr{Cdouble}, PDI_OUT::Cint,
                        "vtx_y"::Cstring, vtx_y::Ptr{Cdouble}, PDI_OUT::Cint,
                        C_NULL::Ptr{Cvoid})::Cint

                        return 1
                        
                    end

                else
                    T_1 = interpolated_temperature(grid, grid.LS[iLS].geoL.projection[II].angle, grid.LS[iLS].geoL.projection[II].point1, concentration_scal, II, periodic_x, periodic_y)
                    dTL = normal_gradient(grid.LS[iLS].geoL.projection[II].d1, T_1, concentration_scal_intfc)
                end
                # grid.V[II] = dTL #+ dTS
                printstyled(color=:cyan, @sprintf "\n v %.2e v from int %.2e %.2e %.2e\n" grid.V[II] dTL*factor T_1 concentration_scal_intfc)
                
                grid.V[II] = dTL*factor 
                v_mean += grid.V[II] #TODO unit use same
                # printstyled(color=:red, @sprintf "\n TODO unit use same \n" )
            
            # elseif num.mass_flux == 2
                # χx = (grid.LS[iLS].geoL.dcap[II,3] .- grid.LS[iLS].geoL.dcap[II,1]) .^ 2
                # χy = (grid.LS[iLS].geoL.dcap[II,4] .- grid.LS[iLS].geoL.dcap[II,2]) .^ 2
                # intfc_length_cell = sqrt(χx + χy)
                # intfc_length += intfc_length_cell

            end #num.mass_flux
        
        end #grid.LS[end].iso[II] != 15.0
    end 
    if average_velocity == 1
        if num.mass_flux == 0
            v_mean = factor * num.sum_mass_flux /intfc_length
        elseif num.mass_flux == 1 
            v_mean = v_mean / num_mixed_cells
            
            if num.current_i == 1 && num.advection_LS_mode == 9
                @test v_mean ≈ 6.889685460499036e-5 atol=1e-12 #6.89e-05

                if abs(v_mean-6.889685460499036e-5) < 1e-12
                @error("error phase-change velocity")
                printstyled(color=:red, @sprintf "\n error velocity\n")
                return 1
                end
            end

        end 

        if v_mean < 0.0
            @error("error phase-change velocity")
            printstyled(color=:red, @sprintf "\n error velocity\n")
            return 1
        end

        @inbounds for II in grid.LS[iLS].MIXED
            grid.V[II] = v_mean
        end 
        
        cfl_tmp = v_mean*num.dt0/grid.dx[1,1]
        printstyled(color=:cyan, @sprintf "\n v_mean %.2e CFL %.2e dt %.2e dx %.2e\n" v_mean cfl_tmp num.dt0 grid.dx[1,1])

        if cfl_tmp > num.CFL
            @error("Error phase-change CFL ")
            return 1
        end
        # else
    #     grid.V ./= intfc_length

    end


    # printstyled(color=:magenta, @sprintf "\n test sign velocity\n" )
    # grid.V.*=-1.0

    # printstyled(color=:cyan, @sprintf "\n flux %.2e factor %.2e intfc_length %.2e\n" num.sum_mass_flux factor intfc_length)

    # printstyled(color=:magenta, @sprintf "\n sum_intfc %.2e sum_intfc/intfc_length %.2e sum all cells %.2e \n" num.sum_mass_flux num.sum_mass_flux/intfc_length sum(mass_flux))

    # printstyled(color=:red, @sprintf "\n test phase-change velocity %.2e intfc_length %.2e πR %.2e\n" sum(mass_flux)*factor/intfc_length intfc_length π*num.R)
    # printstyled(color=:magenta, @sprintf "\n phase-change velocity %.2e intfc_length %.2e πR %.2e\n" num.sum_mass_flux*factor/intfc_length intfc_length π*num.R)

    # printstyled(color=:magenta, @sprintf "\n intfc_length %.2e πR %.2e\n" intfc_length π*num.R)

    if num.extend_field == 0
        i_ext, l_ext, b_ext, r_ext, t_ext = indices_extension(grid, grid.LS[iLS], grid.ind.inside, periodic_x, periodic_y)
        field_extension!(grid, grid.LS[iLS].u, grid.V, i_ext, l_ext, b_ext, r_ext, t_ext, num.NB, periodic_x, periodic_y)
    
    elseif num.extend_field == 1 #constant velocity everywhere
        grid.V .= v_mean
    end


    printstyled(color=:green, @sprintf "\n grid p u v max : %.2e %.2e %.2e\n" maximum(abs.(grid.V[grid.LS[iLS].MIXED])) maximum(abs.(grid_u.V[grid.LS[iLS].MIXED])) maximum(abs.(grid_v.V[grid_v.LS[iLS].MIXED])))

    #plot_electrolysis_velocity!(num, grid, grid.LS[iLS], grid.V, concentration_scalD, grid.LS[iLS].MIXED, periodic_x, periodic_y, concentration_scal_intfc)

    # electrolysis_velocity!(num, grid, grid.LS[iLS], grid.V, concentration_scalD, grid.LS[iLS].MIXED, periodic_x, periodic_y, concentration_scal_intfc,electrolysis_phase_change_case, mass_flux)
    

    # concentration_scalu = zeros(grid_u)
    # concentration_scalv = zeros(grid_v)
    # #interpolate
    # interpolate_scalar!(grid, grid_u, grid_v, reshape(veci(concentration_scalD,grid,1), grid), concentration_scalu, concentration_scalv)
    # # Or give scal directly instead of scalD 
    # # interpolate_scalar!(grid, grid_u, grid_v, reshape(veci(concentration_scal,grid,1), grid), concentration_scalu, concentration_scalv)
    # #TODO H: suppose val on interface for u and v same as for p 
    # electrolysis_velocity!(num, grid_u, grid_u.LS[iLS], grid_u.V, concentration_scalu, grid_u.LS[iLS].MIXED, periodic_x, periodic_y, concentration_scal_intfc)
    # electrolysis_velocity!(num, grid_v, grid_v.LS[iLS], grid_v.V, concentration_scalv, grid_v.LS[iLS].MIXED, periodic_x, periodic_y, concentration_scal_intfc)
    # grid.V[grid.LS[iLS].MIXED] .*= 1. ./ λ
    #H2
    #TODO check order rho1 rho2
    #no surface term

    
    # grid_u.V[grid_u.LS[iLS].MIXED] .*= factor
    # grid_v.V[grid_v.LS[iLS].MIXED] .*= factor

    # grid.V[grid.LS[iLS].MIXED] .*= factor

   
    # if average_velocity
    #     a = mean(grid.V[grid.LS[iLS].MIXED])
    #     grid.V[grid.LS[iLS].MIXED] .= a
    # end

    # grid_u.V .= reshape(vec1(uD,grid_u), grid_u)
    # grid_v.V .= reshape(vec1(vD,grid_v), grid_v)


    # printstyled(color=:green, @sprintf "\n grid p u v max : %.2e %.2e %.2e\n" maximum(abs.(grid.V[grid.LS[iLS].MIXED])) maximum(abs.(grid_u.V[grid.LS[iLS].MIXED])) maximum(abs.(grid_v.V[grid_v.LS[iLS].MIXED])))

    #TODO extension


#     grid_u.V .+.= reshape(veci(uD,grid_u,iLS+1), (grid_u.ny, grid_u.nx))
#     grid_v.V .+.= reshape(veci(vD,grid_v,iLS+1), (grid_v.ny, grid_v.nx))

#     i_u_ext, l_u_ext, b_u_ext, r_u_ext, t_u_ext = indices_extension(grid_u, grid_u.LS[iLS], grid_u.ind.inside, periodic_x, periodic_y)
#     i_v_ext, l_v_ext, b_v_ext, r_v_ext, t_v_ext = indices_extension(grid_v, grid_v.LS[iLS], grid_v.ind.inside, periodic_x, periodic_y)

#     field_extension!(grid_u, grid_u.LS[iLS].u, grid_u.V, i_u_ext, l_u_ext, b_u_ext, r_u_ext, t_u_ext, num.NB, periodic_x, periodic_y)
#     field_extension!(grid_v, grid_v.LS[iLS].u, grid_v.V, i_v_ext, l_v_ext, b_v_ext, r_v_ext, t_v_ext, num.NB, periodic_x, periodic_y)

    # if num.io_pdi>0
        
    #     try
    #         nstep = num.current_i
    #         printstyled(color=:magenta, @sprintf "\n PDI write_iso %.5i \n" num.current_i)

    #         PDI_status = @ccall "libpdi".PDI_multi_expose("write_iso"::Cstring,
    #         "nstep"::Cstring, nstep::Ref{Clonglong}, PDI_OUT::Cint,
    #         "levelset_iso_end"::Cstring, grid.LS[iLSpdi].iso::Ptr{Cdouble}, PDI_OUT::Cint,
    #         C_NULL::Ptr{Cvoid})::Cint

    #     catch error
    #         printstyled(color=:red, @sprintf "\n PDI error \n")
    #         print(error)
    #     end 
    # end #if num.io_pdi>0

    # print("\n sum_mass_flux ",sum_mass_flux)

    return 0
end


"""Interpolate velocity on scalar grid for regular grids for vizualisation

# Arguments
- `grid`: scalar grid
eps parameter TODO 
"""
function interpolate_grid_liquid!( grid::Mesh{Flower.GridCC, Float64, Int64},
    grid_u::Mesh{Flower.GridFCx, Float64, Int64},
    grid_v::Mesh{Flower.GridFCy, Float64, Int64},
    u::Array{Float64, 2},
    v::Array{Float64, 2},
    us::Array{Float64, 2},
    vs::Array{Float64, 2})
    
    #TODO no need to reset?
    # us = p .*0
    # vs = p .*0
    LS_u =grid_u.LS[1]
    LS_v = grid_v.LS[1]
    us .= (
        (u[:,2:end] .* LS_u.geoL.dcap[:,2:end,6] .+ 
        u[:,1:end-1] .* LS_u.geoL.dcap[:,1:end-1,6]) ./ 
        (LS_u.geoL.dcap[:,1:end-1,6] .+ LS_u.geoL.dcap[:,2:end,6] .+ 1e-8 )
    )
    vs .= (
        (v[2:end,:] .* LS_v.geoL.dcap[2:end,:,7] .+ 
        v[1:end-1,:] .* LS_v.geoL.dcap[1:end-1,:,7]) ./
        (LS_v.geoL.dcap[1:end-1,:,7] .+ LS_v.geoL.dcap[2:end,:,7] .+ 1e-8 )
    )

    # phL.p .+= (
    #     (phS.u[:,2:end].^2.0 .* LS_u.geoS.dcap[:,2:end,6] .+ 
    #     phS.u[:,1:end-1].^2.0 .* LS_u.geoS.dcap[:,1:end-1,6]) ./ 
    #     (LS_u.geoS.dcap[:,1:end-1,6] .+ LS_u.geoS.dcap[:,2:end,6] .+ 1e-8 )
    # )
    # phL.p .+= (
    #     (phS.v[2:end,:].^2.0 .* LS_v.geoS.dcap[2:end,:,7] .+ 
    #     phS.v[1:end-1,:].^2.0 .* LS_v.geoS.dcap[1:end-1,:,7]) ./
    #     (LS_v.geoS.dcap[1:end-1,:,7] .+ LS_v.geoS.dcap[2:end,:,7] .+ 1e-8 )
    # )


    # for j = 1:gp.ny
    # for i = 1:gp.nx
    #     us[:,j,i]=(u[:,j,i]+u[:,j,i+1])/2
    #     vs[:,j,i]=(v[:,j,i]+v[:,j+1,i])/2
    # end
    # end
end


"""Interpolate velocity on scalar grid for regular grids for vizualisation

# Arguments
- `grid`: scalar grid
eps parameter TODO 
"""
function interpolate_grid_solid!( grid::Mesh{Flower.GridCC, Float64, Int64},
    grid_u::Mesh{Flower.GridFCx, Float64, Int64},
    grid_v::Mesh{Flower.GridFCy, Float64, Int64},
    u::Array{Float64, 2},
    v::Array{Float64, 2},
    us::Array{Float64, 2},
    vs::Array{Float64, 2})
    
    #TODO no need to reset?
    # us = p .*0
    # vs = p .*0
    LS_u =grid_u.LS[1]
    LS_v = grid_v.LS[1]
    us .= (
        (u[:,2:end] .* LS_u.geoS.dcap[:,2:end,6] .+ 
        u[:,1:end-1] .* LS_u.geoS.dcap[:,1:end-1,6]) ./ 
        (LS_u.geoS.dcap[:,1:end-1,6] .+ LS_u.geoS.dcap[:,2:end,6] .+ 1e-8 )
    )
    vs .= (
        (v[2:end,:] .* LS_v.geoS.dcap[2:end,:,7] .+ 
        v[1:end-1,:] .* LS_v.geoS.dcap[1:end-1,:,7]) ./
        (LS_v.geoS.dcap[1:end-1,:,7] .+ LS_v.geoS.dcap[2:end,:,7] .+ 1e-8 )
    )

    # phL.p .+= (
    #     (phS.u[:,2:end].^2.0 .* LS_u.geoS.dcap[:,2:end,6] .+ 
    #     phS.u[:,1:end-1].^2.0 .* LS_u.geoS.dcap[:,1:end-1,6]) ./ 
    #     (LS_u.geoS.dcap[:,1:end-1,6] .+ LS_u.geoS.dcap[:,2:end,6] .+ 1e-8 )
    # )
    # phL.p .+= (
    #     (phS.v[2:end,:].^2.0 .* LS_v.geoS.dcap[2:end,:,7] .+ 
    #     phS.v[1:end-1,:].^2.0 .* LS_v.geoS.dcap[1:end-1,:,7]) ./
    #     (LS_v.geoS.dcap[1:end-1,:,7] .+ LS_v.geoS.dcap[2:end,:,7] .+ 1e-8 )
    # )


    # for j = 1:gp.ny
    # for i = 1:gp.nx
    #     us[:,j,i]=(u[:,j,i]+u[:,j,i+1])/2
    #     vs[:,j,i]=(v[:,j,i]+v[:,j+1,i])/2
    # end
    # end
end



"""Interpolate velocity on scalar grid for regular grids for vizualisation

# Arguments
- `grid`: scalar grid
eps parameter TODO 
"""
function interpolate_grid_liquid_solid!(num::Numerical{Float64, Int64}, grid::Mesh{Flower.GridCC, Float64, Int64},
    LS_u::Levelset{Float64, Int64},
    LS_v::Levelset{Float64, Int64},
    u::AbstractArray{Float64, 2},
    v::AbstractArray{Float64, 2},
    us::Array{Float64, 2},
    vs::Array{Float64, 2})
    
    #TODO no need to reset?
    us .= 0.0
    vs .= 0.0

    for j in 1:grid.ny
        for i in 1:grid.nx

            us[j,i] = (u[j,i+1] * LS_u.geoL.dcap[j,i+1,6] + u[j,i] * LS_u.geoL.dcap[j,i,6] + u[j,i+1] * LS_u.geoS.dcap[j,i+1,6] + u[j,i] * LS_u.geoS.dcap[j,i,6] ) / (LS_u.geoL.dcap[j,i+1,6] + LS_u.geoL.dcap[j,i,6] + LS_u.geoS.dcap[j,i+1,6] + LS_u.geoS.dcap[j,i,6])

            if us[j,i] !==1.0
                print("\n j",j," i ",i, " us[j,i] ",us[j,i])
            end


            vs[j,i] = (v[j+1,i] * LS_v.geoL.dcap[j+1,i,7] + v[j,i] * LS_v.geoL.dcap[j,i,7] + v[j+1,i] * LS_v.geoS.dcap[j+1,i,7] + v[j,i] * LS_v.geoS.dcap[j,i,7] ) / (LS_v.geoL.dcap[j+1,i,7] + LS_v.geoL.dcap[j,i,7] + LS_v.geoS.dcap[j+1,i,7] + LS_v.geoS.dcap[j,i,7])

            if vs[j,i] !==1.0
                print("\n j",j," i ",i, " us[j,i] ",vs[j,i])
            end

        end
    end 

   
end


"""Interpolate velocity on scalar grid for regular grids for vizualisation

# Arguments
- `grid`: scalar grid
eps parameter TODO 
"""
function interpolate_grid_liquid_2!(num::Numerical{Float64, Int64}, grid::Mesh{Flower.GridCC, Float64, Int64},
    LS_u::Levelset{Float64, Int64},
    LS_v::Levelset{Float64, Int64},
    u::AbstractArray{Float64, 2},
    v::AbstractArray{Float64, 2},
    us::Array{Float64, 2},
    vs::Array{Float64, 2})
    
    #TODO no need to reset?
    us .= 0.0
    vs .= 0.0

    for j in 1:grid.ny
        for i in 1:grid.nx

            if LS_u.geoL.dcap[j,i,5] >num.epsilon_vol
                us[j,i] = (u[j,i+1] * LS_u.geoL.dcap[j,i+1,6] + u[j,i] * LS_u.geoL.dcap[j,i,6]) / (LS_u.geoL.dcap[j,i+1,6] + LS_u.geoL.dcap[j,i,6])

                # vs[j,i] = (v[j+1,i] * LS_v.geoL.dcap[j+1,i,7] + v[j,i] * LS_v.geoL.dcap[j,i,7]) / (LS_v.geoL.dcap[j+1,i,7] + LS_v.geoL.dcap[j,i,7])

            end

            if LS_v.geoL.dcap[j,i,5] >num.epsilon_vol
                # us[j,i] = (u[j,i+1] * LS_u.geoL.dcap[j,i+1,6] + u[j,i] * LS_u.geoL.dcap[j,i,6]) / (LS_u.geoL.dcap[j,i+1,6] + LS_u.geoL.dcap[j,i,6])

                vs[j,i] = (v[j+1,i] * LS_v.geoL.dcap[j+1,i,7] + v[j,i] * LS_v.geoL.dcap[j,i,7]) / (LS_v.geoL.dcap[j+1,i,7] + LS_v.geoL.dcap[j,i,7])

            end

            # if i==64
            #     print("\n j ",j," i ",i, " LS_u.geoL.dcap[j,i,5] ",LS_u.geoL.dcap[j,i,5]," LS_u.geoL.dcap[j,i+1,6] ",LS_u.geoL.dcap[j,i+1,6], " us ",us[j,i]," u ",u[j,i]," u ",u[j,i+1]    )
            # end

            # if (tmp_vec_p[j,i] != 1.0) 
            #     print("\n j",j," i ",i, " tmp_vec_p ",tmp_vec_p[j,i]," tmp_vec_p0 ",tmp_vec_p0[j,i])
            #     print("\n LS_u.geoL.dcap[:,2:end,6] ",gu.LS[1].geoL.dcap[j,i,6]," ",gu.LS[1].geoL.dcap[j,i+1,6])
            # end

        end
    end 

    # us .= (
    #     (u[:,2:end] .* LS_u.geoL.dcap[:,2:end,6] .+ 
    #     u[:,1:end-1] .* LS_u.geoL.dcap[:,1:end-1,6]) ./ 
    #     (LS_u.geoL.dcap[:,1:end-1,6] .+ LS_u.geoL.dcap[:,2:end,6])
    # )
    # vs .= (
    #     (v[2:end,:] .* LS_v.geoL.dcap[2:end,:,7] .+ 
    #     v[1:end-1,:] .* LS_v.geoL.dcap[1:end-1,:,7]) ./
    #     (LS_v.geoL.dcap[1:end-1,:,7] .+ LS_v.geoL.dcap[2:end,:,7])
    # )

    # phL.p .+= (
    #     (phS.u[:,2:end].^2.0 .* LS_u.geoS.dcap[:,2:end,6] .+ 
    #     phS.u[:,1:end-1].^2.0 .* LS_u.geoS.dcap[:,1:end-1,6]) ./ 
    #     (LS_u.geoS.dcap[:,1:end-1,6] .+ LS_u.geoS.dcap[:,2:end,6] .+ 1e-8 )
    # )
    # phL.p .+= (
    #     (phS.v[2:end,:].^2.0 .* LS_v.geoS.dcap[2:end,:,7] .+ 
    #     phS.v[1:end-1,:].^2.0 .* LS_v.geoS.dcap[1:end-1,:,7]) ./
    #     (LS_v.geoS.dcap[1:end-1,:,7] .+ LS_v.geoS.dcap[2:end,:,7] .+ 1e-8 )
    # )


    # for j = 1:gp.ny
    # for i = 1:gp.nx
    #     us[:,j,i]=(u[:,j,i]+u[:,j,i+1])/2
    #     vs[:,j,i]=(v[:,j,i]+v[:,j+1,i])/2
    # end
    # end
end

"""Interpolate velocity on scalar grid for regular grids for vizualisation

# Arguments
- `grid`: scalar grid
eps parameter TODO 
"""
function interpolate_grid_solid_2!(num::Numerical{Float64, Int64}, grid::Mesh{Flower.GridCC, Float64, Int64},
    LS_u::Levelset{Float64, Int64},
    LS_v::Levelset{Float64, Int64},
    u::AbstractArray{Float64, 2},
    v::AbstractArray{Float64, 2},
    us::Array{Float64, 2},
    vs::Array{Float64, 2})
    
    #TODO no need to reset?
    # us .= 0.0
    # vs .= 0.0

    for j in 1:grid.ny
        for i in 1:grid.nx

            if LS_u.geoS.dcap[j,i,5] >num.epsilon_vol
                us[j,i] = (u[j,i+1] * LS_u.geoS.dcap[j,i+1,6] + u[j,i] * LS_u.geoS.dcap[j,i,6]) / (LS_u.geoS.dcap[j,i+1,6] + LS_u.geoS.dcap[j,i,6])

                # vs[j,i] = (v[j+1,i] * LS_v.geoS.dcap[j+1,i,7] + v[j,i] * LS_v.geoS.dcap[j,i,7]) / (LS_v.geoS.dcap[j+1,i,7] + LS_v.geoS.dcap[j,i,7])

            end

            if LS_v.geoS.dcap[j,i,5] >num.epsilon_vol
                # us[j,i] = (u[j,i+1] * LS_u.geoS.dcap[j,i+1,6] + u[j,i] * LS_u.geoS.dcap[j,i,6]) / (LS_u.geoS.dcap[j,i+1,6] + LS_u.geoS.dcap[j,i,6])

                vs[j,i] = (v[j+1,i] * LS_v.geoS.dcap[j+1,i,7] + v[j,i] * LS_v.geoS.dcap[j,i,7]) / (LS_v.geoS.dcap[j+1,i,7] + LS_v.geoS.dcap[j,i,7])

            end

        end
    end 

    
end


"""
Interpolate velocity on scalar grid for regular grids for vizualisation
"""
function interpolate_regular_grid(grid,fwdL)
    
    us = fwdL.p .*0
    vs = fwdL.p .*0

    for j = 1:grid.ny
    for i = 1:grid.nx
        us[:,j,i]=(fwdL.u[:,j,i]+fwdL.u[:,j,i+1])/2
        vs[:,j,i]=(fwdL.v[:,j,i]+fwdL.v[:,j+1,i])/2
    end
    end
    
    return us,vs
end


"""
Interpolate velocity on scalar grid for regular grids for vizualisation
"""
function interpolate_regular_grid(grid,fwdL,u,v)
    
    us = fwdL.p .*0
    vs = fwdL.p .*0

    for j = 1:grid.ny
    for i = 1:grid.nx
        us[:,j,i]=(u[:,j,i]+u[:,j,i+1])/2
        vs[:,j,i]=(v[:,j,i]+v[:,j+1,i])/2
    end
    end
    
    return us,vs
end

# """
# From kinetic_energy
# """
# function scal_magnitude(ph, gp, gu, gv)
#     ph.p = zeros(gp)

#     ph.p .= (
#         (fwdL.u[i,:,2:end].^2.0 .* gu.geoL.dcap[:,2:end,6] .+ 
#         fwdL.u[i,:,1:end-1].^2.0 .* gu.geoL.dcap[:,1:end-1,6]) ./ 
#         (gu.geoL.dcap[:,1:end-1,6] .+ gu.geoL.dcap[:,2:end,6] .+ 1e-8)
#     )
#     ph.p .+= (
#         (fwdL.v[i,2:end,:].^2.0 .* gv.geoL.dcap[2:end,:,7] .+ 
#         fwdL.v[i,1:end-1,:].^2.0 .* gv.geoL.dcap[1:end-1,:,7]) ./
#         (gv.geoL.dcap[1:end-1,:,7] .+ gv.geoL.dcap[2:end,:,7] .+ 1e-8)
#     )
#     ph.p .= sqrt.(ph.p .* gp.geoL.dcap[:,:,5])

# end


# nb_transported_scalars::Int64,
function print_electrolysis_statistics(num::Numerical{Float64, Int64},
    grid::Mesh{Flower.GridCC, Float64, Int64}, phL::Phase{Float64})
    # normscal1L = norm(phL.trans_scal[:,:,1])
    # normscal2L = norm(phL.trans_scal[:,:,2])
    # normscal3L = norm(phL.trans_scal[:,:,3])
    # normphi_eleL = norm(phL.phi_ele)
    # normTL = norm(phL.T)
    # normiL=0.0 #TODO

    

    minscal1L=minscal2L=minscal3L=minphi_eleL=minTL=miniL=0.0
    maxscal1L=maxscal2L=maxscal3L=maxphi_eleL=maxTL=maxiL=0.0
    moyscal1L=moyscal2L=moyscal3L=moyphi_eleL=moyTL=moyiL=0.0


    if num.nb_transported_scalars>0
        minscal1L = minimum(phL.trans_scalD[:,1])
        maxscal1L = maximum(phL.trans_scalD[:,1])
        moyscal1L = mean(phL.trans_scalD[:,1])
    end

    if num.nb_transported_scalars>1
        minscal2L = minimum(phL.trans_scalD[:,2])
        maxscal2L = maximum(phL.trans_scalD[:,2])
        moyscal2L = mean(phL.trans_scalD[:,2])
    
        if num.nb_transported_scalars>2
            minscal3L = minimum(phL.trans_scalD[:,3])
            maxscal3L = maximum(phL.trans_scalD[:,3])
            moyscal3L = mean(phL.trans_scalD[:,3])
        end

    end

    # if heat
    #     minTL = minimum(phL.T)
    #     maxTL = maximum(phL.T)
    #     moyTL = mean(phL.T)
    # end

    minuL = minimum(phL.uD)
    maxuL = maximum(phL.uD)
    moyuL = mean(phL.uD)

    minvL = minimum(phL.vD)
    maxvL = maximum(phL.vD)
    moyvL = mean(phL.vD)

    minpL = minimum(phL.pD)
    maxpL = maximum(phL.pD)
    moypL = mean(phL.pD)

    minTL = minimum(phL.TD)
    maxTL = maximum(phL.TD)
    moyTL = mean(phL.TD)

    # minphi_eleL = minimum(phL.phi_eleD)
    # miniL=minimum(phL.i_current_mag)

    # maxphi_eleL = maximum(phL.phi_eleD)
    # maxiL=maximum(phL.i_current_mag)

    # moyphi_eleL = mean(phL.phi_eleD)
    # moyiL=mean(phL.i_current_mag)

    # print("$(@sprintf("norm(cH2) %.6e", normscal1L))\t$(@sprintf("norm(KOH) %.6e", normscal2L))\t$(@sprintf("norm(H2O) %.6e", normscal3L))\n")
    # print("$(@sprintf("norm(phi_ele) %.6e", normphi_eleL))\t$(@sprintf("norm(T) %.6e", normTL))\t$(@sprintf("norm(i) %.6e", normiL))\n")

    print("$(@sprintf("min(u) %.6e", minuL))\t$(@sprintf("min(v) %.6e", minvL))\t$(@sprintf("min(p) %.6e", minpL))\n")
    print("$(@sprintf("max(u) %.6e", maxuL))\t$(@sprintf("max(v) %.6e", maxvL))\t$(@sprintf("max(p) %.6e", maxpL))\n")
    print("$(@sprintf("moy(u) %.6e", moyuL))\t$(@sprintf("moy(v) %.6e", moyvL))\t$(@sprintf("moy(p) %.6e", moypL))\n")

    print("$(@sprintf("min(cH2) %.6e", minscal1L))\t$(@sprintf("min(KOH) %.6e", minscal2L))\t$(@sprintf("min(H2O) %.6e", minscal3L))\n")
    print("$(@sprintf("max(cH2) %.6e", maxscal1L))\t$(@sprintf("max(KOH) %.6e", maxscal2L))\t$(@sprintf("max(H2O) %.6e", maxscal3L))\n")
    print("$(@sprintf("moy(cH2) %.6e", moyscal1L))\t$(@sprintf("moy(KOH) %.6e", moyscal2L))\t$(@sprintf("moy(H2O) %.6e", moyscal3L))\n")

    # print("$(@sprintf("min(phi_ele) %.6e", minphi_eleL))\t$(@sprintf("min(T) %.6e", minTL))\t$(@sprintf("min(i) %.6e", miniL))\n")
    # print("$(@sprintf("max(phi_ele) %.6e", maxphi_eleL))\t$(@sprintf("max(T) %.6e", maxTL))\t$(@sprintf("max(i) %.6e", maxiL))\n")
    # print("$(@sprintf("moy(phi_ele) %.6e", moyphi_eleL))\t$(@sprintf("moy(T) %.6e", moyTL))\t$(@sprintf("moy(i) %.6e", moyiL))\n")

    # print("mean cH2 interface",mean(veci(phL.trans_scalD[:,1],grid,2)))
    # print("mean cH2 interface",average!(reshape(veci(phL.trans_scalD[:,iscal],grid,2), grid), grid, LS[1].geoL, num))

    # printstyled(color=:green, @sprintf "\n average c %s\n" average!(reshape(veci(phL.trans_scalD[:,iscal],grid,2), grid), grid, LS[1].geoL, num))
    
    # print("nonzero\n")
    # print(nonzero)
    # print("\nmean c(H2) interface \n",mean(nonzero))
    if num.nb_transported_scalars>0
        # nonzero = veci(phL.trans_scalD[:,1],grid,2)[abs.(veci(phL.trans_scalD[:,1],grid,2)) .> 0.0]
        # printstyled(color=:green, @sprintf "\n mean c(H2) interface : %.2e\n" mean(nonzero))

        # iLS = 1
        # nonzero = mean_intfc_non_null(phL.trans_scalD,1,grid,iLS)
        # printstyled(color=:green, @sprintf "\n mean  c(H2) interface : %.2e\n" nonzero)

        for iscal = 1:num.nb_transported_scalars

        printstyled(color=:green, @sprintf "\n iscal  %.2i\n" iscal)

        for iLS in 1:num.nLS+1
            @views nonzero = mean_intfc_non_null_v3( phL.trans_scalD[:,iscal],grid,iLS)
            printstyled(color=:green, @sprintf "\n index %.2i : mean %.2e max %.2e\n" iLS nonzero maximum(veci(phL.trans_scalD[:,iscal],grid,iLS)))
        end
        # printstyled(color=:green, @sprintf "\n mean  c bulk : %.2e %.2e\n" nonzero maximum(veci(phL.trans_scalD[:,iscal],grid,iLS+1)))

        # iLS = 1
        # @views nonzero = mean_intfc_non_null_v2( phL.trans_scalD[:,iscal],grid,iLS)
        # # printstyled(color=:green, @sprintf "\n mean  c interface : %.2e\n" nonzero)
        # printstyled(color=:green, @sprintf "\n mean  c interface : %.2e %.2e\n" nonzero maximum(veci(phL.trans_scalD[:,iscal],grid,iLS+1)))


        # iLS = 2
        # @views nonzero = mean_intfc_non_null_v2( phL.trans_scalD[:,iscal], grid,iLS)
        # # printstyled(color=:green, @sprintf "\n mean  c interface : %.2e\n" nonzero)
        # printstyled(color=:green, @sprintf "\n mean  c interface : %.2e %.2e\n" nonzero maximum(veci(phL.trans_scalD[:,iscal],grid,iLS+1)))

        printstyled(color=:green, @sprintf "\n max c borders : %.2e\n" maximum(vecb(phL.trans_scalD[:,iscal],grid)))



        end

        for iLS in 1:num.nLS
            nonzero = mean_intfc_non_null_v3(phL.phi_eleD,grid,iLS)
            printstyled(color=:green, @sprintf "\n mean  c iLS %.2i  : %.2e %.2e\n" iLS nonzero maximum(veci(phL.phi_eleD,grid,iLS+1)))
        end

        # iLS = 0
        # nonzero = mean_intfc_non_null_v2(phL.phi_eleD,grid,iLS)
        # printstyled(color=:green, @sprintf "\n mean  phi bulk : %.2e\n" nonzero)

        # iLS = 1
        # nonzero = mean_intfc_non_null_v2(phL.phi_eleD,grid,iLS)
        # # printstyled(color=:green, @sprintf "\n mean  phi interface : %.2e\n" nonzero)
        # printstyled(color=:green, @sprintf "\n mean  c interface : %.2e %.2e\n" nonzero maximum(veci(phL.phi_eleD,grid,iLS+1)))

        # iLS = 2
        # nonzero = mean_intfc_non_null_v2(phL.phi_eleD,grid,iLS)
        # # printstyled(color=:green, @sprintf "\n mean  phi interface : %.2e\n" nonzero)
        # printstyled(color=:green, @sprintf "\n mean  c interface : %.2e %.2e\n" nonzero maximum(veci(phL.phi_eleD,grid,iLS+1)))
        
        printstyled(color=:green, @sprintf "\n max c borders : %.2e\n" maximum(vecb(phL.phi_eleD,grid)))

    end

end

function print_electrolysis_statistics_bulk(nb_transported_scalars::Int64,
    grid::Mesh{Flower.GridCC, Float64, Int64}, phL::Phase{Float64})
    # normscal1L = norm(phL.trans_scal[:,:,1])
    # normscal2L = norm(phL.trans_scal[:,:,2])
    # normscal3L = norm(phL.trans_scal[:,:,3])
    # normphi_eleL = norm(phL.phi_ele)
    # normTL = norm(phL.T)
    # normiL=0.0 #TODO

    

    minscal1L=minscal2L=minscal3L=minphi_eleL=minTL=miniL=0.0
    maxscal1L=maxscal2L=maxscal3L=maxphi_eleL=maxTL=maxiL=0.0
    moyscal1L=moyscal2L=moyscal3L=moyphi_eleL=moyTL=moyiL=0.0


    if nb_transported_scalars>0
        minscal1L = minimum(phL.trans_scal[:,:,1])
        maxscal1L = maximum(phL.trans_scal[:,:,1])
        moyscal1L = mean(phL.trans_scal[:,:,1])
    end

    if nb_transported_scalars>1
        minscal2L = minimum(phL.trans_scal[:,:,2])
        maxscal2L = maximum(phL.trans_scal[:,:,2])
        moyscal2L = mean(phL.trans_scal[:,:,2])
    
        if nb_transported_scalars>2
            minscal3L = minimum(phL.trans_scal[:,:,3])
            maxscal3L = maximum(phL.trans_scal[:,:,3])
            moyscal3L = mean(phL.trans_scal[:,:,3])
        end

    end

    # if heat
    #     minTL = minimum(phL.T)
    #     maxTL = maximum(phL.T)
    #     moyTL = mean(phL.T)
    # end

    minuL = minimum(phL.u)
    maxuL = maximum(phL.u)
    moyuL = mean(phL.u)

    minvL = minimum(phL.v)
    maxvL = maximum(phL.v)
    moyvL = mean(phL.v)

    minpL = minimum(phL.p)
    maxpL = maximum(phL.p)
    moypL = mean(phL.p)

    minTL = minimum(phL.T)
    maxTL = maximum(phL.T)
    moyTL = mean(phL.T)

    minphi_eleL = minimum(phL.phi_ele)
    miniL=minimum(phL.i_current_mag)

    maxphi_eleL = maximum(phL.phi_ele)
    maxiL=maximum(phL.i_current_mag)

    moyphi_eleL = mean(phL.phi_ele)
    moyiL=mean(phL.i_current_mag)

    # print("$(@sprintf("norm(cH2) %.6e", normscal1L))\t$(@sprintf("norm(KOH) %.6e", normscal2L))\t$(@sprintf("norm(H2O) %.6e", normscal3L))\n")
    # print("$(@sprintf("norm(phi_ele) %.6e", normphi_eleL))\t$(@sprintf("norm(T) %.6e", normTL))\t$(@sprintf("norm(i) %.6e", normiL))\n")

    print("$(@sprintf("min(u) %.6e", minuL))\t$(@sprintf("min(v) %.6e", minvL))\t$(@sprintf("min(p) %.6e", minpL))\n")
    print("$(@sprintf("max(u) %.6e", maxuL))\t$(@sprintf("max(v) %.6e", maxvL))\t$(@sprintf("max(p) %.6e", maxpL))\n")
    print("$(@sprintf("moy(u) %.6e", moyuL))\t$(@sprintf("moy(v) %.6e", moyvL))\t$(@sprintf("moy(p) %.6e", moypL))\n")

    print("$(@sprintf("min(cH2) %.6e", minscal1L))\t$(@sprintf("min(KOH) %.6e", minscal2L))\t$(@sprintf("min(H2O) %.6e", minscal3L))\n")
    print("$(@sprintf("max(cH2) %.6e", maxscal1L))\t$(@sprintf("max(KOH) %.6e", maxscal2L))\t$(@sprintf("max(H2O) %.6e", maxscal3L))\n")
    print("$(@sprintf("moy(cH2) %.6e", moyscal1L))\t$(@sprintf("moy(KOH) %.6e", moyscal2L))\t$(@sprintf("moy(H2O) %.6e", moyscal3L))\n")

    print("$(@sprintf("min(phi_ele) %.6e", minphi_eleL))\t$(@sprintf("min(T) %.6e", minTL))\t$(@sprintf("min(i) %.6e", miniL))\n")
    print("$(@sprintf("max(phi_ele) %.6e", maxphi_eleL))\t$(@sprintf("max(T) %.6e", maxTL))\t$(@sprintf("max(i) %.6e", maxiL))\n")
    print("$(@sprintf("moy(phi_ele) %.6e", moyphi_eleL))\t$(@sprintf("moy(T) %.6e", moyTL))\t$(@sprintf("moy(i) %.6e", moyiL))\n")
    # print("mean cH2 interface",mean(veci(phL.trans_scalD[:,1],grid,2)))
    # print("mean cH2 interface",average!(reshape(veci(phL.trans_scalD[:,iscal],grid,2), grid), grid, LS[1].geoL, num))

    # printstyled(color=:green, @sprintf "\n average c %s\n" average!(reshape(veci(phL.trans_scalD[:,iscal],grid,2), grid), grid, LS[1].geoL, num))
    
    # print("nonzero\n")
    # print(nonzero)
    # print("\nmean c(H2) interface \n",mean(nonzero))
    if nb_transported_scalars>0
        # nonzero = veci(phL.trans_scalD[:,1],grid,2)[abs.(veci(phL.trans_scalD[:,1],grid,2)) .> 0.0]
        # printstyled(color=:green, @sprintf "\n mean c(H2) interface : %.2e\n" mean(nonzero))

        iLS = 1
        nonzero = mean_intfc_non_null(phL.trans_scalD,1,grid,iLS)
        printstyled(color=:green, @sprintf "\n mean  c(H2) interface : %.2e\n" nonzero)
    end

end


# Gradient of pressure, eq. 17 in 
#"A Conservative Cartesian Cut-Cell Method for Mixed Boundary Conditions and the Incompressible Navier-Stokes Equations on Staggered Meshes"
#From navier_stokes_coupled.jl
# ∇ϕ_x = opC_u.AxT * opC_u.Rx * vec(phi_ele) .+ opC_u.Gx_b * vecb(phi_eleD,grid)
# ∇ϕ_y = opC_v.AyT * opC_v.Ry * vec(phi_ele) .+ opC_v.Gy_b * vecb(phi_eleD,grid)
# for iLS in 1:nLS
#     ∇ϕ_x .+= opC_u.Gx[iLS] * veci(phi_eleD,grid,iLS+1)
#     ∇ϕ_y .+= opC_v.Gy[iLS] * veci(phi_eleD,grid,iLS+1)
# end

# function fill_reshape(a,b)

#     # print("\n size ",size(a)," ",size(a,1)," ",size(a)[1])
#     for j in 1:size(a,1)
#         for i in 1:size(a,2)
#             a[j,i] = b[(j-1)*size(a,1)+i]
#         end
#     end
# end

# function fill_reshapenew(a,b,nx,ny)
#     for j in 1:ny
#         for i in 1:nx
#             a[j,i] = b[(j-1)*ny+i]
#         end
#     end
# end


# function fill_reshape2(a,b)

#     # print("\n size ",size(a)," ",size(a,1)," ",size(a)[1])
#     size1,size2 = size(a)
#     # size2 = size(a,2)
#     for j in 1:size1
#         for i in 1:size2
#             a[j,i] = b[(j-1)*size1+i]
#         end
#     end
# end

# function fill_reshape2(a,b)

#     # print("\n size ",size(a)," ",size(a,1)," ",size(a)[1])
#     size1,size2 = size(a)
#     # size2 = size(a,2)
#     for j in 1:size1 #with threads
#         @views a[j,:] = b[(j-1)*size1+1:j*size1]
#         # for i in 1:size2
#         #     a[j,i] = b[(j-1)*size1+i]
#         # end
#     end
# end

"""
  Compute norm of gradient for exchange current 
    phL.i_current_mag is interpolated
    gradient is cell-averaged
"""
function compute_grad_phi_ele!(num::Numerical{Float64, Int64},
    grid::Mesh{Flower.GridCC, Float64, Int64},
    grid_u::Mesh{Flower.GridFCx, Float64, Int64},
    grid_v::Mesh{Flower.GridFCy, Float64, Int64},
    LS_u::Levelset{Float64, Int64},
    LS_v::Levelset{Float64, Int64},
    phL::Phase{Float64},
    phS::Phase{Float64}, 
    opC_pL::Operators{Float64, Int64}, 
    opC_pS::Operators{Float64, Int64},
    elec_cond::Array{Float64, 2},
    tmp_vec_u::Array{Float64, 2},
    tmp_vec_v::Array{Float64, 2},
    tmp_vec_p::Array{Float64, 2},
    tmp_vec_p0::Array{Float64, 2},
    tmp_vec_p1::Array{Float64, 2},
    )
    
    @unpack nLS = num

    # LS =gp.LS[1]

    #TODO different way to do it?

    #Liquid phase
    @unpack phi_eleD = phL
    opC_p = opC_pL


    # tmp_vec_p1 .= 0.0

    tmp_vec_p1 = zeros(grid.ny,grid.nx)

    for j in 1:grid.ny
        for i in 1:grid.nx
            tmp_vec_p1[j,i] = 0.0
        end
    end

    if maximum(tmp_vec_p1)>0
        printstyled(color=:red, @sprintf "\n tmp_vec_p1 max %.2e \n" maximum(tmp_vec_p1))
        throw(DivideError())
    end

    # mass_flux_vec1   = opC_p.HxT[iLStmp] * opC_p.iMx * opC_p.Bx * vec1(scalD,grid) .+ opC_p.HyT[iLStmp] * opC_p.iMy * opC_p.By * vec1(scalD,grid)
    # mass_flux_vecb   = opC_p.HxT[iLStmp] * opC_p.iMx_b * opC_p.Hx_b * vecb(scalD,grid) .+ opC_p.HyT[iLStmp] *  opC_p.iMy_b * opC_p.Hy_b * vecb(scalD,grid)

    # for iLS in 1:nLS
    #     mass_flux_veci .+= opC_p.HxT[iLS] * opC_p.iMx * opC_p.Hx[iLS] * veci(scalD,grid,iLS+1)
    #     mass_flux_veci .+= opC_p.HyT[iLS] * opC_p.iMy * opC_p.Hy[iLS] * veci(scalD,grid,iLS+1)
    # end

 
    # mass_flux_vec1_2 .= reshape(mass_flux_vec1,grid)
    # mass_flux_vecb_2 .= reshape(mass_flux_vecb,grid)
    # mass_flux_veci_2 .= reshape(mass_flux_veci,grid)

    # mass_flux .= mass_flux_vec1_2 .+ mass_flux_vecb_2 .+ mass_flux_veci_2


    ∇ϕ_x = opC_p.iMx * opC_p.Bx * vec1(phi_eleD,grid) .+ opC_p.iMx_b * opC_p.Hx_b * vecb(phi_eleD,grid)
    ∇ϕ_y = opC_p.iMy * opC_p.By * vec1(phi_eleD,grid) .+ opC_p.iMy_b * opC_p.Hy_b * vecb(phi_eleD,grid)

    for iLS in 1:nLS
        ∇ϕ_x .+= opC_p.iMx * opC_p.Hx[iLS] * veci(phi_eleD,grid,iLS+1)
        ∇ϕ_y .+= opC_p.iMy * opC_p.Hy[iLS] * veci(phi_eleD,grid,iLS+1)
    end

    tmp_vec_u .= reshape(veci(∇ϕ_x,grid_u,1), grid_u)
    tmp_vec_v .= reshape(veci(∇ϕ_y,grid_v,1), grid_v)

   

    # tmp_vec_u .= grd_x
    # tmp_vec_v .= grd_y

    #store in us, vs instead of Eus, Evs
    # interpolate_grid_liquid!(grid,grid_u,grid_v,tmp_vec_u, tmp_vec_v ,tmp_vec_p,tmp_vec_p0)

    interpolate_grid_liquid_2!(num, grid, LS_u, LS_v, tmp_vec_u, tmp_vec_v, tmp_vec_p, tmp_vec_p0)

    # to expose the gradient of phi
    # @ccall "libpdi".PDI_multi_expose("write_data_elec_ix_iy_grad"::Cstring,
    # "grad_phi_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
    # "grad_phi_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,  
    # # "i_current_mag"::Cstring, phL.i_current_mag::Ptr{Cdouble}, PDI_OUT::Cint,
    # # "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
    # C_NULL::Ptr{Cvoid})::Cint

    # TODO 
    tmp_vec_p .*= -elec_cond # i=-κ∇ϕ here magnitude
    tmp_vec_p0 .*= -elec_cond # i=-κ∇ϕ here magnitude




    @ccall "libpdi".PDI_multi_expose("write_data_elec_ix_iy"::Cstring,
    "i_current_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
    "i_current_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,  
    # "i_current_mag"::Cstring, phL.i_current_mag::Ptr{Cdouble}, PDI_OUT::Cint,
    "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
    C_NULL::Ptr{Cvoid})::Cint


    # @ccall "libpdi".PDI_multi_expose("solve_poisson"::Cstring,
    # "i_current_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
    # "i_current_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,  
    # # "i_current_mag"::Cstring, phL.i_current_mag::Ptr{Cdouble}, PDI_OUT::Cint,
    # "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
    # "elec_cond_1D"::Cstring, elec_condD::Ptr{Cdouble}, PDI_OUT::Cint,  
    # # "BC_phi_ele_left"::Cstring, BC_phi_ele.left.val::Ptr{Cdouble}, PDI_OUT::Cint,  
    # # "grad_phi_ele_u"::Cstring, tmp_vec_u::Ptr{Cdouble}, PDI_OUT::Cint,  
    # C_NULL::Ptr{Cvoid})::Cint

    # tmp_vec_p1 .= (
    #     (tmp_vec_u[:,2:end].^2.0 .* LS_u.geoL.dcap[:,2:end,6] .+ 
    #     tmp_vec_u[:,1:end-1].^2.0 .* LS_u.geoL.dcap[:,1:end-1,6]) ./ 
    #     (LS_u.geoL.dcap[:,1:end-1,6] .+ LS_u.geoL.dcap[:,2:end,6])
    # )
    # tmp_vec_p1 .+= (
    #     (tmp_vec_v[2:end,:].^2.0 .* LS_v.geoL.dcap[2:end,:,7] .+ 
    #     tmp_vec_v[1:end-1,:].^2.0 .* LS_v.geoL.dcap[1:end-1,:,7]) ./
    #     (LS_v.geoL.dcap[1:end-1,:,7] .+ LS_v.geoL.dcap[2:end,:,7])
    # )
    if num.average_liquid_solid == 0 #solid contribution set to zeoro

        # printstyled(color=:red, @sprintf "\n tmp_vec_p1 25 49 %.2e \n" tmp_vec_p1[49,25])

    

        for j in 1:grid.ny
            for i in 1:grid.nx
                # print("\n i j tmp_vec_p1[j,i] ",i, " ",j," ",tmp_vec_p1[j,i])
                # printstyled(color=:red, @sprintf "\n tmp_vec_p1 25 49 %.2e \n" tmp_vec_p1[49,25])

                if tmp_vec_p1[j,i] !=0.0
                    printstyled(color=:red, @sprintf "\n  i %.3i j  %.3i tmp_vec_p1 %.2e \n" i j tmp_vec_p1[j,i])
                    return
                end

                if LS_u.geoL.dcap[j,i,5] >num.epsilon_vol
                    tmp_vec_p1[j,i] += (tmp_vec_u[j,i+1]^2.0 * LS_u.geoL.dcap[j,i+1,6] + tmp_vec_u[j,i]^2.0 * LS_u.geoL.dcap[j,i,6]) / (LS_u.geoL.dcap[j,i+1,6] + LS_u.geoL.dcap[j,i,6] + LS_u.geoS.dcap[j,i+1,6] + LS_u.geoS.dcap[j,i,6])

                    tmp_vec_p1[j,i] += (tmp_vec_v[j+1,i]^2.0 * LS_v.geoL.dcap[j+1,i,7] + tmp_vec_v[j,i]^2.0 * LS_v.geoL.dcap[j,i,7]) / (LS_v.geoL.dcap[j+1,i,7] + LS_v.geoL.dcap[j,i,7] + LS_v.geoS.dcap[j+1,i,7] + LS_v.geoS.dcap[j,i,7])

                    try
                        tmp_vec_p1[j,i] = sqrt(tmp_vec_p1[j,i])
                    catch e
                        print("\n i j liq ",i," ",j," ", tmp_vec_p1[j,i]," ",LS_u.geoL.dcap[j,i+1,6] , " ", LS_u.geoL.dcap[j,i,6] ," ", LS_u.geoS.dcap[j,i+1,6] ," ", LS_u.geoS.dcap[j,i,6]," ",LS_u.geoL.dcap[j,i+1,6] + LS_u.geoL.dcap[j,i,6] + LS_u.geoS.dcap[j,i+1,6] + LS_u.geoS.dcap[j,i,6] ," ",tmp_vec_u[j,i+1]^2.0 ," ", tmp_vec_u[j,i]^2.0," ",tmp_vec_v[j+1,i]^2.0 ," ", tmp_vec_v[j,i]^2.0)
                        print("\n i j liq ",i," ",j," ", tmp_vec_p1[j,i]," ",LS_v.geoL.dcap[j+1,i,7] , " ", LS_v.geoL.dcap[j,i,7] ," ", LS_v.geoL.dcap[j+1,i,7] , " ", LS_v.geoL.dcap[j,i,7])

                        print(e)
                        return
                    end

                    # if (LS_u.geoL.dcap[j,i+1,6] + LS_u.geoL.dcap[j,i,6] 
                        
                    #     + LS_u.geoS.dcap[j,i+1,6] + )

                end

            end
        end 

    end # if num.average_liquid_solid == 0


    # if num.average_liquid_solid == 1 #TODO

    #     #Solid phase
    #     @unpack phi_eleD = phS

    #     opC_p = opC_pS

    #     ∇ϕ_x = opC_p.iMx * opC_p.Bx * vec1(phi_eleD,grid) .+ opC_p.iMx_b * opC_p.Hx_b * vecb(phi_eleD,grid)
    #     ∇ϕ_y = opC_p.iMy * opC_p.By * vec1(phi_eleD,grid) .+ opC_p.iMy_b * opC_p.Hy_b * vecb(phi_eleD,grid)

    #     for iLS in 1:nLS
    #         ∇ϕ_x .+= opC_p.iMx * opC_p.Hx[iLS] * veci(phi_eleD,grid,iLS+1)
    #         ∇ϕ_y .+= opC_p.iMy * opC_p.Hy[iLS] * veci(phi_eleD,grid,iLS+1)
    #     end

    #     # grd_x .= reshape(veci(∇ϕ_x,grid_u,1), grid_u)
    #     # grd_y .= reshape(veci(∇ϕ_y,grid_v,1), grid_v)

    #     tmp_vec_u .= reshape(veci(∇ϕ_x,grid_u,1), grid_u)
    #     tmp_vec_v .= reshape(veci(∇ϕ_y,grid_v,1), grid_v)

    #     #Output everything to phL

    #     # phS.i_current_mag .= (
    #     #     (grd_x[:,2:end].^2.0 .* LS_u.geoS.dcap[:,2:end,6] .+ 
    #     #     grd_x[:,1:end-1].^2.0 .* LS_u.geoS.dcap[:,1:end-1,6]) ./ 
    #     #     (LS_u.geoS.dcap[:,1:end-1,6] .+ LS_u.geoS.dcap[:,2:end,6] + eps_den )
    #     # )
    #     # phS.i_current_mag .+= (
    #     #     (grd_y[2:end,:].^2.0 .* LS_v.geoS.dcap[2:end,:,7] .+ 
    #     #     ph.v[1:end-1,:].^2.0 .* LS_v.geoS.dcap[1:end-1,:,7]) ./
    #     #     (LS_v.geoS.dcap[1:end-1,:,7] .+ LS_v.geoS.dcap[2:end,:,7] + eps_den )
    #     # )

    #     # tmp_vec_p1 .+= (
    #     #     (tmp_vec_u[:,2:end].^2.0 .* LS_u.geoS.dcap[:,2:end,6] .+ 
    #     #     tmp_vec_u[:,1:end-1].^2.0 .* LS_u.geoS.dcap[:,1:end-1,6]) ./ 
    #     #     (LS_u.geoS.dcap[:,1:end-1,6] .+ LS_u.geoS.dcap[:,2:end,6])
    #     # )
    #     # # #Store also S value or reset 0
    #     # tmp_vec_p1 .+= (
    #     #     (tmp_vec_v[2:end,:].^2.0 .* LS_v.geoS.dcap[2:end,:,7] .+ 
    #     #     tmp_vec_v[1:end-1,:].^2.0 .* LS_v.geoS.dcap[1:end-1,:,7]) ./
    #     #     (LS_v.geoS.dcap[1:end-1,:,7] .+ LS_v.geoS.dcap[2:end,:,7] )
    #     # )
        
    #     # for j in 1:grid.ny
    #     #     for i in 1:grid.nx

    #     #         if LS_u.geoS.dcap[j,i,5] >num.epsilon_vol
    #     #             print("\n test solid ",tmp_vec_p1[j,i]," ", (tmp_vec_u[j,i+1]^2.0 * LS_u.geoS.dcap[j,i+1,6] + tmp_vec_u[j,i]^2.0 * LS_u.geoS.dcap[j,i,6]) / (LS_u.geoS.dcap[j,i+1,6] + LS_u.geoS.dcap[j,i,6])," ", (tmp_vec_v[j+1,i]^2.0 * LS_v.geoS.dcap[j+1,i,7] + tmp_vec_v[j,i]^2.0 * LS_v.geoS.dcap[j,i,7]) / (LS_v.geoS.dcap[j+1,i,7] + LS_v.geoS.dcap[j,i,7]))
    #     #             tmp_vec_p1[j,i] += (tmp_vec_u[j,i+1]^2.0 * LS_u.geoS.dcap[j,i+1,6] + tmp_vec_u[j,i]^2.0 * LS_u.geoS.dcap[j,i,6]) / (LS_u.geoS.dcap[j,i+1,6] + LS_u.geoS.dcap[j,i,6])

    #     #             tmp_vec_p1[j,i] += (tmp_vec_v[j+1,i]^2.0 * LS_v.geoS.dcap[j+1,i,7] + tmp_vec_v[j,i]^2.0 * LS_v.geoS.dcap[j,i,7]) / (LS_v.geoS.dcap[j+1,i,7] + LS_v.geoS.dcap[j,i,7])
                    
    #     #             try
    #     #                 tmp_vec_p1[j,i] = sqrt(tmp_vec_p1[j,i])
    #     #             catch e
    #     #                 print("\n i j sol ",i," ",j," ", tmp_vec_p1[j,i])
    #     #                 print(e)
    #     #             end

    #     #         end

    #     #     end
    #     # end 


    # end

    # tmp_vec_p1 .= sqrt.(tmp_vec_p1)
    # phS.i_current_mag .= sqrt.(phS.i_current_mag)

    tmp_vec_p1 .*= elec_cond # i=-κ∇ϕ here magnitude


    #TODO Eu Ev S
    # phS.Eu .= grd_x
    # phS.Ev .= grd_y





    @ccall "libpdi".PDI_multi_expose("write_data_elec_imag"::Cstring,
    # "i_current_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
    # "i_current_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,  
    "i_current_mag"::Cstring, tmp_vec_p1::Ptr{Cdouble}, PDI_OUT::Cint,
    # "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
    C_NULL::Ptr{Cvoid})::Cint


    # minphi_eleL = minimum(phL.phi_eleD)
    # miniL=minimum(phL.i_current_mag)

    # maxphi_eleL = maximum(phL.phi_eleD)
    # maxiL=maximum(phL.i_current_mag)

    # moyphi_eleL = mean(phL.phi_eleD)
    # moyiL=mean(phL.i_current_mag)

    # print("$(@sprintf("norm(cH2) %.6e", normscal1L))\t$(@sprintf("norm(KOH) %.6e", normscal2L))\t$(@sprintf("norm(H2O) %.6e", normscal3L))\n")
    # print("$(@sprintf("norm(phi_ele) %.6e", normphi_eleL))\t$(@sprintf("norm(T) %.6e", normTL))\t$(@sprintf("norm(i) %.6e", normiL))\n")

    # print("$(@sprintf("min(phi_ele) %.6e", minphi_eleL))\t$(@sprintf("min(T) %.6e", minTL))\t$(@sprintf("min(i) %.6e", miniL))\n")
    # print("$(@sprintf("max(phi_ele) %.6e", maxphi_eleL))\t$(@sprintf("max(T) %.6e", maxTL))\t$(@sprintf("max(i) %.6e", maxiL))\n")
    # print("$(@sprintf("moy(phi_ele) %.6e", moyphi_eleL))\t$(@sprintf("moy(T) %.6e", moyTL))\t$(@sprintf("moy(i) %.6e", moyiL))\n")


end



"""
  Compute norm of gradient for exchange current 
"""
function compute_grad_phi_ele_liquid_solid!(num::Numerical{Float64, Int64},
    grid::Mesh{Flower.GridCC, Float64, Int64},
    grid_u::Mesh{Flower.GridFCx, Float64, Int64},
    grid_v::Mesh{Flower.GridFCy, Float64, Int64},
    phL::Phase{Float64},
    phS::Phase{Float64}, 
    opC_pL::Operators{Float64, Int64}, 
    opC_pS::Operators{Float64, Int64})
    
    @unpack nLS = num

    LS_u =grid_u.LS[1]
    LS_v = grid_v.LS[1]
    # LS =gp.LS[1]

    eps_den = 1e-8
    #TODO different way to do it?

    #Liquid phase
    @unpack phi_eleD = phL
    opC_p = opC_pL


    # mass_flux_vec1   = opC_p.HxT[iLStmp] * opC_p.iMx * opC_p.Bx * vec1(scalD,grid) .+ opC_p.HyT[iLStmp] * opC_p.iMy * opC_p.By * vec1(scalD,grid)
    # mass_flux_vecb   = opC_p.HxT[iLStmp] * opC_p.iMx_b * opC_p.Hx_b * vecb(scalD,grid) .+ opC_p.HyT[iLStmp] *  opC_p.iMy_b * opC_p.Hy_b * vecb(scalD,grid)

    # for iLS in 1:nLS
    #     mass_flux_veci .+= opC_p.HxT[iLS] * opC_p.iMx * opC_p.Hx[iLS] * veci(scalD,grid,iLS+1)
    #     mass_flux_veci .+= opC_p.HyT[iLS] * opC_p.iMy * opC_p.Hy[iLS] * veci(scalD,grid,iLS+1)
    # end

 
    # mass_flux_vec1_2 .= reshape(mass_flux_vec1,grid)
    # mass_flux_vecb_2 .= reshape(mass_flux_vecb,grid)
    # mass_flux_veci_2 .= reshape(mass_flux_veci,grid)

    # mass_flux .= mass_flux_vec1_2 .+ mass_flux_vecb_2 .+ mass_flux_veci_2


    ∇ϕ_x = opC_p.iMx * opC_p.Bx * vec1(phi_eleD,grid) .+ opC_p.iMx_b * opC_p.Hx_b * vecb(phi_eleD,grid)
    ∇ϕ_y = opC_p.iMy * opC_p.By * vec1(phi_eleD,grid) .+ opC_p.iMy_b * opC_p.Hy_b * vecb(phi_eleD,grid)

    for iLS in 1:nLS
        ∇ϕ_x .+= opC_p.iMx * opC_p.Hx[iLS] * veci(phi_eleD,grid,iLS+1)
        ∇ϕ_y .+= opC_p.iMy * opC_p.Hy[iLS] * veci(phi_eleD,grid,iLS+1)
    end

    grd_x .= reshape(veci(∇ϕ_x,grid_u,1), grid_u)
    grd_y .= reshape(veci(∇ϕ_y,grid_v,1), grid_v)

    phL.i_current_mag .= (
        (grd_x[:,2:end].^2.0 .* LS_u.geoL.dcap[:,2:end,6] .+ 
        grd_x[:,1:end-1].^2.0 .* LS_u.geoL.dcap[:,1:end-1,6]) ./ 
        (LS_u.geoL.dcap[:,1:end-1,6] .+ LS_u.geoL.dcap[:,2:end,6] .+ eps_den )
    )
    phL.i_current_mag .+= (
        (grd_y[2:end,:].^2.0 .* LS_v.geoL.dcap[2:end,:,7] .+ 
        grd_y[1:end-1,:].^2.0 .* LS_v.geoL.dcap[1:end-1,:,7]) ./
        (LS_v.geoL.dcap[1:end-1,:,7] .+ LS_v.geoL.dcap[2:end,:,7] .+ eps_den )
    )

    phL.Eu .= grd_x
    phL.Ev .= grd_y

    #Solid phase
    @unpack phi_eleD = phS

    opC_p = opC_pS

    ∇ϕ_x = opC_p.iMx * opC_p.Bx * vec1(phi_eleD,grid) .+ opC_p.iMx_b * opC_p.Hx_b * vecb(phi_eleD,grid)
    ∇ϕ_y = opC_p.iMy * opC_p.By * vec1(phi_eleD,grid) .+ opC_p.iMy_b * opC_p.Hy_b * vecb(phi_eleD,grid)

    for iLS in 1:nLS
        ∇ϕ_x .+= opC_p.iMx * opC_p.Hx[iLS] * veci(phi_eleD,grid,iLS+1)
        ∇ϕ_y .+= opC_p.iMy * opC_p.Hy[iLS] * veci(phi_eleD,grid,iLS+1)
    end

    grd_x .= reshape(veci(∇ϕ_x,grid_u,1), grid_u)
    grd_y .= reshape(veci(∇ϕ_y,grid_v,1), grid_v)

    #Output everything to phL

    # phS.i_current_mag .= (
    #     (grd_x[:,2:end].^2.0 .* LS_u.geoS.dcap[:,2:end,6] .+ 
    #     grd_x[:,1:end-1].^2.0 .* LS_u.geoS.dcap[:,1:end-1,6]) ./ 
    #     (LS_u.geoS.dcap[:,1:end-1,6] .+ LS_u.geoS.dcap[:,2:end,6] + eps_den )
    # )
    # phS.i_current_mag .+= (
    #     (grd_y[2:end,:].^2.0 .* LS_v.geoS.dcap[2:end,:,7] .+ 
    #     ph.v[1:end-1,:].^2.0 .* LS_v.geoS.dcap[1:end-1,:,7]) ./
    #     (LS_v.geoS.dcap[1:end-1,:,7] .+ LS_v.geoS.dcap[2:end,:,7] + eps_den )
    # )

    phL.i_current_mag .+= (
        (grd_x[:,2:end].^2.0 .* LS_u.geoS.dcap[:,2:end,6] .+ 
        grd_x[:,1:end-1].^2.0 .* LS_u.geoS.dcap[:,1:end-1,6]) ./ 
        (LS_u.geoS.dcap[:,1:end-1,6] .+ LS_u.geoS.dcap[:,2:end,6] .+ eps_den )
    )
    # #Store also S value or reset 0
    # phL.i_current_mag .+= (
    #     (grd_y[2:end,:].^2.0 .* LS_v.geoS.dcap[2:end,:,7] .+ 
    #     grd_y[1:end-1,:].^2.0 .* LS_v.geoS.dcap[1:end-1,:,7]) ./
    #     (LS_v.geoS.dcap[1:end-1,:,7] .+ LS_v.geoS.dcap[2:end,:,7] .+ eps_den )
    # )

    phL.i_current_mag .= sqrt.(phL.i_current_mag)
    # phS.i_current_mag .= sqrt.(phS.i_current_mag)

    phS.Eu .= grd_x
    phS.Ev .= grd_y


end

"""
  Compute norm of gradient of variable pD (located on p-grid)

# Gradient of pressure, eq. 17 in 
#"A Conservative Cartesian Cut-Cell Method for Mixed Boundary Conditions and the Incompressible Navier-Stokes Equations on Staggered Meshes"
#From navier_stokes_coupled.jl
```julia
∇ϕ_x = opC_u.AxT * opC_u.Rx * vec(phi_ele) .+ opC_u.Gx_b * vecb(phi_eleD,grid)
∇ϕ_y = opC_v.AyT * opC_v.Ry * vec(phi_ele) .+ opC_v.Gy_b * vecb(phi_eleD,grid)
for iLS in 1:nLS
    ∇ϕ_x .+= opC_u.Gx[iLS] * veci(phi_eleD,grid,iLS+1)
    ∇ϕ_y .+= opC_v.Gy[iLS] * veci(phi_eleD,grid,iLS+1)
end
```
"""
function compute_grad_p_2!(num,grid, grid_u, grid_v, pD, opC_p,opC_u,opC_v,
    # grd_x,grd_y,
    )
    
    grd_x = zeros(grid)
    grd_y = zeros(grid)
    ∇ϕ_x = fnzeros(grid,num)
    ∇ϕ_y = fnzeros(grid,num)


    #Liquid phase
    # @unpack pD = phL

    ∇ϕ_x .= opC_p.iMx * opC_p.Bx * vec1(pD,grid) .+ opC_p.iMx_b * opC_p.Hx_b * vecb(pD,grid)
    ∇ϕ_y .= opC_p.iMy * opC_p.By * vec1(pD,grid) .+ opC_p.iMy_b * opC_p.Hy_b * vecb(pD,grid)

    for iLS in 1:num.nLS
        ∇ϕ_x .+= opC_p.iMx * opC_p.Hx[iLS] * veci(pD,grid,iLS+1)
        ∇ϕ_y .+= opC_p.iMy * opC_p.Hy[iLS] * veci(pD,grid,iLS+1)
    end

    grd_x .= reshape(veci(∇ϕ_x,grid_u,1), grid_u)
    grd_y .= reshape(veci(∇ϕ_y,grid_v,1), grid_v)

    printstyled(color=:red, @sprintf "\n grad min max x %.2e %.2e y %.2e %.2e\n" minimum(grd_x) maximum(grd_x) minimum(grd_y) maximum(grd_y))

    grd_x .= 0.0
    grd_y .= 0.0
    ∇ϕ_x .= 0.0
    ∇ϕ_y .= 0.0

    ∇ϕ_x .= opC_u.AxT * opC_u.Rx * vec1(pD,grid) .+ opC_u.Gx_b * vecb(pD,grid)
    ∇ϕ_y .= opC_v.AyT * opC_v.Ry * vec1(pD,grid) .+ opC_v.Gy_b * vecb(pD,grid)
    for iLS in 1:num.nLS
        ∇ϕ_x .+= opC_u.Gx[iLS] * veci(pD,grid,iLS+1)
        ∇ϕ_y .+= opC_v.Gy[iLS] * veci(pD,grid,iLS+1)
    end

    iMu = Diagonal(inv_weight_eps2.(num.epsilon_mode,num.epsilon_vol,opC_u.M.diag))
    iMv = Diagonal(inv_weight_eps2.(num.epsilon_mode,num.epsilon_vol,opC_v.M.diag))
    ∇ϕ_x .= iMu * ∇ϕ_x
    ∇ϕ_y .= iMv * ∇ϕ_y
    
    grd_x .= reshape(veci(∇ϕ_x,grid_u,1), grid_u)
    grd_y .= reshape(veci(∇ϕ_y,grid_v,1), grid_v)

    printstyled(color=:magenta, @sprintf "\n grad min max x %.2e %.2e y %.2e %.2e\n" minimum(grd_x) maximum(grd_x) minimum(grd_y) maximum(grd_y))
end

"""
  Compute norm of gradient of variable pD (located on p-grid)

# Gradient of pressure, eq. 17 in 
#"A Conservative Cartesian Cut-Cell Method for Mixed Boundary Conditions and the Incompressible Navier-Stokes Equations on Staggered Meshes"
#From navier_stokes_coupled.jl
```julia
∇ϕ_x = opC_u.AxT * opC_u.Rx * vec(phi_ele) .+ opC_u.Gx_b * vecb(phi_eleD,grid)
∇ϕ_y = opC_v.AyT * opC_v.Ry * vec(phi_ele) .+ opC_v.Gy_b * vecb(phi_eleD,grid)
for iLS in 1:nLS
    ∇ϕ_x .+= opC_u.Gx[iLS] * veci(phi_eleD,grid,iLS+1)
    ∇ϕ_y .+= opC_v.Gy[iLS] * veci(phi_eleD,grid,iLS+1)
end
```
"""
function compute_grad_p!(num,grid, grid_u, grid_v, pD, opC_p,opC_u,opC_v)
    

    #Liquid phase
    # @unpack pD = phL

    ∇ϕ_x = opC_p.iMx * opC_p.Bx * vec1(pD,grid) .+ opC_p.iMx_b * opC_p.Hx_b * vecb(pD,grid)
    ∇ϕ_y = opC_p.iMy * opC_p.By * vec1(pD,grid) .+ opC_p.iMy_b * opC_p.Hy_b * vecb(pD,grid)

    for iLS in 1:num.nLS
        ∇ϕ_x .+= opC_p.iMx * opC_p.Hx[iLS] * veci(pD,grid,iLS+1)
        ∇ϕ_y .+= opC_p.iMy * opC_p.Hy[iLS] * veci(pD,grid,iLS+1)
    end

    grd_x .= reshape(veci(∇ϕ_x,grid_u,1), grid_u)
    grd_y .= reshape(veci(∇ϕ_y,grid_v,1), grid_v)

    printstyled(color=:red, @sprintf "\n grad min max x %.2e %.2e y %.2e %.2e\n" minimum(grd_x) maximum(grd_x) minimum(grd_y) maximum(grd_y))


    ∇ϕ_x = opC_u.AxT * opC_u.Rx * vec1(pD,grid) .+ opC_u.Gx_b * vecb(pD,grid)
    ∇ϕ_y = opC_v.AyT * opC_v.Ry * vec1(pD,grid) .+ opC_v.Gy_b * vecb(pD,grid)
    for iLS in 1:num.nLS
        ∇ϕ_x .+= opC_u.Gx[iLS] * veci(pD,grid,iLS+1)
        ∇ϕ_y .+= opC_v.Gy[iLS] * veci(pD,grid,iLS+1)
    end

    iMu = Diagonal(inv_weight_eps2.(num.epsilon_mode,num.epsilon_vol,opC_u.M.diag))
    iMv = Diagonal(inv_weight_eps2.(num.epsilon_mode,num.epsilon_vol,opC_v.M.diag))
    ∇ϕ_x = iMu * ∇ϕ_x
    ∇ϕ_y = iMv * ∇ϕ_y
    
    grd_x .= reshape(veci(∇ϕ_x,grid_u,1), grid_u)
    grd_y .= reshape(veci(∇ϕ_y,grid_v,1), grid_v)

    printstyled(color=:magenta, @sprintf "\n grad min max x %.2e %.2e y %.2e %.2e\n" minimum(grd_x) maximum(grd_x) minimum(grd_y) maximum(grd_y))
end


"""
Adapt timestep based on velocity, viscosity, ...
# [Kang 2000, “A Boundary Condition Capturing Method for Multiphase Incompressible Flow”](https://doi.org/10.1023/A:1011178417620) 
"""
function adapt_timestep!(num, phL, phS, grid_u, grid_v,adapt_timestep_mode)
    @unpack rho1,rho2,mu1,mu2,grav_x,grav_y,CFL=num

    #With grid_u
    min_spacing_x = minimum(grid_u.dx)
    min_spacing_y = minimum(grid_u.dy)

    # min_spacing_xyz = min(min_spacing_x,min_spacing_y)

    # τ = min(CFL*Δ^2*Re, CFL*Δ/max(abs.(V)..., abs.(phL.u)..., abs.(phL.v)..., abs.(phS.u)..., abs.(phS.v)...))

    #TODO V ?
    # vel = max(abs.(V)..., abs.(phL.u)..., abs.(phL.v)..., abs.(phS.u)..., abs.(phS.v)...)

    if adapt_timestep_mode==1
        c_conv = max(abs.(phL.u)./grid_u.dx..., abs.(phL.v)./grid_v.dy..., abs.(phS.u)./grid_u.dx..., abs.(phS.v)./grid_v.dy...)

        c_visc = max(mu1/rho1, mu2/rho2) * (2.0/min_spacing_x^2 + 2.0/min_spacing_y^2 ) 
        
        c_grav = sqrt(max(abs(grav_x)/min_spacing_x, abs(grav_y)/min_spacing_y))

        #TODO capillary timestep
        c_surf = 0.0

        # c_surf = sigma*kappa / (min(rho1,rho2)*min_spacing_xyz^2)

        c_tot = (c_conv+c_visc)+ sqrt((c_conv+c_visc)^2+4*c_grav^2+4*c_surf^2 )

        return 2*CFL/c_tot
    elseif adapt_timestep_mode==2
        c_conv = max(abs.(phL.u)./grid_u.dx..., abs.(phL.v)./grid_v.dy..., abs.(phS.u)./grid_u.dx..., abs.(phS.v)./grid_v.dy...)

        return CFL/c_conv
    elseif adapt_timestep_mode==3
        return num.dt0
    end

    # printstyled(color=:green, @sprintf "\n rho1 %.2e rho2 %.2e mu1 %.2e mu2 %.2e\n" rho1 rho2 mu1 mu2) 
    # printstyled(color=:green, @sprintf "\n dx %.2e dy %.2e \n" min_spacing_x min_spacing_x) 

    
    # max(mu1/rho1, mu2/rho2) * (2.0/min_spacing_x^2 + 2.0/min_spacing_y^2 )

    #  printstyled(color=:green, @sprintf "\n c_conv %.2e c_visc %.2e c_grad %.2e\n" c_conv c_visc c_grav)
    #  printstyled(color=:green, @sprintf "\n CFL : %.2e dt : %.2e\n" CFL num.τ)
end




# if ls_advection
#     printstyled(color=:red, @sprintf "\n update operators solve_poisson_variable_coeff!\n")
    
#     # update_all_ls_data(num, grid, grid_u, grid_v, BC_int, periodic_x, periodic_y, false)

#     # Mass matrices
#     M.diag .= vec(grid.LS[1].geoL.dcap[:,:,5])
#     Mx = zeros(grid.ny,grid.nx+1)
#     for II in grid.ind.all_indices
#         Mx[II] = grid.LS[1].geoL.dcap[II,8]
#     end
#     for II in grid.ind.b_right[1]
#         Mx[δx⁺(II)] = grid.LS[1].geoL.dcap[II,10]
#     end
#     My = zeros(grid.ny+1,grid.nx)
#     for II in grid.ind.all_indices
#         My[II] = grid.LS[1].geoL.dcap[II,9]
#     end
#     for II in grid.ind.b_top[1]
#         My[δy⁺(II)] = grid.LS[1].geoL.dcap[II,11]
#     end
    
#     # iMx.diag .= 1. ./ (vec(Mx) .+ eps(0.01))
#     # iMy.diag .= 1. ./ (vec(My) .+ eps(0.01))
#     # iMx = Diagonal(inv_weight_eps2.(num.epsilon_mode,num.epsilon_vol,vec(Mx)))
#     # iMy = Diagonal(inv_weight_eps2.(num.epsilon_mode,num.epsilon_vol,vec(My)))
    
#     iMx.diag .= inv_weight_eps2.(num.epsilon_mode,num.epsilon_vol,vec(Mx))
#     iMy.diag .= inv_weight_eps2.(num.epsilon_mode,num.epsilon_vol,vec(My))


#     # Discrete gradient and divergence operators
#     divergence_B!(BxT, ByT, grid.LS[1].geoL.dcap, ny, grid.ind.all_indices)
#     mat_assign!(Bx, sparse(-BxT'))
#     mat_assign!(By, sparse(-ByT'))

#     # Matrices for interior BCs
#     for iLS in 1:num.nLS
#         bc_matrix!(grid, Hx[iLS], Hy[iLS], grid.LS[1].geoL.dcap, grid.LS[1].geoL.dcap, ny, grid.ind.all_indices)

#         mat_assign_T!(HxT[iLS], sparse(Hx[iLS]'))
#         mat_assign_T!(HyT[iLS], sparse(Hy[iLS]'))

#         periodic_bcs!(grid, Bx, By, Hx[iLS], Hy[iLS], periodic_x, periodic_y)

#         χx = (grid.LS[1].geoL.dcap[:,:,3] .- grid.LS[1].geoL.dcap[:,:,1]) .^ 2
#         χy = (grid.LS[1].geoL.dcap[:,:,4] .- grid.LS[1].geoL.dcap[:,:,2]) .^ 2
#         χ[iLS].diag .= sqrt.(vec(χx .+ χy))
#     end
#     mat_assign!(BxT, sparse(-Bx'))
#     mat_assign!(ByT, sparse(-By'))

#     # Matrices for borders BCs
#     set_boundary_indicator!(grid, grid.LS[1].geoL, grid.LS[1].geoL, op)
#     mass_matrix_borders!(num, grid.ind, op.iMx_b, op.iMy_b, op.iMx_bd, op.iMy_bd, grid.LS[1].geoL.dcap, ny)
#     bc_matrix_borders!(grid, grid.ind, op.Hx_b, op.Hy_b, grid.LS[1].geoL.dcap)
#     mat_assign_T!(op.HxT_b, sparse(op.Hx_b'))
#     mat_assign_T!(op.HyT_b, sparse(op.Hy_b'))
#     periodic_bcs_borders!(grid, op.Hx_b, op.Hy_b, periodic_x, periodic_y)
# end

"""
    test_filer_concentration!
    filter concentration to see to what extent abnormal electrolyte concentration (due to small cells) is responsible for abnormal potential

"""
function test_filter_concentration!(num,grid,scalD,scal0)

    for i in 1:num._nLS
        print("\n iLS ",i)

        for II in grid.ind.all_indices
                          
            pII = lexicographic(II, grid.ny)

            scal = veci(scalD, grid,i)[pII]

            if (abs(scal-scal0)/scal0 >1e-2) && (grid.LS[end].geoL.cap[II,5] > num.ϵ) && (scal >0.0) 

                # if i ==1 || (grid.LS[i-1].geoL.cap[II,5] > num.ϵ) #TODO intfc value detect non null

                print("\n II ",II," s ",scal," s0 ",scal0," ratio ", abs(scal-scal0)/scal0)
                print("\n replace scal by scal0 (test)")
                scal = scal0
            end

        end
    end

    #TODO vecb

    # boundaries_idx = [grid.ind.b_left[1], grid.ind.b_bottom[1], grid.ind.b_right[1], grid.ind.b_top[1]]



    # left2 = vcat(grid.ind.all_indices[2,2], grid.ind.all_indices[2:end-1,2], grid.ind.all_indices[end-1,2])
    # bottom2 = vcat(grid.ind.all_indices[2,2], grid.ind.all_indices[2,2:end-1], grid.ind.all_indices[2,end-1])
    # right2 = vcat(grid.ind.all_indices[2,end-1], grid.ind.all_indices[2:end-1,end-1], grid.ind.all_indices[end-1,end-1])
    # top2 = vcat(grid.ind.all_indices[end-1,2], grid.ind.all_indices[end-1,2:end-1], grid.ind.all_indices[end-1,end-1])
    # boundaries2 = [left2, bottom2, right2, top2]

    # left3 = vcat(grid.ind.all_indices[3,3], grid.ind.all_indices[2:end-1,3], grid.ind.all_indices[end-2,3])
    # bottom3 = vcat(grid.ind.all_indices[3,3], grid.ind.all_indices[3,2:end-1], grid.ind.all_indices[3,end-2])
    # right3 = vcat(grid.ind.all_indices[3,end-2], grid.ind.all_indices[2:end-1,end-2], grid.ind.all_indices[end-2,end-2])
    # top3 = vcat(grid.ind.all_indices[end-2,3], grid.ind.all_indices[end-2,2:end-1], grid.ind.all_indices[end-2,end-2])
    # boundaries3 = [left3, bottom3, right3, top3]

    # boundaries_t = [grid.ind.left, grid.ind.bottom, grid.ind.right, grid.ind.top]


    # for (i, (idx, idx2, idx3, xy)) in enumerate(zip(boundaries_idx, boundaries2, boundaries3, direction))
    #     pks, _ = findminima(abs.(u[idx]))
    #     # if is_neumann(boundaries_t[i])

    #     #     for (II, JJ) in zip(idx, idx2)
    #     #         pII = lexicographic(II, grid.ny)
    #     #         pJJ = lexicographic(JJ, grid.ny)

    #     #         A[pII,:] .= 0.0
    #     #         A[pII,pII] = 1.0
    #     #         A[pII,pJJ] = -1.0
    #     #         B[pII,:] .= 0.0
    #     #     end


    # print("left ",grid.ind.b_left)
    #TODO borders

# vecb_L(a,g::G) where {G<:Grid} = @view vecb(a, g)[1:g.ny]
# vecb_B(a,g::G) where {G<:Grid} = @view vecb(a, g)[g.ny+1:g.ny+g.nx]
# vecb_R(a,g::G) where {G<:Grid} = @view vecb(a, g)[g.ny+g.nx+1:2*g.ny+g.nx]
# vecb_T(a,g::G) where {G<:Grid} = @view vecb(a, g)[2*g.ny+g.nx+1:2*g.ny+2*g.nx]

    # left
    for j in 1:grid.ny
        II = CartesianIndex(j,1)
        scal = vecb_L(scalD, grid)[j]
        if abs(scal-scal0)/scal0 >1e-2
            print("\n II ",II," s ",scal," s0 ",scal0," ratio ", abs(scal-scal0)/scal0)
            print("\n replace scal by scal0 (test)")
            scal = scal0
        end
    end

    # bottom
    for i in 1:grid.nx
        II = CartesianIndex(1,i)
        scal = vecb_B(scalD, grid)[i]
        if abs(scal-scal0)/scal0 >1e-2
            print("\n II ",II," s ",scal," s0 ",scal0," ratio ", abs(scal-scal0)/scal0)
            print("\n replace scal by scal0 (test)")
            scal = scal0
        end
    end

    # right
    for j in 1:grid.ny
        II = CartesianIndex(j,grid.nx)
        scal = vecb_R(scalD, grid)[j]
        if abs(scal-scal0)/scal0 >1e-2
            print("\n II ",II," s ",scal," s0 ",scal0," ratio ", abs(scal-scal0)/scal0)
            print("\n replace scal by scal0 (test)")
            scal = scal0
        end
    end

    # top
    for i in 1:grid.nx
        II = CartesianIndex(grid.ny,i)
        scal = vecb_T(scalD, grid)[i]
        if abs(scal-scal0)/scal0 >1e-2
            print("\n II ",II," s ",scal," s0 ",scal0," ratio ", abs(scal-scal0)/scal0)
            print("\n replace scal by scal0 (test)")
            scal = scal0
        end
    end


end





"""
Prepare boundary conditions

!!! danger "Sign of BC"
    * BC.left.val .* ones(ny)
    * BC.bottom.val .* ones(nx)
    * BC.right.val .* ones(ny)
    * BC.bottom.val .* ones(nx)
    In the current implementation, the sign needs to be checked.
    Dirichlet: a1 = -1
    Neumann: b =1
    Robin: a1 = -1 b = 1
    a0 : BC value 
"""
function set_borders_poisson!(grid, cl, u, a0, a1, b, BC, n_ext)
    @unpack nx, ny, x, y, dx, dy, ind = grid

    idx = 1:ny
    @inbounds a0[idx] .= BC.left.val .* ones(ny)
    if is_dirichlet(BC.left)
        @inbounds a1[idx] .= 1.0
        @inbounds b[idx] .= 0.0
    elseif is_neumann(BC.left)
        @inbounds a1[idx] .= 0.0
        @inbounds b[idx] .= 1.0
    elseif is_robin(BC.left)
        @inbounds a1[idx] .= 1.0
        @inbounds b[idx] .= 1.0
    elseif is_periodic(BC.left)
        nothing
    elseif is_navier(BC.left)
        @inbounds a1[idx] .= 1.0
        @inbounds b[idx] .= BC.left.λ
    elseif is_navier_cl(BC.left)
        idx_cl = intersect(CartesianIndices((idx,1)), cl)
        idx_no = setdiff(CartesianIndices((idx,1)), cl)
        ϵb = zeros(grid)
        ϵb[idx_cl] .= n_ext .* dx[idx_cl]

        @inbounds @threads for II in idx_cl
            @inbounds a1[II[1]] = 1.0
            @inbounds b[II[1]] = BC.left.λ * bell_function2(u[II], ϵb[II])
        end
        @inbounds @threads for II in idx_no
            @inbounds a1[II[1]] = 1.0
            @inbounds b[II[1]] = 0.0
        end
    elseif is_gnbc(BC.left)
        # idx_cl = intersect(CartesianIndices((idx,1)), cl)
        # idx_no = setdiff(CartesianIndices((idx,1)), cl)

        # bell = bell_function(grid, BC.ϵ)
        # @inbounds a0[idx] .+= bell[idx] .* BC.σ ./ BC.μ .* (cos(θd) .- cos(BC.θe))
        # @inbounds @threads for II in idx_cl
        #     @inbounds a1[II[1]] = 1.0
        #     @inbounds b[II[1]] = BC.left.λ
        # end
        # @inbounds @threads for II in idx_no
        #     @inbounds a1[II[1]] = 1.0
        #     @inbounds b[II[1]] = 0.0
        # end
        nothing
    else
        @error ("Not implemented yet")
    end

    _idx = ny+1:ny+nx
    idx = 1:nx
    @inbounds a0[_idx] .= BC.bottom.val .* ones(nx)
    if is_dirichlet(BC.bottom)
        @inbounds a1[_idx] .= 1.0
        @inbounds b[_idx] .= 0.0
    elseif is_neumann(BC.bottom)
        @inbounds a1[_idx] .= 0.0
        @inbounds b[_idx] .= 1.0
    elseif is_robin(BC.bottom)
        @inbounds a1[_idx] .= 1.0
        @inbounds b[_idx] .= 1.0
    elseif is_periodic(BC.bottom)
        nothing
    elseif is_navier(BC.bottom)
        @inbounds a1[_idx] .= 1.0
        @inbounds b[_idx] .= BC.bottom.λ
    elseif is_navier_cl(BC.bottom)
        idx_cl = intersect(CartesianIndices((1,idx)), cl)
        idx_no = setdiff(CartesianIndices((1,idx)), cl)
        ϵb = zeros(grid)
        ϵb[idx_cl] .= n_ext .* dy[idx_cl]

        @inbounds @threads for II in idx_cl
            @inbounds a1[II[2]+ny] = 1.0
            @inbounds b[II[2]+ny] = BC.bottom.λ * bell_function2(u[II], ϵb[II])
        end
        @inbounds @threads for II in idx_no
            @inbounds a1[II[2]+ny] = 1.0
            @inbounds b[II[2]+ny] = 0.0
        end
    elseif is_gnbc(BC.bottom)
        # @inbounds a0[idx] .+= bell[idx] .* BC.σ ./ BC.μ .* (cos(θd) .- cos(BC.θe))
        # @inbounds a1[intersect(idx, cl)] .= 1.0
        # @inbounds a1[setdiff(idx, cl)] .= 1.0
        # @inbounds b[intersect(idx, cl)] .= BC.bottom.λ
        # @inbounds b[setdiff(idx, cl)] .= 0.0
        nothing
    else
        @error ("Not implemented yet")
    end

    _idx = ny+nx+1:2*ny+nx
    idx = 1:ny
    @inbounds a0[_idx] .= BC.right.val .* ones(ny)
    if is_dirichlet(BC.right)
        @inbounds a1[_idx] .= 1.0
        @inbounds b[_idx] .= 0.0
    elseif is_neumann(BC.right)
        @inbounds a1[_idx] .= 0.0
        @inbounds b[_idx] .= 1.0
    elseif is_robin(BC.right)
        @inbounds a1[_idx] .= 1.0
        @inbounds b[_idx] .= 1.0
    elseif is_periodic(BC.right)
        nothing
    elseif is_navier(BC.right)
        @inbounds a1[_idx] .= 1.0
        @inbounds b[_idx] .= BC.right.λ
    elseif is_navier_cl(BC.right)
        idx_cl = intersect(CartesianIndices((idx,nx)), cl)
        idx_no = setdiff(CartesianIndices((idx,nx)), cl)
        ϵb = zeros(grid)
        ϵb[idx_cl] .= n_ext .* dx[idx_cl]

        @inbounds @threads for II in idx_cl
            @inbounds a1[II[1]+ny+nx] = 1.0
            @inbounds b[II[1]+ny+nx] = BC.right.λ * bell_function2(u[II], ϵb[II])
        end
        @inbounds @threads for II in idx_no
            @inbounds a1[II[1]+ny+nx] = 1.0
            @inbounds b[II[1]+ny+nx] = 0.0
        end
    elseif is_gnbc(BC.right)
        # @inbounds a0[idx] .+= bell[idx] .* BC.σ ./ BC.μ .* (cos(θd) .- cos(BC.θe))
        # @inbounds a1[intersect(idx, cl)] .= 1.0
        # @inbounds a1[setdiff(idx, cl)] .= 1.0
        # @inbounds b[intersect(idx, cl)] .= BC.right.λ
        # @inbounds b[setdiff(idx, cl)] .= 0.0
        nothing
    else
        @error ("Not implemented yet")
    end

    _idx = 2*ny+nx+1:2*ny+2*nx
    idx = 1:nx
    @inbounds a0[_idx] .= BC.top.val .* ones(nx)
    if is_dirichlet(BC.top)
        @inbounds a1[_idx] .= 1.0
        @inbounds b[_idx] .= 0.0
    elseif is_neumann(BC.top)
        @inbounds a1[_idx] .= 0.0
        @inbounds b[_idx] .= 1.0
    elseif is_robin(BC.top)
        @inbounds a1[_idx] .= 1.0
        @inbounds b[_idx] .= 1.0
    elseif is_periodic(BC.top)
        nothing
    elseif is_navier(BC.top)
        @inbounds a1[_idx] .= 1.0
        @inbounds b[_idx] .= BC.top.λ
    elseif is_navier_cl(BC.top)
        idx_cl = intersect(CartesianIndices((1,idx)), cl)
        idx_no = setdiff(CartesianIndices((1,idx)), cl)
        ϵb = zeros(grid)
        ϵb[idx_cl] .= n_ext .* dy[idx_cl]

        @inbounds @threads for II in idx_cl
            @inbounds a1[II[2]+2*ny+nx] = 1.0
            @inbounds b[II[2]+2*ny+nx] = BC.top.λ * bell_function2(u[II], ϵb[II])
        end
        @inbounds @threads for II in idx_no
            @inbounds a1[II[2]+2*ny+nx] = 1.0
            @inbounds b[II[2]+2*ny+nx] = 0.0
        end
    elseif is_gnbc(BC.top)
        # @inbounds a0[idx] .+= bell[idx] .* BC.σ ./ BC.μ .* (cos(θd) .- cos(BC.θe))
        # @inbounds a1[intersect(idx, cl)] .= 1.0
        # @inbounds a1[setdiff(idx, cl)] .= 1.0
        # @inbounds b[intersect(idx, cl)] .= BC.top.λ
        # @inbounds b[setdiff(idx, cl)] .= 0.0
        nothing
    else
        @error ("Not implemented yet")
    end

    return nothing
end



@doc raw"""
    solve_poisson(bc_type, num, grid, a0, opC, opC_u, opC_v, A, L, bc_L, bc_L_b, BC, ls_advection)

Solves the Poisson equation for a given set of boundary conditions and grid configuration.

### Parameters
- `bc_type`: A vector specifying the type of boundary conditions for interfaces.
- `num`: A structure containing numerical parameters for the simulation.
- `grid`: A structure defining the grid configuration.
- `a0`: Initial values for the solution.
- `opC`: Operator configuration containing various grid and boundary parameters.
- `opC_u`: Operator configuration for the u-component.
- `opC_v`: Operator configuration for the v-component.
- `A`: The matrix to be populated with the Poisson equation coefficients.
- `L`: The right-hand side matrix for the Poisson equation.
- `bc_L`: Boundary conditions for the left side of the domain.
- `bc_L_b`: Boundary conditions for the left side of the outer boundaries.
- `BC`: boundary conditions for the borders of the domain.
- `ls_advection`: Boolean flag indicating whether advection terms are included.

### Returns
- `rhs`: The right-hand side vector for the Poisson equation.

"""
function solve_poisson(
    bc_type, num, grid, a0, opC, opC_u, opC_v,
    A, L, bc_L, bc_L_b, BC,
    ls_advection)
    @unpack Bx, By, Hx, Hy, HxT, HyT, χ, M, iMx, iMy, Hx_b, Hy_b, HxT_b, HyT_b, iMx_b, iMy_b, iMx_bd, iMy_bd, χ_b = opC

    ni = grid.nx * grid.ny
    nb = 2 * grid.nx + 2 * grid.ny

    rhs = fnzeros(grid, num)

    a0_b = zeros(nb)
    _a1_b = zeros(nb)
    _b_b = zeros(nb)
    for iLS in 1:num.nLS
        set_borders_poisson!(grid, grid.LS[iLS].cl, grid.LS[iLS].u, a0_b, _a1_b, _b_b, BC, num.n_ext_cl)
    end
    a1_b = Diagonal(vec(_a1_b))
    b_b = Diagonal(vec(_b_b))

    if ls_advection
        # Poisson equation
        A[1:ni,1:ni] = pad(L, -4.0)
        A[1:ni,end-nb+1:end] = bc_L_b

        # Boundary conditions for outer boundaries
        A[end-nb+1:end,1:ni] = b_b * (HxT_b * iMx_b' * Bx .+ HyT_b * iMy_b' * By) 
        A[end-nb+1:end,end-nb+1:end] = pad(b_b * (HxT_b * iMx_bd * Hx_b .+ HyT_b * iMy_bd * Hy_b) .+ χ_b * a1_b,- 4.0)
    end

    for iLS in 1:num.nLS
        if ls_advection
            if is_dirichlet(bc_type[iLS])
                __a1 = 1.0
                __a2 = 0.0
                __b = 0.0
            elseif is_neumann(bc_type[iLS])
                __a1 = 0.0
                __a2 = 0.0
                __b = 1.0
            elseif is_robin(bc_type[iLS])
                __a1 = 1.0
                __a2 = 0.0
                __b = 1.0
            elseif is_fs(bc_type[iLS])
                __a1 = 0.0
                __a2 = 1.0
                __b = 0.0
            elseif is_wall_no_slip(bc_type[iLS])
                __a1 = 0.0
                __a2 = 0.0
                __b = 1.0
            elseif is_navier(bc_type[iLS])
                __a1 = 0.0
                __a2 = 0.0
                __b = 1.0
            elseif is_navier_cl(bc_type[iLS])
                __a1 = 0.0
                __a2 = 0.0
                __b = 1.0
            else
                __a1 = 0.0
                __a2 = 0.0
                __b = 1.0
            end
    
            _a1 = ones(grid) .* __a1
            a1 = Diagonal(vec(_a1))
            _a2 = ones(grid) .* __a2
            a2 = Diagonal(vec(_a2))
            _b = ones(grid) .* __b
            b = Diagonal(vec(_b))

            fs_mat = HxT[iLS] * Hx[iLS] .+ HyT[iLS] * Hy[iLS]

            sb = iLS*ni+1:(iLS+1)*ni
            
            # Poisson equation
            A[1:ni,sb] = bc_L[iLS]
            # Boundary conditions for inner boundaries
            A[sb,1:ni] = b * (HxT[iLS] * iMx * Bx .+ HyT[iLS] * iMy * By)
            # Contribution to Neumann BC from other boundaries
            for i in 1:num.nLS
                if i != iLS
                    A[sb,i*ni+1:(i+1)*ni] = b * (HxT[iLS] * iMx * Hx[i] .+ HyT[iLS] * iMy * Hy[i])
                end
            end
            A[sb,sb] = pad(
                b *(HxT[iLS] * iMx * Hx[iLS] .+ HyT[iLS] * iMy * Hy[iLS]) .+ χ[iLS] * a1 .-
                a2 * Diagonal(diag(fs_mat)), -4.0
            )
            #TODO was +b... here in set_poisson 
            A[sb,end-nb+1:end] = b * (HxT[iLS] * iMx_b * Hx_b .+ HyT[iLS] * iMy_b * Hy_b) #why was it different in set_poisson ?
            # Boundary conditions for outer boundaries
            A[end-nb+1:end,sb] = b_b * (HxT_b * iMx_b' * Hx[iLS] .+ HyT_b * iMy_b' * Hy[iLS])
        end

        veci(rhs,grid,iLS+1) .= +χ[iLS] * vec(a0[iLS]) #was - in set_poisson
    end

    vecb(rhs,grid) .= +χ_b * vec(a0_b) #was - in set_poisson
    
    return rhs
end


@doc raw"""
    set_poisson_variable_coeff(bc_type, num, grid, a0, opC, opC_u, opC_v, A, L, bc_L, bc_L_b, BC, ls_advection)

Solves the Poisson equation for a given set of boundary conditions and grid configuration.

### Parameters
- `bc_type`: A vector specifying the type of boundary conditions for interfaces.
- `num`: A structure containing numerical parameters for the simulation.
- `grid`: A structure defining the grid configuration.
- `a0`: Initial values for the solution.
- `opC`: Operator configuration containing various grid and boundary parameters.
- `opC_u`: Operator configuration for the u-component.
- `opC_v`: Operator configuration for the v-component.
- `A`: The matrix to be populated with the Poisson equation coefficients.
- `L`: The right-hand side matrix for the Poisson equation.
- `bc_L`: Boundary conditions for the left side of the domain.
- `bc_L_b`: Boundary conditions for the left side of the outer boundaries.
- `BC`: boundary conditions for the borders of the domain.
- `ls_advection`: Boolean flag indicating whether advection terms are included.

### Returns
- `rhs`: The right-hand side vector for the Poisson equation.

"""
function solve_poisson_variable_coeff!(num::Numerical{Float64, Int64},
    grid::Mesh{Flower.GridCC, Float64, Int64},
    grid_u::Mesh{Flower.GridFCx, Float64, Int64},
    grid_v::Mesh{Flower.GridFCy, Float64, Int64},
    opC::Operators{Float64, Int64},
    A::SparseMatrixCSC{Float64, Int64},
    rhs::Array{Float64, 1},
    a0::Array{Float64, 2},
    a1::SparseMatrixCSC{Float64, Int64},
    BC::BoundariesInt,
    ph::Phase{Float64},
    elec_cond,
    coeffD::Array{Float64, 1},
    coeffDu::Array{Float64, 2},
    coeffDv::Array{Float64, 2},
    i_butler::Array{Float64, 1},
    ls_advection::Bool,
    heat::Bool)

    @unpack Bx, By, Hx, Hy, HxT, HyT, χ, M, iMx, iMy, Hx_b, Hy_b, HxT_b, HyT_b, iMx_b, iMy_b, iMx_bd, iMy_bd, χ_b = opC
    @unpack BxT, ByT,tmp_x, tmp_y = opC

    ni = grid.nx * grid.ny
    nb = 2 * grid.nx + 2 * grid.ny

    #TODO reset zero
    rhs .= 0.0
    coeffDu .= 0.0
    coeffDv .= 0.0
    A .= 0.0
    a0 .= 0.0

    a0_b = zeros(nb)
    _a1_b = zeros(nb)
    _b_b = zeros(nb)
    for iLS in 1:num.nLS
        set_borders!(grid, grid.LS[iLS].cl, grid.LS[iLS].u, a0_b, _a1_b, _b_b, BC, num.n_ext_cl)
    end
    a1_b = Diagonal(vec(_a1_b))
    b_b = Diagonal(vec(_b_b))

    #interpolate coefficient
    coeffD_borders = vecb(coeffD,grid)
    interpolate_scalar!(grid, grid_u, grid_v, reshape(veci(coeffD,grid,1), grid), coeffDu, coeffDv)

    # print("\n coeff ",minimum(coeffDu)," ",maximum(coeffDu)," ",minimum(coeffDv)," ",maximum(coeffDv)," ",minimum(coeffD_borders)," ",maximum(coeffD_borders))

    coeffDx_bulk = veci(coeffDu,grid_u)
    coeffDy_bulk = veci(coeffDv,grid_v)

    # mat_coeffDx = Diagonal(vec(coeffDx_bulk)) # coeffDx_bulk is a 2d matrix with shape (grid_u.ny, grid_u.nx), multiplies Bx
    # mat_coeffDy = Diagonal(vec(coeffDy_bulk)) # coeffDx_bulk is a 2d matrix with shape (grid_v.ny, grid_v.nx), multiplies By

    # print("\n sizes",size(coeffDx_bulk)) #
    # print("\n sizes coeffDy_bulk ",size(coeffDy_bulk)) #

    mat_coeffDx = Diagonal(vec(coeffDu)) # coeffDx_bulk is a 2d matrix with shape (grid_u.ny, grid_u.nx), multiplies Bx
    mat_coeffDy = Diagonal(vec(coeffDv)) # coeffDx_bulk is a 2d matrix with shape (grid_v.ny, grid_v.nx), multiplies By

    # Laplacian: bulk 
    # L = BxT * iMx * Bx
    # L is of size nx*ny
    # tmp_x ((nx+1)*ny , nx*ny)
    # BxT (nx*ny, (nx+1)*ny )
    # Bx ((nx+1)*ny, nx*ny)
    # iMx ((nx+1)*ny, (nx+1)*ny )
    # iMx_b (nx+1)*ny,2*nx +2*ny
    # Hx_b  2*nx +2*ny,2*nx +2*ny

    # mat_coeffDx_b should be of size ((nx+1)*ny,(nx+1)*ny)

    old_method = false
    # old_method = true

    if old_method #(wrong)

        # Laplacian
        # diag so ok if not order mat_coeffDx iMx Bx ? No, wrong equation
        mul!(tmp_x, iMx, mat_coeffDx * Bx)
        L = BxT * tmp_x
        mul!(tmp_y, iMy, mat_coeffDy * By)
        L = L .+ ByT * tmp_y


        # mul!(tmp_x, mat_coeffDx * iMx, Bx)
        # L = BxT * tmp_x
        # mul!(tmp_y, mat_coeffDy * iMy, By)
        # L = L .+ ByT * tmp_y

        mat_coeffDx_b = Diagonal(vec(coeffD_borders)) # is a 1d vector with shape (2grid.ny + 2grid.nx), multiplies Hx_b and Hy_b
        
        #Boundary for Laplacian
        bc_L_b = (BxT * iMx_b * mat_coeffDx_b *Hx_b .+ ByT * iMy_b * mat_coeffDx_b *Hy_b)

    else
        
        #TODO wall control volume kappa for bulk vs bulk control volume kappa different  
        mul!(tmp_x, mat_coeffDx * iMx, Bx)
        L = BxT * tmp_x
        mul!(tmp_y, mat_coeffDy * iMy, By)
        L = L .+ ByT * tmp_y

        coeffDu_border = copy(coeffDu)
        coeffDv_border = copy(coeffDv)

        # Interpolate conductivity at center of control volumes for potential gradient at the border
        # interpolate_scalar_to_staggered_u_v_grids_at_border!(num,grid,coeffD,coeffDu,coeffDv)

        interpolate_scalar_to_staggered_u_v_grids_at_border!(num,grid,coeffD,coeffDu_border,coeffDv_border)

        coeffDx_border = veci(coeffDu_border,grid_u)
        coeffDy_border = veci(coeffDv_border,grid_v)

        mat_coeffDx_b = Diagonal(vec(coeffDx_border)) 
        mat_coeffDy_b = Diagonal(vec(coeffDy_border))

        #Boundary for Laplacian
        bc_L_b = (BxT * mat_coeffDx_b * iMx_b * Hx_b .+ ByT * mat_coeffDy_b * iMy_b  * Hy_b)

    end 

    #border


      
    if ls_advection
        # Poisson equation
        A[1:ni,1:ni] = pad(L, 4.0)
        A[1:ni,end-nb+1:end] = bc_L_b

        # Boundary conditions for outer boundaries
        # A[end-nb+1:end,1:ni] = b_b * (HxT_b * iMx_b' * mat_coeffDx * Bx .+ HyT_b * iMy_b' * mat_coeffDy * By)
        # A[end-nb+1:end,end-nb+1:end] = pad(b_b * (HxT_b * iMx_bd * mat_coeffDx_b * Hx_b .+ HyT_b * iMy_bd * mat_coeffDx_b * Hy_b) .+ χ_b * a1_b, -4.0)
        
        printstyled(color=:red, @sprintf "\n test coeff poisson TODO for iLS too:\n")

        A[end-nb+1:end,1:ni] = b_b * (HxT_b * iMx_b' * Bx .+ HyT_b * iMy_b' * By)
        A[end-nb+1:end,end-nb+1:end] = pad(b_b * (HxT_b * iMx_bd  * Hx_b .+ HyT_b * iMy_bd * Hy_b) .+ χ_b * a1_b, -4.0)
        

    end

    for iLS in 1:num.nLS

        a0 .= 0.0 #reset

        if ls_advection
            if is_dirichlet(BC.LS[iLS])
                __a0 = BC.LS[iLS].val
                __a1 = -1.0
                __a2 = 0.0
                __b = 0.0
            elseif is_neumann(BC.LS[iLS])
                __a0 = BC.LS[iLS].val
                __a1 = 0.0
                __a2 = 0.0
                __b = 1.0
            elseif is_robin(BC.LS[iLS])
                __a0 = BC.LS[iLS].val
                __a1 = -1.0
                __a2 = 0.0
                __b = 1.0
            elseif is_fs(BC.LS[iLS])
                print("error not implemented set_poisson_variable_coeff",BC.LS[iLS])
                @error ("error set_poisson_variable_coeff")

                __a1 = 0.0
                __a2 = 1.0
                __b = 0.0
            elseif is_wall_no_slip(BC.LS[iLS])
                print("error not implemented set_poisson_variable_coeff",BC.LS[iLS])
                @error ("error set_poisson_variable_coeff")

                __a1 = 0.0
                __a2 = 0.0
                __b = 1.0
            elseif is_navier(BC.LS[iLS])
                print("error not implemented set_poisson_variable_coeff",BC.LS[iLS])
                @error ("error set_poisson_variable_coeff")

                __a1 = 0.0
                __a2 = 0.0
                __b = 1.0
            elseif is_navier_cl(BC.LS[iLS])
                print("error not implemented set_poisson_variable_coeff",BC.LS[iLS])
                @error ("error set_poisson_variable_coeff")

                __a1 = 0.0
                __a2 = 0.0
                __b = 1.0
            else
                __a1 = 0.0
                __a2 = 0.0
                __b = 1.0
            end
    
            if num.nLS == 1
                # Flags with BCs
                a0 .= __a0
                # a0 = ones(grid) .* __a0
            
            else 
                iLS_elec = 2
                
                if iLS == iLS_elec 

                    #use Butler-Volmer, supposing the interfacial potential is acceptable and phi = phi_ele1 in metal 
                    # for conductivity, use interfacial value or bulk in corresponding cell

                    # TODO -(-i/kappa) in Flower ? so i_butler not -i_butler

                    # if num.bulk_conductivity == 0
                        
                    # elseif num.bulk_conductivity == 1
                    #     @error ("error elseif num.bulk_conductivity == 1")

                    if num.bulk_conductivity == 2
                        # Recommended as long as cell merging not implemented:

                        #TODO remove reshape and use a mapping 
                        # butler_volmer_no_concentration_potential_Neumann!.(num,
                        # reshape(veci(ph.phi_eleD, grid,iLS+1),grid),
                        # reshape(veci(ph.trans_scalD[:,2],grid,iLS+1),grid),
                        # num.temperature0,
                        # a0)
                        # a0 .= butler_volmer_no_concentration_potential_Neumann.(num,
                        # reshape(veci(ph.phi_eleD, grid,iLS+1),grid),
                        # reshape(veci(ph.trans_scalD[:,2],grid,iLS+1),grid),
                        # num.temperature0)

                        for II in grid.LS[iLS].MIXED
                          

                            # a0[II] .= butler_volmer_no_concentration_potential_Neumann.(num,
                            # reshape(veci(ph.phi_eleD, grid,iLS+1),grid),
                            # reshape(veci(ph.trans_scalD[:,2],grid,iLS+1),grid),
                            # num.temperature0)

                            pII = lexicographic(II, grid.ny)

                            # if grid.LS[iLS].geoL.cap[II,5] < num.ϵ
                       
                            if grid.LS[end].geoL.cap[II,5] > num.ϵ #TODO clearer eps

                                a0[II] = butler_volmer_no_concentration_potential_Neumann.(num,
                                veci(ph.phi_eleD, grid,iLS+1)[pII],
                                veci(ph.trans_scalD[:,2],grid,iLS+1)[pII],
                                num.temperature0)

                                if veci(ph.trans_scalD[:,2],grid,iLS+1)[pII] < num.epsilon_concentration[2]
                                    a0[II] = butler_volmer_no_concentration_potential_Neumann.(num,
                                    # reshape(veci(ph.phi_eleD, grid,iLS+1),grid),
                                    veci(ph.phi_eleD, grid,iLS+1)[pII],
                                    ph.trans_scal[II,2],
                                    num.temperature0)
                                end

                                # print("\n II",II,"BC ", BC.LS[iLS].val)
                                # printstyled(color=:red, @sprintf "\n Butler %.2e %.2e \n" a0[II] reshape(veci(ph.trans_scalD[:,2],grid,iLS+1),grid)[II])


                                # a0[II] = butler_volmer_no_concentration_potential_Neumann.(num,
                                # reshape(veci(ph.phi_eleD, grid,iLS+1),grid)[II],
                                # reshape(veci(ph.trans_scalD[:,2],grid,iLS+1),grid)[II],
                                # num.temperature0)

                            else
                                print("\n grid.LS[end].geoL.cap[II,5] > num.ϵ ", (grid.LS[end].geoL.cap[II,5] > num.ϵ),(ph.trans_scal[II,2]<num.ϵ), " ",ph.trans_scal[II,2], " ",num.ϵ)
                                #use bulk conductivity of mixed cell
                                # butler_volmer_no_concentration_potential_Neumann!.(num,
                                # reshape(veci(ph.phi_eleD, grid,iLS+1),grid),
                                # ph.trans_scal[II,2],
                                # num.temperature0,
                                # a0[II]) #TODO if temperature solved temperature[II]

                                if ph.trans_scal[II,2]<num.epsilon_concentration[2] #inside bubble, do not solve, fill with 1 since concentration=0
                                    a0[II] = 1.0 
                                    print("\n zero scal II")
                                else
                                    a0[II] = butler_volmer_no_concentration_potential_Neumann.(num,
                                    veci(ph.phi_eleD, grid,iLS+1)[pII],
                                    ph.trans_scal[II,2],
                                    num.temperature0) #TODO if temperature solved temperature[II]
                                end #ph.trans_scal[II,2]<num.ϵ


                            end #grid.LS[end].geoL.cap[II,5] > num.ϵ: liquid cell

                            print("\n II ",II)
                            # printstyled(color=:red, @sprintf "\n grid.LS[end].geoL.cap[II,5] %.2e grid.LS[iLS].geoL.cap[II,1] %.2e grid.LS[iLS].geoL.cap[II,5] %.2e scal %.2e\n" grid.LS[end].geoL.cap[II,5] grid.LS[iLS].geoL.cap[II,1] grid.LS[iLS].geoL.cap[II,5] ph.trans_scal[II,2] a0[II])
                            printstyled(color=:red, @sprintf "\n grid.LS[end].geoL.cap[II,5] %.2e scaD %.2e scal %.2e phiD %.2e phi %.2e a0 %.2e \n" grid.LS[end].geoL.cap[II,5] veci(ph.trans_scalD[:,2],grid,iLS+1)[pII] ph.trans_scal[II,2] veci(ph.phi_eleD, grid,iLS+1)[pII] ph.phi_ele[II] a0[II])


                            # TODO
                            #Remove Nan when dividing by conductivity which may be null
                            # kill_dead_bc_left_wall!(vecb(elec_condD,grid), grid, iLS,1.0)
                                #Remove Nan when dividing by conductivity which may be null
                        end   #for
                    
                    else #bulk_conductivity!=2
                        a0 .= __a0

                    end #bulk_conductivity

                else #ilS==iLS_elec
                    a0 .= __a0
                end #ilS==iLS_elec

            end #num.nLS == 1


            # _a1 = ones(grid) .* __a1
            # a1 = Diagonal(vec(_a1))
            a1.nzval .= __a1

            _b = ones(grid) .* __b
            b = Diagonal(vec(_b))

            # a0_b = zeros(nb)
            # _a1_b = zeros(nb)
            # _b_b = zeros(nb)
            # set_borders!(grid, grid.LS[1].cl, grid.LS[1].u, a0_b, _a1_b, _b_b, BC, num.n_ext_cl)
            # a1_b = Diagonal(vec(_a1_b))
            # b_b = Diagonal(vec(_b_b))

           
            _a2 = ones(grid) .* __a2
            a2 = Diagonal(vec(_a2))
         
            fs_mat = HxT[iLS] * Hx[iLS] .+ HyT[iLS] * Hy[iLS]

            sb = iLS*ni+1:(iLS+1)*ni

            #interpolate conductivity coefficient for interface term
            #TODO multiple scalars : better to use div grad...?
            # coeffDu0 .= 0.0
            # coeffDv0 .= 0.0
            # interpolate_scalar!(grid, grid_u, grid_v, reshape(veci(coeffD,grid,iLS+1), grid), coeffDu0, coeffDv0)
            # mat_coeffDx_i = Diagonal(vec(coeffDu0)) # coeffDu is a 2d matrix with shape (grid_u.ny, grid_u.nx), multiplies Hx
            # mat_coeffDy_i = Diagonal(vec(coeffDv0)) # coeffDu is a 2d matrix with shape (grid_v.ny, grid_v.nx), multiplies Hy

            # print("\n test coeff")
            coeffDu .= 0.0
            coeffDv .= 0.0
            # TODO will not interpolate correctly, need to extend ...
            interpolate_scalar!(grid, grid_u, grid_v, reshape(veci(coeffD,grid,iLS+1), grid), coeffDu, coeffDv)
            mat_coeffDx_i = Diagonal(vec(coeffDu)) # coeffDu is a 2d matrix with shape (grid_u.ny, grid_u.nx), multiplies Hx
            mat_coeffDy_i = Diagonal(vec(coeffDv)) # coeffDu is a 2d matrix with shape (grid_v.ny, grid_v.nx), multiplies Hy


            
            # Poisson equation
            #Boundary for Laplacian from iLS
            # A[1:ni,sb] = BxT * iMx * mat_coeffDx_i *Hx[iLS] .+ ByT * iMy * mat_coeffDy_i *Hy[iLS]
            A[1:ni,sb] = BxT * mat_coeffDx_i * iMx * Hx[iLS] .+ ByT * mat_coeffDy_i * iMy * Hy[iLS]

            # Boundary conditions for inner boundaries
            A[sb,1:ni] = b * (HxT[iLS] * iMx  * Bx .+ HyT[iLS] * iMy * By) #or vec1
            # Contribution to Neumann BC from other boundaries
            for i in 1:num.nLS
                if i != iLS
                    A[sb,i*ni+1:(i+1)*ni] = b * (HxT[iLS] * iMx  * Hx[i] .+ HyT[iLS] * iMy * Hy[i])
                end
            end
            A[sb,sb] = pad(
                b * (HxT[iLS] * iMx * Hx[iLS] .+ HyT[iLS] * iMy * Hy[iLS]) .+ χ[iLS] * a1 .+
                a2 * Diagonal(diag(fs_mat)), -4.0
            )
            A[sb,end-nb+1:end] = b * (HxT[iLS] * iMx_b * Hx_b .+ HyT[iLS] * iMy_b  * Hy_b)
            # Boundary conditions for outer boundaries
            A[end-nb+1:end,sb] = b_b * (HxT_b * iMx_b' * Hx[iLS] .+ HyT_b * iMy_b' * Hy[iLS])
        end #ls_advection

        veci(rhs,grid,iLS+1) .= χ[iLS] * vec(a0) #vec(a0[iLS])

        

        # printstyled(color=:red, @sprintf "\n veci(rhs,grid,iLS+1) %.2i %.2e %.2e \n" iLS maximum(abs.(veci(rhs,grid,iLS+1))) maximum(abs.(BC.LS[iLS].val)))
        # print("\n a0 max  ", maximum(a0)," min ",minimum(a0))

    end #for iLS in 1:num.nLS

    vecb(rhs,grid) .= χ_b * vec(a0_b)

    printstyled(color=:red, @sprintf "\n vecb(rhs,grid) %.2e %.2e \n" maximum(abs.(vecb(rhs,grid))) maximum(abs.(BC.left.val)))

    # b_phi_ele = zeros(grid)
    # veci(rhs_scal,grid,1) .+= op.opC_pL.M * vec(b_phi_ele)

    if num.null_space == 0
        @time @inbounds @threads for i in 1:A.m
            @inbounds A[i,i] += 1e-10
        end
    end
    
    #region Solve Ax = b
    if num.electrical_potential_nonlinear_solver == 0

        #region Successive substitution

        if num.solver == 0
            @time ph.phi_eleD .= A \ rhs
        elseif num.solver == 1
            # using MUMPS, MPI, SparseArrays, LinearAlgebra
            MPI.Init()
            # A = sprand(10, 10, 0.2) + I
            # rhs = rand(10)
            # x = MUMPS.solve(A, rhs)
            # norm(x - A \ rhs) / norm(x)
            # ph.phi_eleD .= x
            ph.phi_eleD .= MUMPS.solve(A, rhs)
            MPI.Finalize()
        
        elseif num.solver == 2
            # diagA = A.diag

            d = collect(diag(A))
            for i in eachindex(d)
                if iszero(d[i])
                    print("\n diag i ", i,"\n d[i] ",d[i])
                end
                # print("\n diag i ", i,"\n d[i] ",d[i])
                d[i] = ifelse(iszero(d[i]), one(d[i]), 1/d[i])
                # d[i] = ifelse(iszero(d[i]), a*one(d[i]), zero(d[i]))
                
            end
            # A + Diagonal(d)

            # diagA = diag(A,0)
            invdiagA= Diagonal(d)
            newA = invdiagA * A
            newrhs=  invdiagA * rhs
            # newA= inv(diagA) * A 
            # newrhs=  inv(diagA) * rhs
            @time ph.phi_eleD .= newA \ newrhs

            d = collect(diag(newA))
            for i in eachindex(d)
                # if iszero(d[i])
                #     print("\n diag i ", i,"\n d[i] ",d[i])
                # end
                print("\n diag i ", i,"\n d[i] ",d[i])
                # d[i] = ifelse(iszero(d[i]), one(d[i]), 1/d[i])
                # d[i] = ifelse(iszero(d[i]), a*one(d[i]), zero(d[i]))
                
            end

        end

        print("\n rhs max  ", maximum(rhs)," min",minimum(rhs))

        print("\n norm 2 ", norm(A*ph.phi_eleD)," rhs ",norm(rhs))

        if norm(rhs) >0.0
            print("\n norm 2 ", norm(A*ph.phi_eleD -rhs)/norm(rhs))
        end

        print("\n rhs ", minimum(vecb_L(rhs,grid))," rhs ",maximum(vecb_L(rhs,grid)))



        #TODO reevaluate F=Ax-b 
        printstyled(color=:red, @sprintf "\n Residual" )


        rhs_updated = fnzeros(grid,num)

        for iLS in 1:num.nLS
            veci(rhs_updated,grid,iLS+1) .= χ[iLS] * vec(a0) #vec(a0[iLS])
            # printstyled(color=:red, @sprintf "\n veci(rhs,grid,iLS+1) %.2i %.2e %.2e \n" iLS maximum(abs.(veci(rhs,grid,iLS+1))) maximum(abs.(BC.LS[iLS].val)))
            # print("\n a0 max  ", maximum(a0)," min ",minimum(a0))
        end #for iLS in 1:num.nLS

        #TODO reevaluate BC
        update_electrical_current_from_Butler_Volmer!(num,grid,heat,ph.phi_eleD,i_butler)

        # print("\n i_butler ",i_butler)
        update_BC_electrical_potential!(num,grid,BC,elec_cond,coeffD,i_butler)

        a0_b = zeros(nb)
        _a1_b = zeros(nb)
        _b_b = zeros(nb)
        for iLS in 1:num.nLS
            set_borders!(grid, grid.LS[iLS].cl, grid.LS[iLS].u, a0_b, _a1_b, _b_b, BC, num.n_ext_cl)
        end

        vecb(rhs_updated,grid) .= χ_b * vec(a0_b)


        F_residual = A*ph.phi_eleD -rhs_updated

    
        

        if num.io_pdi>0
            iLSpdi = 1
            # dcap_1 for Wall capacity (left)
            # II = ind.b_left[1][i]
            # opC.χ_b[i, i] = geo.dcap[II,1]
            try
                # in YAML file: save only if iscal ==1 for example
                PDI_status = @ccall "libpdi".PDI_multi_expose("check_electrical_potential_convergence"::Cstring,
                "residual_1D"::Cstring, F_residual::Ptr{Cdouble}, PDI_OUT::Cint,
                "phi_ele_1D"::Cstring, ph.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
                "elec_cond_1D"::Cstring, coeffD::Ptr{Cdouble}, PDI_OUT::Cint, 
                "rhs_1D"::Cstring, rhs_updated::Ptr{Cdouble}, PDI_OUT::Cint, 
                "dcap_1"::Cstring, grid.LS[iLSpdi].geoL.dcap[:,:,1]::Ptr{Cdouble}, PDI_OUT::Cint,
                C_NULL::Ptr{Cvoid})::Cint
            catch error
                printstyled(color=:red, @sprintf "\n PDI error \n")
                print(error)
                printstyled(color=:red, @sprintf "\n PDI error \n")
            end
        end #if io_pdi


        print("\n rhs ", minimum(vecb_L(rhs_updated,grid))," rhs ",maximum(vecb_L(rhs_updated,grid)))
    
        #endregion Successive substitution

    elseif num.electrical_potential_nonlinear_solver == 1 #Newton-Raphson

        #region Newton-Raphson

        # print("\n Newton-Raphson")

        # print("\n rhs max  ", maximum(rhs)," min",minimum(rhs))

        # print("\n norm 2 ", norm(A*ph.phi_eleD)," rhs ",norm(rhs))

        # if norm(rhs) >0.0
        #     print("\n norm 2 ", norm(A*ph.phi_eleD -rhs)/norm(rhs))
        # end

        # print("\n rhs ", minimum(vecb_L(rhs,grid))," rhs ",maximum(vecb_L(rhs,grid)))



        #TODO reevaluate F=Ax-b 
        # printstyled(color=:red, @sprintf "\n Residual" )


        rhs_updated = fnzeros(grid,num)

        for iLS in 1:num.nLS
            veci(rhs_updated,grid,iLS+1) .= χ[iLS] * vec(a0) #vec(a0[iLS])
            printstyled(color=:red, @sprintf "\n veci(rhs,grid,iLS+1) %.2i %.2e %.2e \n" iLS maximum(abs.(veci(rhs,grid,iLS+1))) maximum(abs.(BC.LS[iLS].val)))
            print("\n a0 max  ", maximum(a0)," min ",minimum(a0))
        end #for iLS in 1:num.nLS

        #TODO reevaluate BC
        update_electrical_current_from_Butler_Volmer!(num,grid,heat,ph.phi_eleD,i_butler)

        print("\n i_butler ",i_butler)
        update_BC_electrical_potential!(num,grid,BC,elec_cond,coeffD,i_butler)

        a0_b = zeros(nb)
        _a1_b = zeros(nb)
        _b_b = zeros(nb)
        for iLS in 1:num.nLS
            set_borders!(grid, grid.LS[iLS].cl, grid.LS[iLS].u, a0_b, _a1_b, _b_b, BC, num.n_ext_cl)
        end

        vecb(rhs_updated,grid) .= χ_b * vec(a0_b)


        F_residual = A*ph.phi_eleD -rhs_updated

    
        #Wall capacity (left)
        # II = ind.b_left[1][i]
        # opC.χ_b[i, i] = geo.dcap[II,1]

        if num.io_pdi>0
            iLSpdi = 1
            try
                # in YAML file: save only if iscal ==1 for example
                PDI_status = @ccall "libpdi".PDI_multi_expose("check_electrical_potential_convergence"::Cstring,
                "residual_1D"::Cstring, F_residual::Ptr{Cdouble}, PDI_OUT::Cint,
                "phi_ele_1D"::Cstring, ph.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
                "elec_cond_1D"::Cstring, coeffD::Ptr{Cdouble}, PDI_OUT::Cint, 
                "rhs_1D"::Cstring, rhs_updated::Ptr{Cdouble}, PDI_OUT::Cint, 
                "dcap_1"::Cstring, grid.LS[iLSpdi].geoL.dcap[:,:,1]::Ptr{Cdouble}, PDI_OUT::Cint,
                C_NULL::Ptr{Cvoid})::Cint
            catch error
                printstyled(color=:red, @sprintf "\n PDI error \n")
                print(error)
                printstyled(color=:red, @sprintf "\n PDI error \n")
            end
        end #if io_pdi


        print("\n rhs ", minimum(vecb_L(rhs_updated,grid))," rhs ",maximum(vecb_L(rhs_updated,grid)))


        Jacobian = copy(A)

        #TODO reevaluate BC
        i_butler_derivative = zeros(grid.ny)
        jacobian_Butler = zeros(grid.ny)
        update_derivative_electrical_current_from_Butler_Volmer!(num,grid,heat,ph.phi_eleD,i_butler_derivative)

        # print("\n i_butler_derivative ",i_butler_derivative)

        update_BC_derivative_electrical_potential!(num,grid,jacobian_Butler,elec_cond,coeffD,i_butler_derivative)

        # derivative_Butler = - 1.0/elec_cond * fact * 2*np.cosh(fact * (-0.6 - U[0]))

        # BC Jacobian
        ny = grid.ny
        jacobian_bc = zeros(nb)
        jacobian_bc[1:ny] = jacobian_Butler
        # Jacobian[end-nb+1:end,end-nb+1:end] .+= χ_b * jacobian_bc

        # Jacobian[end-nb+1:end,end-nb+1:end] .-= Diagonal(jacobian_bc)

        # print("\n jacobian_bc ",jacobian_bc)
        # print("\n Jacobian[end-nb+1:end,end-nb+1:end]",Jacobian[end-nb+1,:])

        # Jacobian[end-nb+1:end,end-nb+1:end] .-= χ_b * vec(jacobian_bc)
        Jacobian[end-nb+1:end,end-nb+1:end] .-= χ_b * Diagonal(vec(jacobian_bc))


        # print("\n Jacobian[end-nb+1:end,end-nb+1:end]",Jacobian[end-nb+1,:])


        # Jacobian[end-nb+1:end,end-nb+1:end] .+= jacobian_Butler
        # Jacobian[end-nb+1:end-nb+ny,end-nb+1:end-nb+ny] .+= χ_b * #Diagonal(jacobian_Butler) #TODO dx ? diag ?


        phi_increment = Jacobian \ (-F_residual)

        # print("\n ph.phi_eleD border",vecb(ph.phi_eleD,grid))

        # print("\n phi_increment border",vecb(phi_increment,grid))

        # print("\n phi_increment bulk", minimum(veci(phi_increment,grid)), maximum(veci(phi_increment,grid)))


        ph.phi_eleD .+= phi_increment

        

        # print("\n ph.phi_eleD border",vecb(ph.phi_eleD,grid))



        # printstyled(color=:red, @sprintf "\n Residual" )

        #endregion Newton-Raphson

    end
    #endregion

    #TODO or use mul!(rhs_scal, BTL, phL.TD, 1.0, 1.0) like in :

    # kill_dead_cells!(phL.T, grid, grid.LS[1].geoL)
    # veci(phL.TD,grid,1) .= vec(phL.T)
    # rhs = set_heat!(
    #     BC_int[1], num, grid, op.opC_TL, grid.LS[1].geoL, phL, num.θd, BC_TL, grid.LS[1].MIXED, grid.LS[1].geoL.projection,
    #     ATL, BTL,
    #     op.opL, grid_u, grid_u.LS[1].geoL, grid_v, grid_v.LS[1].geoL,
    #     periodic_x, periodic_y, heat_convection, advection, BC_int
    # )
    # mul!(rhs, BTL, phL.TD, 1.0, 1.0)

    # phL.TD .= ATL \ rhs
    # phL.T .= reshape(veci(phL.TD,grid,1), grid)


    # phL.phi_eleD .= Ascal \ rhs_scal

    ph.phi_ele .= reshape(veci(ph.phi_eleD,grid,1), grid)


    if num.io_pdi>0
        try
            # printstyled(color=:magenta, @sprintf "\n PDI write_electrical_potential %.5i \n" num.current_i)
            #in YAML file: save only if iscal ==1 for example
            PDI_status = @ccall "libpdi".PDI_multi_expose("write_electrical_potential"::Cstring,
            # "iscal"::Cstring, iscal::Ref{Clonglong}, PDI_OUT::Cint,
            "rhs_1D"::Cstring, rhs::Ptr{Cdouble}, PDI_OUT::Cint,
            "phi_ele_1D"::Cstring, ph.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
            # "trans_scal_1DT"::Cstring, phL.trans_scalD'::Ptr{Cdouble}, PDI_OUT::Cint,
            C_NULL::Ptr{Cvoid})::Cint
        catch error
            printstyled(color=:red, @sprintf "\n PDI error \n")
            print(error)
            printstyled(color=:red, @sprintf "\n PDI error \n")
        end
    end #if io_pdi
    
end #set_poisson_variable_coeff

# function laplacian_bc_variable_coeff(opC, nLS, grid, coeffD)
#     @unpack BxT, ByT, Hx, Hy, iMx, iMy, Hx_b, Hy_b, iMx_b, iMy_b = opC

#     coeffD_borders = vecb(coeffD,grid)

#     mat_coeffDx_i = Diagonal(vec(coeffDu)) # coeffDu is a 2d matrix with shape (grid_u.ny, grid_u.nx), multiplies Hx
#     mat_coeffDy_i = Diagonal(vec(coeffDv)) # coeffDu is a 2d matrix with shape (grid_v.ny, grid_v.nx), multiplies Hy
#     mat_coeffDx_b = Diagonal(vec(coeffD_borders)) # coeffDu is a 1d vector with shape (2grid.ny + 2grid.nx), multiplies Hx_b and Hy_b


#     bc_L = []
#     for iLS in 1:nLS
#         push!(bc_L, BxT * iMx * mat_coeffDx_i *Hx[iLS] .+ ByT * iMy * mat_coeffDy_i *Hy[iLS])
#     end

#     bc_L_b = (BxT * iMx_b * mat_coeffDx_b *Hx_b .+ ByT * iMy_b * mat_coeffDx_b *Hy_b)

#     return bc_L, bc_L_b
# end

# function laplacian_variable_coeff(opC,coeffD)
#     @unpack Bx, By, BxT, ByT, iMx, iMy, tmp_x, tmp_y = opC

#     mat_coeffDx = Diagonal(vec(coeffDx_bulk)) # coeffDx_bulk is a 2d matrix with shape (grid_u.ny, grid_u.nx), multiplies Bx
#     mat_coeffDy = Diagonal(vec(coeffDy_bulk)) # coeffDx_bulk is a 2d matrix with shape (grid_v.ny, grid_v.nx), multiplies By

#     mul!(tmp_x, iMx, mat_coeffDx, Bx)
#     L = BxT * tmp_x
#     mul!(tmp_y, iMy, mat_coeffDy, By)
#     L = L .+ ByT * tmp_y

#     return L
# end
    


"""
    BC_LS_test!(num, cl, grid, A, B, rhs, BC)

Prints the contact angle, from BC_LS!
"""
function BC_LS_test!(grid, u, A, B, rhs, BC)
    @unpack x, y, nx, ny, dx, dy, ind = grid
    @unpack all_indices, b_left, b_bottom, b_right, b_top = ind
    @unpack left, bottom, right, top = BC

    π2 = π / 2.0

    boundaries_idx = [b_left[1], b_bottom[1], b_right[1], b_top[1]]

    left2 = vcat(all_indices[2,2], all_indices[2:end-1,2], all_indices[end-1,2])
    bottom2 = vcat(all_indices[2,2], all_indices[2,2:end-1], all_indices[2,end-1])
    right2 = vcat(all_indices[2,end-1], all_indices[2:end-1,end-1], all_indices[end-1,end-1])
    top2 = vcat(all_indices[end-1,2], all_indices[end-1,2:end-1], all_indices[end-1,end-1])
    boundaries2 = [left2, bottom2, right2, top2]

    left3 = vcat(all_indices[3,3], all_indices[2:end-1,3], all_indices[end-2,3])
    bottom3 = vcat(all_indices[3,3], all_indices[3,2:end-1], all_indices[3,end-2])
    right3 = vcat(all_indices[3,end-2], all_indices[2:end-1,end-2], all_indices[end-2,end-2])
    top3 = vcat(all_indices[end-2,3], all_indices[end-2,2:end-1], all_indices[end-2,end-2])
    boundaries3 = [left3, bottom3, right3, top3]

    boundaries_t = [left, bottom, right, top]

    direction = [y, x, y, x]

    for (i, (idx, idx2, idx3, xy)) in enumerate(zip(boundaries_idx, boundaries2, boundaries3, direction))
        pks, _ = findminima(abs.(u[idx]))
        # if is_neumann(boundaries_t[i])

        #     for (II, JJ) in zip(idx, idx2)
        #         pII = lexicographic(II, grid.ny)
        #         pJJ = lexicographic(JJ, grid.ny)

        #         A[pII,:] .= 0.0
        #         A[pII,pII] = 1.0
        #         A[pII,pJJ] = -1.0
        #         B[pII,:] .= 0.0
        #     end
        if is_neumann_cl(boundaries_t[i]) && maximum(u[idx]) > 0.0 && minimum(u[idx]) < 0.0 && length(pks) >= 2
            pks1 = idx[pks[1]]
            pkse = idx[pks[end]]

            # Gradually update the contact angle
            Δθe = 90.0 * π / 180

        

            # Find current contact angle
            dist = sqrt((x[idx2[pks[1]]] - x[pks1])^2 + (y[idx2[pks[1]]] - y[pks1])^2)
            old = u[pks1] - u[idx2[pks[1]]]
            # Levelset difference between two consecutive points might be bigger
            # than the distance between them if it's not reinitialized often enough
            if abs(old) > dist
                old = sign(old) * dist
            end
            θe_old = acos(old / dist)

            # printstyled(color=:green, @sprintf "\n θe_old : %.2e \n" θe_old)
            printstyled(color=:green, @sprintf "\n θe_old : %.2e °\n" θe_old*180.0/π)


            return x[idx2[pks]],y[idx2[pks]],x[pks1],y[pks1]

            # # Compute new contact angle
            # if abs(boundaries_t[i].θe - θe_old) > Δθe
            #     θe = θe_old + sign(boundaries_t[i].θe - θe_old) * Δθe
            # else
            #     θe = boundaries_t[i].θe
            # end

            # # distance between the center of the drop and the contact line
            # d = abs(xy[pks1] + u[pks1] - (xy[pkse] + u[pkse])) / 2.0

            # for (II, JJ) in zip(idx[2:end-1], idx2[2:end-1])
            #     pII = lexicographic(II, grid.ny)
            #     pJJ = lexicographic(JJ, grid.ny)

            #     A[pII,:] .= 0.0
            #     A[pII,pII] = 1.0
            #     A[pII,pJJ] = -1.0
            #     B[pII,:] .= 0.0

            #     # Compute levelset angle at a distance u[II] from the contact line
            #     if θe < π2
            #         newθ = atan(tan(θe) * (1.0 - u[II] / d))
            #     else
            #         newθ = π - atan(tan(π - θe) * (1.0 - u[II] / d))
            #     end

            #     rhs[pII] = dist * cos(newθ)
            # end
        # elseif is_neumann_cl(boundaries_t[i]) || is_neumann_inh(boundaries_t[i])
        #     for (II, JJ, KK) in zip(idx, idx2, idx3)
        #         pII = lexicographic(II, grid.ny)
        #         pJJ = lexicographic(JJ, grid.ny)

        #         A[pII,:] .= 0.0
        #         A[pII,pII] = 1.0
        #         A[pII,pJJ] = -1.0
        #         B[pII,:] .= 0.0

        #         rhs[pII] = u[JJ] - u[KK]
        #     end
        end
    end

    return nothing
end

function kill_dead_cells_val!(T::Vector, grid, geo,val)
    @unpack ny, ind = grid

    @inbounds @threads for II in ind.all_indices
        pII = lexicographic(II, ny)
     
        if geo.cap[II,5] < 1e-12
            # print("\n replace ",pII," ", geo.cap[II,5], " ",T[pII], " ",val)
            T[pII] = val
        end
    end
end

function kill_dead_cells_val!(S::SubArray{T,N,P,I,L}, grid, geo,val) where {T,N,P<:Vector{T},I,L}
    @unpack ny, ind = grid

    @inbounds @threads for II in ind.all_indices
        pII = lexicographic(II, ny)
        if geo.cap[II,5] < 1e-12
            # print("\n replace ",pII," ", geo.cap[II,5], " ",S[pII], " ",val)
            S[pII] = val
        end
    end
end

function kill_dead_cells_val!(S::SubArray{T,N,P,I,L}, grid, geo,val) where {T,N,P<:Array{T,3},I,L}
    @unpack ind = grid
    # print("kill dead cells mat")
    @inbounds @threads for II in ind.all_indices
        if geo.cap[II,5] < 1e-12
            # print(II, S[II])
            S[II] = val
            # print("v2",S[II])
        end
    end
end



function kill_dead_cells_val_wall!(S::SubArray{T,N,P,I,L}, grid, geo,val) where {T,N,P<:Vector{T},I,L}
    @unpack ny, ind = grid
    # In ind.b_right, the first vector is the boundary and the second one are the cells right next to them
    @inbounds @threads for II in ind.b_left[1]
        pII = lexicographic(II, ny)
        if geo.cap[II,5] < 1e-12
            print("\n replace ",pII," ", geo.cap[II,5], " ",S[pII], " ",val)
            S[pII] = val
        end
    end

    @inbounds @threads for II in ind.b_bottom[1]
        pII = lexicographic(II, ny)
        if geo.cap[II,5] < 1e-12
            print("\n replace ",pII," ", geo.cap[II,5], " ",S[pII], " ",val)
            S[pII] = val
        end
    end

    @inbounds @threads for II in ind.b_right[1]
        pII = lexicographic(II, ny)
        if geo.cap[II,5] < 1e-12
            print("\n replace ",pII," ", geo.cap[II,5], " ",S[pII], " ",val)
            S[pII] = val
        end
    end

    @inbounds @threads for II in ind.b_top[1]
        pII = lexicographic(II, ny)
        if geo.cap[II,5] < 1e-12
            print("\n replace ",pII," ", geo.cap[II,5], " ",S[pII], " ",val)
            S[pII] = val
        end
    end

end

# function kill_dead_bc_left_wall!(BC,grid,iLS,val)
#     for i = 1:grid.ny
#         II = CartesianIndex(i,1)
#         if grid.LS[iLS].geoL.cap[II,5] < 1e-12
#             BC.left.val[i] = val
#         end
#     end
# end


"""
    returns x,y,f,connectivities,vtx_index,num_seg which can be used to plot interfacial value
    supports on interface 
    field_index is the index of the interfacial value: ex 2 for LS 1 since field_index = 1 is the index of the bulk value
"""
function convert_interfacial_D_to_segments(num,gp,field_D,iLS,field_index)

    # x = zeros(0)
    # y = zeros(0)
    # f = zeros(0)
    x = Float64[]
    y = Float64[]
    f = Float64[]

    connectivities=Int64[]
    ij_index=-ones(Int64,(gp.ny, gp.nx))
    vtx_index = 0 #Integer(0) 

 


    for II in gp.LS[iLS].MIXED
        push!( x, gp.LS[1].mid_point[II].x * gp.dx[II] + gp.x[II] )
        push!( y, gp.LS[1].mid_point[II].y * gp.dy[II] + gp.y[II] )

        # push!( f, field[II] )

        # a[g.ny*g.nx*(p-1)+1:g.ny*g.nx*p]
        
        pII = lexicographic(II, gp.ny)

        push!( f, field_D[gp.ny*gp.nx*(field_index-1) + pII] )



        ij_index[II] = vtx_index

        # print("\n vtx ", vtx_index," ", II, " ",gp.x[II]-gp.dx[II]/2," ",gp.x[II]+gp.dx[II]/2," ",x," ",gp.y[II]-gp.dy[II]/2,
        # " ",gp.y[II]+gp.dy[II]/2,y," ",gp.x[II]," ",gp.y[II])


        #TODO store connvectivities or similar so that after: only need to filter and store non-null values ? 
        #no need to do the work n times for each scalar

        vtx_index +=1
    end

    num_seg = 0
    for j in 1:gp.ny-1
        for i in 1:gp.nx-1
            # II = CartesianIndex(i,j)
            if ij_index[j,i] >=0

                if ij_index[j+1,i]>=0  
                    push!(connectivities, ij_index[j,i])
                    push!(connectivities,ij_index[j+1,i])
                    num_seg+=1
                    # print("\n seg ", i," ",j," ",ij_index[j,i], " ", ij_index[j+1,i])
                end

                if ij_index[j,i+1]>=0  
                    push!(connectivities, ij_index[j,i])
                    push!(connectivities,ij_index[j,i+1])
                    num_seg+=1
                    # print("\n seg ", i," ",j," ",ij_index[j,i], " ", ij_index[j,i+1])
                end

            end

        end
    end

    # num_vtx = vtx_index already incremented vtx_index +=1
    
    return x,y,f,connectivities,vtx_index,num_seg

end





function contact_angle_advancing_receding(num,grid,grid_u, grid_v, iLS, II)

    # From code in run.jl

    # # TODO can be done witu u and v component ?
    # normalx = cos.(grid_u.LS[iLS].α)
    # normaly = sin.(grid_v.LS[iLS].α)

    # # grid_u.V .*= normalx
    # # grid_v.V .*= normaly

    # advancing_receding = normalx[II] * grid_u.V[II] + normaly[II] * grid_v.V[II] #TODO check dim , bounds


    #With interpolation

    normalx = cos.(grid.LS[iLS].α)
    normaly = sin.(grid.LS[iLS].α)

    cap1 = grid_u.LS[iLS].geoL.cap[II,5]
    cap3 = grid_u.LS[iLS].geoL.cap[δx⁺(II),5]
    tmpVx = (grid_u.V[II] * cap1 + grid_u.V[δx⁺(II)] * cap3) * inv_weight_eps(num,cap1 + cap3)
    
    cap2 = grid_v.LS[iLS].geoL.cap[II,5]
    cap4 = grid_v.LS[iLS].geoL.cap[δy⁺(II),5]
    tmpVy = (grid_v.V[II] * cap2 + grid_v.V[δy⁺(II)] * cap4) * inv_weight_eps(num,cap2 + cap4)

    advancing_receding = normalx[II] * tmpVx + normaly[II] * tmpVy

    return advancing_receding



    # V .= 0.0
    # cap1 = grid_u.LS[iLS].geoL.cap[II,5]
    # cap3 = grid_u.LS[iLS].geoL.cap[δx⁺(II),5]
    # tmpVx = (grid_u.V[II] * cap1 + grid_u.V[δx⁺(II)] * cap3) * inv_weight_eps(num,cap1 + cap3)
    
    # cap2 = grid_v.LS[iLS].geoL.cap[II,5]
    # cap4 = grid_v.LS[iLS].geoL.cap[δy⁺(II),5]
    # tmpVy = (grid_v.V[II] * cap2 + grid_v.V[δy⁺(II)] * cap4) * inv_weight_eps(num,cap2 + cap4)
    
    # tmpV = sqrt(tmpVx^2 + tmpVy^2)
    # β = atan(tmpVy, tmpVx)

    # printstyled(color=:magenta, @sprintf "\n beta β: %.2e : %.2e : %.2e" β grid.LS[iLS].α[II] )

    # if grid.LS[iLS].α[II] > 0.0 && β < 0.0
    #     β += 2π
    # end
    # if grid.LS[iLS].α[II] < 0.0 && β > 0.0
    #     β -= 2π
    # end

    # V[II] = tmpV * cos(β - grid.LS[iLS].α[II])

    
    
    
    # Project velocities to the normal and use advecion scheme for advection just
    # in the normal direction
    # tmpVx = zeros(grid)
    # tmpVy = zeros(grid)
    # V .= 0.0
    # @inbounds @threads for II in grid.LS[iLS].MIXED
    #     cap1 = grid_u.LS[iLS].geoL.cap[II,5]
    #     cap3 = grid_u.LS[iLS].geoL.cap[δx⁺(II),5]
    #     # tmpVx[II] = (grid_u.V[II] * cap1 + grid_u.V[δx⁺(II)] * cap3) / (cap1 + cap3 + eps(0.01))
    #     tmpVx[II] = (grid_u.V[II] * cap1 + grid_u.V[δx⁺(II)] * cap3) * inv_weight_eps(num,cap1 + cap3)
        

    #     cap2 = grid_v.LS[iLS].geoL.cap[II,5]
    #     cap4 = grid_v.LS[iLS].geoL.cap[δy⁺(II),5]
    #     # tmpVy[II] = (grid_v.V[II] * cap2 + grid_v.V[δy⁺(II)] * cap4) / (cap2 + cap4 + eps(0.01))
    #     tmpVy[II] = (grid_v.V[II] * cap2 + grid_v.V[δy⁺(II)] * cap4) * inv_weight_eps(num,cap2 + cap4)
        

    #     tmpV = sqrt(tmpVx[II]^2 + tmpVy[II]^2)
    #     β = atan(tmpVy[II], tmpVx[II])
    #     if grid.LS[iLS].α[II] > 0.0 && β < 0.0
    #         β += 2π
    #     end
    #     if grid.LS[iLS].α[II] < 0.0 && β > 0.0
    #         β -= 2π
    #     end

    #     V[II] = tmpV * cos(β - grid.LS[iLS].α[II])
    # end

    # i_ext, l_ext, b_ext, r_ext, t_ext = indices_extension(grid, grid.LS[iLS], grid.ind.inside, periodic_x, periodic_y)
    # field_extension!(grid, grid.LS[iLS].u, V, i_ext, l_ext, b_ext, r_ext, t_ext, num.NB, periodic_x, periodic_y)

    # rhs_LS .= 0.0
    # IIOE_normal!(grid, LS[iLS].A, LS[iLS].B, LS[iLS].u, V, CFL_sc, periodic_x, periodic_y)
    # BC_LS!(grid, LS[iLS].u, LS[iLS].A, LS[iLS].B, rhs_LS, BC_u)
    # BC_LS_interior!(num, grid, iLS, LS[iLS].A, LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
    # LS[iLS].u .= reshape(gmres(LS[iLS].A, LS[iLS].B * vec(LS[iLS].u) .+ rhs_LS), grid)

    # Impose contact angle if a wall is present
    # rhs_LS .= 0.0
    # LS[iLS].A.nzval .= 0.0
    # LS[iLS].B.nzval .= 0.0
    # for II in grid.ind.all_indices
    #     pII = lexicographic(II, grid.ny)
    #     LS[iLS].A[pII,pII] = 1.0
    #     LS[iLS].B[pII,pII] = 1.0
    # end
    # BC_LS_interior!(num, grid, iLS, LS[iLS].A, LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
    # LS[iLS].u .= reshape(gmres(LS[iLS].A, LS[iLS].B * vec(LS[iLS].u) .+ rhs_LS), grid)

end

function check_positivity_of_capacities(op,rhs_scal)

    min_cap = minimum(op.opC_uL.M)

    if min_cap<0.0
        printstyled(color=:red, @sprintf "\n [BUG] Sign capacity u\n" )
    end

    min_cap = minimum(op.opC_vL.M)

    if min_cap<0.0
        printstyled(color=:red, @sprintf "\n [BUG] Sign capacity v\n" )
    end

    # rhs_vec1 = vec1,rhs_scal

end


"""

"""
function compute_divergence!(num::Numerical{Float64, Int64},
    grid,
    # ::Mesh{Flower.GridCC, Float64, Int64},
    # grid_u::Mesh{Flower.GridFCx, Float64, Int64},
    # grid_v::Mesh{Flower.GridFCy, Float64, Int64},
    op::DiscreteOperators{Float64, Int64},
    A::SparseMatrixCSC{Float64, Int64},
    # rhs::Array{Float64, 1},
    # a0::Array{Float64, 2},
    tmp_vec_1D::Array{Float64, 1},
    tmp_vec_1D_2::Array{Float64, 1},
    Lv, 
    bc_Lv, 
    bc_Lv_b
    # vec_1D::Array{Float64, 1},
    # ls_advection::Bool
    )

    @unpack Bx, By, Hx, Hy, HxT, HyT, χ, M, iMx, iMy, Hx_b, Hy_b, HxT_b, HyT_b, iMx_b, iMy_b, iMx_bd, iMy_bd, χ_b = op.opC_vL
    @unpack BxT, ByT,tmp_x, tmp_y = op.opC_vL

    ni = grid.nx * grid.ny
    nb = 2 * grid.nx + 2 * grid.ny

    #TODO reset zero
    # rhs .= 0.0
    tmp_vec_1D .= 0.0
    tmp_vec_1D_2 .= 0.0
    A .= 0.0
    # a0 .= 0.0

    x_centroid = grid.x .+ getproperty.(grid.LS[1].geoL.centroid, :x) .* grid.dx
    y_centroid = grid.y .+ getproperty.(grid.LS[1].geoL.centroid, :y) .* grid.dy

    print("\n size ",size(x_centroid))
    print("\n size ",size(y_centroid))

    print("\n size op.opC_vL.M",size(op.opC_vL.M))


    print("\n x_centroid ",x_centroid[1,:])

    veci(tmp_vec_1D,grid,1)  .= vec(x_centroid)
    vecb_L(tmp_vec_1D,grid)  .= grid.x[:,1]
    vecb_R(tmp_vec_1D,grid)  .= grid.x[:,end]

    veci(tmp_vec_1D_2,grid,1) .= vec(y_centroid)
    vecb_B(tmp_vec_1D_2,grid) .= grid.y[1,:]
    vecb_T(tmp_vec_1D_2,grid) .= grid.y[end,:]

    for iLS in 1:num.nLS
        x_bc = grid.x .+ getproperty.(grid.LS[iLS].mid_point, :x) .* grid.dx
        y_bc = grid.y .+ getproperty.(grid.LS[iLS].mid_point, :y) .* grid.dy
        veci(tmp_vec_1D,grid,iLS+1) .= vec(x_bc)
        veci(tmp_vec_1D_2,grid,iLS+1) .= vec(y_bc)
    end


    # PDI_status = @ccall "libpdi".PDI_multi_expose("mesh"::Cstring, 
    # # "grad_x_coord"::Cstring, tmp_vec_u::Ptr{Cdouble}, PDI_OUT::Cint,   
    # # "grad_y_coord"::Cstring, tmp_vec_v::Ptr{Cdouble}, PDI_OUT::Cint,   
    # "mesh_x_1D"::Cstring, tmp_vec_1D::Ptr{Cdouble}, PDI_OUT::Cint,   
    # "mesh_y_1D"::Cstring, tmp_vec_1D_2::Ptr{Cdouble}, PDI_OUT::Cint,      
    # C_NULL::Ptr{Cvoid})::Cint

    #TODO sparse

    # compute_grad_T_x_array!(num.nLS, grid, grid_u, op.opC_pL, tmp_vec_u, tmp_vec_1D)
    # compute_grad_T_y_array!(num.nLS, grid, grid_v, op.opC_pL, tmp_vec_v, tmp_vec_1D_2)

    tmp_vec_1D .^= 2
    tmp_vec_1D_2 .^= 2



    # mul!(tmp_x, iMx, Bx)
    # L = BxT * tmp_x
    # mul!(tmp_y, iMy, By)
    # L = L .+ ByT * tmp_y

    # #Boundary for Laplacian
    # bc_L_b = (BxT * iMx_b * Hx_b .+ ByT * iMy_b * Hy_b)


    # if ls_advection
        # Poisson equation
    A[1:ni,1:ni] = Lv #pad(Lv, 4.0)
    A[1:ni,end-nb+1:end] = bc_Lv_b

    # # Boundary conditions for outer boundaries
    # A[end-nb+1:end,1:ni] = (HxT_b * iMx_b' * Bx .+ HyT_b * iMy_b' * By) #b_b * 
    # A[end-nb+1:end,end-nb+1:end] = pad((HxT_b * iMx_bd  * Hx_b .+ HyT_b * iMy_bd * Hy_b), -4.0) #b_b * 
    
    # end

    for iLS in 1:num.nLS

        sb = iLS*ni+1:(iLS+1)*ni


        # Poisson equation
        #Boundary for Laplacian from iLS
        A[1:ni,sb] = bc_Lv[iLS] #BxT * iMx * Hx[iLS] .+ ByT * iMy * Hy[iLS]

        # # Boundary conditions for inner boundaries
        # A[sb,1:ni] = (HxT[iLS] * iMx  * Bx .+ HyT[iLS] * iMy * By) #or vec1  b *
        # Contribution to Neumann BC from other boundaries
        # for i in 1:num.nLS
        #     if i != iLS
        #         A[sb,i*ni+1:(i+1)*ni] = (HxT[iLS] * iMx  * Hx[i] .+ HyT[iLS] * iMy * Hy[i]) #b * 
        #     end
        # end
        # A[sb,sb] = pad(
        #         (HxT[iLS] * iMx * Hx[iLS] .+ HyT[iLS] * iMy * Hy[iLS]), -4.0 #b *
        # )
        # A[sb,end-nb+1:end] = (HxT[iLS] * iMx_b * Hx_b .+ HyT[iLS] * iMy_b  * Hy_b) #b *
        # # Boundary conditions for outer boundaries
        # A[end-nb+1:end,sb] = (HxT_b * iMx_b' * Hx[iLS] .+ HyT_b * iMy_b' * Hy[iLS]) #b_b *

    end #for iLS in 1:num.nLS

    # if num.null_space == 0
    #     @time @inbounds @threads for i in 1:A.m
    #         @inbounds A[i,i] += 1e-10
    #     end
    # end

    print("\n size tmp_vec_1D ",size(tmp_vec_1D))

    print("\n size tmp_vec_1D ",size(A))

    
   rhs = A*tmp_vec_1D/2.0 #2V

   print("\n min",minimum(rhs)," max ",maximum(rhs))

   rhs_vec1 = vec1(vec(rhs),grid)

   print("\n rhs_vec1 before end",rhs_vec1[1,1])

   print("\n size op ",size(op.opC_vL.M))

   print("\n rhs_vec1 ",rhs_vec1[1,1])
   
   print("\n op ",op.opC_vL.M[1,1])

   print("\n ratio ",rhs_vec1[1,1]/op.opC_vL.M[1,1])
   
   print("\n size op ",size(op.opC_vL.M))
   print("\n size op diag ",size(op.opC_vL.M.diag))

   print("\n size rhs_vec1 ",size(rhs_vec1))


   op.opC_vL.M.diag .= rhs_vec1

   check_positivity_of_capacities(op,rhs_vec1)

    
end     #divergence


"""
Based on num.bulk_conductivity:
* 0 conductivity computed from wall concentration
* 1 conductivity computed from bulk concentration
* 2 conductivity computed from wall concentration and bulk concentration
* 3 conductivity computed from wall concentration and bulk concentration
"""
function update_BC_electrical_potential!(num,grid,BC_phi_ele,elec_cond,elec_condD,i_butler)

    # BC LS 2 in set_poisson directly
    
    #use Butler-Volmer, supposing the interfacial potential is acceptable and phi = phi_ele1 in metal 
    # for conductivity, use interfacial value or bulk in corresponding cell

    # TODO -(-i/kappa) in Flower ? so i_butler not -i_butler
    # For small cells
    if num.bulk_conductivity == 0
        BC_phi_ele.left.val .= i_butler./vecb_L(elec_condD, grid)

    elseif num.bulk_conductivity == 1
        # Recommended as long as cell merging not implemented:
        # Due to small cells, we may have slivers/small cells at the left wall, then the divergence term is small,
        # which produces higher concentration in front of the contact line
        BC_phi_ele.left.val .= i_butler./elec_cond[:,1]


    elseif num.bulk_conductivity == 2 || num.bulk_conductivity == 3
        BC_phi_ele.left.val .= i_butler./vecb_L(elec_condD, grid)

        iLS = 1 #TODO end ? if several grid.LS ?
        for j in 1:grid.ny
            II = CartesianIndex(j,1)
            if grid.LS[iLS].geoL.cap[II,5] < num.ϵ
                BC_phi_ele.left.val[j] = i_butler[j]/elec_cond[j,1] 
            end
        end
    end

end

"""
If the heat equation is solved, it uses the computed temperature, otherwise num.temperature0 is used.
"""
function update_electrical_current_from_Butler_Volmer!(num,grid,heat,phi_eleD,i_butler;T=nothing)

    #TODO dev multiple levelsets
    if heat
        i_butler .= butler_volmer_no_concentration.(num.alpha_a,num.alpha_c,num.Faraday,num.i0,vecb_L(phi_eleD, grid),
        num.phi_ele1,num.Ru,T)
    else
        if num.nLS == 1
            i_butler .= butler_volmer_no_concentration.(num.alpha_a,num.alpha_c,num.Faraday,num.i0,vecb_L(phi_eleD, grid),
            num.phi_ele1,num.Ru,num.temperature0)
        # else
            #imposed by LS 2
            # iLS_elec = 2
            # i_butler = butler_volmer_no_concentration.(num.alpha_a,num.alpha_c,num.Faraday,num.i0,veci(phL.phi_eleD, grid,iLS_elec+1),
            # num.phi_ele1,num.Ru,num.temperature0)
        end
    end   

end



"""
Based on num.bulk_conductivity:
* 0 conductivity computed from wall concentration
* 1 conductivity computed from bulk concentration
* 2 conductivity computed from wall concentration and bulk concentration
* 3 conductivity computed from wall concentration and bulk concentration
"""
function update_BC_derivative_electrical_potential!(num,grid,jacobian_Butler,elec_cond,elec_condD,i_butler_derivative)

    # BC LS 2 in set_poisson directly
    
    #use Butler-Volmer, supposing the interfacial potential is acceptable and phi = phi_ele1 in metal 
    # for conductivity, use interfacial value or bulk in corresponding cell

    # TODO -(-i/kappa) in Flower ? so i_butler not -i_butler
    # For small cells
    if num.bulk_conductivity == 0
        jacobian_Butler .= i_butler_derivative ./ vecb_L(elec_condD, grid)

    elseif num.bulk_conductivity == 1
        # Recommended as long as cell merging not implemented:
        # Due to small cells, we may have slivers/small cells at the left wall, then the divergence term is small,
        # which produces higher concentration in front of the contact line
        jacobian_Butler .= i_butler_derivative./elec_cond[:,1]


    elseif num.bulk_conductivity == 2 || num.bulk_conductivity == 3
        jacobian_Butler .= i_butler_derivative./vecb_L(elec_condD, grid)

        iLS = 1 #TODO end ? if several grid.LS ?
        for j in 1:grid.ny
            II = CartesianIndex(j,1)
            if grid.LS[iLS].geoL.cap[II,5] < num.ϵ
                jacobian_Butler[j] = i_butler_derivative[j]/elec_cond[j,1] 
            end
        end
    end

end

"""
If the heat equation is solved, it uses the computed temperature, otherwise num.temperature0 is used.
"""
function update_derivative_electrical_current_from_Butler_Volmer!(num,grid,heat,phi_eleD,i_butler_derivative;T=nothing)

    #TODO dev multiple levelsets
    if heat
        i_butler_derivative .= derivative_butler_volmer_no_concentration.(num.alpha_a,num.alpha_c,num.Faraday,num.i0,vecb_L(phi_eleD, grid),
        num.phi_ele1,num.Ru,T)
    else
        if num.nLS == 1
            i_butler_derivative .= derivative_butler_volmer_no_concentration.(num.alpha_a,num.alpha_c,num.Faraday,num.i0,vecb_L(phi_eleD, grid),
            num.phi_ele1,num.Ru,num.temperature0)
        # else
            #imposed by LS 2
            # iLS_elec = 2
            # i_butler = butler_volmer_no_concentration.(num.alpha_a,num.alpha_c,num.Faraday,num.i0,veci(phL.phi_eleD, grid,iLS_elec+1),
            # num.phi_ele1,num.Ru,num.temperature0)
        end
    end   

end


"""
Interpolate conductivity at center of control volumes for potential gradient at the border
"""
function interpolate_scalar_to_staggered_u_v_grids_at_border!(num,grid,coeffD,coeffDu,coeffDv)

    coeffDu .= 0.0
    coeffDv .= 0.0

    # or use lexicographic

    @inbounds @threads for II in grid.ind.b_left[1] #[2:end-1]
    # @inbounds @threads for II in grid_u.ind.b_left[1] #[2:end-1]
        coeffDu[II] = (vecb_L(coeffD,grid)[II[1]]+reshape(veci(coeffD,grid),grid)[II])/2.0 # veci(coeff,grid) or elec_cond
        # print("\n left, II ",II, coeffDu[II])
    end

    @inbounds @threads for II in grid.ind.b_right[1] #[2:end-1]
        coeffDu[ δx⁺(II) ] = (vecb_R(coeffD,grid)[II[1]]+reshape(veci(coeffD,grid),grid)[II])/2.0
        # print("\n right, II ",II, coeffDu[II])

    end

    @inbounds @threads for II in grid.ind.b_bottom[1] #[2:end-1]
        coeffDv[II] = (vecb_B(coeffD,grid)[II[2]]+reshape(veci(coeffD,grid),grid)[II])/2.0
    end

    @inbounds @threads for II in grid.ind.b_top[1] #[2:end-1]
        coeffDv[δy⁺(II)] = (vecb_T(coeffD,grid)[II[2]]+reshape(veci(coeffD,grid),grid)[II])/2.0
    end

    if num.io_pdi>0
        try
            PDI_status = @ccall "libpdi".PDI_multi_expose("write_electrical_conductivity"::Cstring,
            "elec_cond_1D"::Cstring, coeffD::Ptr{Cdouble}, PDI_OUT::Cint,
            "conductivity_u"::Cstring, coeffDu::Ptr{Cdouble}, PDI_OUT::Cint,
            "conductivity_v"::Cstring, coeffDv::Ptr{Cdouble}, PDI_OUT::Cint,   
            C_NULL::Ptr{Cvoid})::Cint
        catch error
            printstyled(color=:red, @sprintf "\n PDI error \n")
            print(error)
            printstyled(color=:red, @sprintf "\n PDI error \n")
        end
    end #if io_pdi

    return coeffDu, coeffDv

end


"""
Interpolate conductivity at center of control volumes for potential gradient at the border
For test and prints
"""
function interpolate_scalar_to_staggered_u_v_grids_at_border_test!(num,grid,coeffD,coeffDu,coeffDv)

    coeffDu .= 0.0
    coeffDv .= 0.0

    # or use lexicographic

    @inbounds @threads for II in grid.ind.b_left[1] #[2:end-1]
    # @inbounds @threads for II in grid_u.ind.b_left[1] #[2:end-1]
        coeffDu[II] = (vecb_L(coeffD,grid)[II[1]]+reshape(veci(coeffD,grid),grid)[II])/2.0 # veci(coeff,grid) or elec_cond
        print("\n left, II ",II, coeffDu[II])
    end

    @inbounds @threads for II in grid.ind.b_right[1] #[2:end-1]
        coeffDu[ δx⁺(II) ] = (vecb_R(coeffD,grid)[II[1]]+reshape(veci(coeffD,grid),grid)[II])/2.0
        print("\n right, II ",II, coeffDu[II])

    end

    @inbounds @threads for II in grid.ind.b_bottom[1] #[2:end-1]
        coeffDv[II] = (vecb_B(coeffD,grid)[II[2]]+reshape(veci(coeffD,grid),grid)[II])/2.0
    end

    @inbounds @threads for II in grid.ind.b_top[1] #[2:end-1]
        coeffDv[δy⁺(II)] = (vecb_T(coeffD,grid)[II[2]]+reshape(veci(coeffD,grid),grid)[II])/2.0
    end

    if num.io_pdi>0
        try
            PDI_status = @ccall "libpdi".PDI_multi_expose("write_electrical_conductivity"::Cstring,
            "elec_cond_1D"::Cstring, coeffD::Ptr{Cdouble}, PDI_OUT::Cint,
            "conductivity_u"::Cstring, coeffDu::Ptr{Cdouble}, PDI_OUT::Cint,
            "conductivity_v"::Cstring, coeffDv::Ptr{Cdouble}, PDI_OUT::Cint,   
            C_NULL::Ptr{Cvoid})::Cint
        catch error
            printstyled(color=:red, @sprintf "\n PDI error \n")
            print(error)
            printstyled(color=:red, @sprintf "\n PDI error \n")
        end
    end #if io_pdi

    return coeffDu, coeffDv

end


"""
Main loop for the resolution of the Poisson equation
"""
function solve_poisson_loop!(num::Numerical{Float64, Int64},
                            grid::Mesh{Flower.GridCC, Float64, Int64},
                            grid_u::Mesh{Flower.GridFCx, Float64, Int64},
                            grid_v::Mesh{Flower.GridFCy, Float64, Int64},
                            op::DiscreteOperators{Float64, Int64},
                            Ascal::SparseMatrixCSC{Float64, Int64},
                            rhs_scal::Array{Float64, 1},
                            tmp_vec_p::Array{Float64, 2},
                            tmp_vec_p0::Array{Float64, 2},
                            tmp_vec_p1::Array{Float64, 2},
                            a1_p::SparseMatrixCSC{Float64, Int64},
                            BC_phi_ele::BoundariesInt,
                            phL::Phase{Float64},
                            phS::Phase{Float64},
                            elec_cond::Array{Float64, 2},
                            elec_condD::Array{Float64, 1},
                            tmp_vec_u::Array{Float64, 2},
                            tmp_vec_v::Array{Float64, 2},
                            i_butler::Array{Float64, 1},
                            ls_advection::Bool,
                            heat::Bool)
            
    # Electroneutrality assumption : rhs in bulk: 0
    a0_p = [] 
    for i in 1:num.nLS
        push!(a0_p, zeros(grid))
    end

    #region Update conductivity
    update_electrical_conductivity!(num,grid,elec_cond,elec_condD)
    #endregion Update conductivity


    #Update Butler-Volmer Boundary Condition with new potential 
    if occursin("Butler",num.electrolysis_reaction) && num.nLS == 1

        printstyled(color=:red, @sprintf "\n Recomputing Butler \n" )

        #region Update current
        if num.electrolysis_reaction == "Butler_no_concentration"                
            update_electrical_current_from_Butler_Volmer!(num,grid,heat,phL.phi_eleD,i_butler;phL.T)
        end
        #endregion Update current

        update_BC_electrical_potential!(num,grid,BC_phi_ele,elec_cond,elec_condD,i_butler)


        # if heat
        #     BC_phi_ele.left.val = -butler_volmer_no_concentration.(num.alpha_a,num.alpha_c,num.Faraday,num.i0,phL.phi_ele[:,1],num.phi_ele1,num.Ru,phL.T)./elec_cond[:,1]
        # else
        #     BC_phi_ele.left.val = -butler_volmer_no_concentration.(num.alpha_a,num.alpha_c,num.Faraday,num.i0,phL.phi_ele[:,1],num.phi_ele1,num.Ru,num.temperature0)./elec_cond[:,1]
            
        #     # for iscal=1:num.nb_transported_scalars
        #     #     BC_trans_scal[iscal].left.val = butler_volmer_no_concentration.(num.alpha_a,num.alpha_c,num.Faraday,num.i0,phL.phi_ele[:,1],num.phi_ele1,num.Ru,num.temperature0)./(2*num.Faraday*num.diffusion_coeff[iscal])
        #     #     if iscal==1 || iscal==2
        #     #         BC_trans_scal[iscal].left.val .*=-1 #H2O
        #     #     end
        #     # end
        # end    

    # elseif num.electrolysis_reaction == ""
    #     # BC_phi_ele.left.val = -butler_volmer_concentration.(num.alpha_a,num.alpha_c,num.Faraday,num.i0,phL.phi_ele[:,1],num.phi_ele1,num.Ru,num.temperature0)./elec_cond
    
        

        # TODO 
        #Remove Nan when dividing by conductivity which may be null

        # TODO bug 1                           

        for iLS in 1:num.nLS
            # kill_dead_bc_left_wall!(vecb(elec_condD,grid), grid, iLS,1.0)
            for i = 1:grid.ny
                # print("vecb cap",vecb_L(grid.LS[iLS].geoL.cap[:,5],grid))
                
                # II = CartesianIndex(i,1)
                # II = grid.ind.b_left[1][i]
                # opC.χ_b[i, i] = geo.dcap[II,1]
                # TODO not cleat why zero: grid.LS[iLS].geoL.cap[II,1]
                #TODO cf update LS convection not convection where something is overwritten
                # wall_liquid_height = grid.LS[iLS].geoL.cap[II,1]
                wall_liquid_height = op.opC_pL.χ_b[i, i]
                if wall_liquid_height < 1e-12
                    BC_phi_ele.left.val[i] = 1.0
                    print("\n bug BC_phi_ele.left.val[i] ",II," ",grid.LS[iLS].geoL.cap[II,:])
                    # print("\n opC.χ_b[i, i] ",op.opC_pL.χ_b[i, i])
                end
            end
        end



    end #if occursin("Butler",num.electrolysis_reaction)


    #TODO nLS
    #TODO kill_dead_cells! ?
    kill_dead_cells!(phL.phi_ele, grid, grid.LS[1].geoL)
    veci(phL.phi_eleD,grid,1) .= vec(phL.phi_ele)
    

    # Store current potential (iteration k)
    phi_eleD_0 = copy(phL.phi_eleD)


    # iterate (non-linear BC with Butler) 
    for poisson_iter=1:num.electrical_potential_max_iter

        # printstyled(color=:orange, @sprintf "\n poisson iter %.2i \n" poisson_iter)

        compute_grad_phi_ele!(num, grid, grid_u, grid_v, grid_u.LS[end], grid_v.LS[end], phL, phS, op.opC_pL, op.opC_pS, 
        elec_cond,tmp_vec_u,tmp_vec_v,tmp_vec_p,tmp_vec_p0,tmp_vec_p1) #TODO current

    
        residual_electrical_potential = maximum(abs.(-tmp_vec_p[div(grid.ny,2),:].+butler_volmer_no_concentration.(num.alpha_a,num.alpha_c,num.Faraday,num.i0,vecb_L(phL.phi_eleD, grid),
                    num.phi_ele1,num.Ru,num.temperature0)))

        # Absolute variation 
        variation_electrical_potential = maximum(abs.(phi_eleD_0 - phL.phi_eleD))           

        @ccall "libpdi".PDI_multi_expose("check_electrical_potential"::Cstring,
        "poisson_iter"::Cstring, poisson_iter ::Ref{Clonglong}, PDI_OUT::Cint,
        "i_current_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
        "i_current_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,  
        "i_current_mag"::Cstring, tmp_vec_p1::Ptr{Cdouble}, PDI_OUT::Cint,
        "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
        "elec_cond_1D"::Cstring, elec_condD::Ptr{Cdouble}, PDI_OUT::Cint,  
        "BC_phi_ele_left"::Cstring, BC_phi_ele.left.val::Ptr{Cdouble}, PDI_OUT::Cint,  
        "levelset_p"::Cstring, grid.LS[num.index_levelset_pdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
        # "levelset_p"::Cstring, grid.LS[1].u::Ptr{Cdouble}, PDI_OUT::Cint,
        # "levelset_p"::Cstring, grid.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
        "residual_electrical_potential"::Cstring, residual_electrical_potential ::Ref{Cdouble}, PDI_OUT::Cint,
        "variation_electrical_potential"::Cstring, variation_electrical_potential ::Ref{Cdouble}, PDI_OUT::Cint,
        # "grad_phi_ele_u"::Cstring, tmp_vec_u::Ptr{Cdouble}, PDI_OUT::Cint,  
        C_NULL::Ptr{Cvoid})::Cint

        electrical_potential_converged = (
            (residual_electrical_potential  < num.electrical_potential_relative_residual) &&
            (variation_electrical_potential < num.electrical_potential_residual)
        )

        if electrical_potential_converged
            # printstyled(color=:orange, @sprintf "\n End Poisson loop \n")
            break
        end

        # Store current potential (iteration k)
        phi_eleD_0 = copy(phL.phi_eleD)
        
        # printstyled(color=:orange, @sprintf "\n grad poisson iter %.2i \n" poisson_iter)

        # print("\n grad ", tmp_vec_u[div(grid_u.ny,2),:]," \n")

        # @ccall "libpdi".PDI_multi_expose("solve_poisson"::Cstring,
        # # "i_current_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
        # # "i_current_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,  
        # # "i_current_mag"::Cstring, phL.i_current_mag::Ptr{Cdouble}, PDI_OUT::Cint,
        # "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
        # "elec_cond_1D"::Cstring, elec_condD::Ptr{Cdouble}, PDI_OUT::Cint,  
        # "BC_phi_ele_left"::Cstring, BC_phi_ele.left.val::Ptr{Cdouble}, PDI_OUT::Cint,  
        # # "grad_phi_ele_u"::Cstring, tmp_vec_u::Ptr{Cdouble}, PDI_OUT::Cint,  
        # C_NULL::Ptr{Cvoid})::Cint

    
        if num.electrolysis_reaction == "Butler_no_concentration"

            # if num.poisson_newton ==1
            #     vecb_L(phL.phi_eleD, grid) = vecb_L(phL.phi_eleD, grid) - (partial...+)/deriv
            # end


            #TODO dev multiple levelsets
            if heat
                i_butler = butler_volmer_no_concentration.(num.alpha_a,num.alpha_c,num.Faraday,num.i0,vecb_L(phL.phi_eleD, grid),
                num.phi_ele1,num.Ru,phL.T)
            else
                if num.nLS == 1
                    i_butler = butler_volmer_no_concentration.(num.alpha_a,num.alpha_c,num.Faraday,num.i0,vecb_L(phL.phi_eleD, grid),
                    num.phi_ele1,num.Ru,num.temperature0)
                # else
                    #imposed by LS 2
                    # iLS_elec = 2
                    # i_butler = butler_volmer_no_concentration.(num.alpha_a,num.alpha_c,num.Faraday,num.i0,veci(phL.phi_eleD, grid,iLS_elec+1),
                    # num.phi_ele1,num.Ru,num.temperature0)
                end
                    
            end   

            if poisson_iter>1
                if num.bulk_conductivity == 0
                    BC_phi_ele.left.val .= i_butler./vecb_L(elec_condD, grid)

                elseif num.bulk_conductivity == 1
                    # Recommended as long as cell merging not implemented:
                    # Due to small cells, we may have slivers/small cells at the left wall, then the divergence term is small,
                    # which produces higher concentration in front of the contact line
                    BC_phi_ele.left.val .= i_butler./elec_cond[:,1]


                elseif num.bulk_conductivity == 2 || num.bulk_conductivity == 3
                    BC_phi_ele.left.val .= i_butler./vecb_L(elec_condD, grid)

                    iLS = 1 #TODO end ? if several grid.LS ?
                    for j in 1:grid.ny
                        II = CartesianIndex(j,1)
                        if grid.LS[iLS].geoL.cap[II,5] < num.ϵ
                            BC_phi_ele.left.val[j] = i_butler[j]/elec_cond[j,1] 
                        end
                    end
                    
                # if num.bulk_conductivity == 3
                #     elec_condD .= compute_ele_cond.(num.Faraday,num.diffusion_coeff[num.index_electrolyte],num.Ru, num.temperature0, num.concentration0[num.index_electrolyte])
                #     elec_cond .= reshape(vec1(elec_condD,grid),grid)
                # end
                
                end
            end 
            # print("\n BC_phi_ele",BC_phi_ele,"\n")

            
        end

        print("\n BC_phi_ele ",BC_phi_ele)

        solve_poisson_variable_coeff!(num, 
        grid, 
        grid_u, 
        grid_v, 
        op.opC_pL,
        Ascal, 
        rhs_scal,
        tmp_vec_p, #a0
        a1_p,
        BC_phi_ele,
        phL,    
        elec_cond,                    
        elec_condD,
        tmp_vec_u,
        tmp_vec_v,
        # tmp_vec_u0,
        # tmp_vec_v0,
        i_butler,
        ls_advection,
        heat)


        @ccall "libpdi".PDI_multi_expose("solve_poisson"::Cstring,
        # "i_current_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
        # "i_current_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,  
        # "i_current_mag"::Cstring, phL.i_current_mag::Ptr{Cdouble}, PDI_OUT::Cint,
        "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
        "elec_cond_1D"::Cstring, elec_condD::Ptr{Cdouble}, PDI_OUT::Cint,
        "BC_phi_ele_left"::Cstring, BC_phi_ele.left.val::Ptr{Cdouble}, PDI_OUT::Cint,  
        C_NULL::Ptr{Cvoid})::Cint

        #TODO or linearize 

        #TODO compute grad

        if num.electrical_potential>0
            compute_grad_phi_ele!(num, grid, grid_u, grid_v, grid_u.LS[end], grid_v.LS[end], phL, phS, op.opC_pL, op.opC_pS, 
            elec_cond,tmp_vec_u,tmp_vec_v,tmp_vec_p,tmp_vec_p0,tmp_vec_p1) #TODO current
            
            # printstyled(color=:orange, @sprintf "\n grad poisson iter %.2i \n" poisson_iter)

            # print("\n grad ", tmp_vec_u[div(grid_u.ny,2),:]," \n")
            
            # print("\n grad ", tmp_vec_u[div(grid_u.ny,2),1]," \n")
            # print("\n BC_phi_ele ", BC_phi_ele.left.val[div(grid_u.ny,2)]," \n")
            # print("\n i_butler ", i_butler," \n")

        end

    end #for loop Poisson

    # printstyled(color=:cyan, @sprintf "\n after solve_poisson_variable_coeff! \n")
    # print_electrolysis_statistics(num,grid,phL)

    PDI_status = @ccall "libpdi".PDI_multi_expose("print_variables"::Cstring,
        "nstep"::Cstring, num.current_i ::Ref{Clonglong}, PDI_OUT::Cint,
        "time"::Cstring, num.time::Ref{Cdouble}, PDI_OUT::Cint,
        "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
        "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
        "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
        "levelset_p"::Cstring, grid.LS[num.index_levelset_pdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
        "levelset_u"::Cstring, grid_u.LS[num.index_levelset_pdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
        "levelset_v"::Cstring, grid_v.LS[num.index_levelset_pdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
        # "levelset_p_wall"::Cstring, LStable::Ptr{Cdouble}, PDI_OUT::Cint,
        "trans_scal_1DT"::Cstring, phL.trans_scalD'::Ptr{Cdouble}, PDI_OUT::Cint,
        "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
        # "i_current_x"::Cstring, Eus::Ptr{Cdouble}, PDI_OUT::Cint,   
        # "i_current_y"::Cstring, Evs::Ptr{Cdouble}, PDI_OUT::Cint,   
        # "velocity_x"::Cstring, us::Ptr{Cdouble}, PDI_OUT::Cint,   
        # "velocity_y"::Cstring, vs::Ptr{Cdouble}, PDI_OUT::Cint,      
        # "radius"::Cstring, current_radius::Ref{Cdouble}, PDI_OUT::Cint,  
        # "intfc_vtx_num"::Cstring, intfc_vtx_num::Ref{Clonglong}, PDI_OUT::Cint, 
        # "intfc_seg_num"::Cstring, intfc_seg_num::Ref{Clonglong}, PDI_OUT::Cint, 
        # "intfc_vtx_x"::Cstring, intfc_vtx_x::Ptr{Cdouble}, PDI_OUT::Cint,
        # "intfc_vtx_y"::Cstring, intfc_vtx_y::Ptr{Cdouble}, PDI_OUT::Cint,
        # "intfc_vtx_field"::Cstring, intfc_vtx_field::Ptr{Cdouble}, PDI_OUT::Cint,
        # "intfc_vtx_connectivities"::Cstring, intfc_vtx_connectivities::Ptr{Clonglong}, PDI_OUT::Cint,
        C_NULL::Ptr{Cvoid})::Cint

    if any(isnan, phL.phi_eleD)
        print("\n phL.uD: ",any(isnan, phL.uD) , "\n phL.vD: ",any(isnan, phL.vD) , "\n phL.TD: ",any(isnan, phL.TD) , "\n phS.uD: ",any(isnan, phS.uD) , "\n phS.vD: ",any(isnan, phS.vD) , "\n phS.TD: ",any(isnan, phS.TD) ,
        "\n phL.trans_scalD: ",any(isnan, phL.trans_scalD) , "\n phL.phi_eleD: ",any(isnan, phL.phi_eleD) ,
        "\n phL.u: ",norm(phL.u) > 1e8 , "\n phS.u: ",norm(phS.u) > 1e8 , "\n phL.T: ",norm(phL.T) > 1e8 , "\n phS.T: ",norm(phS.T) > 1e8 , "\n phL.trans_scal: ",norm(phL.trans_scal) > 1e8 , "\n phL.phi_ele: ",norm(phL.phi_ele) > 1e8)

        print("\n phL.phi_eleD: ",any(isnan, phL.phi_eleD),"\n phL.phi_ele: ",any(isnan, phL.phi_ele),"\n")

        print("\n Ascal: ",any(isnan, Ascal),"\n rhs_scal: ",any(isnan, rhs_scal),"\n")

        print("\n \n vecb_L",vecb_L(phL.phi_eleD[:,1], grid))

    end


    # TODO compute magnitude of exchange current
    # gradient!(::Neumann, Ox, Oy, Bx, By, HNx, HNy, Divx, Divy, dcap, num.n, BC, all_indices, b_left_u, b_bottom_v, b_right_u, b_top_v, b_left_p, b_bottom_p, b_right_p, b_top_p)
    # TODO add post-treatment variables

    #TODO update BC concentration

    
    #region conductivity
    update_electrical_conductivity!(num,grid,elec_cond,elec_condD)
    #endregion conductivity

    if num.electrical_potential>0
        compute_grad_phi_ele!(num, grid, grid_u, grid_v, grid_u.LS[end], grid_v.LS[end], phL, phS, op.opC_pL, op.opC_pS, 
        elec_cond,tmp_vec_u,tmp_vec_v,tmp_vec_p,tmp_vec_p0,tmp_vec_p1) #TODO current
    end

    # scal_magnitude

    # phL.i_current_mag .*= elec_cond # i=-κ∇ϕ here magnitude

    printstyled(color=:green, @sprintf "\n test grad")

    # compute_grad_p!(num,grid, grid_u, grid_v, phL.phi_eleD, op.opC_pL, op.opC_uL, op.opC_vL)


    # #store in us, vs instead of Eus, Evs
    # interpolate_grid_liquid!(grid,grid_u,grid_v,phL.Eu, phL.Ev,tmp_vec_p,tmp_vec_p0)

    # @ccall "libpdi".PDI_multi_expose("write_data_elec"::Cstring,
    # "i_current_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
    # "i_current_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,  
    # "i_current_mag"::Cstring, phL.i_current_mag::Ptr{Cdouble}, PDI_OUT::Cint,
    # "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
    # C_NULL::Ptr{Cvoid})::Cint

end


"""
update electrical conductivity, using temperature array if it is solved, or homogeneous temperature, 
depending on concentration (solved or homogeeneous concentration)
"""
function update_electrical_conductivity!(num,grid,elec_cond,elec_condD)
    # Constant electrical conductivity assumption
    #TODO electrical conductivity depends on concentration
    #iKOH index of KOH 
    # kappa_ele=2*num.Faraday^2*num.concentration0[iKOH]*num.diffusion_coeff[iKOH]/(num.Ru*T)
    # elec_cond=2*num.Faraday^2*trans_scal[iKOH]*num.diffusion_coeff[iKOH]/(num.Ru*T)

    # #print(@sprintf "TODO elec cond and boundary conditions need to be updated for potential\n")

    if num.nb_transported_scalars>1
        if heat 
            elec_condD .= compute_ele_cond.(num.Faraday,num.diffusion_coeff[num.index_electrolyte],num.Ru, phL.TD, phL.trans_scalD[:,num.index_electrolyte])
            elec_cond .= reshape(vec1(elec_condD,grid),grid)
            # elec_cond .= compute_ele_cond.(num.Faraday,num.diffusion_coeff[num.index_electrolyte],num.Ru, phL.T, phL.trans_scal)
            # elec_cond = 2*num.Faraday^2 .*phL.trans_scal[:,:,2].*num.diffusion_coeff[2]./(num.Ru.*phL.T) #phL.T
        else
            elec_condD .= compute_ele_cond.(num.Faraday,num.diffusion_coeff[num.index_electrolyte],num.Ru, num.temperature0, phL.trans_scalD[:,num.index_electrolyte])
            elec_cond .= reshape(vec1(elec_condD,grid),grid)
            # elec_cond .= compute_ele_cond.(num.Faraday,num.diffusion_coeff[num.index_electrolyte],num.Ru, num.temperature0, phL.trans_scal)
            # elec_cond = 2*num.Faraday^2 .*phL.trans_scal[:,:,2].*num.diffusion_coeff[2]./(num.Ru*num.temperature0) 
            
            if num.bulk_conductivity == 3 #test homogeneous conductiviy
                elec_condD .= compute_ele_cond.(num.Faraday,num.diffusion_coeff[num.index_electrolyte],num.Ru, num.temperature0, num.concentration0[num.index_electrolyte])
                elec_cond .= reshape(vec1(elec_condD,grid),grid)
            end
        
        end
    else
        elec_condD .= compute_ele_cond.(num.Faraday,num.diffusion_coeff[num.index_electrolyte],num.Ru, num.temperature0, num.concentration0[num.index_electrolyte])
        elec_cond .= reshape(vec1(elec_condD,grid),grid)
    end

    # printstyled(color=:red, @sprintf "\n test conductivity")
    # # elec_condD .= 2*num.Faraday^2 .*num.concentration0[2].*num.diffusion_coeff[2]./(num.Ru.*num.temperature0)
    # test_filter_concentration!(num,grid,phL.trans_scalD[:,2],num.concentration0[2])

    # elec_condD = 2*num.Faraday^2 .*phL.trans_scalD[:,2].*num.diffusion_coeff[2]./(num.Ru.*num.temperature0)

    # # TODO icurrent mag and replace huge val scal
    
    # print("\n test i ",-2*num.Faraday^2*num.concentration0[2]*num.diffusion_coeff[2]./(num.Ru.*num.temperature0))*(num.phi0-num.phi1)/(1e-4)
    # print("\n test i ",-2*(num.Faraday^2)*num.concentration0[2]*num.diffusion_coeff[2]./(num.Ru.*num.temperature0))*(num.phi0-num.phi1)/(1e-4)
    # printstyled(color=:red, @sprintf "\n test conductivity")

    # if num.nb_transported_scalars>1
    #     if heat 
    #         elec_condD .= compute_ele_cond.(num.Faraday,num.diffusion_coeff[num.index_electrolyte],num.Ru, phL.TD, phL.trans_scalD[:,num.index_electrolyte])
    #         elec_cond .= reshape(vec1(elec_condD,grid),grid)
    #         # elec_cond .= compute_ele_cond.(num.Faraday,num.diffusion_coeff[num.index_electrolyte],num.Ru, phL.T, phL.trans_scal)
    #         # elec_cond = 2*num.Faraday^2 .*phL.trans_scal[:,:,2].*num.diffusion_coeff[2]./(num.Ru.*phL.T) #phL.T
    #     else
    #         elec_condD .= compute_ele_cond.(num.Faraday,num.diffusion_coeff[num.index_electrolyte],num.Ru, num.temperature0, phL.trans_scalD[:,num.index_electrolyte])
    #         elec_cond .= reshape(vec1(elec_condD,grid),grid)
    #         # elec_cond .= compute_ele_cond.(num.Faraday,num.diffusion_coeff[num.index_electrolyte],num.Ru, num.temperature0, phL.trans_scal)
    #         # elec_cond = 2*num.Faraday^2 .*phL.trans_scal[:,:,2].*num.diffusion_coeff[2]./(num.Ru*num.temperature0) 
        
    #         if num.bulk_conductivity == 3
    #             elec_condD .= compute_ele_cond.(num.Faraday,num.diffusion_coeff[num.index_electrolyte],num.Ru, num.temperature0, num.concentration0[num.index_electrolyte])
    #             elec_cond .= reshape(vec1(elec_condD,grid),grid)
    #         end
    #     end
    # end

    # # if electrolysis && num.nb_transported_scalars>1
    # #     if heat 
    # #         elec_cond = 2*num.Faraday^2 .*phL.trans_scal[:,:,2].*num.diffusion_coeff[2]./(num.Ru.*phL.T) #phL.T
    # #     else
    # #         elec_cond = 2*num.Faraday^2 .*phL.trans_scal[:,:,2].*num.diffusion_coeff[2]./(num.Ru*num.temperature0) 
    # #     end
    # # else 
    # #     elec_cond = ones(grid)
    # #     printstyled(color=:green, @sprintf "\n conductivity one")

    # # end 

end