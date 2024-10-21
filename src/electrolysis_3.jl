
# """
# print_sizes

# Debug sizes of expression
# """
# macro print_sizes(e::Expr)
#     print
# 	return esc( :( ones(size($e)[1]) .+ $e) )
# end


# function print_sizes(prog)
#     ex1 = Meta.parse(prog)
#     for i=1:size(eval(ex1.args))[1]
#         print("\n i ",i," ",ex1.args[i]," ",typeof(ex1.args[i]))

#         try
#             print(size(eval(ex1.args[i])))

#         catch e
#             # print("\n i ",i," ",ex1.args[i]," ",typeof(ex1.args[i]))

#         end

#         if typeof(ex1.args[i]) isa Symbol
#             print("\n i ",i," ",ex1.args[i]," ",typeof(ex1.args[i]))
#         else
#             print("\n i ",i," ",ex1.args[i]," ",typeof(ex1.args[i]))
#             print(size(eval(ex1.args[i])))
#         end
#         # try
#         #     print(eval(ex1.args[i]))
#         #     print(size(eval(ex1.args[i])))
#         # catch e
#         # end
#     end
# end

# function print_sizes(prog)
#     ex1 = Meta.parse(prog)
#     for i=1:size(eval(ex1.args))[1]
#         print("\n i ",i," ",ex1.args[i]," ",typeof(ex1.args[i])," ",(typeof(ex1.args[i]) <: Symbol))
       
#         if (typeof(ex1.args[i]) <: Symbol)
#             print("\n symbol")
#         else
#             print(eval(ex1.args[i]))
#             print(typeof(eval(ex1.args[i])))
#             print(size(eval(ex1.args[i])))
#         end
#         # try
#         #     print(test_macro(ex1.args[i]))
#         # catch e
#         # end

#         # try
#         #     print(test_macro2(ex1.args[i]))
#         # catch e
#         # end
       
#         # try
#         #     print(eval(ex1.args[i]))
#         # catch e
#         # end
#         # try
#         #     print(size(eval(ex1.args[i])))
#         # catch e
#         # end
#     end
# end

# macro test_macro(e::Expr)
#     return esc( :( size($e) ) )
# 	# return esc(size($e))
# end
# macro test_macro2(e::Expr)
#     return esc( :( size(eval(eval(:e))) ) )
#     # return esc( :( size($$e) ) )
# 	# return esc(size($e))
# end

"""    
set_poisson_variable_coeff!

# Arguments
- bc_type: BC for interface, num, grid, 
- a0, 
- opC, 
- pC_v,
- A, 
- BC: BC for wall
- ls_advection

# divg(kappa mathbf{q} )=kappa (divg(mathbf{q}))+ mathbf{q}⋅grad(kappa)

"""
function set_poisson_variable_coeff_no_interpolation!(num::Numerical{Float64, Int64},
    grid::Mesh{Flower.GridCC, Float64, Int64},
    grid_u::Mesh{Flower.GridFCx, Float64, Int64},
    grid_v::Mesh{Flower.GridFCy, Float64, Int64},
    # op::Operators{Float64, Int64},
    opC::Operators{Float64, Int64},
    A::SparseMatrixCSC{Float64, Int64},
    rhs::Array{Float64, 1},
    a0::Array{Float64, 2},
    a1::SparseMatrixCSC{Float64, Int64},
    BC::BoundariesInt,
    ph::Phase{Float64},
    coeffD,
    ls_advection::Bool)
    @unpack Bx, By, Hx, Hy, HxT, HyT, χ, M, iMx, iMy, Hx_b, Hy_b, HxT_b, HyT_b, iMx_b, iMy_b, iMx_bd, iMy_bd, χ_b = opC

    @unpack BxT, ByT,tmp_x, tmp_y = opC


    ni = grid.nx * grid.ny
    nb = 2 * grid.nx + 2 * grid.ny

    #TODO reset zero
    rhs .= 0.0

    A .= 0.0

    a0_b = zeros(nb)
    _a1_b = zeros(nb)
    _b_b = zeros(nb)
    for iLS in 1:num.nLS
        set_borders!(grid, grid.LS[iLS].cl, grid.LS[iLS].u, a0_b, _a1_b, _b_b, BC, num.n_ext_cl)
    end
    a1_b = Diagonal(vec(_a1_b))
    b_b = Diagonal(vec(_b_b))

   
    #Laplacian
    mul!(tmp_x, iMx, Bx)
    L = BxT * tmp_x
    mul!(tmp_y, iMy, By)
    L = L .+ ByT * tmp_y

    
    # "
    # cf
    # ∇ϕ_x = opC_p.iMx * opC_p.Bx * vec1(phi_eleD,grid) .+ opC_p.iMx_b * opC_p.Hx_b * vecb(phi_eleD,grid)
    # ∇ϕ_y = opC_p.iMy * opC_p.By * vec1(phi_eleD,grid) .+ opC_p.iMy_b * opC_p.Hy_b * vecb(phi_eleD,grid)

    # for iLS in 1:nLS
    #     ∇ϕ_x .+= opC_p.iMx * opC_p.Hx[iLS] * veci(phi_eleD,grid,iLS+1)
    #     ∇ϕ_y .+= opC_p.iMy * opC_p.Hy[iLS] * veci(phi_eleD,grid,iLS+1)
    # end
    # "

 

    print("\n Sizes")

    # print(size(L),"\n")
    print(size(vec1(coeffD,grid)),"\n")
    # print(size(HxT_b),"\n")
    # print(size(iMx_b'),"\n")
    # print(size(vec1(coeffD,grid)),"\n")
    # print(size(Bx),"\n")
    # print(size(HyT_b),"\n")
    # print(size(iMy_b'),"\n")
    # print(size(By),"\n")


    # L = vec1(coeffD,grid) * L
    
    # print(L)
    # print(size(Diagonal(vec1(coeffD,grid))),"\n")


    L = Diagonal(vec1(coeffD,grid)) * L 

    print("\n L",L[1,:])
    print("\n L",L[:,1])

    
    # L = BxT * iMx * Bx
    # L size nx*ny
    # tmp_x ((nx+1)*ny ?, nx*ny)
    # BxT (nx*ny, (nx+1)*ny ?)
    # Bx ((nx+1)*ny ?, nx*ny)
    # iMx ((nx+1)*ny ?, (nx+1)*ny ?)

    Ltmp = vec1(coeffD,grid) * Bx' * iMx' * iMx * Bx

    print("\n Ltmp ",size(Ltmp))



    print("\n Sizes")

    print("\n L ",size(L))
    print("\n tmp_x ",size(tmp_x))
    print("\n BxT ",size(BxT))
    print("\n Bx ",size(Bx))
    print("\n opC.iMx  ",size(opC.iMx ))
    print("\n opC.iMy  ",size(opC.iMy ))
    print("\n tmp_y  ",size(tmp_y))
    print("\n vec1(coeffD,grid) ",size(vec1(coeffD,grid)))
    print("\n By ",size(By))



    print("\n grad x",size(opC.iMx * opC.Bx * vec1(coeffD,grid)))

    # print(size((opC.iMx * opC.Bx * vec1(coeffD,grid) .* tmp_x)),"\n")

    # print(size(dot(opC.iMx * opC.Bx * vec1(coeffD,grid) , tmp_x)),"\n")

    # print_sizes("opC.iMx * opC.Bx * vec1(coeffD,grid)")
    # print_sizes(" opC.BxT * opC.iMx' *opC.iMx * opC.Bx * vec1(coeffD,grid)")
    # print_sizes("dot(opC.iMx * opC.Bx * vec1(coeffD,grid) , tmp_x) + dot(opC.iMy * opC.By * vec1(coeffD,grid) , tmp_y)")

    # mul!(tmp_x, iMx, Bx)

    # phi -->(∇ϕ_x,∇ϕ_y)
    # dot product ∇ϕ_x,∇ϕ_y -->sum

    # L +=  (opC.iMx * opC.Bx * vec1(coeffD,grid) .* tmp_x) + (opC.iMy * opC.By .* vec1(coeffD,grid) * tmp_y)

    L +=  dot(opC.iMx * opC.Bx * vec1(coeffD,grid) , tmp_x) + dot(opC.iMy * opC.By * vec1(coeffD,grid) , tmp_y)

    #Boundary for Laplacian
    bc_L_b = (BxT * iMx_b * Hx_b .+ ByT * iMy_b *Hy_b)

    bc_L_b = vecb(coeffD,grid) * bc_L_b + opC.iMx_b * opC.Hx_b * vecb(coeffD,grid) * iMx_b * Hx_b .+ opC.iMy_b * opC.Hy_b * vecb(coeffD,grid) * iMy_b *Hy_b

    #TODO still in dev
      
    if ls_advection
        # Poisson equation
        A[1:ni,1:ni] = pad(L, -4.0)
        A[1:ni,end-nb+1:end] = bc_L_b

        # Boundary conditions for outer boundaries
        A[end-nb+1:end,1:ni] = vec1(coeffD,grid) * (-b_b) * (HxT_b * iMx_b' * Bx .+ HyT_b * iMy_b' * By)
        A[end-nb+1:end,end-nb+1:end] = -pad(vecb(coeffD,grid) * b_b * (HxT_b * iMx_bd * Hx_b .+ HyT_b * iMy_bd * Hy_b) .- χ_b * a1_b, 4.0)
    end

    for iLS in 1:num.nLS
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

                #TODO BC LS 2 in set_poisson directly
                #  BC_phi_ele.LS[2].val .= i_butler./vecb_L(elec_condD, grid)

                # #TODO BC LS 2
                # BC_phi_ele.LS[2].val .= i_butler./elec_cond[:,1]

                #use Butler-Volmer, supposing the interfacial potential is acceptable and phi = phi_ele1 in metal 
                # for conductivity, use interfacial value or bulk in corresponding cell

                # TODO -(-i/kappa) in Flower ? so i_butler not -i_butler
                # For small cells

                
                # iLS_elec = 2
                # or use iLS and...
                # if 0 Neumann do not need complicated
                #how to write efficiently BC

               


                if num.bulk_conductivity == 0
                    
                    # for II in grid.LS[iLS].MIXED

                        butler_volmer_no_concentration_potential_Neumann!.(num,
                        reshape(veci(ph.phi_eleD, grid,iLS+1),grid),
                        reshape(veci(ph.trans_scalD[:,2],grid,iLS+1),grid),
                        num.temperature0,
                        a0[II])
                    # end

                elseif num.bulk_conductivity == 1
                    # # Recommended as long as cell merging not implemented:
                    # # Due to small cells, we may have slivers/small cells at the left wall, then the divergence term is small,
                    # # which produces higher concentration in front of the contact line
                    # a0 .= i_butler./elec_cond[:,1]
                    @error ("error elseif num.bulk_conductivity == 1")

                elseif num.bulk_conductivity == 2
                    #TODO remove reshape and use a mapping 
                    butler_volmer_no_concentration_potential_Neumann!.(num,
                    reshape(veci(ph.phi_eleD, grid,iLS+1),grid),
                    reshape(veci(ph.trans_scalD[:,2],grid,iLS+1),grid),
                    num.temperature0,
                    a0)

                    for II in grid.LS[iLS].MIXED
                        if grid.LS[iLS].geoL.cap[II,5] < num.ϵ #volume
                            #use bulk conductivity of mixed cell
                            butler_volmer_no_concentration_potential_Neumann!.(num,
                            reshape(veci(ph.phi_eleD, grid,iLS+1),grid),
                            ph.trans_scal[II,2],
                            num.temperature0,
                            a0[II]) #TODO if temperature solved temperature[II]

                            # a0[II]  = butler_volmer_no_concentration.(num.alpha_a,num.alpha_c,num.Faraday,num.i0,veci(ph.phi_eleD, grid,iLS+1),
                            # num.phi_ele1,num.Ru,num.temperature0)./ph.elec_cond[II]
                        end

                        if grid.LS[iLS].geoL.cap[II,1] < num.ϵ #here length so volume approx num.ϵ^2
                            a0[II] = 1.0
                        end

                        # TODO
                        #Remove Nan when dividing by conductivity which may be null
                        # kill_dead_bc_left_wall!(vecb(elec_condD,grid), grid, iLS,1.0)

                    end

     

                end #bulk_conductivity

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
            interpolate_scalar!(grid, grid_u, grid_v, reshape(veci(coeffD,grid,iLS+1), grid), coeffDu, coeffDv)
            mat_coeffDx_i = Diagonal(vec(coeffDu)) # coeffDu is a 2d matrix with shape (grid_u.ny, grid_u.nx), multiplies Hx
            mat_coeffDy_i = Diagonal(vec(coeffDv)) # coeffDu is a 2d matrix with shape (grid_v.ny, grid_v.nx), multiplies Hy

            
            # Poisson equation
            #Boundary for Laplacian from iLS
            A[1:ni,sb] = veci(coeffD,grid,iLS+1) * (BxT * iMx * Hx[iLS] .+ ByT * iMy * Hy[iLS])

            # Boundary conditions for inner boundaries
            A[sb,1:ni] = veci(coeffD,grid,iLS+1) * (-b) * (HxT[iLS] * iMx * Bx .+ HyT[iLS] * iMy * By) #or vec1
            # Contribution to Neumann BC from other boundaries
            for i in 1:num.nLS
                if i != iLS
                    A[sb,i*ni+1:(i+1)*ni] = veci(coeffD,grid,iLS+1) * (-b) * (HxT[iLS] * iMx * Hx[i] .+ HyT[iLS] * iMy * Hy[i])
                end
            end
            A[sb,sb] = -pad(
                veci(coeffD,grid,iLS+1) * b * (HxT[iLS] * iMx * Hx[iLS] .+ HyT[iLS] * iMy * Hy[iLS]) .- χ[iLS] * a1 .+
                a2 * Diagonal(diag(fs_mat)), 4.0
            )
            A[sb,end-nb+1:end] = vecb(coeffD,grid) * b * (HxT[iLS] * iMx_b * Hx_b .+ HyT[iLS] * iMy_b * Hy_b)
            # Boundary conditions for outer boundaries
            A[end-nb+1:end,sb] = veci(coeffD,grid,iLS+1) * (-b_b) * (HxT_b * iMx_b' * Hx[iLS] .+ HyT_b * iMy_b' * Hy[iLS])
        end

        veci(rhs,grid,iLS+1) .= -χ[iLS] * vec(a0) #vec(a0[iLS])

        #TODO print variable resolution
        printstyled(color=:red, @sprintf "\n veci(rhs,grid,iLS+1) %.2i %.2e %.2e \n" iLS maximum(abs.(veci(rhs,grid,iLS+1))) maximum(abs.(BC.LS[iLS].val)))

    
    end #for iLS in 1:num.nLS

    vecb(rhs,grid) .= -χ_b * vec(a0_b)

    printstyled(color=:red, @sprintf "\n vecb(rhs,grid) %.2e %.2e \n" maximum(abs.(vecb(rhs,grid))) maximum(abs.(BC.left.val)))

    # b_phi_ele = zeros(grid)
    # veci(rhs_scal,grid,1) .+= op.opC_pL.M * vec(b_phi_ele)

    if num.null_space == 0
        @time @inbounds @threads for i in 1:A.m
            @inbounds A[i,i] += 1e-10
        end
    end
    
    @time ph.phi_eleD .= A \ rhs

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
    
end