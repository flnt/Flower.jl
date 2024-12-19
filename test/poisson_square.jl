# using Revise
using Flower

using Polynomials

#From poisson.jl (Quiros)

# Careful here for SPD -

function f(x, y)
    # return cos(π^2 * x * y)*sin(π^2 * x * y)
    # return -10 * (x-2.5) #L=2.5 from 0 to L
    return -10 * (1.25-x) #L=2.5 from -L/2 to L/2
end

function ∇fx(x, y)
    # return π^2 * y * (cos(π^2 * x * y)^2 - sin(π^2 * x * y)^2)
    return 10
end

function minusf(x, y)
    # return cos(π^2 * x * y)*sin(π^2 * x * y)
    # return -10 * (x-2.5) #L=2.5 from 0 to L
    return 10 * (1.25-x) #L=2.5 from -L/2 to L/2
end


function minus∇fx(x, y)
    return -10
end

function ∇fy(x, y)
    # return π^2 * x * (cos(π^2 * x * y)^2 - sin(π^2 * x * y)^2)
    return 0.0
end

function Δf(x, y)
    # return -4 * (π^2)^2 *sin(π^2 * x*y)*cos(π^2 * x*y)*(x^2+y^2)
    return 0.0
end

mutable struct conv_regression
    yreg::Array{Float64, 1}
    coef0::Float64
    coef1::Float64
end

function regression(x, y, x_reg)
    coeffs = fit(log.(x), log.(y), 1).coeffs
    y_reg = exp.(coeffs[2].*log.(x_reg) .+ coeffs[1])
    result = conv_regression(y_reg, -coeffs[1], -coeffs[2])

    return result
end

function dirichlet_bcs!(gp, D)
    @unpack x, y, dx, dy, LS, ind = gp

    @inbounds @threads for II in ind.inside
        x_bc = LS[1].mid_point[II].x * dx[II] + x[II]
        y_bc = LS[1].mid_point[II].y * dy[II] + y[II]
        D[II] = f(x_bc, y_bc)
    end
end

function neumann_bcs!(gp, N)
    @unpack x, y, dx, dy, LS, ind = gp

    @inbounds @threads for II in ind.inside
        x_bc = LS[1].mid_point[II].x * dx[II] + x[II]
        y_bc = LS[1].mid_point[II].y * dy[II] + y[II]
        Nx = ∇fx(x_bc, y_bc)
        Ny = ∇fy(x_bc, y_bc)

        N[II] = Nx * cos(LS[1].α[II]+π) + Ny * sin(LS[1].α[II]+π)
    end

    return nothing
end

function robin_bcs!(gp, R)
    @unpack x, y, dx, dy, LS, ind = gp

    @inbounds @threads for II in ind.inside
        x_bc = LS[1].mid_point[II].x * dx[II] + x[II]
        y_bc = LS[1].mid_point[II].y * dy[II] + y[II]
        Nx = ∇fx(x_bc, y_bc)
        Ny = ∇fy(x_bc, y_bc)

        R[II] = Nx * cos(LS[1].α[II]+π) + Ny * sin(LS[1].α[II]+π) + f(x_bc, y_bc)
    end

    

    return nothing
end

#tolerance for tests
test_tolerance = 1.e-14

test_tolerance_solution_absolute = 1.e-10 #e-12


#verbosity
verbosity = 0 #for more details, set to 1
verbosity = 1 #for more details, set to 1


ϵ = 0.001
bc = rob

bc = dir

if is_dirichlet(bc)
    bc_str = "dir"
elseif is_neumann(bc)
    bc_str = "neu"
elseif is_robin(bc)
    bc_str = "rob"
end

L0 = 2.5

# n_cases = 4 

# st_case = 6 #starting at n=64
# st_case = 3 #starting at n=8

n_cases = 3
st_case = 5 #starting at n=32


npts = 2 .^ [st_case:(st_case+n_cases-1)...]

npts = [32]


print("\n number of points ", npts, "\n")



l1 = zeros(n_cases)
l2 = zeros(n_cases)
loo = zeros(n_cases)
l1_mixed = zeros(n_cases)
l2_mixed = zeros(n_cases)
loo_mixed = zeros(n_cases)
l1_full = zeros(n_cases)
l2_full = zeros(n_cases)
loo_full = zeros(n_cases)


# Convergence study loop
for (i,n) in enumerate(npts)


    x = collect(LinRange(-L0 / 2.0, L0 / 2.0, n + 1))
    y = collect(LinRange(-L0 / 2.0, L0 / 2.0, n + 1))

    # x_old = LinRange(-L0 / 2.0, L0 / 2.0, n + 1)
    # y_old = LinRange(-L0 / 2.0, L0 / 2.0, n + 1)

    # print("\n x ",x," \n")
    # print("\n x ",collect(x)," \n")

    # print("\n x ",x_old," \n")
    # print("\n x ",collect(x_old)," \n")



    num = Numerical(
        case = "Cylinder",
        x = x,
        y = y,
        CFL = 1.0,
        max_iterations = 0,
        R = 1.0,
        ϵ = ϵ,
        verbosity = verbosity,
    )

    gp, gu, gv = init_meshes(num)
    op, phS, phL = init_fields(num, gp, gu, gv)
    # gp.LS[1].u .*= -1.0 #solve inside circle by setting sign of levelset

    # gp.LS[1].u .= 0.0 #deactivate interface
    gp.LS[1].u .= 1.0 #deactivate interface

    # Signs in divergence theorem
    sign_left = -1.0 #n \cdot e_x = -1
    sign_right = 1.0 #n \cdot e_x = 1
    sign_bottom = -1.0 #n \cdot e_y = -1
    sign_top = 1.0 #n \cdot e_y = 1

    sign_left = 1.0
    sign_right = 1.0
    sign_bottom = 1.0
    sign_top = 1.0

    sign_Poisson = -1.0 #because we solve -laplacian p = f

    run_forward!(num, gp, gu, gv, op, phS, phL)
    BC = Boundaries(left = Neumann(val = sign_Poisson * 10.0 * sign_left), # -lapl p = f for SPD so -10
    bottom = Neumann(val = 0.0 * sign_bottom),
    right = Dirichlet(val = 0.0 * sign_right),
    top = Neumann(val = 0.0 * sign_top)) #border boundaries

    BC_int = [bc]

    θd = zeros(gp)
    if is_dirichlet(bc)
        dirichlet_bcs!(gp, θd)
    elseif is_neumann(bc)
        neumann_bcs!(gp, θd)
    elseif is_robin(bc)
        robin_bcs!(gp, θd)
    end

    ni = gp.nx * gp.ny
    nb = 2 * gp.nx + 2 * gp.ny
    nt = (num.nLS + 1) * ni + nb

    if num.verbosity > 0
        printstyled(color=:green, @sprintf "\n ni %.3i nb %.3i nt %.3i \n" ni nb nt)
    end

    A = spzeros(nt, nt)

    periodic_x = false
    periodic_y = false

    update_all_ls_data(num, gp, gu, gv, BC_int, periodic_x, periodic_y, false)


    x_centroid = gp.x .+ getproperty.(gp.LS[1].geoL.centroid, :x) .* gp.dx
    y_centroid = gp.y .+ getproperty.(gp.LS[1].geoL.centroid, :y) .* gp.dy

    laps = set_matrices!(num, gp, [gp.LS[1].geoL], gu, [gu.LS[1].geoL], gv, [gv.LS[1].geoL], op.opC_pL, op.opC_uL, op.opC_vL, periodic_x, periodic_y)
    Lp, bc_Lp, bc_Lp_b, Lu, bc_Lu, bc_Lu_b, Lv, bc_Lv, bc_Lv_b = laps
    a0_p = [θd]


    rhs = set_poisson(BC_int, num, gp, a0_p, op.opC_pL, op.opC_uL, op.opC_vL, A, Lp, bc_Lp, bc_Lp_b, BC, true)

    


    # Test with tolerance 
    @testset "BC" begin

        @testset "Left border rhs" begin
            @test vecb_L(rhs,gp)./gp.dy[:,1] ≈ sign_Poisson * 10 * sign_left .* ones(gp.ny) atol = test_tolerance #-10: -BC because left face (divergence theorem)
        end

        @testset "Right border rhs" begin
            @test vecb_R(rhs,gp)./gp.dx[:,end] ≈ 0 * sign_right .* ones(gp.ny) atol = test_tolerance
        end

        @testset "Bottom border rhs" begin
            @test vecb_B(rhs,gp)./gp.dx[1,:] ≈ 0 * sign_bottom .* ones(gp.nx) atol = test_tolerance
        end

        @testset "Top border rhs" begin
            @test vecb_T(rhs,gp)./gp.dx[end,:] ≈ 0 * sign_bottom .* ones(gp.nx) atol = test_tolerance
        end
    end # BC


    # setting f^\\omega to the exact values of the Laplacian calculated
    # at the cell centroids (x^\\omega i j , y^\\omega i j ) 
    # and gγ to the respective analytical value at the boundary centers (xγ i j , yγ i j ).

    b = Δf.(
        x_centroid,
        y_centroid
    )
    veci(rhs,gp,1) .+= op.opC_pL.M * vec(b)

    res = zeros(size(rhs))

    #testing at (j,j): away from the wall and (j,1): at the wall


    # testing at corner
    j = 1
    II = CartesianIndex(j,j)
    tmp = lexicographic(II,gp.ny)






    i_corner_vecb_left = (num.nLS + 1) * ni + j
    
    i_corner_vecb_bottom = (num.nLS + 1) * ni + gp.ny + j # 1 = j = i
 
    if num.verbosity > 0
        print("\n A ",A[tmp,:],"\n")

        # print("\n i_corner_vecb_left ",i_corner_vecb_left)

        printstyled(color=:green, @sprintf "\n neighbours : %.3i %.3i %.3i %.3i %.3i\n" tmp tmp+1 tmp+gp.ny i_corner_vecb_bottom i_corner_vecb_left )

        print("\n A[i_corner_vecb_left,:]", A[i_corner_vecb_left,:], "\n")
        print("\n A[i_corner_vecb_bottom,:]", A[i_corner_vecb_bottom,:], "\n")
    end

    @testset "Corner bulk" begin

        @testset "A[j,j]" begin
            @test A[tmp,tmp] ≈ -6.0 atol = test_tolerance
        end

        @testset "A[j,i-1]" begin
            @test A[tmp,i_corner_vecb_left] ≈ 2.0 atol = test_tolerance
        end

        @testset "A[j-1,i]" begin
            @test A[tmp,i_corner_vecb_bottom] ≈ 2.0 atol = test_tolerance
        end

        @testset "A[j+1,i]" begin
            @test A[tmp,tmp+1] ≈ 1.0 atol = test_tolerance
        end

        @testset "A[j,i+1]" begin
            @test A[tmp,tmp+gp.ny] ≈ 1.0 atol = test_tolerance
        end

    end # Corner

    @testset "Corner border" begin

        @testset "A[i_corner_vecb_left,tmp]" begin
            @test A[i_corner_vecb_left,tmp] ≈ 2.0 atol = test_tolerance
        end

        @testset "A[i_corner_vecb_left,i_corner_vecb_left]" begin
            @test A[i_corner_vecb_left,i_corner_vecb_left] ≈ -2.0 atol = test_tolerance
        end

        @testset "A[i_corner_vecb_bottom,tmp]" begin
            @test A[i_corner_vecb_bottom,tmp] ≈ 2.0 atol = test_tolerance
        end

        @testset "A[i_corner_vecb_bottom,i_corner_vecb_bottom]" begin
            @test A[i_corner_vecb_bottom,i_corner_vecb_bottom] ≈ -2.0 atol = test_tolerance
        end


    end # Corner

    #Away from corner
    j = div(n,2)

    if num.verbosity > 0
        printstyled(color=:green, @sprintf "\n div(n,2) : %.3i \n" j)
    end

    # A   
    # [15  ]  =  1.0
    # [16  ]  =  -5.0
    # [17  ]  =  1.0
    # [48  ]  =  1.0
    # [1008]  =  2.0
    # [2064]  =  1.0
    # [2128]  =  -1.0
    
    # Test matrix coefficients at the wall 

    IIwall = CartesianIndex(j,1)
    tmp = lexicographic(IIwall,gp.ny)

    i_corner_vecb_left = (num.nLS + 1) * ni + j

    if num.verbosity > 0

        print("\n A ",A[tmp,:],"\n")

        # print("\n A ",A[tmp,tmp],"\n")
        # print("\n A ",A[tmp,tmp-1],"\n")
        # print("\n A ",A[tmp,tmp+1],"\n")
        # print("\n A ",A[tmp,tmp+gp.ny],"\n")

        # print("\n A intfc ",A[tmp,ni+1:2*ni],"\n")

        # print("\n A border ",A[tmp,end-nb+1:end],"\n")

        # testn = gp.ny -tmp+1

        # print("\n A ",A[end-nb+testn,:],"\n")

        # print("\n A[end-nb+testn,1:ni]", A[end-nb+testn,1:ni], "\n")
        # print("\n A[end-nb+testn,ni+1:2*ni]", A[end-nb+testn,ni+1:2*ni], "\n")
        # print("\n A[end-nb+testn,end-nb+1:end]", A[end-nb+testn,end-nb+1:end], "\n")

        print("\n At the wall, not at corner \n")
        printstyled(color=:green, @sprintf "\n neighbours : %.3i %.3i %.3i %.3i %.3i\n" tmp tmp+1 tmp+gp.ny tmp-1 i_corner_vecb_left )
    end

    @testset "At the wall" begin

        @testset "A[j,j]" begin
            @test A[tmp,tmp] ≈ -5.0 atol = test_tolerance
        end

        @testset "A[j,i-1]" begin
            @test A[tmp,i_corner_vecb_left] ≈ 2.0 atol = test_tolerance
        end

        @testset "A[j,j+1]" begin
            @test A[tmp,tmp+1] ≈ 1.0 atol = test_tolerance
        end

        @testset "A[j,j+gp.ny]" begin
            @test A[tmp,tmp+gp.ny] ≈ 1.0 atol = test_tolerance
        end

        # @testset "A[j,j+ni]" begin
        #     @test A[tmp,tmp+ni] ≈ 1.0 atol = test_tolerance
        # end

    end #away from the wall



    # Test matrix coefficients away from the wall
    II = CartesianIndex(j,j)
    tmp = lexicographic(II,gp.ny)

    if num.verbosity > 0
        print("\n A ",A[tmp,:],"\n")
    end

    @testset "Away from the wall" begin

        @testset "A[j,j]" begin
            @test A[tmp,tmp] ≈ -4.0 atol = test_tolerance
        end

        @testset "A[j,i-1]" begin
            @test A[tmp,tmp-1] ≈ 1.0 atol = test_tolerance
        end

        @testset "A[j,j+1]" begin
            @test A[tmp,tmp+1] ≈ 1.0 atol = test_tolerance
        end

        @testset "A[j,j-gp.ny]" begin
            @test A[tmp,tmp-gp.ny] ≈ 1.0 atol = test_tolerance
        end
        @testset "A[j,j+gp.ny]" begin
            @test A[tmp,tmp+gp.ny] ≈ 1.0 atol = test_tolerance
        end

    end #away from the wall

    # Atest = copy(A)

    # @inbounds @threads for i in 1:A.m
    #     @inbounds Atest[i,i] += 1e-10
    # end

    # if num.verbosity > 0
    #     print("\n Away from the wall, after 1e-10 added \n")
    #     print("\n A ",Atest[tmp,:],"\n")
    # end

    # #Test matrix coefficients away from the wall
    # @testset "Away from the wall, after 1e-10 added" begin

    #     @testset "A[j,j]" begin
    #         @test Atest[tmp,tmp] ≈ -4.0 atol = test_tolerance
    #     end

    #     @testset "A[j,i-1]" begin
    #         @test Atest[tmp,tmp-1] ≈ 1.0 atol = test_tolerance
    #     end

    #     @testset "A[j,j+1]" begin
    #         @test Atest[tmp,tmp+1] ≈ 1.0 atol = test_tolerance
    #     end

    #     @testset "A[j,j-gp.ny]" begin
    #         @test Atest[tmp,tmp-gp.ny] ≈ 1.0 atol = test_tolerance
    #     end
    #     @testset "A[j,j+gp.ny]" begin
    #         @test Atest[tmp,tmp+gp.ny] ≈ 1.0 atol = test_tolerance
    #     end

    # end #away from the wall

    @inbounds @threads for i in 1:A.m
        if A[i,i] == 0
            @inbounds A[i,i] += 1e-10
        end
    end

    @testset "Away from the wall, new" begin

        @testset "A[j,j]" begin
            @test A[tmp,tmp] ≈ -4.0 atol = test_tolerance
        end

        @testset "A[j,i-1]" begin
            @test A[tmp,tmp-1] ≈ 1.0 atol = test_tolerance
        end

        @testset "A[j,j+1]" begin
            @test A[tmp,tmp+1] ≈ 1.0 atol = test_tolerance
        end

        @testset "A[j,j-gp.ny]" begin
            @test A[tmp,tmp-gp.ny] ≈ 1.0 atol = test_tolerance
        end
        @testset "A[j,j+gp.ny]" begin
            @test A[tmp,tmp+gp.ny] ≈ 1.0 atol = test_tolerance
        end

    end #away from the wall
  


    @time res .= A \ rhs

    T = reshape(vec1(res, gp), gp) #bulk field
    D = reshape(vec2(res, gp), gp) #interfacial field
    Tana = f.(
        x_centroid,
        y_centroid
    )

    print("\n test val at centroid \n")


    print("\n T ",T[j,:],"\n")
    print("\n Tana ",Tana[j,:],"\n")
    print("\n D ",D[j,:],"\n")

   
 
    # testvec = A * res
    # print("\n A * res ", testvec[j],"\n")

    # print("\n A * res ", testvec[tmp],"\n")

    #Test interface
    @testset "Interface deactivated" begin
        @test maximum(abs.(D)) ≈ 0.0 atol = test_tolerance
    end

    
    # x_bc = gp.x .+ getproperty.(gp.LS[1].mid_point, :x) .* gp.dx
    # y_bc = gp.y .+ getproperty.(gp.LS[1].mid_point, :y) .* gp.dy
    
    # #Initialize bulk value
    # vec1(phL.TD,gp) .= vec(ftest.(x_centroid,y_centroid))
    # #Initialize interfacial value
    # vec2(phL.TD,gp) .= vec(ftest.(x_bc,y_bc))
    
    # #Initialize border value
    
    # vecb_L(phL.TD, gp) .= ftest.(gp.x[:,1] .- gp.x[:,1] ./ 2.0, gp.y[:,1])
    # vecb_B(phL.TD, gp) .= ftest.(gp.x[1,:], gp.y[1,:] .- gp.y[1,:] ./ 2.0)
    # vecb_R(phL.TD, gp) .= ftest.(gp.x[:,end] .+ gp.x[:,1] ./ 2.0, gp.y[:,1])
    # vecb_T(phL.TD, gp) .= ftest.(gp.x[1,:], gp.y[end,:] .+ gp.y[1,:] ./ 2.0)
    
   




    x_bc_left = gp.x[:,1] .- gp.dx[:,1] ./ 2.0

    y_bc_bottom = gp.y[1,:] .- gp.dy[1,:] ./ 2.0

    y_bc_top = gp.y[end,:] .+ gp.dy[end,:] ./ 2.0

    x_bc_right = gp.x[:,end] .+ gp.dx[:,end] ./ 2.0

    if num.verbosity > 0
        print("\n x_bc_left ",x_bc_left,"\n")

        print("\n x_bc_right ",x_bc_right,"\n")

        print("\n y_bc_bottom ",y_bc_bottom,"\n")

        print("\n y_bc_top ",y_bc_top,"\n")
    end

    #analytical solution at left border
    analytical_left = f.(x_bc_left,gp.y[:,1])

    analytical_right= f.(x_bc_right,gp.y[:,end])

    analytical_bottom= f.(gp.x[1,:],y_bc_bottom)

    analytical_top= f.(gp.x[end,:],y_bc_top)

    if num.verbosity > 0

        print("\n analytical_left ",analytical_left,"\n")

        print("\n analytical_right ",analytical_right,"\n")

        print("\n analytical_bottom ",analytical_bottom,"\n")

        print("\n analytical_top ",analytical_top,"\n")

    end

    # Test with tolerance 
    @testset "Left border" begin
        @test analytical_left ≈ vecb_L(res,gp) atol = test_tolerance_solution_absolute
    end

    @testset "Right border" begin
        @test analytical_right ≈ vecb_R(res,gp) atol = test_tolerance_solution_absolute
    end

    @testset "Bottom border" begin
        @test analytical_bottom ≈ vecb_B(res,gp) atol = test_tolerance_solution_absolute
    end

    @testset "Top border" begin
        @test analytical_top ≈ vecb_T(res,gp) atol = test_tolerance_solution_absolute
    end





    print("\n x_centroid ",x_centroid[j,:],"\n")

    printstyled(color=:green, @sprintf "\n gp.x : %.2e %.2e %.2e \n" gp.x[j,1] getproperty.(gp.LS[1].geoL.centroid, :x)[j,1] gp.dx[j,1])

    printstyled(color=:green, @sprintf "\n gp.x : %.2e %.2e %.2e %.2e %.2e \n" gp.x[j,1] x[1] gp.dx[j,1] gp.dx[j,1]/2 x[1]+gp.dx[j,1]/2)

    printstyled(color=:green, @sprintf "\n gp.x : %.2e %.2e \n" gp.x[j,2] x[1]+gp.dx[j,1]*3/2)

    print("\n Tana vecb ",Tana[j,:],"\n")




    @inline function normf_2(field, pos, cap, h)
        AVG = 0.
        RMS = 0.
        VOLUME = 0.
        MAX = 0.
        dv = 0.
        @inbounds for II in pos
            v = abs(field[II])
            dv = cap[II]*h^2
            if dv > 0.
                VOLUME += dv
                AVG += dv*v
                RMS += dv*v^2
                if (v > MAX) MAX = v end
            end
        end
        return AVG/VOLUME, sqrt(RMS/VOLUME), MAX
    end
    
    @inline function normf_2(field, pos, h)
        AVG = 0.
        RMS = 0.
        VOLUME = 0.
        MAX = 0.
        @inbounds for II in pos
            v = abs(field[II])
            VOLUME += h^2
            AVG += (h^2)*v
            RMS += (h^2)*v^2
            if (v > MAX) MAX = v end
        end
        return AVG/VOLUME, sqrt(RMS/VOLUME), MAX
    end



    for II in gp.ind.all_indices
        if gp.LS[1].geoL.cap[II,5] < 1e-12
            Tana[II] = 0.0
            T[II] = 0.0
        end
    end

    vlim = maximum(abs.(Tana))
    err = abs.(Tana .- T)

    LIQUID = gp.ind.all_indices[gp.LS[1].geoL.cap[:,:,5] .> (1-1e-16)]
    MIXED = gp.ind.all_indices[gp.LS[1].geoL.cap[:,:,5] .<= (1-1e-16) .&& gp.LS[1].geoL.cap[:,:,5] .> 1e-16]

    println("$(length(MIXED) * 100 / (length(LIQUID) + length(MIXED)))% of mixed cells")

    norm_all = normf(err, vcat(LIQUID, MIXED), gp.LS[1].geoL.cap[:,:,5], num.Δ)
    println("ALL: $norm_all")
    norm_mixed = normf(err, MIXED, gp.LS[1].geoL.cap[:,:,5], num.Δ)
    println("MIXED: $norm_mixed")
    norm_full = normf(err, LIQUID, gp.LS[1].geoL.cap[:,:,5], num.Δ)
    println("FULL: $norm_full")

    print("\n vecb ",vecb(res,gp))



    l1[i] = norm_all[1]
    l2[i] = norm_all[2]
    loo[i] = norm_all[3]

    l1_mixed[i] = norm_mixed[1]
    l2_mixed[i] = norm_mixed[2]
    loo_mixed[i] = norm_mixed[3]

    l1_full[i] = norm_full[1]
    l2_full[i] = norm_full[2]
    loo_full[i] = norm_full[3]

    for II in gp.ind.all_indices
        if abs.(Tana[II]) < 1e-16
            Tana[II] = NaN
            T[II] = NaN
            err[II] = NaN
        end
        if abs.(D[II]) < 1e-16
            D[II] = NaN
        end
    end
end #convergence


x_reg = 2^st_case:2^(st_case+n_cases-1)
y_tcks = [1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7]

conv_l1 = regression(npts, l1, x_reg)
conv_l2 = regression(npts, l2, x_reg)
conv_loo = regression(npts, loo, x_reg)

conv_l1_mixed = regression(npts, l1_mixed, x_reg)
conv_l2_mixed = regression(npts, l2_mixed, x_reg)
conv_loo_mixed = regression(npts, loo_mixed, x_reg)

conv_l1_full = regression(npts, l1_full, x_reg)
conv_l2_full = regression(npts, l2_full, x_reg)
conv_loo_full = regression(npts, loo_full, x_reg)


print("\n conv_l1 ",conv_l1)
print("\n conv_l2 ",conv_l2)
print("\n conv_loo ",conv_loo)