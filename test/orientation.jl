# using Revise
using Flower

#From poisson.jl

function compute_grad_T_x_T_y_array!(num_LS, grid, grid_u, grid_v, opC_p, grad_x, grad_y, TD)
    
    ∇ϕ_x = opC_p.iMx * opC_p.Bx * vec1(TD,grid) .+ opC_p.iMx_b * opC_p.Hx_b * vecb(TD,grid)
    ∇ϕ_y = opC_p.iMy * opC_p.By * vec1(TD,grid) .+ opC_p.iMy_b * opC_p.Hy_b * vecb(TD,grid)

    for iLS in 1:num_LS
        ∇ϕ_x .+= opC_p.iMx * opC_p.Hx[iLS] * veci(TD,grid,iLS+1)
        ∇ϕ_y .+= opC_p.iMy * opC_p.Hy[iLS] * veci(TD,grid,iLS+1)
    end

    grad_x .= reshape(veci(∇ϕ_x,grid_u,1), grid_u)
    grad_y .= reshape(veci(∇ϕ_y,grid_v,1), grid_v)

end

function compute_grad_T_x_array!(num_LS, grid, grid_u, opC_p, grad_x, TD)
    
    ∇ϕ_x = opC_p.iMx * opC_p.Bx * vec1(TD,grid) .+ opC_p.iMx_b * opC_p.Hx_b * vecb(TD,grid)

    for iLS in 1:num_LS
        ∇ϕ_x .+= opC_p.iMx * opC_p.Hx[iLS] * veci(TD,grid,iLS+1)
    end

    grad_x .= reshape(veci(∇ϕ_x,grid_u,1), grid_u)

end

function compute_grad_T_y_array!(num_LS, grid, grid_v, opC_p, grad_y, TD)
    
    ∇ϕ_y = opC_p.iMy * opC_p.By * vec1(TD,grid) .+ opC_p.iMy_b * opC_p.Hy_b * vecb(TD,grid)

    for iLS in 1:num_LS
        ∇ϕ_y .+= opC_p.iMy * opC_p.Hy[iLS] * veci(TD,grid,iLS+1)
    end

    grad_y .= reshape(veci(∇ϕ_y,grid_v,1), grid_v)

end


function ftest(x,y)
    return x
end

function ftest_1(x,y)
    return 1
end



L0 = 1e-4 

L0 = 4.0

n = 64

x = LinRange(0, L0, n+1)
y = LinRange(0, L0, n+1)


ϵ = 0.001

epsilon=ϵ

radius=2.5e-5 

xcoord = 5e-5
ycoord=5e-5

xcoord = -xcoord
ycoord = -ycoord


num = Numerical(
    case = "Cylinder",
    x = x,
    y = y,
    CFL = 1.0,
    max_iterations = 0,
    shifted = xcoord,
    shifted_y = ycoord,
    R = radius,
    ϵ = ϵ,
)

gp, gu, gv = init_meshes(num)
op, phS, phL = init_fields(num, gp, gu, gv)


figname0=""

gp.LS[1].u .= 1.0 #deactivate interface

figname0="no_intfc"



@time run_forward!(
    num, gp, gu, gv, op, phS, phL;
    BC_uL = Boundaries(
        left = Dirichlet(),
        right = Dirichlet(),
        bottom = Dirichlet(),
        # top=Neumann(val=0.0),
        top=Dirichlet(),
    ),
    BC_vL = Boundaries(
        left = Dirichlet(),
        right=Dirichlet(),
        bottom = Dirichlet(),
        # top=Neumann(val=0.0),
        top=Dirichlet(),
    ),
    BC_pL = Boundaries(
        left  = Neumann(),
        right = Neumann(),
        bottom=Neumann(),
        top= Dirichlet(),
    ),
    BC_int = [WallNoSlip()],

    BC_TL = Boundaries(
    left = Dirichlet(),
    right=Dirichlet(),
    bottom = Dirichlet(),
    top = Dirichlet(),
    # left = Neumann(),
    # right=Neumann(),
    # bottom = Neumann(),
    # top = Neumann(),
    ),
    time_scheme = CN, #or FE?
    navier_stokes = true,
    ns_liquid_phase = true,
    verbose = true,
    show_every = 1,
    ns_advection = false,
)

#Initialize liquid phase
x_centroid = gp.x .+ getproperty.(gp.LS[1].geoL.centroid, :x) .* gp.dx
y_centroid = gp.y .+ getproperty.(gp.LS[1].geoL.centroid, :y) .* gp.dy

x_bc = gp.x .+ getproperty.(gp.LS[1].mid_point, :x) .* gp.dx
y_bc = gp.y .+ getproperty.(gp.LS[1].mid_point, :y) .* gp.dy

#Initialize bulk value
vec1(phL.TD,gp) .= vec(ftest.(x_centroid,y_centroid))
#Initialize interfacial value
vec2(phL.TD,gp) .= vec(ftest.(x_bc,y_bc))


vecb_L(phL.TD, gp) .= ftest.(gp.x[:,1] .- gp.x[:,1] ./ 2.0, gp.y[:,1])
vecb_B(phL.TD, gp) .= ftest.(gp.x[1,:], gp.y[1,:] .- gp.y[1,:] ./ 2.0)
vecb_R(phL.TD, gp) .= ftest.(gp.x[:,end] .+ gp.x[:,1] ./ 2.0, gp.y[:,1])
vecb_T(phL.TD, gp) .= ftest.(gp.x[1,:], gp.y[end,:] .+ gp.y[1,:] ./ 2.0)


#Initialize solid phase


x_centroid = gp.x .+ getproperty.(gp.LS[1].geoS.centroid, :x) .* gp.dx
y_centroid = gp.y .+ getproperty.(gp.LS[1].geoS.centroid, :y) .* gp.dy

x_bc = gp.x .+ getproperty.(gp.LS[1].mid_point, :x) .* gp.dx
y_bc = gp.y .+ getproperty.(gp.LS[1].mid_point, :y) .* gp.dy

vec1(phS.TD,gp) .= vec(ftest.(x_centroid,y_centroid))
vec2(phS.TD,gp) .= vec(ftest.(x_bc,y_bc))


vecb_L(phS.TD, gp) .= ftest.(gp.x[:,1] .- gp.x[:,1] ./ 2.0, gp.y[:,1])
vecb_B(phS.TD, gp) .= ftest.(gp.x[1,:], gp.y[1,:] .- gp.y[1,:] ./ 2.0)
vecb_R(phS.TD, gp) .= ftest.(gp.x[:,end] .+ gp.x[:,1] ./ 2.0, gp.y[:,1])
vecb_T(phS.TD, gp) .= ftest.(gp.x[1,:], gp.y[end,:] .+ gp.y[1,:] ./ 2.0)

phL.T .= reshape(veci(phL.TD,gp,1), gp)
phS.T .= reshape(veci(phS.TD,gp,1), gp)

# phL.T .= ftest.(gp.x,gp.y)

# phS.T .= ftest.(gp.x,gp.y)

LS   = gp.LS
LS_u = gu.LS
LS_v = gv.LS


# compute_grad_T_x!(num,gp, gu, phL, op.opC_pL)
# compute_grad_T_x!(num,gp, gu, phS, op.opC_pS)


grad_x = zeros(gu)
grad_y = zeros(gv)

compute_grad_T_x_T_y_array!(num.nLS, gp, gu, gv, op.opC_pL, grad_x, grad_y, phL.TD)

# compute_grad_T_y_array!(num_LS, gp, gv, op.opC_pL, grad_y, TD)


# grad_analytical = ftest_1.(
#     x_centroid,
#     y_centroid
# )

x_centroid_u = gu.x .+ getproperty.(gu.LS[1].geoS.centroid, :x) .* gu.dx
y_centroid_u = gu.y .+ getproperty.(gu.LS[1].geoS.centroid, :y) .* gu.dy

grad_analytical = ftest_1.(
    x_centroid_u,
    y_centroid_u
)


#tolerance for tests
test_tolerance = 1.e-13 #1.e-14


j = div(n,2)


print("\n grad_analytical ",grad_analytical[j,:],"\n")
print("\n D ",grad_x[j,:],"\n")

@testset "Simple gradient test x" begin
    @test grad_x[j,:] ≈ 1.0 atol=test_tolerance
end

@testset "Simple gradient test y" begin
    @test grad_y[j,:] ≈ 0.0 atol=test_tolerance
end

mass_flux_vec1 = fzeros(gp)
mass_flux_vecb = fzeros(gp)
mass_flux_veci = fzeros(gp)
mass_flux = zeros(gp)



tmp_vec_p = zeros(gp)
tmp_vec_p0 = zeros(gp)
tmp_vec_p1 = zeros(gp)
@testset "Phase change: mass flux" begin

@testset "Phase change: mass flux" begin
    phL.trans_scalD[:,1] .= 1.0 
    mass_flux_vec1 = fzeros(gp)
    mass_flux_vecb = fzeros(gp)
    mass_flux_veci = fzeros(gp)
    mass_flux = zeros(gp)
    # tmp_vec_p = zeros(gp)
    # tmp_vec_p0 = zeros(gp)
    # tmp_vec_p1 = zeros(gp)

  
integrate_mass_flux_over_interface_2_no_writing(num,gp,op.opC_pL,phL.trans_scalD[:,1],mass_flux_vec1,mass_flux_vecb,mass_flux_veci,tmp_vec_p,tmp_vec_p0,tmp_vec_p1,mass_flux)
# @test sum(mass_flux) == 0 
@test sum(mass_flux) ≈ 0 atol=test_tolerance
end #"Phase change: mass flux" begin

@testset "Phase change: mass flux old" begin
    phL.trans_scalD[:,1] .= 1.0 
    mass_flux_vec1 = fzeros(gp)
    mass_flux_vecb = fzeros(gp)
    mass_flux_veci = fzeros(gp)
    mass_flux = zeros(gp)
    # tmp_vec_p = zeros(gp)
    # tmp_vec_p0 = zeros(gp)
    # tmp_vec_p1 = zeros(gp)

integrate_mass_flux_over_interface_no_writing(num,gp,op.opC_pL,phL.trans_scalD[:,1],mass_flux_vec1,mass_flux_vecb,mass_flux_veci,tmp_vec_p,tmp_vec_p0,tmp_vec_p1,mass_flux)
# @test sum(mass_flux) == 0 
@test sum(mass_flux) ≈ 0 atol=test_tolerance
end #"Phase change: mass flux" begin

end #phase change

@testset "Interpolation" begin
gp.V .= 1.0
interpolate_scalar!(gp, gu, gv, gp.V, gu.V, gv.V)

@test minimum(gu.V) ≈ 1.0 atol=test_tolerance
@test maximum(gu.V) ≈ 1.0 atol=test_tolerance
@test minimum(gv.V) ≈ 1.0 atol=test_tolerance
@test maximum(gv.V) ≈ 1.0 atol=test_tolerance

# tmp_vec_p = zeros(gp) 
# tmp_vec_p0 = zeros(gp) 

tmp_vec_u = zeros(gu) 
tmp_vec_v = zeros(gv) 

# tmp_vec_u .= 1.0 
# tmp_vec_v .= 1.0
# interpolate_grid_liquid!(gp,gu,gv,tmp_vec_u, tmp_vec_v,tmp_vec_p,tmp_vec_p0)
# interpolate_grid_solid!(gp,gu,gv,tmp_vec_u, tmp_vec_v,tmp_vec_p,tmp_vec_p0)

# @test minimum(tmp_vec_p) ≈ 1.0 atol=test_tolerance
# @test maximum(tmp_vec_p) ≈ 1.0 atol=test_tolerance
# @test minimum(tmp_vec_p0) ≈ 1.0 atol=test_tolerance
# @test maximum(tmp_vec_p0) ≈ 1.0 atol=test_tolerance

# tmp_vec_u .= 1.0 
# tmp_vec_v .= 1.0
# interpolate_grid_liquid_2!(num,gp,gu.LS[end],gv.LS[end],tmp_vec_u, tmp_vec_v,tmp_vec_p,tmp_vec_p0)

# interpolate_grid_solid_2!(num,gp,gu.LS[end],gv.LS[end],tmp_vec_u, tmp_vec_v,tmp_vec_p,tmp_vec_p0)

# @test minimum(tmp_vec_p) ≈ 1.0 atol=test_tolerance
# @test maximum(tmp_vec_p) ≈ 1.0 atol=test_tolerance
# @test minimum(tmp_vec_p0) ≈ 1.0 atol=test_tolerance
# @test maximum(tmp_vec_p0) ≈ 1.0 atol=test_tolerance

tmp_vec_u .= 1.0 
tmp_vec_v .= 1.0
interpolate_grid_liquid_solid!(num,gp,gu.LS[end],gv.LS[end],tmp_vec_u, tmp_vec_v,tmp_vec_p,tmp_vec_p0)

@test minimum(tmp_vec_p) ≈ 1.0 atol=test_tolerance
@test maximum(tmp_vec_p) ≈ 1.0 atol=test_tolerance
@test minimum(tmp_vec_p0) ≈ 1.0 atol=test_tolerance
@test maximum(tmp_vec_p0) ≈ 1.0 atol=test_tolerance

# LS_u =grid_u.LS[1]
# LS_v = grid_v.LS[1]
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

# u = phL.Eu
# v = phL.Ev


for j in 1:gp.ny
    for i in 1:gp.nx
        if (tmp_vec_p[j,i] != 1.0) 
            print("\n j",j," i ",i, " tmp_vec_p ",tmp_vec_p[j,i]," tmp_vec_p0 ",tmp_vec_p0[j,i])
            print("\n LS_u.geoL.dcap[:,2:end,6] ",gu.LS[1].geoL.dcap[j,i,6]," ",gu.LS[1].geoL.dcap[j,i+1,6])
        end
    end
end

for j in 1:gp.ny
    for i in 1:gp.nx
        if (tmp_vec_p0[j,i] != 1.0) 
            print("\n j",j," i ",i, " tmp_vec_p ",tmp_vec_p[j,i]," tmp_vec_p0 ",tmp_vec_p0[j,i])
            print("\n LS_u.geoL.dcap[:,2:end,6] ",gu.LS[1].geoL.dcap[j,i,6]," ",gu.LS[1].geoL.dcap[j,i+1,6])
        end
    end
end

end






integrate_mass_flux_over_interface(num,gp,op.opC_pL,
phL.TD,
# phL.trans_scalD[:,1],
mass_flux_vec1,
mass_flux_vecb,mass_flux_veci, tmp_vec_p, tmp_vec_p0, tmp_vec_p1, mass_flux,num.index_phase_change)


@testset "Orientation BC derivative left" begin
    @test mass_flux_vec1[1,1] < 0
end

@testset "Orientation BC derivative right" begin
    @test mass_flux_vec1[1,gp.nx] > 0
end

@testset "Gradient x-component" begin
    @test grad_analytical ≈ grad_x atol = test_tolerance
end



compute_grad_T_y!(num,gp, gv, phL, op.opC_pL)
compute_grad_T_y!(num,gp, gv, phS, op.opC_pS)

printstyled(color=:red, @sprintf "\n grad y %f %f %f\n" norm(phL.v) minimum(phL.v) maximum(phL.v))
printstyled(color=:red, @sprintf "\n grad y %f %f %f\n" norm(phS.v) minimum(phS.v) maximum(phS.v))

@testset "Gradient y-component" begin
    @test maximum(abs.(phL.v)) ≈ 0.0 atol = test_tolerance
end
