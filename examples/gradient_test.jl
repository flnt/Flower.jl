using Revise
using Flower

#deprecated
#see gradient in /test


#From poisson.jl

fontsize_theme = Theme(fonts=(;regular="CMU Serif"), fontsize = 60)
set_theme!(fontsize_theme)

prefix="/local/home/pr277828/flower/gradient_test/"

isdir(prefix) || mkdir(prefix)

"""
From plot_grid
"""
function plot_grid_figtest!(fig,ax,
    num, grid;
    linewidth = 0.5,
    limitsx = false,
    limitsy = false,
    hide = false,
    skipx = 0,
    skipy = 0,
    xscale =1.0,
    yscale = 1.0,
    )

    x = grid.x_nodes ./xscale
    y = grid.y_nodes ./yscale
   
    for i = 1:skipx+1:grid.nx+1
        lines!(ones(grid.ny+1) .* x[i], y, linewidth = linewidth, color = :black)
    end
    for i = 1:skipy+1:grid.ny+1
        lines!(x, ones(grid.nx+1) .* y[i], linewidth = linewidth, color = :black)
    end
  
end


"""
From kinetic_energy
"""
function scal_magnitude(phL, phS, gp, gu, gv)
    #TODO eps, den eps
    LS_u =gu.LS[1]
    LS_v = gv.LS[1]
    phL.p .= (
        (phL.u[:,2:end].^2.0 .* LS_u.geoL.dcap[:,2:end,6] .+ 
        phL.u[:,1:end-1].^2.0 .* LS_u.geoL.dcap[:,1:end-1,6]) ./ 
        (LS_u.geoL.dcap[:,1:end-1,6] .+ LS_u.geoL.dcap[:,2:end,6] .+ 1e-8 )
    )
    phL.p .+= (
        (phL.v[2:end,:].^2.0 .* LS_v.geoL.dcap[2:end,:,7] .+ 
        phL.v[1:end-1,:].^2.0 .* LS_v.geoL.dcap[1:end-1,:,7]) ./
        (LS_v.geoL.dcap[1:end-1,:,7] .+ LS_v.geoL.dcap[2:end,:,7] .+ 1e-8 )
    )
    phL.p .+= (
        (phS.u[:,2:end].^2.0 .* LS_u.geoS.dcap[:,2:end,6] .+ 
        phS.u[:,1:end-1].^2.0 .* LS_u.geoS.dcap[:,1:end-1,6]) ./ 
        (LS_u.geoS.dcap[:,1:end-1,6] .+ LS_u.geoS.dcap[:,2:end,6] .+ 1e-8 )
    )
    phL.p .+= (
        (phS.v[2:end,:].^2.0 .* LS_v.geoS.dcap[2:end,:,7] .+ 
        phS.v[1:end-1,:].^2.0 .* LS_v.geoS.dcap[1:end-1,:,7]) ./
        (LS_v.geoS.dcap[1:end-1,:,7] .+ LS_v.geoS.dcap[2:end,:,7] .+ 1e-8 )
    )

    #####################################################################
    # phL.p .= (
    #     (phL.u[:,2:end].^2.0 .* LS_u.geoL.dcap[:,2:end,6] .+ 
    #     phL.u[:,1:end-1].^2.0 .* LS_u.geoL.dcap[:,1:end-1,6]) ./ 
    #     (LS_u.geoL.dcap[:,1:end-1,6] .+ LS_u.geoL.dcap[:,2:end,6] )
    # )
    # phL.p .+= (
    #     (phL.v[2:end,:].^2.0 .* LS_v.geoL.dcap[2:end,:,7] .+ 
    #     phL.v[1:end-1,:].^2.0 .* LS_v.geoL.dcap[1:end-1,:,7]) ./
    #     (LS_v.geoL.dcap[1:end-1,:,7] .+ LS_v.geoL.dcap[2:end,:,7] )
    # )
    # phL.p .+= (
    #     (phS.u[:,2:end].^2.0 .* LS_u.geoS.dcap[:,2:end,6] .+ 
    #     phS.u[:,1:end-1].^2.0 .* LS_u.geoS.dcap[:,1:end-1,6]) ./ 
    #     (LS_u.geoS.dcap[:,1:end-1,6] .+ LS_u.geoS.dcap[:,2:end,6] )
    # )
    # phL.p .+= (
    #     (phS.v[2:end,:].^2.0 .* LS_v.geoS.dcap[2:end,:,7] .+ 
    #     phS.v[1:end-1,:].^2.0 .* LS_v.geoS.dcap[1:end-1,:,7]) ./
    #     (LS_v.geoS.dcap[1:end-1,:,7] .+ LS_v.geoS.dcap[2:end,:,7] )
    # )

    phL.p .= sqrt.(phL.p)
end

"""
From kinetic_energy
"""
function scal_magnitude_L(ph, gp, gu, gv)

    LS_u =gu.LS[1]
    LS_v = gv.LS[1]
    # LS =gp.LS[1]

    ph.p .= (
        (ph.u[:,2:end].^2.0 .* LS_u.geoL.dcap[:,2:end,6] .+ 
        ph.u[:,1:end-1].^2.0 .* LS_u.geoL.dcap[:,1:end-1,6]) ./ 
        (LS_u.geoL.dcap[:,1:end-1,6] .+ LS_u.geoL.dcap[:,2:end,6] 
        # .+ 1e-8
        )
    )
    ph.p .+= (
        (ph.v[2:end,:].^2.0 .* LS_v.geoL.dcap[2:end,:,7] .+ 
        ph.v[1:end-1,:].^2.0 .* LS_v.geoL.dcap[1:end-1,:,7]) ./
        (LS_v.geoL.dcap[1:end-1,:,7] .+ LS_v.geoL.dcap[2:end,:,7] 
        # .+ 1e-8
        )
    )
    ph.p .= sqrt.(ph.p)
    # ph.p .= sqrt.(ph.p .* LS.geoL.dcap[:,:,5])

end


"""
  Compute norm of gradient for exchange T
"""
# function compute_grad_T!(num,grid, ph, opC_p)
    
#     @unpack nLS = num
#     @unpack T, TD = ph

#     ∇ϕ_x = opC_p.iMx * opC_p.Bx * vec1(TD,grid) .+ opC_p.iMx_b * opC_p.Hx_b * vecb(TD,grid)
#     ∇ϕ_y = opC_p.iMy * opC_p.By * vec1(TD,grid) .+ opC_p.iMy_b * opC_p.Hy_b * vecb(TD,grid)

#     for iLS in 1:nLS
#         ∇ϕ_x .+= opC_p.iMx * opC_p.Hx[iLS] * veci(TD,grid,iLS+1)
#         ∇ϕ_y .+= opC_p.iMy * opC_p.Hy[iLS] * veci(TD,grid,iLS+1)
#     end

#     ph.p .= reshape(veci(sqrt.(∇ϕ_x.^2 .+ ∇ϕ_y.^2),grid,1), grid)

#     # Gradient of pressure, eq. 17 in 
#     #"A Conservative Cartesian Cut-Cell Method for Mixed Boundary Conditions and the Incompressible Navier-Stokes Equations on Staggered Meshes"
#     #From navier_stokes_coupled.jl
#     # ∇ϕ_x = opC_u.AxT * opC_u.Rx * vec(T) .+ opC_u.Gx_b * vecb(TD,grid)
#     # ∇ϕ_y = opC_v.AyT * opC_v.Ry * vec(T) .+ opC_v.Gy_b * vecb(TD,grid)
#     # for iLS in 1:nLS
#     #     ∇ϕ_x .+= opC_u.Gx[iLS] * veci(TD,grid,iLS+1)
#     #     ∇ϕ_y .+= opC_v.Gy[iLS] * veci(TD,grid,iLS+1)
#     # end

# end

function compute_grad_T_x!(num,grid,grid_u, ph, opC_p)
    
    @unpack nLS = num
    @unpack T, TD = ph

    ∇ϕ_x = opC_p.iMx * opC_p.Bx * vec1(TD,grid) .+ opC_p.iMx_b * opC_p.Hx_b * vecb(TD,grid)
    ∇ϕ_y = opC_p.iMy * opC_p.By * vec1(TD,grid) .+ opC_p.iMy_b * opC_p.Hy_b * vecb(TD,grid)

    for iLS in 1:nLS
        ∇ϕ_x .+= opC_p.iMx * opC_p.Hx[iLS] * veci(TD,grid,iLS+1)
        ∇ϕ_y .+= opC_p.iMy * opC_p.Hy[iLS] * veci(TD,grid,iLS+1)
    end

    ph.u .= reshape(veci(∇ϕ_x,grid_u,1), grid_u)
end


function compute_grad_T_y!(num,grid, grid_v, ph, opC_p)
    
    @unpack nLS = num
    @unpack T, TD = ph

    ph.v .= 0.0

    ∇ϕ_x = opC_p.iMx * opC_p.Bx * vec1(TD,grid) .+ opC_p.iMx_b * opC_p.Hx_b * vecb(TD,grid)
    ∇ϕ_y = opC_p.iMy * opC_p.By * vec1(TD,grid) .+ opC_p.iMy_b * opC_p.Hy_b * vecb(TD,grid)

    for iLS in 1:nLS
        ∇ϕ_x .+= opC_p.iMx * opC_p.Hx[iLS] * veci(TD,grid,iLS+1)
        ∇ϕ_y .+= opC_p.iMy * opC_p.Hy[iLS] * veci(TD,grid,iLS+1)
    end

    ph.v .= reshape(veci(∇ϕ_y,grid_v,1), grid_v)
end


L0 = 1e-4 
n = 64

x = LinRange(0, L0, n+1)
y = LinRange(0, L0, n+1)


ϵ = 0.001

# ϵ = 0.0

# ϵ = 0.1

# ϵ = 1.0


# ϵ = 0.5

# ϵ = 0.05
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
op, phS, phL, fwd, fwdS, fwdL = init_fields(num, gp, gu, gv)


figname0=""

gp.LS[1].u .= 1.0
figname0="no_intfc"



# run_forward(num, gp, gu, gv, op, phS, phL, fwd, fwdS, fwdL)


@time run_forward(
    num, gp, gu, gv, op, phS, phL, fwd, fwdS, fwdL;
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
    # BC_int = [FreeSurface()], #[Wall()],
    BC_int = [Wall()],
    # BC_int = [FreeSurface()],

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



xlabel = L"x \left(\mu m \right)"
ylabel = L"y \left(\mu m \right)"

xscale = 1e-6
yscale = xscale

# xticks = 0:20:num.L0/xscale
# yticks = 0:20:num.L0/yscale

xticks = 0:20:100
yticks = 0:20:100

velscale = 1e-4 



function ftest(x,y)
    return cos(4*π/L0*x)
end


##########################################################################################
#Initialize liquid phase
##########################################################################################

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


##########################################################################################
#Initialize solid phase
##########################################################################################

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





##########################################################################################
phL.T .= reshape(veci(phL.TD,gp,1), gp)
phS.T .= reshape(veci(phS.TD,gp,1), gp)

# phL.T .= ftest.(gp.x,gp.y)

# phS.T .= ftest.(gp.x,gp.y)



LS   = gp.LS
LS_u = gu.LS
LS_v = gv.LS


@views fwd.T[1,:,:] .= phL.T.*LS[end].geoL.cap[:,:,5] .+ phS.T[:,:].*LS[end].geoS.cap[:,:,5]




compute_grad_T_x!(num,gp, gu, phL, op.opC_pL)
compute_grad_T_x!(num,gp, gu, phS, op.opC_pS)




# @views fwd.u[1,:,:]    .= phL.u.*LS_u[end].geoL.cap[:,:,5] .+ phS.u[:,:].*LS_u[end].geoS.cap[:,:,5]
# @views fwd.T[snap,:,:] .= phL.T.*LS[end].geoL.cap[:,:,5] .+ phS.T.*LS[end].geoS.cap[:,:,5]

# @views fwd.u[1,:,:] .= phL.u

# @views fwd.u[1,:,:] .= phL.u.*LS_u[end].geoL.cap[:,:,5] .+ phS.u[:,:].*LS_u[end].geoS.cap[:,:,5]
@views fwdL.u[1,:,:] .= phL.u.*LS_u[end].geoL.cap[:,:,5] .+ phS.u[:,:].*LS_u[end].geoS.cap[:,:,5]


bool = any(isnan, fwdL.u[1,:,:])
print("Check NaN x ",bool)
if bool
    printstyled(color=:red, @sprintf "\n Check NaN %s \n" bool)
else
    printstyled(color=:green, @sprintf "\n Check NaN %s \n" bool)
end

# phL.u .= phL.u .+ phS.u
# 
# phL.u .= phL.u .+ phS.u


Tscale = 4*π/L0
curent_ticks = -1:0.25:1

# print("T x",phL.u ./Tscale)


T = Figure(size = (1600, 1000))
ax = Axis(T[1,1], aspect = DataAspect(), 
xticks = xticks, yticks = yticks,
)

# co=contourf!(gu.x[1,:]./xscale, gu.y[:,1]./yscale, fwd.u[1,:,:]' ./Tscale, 
co=contourf!(gu.x[1,:]./xscale, gu.y[:,1]./yscale, fwdL.u[1,:,:]' ./Tscale, 
# co=contourf!(gu.x[1,:]./xscale, gu.y[:,1]./yscale, phL.u' ./Tscale, 
# colormap=:dense, 
# colorrange=(0.2, 1.0),
levels = curent_ticks,
)
Colorbar(T[1, 2], co, 
ticks = curent_ticks,
)
# contour!(gu.x[1,:]./xscale, gu.y[:,1]./yscale, gu.LS[1].u', levels = 0:0, color=:red, linewidth = 5);
# contour!(gu.x[1,:]./xscale, gu.y[:,1]./yscale, fwd.u[1,1,:,:]', levels = 0:0, color=:black, linewidth = 5, linestyle=:dot);
limits!(ax, 0.0, num.L0/xscale, 0.0, num.L0/yscale)

plot_grid_figtest!(T,ax,num,gu,xscale=xscale,yscale=yscale)

figname = prefix*figname0*"eps"*string(epsilon)*"T_grad_x"*".pdf"
Makie.save(figname, T)


compute_grad_T_y!(num,gp, gv, phL, op.opC_pL)
compute_grad_T_y!(num,gp, gv, phS, op.opC_pS)


printstyled(color=:red, @sprintf "\n grad y %f %f %f\n" norm(phL.v) minimum(phL.v) maximum(phL.v))
printstyled(color=:red, @sprintf "\n grad y %f %f %f\n" norm(phS.v) minimum(phS.v) maximum(phS.v))


# @views fwd.v[1,:,:] .= phL.v.*LS_v[end].geoL.cap[:,:,5] .+ phS.v[:,:].*LS_v[end].geoS.cap[:,:,5]
@views fwdL.v[1,:,:] .= phL.v.*LS_v[end].geoL.cap[:,:,5] .+ phS.v[:,:].*LS_v[end].geoS.cap[:,:,5]

# @views fwd.v[1,:,:] .= phL.v


bool = any(isnan, fwdL.v[1,:,:])
print("Check NaN y ",bool)
if bool
    printstyled(color=:red, @sprintf "\n Check NaN %s \n" bool)
else
    printstyled(color=:green, @sprintf "\n Check NaN %s \n" bool)
end

##################################################################################

Tscale = 4*π/L0
curent_ticks = -1:0.25:1

# print("T y",phL.v ./Tscale)


T = Figure(size = (1600, 1000))
ax = Axis(T[1,1], aspect = DataAspect(), 
xticks = xticks, yticks = yticks,
)

# co=contourf!(gv.x[1,:]./xscale, gv.y[:,1]./yscale, fwd.v[1,:,:]' ./Tscale, 
co=contourf!(gv.x[1,:]./xscale, gv.y[:,1]./yscale, fwdL.v[1,:,:]' ./Tscale, 
# co=contourf!(gv.x[1,:]./xscale, gv.y[:,1]./yscale, phL.v' ./Tscale, 
# colormap=:dense, 
# colorrange=(0.2, 1.0),
levels = curent_ticks,
)
Colorbar(T[1, 2], co, 
ticks = curent_ticks,
)
# contour!(gv.x[1,:]./xscale, gv.y[:,1]./yscale, gv.LS[1].u', levels = 0:0, color=:red, linewidth = 5);
# contour!(gv.x[1,:]./xscale, gv.y[:,1]./yscale, fwd.u[1,1,:,:]', levels = 0:0, color=:black, linewidth = 5, linestyle=:dot);
limits!(ax, 0.0, num.L0/xscale, 0.0, num.L0/yscale)

plot_grid_figtest!(T,ax,num,gv,xscale=xscale,yscale=yscale)


figname = prefix*figname0*"eps"*string(epsilon)*"T_grad_y"*".pdf"
Makie.save(figname, T)





fgrid = plot_grid(num,gp)

Makie.save(prefix*"grid.pdf", fgrid)



Tscale = 1.0
curent_ticks = -1:0.25:1

T = Figure(size = (1600, 1000))
ax = Axis(T[1,1], aspect = DataAspect(), xticks = xticks, yticks = yticks)

co=contourf!(gp.x[1,:]./xscale, gp.y[:,1]./yscale, fwd.T[1,:,:]' ./Tscale, 
# colormap=:dense, 
# colorrange=(0.2, 1.0),
levels = curent_ticks,
)
Colorbar(T[1, 2], co, 
ticks = curent_ticks,
)
contour!(gp.x[1,:]./xscale, gp.y[:,1]./yscale, gp.LS[1].u', levels = 0:0, color=:red, linewidth = 5);
contour!(gp.x[1,:]./xscale, gp.y[:,1]./yscale, fwd.u[1,1,:,:]', levels = 0:0, color=:black, linewidth = 5, linestyle=:dot);
limits!(ax, 0.0, num.L0/xscale, 0.0, num.L0/yscale)

plot_grid_figtest!(T,ax,num,gp,xscale=xscale,yscale=yscale)

figname = prefix*figname0*"eps"*string(epsilon)*"T"*".pdf"
Makie.save(figname, T)
  


# compute_grad_T!(num,gp, phL, op.opC_pL)

Tscale = 4*π/L0
curent_ticks = -1:0.25:1

# print("T",phL.p ./Tscale)

# scal_magnitude_L(phL,gp,gu,gv)
# scal_magnitude_L(phS,gp,gu,gv)

scal_magnitude(phL,phS,gp,gu,gv)



# @views fwd.p[1,:,:] .= phL.p.*LS[end].geoL.cap[:,:,5] .+ phS.p[:,:].*LS[end].geoS.cap[:,:,5]
# @views fwd.p[1,:,:] .= phL.p


T = Figure(size = (1600, 1000))
ax = Axis(T[1,1], aspect = DataAspect(), xticks = xticks, yticks = yticks)

# co=contourf!(gp.x[1,:]./xscale, gp.y[:,1]./yscale, fwd.p[1,:,:]' ./Tscale, 
co=contourf!(gp.x[1,:]./xscale, gp.y[:,1]./yscale, phL.p' ./Tscale, 
# colormap=:dense, 
# colorrange=(0.2, 1.0),
levels = curent_ticks,
)
Colorbar(T[1, 2], co, 
ticks = curent_ticks,
)
contour!(gp.x[1,:]./xscale, gp.y[:,1]./yscale, gp.LS[1].u', levels = 0:0, color=:red, linewidth = 5);
contour!(gp.x[1,:]./xscale, gp.y[:,1]./yscale, fwd.u[1,1,:,:]', levels = 0:0, color=:black, linewidth = 5, linestyle=:dot);
limits!(ax, 0.0, num.L0/xscale, 0.0, num.L0/yscale)

plot_grid_figtest!(T,ax,num,gp,xscale=xscale,yscale=yscale)

figname = prefix*figname0*"eps"*string(epsilon)*"T_grad_mag"*".pdf"
Makie.save(figname, T)



#  tcks = -num.L0:0.5:num.L0
tcks=xticks

vlim=1.0

 fa = Figure(figure_padding=(0, 50, 50, 50), size = (1600, 1000))
 ax = Axis(fa[1,1], aspect = 1, xlabel=L"x", ylabel=L"y", xticks = tcks, yticks = tcks,
     xgridvisible=false, ygridvisible=false)
 hmap = heatmap!(gu.x[1,:], gu.y[:,1], fwdL.u[1,:,:]'./Tscale, colorrange=(-vlim, vlim))
 cbar = fa[1,2] = Colorbar(fa, hmap, labelpadding=0)
 colsize!(fa.layout, 1, widths(ax.scene.viewport[])[1])
 rowsize!(fa.layout, 1, widths(ax.scene.viewport[])[2])
 resize_to_layout!(fa) 


figname = prefix*figname0*"eps"*string(epsilon)*"T_grad_x_heatmap"*".pdf"
Makie.save(figname, T)



# tcks = -num.L0:0.5:num.L0

vlim=1.0

 fa = Figure(figure_padding=(0, 50, 50, 50), size = (1600, 1000))
 ax = Axis(fa[1,1], aspect = 1, xlabel=L"x", ylabel=L"y", xticks = tcks, yticks = tcks,
     xgridvisible=false, ygridvisible=false)
 hmap = heatmap!(gv.x[1,:], gv.y[:,1], fwdL.v[1,:,:]'./Tscale, colorrange=(-vlim, vlim))
 cbar = fa[1,2] = Colorbar(fa, hmap, labelpadding=0)
 colsize!(fa.layout, 1, widths(ax.scene.viewport[])[1])
 rowsize!(fa.layout, 1, widths(ax.scene.viewport[])[2])
 resize_to_layout!(fa) 

figname = prefix*figname0*"eps"*string(epsilon)*"T_grad_y_heatmap"*".pdf"
Makie.save(figname, T)


fa = Figure(figure_padding=(0, 50, 50, 50), size = (1600, 1000))
ax = Axis(fa[1,1], aspect = 1, xlabel=L"x", ylabel=L"y", xticks = tcks, yticks = tcks,
    xgridvisible=false, ygridvisible=false)
hmap = heatmap!(gp.x[1,:], gp.y[:,1], phL.p'./Tscale, colorrange=(-vlim, vlim))
cbar = fa[1,2] = Colorbar(fa, hmap, labelpadding=0)
colsize!(fa.layout, 1, widths(ax.scene.viewport[])[1])
rowsize!(fa.layout, 1, widths(ax.scene.viewport[])[2])
resize_to_layout!(fa) 


figname = prefix*figname0*"eps"*string(epsilon)*"T_grad_mag_heatmap"*".pdf"
Makie.save(figname, fa)


#TODO TANA error gradient
#TODO extend 