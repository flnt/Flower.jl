using Revise
using Flower

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
    # if isa(limitsx, Tuple{Float64, Float64}) || isa(limitsx, Tuple{Int, Int})
    #     lx = limitsx        
    # else
    #     lx = (min(x...), max(x...))
    # end
    # if isa(limitsy, Tuple{Float64, Float64}) || isa(limitsy, Tuple{Int, Int})
    #     ly = limitsy
    # else
    #     ly = (min(y...), max(y...))
    # end

    # fig = Figure(resolution = (1300, 1000))
    # ax = Axis(fig[1,1], aspect = DataAspect(), xlabel = L"x", ylabel = L"y", xtickalign = 0, 
    #     ytickalign = 0)
    # if hide
    #     hidedecorations!(ax, ticks = false, ticklabels = false)
    # end
    for i = 1:skipx+1:grid.nx+1
        lines!(ones(grid.ny+1) .* x[i], y, linewidth = linewidth, color = :black)
    end
    for i = 1:skipy+1:grid.ny+1
        lines!(x, ones(grid.nx+1) .* y[i], linewidth = linewidth, color = :black)
    end
    # for iLS in 1:num.nLS
    #     contour!(grid.x[1,:], grid.y[:,1], grid.LS[iLS].u', 
    #         levels = [0.0], color = :red, linewidth = 3.0
    #     )
    # end
    # limits!(ax, lx[1], lx[2], ly[1], ly[2])
    # colsize!(fig.layout, 1, widths(ax.scene.viewport[])[1])
    # rowsize!(fig.layout, 1, widths(ax.scene.viewport[])[2])
    # resize_to_layout!(fig)

    # return fig
end


"""
From kinetic_energy
"""
function scal_magnitude(ph, gp, gu, gv)

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
max_iter=1

x = LinRange(0, L0, n+1)
y = LinRange(0, L0, n+1)


ϵ = 0.001


num = Numerical(
    case = "Cylinder",
    x = x,
    y = y,
    CFL = 1.0,
    max_iterations = 0,
    R = 1.0,
    ϵ = ϵ,
)

gp, gu, gv = init_meshes(num)
op, phS, phL, fwd, fwdS, fwdL = init_fields(num, gp, gu, gv)

gp.LS[1].u .= 1.0

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

phL.T .= ftest.(gp.x,gp.y)

phS.T .= ftest.(gp.x,gp.y)


# @views fwd.T[1,:,:] .= phL.T.*LS[end].geoL.cap[:,:,5] .+ phS.T[:,:].*LS[end].geoS.cap[:,:,5]



x_centroid = gp.x .+ getproperty.(gp.LS[1].geoL.centroid, :x) .* gp.dx
y_centroid = gp.y .+ getproperty.(gp.LS[1].geoL.centroid, :y) .* gp.dy

x_bc = gp.x .+ getproperty.(gp.LS[1].mid_point, :x) .* gp.dx
y_bc = gp.y .+ getproperty.(gp.LS[1].mid_point, :y) .* gp.dy

vec1(phL.TD,gp) .= vec(ftest.(x_centroid,y_centroid))
vec2(phL.TD,gp) .= vec(ftest.(x_bc,y_bc))


vecb_L(phL.TD, gp) .= ftest.(gp.x[:,1] .- gp.x[:,1] ./ 2.0, gp.y[:,1])
vecb_B(phL.TD, gp) .= ftest.(gp.x[1,:], gp.y[1,:] .- gp.y[1,:] ./ 2.0)

vecb_R(phL.TD, gp) .= ftest.(gp.x[:,end] .+ gp.x[:,1] ./ 2.0, gp.y[:,1])
vecb_T(phL.TD, gp) .= ftest.(gp.x[1,:], gp.y[end,:] .+ gp.y[1,:] ./ 2.0)

# vecb_L(a,g::G) where {G<:Grid} = @view vecb(a, g)[1:g.ny]
# vecb_B(a,g::G) where {G<:Grid} = @view vecb(a, g)[g.ny+1:g.ny+g.nx]
# vecb_R(a,g::G) where {G<:Grid} = @view vecb(a, g)[g.ny+g.nx+1:2*g.ny+g.nx]
# vecb_T(a,g::G) where {G<:Grid} = @view vecb(a, g)[2*g.ny+g.nx+1:2*g.ny+2*g.nx]


# veci(phL.TD,gp,1) .= vec(phL.T)

# vec1(phL.TD,grid) .= vec(phL.T)
# vec2(phL.TD,grid) .= θd


x_centroid = gp.x .+ getproperty.(gp.LS[1].geoS.centroid, :x) .* gp.dx
y_centroid = gp.y .+ getproperty.(gp.LS[1].geoS.centroid, :y) .* gp.dy

x_bc = gp.x .+ getproperty.(gp.LS[1].mid_point, :x) .* gp.dx
y_bc = gp.y .+ getproperty.(gp.LS[1].mid_point, :y) .* gp.dy

vec1(phS.TD,gp) .= vec(ftest.(x_centroid,y_centroid))
vec2(phS.TD,gp) .= vec(ftest.(x_bc,y_bc))



compute_grad_T_x!(num,gp, gu, phL, op.opC_pL)


Tscale = 4*π/L0
curent_ticks = -1:0.25:1

# print("T x",phL.u ./Tscale)


T = Figure(size = (1600, 1000))
ax = Axis(T[1,1], aspect = DataAspect(), 
xticks = xticks, yticks = yticks,
)

co=contourf!(gu.x[1,:]./xscale, gu.y[:,1]./yscale, phL.u' ./Tscale, 
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


 Makie.save(prefix*"T_grad_x.pdf", T)



compute_grad_T_y!(num,gp, gv, phL, op.opC_pL)

print("Check NaN y ",any(isnan, phL.v))
##################################################################################

Tscale = 4*π/L0
curent_ticks = -1:0.25:1

# print("T y",phL.v ./Tscale)


T = Figure(size = (1600, 1000))
ax = Axis(T[1,1], aspect = DataAspect(), 
xticks = xticks, yticks = yticks,
)

co=contourf!(gv.x[1,:]./xscale, gv.y[:,1]./yscale, phL.v' ./Tscale, 
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


 Makie.save(prefix*"T_grad_y.pdf", T)



fgrid = plot_grid(num,gp)
Makie.save(prefix*"grid.pdf", fgrid)



Tscale = 1.0
curent_ticks = -1:0.25:1

T = Figure(size = (1600, 1000))
ax = Axis(T[1,1], aspect = DataAspect(), xticks = xticks, yticks = yticks)

co=contourf!(gp.x[1,:]./xscale, gp.y[:,1]./yscale, phL.T' ./Tscale, 
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


Makie.save(prefix*"T.pdf", T)
  


# compute_grad_T!(num,gp, phL, op.opC_pL)

Tscale = 4*π/L0
curent_ticks = -1:0.25:1

# print("T",phL.p ./Tscale)

scal_magnitude(phL,gp,gu,gv)

T = Figure(size = (1600, 1000))
ax = Axis(T[1,1], aspect = DataAspect(), xticks = xticks, yticks = yticks)

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

 Makie.save(prefix*"T_grad_mag.pdf", T)