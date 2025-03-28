using Revise
using Flower
using JLD2
using Optim

#System parameters 
Ra = 5e4
St = 1.   
H0 = 0.3
T1 = 0.7
TM = 0.0

#Numerical parameters
ratio = 4
L0 = 1.
n = 16
max_it = 150
save_every = 1
nx = ratio * n
ny = n
x = LinRange(-ratio*L0/2, ratio*L0/2, nx+1)
y = LinRange(-L0/2, L0/2, ny+1)
num = Numerical(
    case = "Planar",
    Re = 1.0,    
    CFL = 0.4,
    x = x,
    y = y,
    max_iterations = max_it,
    u_inf = 0.0,
    θd = TM,
    save_every = save_every,
    NB = 6,
    nb_reinit = 10,
    ϵ = 0.05,
    shift = 0.0,
)

#Optimization parameters
@. model(t, p) = -abs(p[1]) - abs(p[2])*(1 - tanh(t/0.5)^2)
p_desired = [0.3, 2.0]
p_initial = [0., 0.]

γ = [1.0, 1.0, 1.0, 1.0, 1.0]

@. my_gradient(γ, field, indices, x) = -(γ[3]*x + field[indices])
# my_gradient(γ, phS.T, gp.ind.b_top[1], model(gp.x[1,:], p_initial))
@. my_costfunc(γ, temp_des, ls_des, temp_curr, ls_curr) = γ[1]*abs(temp_des - temp_curr) + γ[2]*abs(ls_des - ls_curr)
# my_costfunc(γ, fwdS_des.T[end,:,:], fwd_des.u[1,end,:,:], fwdS.T[end,:,:], fwd.u[1,end,:,:])

# p = curve_fit(model, gp.x[1,:], 100*phS.T[gp.ind.b_top[1]], rand(2))

opt_p = []
opt_S = []
opt_L = []
opt_u = []


#Optimization functions

function fg2!(F, G, x, num, basis, γ, opt_p, opt_S, opt_L, opt_u, fwd_des, fwdS_des, fwdL_des)

    gp, gu, gv = init_meshes(num)
    op, phS, phL, fwd, fwdS, fwdL = init_fields(num, gp, gu, gv)

    local_shift = 0.0001 + num.Δ / 2
    @. gp.LS[1].u = -gp.y - L0/2 + H0 + local_shift
    @. phL.T = T1 - (1. - num.θd)*(gp.y + L0/2) / (H0 + local_shift)
    @. phS.T = num.θd*(gp.y + L0/2 - 1.) / (H0 + local_shift - 1)

    p = curve_fit(basis, gp.x[1,:], x, rand(2))
    @show (p.param)

    try
        @time run_forward(
                num, gp, gu, gv, op, phS, phL, fwd, fwdS, fwdL;
                periodic_x = true,
                BC_TL = Boundaries(
                    bottom = Dirichlet(val = T1),
                    left = Periodic(),
                    right = Periodic(),
                ),
                BC_TS = Boundaries(
                    top = Dirichlet(val = basis(gp.x[1,:], p.param)), 
                    left = Periodic(),
                    right = Periodic(),
                ),
                BC_uL = Boundaries(
                    bottom = Dirichlet(val = num.u_inf),
                    top = Dirichlet(val = num.u_inf),
                    left = Periodic(),
                    right = Periodic(),
                ),
                BC_vL = Boundaries(
                    bottom = Dirichlet(val = num.v_inf),
                    top = Dirichlet(val = num.v_inf),
                    left = Periodic(),
                    right = Periodic(),
                ),
                BC_u = Boundaries(
                    left = Periodic(),
                    right = Periodic(),
                ),
                BC_pL = Boundaries(
                    left = Periodic(),
                    right = Periodic(),
                ),
                BC_int = [Stefan()],
                time_scheme = FE,
                ls_scheme = eno2,
                adaptative_t = false,
                heat = true,
                heat_convection = true,
                heat_liquid_phase = true,
                heat_solid_phase = true,
                navier_stokes = true,
                ns_advection = true,
                ns_liquid_phase = true,
                verbose = false,
                show_every = 10,
                Ra = Ra,
                St = St,
                cutoff_length = 0.9
            )

        op, phS, phL, adj, adjS, adjL = init_fields(num, gp, gu, gv)

        adj = copy(fwd)
        adjS = copy(fwdS)
        adjL = copy(fwdL)

        local_shift = 0.0001 + num.Δ / 2
        @. gp.LS[1].u = fwd.u[1,end,:,:]
        @. phL.T = fwd_des.T[end,:,:] - fwd.T[end,:,:]
        @. phS.T = fwd_des.T[end,:,:] - fwd.T[end,:,:]

        @time run_backward2(
            num, gp, gu, gv, op, phS, phL, adj, adjS, adjL;
            periodic_x = true,
            BC_TL = Boundaries(
                # bottom = Dirichlet(val = T1),
                left = Periodic(),
                right = Periodic(),
            ),
            BC_TS = Boundaries(
                top = Neumann(val = 0), 
                left = Periodic(),
                right = Periodic(),
            ),
            BC_uL = Boundaries(
                bottom = Dirichlet(val = num.u_inf),
                top = Dirichlet(val = num.u_inf),
                left = Periodic(),
                right = Periodic(),
            ),
            BC_vL = Boundaries(
                bottom = Dirichlet(val = num.v_inf),
                top = Dirichlet(val = num.v_inf),
                left = Periodic(),
                right = Periodic(),
            ),
            BC_u = Boundaries(
                left = Periodic(),
                right = Periodic(),
            ),
            BC_pL = Boundaries(
                left = Periodic(),
                right = Periodic(),
            ),
            BC_int = [Stefan()],
            time_scheme = FE,
            ls_scheme = eno2,
            adaptative_t = false,
            heat = true,
            heat_convection = true,
            heat_liquid_phase = true,
            heat_solid_phase = true,
            navier_stokes = true,
            ns_advection = true,
            ns_liquid_phase = true,
            verbose = false,
            show_every = 100,
            Ra = Ra,
            St = St,
            cutoff_length = 0.9
        )
        if G != nothing
            G .= my_gradient(γ, phS.T, gp.ind.b_top[1], x)
            push!(opt_p, p.param)
            push!(opt_S, phS.T)
            push!(opt_L, phL.T)
            push!(opt_u, gp.LS[1].u)
        end
        if F != nothing
        value = sum(my_costfunc(γ, fwdS_des.T[end,:,:], fwd_des.u[1,end,:,:], fwdS.T[end,:,:], fwd.u[1,end,:,:]))
        @show (value)
        return value
        end
    catch
        if G != nothing
            G .= 0*gp.x[1,:]
            # push!(opt_p, p.param)
            # push!(opt_S, phS.T)
            # push!(opt_L, phL.T)
            # push!(opt_u, gp.LS[1].u)
        end
        if F != nothing
            value = 1e10
            @show (value)
            return value
        end
    end
    
  end

function gradient_based_optimization2(x_initial, num, basis, γ, opt_p, opt_S, opt_L, opt_u, fwd_des, fwdS_des, fwdL_des;
    method_opt = ConjugateGradient(),
    opt_iter = 10)

    res = optimize(Optim.only_fg!((F, G, x)->fg2!(F, G, x, num, basis, γ, opt_p, opt_S, opt_L, opt_u, fwd_des, fwdS_des, fwdL_des)), x_initial, method_opt,
    Optim.Options(store_trace = true, show_trace=true, iterations = opt_iter, allow_f_increases = false))

    @show Optim.minimizer(res)

    return res
end


#1 Compute desired solution

gp, gu, gv = init_meshes(num)
op, phS, phL, fwd, fwdS, fwdL = init_fields(num, gp, gu, gv)

local_shift = 0.0001 + num.Δ / 2
@. gp.LS[1].u = -gp.y - L0/2 + H0 + local_shift
@. phL.T = T1 - (1. - num.θd)*(gp.y + L0/2) / (H0 + local_shift)
@. phS.T = num.θd*(gp.y + L0/2 - 1.) / (H0 + local_shift - 1)

fwd_des = copy(fwd)
fwdS_des = copy(fwdS)
fwdL_des = copy(fwdL)

@time run_forward(
        num, gp, gu, gv, op, phS, phL, fwd_des, fwdS_des, fwdL_des;
        periodic_x = true,
        BC_TL = Boundaries(
            bottom = Dirichlet(val = T1),
            left = Periodic(),
            right = Periodic(),
        ),
        BC_TS = Boundaries(
            top = Dirichlet(val = model(gp.x[1,:], p_desired)), 
            left = Periodic(),
            right = Periodic(),
        ),
        BC_uL = Boundaries(
            bottom = Dirichlet(val = num.u_inf),
            top = Dirichlet(val = num.u_inf),
            left = Periodic(),
            right = Periodic(),
        ),
        BC_vL = Boundaries(
            bottom = Dirichlet(val = num.v_inf),
            top = Dirichlet(val = num.v_inf),
            left = Periodic(),
            right = Periodic(),
        ),
        BC_u = Boundaries(
            left = Periodic(),
            right = Periodic(),
        ),
        BC_pL = Boundaries(
            left = Periodic(),
            right = Periodic(),
        ),
        BC_int = [Stefan()],
        time_scheme = FE,
        ls_scheme = eno2,
        adaptative_t = false,
        heat = true,
        heat_convection = true,
        heat_liquid_phase = true,
        heat_solid_phase = true,
        navier_stokes = true,
        ns_advection = true,
        ns_liquid_phase = true,
        verbose = true,
        show_every = 10,
        Ra = Ra,
        St = St,
        cutoff_length = 0.9
    )

Fdes = Figure(size = (1600, 1000))
ax = Axis(Fdes[1,1], aspect = ratio)
contourf!(gp.x[1,:], gp.y[:,1], fwd_des.T[end,:,:]', colormap=:dense, levels = 20)
contour!(gp.x[1,:], gp.y[:,1], fwd_des.u[1,end,:,:]', levels = 0:0, color=:red, linewidth = 5);
resize_to_layout!(Fdes)


#2 Run optimization !

res = gradient_based_optimization2(model(gp.x[1,:], p_initial), num, model, γ, opt_p, opt_S, opt_L, opt_u, fwd_des, fwdS_des, fwdL_des,
    opt_iter = 15,
    method_opt = LBFGS(linesearch = Optim.LineSearches.BackTracking()))

store = zeros(length(res.trace), 2)
for i in axes(store,1)
    store[i, 1] = res.trace[i].iteration
    store[i, 2] = res.trace[i].value
end

f = Figure()
fontsize_theme = Theme(fontsize = 20)
set_theme!(fontsize_theme)
ax = Axis(f[1,1], yscale = log10, xlabel = "Iteration", ylabel = L" J / J_0")

lines!(f[1,1], store[:,1], store[:,2]./store[1,2], color =:black, linewidth = 3)
scatter!(f[1,1], store[:,1], store[:,2]./store[1,2], markersize = 10, color =:black, marker=:rect)

f = current_figure()
# gp, gu, gv = init_meshes(num)
# op, phS, phL, fwd, fwdS, fwdL = init_fields(num, gp, gu, gv)

# local_shift = 0.0001 + num.Δ / 2
# @. gp.LS[1].u = -gp.y - L0/2 + H0 + local_shift
# @. phL.T = T1 - (1. - num.θd)*(gp.y + L0/2) / (H0 + local_shift)
# @. phS.T = num.θd*(gp.y + L0/2 - 1.) / (H0 + local_shift - 1)

# @time run_forward(
#         num, gp, gu, gv, op, phS, phL, fwd, fwdS, fwdL;
#         periodic_x = true,
#         BC_TL = Boundaries(
#             bottom = Dirichlet(val = T1),
#             left = Periodic(),
#             right = Periodic(),
#         ),
#         BC_TS = Boundaries(
#             top = Dirichlet(val = model(gp.x[1,:], p_initial)), 
#             left = Periodic(),
#             right = Periodic(),
#         ),
#         BC_uL = Boundaries(
#             bottom = Dirichlet(val = num.u_inf),
#             top = Dirichlet(val = num.u_inf),
#             left = Periodic(),
#             right = Periodic(),
#         ),
#         BC_vL = Boundaries(
#             bottom = Dirichlet(val = num.v_inf),
#             top = Dirichlet(val = num.v_inf),
#             left = Periodic(),
#             right = Periodic(),
#         ),
#         BC_u = Boundaries(
#             left = Periodic(),
#             right = Periodic(),
#         ),
#         BC_pL = Boundaries(
#             left = Periodic(),
#             right = Periodic(),
#         ),
#         BC_int = [Stefan()],
#         time_scheme = FE,
#         ls_scheme = eno2,
#         adaptative_t = false,
#         heat = true,
#         heat_convection = true,
#         heat_liquid_phase = true,
#         heat_solid_phase = true,
#         navier_stokes = true,
#         ns_advection = true,
#         ns_liquid_phase = true,
#         verbose = true,
#         show_every = 10,
#         Ra = Ra,
#         St = St,
#         cutoff_length = 0.9
#     )

# F1 = Figure(size = (1600, 1000))
# ax = Axis(F1[1,1], aspect = ratio)
# contourf!(gp.x[1,:], gp.y[:,1], fwd.T[end,:,:]', colormap=:dense, levels = 20)
# contour!(gp.x[1,:], gp.y[:,1], fwd.u[1,end,:,:]', levels = 0:0, color=:red, linewidth = 5);
# resize_to_layout!(F1)

# gp, gu, gv = init_meshes(num)
# op, phS, phL, adj, adjS, adjL = init_fields(num, gp, gu, gv)

# local_shift = 0.0001 + num.Δ / 2
# @. gp.LS[1].u = fwd.u[1,end,:,:]
# @. phL.T = fwd_des.T[end,:,:] - fwd.T[end,:,:]
# @. phS.T = fwd_des.T[end,:,:] - fwd.T[end,:,:]

# adj = copy(fwd)
# adjS = copy(fwdS)
# adjL = copy(fwdL)

# Fadj = Figure(size = (1600, 1000))
# ax = Axis(Fadj[1,1], aspect = ratio)
# contourf!(gp.x[1,:], gp.y[:,1], phL.T' + phS.T', colormap=:dense, levels = 20)
# contour!(gp.x[1,:], gp.y[:,1], gp.LS[1].u', levels = 0:0, color=:red, linewidth = 5);
# resize_to_layout!(Fadj)

# @time run_backward2(
#         num, gp, gu, gv, op, phS, phL, adj, adjS, adjL;
#         periodic_x = true,
#         BC_TL = Boundaries(
#             # bottom = Dirichlet(val = T1),
#             left = Periodic(),
#             right = Periodic(),
#         ),
#         BC_TS = Boundaries(
#             top = Neumann(val = 0), 
#             left = Periodic(),
#             right = Periodic(),
#         ),
#         BC_uL = Boundaries(
#             bottom = Dirichlet(val = num.u_inf),
#             top = Dirichlet(val = num.u_inf),
#             left = Periodic(),
#             right = Periodic(),
#         ),
#         BC_vL = Boundaries(
#             bottom = Dirichlet(val = num.v_inf),
#             top = Dirichlet(val = num.v_inf),
#             left = Periodic(),
#             right = Periodic(),
#         ),
#         BC_u = Boundaries(
#             left = Periodic(),
#             right = Periodic(),
#         ),
#         BC_pL = Boundaries(
#             left = Periodic(),
#             right = Periodic(),
#         ),
#         BC_int = [Stefan()],
#         time_scheme = FE,
#         ls_scheme = eno2,
#         adaptative_t = false,
#         heat = true,
#         heat_convection = true,
#         heat_liquid_phase = true,
#         heat_solid_phase = true,
#         navier_stokes = true,
#         ns_advection = true,
#         ns_liquid_phase = true,
#         verbose = true,
#         show_every = 10,
#         Ra = Ra,
#         St = St,
#         cutoff_length = 0.9
#     )

# Fadj2 = Figure(size = (1600, 1000))
# ax = Axis(Fadj2[1,1], aspect = ratio)
# contourf!(gp.x[1,:], gp.y[:,1], phL.T' + phS.T', colormap=:dense, levels = 20)
# contour!(gp.x[1,:], gp.y[:,1], gp.LS[1].u', levels = 0:0, color=:red, linewidth = 5);
# resize_to_layout!(Fadj2)
# # fwd.T[end,gp.ind.b_top[1]]