using Revise
using Flower
using JLD2
using Optim

#System parameters 
St = 1.   
H0 = 0.05
T1 = 0.7
T2 = -0.3
TM = 0.0

#Numerical parameters
ratio = 4
L0 = 1.
n = 32
max_it = 2560 #150 for n = 16 and 600 for n = 32
# TEND = 2.34e-01
save_every = 256
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

for Ra = [1e5, 8e4, 4e4, 1e4, 5e3]

    gp, gu, gv = init_meshes(num)
    op, phS, phL, fwd, fwdS, fwdL = init_fields(num, gp, gu, gv)

    local_shift = 0.0001 + num.Δ / 2
    @. gp.LS[1].u = -gp.y - L0/2 + H0 + local_shift
    @. phL.T = T1 - (1. - num.θd)*(gp.y + L0/2) / (H0 + local_shift)
    @. phS.T = num.θd*(gp.y + L0/2 - 1.) / (H0 + local_shift - 1)


    @time run_forward(
            num, gp, gu, gv, op, phS, phL, fwd, fwdS, fwdL;
            periodic_x = true,
            BC_TL = Boundaries(
                bottom = Dirichlet(val = T1),
                left = Periodic(),
                right = Periodic(),
            ),
            BC_TS = Boundaries(
                top = Dirichlet(val = T2), 
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
    temp = copy(fwd.T[:,:,:])
    ls = copy(fwd.u[1,:,:,:])
    RB = copy(fwd.RB)

    JLD2.@save "/home/tf/Documents/RB_opt/RB_nx_$(nx)_ny_$(ny)_ratio_$(ratio)_tend_$(@sprintf("%.1e", num.τ*max_it))_TM_$(TM)_T1_$(T1)_T2_$(T2)_St_$(St)_Ra_$(@sprintf("%.1e", Ra)).jld2" num gp gu gv phS phL temp ls RB

end
