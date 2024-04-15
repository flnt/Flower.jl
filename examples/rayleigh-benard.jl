using Revise
using Flower
using JLD2

prefix = "/home/tf/Documents/Flower_figures/"

for vRa = [1e6, 1e5, 5e4, 2e4, 1e4]
# vRa = 1e6
    Ra = vRa
    St = 1.0   
    H0 = 0.05

    T1 = 0.7
    T2 = -0.3
    TM = 0.0

    ratio = 4
    L0 = 1.

    if vRa > 1e5
        n = 32
        max_it = 10000
        save_every = 100
    else
        n = 32
        max_it = 10000
        save_every = 100
    end

    nx = ratio * n
    ny = n
    
    x = LinRange(-ratio*L0/2, ratio*L0/2, nx+1)
    y = LinRange(-L0/2, L0/2, ny+1)
    
    num = Numerical(
        case = "Planar",
        Re = 1.0,    
        CFL = 0.1,
        x = x,
        y = y,
        max_iterations = max_it,
        u_inf = 0.0,
        θd = TM,
        save_every = save_every,
        NB = 6,
        nb_reinit = 50,
        ϵ = 0.05,
        shift = 0.0,
    )

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
        adaptative_t = true,
        heat = true,
        heat_convection = true,
        heat_liquid_phase = true,
        heat_solid_phase = true,
        navier_stokes = true,
        ns_advection = true,
        ns_liquid_phase = true,
        verbose = true,
        show_every = 1,
        Ra = Ra,
        St = St,
        cutoff_length = 0.8
    )

    make_video(num, gp, fwd.u, fwd.T; title_prefix=prefix*"T_field_newops_nx_$(nx)_ny_$(ny)_ratio_$(ratio)_maxiter_$(@sprintf("%.1e", max_it))_TM_$(TM)_T1_$(T1)_T2_$(T2)_St_$(St)_Ra_$(@sprintf("%.1e", Ra))",
        title_suffix="", framerate=240)
    make_video(num, gu, fwd.ux, fwdL.u; title_prefix=prefix*"u_field_newops_nx_$(nx)_ny_$(ny)_ratio_$(ratio)_maxiter_$(@sprintf("%.1e", max_it))_TM_$(TM)_T1_$(T1)_T2_$(T2)_St_$(St)_Ra_$(@sprintf("%.1e", Ra))",
        title_suffix="", framerate=240)
    make_video(num, gv, fwd.uy, fwdL.v; title_prefix=prefix*"v_field_newops_nx_$(nx)_ny_$(ny)_ratio_$(ratio)_maxiter_$(@sprintf("%.1e", max_it))_TM_$(TM)_T1_$(T1)_T2_$(T2)_St_$(St)_Ra_$(@sprintf("%.1e", Ra))",
        title_suffix="", framerate=240)

    JLD2.@save "/home/tf/Documents/Flower_figures/newops_nx_$(nx)_ny_$(ny)_ratio_$(ratio)_maxiter_$(@sprintf("%.1e", max_it))_TM_$(TM)_T1_$(T1)_T2_$(T2)_St_$(St)_Ra_$(@sprintf("%.1e", Ra)).jld2" num gp gu gv fwd Ra St
end


