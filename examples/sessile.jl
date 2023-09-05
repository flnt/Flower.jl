using Revise
using Flower

prefix = "/Users/alex/Documents/PhD/Cutcell/New_ops/robin/contact_line/sessile/"

fontsize_theme = Theme(fonts=(;regular="CMU Serif"), fontsize = 30)
set_theme!(fontsize_theme)

Rf(θ, V) = sqrt(V / (θ - sin(θ) * cos(θ)))
RR0(θ) = sqrt(π / (2 * (θ - sin(θ) * cos(θ))))
center(r, θ) = r * cos(π - θ)

h0 = 0.5
L0x = 3.0
L0y = 2.0
n = 96

x = collect(LinRange(-L0x / 2, L0x / 2, n + 1))
y = collect(LinRange(-L0y / 2, 0, n ÷ 3 + 1))

function run_sessile(θe = 90)
    _θe = acos((diff(y)[1] + cos(θe * π / 180) * h0) / h0) * 180 / π
    println("θe = $(_θe)")
    num = Numerical(
        case = "Planar",
        x = x,
        y = y,
        Re = 1.0,
        CFL = 1.0,
        max_iterations = 300,
        u_inf = 0.0,
        v_inf = 0.0,
        shifted = 0.0,
        save_every = 30,
        reinit_every = 20,
        nb_reinit = 2,
        σ = 10.0,
        ϵ = 0.05,
        θe = _θe * π / 180,
        NB = 12,
        subdomains = 2,
        overlaps = 1,
    )

    gp, gu, gv = init_meshes(num)
    op, phS, phL, fwd, fwdS, fwdL = init_fields(num, gp, gu, gv)

    r = 0.5
    gp.u .= sqrt.(gp.x.^2 + (gp.y .+ L0y / 2).^2) - r * ones(gp)
    gp.u .*= -1.0

    phL.u .= 0.0
    phL.v .= 0.0

    @time MIXED, SOLID, LIQUID = run_forward(
        num, gp, gu, gv, op, phS, phL, fwd, fwdS, fwdL;
        BC_uL = Boundaries(
            bottom = Navier_cl(λ = 1e-2),
            top = Dirichlet(),
        ),
        BC_vL = Boundaries(
            bottom = Dirichlet(),
            top = Dirichlet(),
        ),
        BC_pL = Boundaries(
            left = Dirichlet(),
            right = Dirichlet(),
        ),
        BC_u = Boundaries(
            bottom = Neumann(val = gp.dy[1,1] * cos(_θe * π / 180)),
            top = Neumann(val = gp.dy[end,end] * cos(_θe * π / 180)),
        ),
        time_scheme = FE,
        advection = true,
        navier_stokes = true,
        ns_advection = false,
        ns_solid_phase = false,
        ns_liquid_phase = true,
        free_surface = true,
        verbose = true,
        show_every = 10,
        save_length = true,
    )

    V0 = 0.5 * π * 0.5^2
    Vf = volume(gp.geoL)
    Vratio = Vf / V0

    mean_rad = 1 / abs(mean(gp.κ[gp.ind.MIXED[5:end-5]]))
    RR0_sim = mean_rad / 0.5
    RR0_teo = RR0(θe * π / 180)

    println("Vratio = $(Vratio)")
    println("mean rad = $(mean_rad)")
    println("RR0_sim = $(RR0_sim)")
    println("RR0_teo = $(RR0_teo)")

    suffix = "$(θe)deg_$(num.max_iterations)_$(n)_reinit$(num.reinit_every)_nb$(num.nb_reinit)"
    file = suffix*".jld2"
    save_field(prefix*file, num, gp, phL, fwdL, fwd)

    tcks = -num.L0/2:0.5:num.L0
    lim = (num.L0 + num.Δ) / 2
    lim = 1.0

    fu = Figure(resolution = (1600, 1000))
    colsize!(fu.layout, 1, Aspect(1, 1.0))
    ax = Axis(fu[1,1], aspect=DataAspect(), xlabel=L"x", ylabel=L"y",
        xtickalign=0,  ytickalign=0, yticks = tcks)
    heatmap!(gu.x[1,:], gu.y[:,1], phL.u')
    contour!(gu.x[1,:], gu.y[:,1], gu.u', levels = 0:0, color=:red, linewidth = 3);
    resize_to_layout!(fu)

    fv = Figure(resolution = (1600, 1000))
    colsize!(fv.layout, 1, Aspect(1, 1.0))
    ax = Axis(fv[1,1], aspect=DataAspect(), xlabel=L"x", ylabel=L"y",
        xtickalign=0,  ytickalign=0, yticks = tcks)
    heatmap!(gv.x[1,:], gv.y[:,1], phL.v')
    contour!(gv.x[1,:], gv.y[:,1], gv.u', levels = 0:0, color=:red, linewidth = 3);
    resize_to_layout!(fv)

    fp = Figure(resolution = (1600, 1000))
    colsize!(fp.layout, 1, Aspect(1, 1.0))
    ax  = Axis(fp[1,1], aspect=DataAspect(), xlabel=L"x", ylabel=L"y",
                xtickalign=0,  ytickalign=0, yticks = tcks)
    heatmap!(gp.x[1,:], gp.y[:,1], fwdL.p[end,:,:]')
    contour!(gp.x[1,:], gp.y[:,1], fwd.u[1,:,:]', levels = 0:0, color=:red, linewidth = 2);
    resize_to_layout!(fp)

    fk = Figure(resolution = (1600, 1000))
    colsize!(fk.layout, 1, Aspect(1, 1.0))
    ax = Axis(fk[1,1], aspect=DataAspect(), xlabel=L"x", ylabel=L"y",
        xtickalign=0,  ytickalign=0, yticks = tcks)
    hmap = heatmap!(gp.x[1,:], gp.y[:,1], gp.κ')#, colorrange=(-1.5, 1.5))
    contour!(gp.x[1,:], gp.y[:,1], gp.u', levels = 0:0, color=:red, linewidth = 3);
    cbar = fk[1,2] = Colorbar(fk, hmap, labelpadding=0)
    limits!(ax, num.x[1], num.x[end], num.y[1], num.y[end])
    resize_to_layout!(fk)

    fLS = Figure(resolution = (1600, 1000))
    colsize!(fLS.layout, 1, Aspect(1, 1.0))
    ax = Axis(fLS[1,1], aspect=DataAspect(), xlabel=L"x", ylabel=L"y",
        xtickalign=0,  ytickalign=0, yticks = tcks)
    arc!(Point2f(0, -1.0 + center(Rf(θe * π / 180, V0), θe * π / 180)), Rf(θe * π / 180, V0), -π, π,
        linewidth = 3, linestyle = :dash, color = :red, label = "Teo")
    contour!(gp.x[1,:], gp.y[:,1], gp.u', levels = 0:0,
        color = :black, linewidth = 3, label = "Sim");
    fLS[1,2] = Legend(fLS, ax, framevisible = false)
    limits!(ax, num.x[1], num.x[end], num.y[1], num.y[end])
    resize_to_layout!(fLS)

    fLS0 = Figure(resolution = (1600, 1000))
    colsize!(fLS0.layout, 1, Aspect(1, 1.0))
    ax = Axis(fLS0[1,1], aspect=DataAspect(), xlabel=L"x", ylabel=L"y",
        xtickalign=0,  ytickalign=0, yticks = tcks)
    heatmap!(gp.x[1,:], gp.y[:,1], fwd.u[1,:,:]')
    contour!(gp.x[1,:], gp.y[:,1], fwd.u[1,:,:]', levels = 0:0, color=:red, linewidth = 3);
    limits!(ax, num.x[1], num.x[end], num.y[1], num.y[end])
    resize_to_layout!(fLS0)


    limx = num.x[end]
    limy0 = num.y[1]
    limye = num.y[end]
    make_video(gu, fwd.ux, fwdL.u; title_prefix=prefix*"u_field_",
            title_suffix=suffix, framerate=1000÷num.save_every, limitsx=(-limx-num.Δ/2,limx+num.Δ/2), limitsy=(limy0,limye))
    make_video(gv, fwd.uy, fwdL.v; title_prefix=prefix*"v_field_",
            title_suffix=suffix, framerate=1000÷num.save_every, limitsx=(-limx,limx), limitsy=(limy0-num.Δ/2,limye+num.Δ/2))
    make_video(gp, fwd.u, fwdL.p; title_prefix=prefix*"p_field_",
            title_suffix=suffix, framerate=1000÷num.save_every, limitsx=(-limx,limx), limitsy=(limy0,limye))
    make_video(gp, fwd.u, fwd.κ; title_prefix=prefix*"k_field_",
            title_suffix=suffix, framerate=1000÷num.save_every, limitsx=(-limx,limx), limitsy=(limy0,limye))
    make_video(gp, fwd.u; title_prefix=prefix*"none_",
            title_suffix=suffix, framerate=1000÷num.save_every, limitsx=(-limx,limx), limitsy=(limy0,limye))

    return fLS
end

# for θe in 30:15:165
#     run_sessile(θe)
# end

fLS = run_sessile(60)