using Revise
using Flower

prefix = "/Users/alex/Documents/PhD/Cutcell/New_ops/robin/contact_line/sessile_2LS_bumpy/"
# prefix = "/Users/alex/Documents/PhD/Cutcell/New_ops/robin/contact_line/sessile_2LS/"

fontsize_theme = Theme(fonts=(;regular="CMU Serif"), fontsize = 50)
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

θe = 120
# function run_sessile(θe = 90)
    # if θe < 40
    #     max_its = 35000
    #     n_ext = 10
    #     CFL = 0.5
    # elseif θe < 100
    #     max_its = 15000
    #     n_ext = 10
    #     CFL = 0.5
    # else
    #     max_its = 5000
    #     max_its = 1
    #     n_ext = 10
    #     CFL = 0.5
    # end
    max_its = 500
    if max_its <= 100
        save_every = 1
    else
        save_every = max_its÷100
    end
    n_ext = 10
    CFL = 0.5
    num = Numerical(
        case = "Planar",
        x = x,
        y = y,
        Re = 1.0,
        CFL = CFL,
        max_iterations = max_its,
        u_inf = 0.0,
        v_inf = 0.0,
        save_every = save_every,
        reinit_every = 1,
        nb_reinit = 10,
        δreinit = 0.0,
        σ = 1.0,
        ϵ = 0.05,
        n_ext_cl = n_ext,
        NB = 24,
        nLS = 2
    )

    gp, gu, gv = init_meshes(num)
    op, phS, phL, fwd, fwdS, fwdL = init_fields(num, gp, gu, gv)

    r = 0.5
    gp.LS[1].u .= sqrt.(gp.x.^2 + (gp.y .+ L0y / 2).^2) - r * ones(gp)
    gp.LS[1].u .*= -1.0

    gp.LS[2].u .= gp.y .+ 0.85  .+ 0.1 .* cos.(2.0 .*π .* gp.x)

    phL.u .= 0.0
    phL.v .= 0.0

    @time run_forward(
        num, gp, gu, gv, op, phS, phL, fwd, fwdS, fwdL;
        BC_uL = Boundaries(
            bottom = Dirichlet(),
            top = Dirichlet(),
        ),
        BC_vL = Boundaries(
            bottom = Dirichlet(),
            top = Dirichlet(),
        ),
        BC_pL = Boundaries(),
        BC_u = Boundaries(
            bottom = Neumann_inh(),
            top = Neumann_inh(),
            left = Neumann_inh(),
            right = Neumann_inh()
        ),
        BC_int = [FreeSurface(), Wall(θe = θe * π / 180)],
        time_scheme = FE,
        auto_reinit = true,
        navier_stokes = true,
        ns_advection = false,
        ns_liquid_phase = true,
        verbose = true,
        show_every = 1,
        save_length = true,
    )

    V0 = 0.5 * π * 0.5^2
    Vf = volume(gp.LS[1].geoL)
    Vratio = Vf / V0

    mean_rad = 1 / abs(mean(gp.LS[1].κ[gp.LS[1].MIXED[5:end-5]]))
    RR0_sim = mean_rad / 0.5
    RR0_teo = RR0(θe * π / 180)

    println("Vratio = $(Vratio)")
    println("mean rad = $(mean_rad)")
    println("RR0_sim = $(RR0_sim)")
    println("RR0_teo = $(RR0_teo)")

    suffix = "$(θe)deg_$(max_its)_$(n)_reinit$(num.reinit_every)_nb$(num.nb_reinit)_δ$(num.δreinit)"
    # suffix = "$(θe)deg_$(num.max_iterations)_$(n)_reinit$(num.reinit_every)_nb$(num.nb_reinit)"
    file = suffix*".jld2"
    # save_field(prefix*file, num, gp, phL, fwdL, fwd)

    tcks = -num.L0/2:0.5:num.L0
    lim = (num.L0 + num.Δ) / 2
    lim = 1.0

    fu = Figure(size = (1600, 1000))
    ax = Axis(fu[1,1], aspect=DataAspect(), xlabel=L"x", ylabel=L"y",
        xtickalign=0,  ytickalign=0, yticks = tcks)
    hmap = heatmap!(gu.x[1,:], gu.y[:,1], reshape(vec1(phL.uD, gu), gu)')
    for iLS in 1:num.nLS
        contour!(gu.x[1,:], gu.y[:,1], gu.LS[iLS].u', levels = 0:0, color=:red, linewidth = 3);
    end
    scatter!([gu.x[id]], [gu.y[id]], color=:red)
    # cbar = fu[1,2] = Colorbar(fu, hmap, labelpadding=0)
    limits!(ax, -1, 1, -1, 0)
    colsize!(fu.layout, 1, widths(ax.scene.viewport[])[1])
    rowsize!(fu.layout, 1, widths(ax.scene.viewport[])[2])
    resize_to_layout!(fu)

    fv = Figure(size = (1600, 1000))
    ax = Axis(fv[1,1], aspect=DataAspect(), xlabel=L"x", ylabel=L"y",
        xtickalign=0,  ytickalign=0, yticks = tcks)
    hmap = heatmap!(gv.x[1,:], gv.y[:,1], reshape(vec1(phL.vD, gv), gv)')
    for iLS in 1:num.nLS
        contour!(gv.x[1,:], gv.y[:,1], gv.LS[iLS].u', levels = 0:0, color=:red, linewidth = 3);
    end
    # cbar = fv[1,2] = Colorbar(fv, hmap, labelpadding=0)
    limits!(ax, -1, 1, -1, 0)
    colsize!(fv.layout, 1, widths(ax.scene.viewport[])[1])
    rowsize!(fv.layout, 1, widths(ax.scene.viewport[])[2])
    resize_to_layout!(fv)

    fp = Figure(size = (1600, 1000))
    ax  = Axis(fp[1,1], aspect=DataAspect(), xlabel=L"x", ylabel=L"y",
                xtickalign=0,  ytickalign=0, yticks = tcks)
    hmap = heatmap!(gp.x[1,:], gp.y[:,1], reshape(vec1(fwdL.pD[end,:], gp), gp)')
    for iLS in 1:num.nLS
        contour!(gp.x[1,:], gp.y[:,1], gp.LS[iLS].u', levels = 0:0, color=:red, linewidth = 3);
    end
    cbar = fp[1,2] = Colorbar(fp, hmap, labelpadding=0)
    limits!(ax, -1, 1, -1, 0)
    colsize!(fp.layout, 1, widths(ax.scene.viewport[])[1])
    rowsize!(fp.layout, 1, widths(ax.scene.viewport[])[2])
    resize_to_layout!(fp)

    fk = Figure(size = (1600, 1000))
    ax = Axis(fk[1,1], aspect=DataAspect(), xlabel=L"x", ylabel=L"y",
        xtickalign=0,  ytickalign=0, yticks = tcks)
    hmap = heatmap!(gp.x[1,:], gp.y[:,1], gp.LS[1].κ')
    contour!(gp.x[1,:], gp.y[:,1], gp.LS[1].u', levels = 0:0, color=:red, linewidth = 3);
    limits!(ax, -1.0, 1.0, -1.0, 0.0)
    colsize!(fk.layout, 1, widths(ax.scene.viewport[])[1])
    rowsize!(fk.layout, 1, widths(ax.scene.viewport[])[2])
    resize_to_layout!(fk)

    fLS = Figure(size = (1600, 1000))
    ax = Axis(fLS[1,1], aspect=DataAspect(), xlabel=L"x", ylabel=L"y",
        xtickalign=0,  ytickalign=0, yticks = tcks)
    heatmap!(gp.x[1,:], gp.y[:,1], gp.LS[1].u')
    contour!(gp.x[1,:], gp.y[:,1], gp.LS[1].u', levels = 0:0, color = :red, linewidth = 3);
    contour!(gp.x[1,:], gp.y[:,1], gp.LS[2].u', levels = 0:0, color = :red, linewidth = 3);
    limits!(ax, num.x[1]+num.Δ, num.x[end]-num.Δ, num.y[1]+num.Δ, num.y[end]-num.Δ)
    # limits!(ax, -0.5, 0.5, -0.96, 0.0)
    colsize!(fLS.layout, 1, widths(ax.scene.viewport[])[1])
    rowsize!(fLS.layout, 1, widths(ax.scene.viewport[])[2])
    resize_to_layout!(fLS)

    fLS0 = Figure(size = (1600, 1000))
    colsize!(fLS0.layout, 1, Aspect(1, 1.0))
    ax = Axis(fLS0[1,1], aspect=DataAspect(), xlabel=L"x", ylabel=L"y",
        xtickalign=0,  ytickalign=0, yticks = tcks)
    heatmap!(gp.x[1,:], gp.y[:,1], fwd.u[3,1,:,:]')
    contour!(gp.x[1,:], gp.y[:,1], fwd.u[3,1,:,:]', levels = 0:0, color=:red, linewidth = 3);
    limits!(ax, num.x[1], num.x[end], num.y[1], num.y[end])
    colsize!(fLS0.layout, 1, widths(ax.scene.viewport[])[1])
    rowsize!(fLS0.layout, 1, widths(ax.scene.viewport[])[2])
    resize_to_layout!(fLS0)


    limx = num.x[end]
    limy0 = num.y[1]
    limye = num.y[end]
    make_video(num, gu, fwd.ux, fwdL.u; title_prefix=prefix*"u_field_",
            title_suffix=suffix, framerate=1000÷num.save_every, limitsx=(-limx-num.Δ/2,limx+num.Δ/2), limitsy=(limy0,limye))
    make_video(num, gv, fwd.uy, fwdL.v; title_prefix=prefix*"v_field_",
            title_suffix=suffix, framerate=1000÷num.save_every, limitsx=(-limx,limx), limitsy=(limy0-num.Δ/2,limye+num.Δ/2))
    make_video(num, gp, fwd.u, fwdL.p; title_prefix=prefix*"p_field_",
            title_suffix=suffix, framerate=1000÷num.save_every, limitsx=(-limx,limx), limitsy=(limy0,limye))
    make_video(num, gp, fwd.u, fwd.κ[1,:,:,:]; title_prefix=prefix*"k_field_",
            title_suffix=suffix, framerate=1000÷num.save_every, limitsx=(-limx,limx), limitsy=(limy0,limye))
    make_video(num, gp, fwd.u; title_prefix=prefix*"none_",
            title_suffix=suffix, framerate=1000÷num.save_every, limitsx=(-limx,limx), limitsy=(limy0,limye));

    # return fLS, fu, fv, phL, gu, gp, fk, fwdL
# end

# fLS, fu, fv, phL, gu, gp, fk, fwdL = run_sessile(160);