using Revise
using Flower

prefix = "/Users/alex/Documents/PhD/EnVar/drop/sims"
case = ""

fontsize_theme = Theme(fonts=(;regular="CMU Serif"), fontsize = 50)
set_theme!(fontsize_theme)

L0x = 4.0
L0y = 6.0
n = 96+1

y0 = 2.25

x = collect(LinRange(-L0x / 2, L0x / 2, n + 1))
dx = diff(x)[1]
y = collect(-L0y/2:dx:L0y/2+dx)
# y = collect(LinRange(-L0y / 2, L0y / 2, n + 1))

# σH = 1e-4
# σP = 1e-3
# nEns = 1
# nits = 50
# β = 0.9
# opt_prefix = "/Users/alex/Documents/PhD/EnVar/drop/optims"
# case = "/nEns$(nEns)_nits$(nits)_σH$(σH)_σP$(σP)_β$(β)"
# data = load(opt_prefix * case * "/data.jld2")
# A = data["x0_hist"][end,:]

θe = 90
θe2 = 135
A = zeros(8)
A[1] = 0.1
# function run_sessile(θe = 90)
    max_its = 2500
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
        δreinit = 0.65,
        g = 20.0,
        β = 0.0,
        σ = 1.3,
        ϵ = 0.01,
        ϵwall = 0.05,
        n_ext_cl = n_ext,
        NB = 24,
        nLS = 2,
    )

    gp, gu, gv = init_meshes(num)
    op, phS, phL, fwd, fwdS, fwdL = init_fields(num, gp, gu, gv)

    r = 0.5
    gp.LS[1].u .= sqrt.(gp.x.^2 + (gp.y .- y0).^2) - r * ones(gp)
    gp.LS[1].u .*= -1.0

    wall_shape = (
        A[1] .* cos.(2π .* gp.x) .+ A[2] .* sin.(2π .* gp.x) .+
        A[3] .* cos.(2π .* gp.x) .^ 2 .+ A[4] .* sin.(2π .* gp.x) .^ 2 .+
        A[5] .* cos.(2π .* gp.x) .^ 3 .+ A[6] .* sin.(2π .* gp.x) .^ 3 .+
        A[7] .* cos.(2π .* gp.x) .^ 4 .+ A[8] .* sin.(2π .* gp.x) .^ 4
    )

    # u1 = sqrt.((gp.x .+ 0.49).^2 + (gp.y .+ 0.2).^2) - r * ones(gp)
    # u2 = sqrt.((gp.x .- 0.49).^2 + (gp.y .+ 0.2).^2) - r * ones(gp)
    # gp.LS[2].u .= combine_levelsets(gp, u1, u2)

    gp.LS[2].u .= gp.y .- (L0y/2 .- 0.75) .+ wall_shape
    gp.LS[2].u .*= -1.0

    phL.u .= 0.0
    phL.v .= 0.0

    @time vol_drop,V0L = run_forward(
        num, gp, gu, gv, op, phS, phL, fwd, fwdS, fwdL;
        BC_uL = Boundaries(bottom = Navier_cl(λ = 1e-2),),
        BC_vL = Boundaries(bottom = Dirichlet(),),
        BC_pL = Boundaries(),
        BC_u = Boundaries(
            bottom = Neumann_cl(θe = θe2 * π / 180),
            top = Neumann(),
            left = Neumann_inh(),
            right = Neumann_inh()
        ),
        BC_int = [FreeSurface(), WallNoSlip(θe = θe * π / 180)],
        time_scheme = FE,
        auto_reinit = true,
        navier_stokes = true,
        ns_advection = false,
        ns_liquid_phase = true,
        verbose = true,
        show_every = 1,
        save_length = true,
        breakup = true,
    )

    suffix = "/$(θe)deg_$(num.current_i)_$(n)_reinit$(num.reinit_every)_nb$(num.nb_reinit)_δ$(num.δreinit)_σ$(num.σ)_hart_breakup"
    file = suffix*".jld2"
    # save_field(prefix*case*file, num, gp, phL, fwdL, fwd)

    tcks = -L0x/2:0.5:L0x
    lim = (L0x + num.Δ) / 2
    lim = 1.0

    fu = Figure(figure_pading = (0, 50, 0, 50), size = (1600, 1000))
    ax = Axis(fu[1,1], aspect=DataAspect(), xlabel=L"x", ylabel=L"y",
        xtickalign=0,  ytickalign=0, yticks = tcks)
    hmap = heatmap!(gu.x[1,:], gu.y[:,1], reshape(vec1(phL.uD, gu), gu)')
    for iLS in 1:num.nLS
        contour!(gu.x[1,:], gu.y[:,1], gu.LS[iLS].u', levels = 0:0, color=:red, linewidth = 3);
    end
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
    for iLS in 1:num.nLS
        contour!(gp.x[1,:], gp.y[:,1], gp.LS[iLS].u', levels = 0:0, color=:red, linewidth = 3);
    end
    limits!(ax, -1.0, 1.0, -1.0, 0.0)
    colsize!(fk.layout, 1, widths(ax.scene.viewport[])[1])
    rowsize!(fk.layout, 1, widths(ax.scene.viewport[])[2])
    resize_to_layout!(fk)

    fV = Figure(size = (1600, 1000))
    ax = Axis(fV[1,1], aspect=DataAspect(), xlabel=L"x", ylabel=L"y",
        xtickalign=0,  ytickalign=0, yticks = tcks)
    hmap = heatmap!(gp.x[1,:], gp.y[:,1], gp.V')
    contour!(gp.x[1,:], gp.y[:,1], gp.LS[1].u', levels = 0:0, color=:red, linewidth = 3);
    for iLS in 1:num.nLS
        contour!(gp.x[1,:], gp.y[:,1], gp.LS[iLS].u', levels = 0:0, color=:red, linewidth = 3);
    end
    limits!(ax, -1.0, 1.0, -1.0, 0.0)
    colsize!(fV.layout, 1, widths(ax.scene.viewport[])[1])
    rowsize!(fV.layout, 1, widths(ax.scene.viewport[])[2])
    resize_to_layout!(fV)

    fLS = Figure(size = (1600, 1000), figure_padding = (0, 40, 0, 40))
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

    fLS0 = Figure(size = (1600, 1000), figure_padding = (0, 40, 0, 40))
    ax = Axis(fLS0[1,1], aspect=DataAspect(), xlabel=L"x", ylabel=L"y",
        xtickalign=0,  ytickalign=0, yticks = tcks)
    heatmap!(gp.x[1,:], gp.y[:,1], fwd.u[3,1,:,:]')
    contour!(gp.x[1,:], gp.y[:,1], fwd.u[3,1,:,:]', levels = 0:0, color=:red, linewidth = 3);
    limits!(ax, num.x[1], num.x[end], num.y[1], num.y[end])
    colsize!(fLS0.layout, 1, widths(ax.scene.viewport[])[1])
    rowsize!(fLS0.layout, 1, widths(ax.scene.viewport[])[2])
    resize_to_layout!(fLS0)

    suffix = "$(θe)deg_$(num.current_i)_$(n)_reinit$(num.reinit_every)_nb$(num.nb_reinit)_δ$(num.δreinit)_σ$(num.σ)_hart_breakup"

    limx = num.x[end]
    limy0 = num.y[1]
    limye = num.y[end]
    # make_video(num, gv, fwd.uy, fwdL.v; title_prefix=prefix*"v_field_",
    #         title_suffix=suffix, framerate=1000÷num.save_every, limitsx=(-limx,limx), limitsy=(limy0-num.Δ/2,limye+num.Δ/2))
    # make_video(num, gp, fwd.u, fwdL.p; title_prefix=prefix*"p_field_",
    #         title_suffix=suffix, framerate=1000÷num.save_every, limitsx=(-limx,limx), limitsy=(limy0,limye), plot_levelsets=true)
    # make_video(num, gp, fwd.u, fwd.κ[1,:,:,:]; title_prefix=prefix*"k_field_",
    #         title_suffix=suffix, framerate=1000÷num.save_every, limitsx=(-limx,limx), limitsy=(limy0,limye), plot_levelsets=true)
    try
        make_video(num, gu, fwd.ux, fwdL.u; title_prefix=prefix*case*"/u_field_",
            title_suffix=suffix, framerate=1000÷num.save_every, limitsx=(-limx-num.Δ/2,limx+num.Δ/2),
            limitsy=(limy0,limye))
        # make_video(num, gv, fwd.vx, fwdL.v; title_prefix=prefix*"v_field_",
        #     title_suffix=suffix, framerate=1000÷num.save_every, limitsx=(-limx,limx),
        #     limitsy=(limy0-num.Δ/2,limye+num.Δ/2))
        # make_video(num, gp, fwd.u, fwdL.p; title_prefix=prefix*"p_field_",
        #     title_suffix=suffix, framerate=1000÷num.save_every, limitsx=(-limx,limx),
        #     limitsy=(limy0,limye))
        make_video(num, gp, fwd.u; title_prefix=prefix*case*"/none_",
            title_suffix=suffix, framerate=1000÷num.save_every, limitsx=(-limx,limx),
            limitsy=(limy0,limye), plot_levelsets=true);
        make_video(num, gp, fwd.u; title_prefix=prefix*case*"/",
            title_suffix=suffix, framerate=1000÷num.save_every, limitsx=(-limx,limx),
            limitsy=(limy0,limye), liquid=[1], solid=[2], hide_decors=true, plot_bottom=true);
    catch e
        last_snap = (num.current_i-1)÷num.save_every
        make_video(num, gu, fwd.ux, fwdL.u; title_prefix=prefix*case*"/u_field_",
            title_suffix=suffix, framerate=1000÷num.save_every, limitsx=(-limx-num.Δ/2,limx+num.Δ/2),
            limitsy=(limy0,limye), stepf = last_snap)
        # make_video(num, gv, fwd.vx, fwdL.v; title_prefix=prefix*"v_field_",
        #     title_suffix=suffix, framerate=1000÷num.save_every, limitsx=(-limx,limx),
        #     limitsy=(limy0-num.Δ/2,limye+num.Δ/2), stepf = last_snap)
        # make_video(num, gp, fwd.u, fwdL.p; title_prefix=prefix*"p_field_",
        #     title_suffix=suffix, framerate=1000÷num.save_every, limitsx=(-limx,limx),
        #     limitsy=(limy0,limye), stepf = last_snap)
        make_video(num, gp, fwd.u; title_prefix=prefix*case*"/none_",
            title_suffix=suffix, framerate=1000÷num.save_every, limitsx=(-limx,limx),
            limitsy=(limy0,limye), plot_levelsets=true, stepf = last_snap);
        make_video(num, gp, fwd.u; title_prefix=prefix*case*"/",
            title_suffix=suffix, framerate=1000÷num.save_every, limitsx=(-limx,limx),
            limitsy=(limy0,limye), stepf = last_snap, liquid=[1], solid=[2], hide_decors=true, plot_bottom=true);
    end

    # return fLS, fu, fv, phL, gu, gp, fk, fwdL
# end

# fLS, fu, fv, phL, gu, gp, fk, fwdL = run_sessile(160);

# fLS = Figure(size = (1600, 1000), figure_padding = (0, 20, 0, 20))
# ax = Axis(fLS[1,1], aspect=DataAspect(), xlabel=L"x", ylabel=L"y",
#     xtickalign=0,  ytickalign=0)#, yticks = tcks)
# contourf!(gp.x[1,:], gp.y[:,1], fwd.u[1,50,:,:]', levels = 0:0, color = :blue, linewidth = 3, extendhigh = :blue);
# contourf!(gp.x[1,:], gp.y[:,1], gp.LS[2].u', levels = 0:0, color = :gray, linewidth = 3, extendlow = :gray);
# limits!(ax, num.x[1]+num.Δ, num.x[end]-num.Δ, num.y[1]+num.Δ, num.y[end]-num.Δ)
# # limits!(ax, -0.5, 0.5, -0.96, 0.0)
# colsize!(fLS.layout, 1, widths(ax.scene.viewport[])[1])
# rowsize!(fLS.layout, 1, widths(ax.scene.viewport[])[2])
# resize_to_layout!(fLS)
# save(prefix * "/drop_optim2.pdf", fLS)


# σH = 1e-4
# σP = 1e-3
# nEns = 1
# nits = 50
# β = 0.9
# opt_prefix = "/Users/alex/Documents/PhD/EnVar/drop/optims"
# case = "/nEns$(nEns)_nits$(nits)_σH$(σH)_σP$(σP)_β$(β)"
# data = load(opt_prefix * case * "/data.jld2")


# fig = Figure(size = (1600, 1000), figure_padding = (0, 40, 0, 40))
# ax = Axis(fig[1,1], aspect=1.3,# xscale=log10, yscale=log10,
#     xlabel=L"\text{its}", ylabel=L"f(x)",
#     xticks = 0:10:100, yticks = 0:0.1:2
# )
# lines!(1:length(data["f0_hist"]), data["y0_hist"][:,1], linewidth=5)
# colsize!(fig.layout, 1, widths(ax.scene.viewport[])[1])
# rowsize!(fig.layout, 1, widths(ax.scene.viewport[])[2])
# resize_to_layout!(fig)
# save(opt_prefix * case *  "/drop_evol.pdf", fig)