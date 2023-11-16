using Revise
using Flower

prefix = "/Users/alex/Documents/PhD/Cutcell/New_ops/robin/conservation/channel/"

fontsize_theme = Theme(fonts=(;regular="CMU Serif"), fontsize = 50)
set_theme!(fontsize_theme)

L0x = 1.2
L0y = 1.2
n = 256
nx = n
ny = n

x = collect(LinRange(-L0x / 2.0, L0x / 2.0, nx + 1))
y = collect(LinRange(-L0y / 2.0, L0y / 2.0, ny + 1))

max_its = 1000
num = Numerical(
    case = "Cylinder",
    x = x,
    y = y,
    R = 0.5,
    Re = 1e3,
    CFL = 0.45,
    max_iterations = max_its,
    u_inf = 0.0,
    v_inf = 0.0,
    uD = 0.0,
    vD = 0.0,
    shifted = 0.0,
    save_every = 10,
    reinit_every = max_its,
    nb_reinit = 0,
    ϵ = 0.0,
    NB = 12,
    subdomains = 2,
    overlaps = 1,
)

gp, gu, gv = init_meshes(num)
op, phS, phL, fwd, fwdS, fwdL = init_fields(num, gp, gu, gv)
gp.u .*= -1.0

phL.u .= -0.5 .* gu.y
phL.v .= 0.5 .* gv.x

@time MIXED, SOLID, LIQUID = run_forward(
    num, gp, gu, gv, op, phS, phL, fwd, fwdS, fwdL;
    periodic_x = true,
    periodic_y = true,
    BC_uL = Boundaries(
        bottom = Periodic(),
        top = Periodic(),
        right = Periodic(),
        left = Periodic()
    ),
    BC_vL = Boundaries(
        bottom = Periodic(),
        top = Periodic(),
        right = Periodic(),
        left = Periodic()
    ),
    BC_pL = Boundaries(
        bottom = Periodic(),
        top = Periodic(),
        right = Periodic(),
        left = Periodic()
    ),
    BC_u = Boundaries(
        bottom = Periodic(),
        top = Periodic(),
        right = Periodic(),
        left = Periodic()
    ),
    time_scheme = CN,
    navier_stokes = true,
    ns_advection = true,
    ns_liquid_phase = true,
    verbose = true,
    show_every = 1,
)

suffix = "n$(n)_its$(num.max_iterations)_ϵ$(num.ϵ)"

# file = suffix*".jld2"
# save_field(prefix*file, num, gp, phL, fwdL, fwd)

function momentum(fwdL, gu, gv)
    its = size(fwdL.u)[1]
    mom_u = zeros(its)
    mom_v = zeros(its)
    _mom_u = zeros(gu)
    _mom_v = zeros(gv)

    for i in 1:its
        _mom_u .= fwdL.u[i,:,:] .* gu.geoL.dcap[:,:,5]
        _mom_v .= fwdL.v[i,:,:] .* gv.geoL.dcap[:,:,5]
        mom_u[i] = sum(_mom_u)
        mom_v[i] = sum(_mom_v)
    end

    return mom_u, mom_v
end

function kinetic_energy(fwdL, gp, gu, gv)
    its = size(fwdL.u)[1]
    energy = zeros(its)
    _energy = zeros(gp)

    for i in 1:its
        _energy .= (
            (fwdL.u[i,:,2:end].^2.0 .* gu.geoL.dcap[:,2:end,6] .+ 
            fwdL.u[i,:,1:end-1].^2.0 .* gu.geoL.dcap[:,1:end-1,6]) ./ 
            (gu.geoL.dcap[:,1:end-1,6] .+ gu.geoL.dcap[:,2:end,6] .+ 1e-8)
        )
        _energy .+= (
            (fwdL.v[i,2:end,:].^2.0 .* gv.geoL.dcap[2:end,:,7] .+ 
            fwdL.v[i,1:end-1,:].^2.0 .* gv.geoL.dcap[1:end-1,:,7]) ./
            (gv.geoL.dcap[1:end-1,:,7] .+ gv.geoL.dcap[2:end,:,7] .+ 1e-8)
        )
        energy[i] = sum(_energy .* gp.geoL.dcap[:,:,5])
    end

    return energy ./ energy[2]
end

# data = load(prefix*"n128_its1000_Re100.0_ϵ0.05_CN_2.jld2")
# _fwdL = data["fwdPh"]
# _mom_u, _mom_v = momentum(_fwdL, gu, gv)
# _kin_energy = kinetic_energy(_fwdL, gp, gu, gv)

mom_u, mom_v = momentum(fwdL, gu, gv)
kin_energy = kinetic_energy(fwdL, gp, gu, gv)

fmu = Figure(resolution = (1200, 1000))
ax  = Axis(fmu[1,1], aspect=1, xlabel=L"t", ylabel=L"\text{Mom} _ x")
lines!(fwd.t[2:end].-fwd.t[2], mom_u[2:end], linewidth=3)
colsize!(fmu.layout, 1, Auto(1))
resize_to_layout!(fmu)
Makie.save(prefix*suffix*"_momu.pdf", fmu)

fmv = Figure(resolution = (1200, 1000))
ax  = Axis(fmv[1,1], aspect=1, xlabel=L"t", ylabel=L"\text{Mom} _ y")
lines!(fwd.t[2:end].-fwd.t[2], mom_v[2:end], linewidth=3)
colsize!(fmv.layout, 1, Auto(1))
resize_to_layout!(fmv)
Makie.save(prefix*suffix*"_momv.pdf", fmv)

fe = Figure(resolution = (1200, 1000))
ax  = Axis(fe[1,1], aspect=1, xlabel=L"t", ylabel=L"K / K _ 0")
lines!(fwd.t[2:end].-fwd.t[2], kin_energy[2:end], linewidth=3)
ylims!(ax, 1-2e-7, 1+2e-7)
colsize!(fe.layout, 1, Auto(1))
resize_to_layout!(fe)
Makie.save(prefix*suffix*"_kin.pdf", fe)

nothing