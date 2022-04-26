function cost_functional(u, u_desired, T_all, T_desired, inside, MIXED, BC_indices, γ)
    s1 = 0.
    s2 = 0.
    s3 = 0.
    @inbounds for II in CartesianIndices(u)
        if II ∈ inside
            if II ∉ MIXED
                s1 += (T_all[II] - T_desired[II])^2
            else
                s2 += (u[II] - u_desired[II])^2
            end
        elseif II ∈ BC_indices
            s3 += (T_all[II] - T_desired[II])^2
        end
    end
    return γ[1]*s1/(2*(length(inside)-length(MIXED))) + γ[2]*s2/(2*length(MIXED)) + γ[3]*s3/(2*length(BC_indices))
end

function fg!(F, G, x, des, opt, num, idx, initial_levelset, basis)
  tmp, fwd = init_fields(num, idx)
  fwd.u .= initial_levelset
  xdata = num.H[opt.ind];
  p = curve_fit(basis, xdata, x, zeros(4))
  boundary_values = basis(num.H, p.param)

  MIXED, SOLID, LIQUID = run_forward(num, idx, tmp, fwd,
      BC_TL = Boundaries(top = Boundary(t = dir, f = dirichlet, val = boundary_values),
      bottom = Boundary(t = dir, f = dirichlet, val = boundary_values[end:-1:1]),
      right = Boundary(t = dir, f = dirichlet, val = boundary_values),
      left = Boundary(t = dir, f = dirichlet, val = boundary_values[end:-1:1])),
      stefan = true,
      heat = true,
      liquid_phase = true,
      solid_phase = true,
      verbose = false,
      advection = true
      );
  s = similar(fwd.u)
  adj = my_Adjoint(s, fwd.u, opt.γ[1].*(des.TS - fwd.TS), opt.γ[1].*(des.TL - fwd.TL), s, s, s, s)
  tmp, fwd_ = init_fields(num, idx)
  run_backward(num, idx, tmp, fwd, adj,
      stefan = true,
      heat = true,
      liquid_phase = true,
      solid_phase = true,
      verbose = false,
      advection = true
      );
  if G != nothing
      G .= gradient(adj.TL, opt, x)
      @show (p.param)
      push!(opt.p, p.param)
      push!(opt.TLsave, fwd.TL)
      push!(opt.TSsave, fwd.TS)
      push!(opt.usave, fwd.usave)
  end
  if F != nothing
    value = cost_functional(fwd, des, opt, idx, num, MIXED)
    @show (value)
    return value
  end
end

function gradient_based_optimization(x_desired, x_initial, opt, num, idx, initial_levelset, basis;
    method_opt = ConjugateGradient(),
    opt_iter = 10)

    tmp, fwd = init_fields(num, idx)
    fwd.u .= initial_levelset
    xdata = num.H[opt.ind];
    p = curve_fit(basis, xdata, x_desired, zeros(4))
    boundary_values = basis(num.H, p.param)
    @show (p.param)
    @time MIXED, SOLID, LIQUID = run_forward(num, idx, tmp, fwd,
        BC_TL = Boundaries(top = Boundary(t = dir, f = dirichlet, val = boundary_values),
        bottom = Boundary(t = dir, f = dirichlet, val = boundary_values[end:-1:1]),
        right = Boundary(t = dir, f = dirichlet, val = boundary_values),
        left = Boundary(t = dir, f = dirichlet, val = boundary_values[end:-1:1])),
        stefan = true,
        heat = true,
        liquid_phase = true,
        solid_phase = true,
        verbose = false,
        advection = true,
        show_every = 20
        );
    des = Desired(x_desired, fwd.u, fwd.usave, fwd.TL, fwd.TS)
    res = optimize(Optim.only_fg!((F, G, x)->fg!(F, G, x, des, opt, num, idx, initial_levelset, basis)), x_initial, method_opt,
    Optim.Options(store_trace = true, show_trace=true, iterations = opt_iter, allow_f_increases = true))

    @show Optim.minimizer(res)

    return res, des
end

function cost_functional(fwd, des, opt, idx, num, MIXED)
    s1 = 0.; s2 = 0.;
    @inbounds for II in idx.inside
        Tall = fwd.TL[II] + fwd.TS[II]
        Tall_d = des.TL[II] + des.TS[II]
        s1 += (Tall - Tall_d)^2 * num.Δ^2
        if II ∈ MIXED
            s2 += (fwd.u[II] - des.u[II])^2 * num.Δ^2
        end
    end
    return opt.γ[1]*s1/2 + opt.γ[2]*s2/2
end

@. gradient(field, opt, x) = -(opt.γ[3]*x + field[opt.bc_indices])
