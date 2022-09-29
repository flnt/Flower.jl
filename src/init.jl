function Indices(nx, ny)
    T = zeros(ny, nx)

    all_indices = CartesianIndices(T)
    inside = CartesianIndices(OffsetArray(T[2:end-1,2:end-1], CartesianIndex(2, 2):CartesianIndex(ny-1, nx-1)))

    LEFT = filter(x->x.I[2]==1, all_indices)
    LEFT2 = filter(x->x.I[2]==2, all_indices)
    BOTTOM = filter(x->x.I[1]==1, all_indices)
    BOTTOM2 = filter(x->x.I[1]==2, all_indices)
    RIGHT = filter(x->x.I[2]==nx, all_indices)
    RIGHT2 = filter(x->x.I[2]==nx-1, all_indices)
    TOP = filter(x->x.I[1]==ny, all_indices)
    TOP2 = filter(x->x.I[1]==ny-1, all_indices)

    return Indices(all_indices, inside, (LEFT, RIGHT2), (RIGHT, LEFT2), (BOTTOM, TOP2),
                (TOP, BOTTOM2), (LEFT, LEFT2), (BOTTOM, BOTTOM2), (RIGHT, RIGHT2), (TOP, TOP2))
end

function Mesh(grid, x_nodes, y_nodes)
    nx = length(x_nodes) - 1
    ny = length(y_nodes) - 1

    _dx = diff(x_nodes)
    _dy = diff(y_nodes)

    x = Matrix(transpose([x_nodes[i] + _dx[i]/2 for i = 1:nx, j = 1:ny]))
    y = [y_nodes[i] + _dy[i]/2 for i = 1:ny, j = 1:nx]

    dx = Matrix(transpose(repeat(_dx, 1, ny)))
    dy = repeat(_dy, 1, nx)

    ind = Indices(nx, ny)

    u = zeros(ny, nx)
    iso = zeros(ny, nx)
    faces = zeros(ny, nx, 4)

    SOL = zeros(ny, nx, 11)
    LIQ = ones(ny, nx, 11)
    LIQ[:,:,8:end] .*= 0.5

    dSOL = zeros(ny, nx, 11)
    dLIQ = ones(ny, nx, 11)

    liq_projection = Array{Gradient{Float64}, 2}(undef, ny, nx)
    sol_projection = Array{Gradient{Float64}, 2}(undef, ny, nx)

    sol_centroid = Array{Point{Float64}, 2}(undef, ny, nx)
    liq_centroid = Array{Point{Float64}, 2}(undef, ny, nx)
    mid_point = Array{Point{Float64}, 2}(undef, ny, nx)
    cut_points = Array{Vector{Point{Float64}}, 2}(undef, ny, nx)
    @inbounds @threads for II in eachindex(sol_centroid)
        sol_centroid[II] = Point(0.0, 0.0)
        liq_centroid[II] = Point(0.0, 0.0)
        mid_point[II] = Point(0.0, 0.0)
        cut_points[II] = [Point(0.0, 0.0), Point(0.0, 0.0)]
    end

    emptiedS = zeros(Bool, ny, nx)
    emptiedL = zeros(Bool, ny, nx)

    geoS = GeometricInfo(SOL, dSOL, sol_projection, sol_centroid, emptiedS)
    geoL = GeometricInfo(LIQ, dLIQ, liq_projection, liq_centroid, emptiedL)

    α = zeros(ny, nx)
    α .= NaN
    κ = zeros(ny, nx)
    V = zeros(ny, nx)

    ii = collect(i for i = 1:nx*ny)
    iw = collect(i for i = ny+1:nx*ny)
    is = collect(i for i = 2:nx*ny)
    iN = collect(i for i = 1:nx*ny-1)
    ie = collect(i for i = 1:nx*ny-ny)
    iwp = collect(i for i = 1:ny)
    isp = collect(i for i = 1:ny:nx*ny)
    inp = collect(i for i = ny:ny:nx*ny)
    iep = collect(i for i = nx*ny-ny+1:nx*ny)

    II = vcat(ii,iw,is,iN,ie,iwp,isp,inp,iep)

    jj = collect(i for i = 1:nx*ny)
    jw = collect(i for i = 1:nx*ny-ny)
    js = collect(i for i = 1:nx*ny-1)
    jn = collect(i for i = 2:nx*ny)
    je = collect(i for i = ny+1:nx*ny)
    jwp = collect(i for i = nx*ny-ny+1:nx*ny)
    jsp = collect(i for i = ny:ny:nx*ny)
    jnp = collect(i for i = 1:ny:nx*ny)
    jep = collect(i for i = 1:ny)

    JJ = vcat(jj,jw,js,jn,je,jwp,jsp,jnp,jep)

    a = ones(length(jj))
    b = zeros(length(jw)+length(js)+length(jn)+length(je))
    c = zeros(length(jwp)+length(jsp)+length(jnp)+length(jep))

    LSA = sparse(II,JJ,vcat(a,b,c))
    LSB = sparse(II,JJ,vcat(a,b,c))

    return Mesh{grid,Float64,Int64}(x_nodes, y_nodes, x, y, nx, ny, dx, dy, ind, u, iso, faces, 
                geoS, geoL, mid_point, cut_points, α, κ, V, LSA, LSB)
end

function init_meshes(num::NumericalParameters)
    mesh_cc = Mesh(GridCC, num.x, num.y)

    xx = vcat(num.x[1] - mesh_cc.dx[1,1]/2, num.x[1:end-1] .+ mesh_cc.dx[1,:]/2, num.x[end] + mesh_cc.dx[1,end]/2)
    mesh_stx = Mesh(GridFCx, xx, num.y)

    yy = vcat(num.y[1] - mesh_cc.dy[1,1]/2, num.y[1:end-1] .+ mesh_cc.dy[:,1]/2, num.y[end] + mesh_cc.dy[end,1]/2)
    mesh_sty = Mesh(GridFCy, num.x, yy)

    return (mesh_cc, mesh_stx, mesh_sty)
end

function init_fields(num::NumericalParameters, grid, grid_u, grid_v)
    @unpack τ, N, T_inf, u_inf, v_inf, A, R, L0, Δ, shifted, max_iterations, save_every, CFL, x_airfoil, y_airfoil = num
    @unpack x, y, nx, ny, u, ind = grid

    SCUTT = zeros(nx*ny)
    LCUTT = zeros(nx*ny)

    SCUTp = zeros(nx*ny)
    LCUTp = zeros(nx*ny)

    SCUTu = zeros(grid_u.nx*grid_u.ny)
    LCUTu = zeros(grid_u.nx*grid_u.ny)

    SCUTv = zeros(grid_v.nx*grid_v.ny)
    LCUTv = zeros(grid_v.nx*grid_v.ny)

    SCUTDx = zeros(nx*ny)
    SCUTDy = zeros(nx*ny)
    LCUTDx = zeros(nx*ny)
    LCUTDy = zeros(nx*ny)

    SCUTCT = zeros(nx*ny)
    LCUTCT = zeros(nx*ny)

    SCUTGxT = zeros(grid_u.nx*grid_u.ny)
    LCUTGxT = zeros(grid_u.nx*grid_u.ny)
    SCUTGyT = zeros(grid_v.nx*grid_v.ny)
    LCUTGyT = zeros(grid_v.nx*grid_v.ny)

    SCUTGxp = zeros(grid_u.nx*grid_u.ny)
    LCUTGxp = zeros(grid_u.nx*grid_u.ny)
    SCUTGyp = zeros(grid_v.nx*grid_v.ny)
    LCUTGyp = zeros(grid_v.nx*grid_v.ny)

    SCUTGxϕ = zeros(grid_u.nx*grid_u.ny)
    LCUTGxϕ = zeros(grid_u.nx*grid_u.ny)
    SCUTGyϕ = zeros(grid_v.nx*grid_v.ny)
    LCUTGyϕ = zeros(grid_v.nx*grid_v.ny)

    SCUTCu = zeros(grid_u.nx*grid_u.ny)
    LCUTCu = zeros(grid_u.nx*grid_u.ny)
    SCUTCv = zeros(grid_v.nx*grid_v.ny)
    LCUTCv = zeros(grid_v.nx*grid_v.ny)

    # Star stencil (p and T grid)
    ii = collect(i for i = 1:nx*ny)
    iw = collect(i for i = ny+1:nx*ny)
    is = collect(i for i = 2:nx*ny)
    iN = collect(i for i = 1:nx*ny-1)
    ie = collect(i for i = 1:nx*ny-ny)
    iwp = collect(i for i = 1:ny)
    isp = collect(i for i = 1:ny:nx*ny)
    inp = collect(i for i = ny:ny:nx*ny)
    iep = collect(i for i = nx*ny-ny+1:nx*ny)

    II = vcat(ii,iw,is,iN,ie,iwp,isp,inp,iep)

    jj = collect(i for i = 1:nx*ny)
    jw = collect(i for i = 1:nx*ny-ny)
    js = collect(i for i = 1:nx*ny-1)
    jn = collect(i for i = 2:nx*ny)
    je = collect(i for i = ny+1:nx*ny)
    jwp = collect(i for i = nx*ny-ny+1:nx*ny)
    jsp = collect(i for i = ny:ny:nx*ny)
    jnp = collect(i for i = 1:ny:nx*ny)
    jep = collect(i for i = 1:ny)

    JJ = vcat(jj,jw,js,jn,je,jwp,jsp,jnp,jep)

    a = ones(length(jj))
    _a = zeros(length(jj))
    b = zeros(length(jw)+length(js)+length(jn)+length(je))
    c = zeros(length(jwp)+length(jsp)+length(jnp)+length(jep))

    LTS = sparse(II,JJ,vcat(a,b,c))
    LTL = sparse(II,JJ,vcat(a,b,c))
    LpS = sparse(II,JJ,vcat(a,b,c))
    LpL = sparse(II,JJ,vcat(a,b,c))
    AS = sparse(II,JJ,vcat(a,b,c))
    AL = sparse(II,JJ,vcat(a,b,c))
    BS = sparse(II,JJ,vcat(a,b,c))
    BL = sparse(II,JJ,vcat(a,b,c))
    ApS = sparse(II,JJ,vcat(a,b,c))
    ApL = sparse(II,JJ,vcat(a,b,c))
    CTS = sparse(II,JJ,vcat(_a,b,c))
    CTL = sparse(II,JJ,vcat(_a,b,c))

    # Star stencil (u grid)
    ii = collect(i for i = 1:grid_u.nx*grid_u.ny)
    iw = collect(i for i = grid_u.ny+1:grid_u.nx*grid_u.ny)
    is = collect(i for i = 2:grid_u.nx*grid_u.ny)
    iN = collect(i for i = 1:grid_u.nx*grid_u.ny-1)
    ie = collect(i for i = 1:grid_u.nx*grid_u.ny-grid_u.ny)
    iwp = collect(i for i = 1:grid_u.ny)
    isp = collect(i for i = 1:grid_u.ny:grid_u.nx*grid_u.ny)
    inp = collect(i for i = grid_u.ny:grid_u.ny:grid_u.nx*grid_u.ny)
    iep = collect(i for i = grid_u.nx*grid_u.ny-grid_u.ny+1:grid_u.nx*grid_u.ny)

    II = vcat(ii,iw,is,iN,ie,iwp,isp,inp,iep)

    jj = collect(i for i = 1:grid_u.nx*grid_u.ny)
    jw = collect(i for i = 1:grid_u.nx*grid_u.ny-grid_u.ny)
    js = collect(i for i = 1:grid_u.nx*grid_u.ny-1)
    jn = collect(i for i = 2:grid_u.nx*grid_u.ny)
    je = collect(i for i = grid_u.ny+1:grid_u.nx*grid_u.ny)
    jwp = collect(i for i = grid_u.nx*grid_u.ny-grid_u.ny+1:grid_u.nx*grid_u.ny)
    jsp = collect(i for i = grid_u.ny:grid_u.ny:grid_u.nx*grid_u.ny)
    jnp = collect(i for i = 1:grid_u.ny:grid_u.nx*grid_u.ny)
    jep = collect(i for i = 1:grid_u.ny)

    JJ = vcat(jj,jw,js,jn,je,jwp,jsp,jnp,jep)

    a = ones(length(jj))
    b = zeros(length(jw)+length(js)+length(jn)+length(je))
    c = zeros(length(jwp)+length(jsp)+length(jnp)+length(jep))

    LuS = sparse(II,JJ,vcat(a,b,c))
    LuL = sparse(II,JJ,vcat(a,b,c))
    AuS = sparse(II,JJ,vcat(a,b,c))
    AuL = sparse(II,JJ,vcat(a,b,c))
    CuS = sparse(II,JJ,vcat(a,b,c))
    CuL = sparse(II,JJ,vcat(a,b,c))

    # Star stencil (v grid)
    ii = collect(i for i = 1:grid_v.nx*grid_v.ny)
    iw = collect(i for i = grid_v.ny+1:grid_v.nx*grid_v.ny)
    is = collect(i for i = 2:grid_v.nx*grid_v.ny)
    iN = collect(i for i = 1:grid_v.nx*grid_v.ny-1)
    ie = collect(i for i = 1:grid_v.nx*grid_v.ny-grid_v.ny)
    iwp = collect(i for i = 1:grid_v.ny)
    isp = collect(i for i = 1:grid_v.ny:grid_v.nx*grid_v.ny)
    inp = collect(i for i = grid_v.ny:grid_v.ny:grid_v.nx*grid_v.ny)
    iep = collect(i for i = grid_v.nx*grid_v.ny-grid_v.ny+1:grid_v.nx*grid_v.ny)

    II = vcat(ii,iw,is,iN,ie,iwp,isp,inp,iep)

    jj = collect(i for i = 1:grid_v.nx*grid_v.ny)
    jw = collect(i for i = 1:grid_v.nx*grid_v.ny-grid_v.ny)
    js = collect(i for i = 1:grid_v.nx*grid_v.ny-1)
    jn = collect(i for i = 2:grid_v.nx*grid_v.ny)
    je = collect(i for i = grid_v.ny+1:grid_v.nx*grid_v.ny)
    jwp = collect(i for i = grid_v.nx*grid_v.ny-grid_v.ny+1:grid_v.nx*grid_v.ny)
    jsp = collect(i for i = grid_v.ny:grid_v.ny:grid_v.nx*grid_v.ny)
    jnp = collect(i for i = 1:grid_v.ny:grid_v.nx*grid_v.ny)
    jep = collect(i for i = 1:grid_v.ny)

    JJ = vcat(jj,jw,js,jn,je,jwp,jsp,jnp,jep)

    a = ones(length(jj))
    b = zeros(length(jw)+length(js)+length(jn)+length(je))
    c = zeros(length(jwp)+length(jsp)+length(jnp)+length(jep))

    LvS = sparse(II,JJ,vcat(a,b,c))
    LvL = sparse(II,JJ,vcat(a,b,c))
    AvS = sparse(II,JJ,vcat(a,b,c))
    AvL = sparse(II,JJ,vcat(a,b,c))
    CvS = sparse(II,JJ,vcat(a,b,c))
    CvL = sparse(II,JJ,vcat(a,b,c))

    # 2 points stencil (p grid to u grid)
    iw = collect(i for i = grid_u.ny+1:grid_u.nx*grid_u.ny)
    ie = collect(i for i = 1:grid_u.nx*grid_u.ny-grid_u.ny)
    iwp = collect(i for i = 1:grid_u.ny)
    iep = collect(i for i = grid_u.nx*grid_u.ny-grid_u.ny+1:grid_u.nx*grid_u.ny)

    II = vcat(iw,ie,iwp,iep)

    jw = collect(i for i = 1:nx*ny)
    je = collect(i for i = 1:nx*ny)
    jwp = collect(i for i = nx*ny-ny+1:nx*ny)
    jep = collect(i for i = 1:ny)

    JJ = vcat(jw,je,jwp,jep)

    a = zeros(length(jw)+length(je)+length(jwp)+length(jep))

    GxpS = sparse(II,JJ,a)
    GxpL = sparse(II,JJ,a)
    GxTS = sparse(II,JJ,a)
    GxTL = sparse(II,JJ,a)
    GxϕS = sparse(II,JJ,a)
    GxϕL = sparse(II,JJ,a)

    # 2 points stencil (p grid to v grid)
    is = collect(i for i = 2:grid_v.ny)
    for j=2:grid_v.nx
        is = vcat(is, collect(i for i = (j-1)*grid_v.ny+2:(j-1)*grid_v.ny+grid_v.ny))
    end
    iN = collect(i for i = 1:grid_v.ny-1)
    for j=2:grid_v.nx
        iN = vcat(iN, collect(i for i = (j-1)*grid_v.ny+1:(j-1)*grid_v.ny+grid_v.ny-1))
    end
    isp = collect(i for i = 1:grid_v.ny:grid_v.nx*grid_v.ny)
    inp = collect(i for i = grid_v.ny:grid_v.ny:grid_v.nx*grid_v.ny)
    
    II = vcat(is,iN,isp,inp)

    js = collect(i for i = 1:nx*ny)
    jn = collect(i for i = 1:nx*ny)
    jsp = collect(i for i = ny:ny:nx*ny)
    jnp = collect(i for i = 1:ny:nx*ny)

    JJ = vcat(js,jn,jsp,jnp)

    a = zeros(length(js)+length(jn)+length(jsp)+length(jnp))

    GypS = sparse(II,JJ,a)
    GypL = sparse(II,JJ,a)
    GyTS = sparse(II,JJ,a)
    GyTL = sparse(II,JJ,a)
    GyϕS = sparse(II,JJ,a)
    GyϕL = sparse(II,JJ,a)

    # 2 points stencil (u grid to p grid)
    iw = collect(i for i = 1:nx*ny)
    ie = collect(i for i = 1:nx*ny)
    II = vcat(iw,ie)

    jw = collect(i for i = 1:grid_u.nx*grid_u.ny-grid_u.ny)
    je = collect(i for i = grid_u.ny+1:grid_u.nx*grid_u.ny)
    JJ = vcat(jw,je)

    a = zeros(length(jw)+length(je))

    DxuS = sparse(II,JJ,a)
    DxuL = sparse(II,JJ,a)
    ftcGxTS = sparse(II,JJ,a)
    ftcGxTL = sparse(II,JJ,a)
    E11 = sparse(II,JJ,a)
    utpS = sparse(II,JJ,a)
    utpL = sparse(II,JJ,a)

    # 2 points stencil (v grid to p grid)
    is = collect(i for i = 1:nx*ny)
    iN = collect(i for i = 1:nx*ny)
    II = vcat(is,iN)

    js = collect(i for i = 1:grid_v.ny-1)
    for j=2:grid_v.nx
        js = vcat(js, collect(i for i = (j-1)*grid_v.ny+1:(j-1)*grid_v.ny+grid_v.ny-1))
    end
    jn = collect(i for i = 2:grid_v.ny)
    for j=2:grid_v.nx
        jn = vcat(jn, collect(i for i = (j-1)*grid_v.ny+2:(j-1)*grid_v.ny+grid_v.ny))
    end
    JJ = vcat(js,jn)

    a = zeros(length(js)+length(jn))

    DyvS = sparse(II,JJ,a)
    DyvL = sparse(II,JJ,a)
    ftcGyTS = sparse(II,JJ,a)
    ftcGyTL = sparse(II,JJ,a)
    E22 = sparse(II,JJ,a)
    vtpS = sparse(II,JJ,a)
    vtpL = sparse(II,JJ,a)

    # 2 points stencil (u grid to p' grid)
    is = collect(i for i = 1:nx*ny)
    iN = collect(i for i = 1:nx*ny)
    iF = [nx*ny]
    II = vcat(is,iN,iF)

    js = collect(i for i = 1:grid_u.nx*grid_u.ny-grid_u.ny)
    jn = collect(i for i = 2:grid_u.nx*grid_u.ny-grid_u.ny+1)
    jf = [grid_u.nx*grid_u.ny]
    JJ = vcat(js,jn,jf)

    a = zeros(length(js)+length(jn)+1)

    E12_x = sparse(II,JJ,a)

    # 2 points stencil (v grid to p' grid)
    iw = collect(i for i = ny+1:nx*ny)
    ie = collect(i for i = 1:nx*ny)
    II = vcat(iw,ie)

    jw = []
    for j=2:grid_v.nx
        jw = vcat(jw, collect(i for i = (j-1)*grid_v.ny-(grid_v.ny-2):(j-1)*grid_v.ny))
    end
    je = []
    for j=1:grid_v.nx
        je = vcat(je, collect(i for i = j*grid_v.ny-(grid_v.ny-2):j*grid_v.ny))
    end
    JJ = vcat(jw,je)

    a = zeros(length(jw)+length(je))

    E12_y = sparse(II,JJ,a)

    TS = zeros(ny, nx)
    TL = zeros(ny, nx)
    pS = zeros(ny, nx)
    pL = zeros(ny, nx)
    ϕS = zeros(ny, nx)
    ϕL = zeros(ny, nx)
    Gxm1S = zeros(grid_u.nx*grid_u.ny)
    Gym1S = zeros(grid_v.nx*grid_v.ny)
    Gxm1L = zeros(grid_u.nx*grid_u.ny)
    Gym1L = zeros(grid_v.nx*grid_v.ny)
    uS = zeros(grid_u.ny, grid_u.nx)
    uL = zeros(grid_u.ny, grid_u.nx)
    vS = zeros(grid_v.ny, grid_v.nx)
    vL = zeros(grid_v.ny, grid_v.nx)
    ucorrS = zeros(grid_u.ny, grid_u.nx)
    ucorrL = zeros(grid_u.ny, grid_u.nx)
    vcorrS = zeros(grid_v.ny, grid_v.nx)
    vcorrL = zeros(grid_v.ny, grid_v.nx)
    Tall = zeros(ny, nx)
    DTS = zeros(ny, nx)
    DTL = zeros(ny, nx)
    DϕS = zeros(ny, nx)
    DϕL = zeros(ny, nx)
    DuS = zeros(grid_u.ny, grid_u.nx)
    DuL = zeros(grid_u.ny, grid_u.nx)
    DvS = zeros(grid_v.ny, grid_v.nx)
    DvL = zeros(grid_v.ny, grid_v.nx)
    # tmpS = zeros(grid_u.ny, grid_u.nx)
    # tmpL = zeros(grid_u.ny, grid_u.nx)
    tmpS = zeros(ny, nx)
    tmpL = zeros(ny, nx)

    n_snaps = iszero(max_iterations%save_every) ? max_iterations÷save_every+1 : max_iterations÷save_every+2
    
    usave = zeros(n_snaps, ny, nx)
    uusave = zeros(n_snaps, grid_u.ny, grid_u.nx)
    uvsave = zeros(n_snaps, grid_v.ny, grid_v.nx)
    TSsave = zeros(n_snaps, ny, nx)
    TLsave = zeros(n_snaps, ny, nx)
    Tsave = zeros(n_snaps, ny, nx)
    psave = zeros(n_snaps, ny, nx)
    ϕsave = zeros(n_snaps, ny, nx)
    Uxsave = zeros(n_snaps, grid_u.ny, grid_u.nx)
    Uysave = zeros(n_snaps, grid_v.ny, grid_v.nx)
    Uxcorrsave = zeros(n_snaps, grid_u.ny, grid_u.nx)
    Uycorrsave = zeros(n_snaps, grid_v.ny, grid_v.nx)
    Vsave = zeros(n_snaps, ny, nx)
    κsave = zeros(n_snaps, ny, nx)
    lengthsave = zeros(n_snaps)
    time = zeros(n_snaps)
    Cd = zeros(n_snaps)
    Cl = zeros(n_snaps)

    TL[:,:] = zeros(ny, nx)
    TS[:,:] = zeros(ny, nx)
    Tall[:,:] = zeros(ny, nx)
    if num.case == "Planar"
        u .= y .+ shifted
    elseif num.case == "Square"
        @. u = max((-y - R),(y - R),(-x - R),(x - R))
    elseif num.case == "Sphere"
        u .= sqrt.((x .+ shifted).^ 2 + y .^ 2) - (R) * ones(ny, nx);
        init_franck!(grid, TL, R, T_inf, 0)

        # su = sqrt.(Xu.^2 .+ Yu.^2)
        # R1 = R + 3.0*Δ
        # uL .= u_inf
        # vL .= v_inf
        # bl = 4.0
        # for II in idxu.all_indices
        #     if su[II] <= R1
        #         uL[II] = 0.0
        #     elseif su[II] > R1
        #         uL[II] = tanh(bl*(su[II]-R1))
        #     end
        # end
    elseif num.case == "Cylinder"
        u .= sqrt.((x .+ shifted).^ 2 + y .^ 2) - (R) * ones(ny, nx);
        init_franck!(grid, TL, R, T_inf, 0)

        su = sqrt.((grid_u.x .+ shifted).^2 .+ grid_u.y.^2)
        R1 = R + 3.0*Δ
        uL .= u_inf
        vL .= v_inf
        bl = 4.0
        for II in grid_u.ind.all_indices
            if su[II] <= R1
                uL[II] = 0.0
            elseif su[II] > R1
                uL[II] = tanh(bl*(su[II]-R1))
            end
        end
    elseif num.case == "Ellipse"
        u .= sqrt.((x .+ shifted) .^ 2 + (2.0 .* y) .^ 2) - (R) * ones(ny, nx);
    elseif num.case == "Mullins"
        u .= y .+ shifted .+ A*sin.(N*pi*x) .+ L0/4;
        init_mullins!(grid, TL, T_inf, 0., A, N, L0/4)
    elseif num.case == "Mullins_cos"
        @. u = y + 0.6L0/2 - 0Δ + A*cos(N*pi*x)
        init_mullins2!(grid, TL, T_inf, 0., A, N, 0.6L0/2)
    elseif num.case == "Crystal"
        u .= -R*ones(ny, nx) + sqrt.(x.^2 + y.^2).*(ones(ny, nx) + A*cos.(N*atan.(x./(y + 1E-30*ones(ny, nx)))))
        init_franck!(grid, TL, R, T_inf, 0)
    elseif num.case == "3Crystals"
        x_c = 0.2
        y_c = 0.2
        for II in CartesianIndices(u)
            x_ = x[II]
            y_ = y[II]
            u1 = -R + sqrt((x_+x_c)^2 + (y_+y_c)^2)*(1 + A*cos(N*atan((x_+x_c)/(y_+y_c))))
            u2 = -R + sqrt((x_-x_c)^2 + (y_+y_c)^2)*(1 + A*cos(N*atan((x_-x_c)/(y_+y_c))))
            u3 = -R + sqrt((x_+x_c)^2 + (y_-y_c)^2)*(1 + A*cos(N*atan((x_+x_c)/(y_-y_c))))
            u4 = -R + sqrt((x_-x_c)^2 + (y_-y_c)^2)*(1 + A*cos(N*atan((x_-x_c)/(y_-y_c))))
            u5 = -R + sqrt((x_)^2 + (y_)^2)*(1 + A*cos(N*atan((x_)/(y_))))
            u[II] = min(u1,u2,u3)
        end
        TL .= T_inf;
    elseif num.case == "Nothing"
        u .= sqrt.(x.^2 + y.^2) .+ 1e-8
    elseif num.case == "Airfoil"
        @inbounds @threads for II in ind.all_indices
            d = 2L0
            count = 0
            @inbounds for (x1, y1, x2, y2) in zip(x_airfoil[1:end-1], y_airfoil[1:end-1], x_airfoil[2:end], y_airfoil[2:end])
                # compute distance from the cell center to the closest point in the solid
                den = (x2-x1)^2 + (y2-y1)^2
                t = -((x1-x[II])*(x2-x1) + (y1-y[II])*(y2-y1)) / den
                if t >= 0 && t <= 1.0
                    d1 = abs((x2-x1)*(y1-y[II]) - (y2-y1)*(x1-x[II])) / sqrt(den)
                    d = min(d1, d)
                else
                    d1 = sqrt((x1-x[II])^2 + (y1-y[II])^2)
                    d2 = sqrt((x2-x[II])^2 + (y2-y[II])^2)
                    d = min(d1, d2, d)
                end
                
                # compute the sign of the distance
                if ((y1 ≤ y[II] && y2 ≥ y[II]) || (y1 > y[II] && y2 < y[II])) && (x1 ≥ x[II] || x2 ≥ x[II])
                    count += 1
                end
            end

            @inbounds u[II] = iseven(count) ? d : -d
        end
        uL .= u_inf
        vL .= v_inf
    elseif num.case == "Jet"
        Δy = max(y...) - min(y...)
        u .= y .- min(y...) .- Δy ./ 2. .+ R .+ A.*cos.(N.*pi.*x)
        @inbounds for i = 0:(ny÷2-1)
            u[end-i,:] = u[i+1,:]
        end
    end

    return (Operators(SCUTT, SCUTp, SCUTu, SCUTv, SCUTDx, SCUTDy, SCUTCT, SCUTGxT, SCUTGyT, SCUTGxp, SCUTGyp, SCUTGxϕ, SCUTGyϕ, SCUTCu, SCUTCv, LTS, LpS, LuS, LvS, AS, BS, GxpS, GypS, GxϕS, GyϕS, DxuS, DyvS, ApS, AuS, AvS, CTS, GxTS, GyTS, ftcGxTS, ftcGyTS, CuS, CvS, E11, E12_x, E12_y, E22, utpS, vtpS),
            Operators(LCUTT, LCUTp, LCUTu, LCUTv, LCUTDx, LCUTDy, LCUTCT, LCUTGxT, LCUTGyT, LCUTGxp, LCUTGyp, LCUTGxϕ, LCUTGyϕ, LCUTCu, LCUTCv, LTL, LpL, LuL, LvL, AL, BL, GxpL, GypL, GxϕL, GyϕL, DxuL, DyvL, ApL, AuL, AvL, CTL, GxTL, GyTL, ftcGxTL, ftcGyTL, CuL, CvL, E11, E12_x, E12_y, E22, utpL, vtpL),
            Phase(TS, pS, ϕS, Gxm1S, Gym1S, uS, vS, ucorrS, vcorrS, DTS, DϕS, DuS, DvS, tmpS),
            Phase(TL, pL, ϕL, Gxm1L, Gym1L, uL, vL, ucorrL, vcorrL, DTL, DϕL, DuL, DvL, tmpL),
            Forward(Tall, usave, uusave, uvsave, TSsave, TLsave, Tsave, psave, ϕsave, Uxsave, Uysave, Uxcorrsave, Uycorrsave, Vsave, κsave, lengthsave, time, Cd, Cl))
end

function init_mullins!(grid, T, V, t, A, N, shift)
    @unpack x, y, nx, ny, ind = grid
    @unpack all_indices = ind
    
    @inbounds for II in all_indices
        ystar = shift + y[II] + A*sin(N*pi*x[II])
        if ystar > V*t
            T[II] = -1+exp(-V*(ystar-V*t))
        end
    end
end

function init_mullins2!(grid, T, V, t, A, N, shift)
    @unpack x, y, nx, ny, ind = grid
    @unpack all_indices = ind

    @inbounds for II in all_indices
        ystar = shift + y[II] + A*cos(N*pi*x[II])
        if ystar > V*t
            T[II] = -1+exp(-V*(ystar-V*t))
        end
    end
end

function init_franck!(grid, temp, R, T_inf, h)
    @unpack x, y, nx, ny, ind = grid
    @unpack all_indices = ind

    @inbounds for II in all_indices
        s = sqrt(x[II]^2 + y[II]^2)
        if s >= (R - h)
            temp[II] = T_inf*(1-(expint(0.25*s^2))/(expint(0.25*R^2)))
        end
    end
end

function init_franck!(grid, temp, R, T_inf, h, t)
    @unpack x, y, nx, ny, ind = grid
    @unpack all_indices = ind

    @inbounds for II in all_indices
        s = sqrt(x[II]^2 + y[II]^2)/√t
        if s >= (R - h)
            temp[II] = T_inf*(1-(expint(0.25*s^2))/(expint(0.25*R^2)))
        end
    end
end

function init_smooth(X, Y)
    u = similar(X)
    @. u = ((X-1)^2 + (Y-1)^2 + 0.1 )*(sqrt(X^2 + Y^2) - 1.0)
    return u
end
