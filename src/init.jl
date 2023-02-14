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

    return Indices(all_indices, inside, (LEFT, RIGHT), (BOTTOM, TOP),
                   (LEFT, LEFT2), (BOTTOM, BOTTOM2), (RIGHT, RIGHT2), (TOP, TOP2))
end

function Mesh(grid, x_nodes, y_nodes, s, o)
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

    freshS = zeros(Bool, ny, nx)
    freshL = zeros(Bool, ny, nx)

    geoS = GeometricInfo(SOL, dSOL, sol_projection, sol_centroid, emptiedS, freshS)
    geoL = GeometricInfo(LIQ, dLIQ, liq_projection, liq_centroid, emptiedL, freshL)

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

    dom = domain((OneTo(ny), OneTo(nx), OneTo(2)))
    subs = (s, s, 1)
    over = (o, o, 0)
    dec = DDM.decompose(dom, subs, over)

    domdec = decomposition(dom, dec)

    pou = uniform(domdec)

    return Mesh{grid,Float64,Int64}(x_nodes, y_nodes, x, y, nx, ny, dx, dy, ind, u, iso, faces, 
                geoS, geoL, mid_point, cut_points, α, κ, V, LSA, LSB, domdec, pou)
end

function init_meshes(num::NumericalParameters)
    mesh_cc = Mesh(GridCC, num.x, num.y, num.subdomains, num.overlaps)

    xx = vcat(num.x[1] - mesh_cc.dx[1,1]/2, num.x[1:end-1] .+ mesh_cc.dx[1,:]/2, num.x[end] + mesh_cc.dx[1,end]/2)
    mesh_stx = Mesh(GridFCx, xx, num.y, num.subdomains, num.overlaps)

    yy = vcat(num.y[1] - mesh_cc.dy[1,1]/2, num.y[1:end-1] .+ mesh_cc.dy[:,1]/2, num.y[end] + mesh_cc.dy[end,1]/2)
    mesh_sty = Mesh(GridFCy, num.x, yy, num.subdomains, num.overlaps)

    return (mesh_cc, mesh_stx, mesh_sty)
end

function init_sparse_Bx(grid)
    @unpack nx, ny = grid
    iw = collect(i for i = ny+1:(nx+1)*ny)
    ie = collect(i for i = 1:nx*ny)
    iwp = collect(i for i = 1:ny)
    iep = collect(i for i = nx*ny+1:(nx+1)*ny)

    II = vcat(iw,ie,iwp,iep)

    jw = collect(i for i = 1:nx*ny)
    je = collect(i for i = 1:nx*ny)
    jwp = collect(i for i = (nx-1)*ny+1:(nx)*ny)
    jep = collect(i for i = 1:ny)

    JJ = vcat(jw,je,jwp,jep)

    a = zeros(length(jw)+length(je)+length(jwp)+length(jep))

    Bx = sparse(II,JJ,a)
end

function init_sparse_By(grid)
    @unpack nx, ny = grid

    is = collect(i for i = 2:(ny+1))
    for j=2:nx
        is = vcat(is, collect(i for i = (j-1)*(ny+1)+2:j*(ny+1)))
    end
    iN = collect(i for i = 1:ny)
    for j=2:nx
        iN = vcat(iN, collect(i for i = (j-1)*(ny+1)+1:(j-1)*(ny+1)+ny))
    end
    isp = collect(i for i = 1:(ny+1):nx*(ny+1))
    inp = collect(i for i = (ny+1):(ny+1):nx*(ny+1))
    
    II = vcat(is,iN,isp,inp)

    js = collect(i for i = 1:nx*ny)
    jn = collect(i for i = 1:nx*ny)
    jsp = collect(i for i = ny:ny:nx*ny)
    jnp = collect(i for i = 1:ny:nx*ny)

    JJ = vcat(js,jn,jsp,jnp)

    a = zeros(length(js)+length(jn)+length(jsp)+length(jnp))

    By = sparse(II,JJ,a)
end

function init_sparse_BxT(grid)
    @unpack nx, ny = grid

    iw = collect(i for i = 1:nx*ny)
    ie = collect(i for i = 1:nx*ny)
    II = vcat(iw,ie)

    jw = collect(i for i = 1:nx*ny)
    je = collect(i for i = ny+1:(nx+1)*ny)
    JJ = vcat(jw,je)

    a = zeros(length(jw)+length(je))

    BxT = sparse(II,JJ,a)    
end

function init_sparse_ByT(grid)
    @unpack nx, ny = grid

    is = collect(i for i = 1:nx*ny)
    iN = collect(i for i = 1:nx*ny)
    II = vcat(is,iN)

    js = collect(i for i = 1:ny)
    for j=2:nx
        js = vcat(js, collect(i for i = (j-1)*(ny+1)+1:j*(ny+1)-1))
    end
    jn = collect(i for i = 2:(ny+1))
    for j=2:nx
        jn = vcat(jn, collect(i for i = (j-1)*(ny+1)+2:j*(ny+1)))
    end
    JJ = vcat(js,jn)

    a = zeros(length(js)+length(jn))

    ByT = sparse(II,JJ,a)
end

function init_fields(num::NumericalParameters, grid, grid_u, grid_v)
    @unpack τ, N, T_inf, u_inf, v_inf, A, R, L0, Δ, shifted, shift_y, max_iterations, save_every, CFL, x_airfoil, y_airfoil = num
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
    GxpS = init_sparse_Bx(grid)
    GxpL = init_sparse_Bx(grid)
    GxTS = init_sparse_Bx(grid)
    GxTL = init_sparse_Bx(grid)
    GxϕS = init_sparse_Bx(grid)
    GxϕL = init_sparse_Bx(grid)

    # Coupled system operators
    Bx_TS = init_sparse_Bx(grid)
    Bx_TL = init_sparse_Bx(grid)
    Hx_TS = init_sparse_Bx(grid)
    Hx_TL = init_sparse_Bx(grid)

    Bx_pS = init_sparse_Bx(grid)
    Bx_pL = init_sparse_Bx(grid)
    Hx_pS = init_sparse_Bx(grid)
    Hx_pL = init_sparse_Bx(grid)
    Gx_S = init_sparse_Bx(grid)
    Gx_L = init_sparse_Bx(grid)

    Bx_uS = init_sparse_Bx(grid_u)
    Bx_uL = init_sparse_Bx(grid_u)
    Hx_uS = init_sparse_Bx(grid_u)
    Hx_uL = init_sparse_Bx(grid_u)

    Bx_vS = init_sparse_Bx(grid_v)
    Bx_vL = init_sparse_Bx(grid_v)
    Hx_vS = init_sparse_Bx(grid_v)
    Hx_vL = init_sparse_Bx(grid_v)

    # 2 points stencil (p grid to v grid)
    GypS = init_sparse_By(grid)
    GypL = init_sparse_By(grid)
    GyTS = init_sparse_By(grid)
    GyTL = init_sparse_By(grid)
    GyϕS = init_sparse_By(grid)
    GyϕL = init_sparse_By(grid)

    # Coupled system operators
    By_TS = init_sparse_By(grid)
    By_TL = init_sparse_By(grid)
    Hy_TS = init_sparse_By(grid)
    Hy_TL = init_sparse_By(grid)

    By_pS = init_sparse_By(grid)
    By_pL = init_sparse_By(grid)
    Hy_pS = init_sparse_By(grid)
    Hy_pL = init_sparse_By(grid)
    Gy_S = init_sparse_By(grid)
    Gy_L = init_sparse_By(grid)

    By_uS = init_sparse_By(grid_u)
    By_uL = init_sparse_By(grid_u)
    Hy_uS = init_sparse_By(grid_u)
    Hy_uL = init_sparse_By(grid_u)

    By_vS = init_sparse_By(grid_v)
    By_vL = init_sparse_By(grid_v)
    Hy_vS = init_sparse_By(grid_v)
    Hy_vL = init_sparse_By(grid_v)

    # 2 points stencil (u grid to p grid)
    DxuS = init_sparse_BxT(grid)
    DxuL = init_sparse_BxT(grid)
    ftcGxTS = init_sparse_BxT(grid)
    ftcGxTL = init_sparse_BxT(grid)
    E11 = init_sparse_BxT(grid)
    utpS = init_sparse_BxT(grid)
    utpL = init_sparse_BxT(grid)

    # Coupled system operators
    AxT_TS = init_sparse_BxT(grid)
    AxT_TL = init_sparse_BxT(grid)
    BxT_TS = init_sparse_BxT(grid)
    BxT_TL = init_sparse_BxT(grid)
    HxT_TS = init_sparse_BxT(grid)
    HxT_TL = init_sparse_BxT(grid)

    AxT_pS = init_sparse_BxT(grid)
    AxT_pL = init_sparse_BxT(grid)
    BxT_pS = init_sparse_BxT(grid)
    BxT_pL = init_sparse_BxT(grid)
    HxT_pS = init_sparse_BxT(grid)
    HxT_pL = init_sparse_BxT(grid)
    GxT_S = init_sparse_BxT(grid)
    GxT_L = init_sparse_BxT(grid)

    AxT_uS = init_sparse_BxT(grid_u)
    AxT_uL = init_sparse_BxT(grid_u)
    BxT_uS = init_sparse_BxT(grid_u)
    BxT_uL = init_sparse_BxT(grid_u)
    HxT_uS = init_sparse_BxT(grid_u)
    HxT_uL = init_sparse_BxT(grid_u)

    AxT_vS = init_sparse_BxT(grid_v)
    AxT_vL = init_sparse_BxT(grid_v)
    BxT_vS = init_sparse_BxT(grid_v)
    BxT_vL = init_sparse_BxT(grid_v)
    HxT_vS = init_sparse_BxT(grid_v)
    HxT_vL = init_sparse_BxT(grid_v)

    # 2 points stencil (v grid to p grid)
    DyvS = init_sparse_ByT(grid)
    DyvL = init_sparse_ByT(grid)
    ftcGyTS = init_sparse_ByT(grid)
    ftcGyTL = init_sparse_ByT(grid)
    E22 = init_sparse_ByT(grid)
    vtpS = init_sparse_ByT(grid)
    vtpL = init_sparse_ByT(grid)

    # Coupled system operators
    AyT_TS = init_sparse_ByT(grid)
    AyT_TL = init_sparse_ByT(grid)
    ByT_TS = init_sparse_ByT(grid)
    ByT_TL = init_sparse_ByT(grid)
    HyT_TS = init_sparse_ByT(grid)
    HyT_TL = init_sparse_ByT(grid)

    AyT_pS = init_sparse_ByT(grid)
    AyT_pL = init_sparse_ByT(grid)
    ByT_pS = init_sparse_ByT(grid)
    ByT_pL = init_sparse_ByT(grid)
    HyT_pS = init_sparse_ByT(grid)
    HyT_pL = init_sparse_ByT(grid)
    GyT_S = init_sparse_ByT(grid)
    GyT_L = init_sparse_ByT(grid)

    AyT_uS = init_sparse_ByT(grid_u)
    AyT_uL = init_sparse_ByT(grid_u)
    ByT_uS = init_sparse_ByT(grid_u)
    ByT_uL = init_sparse_ByT(grid_u)
    HyT_uS = init_sparse_ByT(grid_u)
    HyT_uL = init_sparse_ByT(grid_u)

    AyT_vS = init_sparse_ByT(grid_v)
    AyT_vL = init_sparse_ByT(grid_v)
    ByT_vS = init_sparse_ByT(grid_v)
    ByT_vL = init_sparse_ByT(grid_v)
    HyT_vS = init_sparse_ByT(grid_v)
    HyT_vL = init_sparse_ByT(grid_v)

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

    # Auxiliar matrices from cell centered to 2-staggered
    ii = collect(i for i = ny+1:ny+nx*ny)
    ie = collect(i for i = (nx+1)*ny+1:(nx+2)*ny)
    iwp = collect(i for i = 1:ny)
    iep = collect(i for i = (nx+1)*ny+1:(nx+2)*ny)

    jj = collect(i for i = 1:nx*ny)
    je = collect(i for i = (nx-1)*ny+1:nx*ny)
    jwp = collect(i for i = (nx-1)*ny+1:nx*ny)
    jep = collect(i for i = 1:ny)

    II = vcat(ii,ie,iwp,iep)
    JJ = vcat(jj,je,jwp,jep)
    a = ones(length(jj))
    b = zeros(length(je)+length(jwp)+length(jep))
    Rx = sparse(II,JJ,vcat(a,b))

    ii = []
    for j=1:nx
        ii = vcat(ii, collect(i for i = 2+(j-1)*(ny+2):ny+1+(j-1)*(ny+2)))
    end
    is = []
    for j=1:nx
        is = vcat(is, 1+(j-1)*(ny+2))
    end
    iN = []
    for j=1:nx
        iN = vcat(iN, ny+2+(j-1)*(ny+2))
    end
    isp = collect(i for i = 1:ny+2:(nx-1)*(ny+2)+1)
    inp = collect(i for i = ny+2:ny+2:nx*(ny+2))

    jj = []
    for j=1:nx
        jj = vcat(jj, collect(i for i = 1+(j-1)*ny:ny+(j-1)*ny))
    end
    js = []
    for j=1:nx
        js = vcat(js, 1+(j-1)*ny)
    end
    jn = []
    for j=1:nx
        jn = vcat(jn, ny+(j-1)*ny)
    end
    jsp = collect(i for i = ny:ny:nx*ny)
    jnp = collect(i for i = 1:ny:(nx-1)*ny+1)

    II = vcat(ii,iN,isp,inp)
    JJ = vcat(jj,jn,jsp,jnp)
    a = ones(length(jj))
    b = zeros(length(jn)+length(jsp)+length(jnp))
    Ry = sparse(II,JJ,vcat(a,b))

    tmp_x_TS = copy(Bx_TS)
    tmp_x_TL = copy(Bx_TL)
    tmp_x_pS = copy(Bx_pS)
    tmp_x_pL = copy(Bx_pL)
    tmp_x_uS = copy(Bx_uS)
    tmp_x_uL = copy(Bx_uL)
    tmp_x_vS = copy(Bx_vS)
    tmp_x_vL = copy(Bx_vL)

    tmp_y_TS = copy(By_TS)
    tmp_y_TL = copy(By_TL)
    tmp_y_pS = copy(By_pS)
    tmp_y_pL = copy(By_pL)
    tmp_y_uS = copy(By_uS)
    tmp_y_uL = copy(By_uL)
    tmp_y_vS = copy(By_vS)
    tmp_y_vL = copy(By_vL)

    M_TS = Diagonal(zeros(nx*ny))
    M_TL = Diagonal(zeros(nx*ny))
    iMx_TS = Diagonal(zeros((nx+1)*ny))
    iMx_TL = Diagonal(zeros((nx+1)*ny))
    iMy_TS = Diagonal(zeros(nx*(ny+1)))
    iMy_TL = Diagonal(zeros(nx*(ny+1)))
    χ_TS = Diagonal(zeros(nx*ny))
    χ_TL = Diagonal(zeros(nx*ny))

    M_pS = Diagonal(zeros(nx*ny))
    M_pL = Diagonal(zeros(nx*ny))
    iMx_pS = Diagonal(zeros((nx+1)*ny))
    iMx_pL = Diagonal(zeros((nx+1)*ny))
    iMy_pS = Diagonal(zeros(nx*(ny+1)))
    iMy_pL = Diagonal(zeros(nx*(ny+1)))
    χ_pS = Diagonal(zeros(nx*ny))
    χ_pL = Diagonal(zeros(nx*ny))

    M_uS = Diagonal(zeros(grid_u.nx*grid_u.ny))
    M_uL = Diagonal(zeros(grid_u.nx*grid_u.ny))
    iMx_uS = Diagonal(zeros((grid_u.nx+1)*grid_u.ny))
    iMx_uL = Diagonal(zeros((grid_u.nx+1)*grid_u.ny))
    iMy_uS = Diagonal(zeros(grid_u.nx*(grid_u.ny+1)))
    iMy_uL = Diagonal(zeros(grid_u.nx*(grid_u.ny+1)))
    χ_uS = Diagonal(zeros(grid_u.nx*grid_u.ny))
    χ_uL = Diagonal(zeros(grid_u.nx*grid_u.ny))

    M_vS = Diagonal(zeros(grid_v.nx*grid_v.ny))
    M_vL = Diagonal(zeros(grid_v.nx*grid_v.ny))
    iMx_vS = Diagonal(zeros((grid_v.nx+1)*grid_v.ny))
    iMx_vL = Diagonal(zeros((grid_v.nx+1)*grid_v.ny))
    iMy_vS = Diagonal(zeros(grid_v.nx*(grid_v.ny+1)))
    iMy_vL = Diagonal(zeros(grid_v.nx*(grid_v.ny+1)))
    χ_vS = Diagonal(zeros(grid_v.nx*grid_v.ny))
    χ_vL = Diagonal(zeros(grid_v.nx*grid_v.ny))

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
    TDS = zeros(2*ny*nx)
    TDL = zeros(2*ny*nx)
    pDS = zeros(2*ny*nx)
    pDL = zeros(2*ny*nx)
    ϕDS = zeros(2*ny*nx)
    ϕDL = zeros(2*ny*nx)
    uDS = zeros(2*grid_u.ny*grid_u.nx)
    uDL = zeros(2*grid_u.ny*grid_u.nx)
    vDS = zeros(2*grid_v.ny*grid_v.nx)
    vDL = zeros(2*grid_v.ny*grid_v.nx)
    uvDS = zeros(2*grid_u.ny*grid_u.nx + 2*grid_v.ny*grid_v.nx)
    uvDL = zeros(2*grid_u.ny*grid_u.nx + 2*grid_v.ny*grid_v.nx)
    uvϕDS = zeros(2*grid_u.ny*grid_u.nx + 2*grid_v.ny*grid_v.nx + 2*grid.ny*grid.nx)
    uvϕDL = zeros(2*grid_u.ny*grid_u.nx + 2*grid_v.ny*grid_v.nx + 2*grid.ny*grid.nx)

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
        u .= sqrt.((x .+ shifted).^ 2 + (y .+ shift_y) .^ 2) - (R) * ones(ny, nx);
        init_franck!(grid, TL, R, T_inf, 0)

        su = sqrt.((grid_u.x .+ shifted).^2 .+ (grid_u.y .+ shift_y).^2)
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
        @. u = y + shifted + 0.6L0/2 - 0Δ + A*cos(N*pi*x)
        init_mullins2!(grid, TL, T_inf, 0., A, N, 0.6L0/2)
    elseif num.case == "Drop"
        maxy = max(y...)
        @. u = y - (maxy - R) + R*A*cos(2*pi*x/(x[end]-x[1])) + shifted
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
            OperatorsCoupled(AxT_TS, AyT_TS, Bx_TS, By_TS, BxT_TS, ByT_TS, Hx_TS, Hy_TS, HxT_TS, HyT_TS, tmp_x_TS, tmp_y_TS, M_TS, iMx_TS, iMy_TS, χ_TS, Rx, Ry, GxT_S, GyT_S),
            OperatorsCoupled(AxT_TL, AyT_TL, Bx_TL, By_TL, BxT_TL, ByT_TL, Hx_TL, Hy_TL, HxT_TL, HyT_TL, tmp_x_TL, tmp_y_TL, M_TL, iMx_TL, iMy_TL, χ_TL, Rx, Ry, GxT_L, GyT_L),
            OperatorsCoupled(AxT_pS, AyT_pS, Bx_pS, By_pS, BxT_pS, ByT_pS, Hx_pS, Hy_pS, HxT_pS, HyT_pS, tmp_x_pS, tmp_y_pS, M_pS, iMx_pS, iMy_pS, χ_pS, Rx, Ry, GxT_S, GyT_S),
            OperatorsCoupled(AxT_pL, AyT_pL, Bx_pL, By_pL, BxT_pL, ByT_pL, Hx_pL, Hy_pL, HxT_pL, HyT_pL, tmp_x_pL, tmp_y_pL, M_pL, iMx_pL, iMy_pL, χ_pL, Rx, Ry, GxT_L, GyT_L),
            OperatorsCoupled(AxT_uS, AyT_uS, Bx_uS, By_uS, BxT_uS, ByT_uS, Hx_uS, Hy_uS, HxT_uS, HyT_uS, tmp_x_uS, tmp_y_uS, M_uS, iMx_uS, iMy_uS, χ_uS, Rx, Ry, Gx_S, Gy_S),
            OperatorsCoupled(AxT_uL, AyT_uL, Bx_uL, By_uL, BxT_uL, ByT_uL, Hx_uL, Hy_uL, HxT_uL, HyT_uL, tmp_x_uL, tmp_y_uL, M_uL, iMx_uL, iMy_uL, χ_uL, Rx, Ry, Gx_L, Gy_L),
            OperatorsCoupled(AxT_vS, AyT_vS, Bx_vS, By_vS, BxT_vS, ByT_vS, Hx_vS, Hy_vS, HxT_vS, HyT_vS, tmp_x_vS, tmp_y_vS, M_vS, iMx_vS, iMy_vS, χ_vS, Rx, Ry, Gx_S, Gy_S),
            OperatorsCoupled(AxT_vL, AyT_vL, Bx_vL, By_vL, BxT_vL, ByT_vL, Hx_vL, Hy_vL, HxT_vL, HyT_vL, tmp_x_vL, tmp_y_vL, M_vL, iMx_vL, iMy_vL, χ_vL, Rx, Ry, Gx_L, Gy_L),
            Phase(TS, pS, ϕS, Gxm1S, Gym1S, uS, vS, ucorrS, vcorrS, DTS, DϕS, DuS, DvS, TDS, pDS, ϕDS, uDS, vDS, uvDS, uvϕDS),
            Phase(TL, pL, ϕL, Gxm1L, Gym1L, uL, vL, ucorrL, vcorrL, DTL, DϕL, DuL, DvL, TDL, pDL, ϕDL, uDL, vDL, uvDL, uvϕDL),
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

function init_smooth(X, Y)
    u = similar(X)
    @. u = ((X-1)^2 + (Y-1)^2 + 0.1 )*(sqrt(X^2 + Y^2) - 1.0)
    return u
end
