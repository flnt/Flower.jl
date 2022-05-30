function set_indices(n)
    p = zeros(n, n)
    u = zeros(n, n+1)
    v = zeros(n+1, n)

    all_indices_p = CartesianIndices(p)
    all_indices_u = CartesianIndices(u)
    all_indices_v = CartesianIndices(v)
    inside_p = CartesianIndices(OffsetArray(p[2:end-1,2:end-1], CartesianIndex(2, 2):CartesianIndex(n-1, n-1)))
    inside_u = CartesianIndices(OffsetArray(u[2:end-1,2:end-1], CartesianIndex(2, 2):CartesianIndex(n-1, n)))
    inside_v = CartesianIndices(OffsetArray(v[2:end-1,2:end-1], CartesianIndex(2, 2):CartesianIndex(n, n-1)))

    LEFT_p = filter(x->x.I[2]==1, all_indices_p)
    LEFT2_p = filter(x->x.I[2]==2, all_indices_p)
    BOTTOM_p = filter(x->x.I[1]==1, all_indices_p)
    BOTTOM2_p = filter(x->x.I[1]==2, all_indices_p)
    RIGHT_p = filter(x->x.I[2]==n, all_indices_p)
    RIGHT2_p = filter(x->x.I[2]==n-1, all_indices_p)
    TOP_p = filter(x->x.I[1]==n, all_indices_p)
    TOP2_p = filter(x->x.I[1]==n-1, all_indices_p)

    LEFT_u = filter(x->x.I[2]==1, all_indices_u)
    LEFT2_u = filter(x->x.I[2]==2, all_indices_u)
    BOTTOM_u = filter(x->x.I[1]==1, all_indices_u)
    BOTTOM2_u = filter(x->x.I[1]==2, all_indices_u)
    RIGHT_u = filter(x->x.I[2]==n+1, all_indices_u)
    RIGHT2_u = filter(x->x.I[2]==n, all_indices_u)
    TOP_u = filter(x->x.I[1]==n, all_indices_u)
    TOP2_u = filter(x->x.I[1]==n-1, all_indices_u)

    LEFT_v = filter(x->x.I[2]==1, all_indices_v)
    LEFT2_v = filter(x->x.I[2]==2, all_indices_v)
    BOTTOM_v = filter(x->x.I[1]==1, all_indices_v)
    BOTTOM2_v = filter(x->x.I[1]==2, all_indices_v)
    RIGHT_v = filter(x->x.I[2]==n, all_indices_v)
    RIGHT2_v = filter(x->x.I[2]==n-1, all_indices_v)
    TOP_v = filter(x->x.I[1]==n+1, all_indices_v)
    TOP2_v = filter(x->x.I[1]==n, all_indices_v)

    return (Indices(all_indices_p, inside_p, (LEFT_p, RIGHT2_p), (RIGHT_p, LEFT2_p), (BOTTOM_p, TOP2_p), (TOP_p, BOTTOM2_p), (LEFT_p, LEFT2_p), (BOTTOM_p, BOTTOM2_p), (RIGHT_p, RIGHT2_p), (TOP_p, TOP2_p)),
            Indices(all_indices_u, inside_u, (LEFT_u, RIGHT2_u), (RIGHT_u, LEFT2_u), (BOTTOM_u, TOP2_u), (TOP_u, BOTTOM2_u), (LEFT_u, LEFT2_u), (BOTTOM_u, BOTTOM2_u), (RIGHT_u, RIGHT2_u), (TOP_u, TOP2_u)),
            Indices(all_indices_v, inside_v, (LEFT_v, RIGHT2_v), (RIGHT_v, LEFT2_v), (BOTTOM_v, TOP2_v), (TOP_v, BOTTOM2_v), (LEFT_v, LEFT2_v), (BOTTOM_v, BOTTOM2_v), (RIGHT_v, RIGHT2_v), (TOP_v, TOP2_v)))
end

function init_fields(num::NumericalParameters, idx::NumericalParameters, idxu::NumericalParameters, idxv::NumericalParameters)
    @unpack N, T_inf, u_inf, v_inf, A, R, n, L0, Δ, shifted, X, Y, Xu, Yu, Xv, Yv, H, max_iterations, save_every, CFL = num

    SCUTT = zeros(n^2)
    LCUTT = zeros(n^2)

    SCUTp = zeros(n^2)
    LCUTp = zeros(n^2)

    SCUTu = zeros(n*(n+1))
    LCUTu = zeros(n*(n+1))

    SCUTv = zeros(n*(n+1))
    LCUTv = zeros(n*(n+1))

    SCUTDx = zeros(n^2)
    SCUTDy = zeros(n^2)
    LCUTDx = zeros(n^2)
    LCUTDy = zeros(n^2)

    SCUTCT = zeros(n^2)
    LCUTCT = zeros(n^2)

    SCUTGxT = zeros(n*(n+1))
    LCUTGxT = zeros(n*(n+1))
    SCUTGyT = zeros(n*(n+1))
    LCUTGyT = zeros(n*(n+1))

    SCUTCu = zeros(n*(n+1))
    LCUTCu = zeros(n*(n+1))
    SCUTCv = zeros(n*(n+1))
    LCUTCv = zeros(n*(n+1))

    SOL = zeros(n,n,11)
    LIQ = ones(n,n,11)
    LIQ[:,:,8:end] .*= 0.5

    SOLu = zeros(n,n+1,11)
    LIQu = ones(n,n+1,11)
    LIQu[:,:,8:end] .*= 0.5

    SOLv = zeros(n+1,n,11)
    LIQv = ones(n+1,n,11)
    LIQv[:,:,8:end] .*= 0.5

    sol_projection = Array{Gradient{Float64}, 2}(undef, n, n)
    liq_projection = Array{Gradient{Float64}, 2}(undef, n, n)

    sol_projectionu = Array{Gradient{Float64}, 2}(undef, n, n+1)
    liq_projectionu = Array{Gradient{Float64}, 2}(undef, n, n+1)

    sol_projectionv = Array{Gradient{Float64}, 2}(undef, n+1, n)
    liq_projectionv = Array{Gradient{Float64}, 2}(undef, n+1, n)

    sol_centroid = Array{Point{Float64}, 2}(undef, n, n)
    liq_centroid = Array{Point{Float64}, 2}(undef, n, n)
    mid_point = Array{Point{Float64}, 2}(undef, n, n)

    sol_centroidu = Array{Point{Float64}, 2}(undef, n, n+1)
    liq_centroidu = Array{Point{Float64}, 2}(undef, n, n+1)
    mid_pointu = Array{Point{Float64}, 2}(undef, n, n+1)

    sol_centroidv = Array{Point{Float64}, 2}(undef, n+1, n)
    liq_centroidv = Array{Point{Float64}, 2}(undef, n+1, n)
    mid_pointv = Array{Point{Float64}, 2}(undef, n+1, n)

    @inbounds @threads for II in eachindex(sol_centroid)
        sol_centroid[II] = Point(0.0, 0.0)
        liq_centroid[II] = Point(0.0, 0.0)
        mid_point[II] = Point(0.0, 0.0)
    end
    @inbounds @threads for II in eachindex(sol_centroidu)
        sol_centroidu[II] = Point(0.0, 0.0)
        liq_centroidu[II] = Point(0.0, 0.0)
        mid_pointu[II] = Point(0.0, 0.0)
    end
    @inbounds @threads for II in eachindex(sol_centroidv)
        sol_centroidv[II] = Point(0.0, 0.0)
        liq_centroidv[II] = Point(0.0, 0.0)
        mid_pointv[II] = Point(0.0, 0.0)
    end

    α = zeros(n,n)
    α .= NaN

    αu = zeros(n,n+1)
    αu.= NaN

    αv = zeros(n+1,n)
    αv .= NaN

    ii = collect(i for i = 1:n^2)
    iw = collect(i for i = n+1:n^2)
    is = collect(i for i = 2:n^2)
    iN = collect(i for i = 1:n^2-1)
    ie = collect(i for i = 1:n^2-n)
    iwp = collect(i for i = 1:n)
    isp = collect(i for i = 1:n:n^2)
    inp = collect(i for i = n:n:n^2)
    iep = collect(i for i = n^2-n+1:n^2)

    II = vcat(ii,iw,is,iN,ie,iwp,isp,inp,iep)

    jj = collect(i for i = 1:n^2)
    jw = collect(i for i = 1:n^2-n)
    js = collect(i for i = 1:n^2-1)
    jn = collect(i for i = 2:n^2)
    je = collect(i for i = n+1:n^2)
    jwp = collect(i for i = n^2-n+1:n^2)
    jsp = collect(i for i = n:n:n^2)
    jnp = collect(i for i = 1:n:n^2)
    jep = collect(i for i = 1:n)

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

    LSA = sparse(II,JJ,vcat(a,b,c))
    LSB = sparse(II,JJ,vcat(a,b,c))

    ii = collect(i for i = 1:n*(n+1))
    iw = collect(i for i = n+1:n*(n+1))
    is = collect(i for i = 2:n*(n+1))
    iN = collect(i for i = 1:n*(n+1)-1)
    ie = collect(i for i = 1:n*(n+1)-n)
    iwp = collect(i for i = 1:n)
    isp = collect(i for i = 1:n:n*(n+1))
    inp = collect(i for i = n:n:n*(n+1))
    iep = collect(i for i = n*(n+1)-n+1:n*(n+1))

    II = vcat(ii,iw,is,iN,ie,iwp,isp,inp,iep)

    jj = collect(i for i = 1:n*(n+1))
    jw = collect(i for i = 1:n*(n+1)-n)
    js = collect(i for i = 1:n*(n+1)-1)
    jn = collect(i for i = 2:n*(n+1))
    je = collect(i for i = n+1:n*(n+1))
    jwp = collect(i for i = n*(n+1)-n+1:n*(n+1))
    jsp = collect(i for i = n:n:n*(n+1))
    jnp = collect(i for i = 1:n:n*(n+1))
    jep = collect(i for i = 1:n)

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

    ii = collect(i for i = 1:(n+1)*n)
    iw = collect(i for i = n+2:(n+1)*n)
    is = collect(i for i = 2:(n+1)*n)
    iN = collect(i for i = 1:(n+1)*n-1)
    ie = collect(i for i = 1:(n+1)*n-n-1)
    iwp = collect(i for i = 1:n+1)
    isp = collect(i for i = 1:n+1:(n+1)*n)
    inp = collect(i for i = n+1:n+1:(n+1)*n)
    iep = collect(i for i = (n+1)*n-n:(n+1)*n)

    II = vcat(ii,iw,is,iN,ie,iwp,isp,inp,iep)

    jj = collect(i for i = 1:(n+1)*n)
    jw = collect(i for i = 1:(n+1)*n-n-1)
    js = collect(i for i = 1:(n+1)*n-1)
    jn = collect(i for i = 2:(n+1)*n)
    je = collect(i for i = n+2:(n+1)*n)
    jwp = collect(i for i = (n+1)*n-n:(n+1)*n)
    jsp = collect(i for i = n+1:n+1:(n+1)*n)
    jnp = collect(i for i = 1:n+1:(n+1)*n)
    jep = collect(i for i = 1:n+1)

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

    iw = collect(i for i = n+1:n*(n+1))
    ie = collect(i for i = 1:n*n)
    iwp = collect(i for i = 1:n)
    iep = collect(i for i = n*(n+1)-n+1:n*(n+1))

    II = vcat(iw,ie,iwp,iep)

    jw = collect(i for i = 1:n*n)
    je = collect(i for i = 1:n*n)
    jwp = collect(i for i = n*n-n+1:n*n)
    jep = collect(i for i = 1:n)

    JJ = vcat(jw,je,jwp,jwp)

    a = zeros(length(jw)+length(je)+length(jwp)+length(jep))

    GxpS = sparse(II,JJ,a)
    GxpL = sparse(II,JJ,a)
    GxTS = sparse(II,JJ,a)
    GxTL = sparse(II,JJ,a)

    is = collect(i for i = 2:n+1)
    for j=2:n
        is = vcat(is, collect(i for i = (j-1)*(n+1)+2:(j-1)*(n+1)+n+1))
    end
    iN = collect(i for i = 1:n)
    for j=2:n
        iN = vcat(iN, collect(i for i = (j-1)*(n+1)+1:(j-1)*(n+1)+n))
    end
    isp = collect(i for i = 1:n+1:(n+1)*n)
    inp = collect(i for i = n+1:n+1:(n+1)*n)
    
    II = vcat(is,iN,isp,inp)

    js = collect(i for i = 1:n*n)
    jn = collect(i for i = 1:n*n)
    jsp = collect(i for i = n:n:n*n)
    jnp = collect(i for i = 1:n:n*n)

    JJ = vcat(js,jn,jsp,jnp)

    a = zeros(length(js)+length(jn)+length(jsp)+length(jnp))

    GypS = sparse(II,JJ,a)
    GypL = sparse(II,JJ,a)
    GyTS = sparse(II,JJ,a)
    GyTL = sparse(II,JJ,a)

    iw = collect(i for i = 1:n*n)
    ie = collect(i for i = 1:n*n)
    II = vcat(iw,ie)

    jw = collect(i for i = 1:n*n)
    je = collect(i for i = n+1:n*(n+1))
    JJ = vcat(jw,je)

    a = zeros(length(jw)+length(je))

    DxuS = sparse(II,JJ,a)
    DxuL = sparse(II,JJ,a)
    ftcGxTS = sparse(II,JJ,a)
    ftcGxTL = sparse(II,JJ,a)

    is = collect(i for i = 1:n*n)
    iN = collect(i for i = 1:n*n)
    II = vcat(is,iN)

    js = collect(i for i = 1:n)
    for j=2:n
        js = vcat(js, collect(i for i = (j-1)*(n+1)+1:(j-1)*(n+1)+n))
    end
    jn = collect(i for i = 2:n+1)
    for j=2:n
        jn = vcat(jn, collect(i for i = (j-1)*(n+1)+2:(j-1)*(n+1)+n+1))
    end
    JJ = vcat(js,jn)

    a = zeros(length(js)+length(jn))

    DyvS = sparse(II,JJ,a)
    DyvL = sparse(II,JJ,a)
    ftcGyTS = sparse(II,JJ,a)
    ftcGyTL = sparse(II,JJ,a)

    iso = zeros(n, n)
    isou = zeros(n, n+1)
    isov = zeros(n+1, n)
    u = zeros(n, n)
    uu = zeros(n, n+1)
    uv = zeros(n+1, n)
    TS = zeros(n, n)
    TL = zeros(n, n)
    pS = zeros(n, n)
    pL = zeros(n, n)
    ϕS = zeros(n, n)
    ϕL = zeros(n, n)
    uS = zeros(n, n+1)
    uL = zeros(n, n+1)
    vS = zeros(n+1, n)
    vL = zeros(n+1, n)
    Tall = zeros(n, n)
    DTS = zeros(n, n)
    DTL = zeros(n, n)
    V = zeros(n, n)
    Vu = zeros(n, n+1)
    Vv = zeros(n+1, n)
    κ = zeros(n, n)
    κu = zeros(n, n+1)
    κv = zeros(n+1, n)

    n_snaps = iszero(max_iterations%save_every) ? max_iterations÷save_every+1 : max_iterations÷save_every+2
    
    usave = zeros(n_snaps, n, n)
    uusave = zeros(n_snaps, n, n+1)
    uvsave = zeros(n_snaps, n+1, n)
    TSsave = zeros(n_snaps, n, n)
    TLsave = zeros(n_snaps, n, n)
    Tsave = zeros(n_snaps, n, n)
    psave = zeros(n_snaps, n, n)
    Uxsave = zeros(n_snaps, n, n+1)
    Uysave = zeros(n_snaps, n+1, n)
    Vsave = zeros(n_snaps, n, n)
    κsave = zeros(n_snaps, n, n)
    lengthsave = zeros(n_snaps)

    TL[:,:] = zeros(n, n)
    TS[:,:] = zeros(n, n)
    Tall[:,:] = zeros(n, n)
    if num.case == "Planar"
        u[:,:] = Y .+ shifted
    elseif num.case == "Sphere"
        u[:,:] = sqrt.((X .+ shifted).^ 2 + Y .^ 2) - (R) * ones(n, n);
        init_franck!(TL, R, T_inf, H, n, 0)

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
        u[:,:] = sqrt.((X .+ shifted).^ 2 + Y .^ 2) - (R) * ones(n, n);
        init_franck!(TL, R, T_inf, H, n, 0)

        su = sqrt.((Xu .+ shifted).^2 .+ Yu.^2)
        R1 = R + 3.0*Δ
        uL .= u_inf
        vL .= v_inf
        bl = 4.0
        for II in idxu.all_indices
            if su[II] <= R1
                uL[II] = 0.0
            elseif su[II] > R1
                uL[II] = tanh(bl*(su[II]-R1))
            end
        end
    elseif num.case == "Mullins"
        u[:,:] = Y.+ shifted .+ A*sin.(N*pi*X) .+ L0/4;
        init_mullins!(TL, T_inf, 0., A, N, n, H, L0/4)
    elseif num.case == "Mullins_cos"
        @. u = Y + 0.6L0/2 - 0Δ + A*cos(N*pi*X)
        init_mullins2!(TL, T_inf, 0., A, N, n, H, 0.6L0/2)
    elseif num.case == "Crystal"
        u[:,:] = -R*ones(n, n) + sqrt.(X.^2 + Y.^2).*(ones(n,n) + A*cos.(N*atan.(X./(Y + 1E-30*ones(n, n)))))
        init_franck!(TL, R, T_inf, H, n, 0)
    elseif num.case == "3Crystals"
        x_c = 0.2
        y_c = 0.2
        for II in CartesianIndices(u)
            x_ = X[II]
            y_ = Y[II]
            u1 = -R + sqrt((x_+x_c)^2 + (y_+y_c)^2)*(1 + A*cos(N*atan((x_+x_c)/(y_+y_c))))
            u2 = -R + sqrt((x_-x_c)^2 + (y_+y_c)^2)*(1 + A*cos(N*atan((x_-x_c)/(y_+y_c))))
            u3 = -R + sqrt((x_+x_c)^2 + (y_-y_c)^2)*(1 + A*cos(N*atan((x_+x_c)/(y_-y_c))))
            u4 = -R + sqrt((x_-x_c)^2 + (y_-y_c)^2)*(1 + A*cos(N*atan((x_-x_c)/(y_-y_c))))
            u5 = -R + sqrt((x_)^2 + (y_)^2)*(1 + A*cos(N*atan((x_)/(y_))))
            u[II] = min(u1,u2,u3)
        end
        TL .= T_inf;
    elseif num.case == "Nothing"
        u .= sqrt.(X.^2 + Y.^2) .+ 1e-8
    end

    return TempArrays(SCUTT, LCUTT, SCUTp, LCUTp, SCUTu, LCUTu, SCUTv, LCUTv, SCUTDx, SCUTDy, LCUTDx, LCUTDy, SCUTCT, LCUTCT, SCUTGxT, LCUTGxT, SCUTGyT, LCUTGyT, SCUTCu, LCUTCu, SCUTCv, LCUTCv, SOL, LIQ, SOLu, LIQu, SOLv, LIQv, sol_projection, liq_projection, sol_projectionu, liq_projectionu, sol_projectionv, liq_projectionv, sol_centroid, liq_centroid, mid_point, sol_centroidu, liq_centroidu, mid_pointu, sol_centroidv, liq_centroidv, mid_pointv, α, αu, αv, LTS, LTL, LpS, LpL, LuS, LuL, LvS, LvL, AS, AL, BS, BL, LSA, LSB, GxpS, GxpL, GypS, GypL, DxuS, DxuL, DyvS, DyvL, ApS, ApL, AuS, AuL, AvS, AvL, CTS, CTL, GxTS, GxTL, GyTS, GyTL, ftcGxTS, ftcGxTL, ftcGyTS, ftcGyTL, CuS, CuL, CvS, CvL), Forward(iso, isou, isov, u, uu, uv, TS, TL, pS, pL, ϕS, ϕL, uS, uL, vS, vL, Tall, DTS, DTL, V, Vu, Vv, κ, κu, κv, usave, uusave, uvsave, TSsave, TLsave, Tsave, psave, Uxsave, Uysave, Vsave, κsave, lengthsave)
end

function init_mullins!(T, V, t, A, N, n, H, shift)
    @inbounds for j = 1:n
        @inbounds for i = 1:n
            ystar = shift + H[i] + A*sin(N*pi*(H[j]))
            if ystar > V*t
                T[i, j] = -1+exp(-V*(ystar-V*t))
            end
        end
    end
end

function init_mullins2!(T, V, t, A, N, n, H, shift)
    @inbounds for j = 1:n
        @inbounds for i = 1:n
            ystar = shift + H[i] + A*cos(N*pi*(H[j]))
            if ystar > V*t
                T[i, j] = -1+exp(-V*(ystar-V*t))
            end
        end
    end
end

function init_franck!(temp, R, T_inf, H, n, h)
    @inbounds for j = 1:n
        @inbounds for i = 1:n
            s = sqrt(H[j]^2+H[i]^2)
            if s >= (R - h)
                temp[i, j] = T_inf*(1-(expint(0.25*s^2))/(expint(0.25*R^2)))
            end
        end
    end
end

function init_smooth(X, Y)
    u = similar(X)
    @. u = ((X-1)^2 + (Y-1)^2 + 0.1 )*(sqrt(X^2 + Y^2) - 1.0)
    return u
end
