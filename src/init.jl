function set_indices(n)
    u = zeros(n,n)
    C = vcat(collect(i + j for j = 0:n-3, i = n+2:n:((n)^2-(n)))[:])

    all_indices = CartesianIndices(u)
    inside = CartesianIndices(OffsetArray(u[2:end-1,2:end-1], CartesianIndex(2, 2):CartesianIndex(n-1, n-1)))

    LEFT = filter(x->x.I[2]==1, all_indices)
    LEFT2 = filter(x->x.I[2]==2, all_indices)
    BOTTOM = filter(x->x.I[1]==1, all_indices)
    BOTTOM2 = filter(x->x.I[1]==2, all_indices)
    RIGHT = filter(x->x.I[2]==n, all_indices)
    RIGHT2 = filter(x->x.I[2]==n-1, all_indices)
    TOP = filter(x->x.I[1]==n, all_indices)
    TOP2 = filter(x->x.I[1]==n-1, all_indices)

    return Indices(inside, (LEFT, RIGHT2), (RIGHT, LEFT2), (BOTTOM, TOP2), (TOP, BOTTOM2), (LEFT, LEFT2), (BOTTOM, BOTTOM2), (RIGHT, RIGHT2), (TOP, TOP2))
end

function init_fields(num::NumericalParameters, idx::NumericalParameters)
    @unpack N, T_inf, A, R, n, L0, Δ, shifted, X, Y, H, max_iterations, CFL = num

    SCUT = zeros(n^2)
    LCUT = zeros(n^2)

    SOL = zeros(n,n,11)
    LIQ = zeros(n,n,11)

    sol_projection = Array{Gradient{Float64}, 2}(undef, n, n)
    liq_projection = Array{Gradient{Float64}, 2}(undef, n, n)

    ii = collect(i for i = 1:n^2)
    iw = collect(i for i = n+1:n^2)
    is = collect(i for i = 2:n^2)
    iN = collect(i for i = 1:n^2-1)
    ie = collect(i for i = 1:n^2-n)

    II = vcat(ii,iw,is,iN,ie)

    jj = collect(i for i = 1:n^2)
    jw = collect(i for i = 1:n^2-n)
    js = collect(i for i = 1:n^2-1)
    jn = collect(i for i = 2:n^2)
    je = collect(i for i = n+1:n^2)

    JJ = vcat(jj,jw,js,jn,je)

    a = ones(length(jj))
    b = zeros(length(jw)+length(js)+length(jn)+length(je))

    AS = sparse(II,JJ,vcat(a,b))
    AL = sparse(II,JJ,vcat(a,b))
    BS = sparse(II,JJ,vcat(a,b))
    BL = sparse(II,JJ,vcat(a,b))

    LSA = sparse(II,JJ,vcat(a,b))
    LSB = sparse(II,JJ,vcat(a,b))

    iso = zeros(n, n)
    u = zeros(n, n)
    TS = zeros(n, n)
    TL = zeros(n, n)
    Tall = zeros(n, n)
    V = zeros(n, n)
    κ = zeros(n, n)

    usave = zeros(max_iterations+1, n, n)
    TSsave = zeros(max_iterations+1, n, n)
    TLsave = zeros(max_iterations+1, n, n)
    Tsave = zeros(max_iterations+1, n, n)
    Vsave = zeros(max_iterations+1, n, n)
    κsave = zeros(max_iterations+1, n, n)
    lengthsave = zeros(max_iterations+1)

    TL[:,:] = zeros(n, n)
    TS[:,:] = zeros(n, n)
    Tall[:,:] = zeros(n, n)
    if num.case == "Planar"
        u[:,:] = Y .+ shifted
    elseif num.case == "Sphere"
        u[:,:] = sqrt.((X .+ shifted).^ 2 + Y .^ 2) - (R) * ones(n, n);
        init_franck!(TL, R, T_inf, H, n, 0)
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
    end

    return TempArrays(SCUT, LCUT, LIQ, SOL, sol_projection, liq_projection, AS, AL, BS, BL, LSA, LSB), Forward(iso, u, TS, TL, Tall, V, κ, usave, TSsave, TLsave, Tsave, Vsave, κsave, lengthsave)
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
