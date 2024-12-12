

"""
    computes the cell-averaged gradient (grad_x, grad_y) for scalar TD (stored in 1D)
* grad_x : the useful information is of size (nx-1,ny) in an array of size (nx+1,ny) : grad_x[:,2:end-1]
* grad_y : the useful information is of size (nx,ny-1) in an array of size (nx,ny+1) : grad_y[2:end-1,:]
"""
function compute_grad_T_x_T_y_array!(num_LS, grid, grid_u, grid_v, opC_p, grad_x, grad_y, TD)

    ∇ϕ_x = opC_p.iMx * opC_p.Bx * vec1(TD,grid) .+ opC_p.iMx_b * opC_p.Hx_b * vecb(TD,grid)
    ∇ϕ_y = opC_p.iMy * opC_p.By * vec1(TD,grid) .+ opC_p.iMy_b * opC_p.Hy_b * vecb(TD,grid)

    for iLS in 1:num_LS
        ∇ϕ_x .+= opC_p.iMx * opC_p.Hx[iLS] * veci(TD,grid,iLS+1)
        ∇ϕ_y .+= opC_p.iMy * opC_p.Hy[iLS] * veci(TD,grid,iLS+1)
    end

    grad_x .= reshape(veci(∇ϕ_x,grid_u,1), grid_u)
    grad_y .= reshape(veci(∇ϕ_y,grid_v,1), grid_v)

end


"""
    computes the x component of the cell-averaged gradient: grad_x, for scalar TD (stored in 1D)
* grad_x : the useful information is of size (nx-1,ny) in an array of size (nx+1,ny) : grad_x[:,2:end-1]
"""
function compute_grad_T_x_array!(num_LS, grid, grid_u, opC_p, grad_x, TD)
    
    ∇ϕ_x = opC_p.iMx * opC_p.Bx * vec1(TD,grid) .+ opC_p.iMx_b * opC_p.Hx_b * vecb(TD,grid)

    for iLS in 1:num_LS
        ∇ϕ_x .+= opC_p.iMx * opC_p.Hx[iLS] * veci(TD,grid,iLS+1)
    end

    grad_x .= reshape(veci(∇ϕ_x,grid_u,1), grid_u)

end


"""
    computes the y component of the cell-averaged gradient: grad_y) for scalar TD (stored in 1D)
* grad_y : the useful information is of size (nx,ny-1) in an array of size (nx,ny+1) : grad_y[2:end-1,:]
"""
function compute_grad_T_y_array!(num_LS, grid, grid_v, opC_p, grad_y, TD)
    
    ∇ϕ_y = opC_p.iMy * opC_p.By * vec1(TD,grid) .+ opC_p.iMy_b * opC_p.Hy_b * vecb(TD,grid)

    for iLS in 1:num_LS
        ∇ϕ_y .+= opC_p.iMy * opC_p.Hy[iLS] * veci(TD,grid,iLS+1)
    end

    grad_y .= reshape(veci(∇ϕ_y,grid_v,1), grid_v)

end



"""
    Integrate mass transfer rate along a line (the interface in 2D), in  in mol/m/s:

    interface_id: id of bubble interface we want to integrate
    
In Flower, it seems we can identify a contribution of the divergence coming from the bubble (all terms beginning with:

```math
H^{\\Gamma 1 , T}
```

, but if I understood correctly we can say that even though we can identify a contribution 
from the bubble in the divergence, the gradient is a cell-averaged gradient 
and we cannot split it to keep only the derivative along the normal. 
As opposed to the term 

```math
\\nabla C \\cdot n
```

In one case the gradient includes the effect from the whole cell 
(and we have to divide the part of the integrated divergence by the interface length).
In the other case cf  

[`Johansen's method`](@ref interpolated_temperature),

, the gradient is a derivative along the normal. 
This is the cell-averaged gradient dot the normal integrated along the interface 

to increment mol/m in 2D : use dn=result*dt
to advect: need to go back to surface mass flux so divide by interface length : mol/m^2/s
then velocity m/s is obtained  (*M (kg/mol) and * combination of 1/rho_i ))

need to pay attention to operators of levelsets describing walls (do not compute mass flux at interface outside domain)

with iLS=1 for example for first levelset

```julia
# Interface
opC_p.HxT[iLS] (nx*ny, (nx+1)*ny) TODO (nx+1)*ny or nx*(ny+1)
# Wall
opC_p.iMx_b ((nx+1)*ny, 2*nx+2*ny)
opC_p.Hx_b (2*nx+2*ny, 2*nx+2*ny)
```

computes also in wall described  by LS: you need to select the contribution based on the levelset afterwards
Interface described by LS number one
```julia
opC_p.Hx_b: left: -dy (cell height along y), bottom: 0, right: +dy, top: 0
opC_p.Hy_b: left: 0 , bottom: -dx, right: 0, top: +dx
```

"""
function integrate_mass_flux_over_interface(num::Numerical{Float64, Int64},
    grid::Mesh{Flower.GridCC, Float64, Int64},
    opC_pL::Operators{Float64, Int64}, 
    scalD::AbstractArray{Float64, 1},
    mass_flux_vec1::Array{Float64, 1},
    mass_flux_vecb::Array{Float64, 1}, 
    mass_flux_veci::Array{Float64, 1},
    mass_flux_vec1_2::Array{Float64, 2},
    mass_flux_vecb_2::Array{Float64, 2},
    mass_flux_veci_2::Array{Float64, 2},
    mass_flux::Array{Float64, 2},
    interface_id::Int64,
    )

    opC_p = opC_pL

    # Interface described by LS number one

    #TODO opC_p.HxT[iLStmp] everywhere

    # for iLS in 1:num.nLS
    #     mass_flux_veci .+= opC_p.HxT[iLStmp] * opC_p.iMx * opC_p.Hx[iLS] * veci(scalD,grid,iLS+1)
    #     mass_flux_veci .+= opC_p.HyT[iLStmp] * opC_p.iMy * opC_p.Hy[iLS] * veci(scalD,grid,iLS+1)
    # end

    #size (nx*ny)
    mass_flux_vec1 .= 0.0 
    mass_flux_vecb .= 0.0
    mass_flux_veci .= 0.0
    
    #size (ny,nx)
    mass_flux .= 0.0
    mass_flux_vec1_2 .= 0.0
    mass_flux_vecb_2 .= 0.0
    mass_flux_veci_2 .= 0.0


    mass_flux_vec1   .= opC_p.HxT[interface_id] * opC_p.iMx * opC_p.Bx * vec1(scalD,grid) .+ opC_p.HyT[interface_id] * opC_p.iMy * opC_p.By * vec1(scalD,grid)
    mass_flux_vecb   .= opC_p.HxT[interface_id] * opC_p.iMx_b * opC_p.Hx_b * vecb(scalD,grid) .+ opC_p.HyT[interface_id] *  opC_p.iMy_b * opC_p.Hy_b * vecb(scalD,grid)

    for iLS in 1:num.nLS
        mass_flux_veci .+= opC_p.HxT[interface_id] * opC_p.iMx * opC_p.Hx[iLS] * veci(scalD,grid,iLS+1)
        mass_flux_veci .+= opC_p.HyT[interface_id] * opC_p.iMy * opC_p.Hy[iLS] * veci(scalD,grid,iLS+1)
    end

    # printstyled(color=:red, @sprintf "\n vec1 x y %.2e %.2e \n" sum(opC_p.HxT[iLStmp] * opC_p.iMx * opC_p.Bx * vec1(scalD,grid)) sum(opC_p.HyT[iLStmp] * opC_p.iMy * opC_p.By * vec1(scalD,grid)))
    # printstyled(color=:red, @sprintf "\n vecb x y %.2e %.2e\n" sum(opC_p.HxT[iLStmp] * opC_p.iMx_b * opC_p.Hx_b * vecb(scalD,grid)) sum(opC_p.HyT[iLStmp] *  opC_p.iMy_b * opC_p.Hy_b * vecb(scalD,grid)))

    # print("\n vecb(scalD,grid)", vecb(scalD,grid) ) 

    # print("\n opC_p.iMx_b ", opC_p.iMx_b,"\n" ) 
    # print("\n opC_p.Hx_b ", opC_p.Hx_b ,"\n" ) 


    # print("\n x \n")

    # print("\n test ",opC_p.Hx_b * vecb(scalD,grid))

    # print("\n test ",opC_p.iMx_b * opC_p.Hx_b * vecb(scalD,grid))

    # print("\n y \n")

    # print("\n test ",opC_p.Hy_b * vecb(scalD,grid))

    # print("\n test ", opC_p.iMy_b * opC_p.Hy_b * vecb(scalD,grid) )


    # testvec = ones(2*grid.nx+2*grid.ny)
    # print("\n opC_p.Hx_b * testvec ",opC_p.Hx_b * testvec)
    # #   print("\n test ",opC_p.iMx_b * opC_p.Hx_b * testvec)
    # print("\n test ", opC_p.HxT[iLStmp] * opC_p.iMx_b * opC_p.Hx_b * testvec)

    # print("\n new test ", opC_p.BxT[iLStmp] * opC_p.iMx_b * opC_p.Hx_b * testvec)


    # testvec = ones(2*grid.nx+2*grid.ny)
    # print("\n opC_p.Hy_b * testvec ",opC_p.Hy_b * testvec)
    # #   print("\n test ",opC_p.iMy_b * opC_p.Hy_b * testvec)
    # print("\n test ", opC_p.HyT[iLStmp] * opC_p.iMy_b * opC_p.Hy_b * testvec)

    # print("\n new test ", opC_p.ByT[iLStmp] * opC_p.iMy_b * opC_p.Hy_b * testvec)


  # mass_flux = mass_flux_vec1 .+ mass_flux_vecb .+ mass_flux_veci

    # mass_flux_2 .= reshape(mass_flux,grid)
    mass_flux_vec1_2 .= reshape(mass_flux_vec1,grid)
    mass_flux_vecb_2 .= reshape(mass_flux_vecb,grid)
    mass_flux_veci_2 .= reshape(mass_flux_veci,grid)

    mass_flux .= mass_flux_vec1_2 .+ mass_flux_vecb_2 .+ mass_flux_veci_2

    print("\n sum mass flux all levelsets (walls and interfaces alike) ", sum(mass_flux),"\n ")
    
    print("\n mass_flux ", sum(mass_flux),"\n ")
    print("\n mass_flux_vec1_2 ", sum(mass_flux_vec1_2),"\n ")
    print("\n mass_flux_vecb_2 ", sum(mass_flux_vecb_2),"\n ")
    print("\n mass_flux_veci_2 ", sum(mass_flux_veci_2),"\n ")



    # iplot = 1
    # for jplot in 1:grid.ny
    #     #for iplot in 1:grid.nx
    #     II = CartesianIndex(jplot, iplot) #(id_y, id_x)
    #     pII = lexicographic(II, grid.ny)

    #     if mass_flux_vec1_2[II]>0
    #         printstyled(color=:green, @sprintf "\n j %.5i m %.2e HxT %.2e\n" jplot mass_flux_vec1_2[II] opC_p.HxT[iLStmp][II])
    #         printstyled(color=:red, @sprintf "\n iMx %.10e iMy %.10e \n" opC_p.iMy.diag[pII] opC_p.iMy.diag[pII] )
    #         print("\n B ", II," ",opC_p.Bx[pII,pII]," ",opC_p.BxT[pII,pII])
    #     end

    # end


    # iplot=1
    # jplot = 59
    # II = CartesianIndex(jplot, iplot) #(id_y, id_x)
    # pII = lexicographic(II, grid.ny)
    # printstyled(color=:magenta, @sprintf "\n j %.5i m %.2e HxT %.2e HyT %.2e Bx %.2e By %.2e iMx %.10e iMy %.10e\n" jplot mass_flux_vec1_2[II] opC_p.HxT[iLStmp][II] opC_p.HyT[iLStmp][II] opC_p.Bx[pII,pII] opC_p.By[pII,pII] opC_p.iMy.diag[pII] opC_p.iMy.diag[pII])

    ######################################################################################


    # Matrices for interior BCs
    # for iLS in 1:num.nLS
    #     χx = (geo.dcap[:,:,3] .- geo.dcap[:,:,1]) .^ 2
    #     χy = (geo.dcap[:,:,4] .- geo.dcap[:,:,2]) .^ 2
    #     χ[iLS].diag .= sqrt.(vec(χx .+ χy))
    # end

    # χx = (geo.dcap[:,:,3] .- geo.dcap[:,:,1]) .^ 2
    # χy = (geo.dcap[:,:,4] .- geo.dcap[:,:,2]) .^ 2
    # #     χ[iLS].diag .= sqrt.(vec(χx .+ χy))
    # radial_flux_surf = mass_flux_2 ./ sqrt.(vec(χx .+ χy))
    # # radial_flux_surf = mass_flux_2 ./ χ[1]

    # printstyled(color=:green, @sprintf "\n Radial flux: %.2e \n" radial_flux_surf)

    if num.io_pdi>0
        printstyled(color=:magenta, @sprintf "\n PDI write_mass_flux %.5i \n" num.current_i)
        #nstep needs to be updated beforehand
        @ccall "libpdi".PDI_multi_expose("write_mass_flux"::Cstring,
        "mass_flux"::Cstring, mass_flux::Ptr{Cdouble}, PDI_OUT::Cint,
        "mass_flux_bulk"::Cstring, mass_flux_vec1_2::Ptr{Cdouble}, PDI_OUT::Cint,
        "mass_flux_border"::Cstring, mass_flux_vecb_2::Ptr{Cdouble}, PDI_OUT::Cint,
        "mass_flux_intfc"::Cstring, mass_flux_veci_2::Ptr{Cdouble}, PDI_OUT::Cint,
        C_NULL::Ptr{Cvoid})::Cvoid
    end #if num.io_pdi>0
     
    print("\n sum mass flux ", sum(mass_flux),"\n ")
end


"""
    Integrate mass transfer rate along a line (the interface in 2D), in  in mol/m/s

    to increment mol/m in 2D : use dn=result*dt
    to advect: need to go back to surface mass flux so divide by interface length : mol/m^2/s
    then velocity m/s is obtained  (*M (kg/mol) and * combination of 1/rho_i ))

    need to pay attention to operators of levelsets describing walls (do not compute mass flux at interface outside domain)

    with iLS=1 for example for first levelset

    # Interface
    opC_p.HxT[iLS] (nx*ny, (nx+1)*ny) TODO (nx+1)*ny or nx*(ny+1)
    # Wall
    opC_p.iMx_b ((nx+1)*ny, 2*nx+2*ny)
    opC_p.Hx_b (2*nx+2*ny, 2*nx+2*ny)

    
    computes also in wall described  by LS: you need to select the contribution based on the levelset afterwards
    Interface described by LS number one

    opC_p.Hx_b: left: -dy (cell height along y), bottom: 0, right: +dy, top: 0
    opC_p.Hy_b: left: 0 , bottom: -dx, right: 0, top: +dx


"""
function integrate_mass_flux_over_interface_old(num::Numerical{Float64, Int64},
    grid::Mesh{Flower.GridCC, Float64, Int64},
    opC_pL::Operators{Float64, Int64}, 
    scalD::AbstractArray{Float64, 1},
    mass_flux_vec1::Array{Float64, 1},
    mass_flux_vecb::Array{Float64, 1}, 
    mass_flux_veci::Array{Float64, 1},
    mass_flux::Array{Float64, 2}
    )

    opC_p = opC_pL

    # Interface described by LS number one
    iLStmp=1

    #size (nx*ny)
    mass_flux_vec1 .= 0.0 
    mass_flux_vecb .= 0.0
    mass_flux_veci .= 0.0
    
    #size (ny,nx)
    mass_flux .= 0.0

    mass_flux_vec1   = opC_p.HxT[iLStmp] * opC_p.iMx * opC_p.Bx * vec1(scalD,grid) .+ opC_p.HyT[iLStmp] * opC_p.iMy * opC_p.By * vec1(scalD,grid)
    mass_flux_vecb   = opC_p.HxT[iLStmp] * opC_p.iMx_b * opC_p.Hx_b * vecb(scalD,grid) .+ opC_p.HyT[iLStmp] *  opC_p.iMy_b * opC_p.Hy_b * vecb(scalD,grid)

    for iLS in 1:num.nLS
        mass_flux_veci .+= opC_p.HxT[iLS] * opC_p.iMx * opC_p.Hx[iLS] * veci(scalD,grid,iLS+1)
        mass_flux_veci .+= opC_p.HyT[iLS] * opC_p.iMy * opC_p.Hy[iLS] * veci(scalD,grid,iLS+1)
    end

    # printstyled(color=:red, @sprintf "\n vec1 x y %.2e %.2e \n" sum(opC_p.HxT[iLStmp] * opC_p.iMx * opC_p.Bx * vec1(scalD,grid)) sum(opC_p.HyT[iLStmp] * opC_p.iMy * opC_p.By * vec1(scalD,grid)))
    # printstyled(color=:red, @sprintf "\n vecb x y %.2e %.2e\n" sum(opC_p.HxT[iLStmp] * opC_p.iMx_b * opC_p.Hx_b * vecb(scalD,grid)) sum(opC_p.HyT[iLStmp] *  opC_p.iMy_b * opC_p.Hy_b * vecb(scalD,grid)))

    # print("\n vecb(scalD,grid)", vecb(scalD,grid) ) 

    # print("\n opC_p.iMx_b ", opC_p.iMx_b,"\n" ) 
    # print("\n opC_p.Hx_b ", opC_p.Hx_b ,"\n" ) 


    # print("\n x \n")

    # print("\n test ",opC_p.Hx_b * vecb(scalD,grid))

    # print("\n test ",opC_p.iMx_b * opC_p.Hx_b * vecb(scalD,grid))

    # print("\n y \n")

    # print("\n test ",opC_p.Hy_b * vecb(scalD,grid))

    # print("\n test ", opC_p.iMy_b * opC_p.Hy_b * vecb(scalD,grid) )


    # testvec = ones(2*grid.nx+2*grid.ny)
    # print("\n opC_p.Hx_b * testvec ",opC_p.Hx_b * testvec)
    # #   print("\n test ",opC_p.iMx_b * opC_p.Hx_b * testvec)
    # print("\n test ", opC_p.HxT[iLStmp] * opC_p.iMx_b * opC_p.Hx_b * testvec)

    # print("\n new test ", opC_p.BxT[iLStmp] * opC_p.iMx_b * opC_p.Hx_b * testvec)


    # testvec = ones(2*grid.nx+2*grid.ny)
    # print("\n opC_p.Hy_b * testvec ",opC_p.Hy_b * testvec)
    # #   print("\n test ",opC_p.iMy_b * opC_p.Hy_b * testvec)
    # print("\n test ", opC_p.HyT[iLStmp] * opC_p.iMy_b * opC_p.Hy_b * testvec)

    # print("\n new test ", opC_p.ByT[iLStmp] * opC_p.iMy_b * opC_p.Hy_b * testvec)


  # mass_flux = mass_flux_vec1 .+ mass_flux_vecb .+ mass_flux_veci

    # mass_flux_2 .= reshape(mass_flux,grid)
    mass_flux_vec1_2 .= reshape(mass_flux_vec1,grid)
    mass_flux_vecb_2 .= reshape(mass_flux_vecb,grid)
    mass_flux_veci_2 .= reshape(mass_flux_veci,grid)

    mass_flux .= mass_flux_vec1_2 .+ mass_flux_vecb_2 .+ mass_flux_veci_2

    print("\n sum mass flux all levelsets (walls and interfaces alike) ", sum(mass_flux),"\n ")
    
    print("\n mass_flux ", sum(mass_flux),"\n ")
    print("\n mass_flux_vec1_2 ", sum(mass_flux_vec1_2),"\n ")
    print("\n mass_flux_vecb_2 ", sum(mass_flux_vecb_2),"\n ")
    print("\n mass_flux_veci_2 ", sum(mass_flux_veci_2),"\n ")



    # iplot = 1
    # for jplot in 1:grid.ny
    #     #for iplot in 1:grid.nx
    #     II = CartesianIndex(jplot, iplot) #(id_y, id_x)
    #     pII = lexicographic(II, grid.ny)

    #     if mass_flux_vec1_2[II]>0
    #         printstyled(color=:green, @sprintf "\n j %.5i m %.2e HxT %.2e\n" jplot mass_flux_vec1_2[II] opC_p.HxT[iLStmp][II])
    #         printstyled(color=:red, @sprintf "\n iMx %.10e iMy %.10e \n" opC_p.iMy.diag[pII] opC_p.iMy.diag[pII] )
    #         print("\n B ", II," ",opC_p.Bx[pII,pII]," ",opC_p.BxT[pII,pII])
    #     end

    # end


    # iplot=1
    # jplot = 59
    # II = CartesianIndex(jplot, iplot) #(id_y, id_x)
    # pII = lexicographic(II, grid.ny)
    # printstyled(color=:magenta, @sprintf "\n j %.5i m %.2e HxT %.2e HyT %.2e Bx %.2e By %.2e iMx %.10e iMy %.10e\n" jplot mass_flux_vec1_2[II] opC_p.HxT[iLStmp][II] opC_p.HyT[iLStmp][II] opC_p.Bx[pII,pII] opC_p.By[pII,pII] opC_p.iMy.diag[pII] opC_p.iMy.diag[pII])

    ######################################################################################


    # Matrices for interior BCs
    # for iLS in 1:num.nLS
    #     χx = (geo.dcap[:,:,3] .- geo.dcap[:,:,1]) .^ 2
    #     χy = (geo.dcap[:,:,4] .- geo.dcap[:,:,2]) .^ 2
    #     χ[iLS].diag .= sqrt.(vec(χx .+ χy))
    # end

    # χx = (geo.dcap[:,:,3] .- geo.dcap[:,:,1]) .^ 2
    # χy = (geo.dcap[:,:,4] .- geo.dcap[:,:,2]) .^ 2
    # #     χ[iLS].diag .= sqrt.(vec(χx .+ χy))
    # radial_flux_surf = mass_flux_2 ./ sqrt.(vec(χx .+ χy))
    # # radial_flux_surf = mass_flux_2 ./ χ[1]

    # printstyled(color=:green, @sprintf "\n Radial flux: %.2e \n" radial_flux_surf)

    if num.io_pdi>0
        printstyled(color=:magenta, @sprintf "\n PDI write_mass_flux %.5i \n" num.current_i)
        #nstep needs to be updated beforehand
        @ccall "libpdi".PDI_multi_expose("write_mass_flux"::Cstring,
        "mass_flux"::Cstring, mass_flux::Ptr{Cdouble}, PDI_OUT::Cint,
        "mass_flux_bulk"::Cstring, mass_flux_vec1_2::Ptr{Cdouble}, PDI_OUT::Cint,
        "mass_flux_border"::Cstring, mass_flux_vecb_2::Ptr{Cdouble}, PDI_OUT::Cint,
        "mass_flux_intfc"::Cstring, mass_flux_veci_2::Ptr{Cdouble}, PDI_OUT::Cint,
        C_NULL::Ptr{Cvoid})::Cvoid
    end #if num.io_pdi>0
     
    print("\n sum mass flux ", sum(mass_flux),"\n ")
end


"""
    Integrate mass transfer rate along a line (the interface in 2D), in  in mol/m/s

    to increment mol/m in 2D : use dn=result*dt
    to advect: need to go back to surface mass flux so divide by interface length : mol/m^2/s
    then velocity m/s is obtained  (*M (kg/mol) and * combination of 1/rho_i ))

    need to pay attention to operators of levelsets describing walls (do not compute mass flux at interface outside domain)

    with iLS=1 for example for first levelset

    # Interface
    opC_p.HxT[iLS] (nx*ny, (nx+1)*ny) TODO (nx+1)*ny or nx*(ny+1)
    # Wall
    opC_p.iMx_b ((nx+1)*ny, 2*nx+2*ny)
    opC_p.Hx_b (2*nx+2*ny, 2*nx+2*ny)

    
    computes also in wall described  by LS: you need to select the contribution based on the levelset afterwards
    Interface described by LS number one

    opC_p.Hx_b: left: -dy (cell height along y), bottom: 0, right: +dy, top: 0
    opC_p.Hy_b: left: 0 , bottom: -dx, right: 0, top: +dx


"""
function integrate_mass_flux_over_interface_no_writing(num::Numerical{Float64, Int64},
    grid::Mesh{Flower.GridCC, Float64, Int64},
    opC_pL::Operators{Float64, Int64}, 
    scalD::AbstractArray{Float64, 1},
    mass_flux_vec1::Array{Float64, 1},
    mass_flux_vecb::Array{Float64, 1}, 
    mass_flux_veci::Array{Float64, 1},
    mass_flux_vec1_2::Array{Float64, 2},
    mass_flux_vecb_2::Array{Float64, 2},
    mass_flux_veci_2::Array{Float64, 2},
    mass_flux::Array{Float64, 2}
    )

    opC_p = opC_pL

    # Interface described by LS number one
    iLStmp=1

    #size (nx*ny)
    mass_flux_vec1 .= 0.0 
    mass_flux_vecb .= 0.0
    mass_flux_veci .= 0.0
    
    #size (ny,nx)
    mass_flux .= 0.0
    mass_flux_vec1_2 .= 0.0
    mass_flux_vecb_2 .= 0.0
    mass_flux_veci_2 .= 0.0

    mass_flux_vec1   .= opC_p.HxT[iLStmp] * opC_p.iMx * opC_p.Bx * vec1(scalD,grid) .+ opC_p.HyT[iLStmp] * opC_p.iMy * opC_p.By * vec1(scalD,grid)
    mass_flux_vecb   .= opC_p.HxT[iLStmp] * opC_p.iMx_b * opC_p.Hx_b * vecb(scalD,grid) .+ opC_p.HyT[iLStmp] *  opC_p.iMy_b * opC_p.Hy_b * vecb(scalD,grid)

    for iLS in 1:num.nLS
        mass_flux_veci .+= opC_p.HxT[iLS] * opC_p.iMx * opC_p.Hx[iLS] * veci(scalD,grid,iLS+1)
        mass_flux_veci .+= opC_p.HyT[iLS] * opC_p.iMy * opC_p.Hy[iLS] * veci(scalD,grid,iLS+1)
    end

    # printstyled(color=:red, @sprintf "\n vec1 x y %.2e %.2e \n" sum(opC_p.HxT[iLStmp] * opC_p.iMx * opC_p.Bx * vec1(scalD,grid)) sum(opC_p.HyT[iLStmp] * opC_p.iMy * opC_p.By * vec1(scalD,grid)))
    # printstyled(color=:red, @sprintf "\n vecb x y %.2e %.2e\n" sum(opC_p.HxT[iLStmp] * opC_p.iMx_b * opC_p.Hx_b * vecb(scalD,grid)) sum(opC_p.HyT[iLStmp] *  opC_p.iMy_b * opC_p.Hy_b * vecb(scalD,grid)))

    # print("\n vecb(scalD,grid)", vecb(scalD,grid) ) 

    # print("\n opC_p.iMx_b ", opC_p.iMx_b,"\n" ) 
    # print("\n opC_p.Hx_b ", opC_p.Hx_b ,"\n" ) 


    # print("\n x \n")

    # print("\n test ",opC_p.Hx_b * vecb(scalD,grid))

    # print("\n test ",opC_p.iMx_b * opC_p.Hx_b * vecb(scalD,grid))

    # print("\n y \n")

    # print("\n test ",opC_p.Hy_b * vecb(scalD,grid))

    # print("\n test ", opC_p.iMy_b * opC_p.Hy_b * vecb(scalD,grid) )


    # testvec = ones(2*grid.nx+2*grid.ny)
    # print("\n opC_p.Hx_b * testvec ",opC_p.Hx_b * testvec)
    # #   print("\n test ",opC_p.iMx_b * opC_p.Hx_b * testvec)
    # print("\n test ", opC_p.HxT[iLStmp] * opC_p.iMx_b * opC_p.Hx_b * testvec)

    # print("\n new test ", opC_p.BxT[iLStmp] * opC_p.iMx_b * opC_p.Hx_b * testvec)


    # testvec = ones(2*grid.nx+2*grid.ny)
    # print("\n opC_p.Hy_b * testvec ",opC_p.Hy_b * testvec)
    # #   print("\n test ",opC_p.iMy_b * opC_p.Hy_b * testvec)
    # print("\n test ", opC_p.HyT[iLStmp] * opC_p.iMy_b * opC_p.Hy_b * testvec)

    # print("\n new test ", opC_p.ByT[iLStmp] * opC_p.iMy_b * opC_p.Hy_b * testvec)


  # mass_flux = mass_flux_vec1 .+ mass_flux_vecb .+ mass_flux_veci

    # mass_flux_2 .= reshape(mass_flux,grid)
    mass_flux_vec1_2 .= reshape(mass_flux_vec1,grid)
    mass_flux_vecb_2 .= reshape(mass_flux_vecb,grid)
    mass_flux_veci_2 .= reshape(mass_flux_veci,grid)

    mass_flux .= mass_flux_vec1_2 .+ mass_flux_vecb_2 .+ mass_flux_veci_2

    print("\n sum mass flux all levelsets (walls and interfaces alike) ", sum(mass_flux),"\n ")
    
    print("\n mass_flux ", sum(mass_flux),"\n ")
    print("\n mass_flux_vec1_2 ", sum(mass_flux_vec1_2),"\n ")
    print("\n mass_flux_vecb_2 ", sum(mass_flux_vecb_2),"\n ")
    print("\n mass_flux_veci_2 ", sum(mass_flux_veci_2),"\n ")



    # iplot = 1
    # for jplot in 1:grid.ny
    #     #for iplot in 1:grid.nx
    #     II = CartesianIndex(jplot, iplot) #(id_y, id_x)
    #     pII = lexicographic(II, grid.ny)

    #     if mass_flux_vec1_2[II]>0
    #         printstyled(color=:green, @sprintf "\n j %.5i m %.2e HxT %.2e\n" jplot mass_flux_vec1_2[II] opC_p.HxT[iLStmp][II])
    #         printstyled(color=:red, @sprintf "\n iMx %.10e iMy %.10e \n" opC_p.iMy.diag[pII] opC_p.iMy.diag[pII] )
    #         print("\n B ", II," ",opC_p.Bx[pII,pII]," ",opC_p.BxT[pII,pII])
    #     end

    # end


    # iplot=1
    # jplot = 59
    # II = CartesianIndex(jplot, iplot) #(id_y, id_x)
    # pII = lexicographic(II, grid.ny)
    # printstyled(color=:magenta, @sprintf "\n j %.5i m %.2e HxT %.2e HyT %.2e Bx %.2e By %.2e iMx %.10e iMy %.10e\n" jplot mass_flux_vec1_2[II] opC_p.HxT[iLStmp][II] opC_p.HyT[iLStmp][II] opC_p.Bx[pII,pII] opC_p.By[pII,pII] opC_p.iMy.diag[pII] opC_p.iMy.diag[pII])

    ######################################################################################


    # Matrices for interior BCs
    # for iLS in 1:num.nLS
    #     χx = (geo.dcap[:,:,3] .- geo.dcap[:,:,1]) .^ 2
    #     χy = (geo.dcap[:,:,4] .- geo.dcap[:,:,2]) .^ 2
    #     χ[iLS].diag .= sqrt.(vec(χx .+ χy))
    # end

    # χx = (geo.dcap[:,:,3] .- geo.dcap[:,:,1]) .^ 2
    # χy = (geo.dcap[:,:,4] .- geo.dcap[:,:,2]) .^ 2
    # #     χ[iLS].diag .= sqrt.(vec(χx .+ χy))
    # radial_flux_surf = mass_flux_2 ./ sqrt.(vec(χx .+ χy))
    # # radial_flux_surf = mass_flux_2 ./ χ[1]

    # printstyled(color=:green, @sprintf "\n Radial flux: %.2e \n" radial_flux_surf)
     
    print("\n sum mass flux ", sum(mass_flux),"\n ")
end


"""
    Integrate mass transfer rate along a line (the interface in 2D), in  in mol/m/s

    to increment mol/m in 2D : use dn=result*dt
    to advect: need to go back to surface mass flux so divide by interface length : mol/m^2/s
    then velocity m/s is obtained  (*M (kg/mol) and * combination of 1/rho_i ))

    need to pay attention to operators of levelsets describing walls (do not compute mass flux at interface outside domain)

    with iLS=1 for example for first levelset

    # Interface
    opC_p.HxT[iLS] (nx*ny, (nx+1)*ny) TODO (nx+1)*ny or nx*(ny+1)

    # Wall
    opC_p.iMx_b ((nx+1)*ny, 2*nx+2*ny)
    opC_p.Hx_b (2*nx+2*ny, 2*nx+2*ny)

    computes also in wall described  by LS: you need to select the contribution based on the levelset afterwards
    Interface described by LS number one

    opC_p.Hx_b: left: -dy (cell height along y), bottom: 0, right: +dy, top: 0
    opC_p.Hy_b: left: 0 , bottom: -dx, right: 0, top: +dx
"""
function integrate_mass_flux_over_interface_2(num::Numerical{Float64, Int64},
    grid::Mesh{Flower.GridCC, Float64, Int64},
    opC_pL::Operators{Float64, Int64}, 
    scalD::AbstractArray{Float64, 1},
    mass_flux_vec1::Array{Float64, 1},
    mass_flux_vecb::Array{Float64, 1}, 
    mass_flux_veci::Array{Float64, 1},
    mass_flux::Array{Float64, 2}
    )

    opC_p = opC_pL


    #size (nx*ny)
    mass_flux_vec1 .= 0.0 
    mass_flux_vecb .= 0.0
    mass_flux_veci .= 0.0
    
    #size (ny,nx)
    mass_flux .= 0.0

    mass_flux_vec1   = opC_p.BxT * opC_p.iMx * opC_p.Bx * vec1(scalD,grid) .+ opC_p.ByT * opC_p.iMy * opC_p.By * vec1(scalD,grid)
    mass_flux_vecb   = opC_p.BxT * opC_p.iMx_b * opC_p.Hx_b * vecb(scalD,grid) .+ opC_p.ByT *  opC_p.iMy_b * opC_p.Hy_b * vecb(scalD,grid)

    for iLS in 1:num.nLS
        #TODO
        # mass_flux_veci .+= opC_p.HxT[iLS] * opC_p.iMx * opC_p.Hx[iLS] * veci(scalD,grid,iLS+1)
        # mass_flux_veci .+= opC_p.HyT[iLS] * opC_p.iMy * opC_p.Hy[iLS] * veci(scalD,grid,iLS+1)
        mass_flux_veci .+= opC_p.BxT * opC_p.iMx * opC_p.Hx[iLS] * veci(scalD,grid,iLS+1)
        mass_flux_veci .+= opC_p.ByT * opC_p.iMy * opC_p.Hy[iLS] * veci(scalD,grid,iLS+1)
    end

    # printstyled(color=:red, @sprintf "\n vec1 x y %.2e %.2e \n" sum(opC_p.BxT * opC_p.iMx * opC_p.Bx * vec1(scalD,grid)) sum(opC_p.ByT * opC_p.iMy * opC_p.By * vec1(scalD,grid)))
    # printstyled(color=:red, @sprintf "\n vecb x y %.2e %.2e\n" sum(opC_p.BxT * opC_p.iMx_b * opC_p.Hx_b * vecb(scalD,grid)) sum(opC_p.ByT *  opC_p.iMy_b * opC_p.Hy_b * vecb(scalD,grid)))

    # print("\n vecb(scalD,grid)", vecb(scalD,grid) ) 

    # print("\n opC_p.iMx_b ", opC_p.iMx_b,"\n" ) 
    # print("\n opC_p.Hx_b ", opC_p.Hx_b ,"\n" ) 

    # print("\n x \n")

    # print("\n test ",opC_p.Hx_b * vecb(scalD,grid))

    # print("\n test ",opC_p.iMx_b * opC_p.Hx_b * vecb(scalD,grid))

    # print("\n y \n")

    # print("\n test ",opC_p.Hy_b * vecb(scalD,grid))

    # print("\n test ", opC_p.iMy_b * opC_p.Hy_b * vecb(scalD,grid) )


    # testvec = ones(2*grid.nx+2*grid.ny)
    # print("\n opC_p.Hx_b * testvec ",opC_p.Hx_b * testvec)
    # #   print("\n test ",opC_p.iMx_b * opC_p.Hx_b * testvec)

    # # print("\n new test ", opC_p.BxT[iLStmp] * opC_p.iMx_b * opC_p.Hx_b * testvec)
    # print("\n new test ", opC_p.BxT * opC_p.iMx_b * opC_p.Hx_b * testvec)


    # testvec = ones(2*grid.nx+2*grid.ny)
    # print("\n opC_p.Hy_b * testvec ",opC_p.Hy_b * testvec)
    #   print("\n test ",opC_p.iMy_b * opC_p.Hy_b * testvec)
    # print("\n test ", opC_p.HyT[iLStmp] * opC_p.iMy_b * opC_p.Hy_b * testvec)

    # print("\n new test ", opC_p.ByT * opC_p.iMy_b * opC_p.Hy_b * testvec)


  # mass_flux = mass_flux_vec1 .+ mass_flux_vecb .+ mass_flux_veci

    # mass_flux_2 .= reshape(mass_flux,grid)
    mass_flux_vec1_2 .= reshape(mass_flux_vec1,grid)
    mass_flux_vecb_2 .= reshape(mass_flux_vecb,grid)
    mass_flux_veci_2 .= reshape(mass_flux_veci,grid)

    mass_flux .= mass_flux_vec1_2 .+ mass_flux_vecb_2 .+ mass_flux_veci_2

    print("\n sum mass flux all levelsets (walls and interfaces alike) ", sum(mass_flux),"\n ")
    
    print("\n mass_flux ", sum(mass_flux),"\n ")
    print("\n mass_flux_vec1_2 ", sum(mass_flux_vec1_2),"\n ")
    print("\n mass_flux_vecb_2 ", sum(mass_flux_vecb_2),"\n ")
    print("\n mass_flux_veci_2 ", sum(mass_flux_veci_2),"\n ")



    # iplot = 1
    # for jplot in 1:grid.ny
    #     #for iplot in 1:grid.nx
    #     II = CartesianIndex(jplot, iplot) #(id_y, id_x)
    #     pII = lexicographic(II, grid.ny)

    #     if mass_flux_vec1_2[II]>0
    #         printstyled(color=:green, @sprintf "\n j %.5i m %.2e HxT %.2e\n" jplot mass_flux_vec1_2[II] opC_p.HxT[iLStmp][II])
    #         printstyled(color=:red, @sprintf "\n iMx %.10e iMy %.10e \n" opC_p.iMy.diag[pII] opC_p.iMy.diag[pII] )
    #         print("\n B ", II," ",opC_p.Bx[pII,pII]," ",opC_p.BxT[pII,pII])
    #     end

    # end


    # iplot=1
    # jplot = 59
    # II = CartesianIndex(jplot, iplot) #(id_y, id_x)
    # pII = lexicographic(II, grid.ny)
    # printstyled(color=:magenta, @sprintf "\n j %.5i m %.2e HxT %.2e HyT %.2e Bx %.2e By %.2e iMx %.10e iMy %.10e\n" jplot mass_flux_vec1_2[II] opC_p.HxT[iLStmp][II] opC_p.HyT[iLStmp][II] opC_p.Bx[pII,pII] opC_p.By[pII,pII] opC_p.iMy.diag[pII] opC_p.iMy.diag[pII])

    ######################################################################################


    # Matrices for interior BCs
    # for iLS in 1:num.nLS
    #     χx = (geo.dcap[:,:,3] .- geo.dcap[:,:,1]) .^ 2
    #     χy = (geo.dcap[:,:,4] .- geo.dcap[:,:,2]) .^ 2
    #     χ[iLS].diag .= sqrt.(vec(χx .+ χy))
    # end

    # χx = (geo.dcap[:,:,3] .- geo.dcap[:,:,1]) .^ 2
    # χy = (geo.dcap[:,:,4] .- geo.dcap[:,:,2]) .^ 2
    # #     χ[iLS].diag .= sqrt.(vec(χx .+ χy))
    # radial_flux_surf = mass_flux_2 ./ sqrt.(vec(χx .+ χy))
    # # radial_flux_surf = mass_flux_2 ./ χ[1]

    # printstyled(color=:green, @sprintf "\n Radial flux: %.2e \n" radial_flux_surf)
    
    if num.io_pdi>0
        printstyled(color=:magenta, @sprintf "\n PDI write_mass_flux %.5i \n" num.current_i)
        #nstep needs to be updated beforehand
        @ccall "libpdi".PDI_multi_expose("write_mass_flux"::Cstring,
        "mass_flux"::Cstring, mass_flux::Ptr{Cdouble}, PDI_OUT::Cint,
        "mass_flux_bulk"::Cstring, mass_flux_vec1_2::Ptr{Cdouble}, PDI_OUT::Cint,
        "mass_flux_border"::Cstring, mass_flux_vecb_2::Ptr{Cdouble}, PDI_OUT::Cint,
        "mass_flux_intfc"::Cstring, mass_flux_veci_2::Ptr{Cdouble}, PDI_OUT::Cint,
        C_NULL::Ptr{Cvoid})::Cvoid
    end #if num.io_pdi>0

    print("\n sum mass flux ", sum(mass_flux),"\n ")
end


"""
    Integrate mass transfer rate along a line (the interface in 2D), in  in mol/m/s

    to increment mol/m in 2D : use dn=result*dt
    to advect: need to go back to surface mass flux so divide by interface length : mol/m^2/s
    then velocity m/s is obtained  (*M (kg/mol) and * combination of 1/rho_i ))

    need to pay attention to operators of levelsets describing walls (do not compute mass flux at interface outside domain)

    with iLS=1 for example for first levelset

    # Interface
    opC_p.HxT[iLS] (nx*ny, (nx+1)*ny) TODO (nx+1)*ny or nx*(ny+1)

    # Wall
    opC_p.iMx_b ((nx+1)*ny, 2*nx+2*ny)
    opC_p.Hx_b (2*nx+2*ny, 2*nx+2*ny)

    computes also in wall described  by LS: you need to select the contribution based on the levelset afterwards
    Interface described by LS number one

    opC_p.Hx_b: left: -dy (cell height along y), bottom: 0, right: +dy, top: 0
    opC_p.Hy_b: left: 0 , bottom: -dx, right: 0, top: +dx
"""
function integrate_mass_flux_over_interface_2_no_writing(num::Numerical{Float64, Int64},
    grid::Mesh{Flower.GridCC, Float64, Int64},
    opC_pL::Operators{Float64, Int64}, 
    scalD::AbstractArray{Float64, 1},
    mass_flux_vec1::Array{Float64, 1},
    mass_flux_vecb::Array{Float64, 1}, 
    mass_flux_veci::Array{Float64, 1},
    mass_flux_vec1_2::Array{Float64, 2},
    mass_flux_vecb_2::Array{Float64, 2},
    mass_flux_veci_2::Array{Float64, 2},
    mass_flux::Array{Float64, 2}
    )

    opC_p = opC_pL


    #size (nx*ny)
    mass_flux_vec1 .= 0.0 
    mass_flux_vecb .= 0.0
    mass_flux_veci .= 0.0
    
    #size (ny,nx)
    mass_flux .= 0.0
    mass_flux_vec1_2 .= 0.0
    mass_flux_vecb_2 .= 0.0
    mass_flux_veci_2 .= 0.0

    mass_flux_vec1   .= opC_p.BxT * opC_p.iMx * opC_p.Bx * vec1(scalD,grid) .+ opC_p.ByT * opC_p.iMy * opC_p.By * vec1(scalD,grid)
    mass_flux_vecb   .= opC_p.BxT * opC_p.iMx_b * opC_p.Hx_b * vecb(scalD,grid) .+ opC_p.ByT *  opC_p.iMy_b * opC_p.Hy_b * vecb(scalD,grid)

    for iLS in 1:num.nLS
        mass_flux_veci .+= opC_p.BxT * opC_p.iMx * opC_p.Hx[iLS] * veci(scalD,grid,iLS+1)
        mass_flux_veci .+= opC_p.ByT * opC_p.iMy * opC_p.Hy[iLS] * veci(scalD,grid,iLS+1)
    end

    # printstyled(color=:red, @sprintf "\n vec1 x y %.2e %.2e \n" sum(opC_p.BxT * opC_p.iMx * opC_p.Bx * vec1(scalD,grid)) sum(opC_p.ByT * opC_p.iMy * opC_p.By * vec1(scalD,grid)))
    # printstyled(color=:red, @sprintf "\n vecb x y %.2e %.2e\n" sum(opC_p.BxT * opC_p.iMx_b * opC_p.Hx_b * vecb(scalD,grid)) sum(opC_p.ByT *  opC_p.iMy_b * opC_p.Hy_b * vecb(scalD,grid)))

    # print("\n vecb(scalD,grid)", vecb(scalD,grid) ) 

    # print("\n opC_p.iMx_b ", opC_p.iMx_b,"\n" ) 
    # print("\n opC_p.Hx_b ", opC_p.Hx_b ,"\n" ) 


    # print("\n x \n")

    # print("\n test ",opC_p.Hx_b * vecb(scalD,grid))

    # print("\n test ",opC_p.iMx_b * opC_p.Hx_b * vecb(scalD,grid))

    # print("\n y \n")

    # print("\n test ",opC_p.Hy_b * vecb(scalD,grid))

    # print("\n test ", opC_p.iMy_b * opC_p.Hy_b * vecb(scalD,grid) )


    # testvec = ones(2*grid.nx+2*grid.ny)
    # print("\n opC_p.Hx_b * testvec ",opC_p.Hx_b * testvec)
    # #   print("\n test ",opC_p.iMx_b * opC_p.Hx_b * testvec)

    # # print("\n new test ", opC_p.BxT[iLStmp] * opC_p.iMx_b * opC_p.Hx_b * testvec)
    # print("\n new test ", opC_p.BxT * opC_p.iMx_b * opC_p.Hx_b * testvec)


    # testvec = ones(2*grid.nx+2*grid.ny)
    # print("\n opC_p.Hy_b * testvec ",opC_p.Hy_b * testvec)
    #   print("\n test ",opC_p.iMy_b * opC_p.Hy_b * testvec)
    # print("\n test ", opC_p.HyT[iLStmp] * opC_p.iMy_b * opC_p.Hy_b * testvec)

    # print("\n new test ", opC_p.ByT * opC_p.iMy_b * opC_p.Hy_b * testvec)


  # mass_flux = mass_flux_vec1 .+ mass_flux_vecb .+ mass_flux_veci

    # mass_flux_2 .= reshape(mass_flux,grid)
    mass_flux_vec1_2 .= reshape(mass_flux_vec1,grid)
    mass_flux_vecb_2 .= reshape(mass_flux_vecb,grid)
    mass_flux_veci_2 .= reshape(mass_flux_veci,grid)

    mass_flux .= mass_flux_vec1_2 .+ mass_flux_vecb_2 .+ mass_flux_veci_2

    print("\n sum mass flux all levelsets (walls and interfaces alike) ", sum(mass_flux),"\n ")
    
    print("\n mass_flux ", sum(mass_flux),"\n ")
    print("\n mass_flux_vec1_2 ", sum(mass_flux_vec1_2),"\n ")
    print("\n mass_flux_vecb_2 ", sum(mass_flux_vecb_2),"\n ")
    print("\n mass_flux_veci_2 ", sum(mass_flux_veci_2),"\n ")


    print("\n test mass flux  ", mass_flux[div(grid.ny,2),:],"\n")
    print("\n test mass_flux_vec1_2", mass_flux_vec1_2[div(grid.ny,2),:],"\n")
    print("\n test mass flux b", mass_flux_vecb_2[div(grid.ny,2),:],"\n")
    print("\n test mass_flux_veci_2", mass_flux_veci_2[div(grid.ny,2),:],"\n")

    # iplot = 1
    # for jplot in 1:grid.ny
    #     #for iplot in 1:grid.nx
    #     II = CartesianIndex(jplot, iplot) #(id_y, id_x)
    #     pII = lexicographic(II, grid.ny)

    #     if mass_flux_vec1_2[II]>0
    #         printstyled(color=:green, @sprintf "\n j %.5i m %.2e HxT %.2e\n" jplot mass_flux_vec1_2[II] opC_p.HxT[iLStmp][II])
    #         printstyled(color=:red, @sprintf "\n iMx %.10e iMy %.10e \n" opC_p.iMy.diag[pII] opC_p.iMy.diag[pII] )
    #         print("\n B ", II," ",opC_p.Bx[pII,pII]," ",opC_p.BxT[pII,pII])
    #     end

    # end


    # iplot=1
    # jplot = 59
    # II = CartesianIndex(jplot, iplot) #(id_y, id_x)
    # pII = lexicographic(II, grid.ny)
    # printstyled(color=:magenta, @sprintf "\n j %.5i m %.2e HxT %.2e HyT %.2e Bx %.2e By %.2e iMx %.10e iMy %.10e\n" jplot mass_flux_vec1_2[II] opC_p.HxT[iLStmp][II] opC_p.HyT[iLStmp][II] opC_p.Bx[pII,pII] opC_p.By[pII,pII] opC_p.iMy.diag[pII] opC_p.iMy.diag[pII])

    ######################################################################################


    # Matrices for interior BCs
    # for iLS in 1:num.nLS
    #     χx = (geo.dcap[:,:,3] .- geo.dcap[:,:,1]) .^ 2
    #     χy = (geo.dcap[:,:,4] .- geo.dcap[:,:,2]) .^ 2
    #     χ[iLS].diag .= sqrt.(vec(χx .+ χy))
    # end

    # χx = (geo.dcap[:,:,3] .- geo.dcap[:,:,1]) .^ 2
    # χy = (geo.dcap[:,:,4] .- geo.dcap[:,:,2]) .^ 2
    # #     χ[iLS].diag .= sqrt.(vec(χx .+ χy))
    # radial_flux_surf = mass_flux_2 ./ sqrt.(vec(χx .+ χy))
    # # radial_flux_surf = mass_flux_2 ./ χ[1]

    # printstyled(color=:green, @sprintf "\n Radial flux: %.2e \n" radial_flux_surf)
   
    print("\n sum mass flux ", sum(mass_flux),"\n ")
end