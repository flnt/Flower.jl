# using Revise
using Flower

#From poisson.jl


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
function integrate_mass_flux_over_interface_3(num::Numerical{Float64, Int64},
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
function integrate_mass_flux_over_interface_3_no_writing(num::Numerical{Float64, Int64},
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

function compute_grad_T_x_T_y_array_test!(num_LS, grid, grid_u, grid_v, opC_p, grad_x, grad_y, TD)

    test_tolerance = 1e-14

    verbosity = 0
    # verbosity = 3

    j = div(grid.ny,2)
    i = div(grid.nx,2)

    @testset " capacities bulk" begin
        if verbosity > 2 
            printstyled(color=:green, @sprintf "\n i: %.3i j: %.3i\n" i j) 
        end

        II = CartesianIndex(j,i)
        pII = lexicographic(II,grid.ny)


        if verbosity > 2 
            print("\n iMx ", opC_p.iMx[pII,pII], "\n")
            print("\n iMx ", opC_p.iMx.diag[pII], "\n")
        end

        @test opC_p.iMx[pII,pII] ≈ 1.0/100.0 atol=test_tolerance #volume h_x*h_y = 10*10


        @test opC_p.Bx[pII,pII] ≈ 10.0 atol=test_tolerance #h_y = 10

        if verbosity > 2 
            print("\n Bx ", opC_p.Bx[pII,pII], "\n")
            print("\n BxT ", opC_p.BxT[pII,pII], "\n")

            # print("\n capacities ", grid.LS[1].geoL.cap[II,:], "\n")

            print("\n capacities ", grid.LS[1].geoL.dcap[II,:], "\n")


            # geo_v[end].dcap[II,8]

            # tmp = lexicographic(II,gp.ny)
            # II is a two-dimensional index, and in 1D index is required. 
            # It is recommended to use , it's faster than accessing .
        end

    end

    ∇ϕ_x = opC_p.iMx * opC_p.Bx * vec1(TD,grid) .+ opC_p.iMx_b * opC_p.Hx_b * vecb(TD,grid)
    ∇ϕ_y = opC_p.iMy * opC_p.By * vec1(TD,grid) .+ opC_p.iMy_b * opC_p.Hy_b * vecb(TD,grid)

    for iLS in 1:num_LS
        ∇ϕ_x .+= opC_p.iMx * opC_p.Hx[iLS] * veci(TD,grid,iLS+1)
        ∇ϕ_y .+= opC_p.iMy * opC_p.Hy[iLS] * veci(TD,grid,iLS+1)
    end

    grad_x .= reshape(veci(∇ϕ_x,grid_u,1), grid_u)
    grad_y .= reshape(veci(∇ϕ_y,grid_v,1), grid_v)

    if verbosity > 2 
        printstyled(color=:green, @sprintf "\n grad \n")

        printstyled(color=:green, @sprintf "\n i: %.3i j: %.3i\n" i j)


        print("\n grad_x ", grad_x[j,i], "\n")
        print("\n grad_y ", grad_y[j,i], "\n")


        j = div(grid.ny,2)
        i = 1

        printstyled(color=:green, @sprintf "\n i: %.3i j: %.3i\n" i j)



        print("\n grad_x line", grad_x[j,:], "\n")

        print("\n length grad_x line ", length(grad_x[j,:]), "\n")

        print("\n gp.x line", grid.x[j,:], "\n")
        print("\n gp.dx line", grid.dx[j,:], "\n")
        
        print("\n length gp.x line", length(grid.x[j,:]), "\n")
        print("\n length gp.dx line", length(grid.dx[j,:]), "\n")

        print("\n gu.x line", grid_u.x[j,:], "\n")
        print("\n gu.dx line", grid_u.dx[j,:], "\n")
        
        print("\n length gu.x line", length(grid_u.x[j,:]), "\n")
        print("\n length gu.dx line", length(grid_u.dx[j,:]), "\n")


        print("\n grad_x ", grad_x[j,i], "\n")
        print("\n grad_y ", grad_y[j,i], "\n")

        II = CartesianIndex(j,i)
        pII = lexicographic(II,grid.ny)


        @test opC_p.iMx[pII,pII] ≈ 2.0/100.0 atol=test_tolerance #volume h_x*h_y/2 = 10*10/2
        print("\n iMx ", opC_p.iMx[pII,pII], "\n")
        print("\n iMx ", opC_p.iMx.diag[pII], "\n")

        @test opC_p.Bx[pII,pII] ≈ 10.0 atol=test_tolerance


        print("\n Bx ", opC_p.Bx[pII,pII], "\n")
        print("\n BxT ", opC_p.BxT[pII,pII], "\n")

        print("\n Bx ", opC_p.Bx[pII,:], "\n")
        print("\n Bx ", opC_p.Bx[:,pII], "\n")


        print("\n By ", opC_p.By[pII,pII], "\n")
        print("\n ByT ", opC_p.ByT[pII,pII], "\n")

        # print("\n capacities ", grid.LS[1].geoL.cap[II,:], "\n")

        print("\n capacities ", grid.LS[1].geoL.dcap[II,:], "\n")
    end

end

function compute_grad_T_x_array!(num_LS, grid, grid_u, opC_p, grad_x, TD)
    
    ∇ϕ_x = opC_p.iMx * opC_p.Bx * vec1(TD,grid) .+ opC_p.iMx_b * opC_p.Hx_b * vecb(TD,grid)

    for iLS in 1:num_LS
        ∇ϕ_x .+= opC_p.iMx * opC_p.Hx[iLS] * veci(TD,grid,iLS+1)
    end

    grad_x .= reshape(veci(∇ϕ_x,grid_u,1), grid_u)

end

function compute_grad_T_y_array!(num_LS, grid, grid_v, opC_p, grad_y, TD)
    
    ∇ϕ_y = opC_p.iMy * opC_p.By * vec1(TD,grid) .+ opC_p.iMy_b * opC_p.Hy_b * vecb(TD,grid)

    for iLS in 1:num_LS
        ∇ϕ_y .+= opC_p.iMy * opC_p.Hy[iLS] * veci(TD,grid,iLS+1)
    end

    grad_y .= reshape(veci(∇ϕ_y,grid_v,1), grid_v)

end


function ftest(x,y)
    return x
end

function ftest_1(x,y)
    return 1
end

verbosity = 0

#tolerance for tests
test_tolerance = 1.e-14


L0 = 100

n = 10

x = LinRange(0, L0, n+1)
y = LinRange(0, L0, n+1)


ϵ = 0.001

epsilon=ϵ

radius=2.5e-5 

xcoord = 5e-5
ycoord=5e-5

xcoord = -xcoord
ycoord = -ycoord


num = Numerical(
    case = "Cylinder",
    x = x,
    y = y,
    CFL = 1.0,
    max_iterations = 0,
    shifted = xcoord,
    shifted_y = ycoord,
    R = radius,
    ϵ = ϵ,
    solve_solid=1,
)

gp, gu, gv = init_meshes(num)
op, phS, phL = init_fields(num, gp, gu, gv)

gp.LS[1].u .= 1.0 #deactivate interface


@time run_forward!(num, gp, gu, gv, op, phS, phL; navier_stokes = true)
#navier_stokes = true 
# or:
# run_forward!(num, gp, gu, gv, op, phS, phL)
# update_all_ls_data(num, gp, gu, gv, BC_int, true, true, false)
# laps = set_matrices!(num, gp, [gp.LS[1].geoL], gu, [gu.LS[1].geoL], gv, [gv.LS[1].geoL], op.opC_pL, op.opC_uL, op.opC_vL, true, true)
# Lp, bc_Lp, bc_Lp_b, Lu, bc_Lu, bc_Lu_b, Lv, bc_Lv, bc_Lv_b = laps



mass_flux_vec1 = fzeros(gp)
mass_flux_vecb = fzeros(gp)
mass_flux_veci = fzeros(gp)
mass_flux = zeros(gp)

tmp_vec_p = zeros(gp)
tmp_vec_p0 = zeros(gp)
tmp_vec_p1 = zeros(gp)
@testset "Phase change: mass flux" begin

@testset "Phase change: mass flux" begin
    phL.TD .= 1.0 
    mass_flux_vec1 = fzeros(gp)
    mass_flux_vecb = fzeros(gp)
    mass_flux_veci = fzeros(gp)
    mass_flux = zeros(gp)
    # tmp_vec_p = zeros(gp)
    # tmp_vec_p0 = zeros(gp)
    # tmp_vec_p1 = zeros(gp)

print("\n TODO check capacities \n")
integrate_mass_flux_over_interface_2_no_writing(num,gp,op.opC_pL,phL.TD,mass_flux_vec1,mass_flux_vecb,mass_flux_veci,tmp_vec_p,tmp_vec_p0,tmp_vec_p1,mass_flux)
# @test sum(mass_flux) == 0 
@test sum(mass_flux) ≈ 0 atol=test_tolerance


integrate_mass_flux_over_interface(num,gp,op.opC_pL,
phL.TD,
mass_flux_vec1,
mass_flux_vecb,mass_flux_veci, tmp_vec_p, tmp_vec_p0, tmp_vec_p1, mass_flux,num.index_phase_change)

@test sum(mass_flux) ≈ 0 atol=test_tolerance

end #"Phase change: mass flux" begin

# @testset "Phase change: mass flux old" begin
#     phL.TD .= 1.0 
#     mass_flux_vec1 = fzeros(gp)
#     mass_flux_vecb = fzeros(gp)
#     mass_flux_veci = fzeros(gp)
#     mass_flux = zeros(gp)
#     # tmp_vec_p = zeros(gp)
#     # tmp_vec_p0 = zeros(gp)
#     # tmp_vec_p1 = zeros(gp)

# integrate_mass_flux_over_interface_no_writing(num,gp,op.opC_pL,phL.TD,mass_flux_vec1,mass_flux_vecb,mass_flux_veci,tmp_vec_p,tmp_vec_p0,tmp_vec_p1,mass_flux)
# # @test sum(mass_flux) == 0 
# @test sum(mass_flux) ≈ 0 atol=test_tolerance
# end #"Phase change: mass flux" begin

end #phase change

@testset "Interpolation" begin
gp.V .= 1.0
interpolate_scalar!(gp, gu, gv, gp.V, gu.V, gv.V)

@test minimum(gu.V) ≈ 1.0 atol=test_tolerance
@test maximum(gu.V) ≈ 1.0 atol=test_tolerance
@test minimum(gv.V) ≈ 1.0 atol=test_tolerance
@test maximum(gv.V) ≈ 1.0 atol=test_tolerance

# tmp_vec_p = zeros(gp) 
# tmp_vec_p0 = zeros(gp) 

tmp_vec_u = zeros(gu) 
tmp_vec_v = zeros(gv) 

# tmp_vec_u .= 1.0 
# tmp_vec_v .= 1.0
# interpolate_grid_liquid!(gp,gu,gv,tmp_vec_u, tmp_vec_v,tmp_vec_p,tmp_vec_p0)
# interpolate_grid_solid!(gp,gu,gv,tmp_vec_u, tmp_vec_v,tmp_vec_p,tmp_vec_p0)

# @test minimum(tmp_vec_p) ≈ 1.0 atol=test_tolerance
# @test maximum(tmp_vec_p) ≈ 1.0 atol=test_tolerance
# @test minimum(tmp_vec_p0) ≈ 1.0 atol=test_tolerance
# @test maximum(tmp_vec_p0) ≈ 1.0 atol=test_tolerance

# tmp_vec_u .= 1.0 
# tmp_vec_v .= 1.0
# interpolate_grid_liquid_2!(num,gp,gu.LS[end],gv.LS[end],tmp_vec_u, tmp_vec_v,tmp_vec_p,tmp_vec_p0)

# interpolate_grid_solid_2!(num,gp,gu.LS[end],gv.LS[end],tmp_vec_u, tmp_vec_v,tmp_vec_p,tmp_vec_p0)

# @test minimum(tmp_vec_p) ≈ 1.0 atol=test_tolerance
# @test maximum(tmp_vec_p) ≈ 1.0 atol=test_tolerance
# @test minimum(tmp_vec_p0) ≈ 1.0 atol=test_tolerance
# @test maximum(tmp_vec_p0) ≈ 1.0 atol=test_tolerance

tmp_vec_u .= 1.0 
tmp_vec_v .= 1.0
interpolate_grid_liquid_solid!(num,gp,gu.LS[end],gv.LS[end],tmp_vec_u, tmp_vec_v,tmp_vec_p,tmp_vec_p0)

@test minimum(tmp_vec_p) ≈ 1.0 atol=test_tolerance
@test maximum(tmp_vec_p) ≈ 1.0 atol=test_tolerance
@test minimum(tmp_vec_p0) ≈ 1.0 atol=test_tolerance
@test maximum(tmp_vec_p0) ≈ 1.0 atol=test_tolerance

# LS_u =grid_u.LS[1]
# LS_v = grid_v.LS[1]
# us .= (
#     (u[:,2:end] .* LS_u.geoL.dcap[:,2:end,6] .+ 
#     u[:,1:end-1] .* LS_u.geoL.dcap[:,1:end-1,6]) ./ 
#     (LS_u.geoL.dcap[:,1:end-1,6] .+ LS_u.geoL.dcap[:,2:end,6])
# )
# vs .= (
#     (v[2:end,:] .* LS_v.geoL.dcap[2:end,:,7] .+ 
#     v[1:end-1,:] .* LS_v.geoL.dcap[1:end-1,:,7]) ./
#     (LS_v.geoL.dcap[1:end-1,:,7] .+ LS_v.geoL.dcap[2:end,:,7])
# )

# u = phL.Eu
# v = phL.Ev

if verbosity>2

    for j in 1:gp.ny
        for i in 1:gp.nx
            if (tmp_vec_p[j,i] != 1.0) 
                print("\n j",j," i ",i, " tmp_vec_p ",tmp_vec_p[j,i]," tmp_vec_p0 ",tmp_vec_p0[j,i])
                print("\n LS_u.geoL.dcap[:,2:end,6] ",gu.LS[1].geoL.dcap[j,i,6]," ",gu.LS[1].geoL.dcap[j,i+1,6])
            end
        end
    end

    for j in 1:gp.ny
        for i in 1:gp.nx
            if (tmp_vec_p0[j,i] != 1.0) 
                print("\n j",j," i ",i, " tmp_vec_p ",tmp_vec_p[j,i]," tmp_vec_p0 ",tmp_vec_p0[j,i])
                print("\n LS_u.geoL.dcap[:,2:end,6] ",gu.LS[1].geoL.dcap[j,i,6]," ",gu.LS[1].geoL.dcap[j,i+1,6])
            end
        end
    end

end

end


print("\n Test interpolation linear function ")


j = div(n,2)


@testset "Interpolation: linear" begin

    # Interpolation 
    # Example f(x)=x on scalar grid
    # Scalar
    # gp.V j [5.0, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0, 85.0, 95.0]
    # u 
    # gu.V j [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0]
    # v
    # gv.V j [5.0, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0, 85.0, 95.0]

    print("\n interpolation size ",size(gp.V))
    gp.V .= ftest.(gp.x,gp.y)

    interpolate_scalar!(gp, gu, gv, gp.V, gu.V, gv.V)
    
    # @test minimum(gu.V) ≈ 1.0 atol=test_tolerance
    # @test maximum(gu.V) ≈ 1.0 atol=test_tolerance
    # @test minimum(gv.V) ≈ 1.0 atol=test_tolerance
    # @test maximum(gv.V) ≈ 1.0 atol=test_tolerance

    print("\n gp.V j ", gp.V[j,:])
    print("\n gu.V j ", gu.V[j,:])
    print("\n gv.V j ", gv.V[j,:])

    
    # tmp_vec_p = zeros(gp) 
    # tmp_vec_p0 = zeros(gp) 
    
    # tmp_vec_u = zeros(gu) 
    # tmp_vec_v = zeros(gv) 

end



#region interpolate_scalar_to_staggered_u_v_grids_at_border!
@testset "interpolate_scalar_to_staggered_u_v_grids_at_border!" begin

    printstyled(color=:green, @sprintf "\n interpolate_scalar_to_staggered_u_v_grids_at_border!" )

    @unpack Bx, By, Hx, Hy, HxT, HyT, χ, M, iMx, iMy, Hx_b, Hy_b, HxT_b, HyT_b, iMx_b, iMy_b, iMx_bd, iMy_bd, χ_b = op.opC_pL
    @unpack BxT, ByT,tmp_x, tmp_y = op.opC_pL

    grid = gp
    grid_u = gu
    grid_v = gv
   
    ni = grid.nx * grid.ny
    nb = 2 * grid.nx + 2 * grid.ny

    # #TODO reset zero
    # rhs .= 0.0
    # coeffDu .= 0.0
    # coeffDv .= 0.0
    # A .= 0.0
    # a0 .= 0.0

    # a0_b = zeros(nb)
    # _a1_b = zeros(nb)
    # _b_b = zeros(nb)
    # for iLS in 1:num.nLS
    #     set_borders!(grid, grid.LS[iLS].cl, grid.LS[iLS].u, a0_b, _a1_b, _b_b, BC, num.n_ext_cl)
    # end
    # a1_b = Diagonal(vec(_a1_b))
    # b_b = Diagonal(vec(_b_b))

    coeffD = fnones(grid,num)

    vecb(coeffD,grid) .= 9.0 # we test how the wall conductivity intervenes

    coeffDu = zeros(gu)
    coeffDv = zeros(gv)


    # interpolate coefficient from p grid to u and v grids
    coeffD_borders = vecb(coeffD,grid)

    interpolate_scalar!(grid, grid_u, grid_v, reshape(veci(coeffD,grid,1), grid), coeffDu, coeffDv)

    print("\n coeff ",minimum(coeffDu)," ",maximum(coeffDu)," ",minimum(coeffDv)," ",maximum(coeffDv)," ",minimum(coeffD_borders)," ",maximum(coeffD_borders))

    coeffDx_bulk = veci(coeffDu,grid_u)
    coeffDy_bulk = veci(coeffDv,grid_v)

    # mat_coeffDx = Diagonal(vec(coeffDx_bulk)) # coeffDx_bulk is a 2d matrix with shape (grid_u.ny, grid_u.nx), multiplies Bx
    # mat_coeffDy = Diagonal(vec(coeffDy_bulk)) # coeffDx_bulk is a 2d matrix with shape (grid_v.ny, grid_v.nx), multiplies By

    print("\n sizes",size(coeffDx_bulk)) #
    print("\n sizes coeffDy_bulk ",size(coeffDy_bulk)) #

    mat_coeffDx = Diagonal(vec(coeffDu)) # coeffDx_bulk is a 2d matrix with shape (grid_u.ny, grid_u.nx), multiplies Bx
    mat_coeffDy = Diagonal(vec(coeffDv)) # coeffDx_bulk is a 2d matrix with shape (grid_v.ny, grid_v.nx), multiplies By

   
    ni = gp.nx * gp.ny
    nb = 2 * gp.nx + 2 * gp.ny


    mul!(tmp_x, mat_coeffDx * iMx, Bx)
    L = BxT * tmp_x
    mul!(tmp_y, mat_coeffDy * iMy, By)
    L = L .+ ByT * tmp_y

    coeffDu_border = copy(coeffDu)
    coeffDv_border = copy(coeffDv)

    # Interpolate conductivity at center of control volumes for potential gradient at the border
    # interpolate_scalar_to_staggered_u_v_grids_at_border!(num,grid,coeffD,coeffDu,coeffDv)

    interpolate_scalar_to_staggered_u_v_grids_at_border_test!(num,gp,coeffD,coeffDu_border,coeffDv_border)

    coeffDx_border = veci(coeffDu_border,gu)
    coeffDy_border = veci(coeffDv_border,gv)

    mat_coeffDx_b = Diagonal(vec(coeffDx_border)) 
    mat_coeffDy_b = Diagonal(vec(coeffDy_border))

    # mat_coeffDx_b = Diagonal(vec(coeffDu_border)) 
    # mat_coeffDy_b = Diagonal(vec(coeffDv_border))

    # for j in 1:gp.ny
    # #  for i in 1:gp.nx+1

    # #  end
    # print("\n iMx_b * Hx_b j ", (iMx_b * Hx_b)[j,:])
    # end

    # for j in 1:gp.ny+1
    # #  for i in 1:gp.nx+1

    # #  end
    # print("\n iMy_b * Hy_b j ", (iMy_b * Hy_b)[j,:])
    # end

    # print("\n iMx_b * Hx_b ", iMx_b * Hx_b)

    # print("\n iMy_b * Hy_b ", iMy_b  * Hy_b)

    print("\n iMx_b * Hx_b ", iMx_b * Hx_b * coeffD_borders)

    coeffD_borders_test = copy(coeffD_borders)
    coeffD_borders_test .=1.0

    print("\n iMx_b * Hx_b ", mat_coeffDx_b * iMx_b * Hx_b * coeffD_borders_test)

    print("\n iMy_b * Hy_b ", iMy_b  * Hy_b * coeffD_borders)

    print("\n iMy_b * Hy_b ", mat_coeffDy_b * iMy_b  * Hy_b * coeffD_borders_test) #ones(nb))

    # for j in 1:gp.ny
    #     #  for i in 1:gp.nx+1

    #     #  end
    #     # print("\n iMx_b * Hx_b j ", (iMx_b * Hx_b)[j,:])
    #     print("\n iMx_b * Hx_b ", (iMx_b * Hx_b * coeffD_borders)[j,:])

    #     print("\n iMx_b * Hx_b ", (mat_coeffDx_b * iMx_b * Hx_b * ones(nb))[j,:])
    # end


    # bc_L_b = (BxT * mat_coeffDx_b * iMx_b * Hx_b .+ ByT * mat_coeffDy_b * iMy_b  * Hy_b)

    print("\n iMx_b * Hx_b ", BxT * mat_coeffDx_b * iMx_b * Hx_b * coeffD_borders_test)

    print("\n iMy_b * Hy_b ",  ByT * iMy_b  * Hy_b * coeffD_borders)

    print("\n iMy_b * Hy_b ",  ByT * mat_coeffDy_b * iMy_b  * Hy_b * coeffD_borders_test) #ones(nb))


    test_error = maximum(BxT * mat_coeffDx_b * iMx_b * Hx_b * coeffD_borders_test)

    print("\n test_error ",test_error)

    @test test_error ≈ 10.0 atol=test_tolerance #we should have the average conductivity:
    #(9+1)/2, not 9 multiplied by the Laplacian coeff 2

    mul!(tmp_x, mat_coeffDx * iMx, Bx)
    L = BxT * tmp_x
    # mul!(tmp_y, mat_coeffDy * iMy, By)
    # L = L .+ ByT * tmp_y

    # we check only the x contribution
    L = BxT * mat_coeffDx * iMx * Bx 

    print("\n L ", L[1,:])

    print("\n L ", L[1,1])

    @testset "bulk interpolation" begin

        @test L[1,1] ≈ -3.0 atol=test_tolerance

    end

    print("\n L ", L[1,1+grid.ny])


    print("\n L ", L[10,:])


end

#endregion interpolate_scalar_to_staggered_u_v_grids_at_border!


#region gradient
printstyled(color=:green, @sprintf "\n gradient" )

#Initialize liquid phase
x_centroid = gp.x .+ getproperty.(gp.LS[1].geoL.centroid, :x) .* gp.dx
y_centroid = gp.y .+ getproperty.(gp.LS[1].geoL.centroid, :y) .* gp.dy


if verbosity > 2
    print("\n x_centroid ", x_centroid[div(gp.ny,2),:],"\n")
    print("\n gp.dx ", gp.dx[div(gp.ny,2),:],"\n")
end

x_bc = gp.x .+ getproperty.(gp.LS[1].mid_point, :x) .* gp.dx
y_bc = gp.y .+ getproperty.(gp.LS[1].mid_point, :y) .* gp.dy

#Initialize bulk value
vec1(phL.TD,gp) .= vec(ftest.(x_centroid,y_centroid))
#Initialize interfacial value
vec2(phL.TD,gp) .= vec(ftest.(x_bc,y_bc))


vecb_L(phL.TD, gp) .= ftest.(gp.x[:,1] .- gp.dx[:,1] ./ 2.0, gp.y[:,1])
vecb_B(phL.TD, gp) .= ftest.(gp.x[1,:], gp.y[1,:] .- gp.dy[1,:] ./ 2.0)
vecb_R(phL.TD, gp) .= ftest.(gp.x[:,end] .+ gp.dx[:,1] ./ 2.0, gp.y[:,1])
vecb_T(phL.TD, gp) .= ftest.(gp.x[1,:], gp.y[end,:] .+ gp.dy[1,:] ./ 2.0)


#Initialize solid phase


x_centroid = gp.x .+ getproperty.(gp.LS[1].geoS.centroid, :x) .* gp.dx
y_centroid = gp.y .+ getproperty.(gp.LS[1].geoS.centroid, :y) .* gp.dy

x_bc = gp.x .+ getproperty.(gp.LS[1].mid_point, :x) .* gp.dx
y_bc = gp.y .+ getproperty.(gp.LS[1].mid_point, :y) .* gp.dy

vec1(phS.TD,gp) .= vec(ftest.(x_centroid,y_centroid))
vec2(phS.TD,gp) .= vec(ftest.(x_bc,y_bc))


vecb_L(phS.TD, gp) .= ftest.(gp.x[:,1] .- gp.x[:,1] ./ 2.0, gp.y[:,1])
vecb_B(phS.TD, gp) .= ftest.(gp.x[1,:], gp.y[1,:] .- gp.y[1,:] ./ 2.0)
vecb_R(phS.TD, gp) .= ftest.(gp.x[:,end] .+ gp.x[:,1] ./ 2.0, gp.y[:,1])
vecb_T(phS.TD, gp) .= ftest.(gp.x[1,:], gp.y[end,:] .+ gp.y[1,:] ./ 2.0)

phL.T .= reshape(veci(phL.TD,gp,1), gp)
phS.T .= reshape(veci(phS.TD,gp,1), gp)

# phL.T .= ftest.(gp.x,gp.y)

# phS.T .= ftest.(gp.x,gp.y)

LS   = gp.LS
LS_u = gu.LS
LS_v = gv.LS


# compute_grad_T_x!(num,gp, gu, phL, op.opC_pL)
# compute_grad_T_x!(num,gp, gu, phS, op.opC_pS)


grad_x = zeros(gu)
grad_y = zeros(gv)

# compute_grad_T_x_T_y_array!(num.nLS, gp, gu, gv, op.opC_pL, grad_x, grad_y, phL.TD)
compute_grad_T_x_T_y_array_test!(num.nLS, gp, gu, gv, op.opC_pL, grad_x, grad_y, phL.TD)

# compute_grad_T_y_array!(num_LS, gp, gv, op.opC_pL, grad_y, TD)


# grad_analytical = ftest_1.(
#     x_centroid,
#     y_centroid
# )

x_centroid_u = gu.x .+ getproperty.(gu.LS[1].geoS.centroid, :x) .* gu.dx
y_centroid_u = gu.y .+ getproperty.(gu.LS[1].geoS.centroid, :y) .* gu.dy

print("\nx_centroid_u ",x_centroid_u)
print("\ny_centroid_u ",y_centroid_u)

x_centroid_v = gv.x .+ getproperty.(gv.LS[1].geoS.centroid, :x) .* gv.dx
y_centroid_v = gv.y .+ getproperty.(gv.LS[1].geoS.centroid, :y) .* gv.dy

print("\nx_centroid_v ",x_centroid_v)
print("\ny_centroid_v ",y_centroid_v)

x_centroid_u = gu.x .+ getproperty.(gu.LS[1].geoL.centroid, :x) .* gu.dx
y_centroid_u = gu.y .+ getproperty.(gu.LS[1].geoL.centroid, :y) .* gu.dy

print("\nx_centroid_u ",x_centroid_u)
print("\ny_centroid_u ",y_centroid_u)

x_centroid_v = gv.x .+ getproperty.(gv.LS[1].geoL.centroid, :x) .* gv.dx
y_centroid_v = gv.y .+ getproperty.(gv.LS[1].geoL.centroid, :y) .* gv.dy

print("\nx_centroid_v ",x_centroid_v)
print("\ny_centroid_v ",y_centroid_v)

grad_analytical = ftest_1.(
    x_centroid_u,
    y_centroid_u
)

j = div(n,2)

print("\n grad_analytical ",grad_analytical[j,:],"\n")
print("\n D grad_x[j,:]",grad_x[j,:],"\n")
print("\n D grad_x[j,:]",grad_x[j,:],"\n")

@testset "Simple gradient test x" begin
    print("\n grad_x[j,:]",grad_x[j,:]," len ",length(grad_x[j,:]))
    test_error = maximum(abs.(grad_x[j,:] .- 1.0))
    @test test_error ≈ 0.0 atol=test_tolerance #first and last values are not gradients and correspond to left and right walls
    # @test grad_x[j,1:end-1] .≈ 1.0 atol=test_tolerance #first and last values are not gradients and correspond to left and right walls
end

# print("\n grad_x at x j",grad_x[j,:]," len ",length(grad_x[j,:]))
# print("\n gu.y at x j",gu.y[j,:]," len ",length(gu.y[j,:]))


@testset "Simple gradient test y" begin
    print("\n grad_y[j,:]",grad_y[j,:]," len ",length(grad_y[j,:]))

    test_error = maximum(abs.(grad_y[j,:] .- 0.0))
    @test test_error ≈ 0.0 atol=test_tolerance #first and last values are not gradients and correspond to bottom and top walls
end

#endregion gradient


#region gradient staggered

@testset "Simple gradient staggered" begin

    printstyled(color=:red, @sprintf "\n compute_grad_T_x_T_y_array_u_v_capacities! \n") 

    

    # x_centroid_u = gu.x .+ getproperty.(gu.LS[1].geoS.centroid, :x) .* gu.dx
    # y_centroid_v = gv.y .+ getproperty.(gv.LS[1].geoS.centroid, :y) .* gv.dy
    #geoS ?

    # print("\n x_centroid_u ",x_centroid_u[1,:])
    # print("\n y_centroid_v ",y_centroid_v[:,1])


    grad_x = zeros(gu)
    grad_y = zeros(gv)
    
    # compute_grad_T_x_T_y_array!(num.nLS, gp, gu, gv, op.opC_pL, grad_x, grad_y, phL.TD)
    compute_grad_T_x_T_y_array_u_v_capacities!(num, gp, gu, gv, op.opC_uL,op.opC_vL, grad_x, grad_y, phL.TD)
    
    # compute_grad_T_y_array!(num_LS, gp, gv, op.opC_pL, grad_y, TD)
    
    
    # grad_analytical = ftest_1.(
    #     x_centroid,
    #     y_centroid
    # )
    
    # L
    x_centroid_u = gu.x .+ getproperty.(gu.LS[1].geoL.centroid, :x) .* gu.dx
    y_centroid_u = gu.y .+ getproperty.(gu.LS[1].geoL.centroid, :y) .* gu.dy
    
    grad_analytical = ftest_1.(
        x_centroid_u,
        y_centroid_u
    )
    
    j = div(n,2)
    print("\n length grad_analytical ",length(grad_analytical[j,:]),"\n")

    print("\n grad_analytical ",grad_analytical[j,:],"\n")

    print("\n length grad_x ",length(grad_x[j,:]),"\n")

    print("\n D grad_x[j,:]",grad_x[j,:],"\n")
    print("\n D grad_x[j,:]",grad_x[j,:],"\n")
    
    print("\n D grad_y[:,i]",grad_x[:,j],"\n")

    print("\n vecb_L",vecb_L(phL.TD,gp),"\n")



    @testset "Simple gradient test x" begin
        print("\n grad_x[j,:]",grad_x[j,:]," len ",length(grad_x[j,:]))
        test_error = maximum(abs.(grad_x[j,:] .- 1.0))
        @test test_error ≈ 0.0 atol=test_tolerance #first and last values are not gradients and correspond to left and right walls
        # @test grad_x[j,1:end-1] .≈ 1.0 atol=test_tolerance #first and last values are not gradients and correspond to left and right walls
    end
    
    # print("\n grad_x at x j",grad_x[j,:]," len ",length(grad_x[j,:]))
    # print("\n gu.y at x j",gu.y[j,:]," len ",length(gu.y[j,:]))
    
    
    @testset "Simple gradient test y" begin
        print("\n grad_y[j,:]",grad_y[j,:]," len ",length(grad_y[j,:]))
    
        test_error = maximum(abs.(grad_y[j,:] .- 0.0))
        @test test_error ≈ 0.0 atol=test_tolerance #first and last values are not gradients and correspond to bottom and top walls
    end

end
#endregion gradient staggered


#region laplacian
printstyled(color=:red, @sprintf "\n test laplacian \n") 

@testset "Laplacian" begin
    #Allocations for scalar grid
    ni = gv.nx * gv.ny
    nb = 2 * gv.nx + 2 * gv.ny
    nt = (num.nLS + 1) * ni + nb

    AvL =  spzeros(nt, nt) 
    tmp_vec_1D_v = fnzeros(gv,num)
    tmp_vec_1D_v0 = fnzeros(gv,num)

    BC_int = [WallNoSlip()]

    periodic_x = false
    periodic_y = false

    update_all_ls_data(num, gp, gu, gv, BC_int, periodic_x, periodic_y, false) #not periodic

    # geoL = [grid.LS[iLS].geoL for iLS in 1:num._nLS]
    # geo_uL = [grid_u.LS[iLS].geoL for iLS in 1:num._nLS]
    # geo_vL = [grid_v.LS[iLS].geoL for iLS in 1:num._nLS]
    # Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, Lum1_L, bc_Lum1_L, bc_Lum1_b_L, Lvm1_L, bc_Lvm1_L, bc_Lvm1_b_L = set_matrices!(
    #     num, grid, geoL, grid_u, geo_uL, grid_v, geo_vL,
    #     op.opC_pL, op.opC_uL, op.opC_vL, periodic_x, periodic_y
    # )

    laps = set_matrices!(num, gp, [gp.LS[1].geoL], gu, [gu.LS[1].geoL], gv, [gv.LS[1].geoL], 
    op.opC_pL, op.opC_uL, op.opC_vL, periodic_x, periodic_y) #not periodic
    # Lp, bc_Lp, bc_Lp_b, Lu, bc_Lu, bc_Lu_b, Lv, bc_Lv, bc_Lv_b = laps
    Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, Lum1_L, bc_Lum1_L, bc_Lum1_b_L, Lvm1_L, bc_Lvm1_L, bc_Lvm1_b_L = laps

    compute_divergence_test!(num, 
    # grid, 
    # grid_u, 
    gv, 
    op,
    AvL, 
    # rhs_scal,
    # tmp_vec_p, #a0
    tmp_vec_1D_v0,
    tmp_vec_1D_v,
    Lvm1_L, 
    bc_Lvm1_L, 
    bc_Lvm1_b_L
    # tmp_vec_u0,
    # tmp_vec_v0,
    # tmp_vec_1D,
    # ls_advection
    )

    dx0 = L0/n
    volume0 = dx0^2
    volume02 = volume0/2
    print("\n volume ",volume0," wall ",volume02,"op.opC_vL.M.diag[5,1] ",op.opC_vL.M.diag[5,1])
    # @test rhs_vec1[5,1] ≈ volume02 atol = test_tolerance
    # @test op.opC_vL.M.diag ≈ volume02 atol = test_tolerance
    @test op.opC_vL.M.diag[5,1] ≈ volume0*3/4 atol = test_tolerance
end
#endregion laplacian

#region mask

@testset "mask_1D" begin
    # TODO add bubble wall to check mask
    mask_1D = fnzeros(gp,num)
    compute_mask_1D!(num,gp,mask_1D)
    # mask_1D ≈ 0 atol=test_tolerance
    @test mask_1D ≈ 0 atol=test_tolerance

end 

#endregion mask

#region wall gradient
#TODO

#endregion wall gradient



#region mass flux

printstyled(color=:red, @sprintf "\n testing divergence \n") 

integrate_mass_flux_over_interface_3_no_writing(num,gp,op.opC_pL,
phL.TD,mass_flux_vec1,mass_flux_vecb,mass_flux_veci,tmp_vec_p,tmp_vec_p0,tmp_vec_p1,mass_flux)
# @test sum(mass_flux) == 0 
# @test sum(mass_flux) ≈ 0 atol=test_tolerance


integrate_mass_flux_over_interface_3(num,gp,op.opC_pL,
phL.TD,
mass_flux_vec1,
mass_flux_vecb,mass_flux_veci, tmp_vec_p, tmp_vec_p0, tmp_vec_p1, mass_flux,num.index_phase_change)

# @test sum(mass_flux) ≈ 0 atol=test_tolerance
print("\n num.index_phase_change ",num.index_phase_change,"\n")
printstyled(color=:red, @sprintf "\n testing divergence \n") 



integrate_mass_flux_over_interface(num,gp,op.opC_pL,
phL.TD,
mass_flux_vec1,
mass_flux_vecb,mass_flux_veci, tmp_vec_p, tmp_vec_p0, tmp_vec_p1, mass_flux,num.index_phase_change)




print("\n flux j line ", mass_flux_vec1[j,:], "\n")


#zero because with H interface 
@testset "mass flux 0" begin
    @test mass_flux[1,1] ≈ 0.0 atol = test_tolerance
end

@testset "mass flux 0" begin
    @test mass_flux[1,gp.nx] ≈ 0.0 atol = test_tolerance
end

# @testset "Gradient x-component" begin
#     @test grad_analytical ≈ grad_x atol = test_tolerance
# end

# compute_grad_T_y!(num,gp, gv, phL, op.opC_pL)
# compute_grad_T_y!(num,gp, gv, phS, op.opC_pS)

# printstyled(color=:red, @sprintf "\n grad y %f %f %f\n" norm(phL.v) minimum(phL.v) maximum(phL.v))
# printstyled(color=:red, @sprintf "\n grad y %f %f %f\n" norm(phS.v) minimum(phS.v) maximum(phS.v))

# @testset "Gradient y-component" begin
#     @test maximum(abs.(phL.v)) ≈ 0.0 atol = test_tolerance
# end

#TODO test mass flux 

#endregion mass flux


#region diffusion u v 

    function set_cutcell_matrices_test!(num, grid, geo, geo_p, opC, periodic_x, periodic_y)
        @unpack nx, ny, ind = grid
        @unpack AxT, AyT, Bx, By, BxT, ByT, Hx, Hy, HxT, HyT, M, iMx, iMy, χ = opC

        M.diag .= vec(geo[end].dcap[:,:,5])
        Mx = zeros(ny,nx+1)
        for II in ind.all_indices
            Mx[II] = geo[end].dcap[II,8]
            pII = lexicographic(II, grid.ny)
            iMx.diag[pII] = inv_weight_eps(num,Mx[II])
        end
        for II in ind.b_right[1]
            Mx[δx⁺(II)] = geo[end].dcap[II,10]
            pII = lexicographic(δx⁺(II), grid.ny)
            iMx.diag[pII] = inv_weight_eps(num,Mx[δx⁺(II)])
        end

        print("\n Mx ",size(Mx))
        print("\n Mx ",Mx)



        My = zeros(ny+1,nx)
        for II in ind.all_indices
            My[II] = geo[end].dcap[II,9]
            pII = lexicographic(II, grid.ny + 1)
            iMy.diag[pII] = inv_weight_eps(num,My[II])
        end
        for II in ind.b_top[1]
            My[δy⁺(II)] = geo[end].dcap[II,11]
            pII = lexicographic(δy⁺(II), grid.ny + 1)
            iMy.diag[pII] = inv_weight_eps(num,My[δy⁺(II)])
        end   

        # Discrete gradient and divergence operators
        divergence_A!(grid, AxT, AyT, geo[end].dcap, ny, ind.all_indices, periodic_x, periodic_y)
        divergence_B!(BxT, ByT, geo[end].dcap, ny, ind.all_indices)

        mat_assign!(Bx, sparse(-BxT'))
        mat_assign!(By, sparse(-ByT'))

        # Matrices for BCs
        for iLS in 1:num.nLS
            bc_matrix!(grid, Hx[iLS], Hy[iLS], geo[iLS].dcap, geo_p[iLS].dcap, ny, ind.all_indices)

            mat_assign_T!(HxT[iLS], sparse(Hx[iLS]'))
            mat_assign_T!(HyT[iLS], sparse(Hy[iLS]'))

            periodic_bcs!(grid, Bx, By, Hx[iLS], Hy[iLS], periodic_x, periodic_y)

            χx = (geo[iLS].dcap[:,:,3] .- geo[iLS].dcap[:,:,1]) .^ 2
            χy = (geo[iLS].dcap[:,:,4] .- geo[iLS].dcap[:,:,2]) .^ 2
            χ[iLS].diag .= sqrt.(vec(χx .+ χy))
        end

        mat_assign!(BxT, sparse(-Bx'))
        mat_assign!(ByT, sparse(-By'))

        return nothing
    end
@testset "diffusion_u_v" begin

    printstyled(color=:magenta, @sprintf "\n testing diffusion_u_v \n") 


    @unpack Bx, By, Hx, Hy, HxT, HyT, χ, M, iMx, iMy, Hx_b, Hy_b, HxT_b, HyT_b, iMx_b, iMy_b, iMx_bd, iMy_bd, χ_b = op.opC_pL
    @unpack BxT, ByT,tmp_x, tmp_y = op.opC_pL

    opC_p = op.opC_pL
    opC_u = op.opC_uL
    opC_v = op.opC_vL


    grid = gp
    grid_u = gu
    grid_v = gv
   
    ni = grid.nx * grid.ny
    nb = 2 * grid.nx + 2 * grid.ny

    coeffD = fnones(grid,num)


    #region interpolate at border

    print("\n opC_u.iMx_b", size(opC_u.iMx_b)," opC_u.BxT ",size(opC_u.BxT))
    print("\n opC_u.iMy_b", size(opC_u.iMy_b)," opC_u.ByT ",size(opC_u.ByT))

    print("\n opC_v.iMx_b", size(opC_v.iMx_b)," opC_v.BxT ",size(opC_v.BxT))
    print("\n opC_v.iMy_b", size(opC_v.iMy_b)," opC_v.ByT ",size(opC_v.ByT))


    laps = set_matrices!(num, gp, [gp.LS[1].geoL], gu, [gu.LS[1].geoL], gv, [gv.LS[1].geoL], 
    op.opC_pL, op.opC_uL, op.opC_vL, false, false)
    Lp, bc_Lp, bc_Lp_b, Lu, bc_Lu, bc_Lu_b, Lv, bc_Lv, bc_Lv_b = laps
    
    print("\n bc_Lu_b", size(bc_Lu_b)," bc_Lv_b ",size(bc_Lv_b))

    j = Int(grid.ny /2) 
    i = 1
    II = CartesianIndex(j,i)
    pII = lexicographic(II,grid.ny)

    # test_index = Int(grid.nx /2)
    test_index = pII
    print("\n Lu ", Lu[test_index,:])

    print("\n bc_Lu_b ", bc_Lu_b[test_index,:])

    print("\n bc_Lv_b ", bc_Lv_b[test_index,:])

  

    print("\nprint scal \n")

    print("\n Lp ", Lp[test_index,:])

    print("\n bc_Lp_b ", bc_Lp_b[test_index,:])

    # print("\nprint matrix \n")

    # # display(test_bcL)
    # Base.print_matrix(stdout,test_bcL)


    @testset " Laplacian u" begin
        testval = -(1 + 4/3 + 4) 
        print("\n -(1 + 4/3 + 4) ",testval,Lu[test_index,test_index])
        print("\n test Lu")
        @test Lu[test_index,test_index] ≈ testval atol=test_tolerance # 1 because half cell 0.5
        
        print("\n test bc_Lu_b")

        @test bc_Lu_b[test_index,test_index] ≈ 4.0 atol=test_tolerance # 1 because half cell 0.5

        printstyled(color=:magenta, @sprintf "\n testing bottom of the domain \n") 
        II = CartesianIndex(1,2)
        print("\n II ",II)
        pII = lexicographic(II,grid.ny)
        print("\n Lu ", Lu[pII,:])
        print("\n bc_Lu_b ", bc_Lu_b[pII,:])


        II = CartesianIndex(1,Int(grid.nx /2))
        print("\n II ",II)
        pII = lexicographic(II,grid.ny)
        print("\n Lu ", Lu[pII,:])
        print("\n bc_Lu_b ", bc_Lu_b[pII,:])



        # @test Lu[pII,pII] ≈ testval atol=test_tolerance # 1 because half cell 0.5
        # print("\n test bc_Lu_b")
        # @test bc_Lu_b[pII,pII] ≈ 4.0 atol=test_tolerance # 1 because half cell 0.5
        
    end

    @testset " Laplacian u with coefficient" begin


        viscosity_coeff_for_u_border_x = zeros(grid.ny, grid.nx+2) # copy(coeffDu)
        viscosity_coeff_for_u_border_y = zeros(grid.ny+1, grid.nx+1) # copy(coeffDv)

        viscosity_coeff_for_v_border_x = zeros(grid.ny+1, grid.nx+1) #copy(coeffDu)
        viscosity_coeff_for_v_border_y = zeros(grid.ny+2, grid.nx) #copy(coeffDv)
        print("\n grid_u size ")

        print("\n grid.ind.b_left[1][1] grid.ind.b_right[1][1]",grid_u.ind.b_left[1][1], " ",grid_u.ind.b_right[1][1] )
        print("\n grid.ind.b_bottom[1][1] grid.ind.b_top[1][1]",grid_u.ind.b_bottom[1][1], " ",grid_u.ind.b_top[1][1] )


        printstyled(color=:magenta, @sprintf "\n testing variable coefficient \n") 

        viscosity_coeff_for_du_dx = zeros(grid_u.ny,grid_u.nx+1)
        viscosity_coeff_for_du_dx[:,2:grid_u.nx] .= 1.0 #volume_fraction #grid.LS[end].geoL.cap[:,:,5]

        #TODO contact angle change  viscosity_coeff_for_du_dx[:,1] and at end (not interpolating right now)
        viscosity_coeff_for_du_dx[:,1] = viscosity_coeff_for_du_dx[:,2]
        viscosity_coeff_for_du_dx[:,end] = viscosity_coeff_for_du_dx[:,end-1]

        display(viscosity_coeff_for_du_dx)

        viscosity_coeff_for_u_border_x = viscosity_coeff_for_du_dx
        
        viscosity_coeff_for_u_border_y = zeros(grid.ny+1, grid.nx+1)

        viscosity_coeff_for_u_border_y .= 1.0


        viscosity_coeff_for_v_border_x = ones(grid.ny+1, grid.nx+1)
        viscosity_coeff_for_v_border_y = ones(grid.ny+2, grid.nx)

        # print("\n grid.ind.b_left[1][1] grid.ind.b_right[1][1]",size(coeffDu), " ",size(coeffDv) )

        # Interpolate at the border
        # interpolate_scalar_to_staggered_u_v_grids_at_border!(num,grid_u,coeffD,viscosity_coeff_for_u_border_x,viscosity_coeff_for_u_border_y)

        # display(viscosity_coeff_for_u_border_x)
        # display(viscosity_coeff_for_u_border_y)


        # interpolate_scalar_to_staggered_u_v_grids_at_border!(num,grid_v,coeffD,viscosity_coeff_for_v_border_x,viscosity_coeff_for_v_border_y)

        # viscosity_coeff_for_u_border_x = veci(viscosity_coeff_for_u_border_x,grid_u) #should be of size (ny, nx+2)
        # viscosity_coeff_for_u_border_y = veci(viscosity_coeff_for_u_border_y,grid_v) #should be of size (ny+1, nx+1)

        # viscosity_coeff_for_v_border_x = veci(viscosity_coeff_for_v_border_x,grid_u) #should be of size (ny+1, nx+1)
        # viscosity_coeff_for_v_border_y = veci(viscosity_coeff_for_v_border_y,grid_v) #should be of size (ny+2, nx)

        diag_viscosity_coeff_for_u_border_x = Diagonal(vec(viscosity_coeff_for_u_border_x)) 
        diag_viscosity_coeff_for_u_border_y = Diagonal(vec(viscosity_coeff_for_u_border_y)) 

        diag_viscosity_coeff_for_v_border_x = Diagonal(vec(viscosity_coeff_for_v_border_x)) 
        diag_viscosity_coeff_for_v_border_y = Diagonal(vec(viscosity_coeff_for_v_border_y)) 

        # print("\n size(diag_viscosity_coeff_for_u_border_x) ",size(diag_viscosity_coeff_for_u_border_x), " ",size(opC_u.BxT)," ",size(opC_u.iMx_b))

        bc_Lu_b = (opC_u.BxT * diag_viscosity_coeff_for_u_border_x * opC_u.iMx_b * opC_u.Hx_b .+ opC_u.ByT * diag_viscosity_coeff_for_u_border_y * opC_u.iMy_b * opC_u.Hy_b)

        bc_Lv_b = (opC_v.BxT * diag_viscosity_coeff_for_v_border_x * opC_v.iMx_b * opC_v.Hx_b .+ opC_v.ByT * diag_viscosity_coeff_for_v_border_y * opC_v.iMy_b * opC_v.Hy_b)


        print("\nprint matrix \n")
       
        test_bcL = opC_u.BxT * opC_u.iMx_b * opC_u.Hx_b .+ opC_u.ByT * opC_u.iMy_b * opC_u.Hy_b

        # print("\nprint matrix ", test_bcL[1,:] )
        # print("\nprint matrix ", test_bcL[2,:] )

        testyindex = Int(grid.ny /2)
        
        II = CartesianIndex( testyindex, 1)
        print("\n II ",II)
        pII = lexicographic(II,grid.ny)
        print("\n Lu ", Lu[pII,:])
        print("\n bc_Lu_b ", bc_Lu_b[pII,:])

        II = CartesianIndex( testyindex-1, 1)
        print("\n II ",II)
        pII = lexicographic(II,grid.ny)
        print("\n Lu ", Lu[pII,:])
        print("\n bc_Lu_b ", bc_Lu_b[pII,:])

        II = CartesianIndex( testyindex+1, 1)
        print("\n II ",II)
        pII = lexicographic(II,grid.ny)
        print("\n Lu ", Lu[pII,:])
        print("\n bc_Lu_b ", bc_Lu_b[pII,:])


        viscosity_coeff_for_u_border_x .= 0.0
        
        viscosity_coeff_for_u_border_y .= 0.0

        # II = CartesianIndex( testyindex, 1)
        viscosity_coeff_for_u_border_x[testyindex,1] = 1.0 #TODO why 2, and 4 ? staggered +1 -1

        # for j in 1:size(viscosity_coeff_for_u_border_x,1)
        #     for i in 1:size(viscosity_coeff_for_u_border_x,2)

        #         print("\n j i ",j," ",i," ",j*1000+i)

        #         viscosity_coeff_for_u_border_x[j,i] = (j*1000+i)/4

        #     end

        # end

        # print("\n opC_u.iMx_b * opC_u.Hx_b")
        # display(opC_u.iMx_b * opC_u.Hx_b)

        # print("\n opC_u.BxT")
        # display(opC_u.BxT)
        
        # viscosity_coeff_for_u_border_y .= 0.0

        diag_viscosity_coeff_for_u_border_x = Diagonal(vec(viscosity_coeff_for_u_border_x)) 
        diag_viscosity_coeff_for_u_border_y = Diagonal(vec(viscosity_coeff_for_u_border_y)) 

        bc_Lu_b = (opC_u.BxT * diag_viscosity_coeff_for_u_border_x * opC_u.iMx_b * opC_u.Hx_b .+ opC_u.ByT * diag_viscosity_coeff_for_u_border_y * opC_u.iMy_b * opC_u.Hy_b)

        II = CartesianIndex( testyindex, 1)
        print("\n II ",II)
        pII = lexicographic(II,grid.ny)
        print("\n Lu ", Lu[pII,:])
        print("\n bc_Lu_b ", bc_Lu_b[pII,:])

        II = CartesianIndex( testyindex-1, 1)
        print("\n II ",II)
        pII = lexicographic(II,grid.ny)
        print("\n Lu ", Lu[pII,:])
        print("\n bc_Lu_b ", bc_Lu_b[pII,:])

        II = CartesianIndex( testyindex+1, 1)
        print("\n II ",II)
        pII = lexicographic(II,grid.ny)
        print("\n Lu ", Lu[pII,:])
        print("\n bc_Lu_b ", bc_Lu_b[pII,:])

        ###################################

        testval = -(1 + 4/3 + 4) 
        print("\n -(1 + 4/3 + 4) ",testval,Lu[test_index,test_index])
        print("\n test Lu")
        @test Lu[test_index,test_index] ≈ testval atol=test_tolerance # 1 because half cell 0.5
        
        print("\n test bc_Lu_b")

        @test bc_Lu_b[test_index,test_index] ≈ 4.0 atol=test_tolerance # 1 because half cell 0.5


        printstyled(color=:magenta, @sprintf "\n testing variable coefficient \n") 


    end



    #endregion interpolate at border



    #region scalar_all_nodes



    create_2D_grid = zeros(grid.ny+2,grid.nx+2)

    x_centroid = gp.x .+ getproperty.(gp.LS[1].geoS.centroid, :x) .* gp.dx
    y_centroid = gp.y .+ getproperty.(gp.LS[1].geoS.centroid, :y) .* gp.dy

    create_2D_grid[2:grid.ny+1,2:grid.nx+1] = x_centroid #grid.x

    x_bc_left = gp.x[:,1] .- gp.dx[:,1] ./ 2.0

    y_bc_bottom = gp.y[1,:] .- gp.dy[1,:] ./ 2.0

    y_bc_top = gp.y[end,:] .+ gp.dy[end,:] ./ 2.0

    x_bc_right = gp.x[:,end] .+ gp.dx[:,end] ./ 2.0

    # create_2D_grid[1,2:grid.nx] = create_2D_grid[2,2:grid.nx]

    # create_2D_grid[end,2:grid.nx] = create_2D_grid[end-1,2:grid.nx]

    # display(create_2D_grid)


    create_2D_grid[2:grid.ny+1,1] = x_bc_left

    create_2D_grid[2:grid.ny+1,end] = x_bc_right

    create_2D_grid[1,:] = create_2D_grid[2,:]

    create_2D_grid[end,:] = create_2D_grid[end-1,:]

    printstyled(color=:magenta, @sprintf "\n create_2D_grid all scalar nodes \n") 
    
    # display(create_2D_grid)

    create_2D_grid = create_2D_grid_x(gp)

    display(create_2D_grid)


    #endregion scalar_all_nodes


    #region viscosity_coeff_for_du_dx
    # Mx is the volume of the control volume associated to the gradient, between two nodes
    # for grid_u: staggered in x, between (border + bulk) nodes there are n+2 control volumes 

    viscosity_coeff_for_du_dx = zeros(grid_u.ny,grid_u.nx+1)

    viscosity_coeff_for_du_dx[:,2:grid_u.nx] = grid.x #grid.x*1000 +grid.y

    x_centroid_u = gu.x .+ getproperty.(gu.LS[1].geoL.centroid, :x) .* gu.dx
    y_centroid_u = gu.y .+ getproperty.(gu.LS[1].geoL.centroid, :y) .* gu.dy

    viscosity_coeff_for_du_dx[:,1] = (gu.x[:,1] + x_centroid_u[:,1]) / 2
    viscosity_coeff_for_du_dx[:,grid_u.nx+1] = (x_centroid_u[:,grid_u.nx] + gu.x[:,grid_u.nx]) / 2

    #endregion viscosity_coeff_for_du_dx

    print("\ngridu x ",grid_u.x[1,:])
    print("\ngrid x ",grid.x[1,:])
    print("\nx_centroid_u ",x_centroid_u)
    print("\ny_centroid_u ",y_centroid_u)
    

    printstyled(color=:magenta, @sprintf "\ninterp for du/dx \n")

    display(viscosity_coeff_for_du_dx)

    
    viscosity_coeff_for_du_dx = zeros(grid_u.ny,grid_u.nx+1)
    viscosity_coeff_for_du_dx[:,2:grid_u.nx] = grid.LS[end].geoL.cap[:,:,5]

    viscosity_coeff_for_du_dx[:,1] = viscosity_coeff_for_du_dx[:,2]
    viscosity_coeff_for_du_dx[:,end] = viscosity_coeff_for_du_dx[:,end-1]
    #TODO contact angle change  viscosity_coeff_for_du_dx[:,1] and at end

    display(viscosity_coeff_for_du_dx)


    # print("\n size ",size(viscosity_coeff_for_du_dx))



    # laps = set_matrices!(num, gp, [gp.LS[1].geoL], gu, [gu.LS[1].geoL], gv, [gv.LS[1].geoL], 
    # op.opC_pL, op.opC_uL, op.opC_vL, periodic_x, periodic_y)
    
    # geo = [gp.LS[1].geoL]
    # geo_u = [gu.LS[1].geoL]
    # geo_v = [gv.LS[1].geoL]

    # set_cutcell_matrices_test!(num, grid, geo, geo, opC_p, false, false)

    # set_cutcell_matrices_test!(num, grid_u, geo_u, geo, opC_u, false, false)

    # set_cutcell_matrices_test!(num, grid_v, geo_v, geo, opC_v, false, false)

    # print("\n size(opC_u.Bx) ",size(opC_u.Bx)) #from u (n+1) to grad u (n+2)
    # print("\n size(opC_u.iMx) ",size(opC_u.iMx)) #size grad u 
    # print("\n size(mat_coeffD) ",size(diag_viscosity_coeff_for_du_dx)) #
    # print("\n size(tmp_x) ",size(tmp_x))

    # print("\n size(mat_coeffD * opC_u.iMx) ",size(diag_viscosity_coeff_for_du_dx * opC_u.iMx))

    
    diag_viscosity_coeff_for_du_dx = Diagonal(vec(viscosity_coeff_for_du_dx))

    mul!(opC_u.tmp_x, diag_viscosity_coeff_for_du_dx * opC_u.iMx, opC_u.Bx)

    diffusion_bulk_u = opC_u.BxT * opC_u.tmp_x




    printstyled(color=:magenta, @sprintf "\n viscosity_coeff_for_dv_dy \n") 

    viscosity_coeff_for_du_dy = zeros(grid_u.ny+1,grid_u.nx)

    #region Viscosity coefficient for \frac{\partial v}{\partial y}
    #cf test in orientation.jl
    viscosity_coeff_for_dv_dy = zeros(grid_v.ny+1,grid_v.nx)

    viscosity_coeff_for_dv_dy[2:grid_v.ny,:] = grid.y #volume_fraction #grid.LS[end].geoL.cap[:,:,5]

    # x_centroid_v

    viscosity_coeff_for_dv_dy[1,:] = (gv.y[1,:] + y_centroid_v[1,:]) / 2
    viscosity_coeff_for_dv_dy[grid_v.ny+1,:] = (y_centroid_v[grid_v.ny,:] + gv.y[grid_v.ny,:]) / 2


    display(viscosity_coeff_for_dv_dy)

    viscosity_coeff_for_dv_dy[2:grid_v.ny,:] = grid.LS[end].geoL.cap[:,:,5]

    #TODO contact angle change  viscosity_coeff_for_dv_dy[:,1] and at end (not interpolating right now)
    viscosity_coeff_for_dv_dy[1,:] = viscosity_coeff_for_dv_dy[2,:]
    viscosity_coeff_for_dv_dy[end,:] = viscosity_coeff_for_dv_dy[end-1,:]

    display(viscosity_coeff_for_dv_dy)

    # arithmetic average 
    # viscosity_coeff_for_dv_dy .= 2 * ( (num.mu1 - num.mu2) * viscosity_coeff_for_dv_dy  .+ num.mu2 )

    # # display(viscosity_coeff_for_dv_dy)
    # PDI_status = @ccall "libpdi".PDI_multi_expose("viscosity_coeff_for_dv_dy"::Cstring,
    # "viscosity_coeff_for_dv_dy"::Cstring, viscosity_coeff_for_dv_dy::Ptr{Cdouble}, PDI_OUT::Cint,        
    # C_NULL::Ptr{Cvoid})::Cint

    # diag_viscosity_coeff_for_dv_dy = Diagonal(vec(viscosity_coeff_for_dv_dy))

    #endregion Viscosity coefficient for \frac{\partial v}{\partial y}


   

    mu1 = 1
    mu2 = 100





   


    #region interpolate 
    volume_fraction_full = create_2D_grid_volume_fraction(gp)
    grid_x_full_2D = create_2D_grid_x(gp,true,true)
    grid_y_full_2D = create_2D_grid_y(gp,true,true)

    all_grid_u_nodes_2D_x_for_du_dy_interp = create_2D_grid_x(gu,false,true)
    all_grid_u_nodes_2D_y_for_du_dy_interp = create_2D_grid_y(gu,false,true)

    printstyled(color=:magenta, @sprintf "\n grid_x_full_2D \n") 

    display(grid_x_full_2D)
    
    printstyled(color=:magenta, @sprintf "\n grid_y_full_2D \n") 

    display(grid_y_full_2D)

    printstyled(color=:magenta, @sprintf "\n x for viscosity_coeff_for_du_dy \n") 

    display(all_grid_u_nodes_2D_x_for_du_dy_interp)    
    # display(grid_u.x)    
    # display(x_centroid_u)    

    
    printstyled(color=:magenta, @sprintf "\n all_grid_u_nodes_2D_y_for_du_dy_interp \n") 

    display(all_grid_u_nodes_2D_y_for_du_dy_interp)


    for j in 1:grid_u.ny+1
        for i in 1:grid_u.nx

            print("\nvolume_fraction i ",i," j ",j,"\n")

            du_dy_coord_x = all_grid_u_nodes_2D_x_for_du_dy_interp[j,i] 
            # du_dy_coord_x = x_centroid_u[j,i]  #grid_u.x[j,i] 

            du_dy_coord_y = (all_grid_u_nodes_2D_y_for_du_dy_interp[j,i] + all_grid_u_nodes_2D_y_for_du_dy_interp[j+1,i])/2 

            # x1 = grid.x[j-1,i-1]
            # x2 = grid_u.x[j-1,i]
            # y1 = grid.y[j-1,i-1]
            # y2 = grid.y[j,i-1]

            # x1 = grid_x_full_2D[j-1,i-1]
            # x2 = grid_x_full_2D[j-1,i]
            # y1 = grid_y_full_2D[j-1,i-1]
            # y2 = grid_y_full_2D[j,i-1]

            # Q11 = volume_fraction_full[j-1,i-1]
            # Q12 = volume_fraction_full[j,i-1]
            # Q21 = volume_fraction_full[j-1,i]
            # Q22 = volume_fraction_full[j,i]



            x1 = grid_x_full_2D[j,i]
            x2 = grid_x_full_2D[j,i+1]
            y1 = grid_y_full_2D[j,i]
            y2 = grid_y_full_2D[j+1,i+1]

            Q11 = volume_fraction_full[j,i]
            Q12 = volume_fraction_full[j+1,i]
            Q21 = volume_fraction_full[j,i+1]
            Q22 = volume_fraction_full[j+1,i+1]

            volume_fraction_face = bilinear_interpolation(du_dy_coord_x, du_dy_coord_y, x1, y1, x2, y2, Q11, Q12, Q21, Q22)

            printstyled(color=:green, @sprintf "\n i %.5i j %.5i x %.2e y %.2e x1 %.2e y1 %.2e x2 %.2e y2 %.2e Q11 %.2e Q12 %.2e Q21 %.2e Q22 %.2e\n" i j du_dy_coord_x du_dy_coord_y x1 y1 x2 y2 Q11 Q12 Q21 Q22)

            print("\nvolume_fraction i ",i," j ",j," ",volume_fraction_face)
            viscosity_coeff_for_du_dy[j,i] = harmonic_average_one_fluid(mu1,mu2,volume_fraction_face)

        end
    end
    
    printstyled(color=:cyan, @sprintf "\n viscosity_coeff_for_du_dy\n")

    display(viscosity_coeff_for_du_dy)

    #endregion interpolate 

  
    #region old method of interpolation
   
    old_method = false 

    if old_method
        volume_fraction = grid.LS[end].geoL.cap[:,:,5]

        for j in 2:grid_u.ny
            for i in 2:grid_u.nx-1

                du_dy_coord_x = x_centroid_u[j,i] #grid_u.x[j,i]
                du_dy_coord_y = (y_centroid_u[j-1,i] + y_centroid_u[j,i])/2 #(grid_u.y[j,i] + grid_u.y[j+1,i])/2
            
                # print("\n du_dy_coord_x ",du_dy_coord_x," ",du_dy_coord_y," ",x_centroid_u[j,i]," ",y_centroid_u[j,i])
                # test_cooord = bilinear_interpolation(grid, du_dy_coord_x, du_dy_coord_y,grid.x)
                # print("\ntest x i ",i," j ",j," ",test_cooord," ",du_dy_coord_x)
                # @test test_cooord == du_dy_coord_x
                # test_cooord = bilinear_interpolation(grid, du_dy_coord_x, du_dy_coord_y,grid.y)
                # print("\ntest y i ",i," j ",j," ",test_cooord," ",du_dy_coord_y)
                # @test test_cooord == du_dy_coord_y


                volume_fraction_face = bilinear_interpolation(grid, du_dy_coord_x, du_dy_coord_y,volume_fraction)
                # print("\nvolume_fraction i ",i," j ",j," ",volume_fraction_face)
                viscosity_coeff_for_du_dy[j,i] = harmonic_average_one_fluid(mu1,mu2,volume_fraction_face)
            end
        end

        display(viscosity_coeff_for_du_dy)

        printstyled(color=:magenta, @sprintf "\n last row \n") 


        # # Example usage:
        # x = 1.2
        # y = 2.3
        # x1, y1 = 1.0, 2.0
        # x2, y2 = 3.0, 4.0
        # Q11, Q12, Q21, Q22 = 10.0, 20.0, 30.0, 40.0


        # y2 | Q21 Q22 
        #    |    x
        # y1 | Q11 Q21
        
        #      x1  x2

        # result = bilinear_interpolation(x, y, x1, y1, x2, y2, Q11, Q12, Q21, Q22)
        # println("The interpolated value at ($x, $y) is $result")

        # x = 1.0
        # y = 2.0 
        # result = bilinear_interpolation(x, y, x1, y1, x2, y2, Q11, Q12, Q21, Q22)
        # println("The interpolated value at ($x, $y) is $result")

        # x = 1.0
        # y = 4.0 
        # result = bilinear_interpolation(x, y, x1, y1, x2, y2, Q11, Q12, Q21, Q22)
        # println("The interpolated value at ($x, $y) is $result")

        # x = 3.0
        # y = 2.0 
        # result = bilinear_interpolation(x, y, x1, y1, x2, y2, Q11, Q12, Q21, Q22)
        # println("The interpolated value at ($x, $y) is $result")

        # x = 3.0
        # y = 4.0 
        # result = bilinear_interpolation(x, y, x1, y1, x2, y2, Q11, Q12, Q21, Q22)
        # println("The interpolated value at ($x, $y) is $result")

        print("\n grid.x ",grid.x[1,:])
        print("\n grid_u.x ",grid_u.x[1,:])

        volume_fraction_1D = fnones(grid,num)

        volume_fraction_1D .= 10.0

        # right border
        i = grid_u.nx
        for j in 2:grid_u.ny
            du_dy_coord_x = x_centroid_u[j,i] #grid_u.x[j,i]
            du_dy_coord_y = (y_centroid_u[j-1,i] + y_centroid_u[j,i])/2 #(grid_u.y[j,i] + grid_u.y[j+1,i])/2
            x1 = grid.x[j-1,i-1]
            x2 = grid_u.x[j-1,i]
            y1 = grid.y[j-1,i-1]
            y2 = grid.y[j,i-1]

            Q11 = volume_fraction[j-1,i-1]
            Q12 = volume_fraction[j,i-1]
            Q21 = vecb_R(volume_fraction_1D,grid)[j-1]
            Q22 = vecb_R(volume_fraction_1D,grid)[j]

            volume_fraction_face = bilinear_interpolation(du_dy_coord_x, du_dy_coord_y, x1, y1, x2, y2, Q11, Q12, Q21, Q22)

            printstyled(color=:green, @sprintf "\n i %.5i j %.5i x %.2e y %.2e x1 %.2e y1 %.2e x2 %.2e y2 %.2e Q11 %.2e Q12 %.2e Q21 %.2e Q22 %.2e\n" i j du_dy_coord_x du_dy_coord_y x1 y1 x2 y2 Q11 Q12 Q21 Q22)

            print("\nvolume_fraction i ",i," j ",j," ",volume_fraction_face)
            viscosity_coeff_for_du_dy[j,i] = harmonic_average_one_fluid(mu1,mu2,volume_fraction_face)
            
        end
        

        # left border
        i = 1
        for j in 2:grid_u.ny
            du_dy_coord_x = x_centroid_u[j,i] #grid_u.x[j,i]
            du_dy_coord_y = (y_centroid_u[j-1,i] + y_centroid_u[j,i])/2 #(grid_u.y[j,i] + grid_u.y[j+1,i])/2
            x1 = grid_u.x[j-1,i]
            x2 = grid.x[j-1,i]
            y1 = grid.y[j-1,i]
            y2 = grid.y[j,i]

            Q11 = vecb_L(volume_fraction_1D,grid)[j-1]
            Q12 = vecb_L(volume_fraction_1D,grid)[j]
            Q21 = volume_fraction[j-1,i]
            Q22 = volume_fraction[j,i]

            volume_fraction_face = bilinear_interpolation(du_dy_coord_x, du_dy_coord_y, x1, y1, x2, y2, Q11, Q12, Q21, Q22)

            printstyled(color=:green, @sprintf "\n i %.5i j %.5i x %.2e y %.2e x1 %.2e y1 %.2e x2 %.2e y2 %.2e Q11 %.2e Q12 %.2e Q21 %.2e Q22 %.2e\n" i j du_dy_coord_x du_dy_coord_y x1 y1 x2 y2 Q11 Q12 Q21 Q22)

            print("\nvolume_fraction i ",i," j ",j," ",volume_fraction_face)
            viscosity_coeff_for_du_dy[j,i] = harmonic_average_one_fluid(mu1,mu2,volume_fraction_face)
            
        end

    end #old method
    #endregion old method of interpolation



#     display(viscosity_coeff_for_du_dy)

#     print("\n size viscosity_coeff_for_du_dy ",size(viscosity_coeff_for_du_dy),"\n")


#     diag_viscosity_coeff_for_du_dy = Diagonal(vec(viscosity_coeff_for_du_dy))

#     mul!(opC_u.tmp_y, diag_viscosity_coeff_for_du_dy * opC_u.iMy, opC_u.By)
#     diffusion_bulk_u = L .+ opC_u.ByT * opC_u.tmp_y

#     j = div(grid.ny,2)
#     i = 1

#     II = CartesianIndex(j,i)

#     print("\n diffusion_bulk_u ", diffusion_bulk_u[II])


    printstyled(color=:magenta, @sprintf "\n testing diffusion_u_v \n") 
end
#endregion diffusion u v