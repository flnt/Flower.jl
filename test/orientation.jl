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
)

gp, gu, gv = init_meshes(num)
op, phS, phL = init_fields(num, gp, gu, gv)

gp.LS[1].u .= 1.0 #deactivate interface


@time run_forward!(num, gp, gu, gv, op, phS, phL; navier_stokes = true,)
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


    #interpolate coefficient
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

        @test L[1,1]=-3.0 atol=test_tolerance

    end

    print("\n L ", L[1,1+grid.ny])


    print("\n L ", L[10,:])


end

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


vecb_L(phL.TD, gp) .= ftest.(gp.x[:,1] .- gp.x[:,1] ./ 2.0, gp.y[:,1])
vecb_B(phL.TD, gp) .= ftest.(gp.x[1,:], gp.y[1,:] .- gp.y[1,:] ./ 2.0)
vecb_R(phL.TD, gp) .= ftest.(gp.x[:,end] .+ gp.x[:,1] ./ 2.0, gp.y[:,1])
vecb_T(phL.TD, gp) .= ftest.(gp.x[1,:], gp.y[end,:] .+ gp.y[1,:] ./ 2.0)


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

grad_analytical = ftest_1.(
    x_centroid_u,
    y_centroid_u
)

j = div(n,2)

print("\n grad_analytical ",grad_analytical[j,:],"\n")
print("\n D ",grad_x[j,2:end-1],"\n")

@testset "Simple gradient test x" begin
    print("\n grad_x[j,2:end-1]",grad_x[j,2:end-1]," len ",length(grad_x[j,2:end-1]))
    test_error = maximum(abs.(grad_x[j,2:end-1] .- 1.0))
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