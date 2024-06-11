# module Plotpython

# using PyCall
# using LaTeXStrings
# using PyPlot
# @pyimport matplotlib.animation as anim

# __precompile__() # this module is safe to precompile
# module MyModule
# using PyCall

# const scipy_opt = PyNULL()

# function __init__()
#     copy!(scipy_opt, pyimport_conda("scipy.optimize", "scipy"))
# end

# const anim = PyNULL()

# function __init__()
#     copy!(anim, pyimport_conda("matplotlib.animation", "matplotlib"))
# end

# @pyimport matplotlib.animation as anim


# end


# rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
# rcParams["text.latex.preamble"] = [raw"\usepackage{siunitx}"]

# plt.rc("text", usetex=true)
# # PyCall.PyDict(plt."rcParams")["text.latex.preamble"] = [raw"\usepackage{amsmath}"]
# PyCall.PyDict(plt."rcParams")["text.latex.preamble"] = [raw"\usepackage{siunitx}"]

# Pyplot.rc("text", usetex=true)
# PyCall.PyDict(Pyplot."rcParams")["text.latex.preamble"] = [raw"\usepackage{siunitx}"]


# function python

function strtitlefunc(isnap,fwd)
    # strtitle = @sprintf "t %.2e radius %.2e" fwd.t[i+1] fwd.radius[i+1]
    strtitle = @sprintf "t %.2e (ms) radius %.2e (mu m)" fwd.t[isnap]*1e3 fwd.radius[isnap]*1e6
    return strtitle
end

function plot_python_pdf(itmp,field0,figname,prefix,plot_levelset,isocontour,plot_grid,plot_mode,levels,range,cmap,x_array,y_array,gp,cbarlabel,i0,i1,j0,j1,fwd)

    if itmp<0
        # print("\n testi ", itmp)
        i=-itmp
        # field = field0[1:j1-j0+1,1:i1-i0+1]
        field = field0[j0:j1,i0:i1]
        # i0=1
        # i1=ii1-ii0+1
        # j0=1
        # j1=

    else
        i=itmp
        field = field0[i,j0:j1,i0:i1]
        # i0=ii0
        # i1=ii1
        # j0=jj0
        # j1=jj1
    end
    # PyPlot.rc("text", usetex=true)
    # rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    # rcParams["text.latex.preamble"] = raw"\usepackage{siunitx}"

    x_arr=x_array[i0:i1]
    y_arr=y_array[j0:j1]

    # field = field0[i,i0:i1,j0:j1]


    # print("range",range)
    fig1, ax2 = plt.subplots(layout="constrained")
    # CS = ax2.contourf(x_arr,y_arr,max.((phL.trans_scal[:,:,1] .-c0_H2)./c0_H2,0.0), 10, cmap=cmap)

   
    if plot_mode == "contourf"
        if levels==0
            CS = ax2.contourf(x_arr,y_arr,field, 
            levels=range, #10, 
            cmap=cmap)
        else
            CS = ax2.contourf(x_arr,y_arr,field, 
            levels=levels,
            cmap=cmap)           
        end
    elseif plot_mode == "pcolormesh"
        if levels==0            
            norm = mpl_colors.BoundaryNorm(range, ncolors=cmap.N, clip=true)
            CS = ax2.pcolormesh(x_arr,y_arr,field, cmap=cmap, norm=norm)
        else          
            levels = mpl_tickers.MaxNLocator(nbins=levels).tick_values(minimum(field), maximum(field))
            norm = mpl_colors.BoundaryNorm(levels, ncolors=cmap.N, clip=true)
            CS = ax2.pcolormesh(x_arr,y_arr,field, cmap=cmap, norm=norm)
        end
    end

    
    lcolor= "w" #"k"
    lw=0.5
    ms=0.5
    fontsize=5

    if plot_grid
        # for igrid in i0:i1
        #     ax2.axvline(x_array[igrid],c=lcolor,lw=lw)
        # end
        # for igrid in j0:j1
        #     ax2.axhline(y_array[igrid],c=lcolor,lw=lw)
        # end
        for igrid0 in i0:i1        
            for jgrid0 in j0:j1
                jgrid=jgrid0-j0+1
                igrid=igrid0-i0+1

                # print("\n",x_array[igrid],y_array[jgrid])
                # ax2.scatter(x_array[igrid],y_array[jgrid],
                # c=lcolor,
                # s=ms,
                # )

                # str=@sprintf "%.2e" field0[i,jgrid,igrid]
                # str=@sprintf "%.5e" field0[i,jgrid,igrid]
                str=@sprintf "%.2e" field[jgrid,igrid]


                ax2.annotate(str,(x_array[igrid],y_array[jgrid]),fontsize=fontsize,c=lcolor,ha="center")

            end
        end
    end


    # ax2.set_title("Title")
    # ax2.set_xlabel(L"$x (\mu m)$")
    # ax2.set_ylabel(L"$y (\mu m)$")

    ax2.set_xlabel(raw"$x ( \unit{\um})$")
    ax2.set_ylabel(L"$y (\mu m)$")

    # Make a colorbar for the ContourSet returned by the contourf call.
    cbar = fig1.colorbar(CS)
    cbar.ax.set_ylabel(cbarlabel)
    # Add the contour line levels to the colorbar
    if isocontour
        CS2 = ax2.contour(CS, 
        # levels=CS.levels[::2], 
        # levels=
        colors="r")
        cbar.add_lines(CS2)
    end

    indLS = max(i-1,1) 
   

    if plot_levelset
        # CSlvl = ax2.contourf(x_arr,y_arr,(phL.trans_scal[:,:,1] .-c0_H2)./c0_H2, levels=0.0, cmap=cmap)
        # CS2 = ax2.contour(CSlvl, 
        # # levels=CS.levels[::2], 
        # # levels=
        # colors="r")
        # cbar.add_lines(CS2)
        # CSlvl = ax2.contour(x_arr,y_arr, gp.LS[1].u, [0.0],colors="r")
        # CSlvl = ax2.contour(x_arr,y_arr, fwd.u[1,i,i0:i1,j0:j1], [0.0],colors="r")
     

        CSlvl = ax2.contour(x_arr,y_arr, fwd.u[1,indLS,j0:j1,i0:i1], [0.0],colors="r")

    end

    # strtitle=strtitlefunc(indLS,fwd)
    isnap = indLS
    strtime = @sprintf "%.2e" fwd.t[isnap]*1e3 
    strrad = @sprintf "%.2e" fwd.radius[isnap]*1e6

    plt.title("t "*strtime*raw"$( \unit{\ms})$"*"radius "*strrad*raw"$( \unit{\um})$")

    # if plot_levelset
    #     gp.LS[1].u .= sqrt.((gp.x .+ xcoord).^ 2 + (gp.y .+ ycoord) .^ 2) - (radius) * ones(gp);

    #     CSlvl = ax2.contour(x_arr,y_arr, gp.LS[1].u, [0.0],colors="r")
    # end

    plt.axis("equal")

    str_it = @sprintf "_%.5i" i
    plt.savefig(prefix*figname*str_it*".pdf")
    plt.close(fig1)

end


function plot_python_pdf_full(i,field0,figname,prefix,plot_levelset,isocontour,plot_grid,plot_mode,levels,range,cmap,x_array,y_array,gp,cbarlabel,ii0,ii1,jj0,jj1,fwd,fwdL,xscale)

    # PyPlot.rc("text", usetex=true)
    # rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    # rcParams["text.latex.preamble"] = raw"\usepackage{siunitx}"

    # if itmp<0
    #     i=-itmp
    #     field = field0[1:j1-j0+1,1:i1-i0+1]
    # else
    #     i=itmp
    #     field = field0[i,j0:j1,i0:i1]
    # end

    # field = field0[i,i0:i1,j0:j1]

    i0=ii0
    i1=ii1
    j0=jj0
    j1=jj1

    x_arr=x_array[i0:i1]
    y_arr=y_array[j0:j1]

    vecb_l=false


    if ii0 == 1
        vecb_l=true
        i1+=1

        # print("\n x_arr",x_arr)

        pushfirst!(x_arr,x_array[i0]-0.5*gp.dx[1,1]/xscale)

        # print("\n x_arr",x_arr)
        
        fieldtmp = zeros(gp.ny, gp.nx + 1)
        # print(size(field))
        # print(size(field0))
        # printstyled(color=:green, @sprintf "\n j: %5i %5i %5i %5i\n" i0 i1 j0 j1)

        @views fieldtmp[j0:j1,2:i1] = field0[i,jj0:jj1,ii0:ii1]
        @views fieldtmp[:,1] = vecb_L(fwdL.trans_scalD[i,:,1],gp)


        field = fieldtmp[j0:j1,i0:i1]

        # @views pushfirst!(field[j0:j1,:],vecb_L(fwd.trans_scalD[i,:,1],gp))
    else
        field = field0[i,j0:j1,i0:i1]

    end


    # print("range",range)
    fig1, ax2 = plt.subplots(layout="constrained")
    # CS = ax2.contourf(x_arr,y_arr,max.((phL.trans_scal[:,:,1] .-c0_H2)./c0_H2,0.0), 10, cmap=cmap)

   
    
    if levels==0
        CS = ax2.contourf(x_arr,y_arr,field, 
        levels=range, #10, 
        cmap=cmap)
    else
        # CS = ax2.contourf(x_arr,y_arr,field, 
        # levels=levels,
        # cmap=cmap)
        # # print("levels " ,levels)
        mpl_levels = mpl_tickers.MaxNLocator(nbins=levels).tick_values(minimum(field), maximum(field))
        norm = mpl_colors.BoundaryNorm(mpl_levels, ncolors=cmap.N, clip=true)
        CS = ax2.pcolormesh(x_arr,y_arr,field, cmap=cmap, norm=norm)
        # fig.colorbar(im, ax=ax0)
    end

    
    lcolor= "k" #"w" #"k"
    lw=0.5
    ms=0.5
    fontsize=5

    if plot_grid
        # for igrid in i0:i1
        #     ax2.axvline(x_array[igrid],c=lcolor,lw=lw)
        # end
        # for igrid in j0:j1
        #     ax2.axhline(y_array[igrid],c=lcolor,lw=lw)
        # end


        for igrid0 in i0:i1        
            for jgrid0 in j0:j1
                # print("\n",x_array[igrid],y_array[jgrid])
                # print("\n ",igrid," ",jgrid)
                igrid=igrid0
                jgrid = jgrid0-j0+1

                if igrid==1
                    lcolor= "k"
                else
                    lcolor= "w"
                end
                if igrid%2 == 0
                    va="top"
                else
                    va="bottom"
                end
                ax2.scatter(x_arr[igrid],y_arr[jgrid],
                c=lcolor,
                s=ms,
                )

                
                # str=@sprintf "%.2e" fieldtmp[jgrid0,igrid0]
                str=@sprintf "%.2e %.3i %.3i" fieldtmp[jgrid0,igrid0] igrid0 jgrid0


                ax2.annotate(str,(x_arr[igrid],y_arr[jgrid]),fontsize=fontsize,c=lcolor,ha="center",va=va)

            end
        end
    end


    # ax2.set_title("Title")
    # ax2.set_xlabel(L"$x (\mu m)$")
    # ax2.set_ylabel(L"$y (\mu m)$")

    ax2.set_xlabel(raw"$x ( \unit{\um})$")
    ax2.set_ylabel(L"$y (\mu m)$")

    # Make a colorbar for the ContourSet returned by the contourf call.
    cbar = fig1.colorbar(CS)
    cbar.ax.set_ylabel(cbarlabel)
    # Add the contour line levels to the colorbar
    if isocontour
        CS2 = ax2.contour(CS, 
        # levels=CS.levels[::2], 
        # levels=
        colors="r")
        cbar.add_lines(CS2)
    end

    indLS = max(i-1,1) 
   
    if plot_levelset
        # CSlvl = ax2.contourf(x_arr,y_arr,(phL.trans_scal[:,:,1] .-c0_H2)./c0_H2, levels=0.0, cmap=cmap)
        # CS2 = ax2.contour(CSlvl, 
        # # levels=CS.levels[::2], 
        # # levels=
        # colors="r")
        # cbar.add_lines(CS2)
        # CSlvl = ax2.contour(x_arr,y_arr, gp.LS[1].u, [0.0],colors="r")
        # CSlvl = ax2.contour(x_arr,y_arr, fwd.u[1,i,i0:i1,j0:j1], [0.0],colors="r")
     

        CSlvl = ax2.contour(x_arr,y_arr, fwd.u[1,indLS,j0:j1,i0:i1], [0.0],colors="r")

    end

    # strtitle=strtitlefunc(indLS,fwd)
    isnap = indLS
    strtime = @sprintf "%.2e" fwd.t[isnap]*1e3 
    strrad = @sprintf "%.2e" fwd.radius[isnap]*1e6

    plt.title("t "*strtime*raw"$( \unit{\ms})$"*"radius "*strrad*raw"$( \unit{\um})$")

    # if plot_levelset
    #     gp.LS[1].u .= sqrt.((gp.x .+ xcoord).^ 2 + (gp.y .+ ycoord) .^ 2) - (radius) * ones(gp);

    #     CSlvl = ax2.contour(x_arr,y_arr, gp.LS[1].u, [0.0],colors="r")
    # end

    plt.axis("equal")



    str_it = @sprintf "_%.5i" i
    plt.savefig(prefix*figname*str_it*".pdf")
    plt.close(fig1)

end



function python_movie_zoom(field0,figname,prefix,plot_levelset,isocontour,plot_mode,levels,range,cmap,x_array,y_array,gp,cbarlabel,size_frame,i0,i1,j0,j1,fwd)

    # if step!=0
    #     levels=range(lmin,lmax,step)
    # else
    #     levels=10
    # end
    nlevels=10

    fig1, ax2 = plt.subplots(layout="constrained")

    x_arr=x_array[i0:i1]
    y_arr=y_array[j0:j1]

    # field = field0[:,i0:i1,j0:j1]
    field = field0[:,j0:j1,i0:i1]



    # if levels==0
    #     CS = ax2.contourf(x_arr,y_arr,field[1,:,:], 
    #     levels=range,
    #     cmap=cmap)
    # else
    #     CS = ax2.contourf(x_arr,y_arr,field[1,:,:], 
    #     levels=levels,
    #     cmap=cmap)
    # end

    if plot_mode == "contourf"
        plot_mode_contourf = true
        plot_mode_pcolormesh = false
        if levels==0
            CS = ax2.contourf(x_arr,y_arr,field[1,:,:], 
            levels=range, #10, 
            cmap=cmap)
        else
            CS = ax2.contourf(x_arr,y_arr,field[1,:,:], 
            levels=levels,
            cmap=cmap)           
        end
    elseif plot_mode == "pcolormesh"
        plot_mode_contourf = false
        plot_mode_pcolormesh = true
        if levels==0            
            norm = mpl_colors.BoundaryNorm(range, ncolors=cmap.N, clip=true)
            CS = ax2.pcolormesh(x_arr,y_arr,field[1,:,:], cmap=cmap, norm=norm)
        else          
            mpl_levels = mpl_tickers.MaxNLocator(nbins=levels).tick_values(minimum(field[1,:,:]), maximum(field[1,:,:]))
            norm = mpl_colors.BoundaryNorm(mpl_levels, ncolors=cmap.N, clip=true)
            CS = ax2.pcolormesh(x_arr,y_arr,field[1,:,:], cmap=cmap, norm=norm)
        end
    end
  

    # Make a colorbar for the ContourSet returned by the contourf call.
    cbar = fig1.colorbar(CS)
    cbar.ax.set_ylabel(cbarlabel)

    function make_frame(i)
        # ax1.clear()
        ax2.clear()
        # ax1.imshow(A[:,:,i+1, 1])

        # CS = ax2.contourf(x_arr,y_arr,max.((phL.trans_scal[:,:,1] .-c0_H2)./c0_H2,0.0), 10, cmap=cmap)

        # if levels==0
        #     CS = ax2.contourf(x_arr,y_arr,field[i+1,:,:], 
        #     levels=range,
        #     cmap=cmap)
        # else
        #     CS = ax2.contourf(x_arr,y_arr,field[i+1,:,:], 
        #     levels=levels,
        #     cmap=cmap)
        # end

        if plot_mode_contourf
            if levels==0
                CS = ax2.contourf(x_arr,y_arr,field[i+1,:,:], 
                levels=range, #10, 
                cmap=cmap)
            else
                CS = ax2.contourf(x_arr,y_arr,field[i+1,:,:], 
                levels=levels,
                cmap=cmap)           
            end
        elseif plot_mode_pcolormesh
            if levels==0            
                norm = mpl_colors.BoundaryNorm(range, ncolors=cmap.N, clip=true)
                CS = ax2.pcolormesh(x_arr,y_arr,field[i+1,:,:], cmap=cmap, norm=norm)
            else          
                mpl_levels = mpl_tickers.MaxNLocator(nbins=levels).tick_values(minimum(field[i+1,:,:]), maximum(field[i+1,:,:]))
                norm = mpl_colors.BoundaryNorm(mpl_levels, ncolors=cmap.N, clip=true)
                CS = ax2.pcolormesh(x_arr,y_arr,field[i+1,:,:], cmap=cmap, norm=norm)
            end
        end

        # CS = ax2.contourf(x_arr,y_arr,field[i+1,:,:], 
        # # levels=10, 
        # levels=levels,
        # cmap=cmap)

        indLS= i 
        if i ==0 
            indLS = 1
        end


        strtitle=strtitlefunc(indLS,fwd)
        plt.title(strtitle)

        # CS = ax2.contourf(x_arr,y_arr,(fwd.trans_scal[i+1,:,:,1] .-c0_H2)./c0_H2, 
        # # levels=10, 
        # levels=range(0,1400,step=200),
        # cmap=cmap)

        plt.axis("equal")


        # ax2.imshow(A[:,:,i+1])

        # ax2.set_title("Title")
        ax2.set_xlabel(L"$x (\mu m)$")
        ax2.set_ylabel(L"$y (\mu m)$")

        # Make a colorbar for the ContourSet returned by the contourf call.
        #  cbar = fig1.colorbar(CS)
        #  cbar.ax.set_ylabel(cbarlabel)

        if step!=0
            fig1.colorbar(CS,cax=cbar.ax)
        end

        #https://stackoverflow.com/questions/5180518/duplicated-colorbars-when-creating-an-animation

        if (i==0) #(i+1==1) python starts at 0
            # # Make a colorbar for the ContourSet returned by the contourf call.
            # cbar = fig1.colorbar(CS)
            # cbar.ax.set_ylabel(cbarlabel)
            
            # Add the contour line levels to the colorbar
            if isocontour
                CS2 = ax2.contour(CS, 
                # levels=CS.levels[::2], 
                # levels=
                colors="r")
                cbar.add_lines(CS2)
            end
        end
    
        if plot_levelset
            # CSlvl = ax2.contour(x_arr,y_arr, fwd.u[1,i+1,i0:i1,j0:j1], [0.0],colors="r")
            CSlvl = ax2.contour(x_arr,y_arr, fwd.u[1,indLS,j0:j1,i0:i1], [0.0],colors="r")
            # print("check levelset")
        end


        # plt.savefig(prefix*"H2.pdf")
        # plt.close(fig1)
    end

    myanim = anim.FuncAnimation(fig1, make_frame, frames=size_frame, interval=size_frame, blit=false)

    # myanim.save(prefix*"test.gif")
    myanim.save(prefix*figname*".mp4")

    plt.close("all")

end
###############################################################################################


"""
From Stefan_velocity!
"""
function plot_electrolysis_velocity!(num, grid, LS, V, TL, MIXED, periodic_x, periodic_y, concentration_scal_intfc)
    @unpack geoS, geoL, κ = LS

    # V .= 0
    # @inbounds @threads for II in MIXED

    
    x_array=grid.x[1,:]/num.plot_xscale
    y_array=grid.y[:,1]/num.plot_xscale


    fig1, ax2 = plt.subplots(layout="constrained")
    cmap = plt.cm.viridis

    CS = ax2.contourf(x_array,y_array, reshape(veci(TL,grid,1), grid), 10, cmap=cmap)


    # ax2.set_title("Title")
    ax2.set_xlabel(L"$x (\mu m)$")
    ax2.set_ylabel(L"$y (\mu m)$")

    # Make a colorbar for the ContourSet returned by the contourf call.
    cbar = fig1.colorbar(CS)
    cbar.ax.set_ylabel("concentration")

    # cbar.ax.set_title(L"$10^{-4}$")
    plt.axis("equal")    

    @inbounds for II in MIXED

        θ_d = concentration_scal_intfc

        # dTS = 0.
        dTL = 0.
        if geoL.projection[II].flag
            T_1, T_2 = interpolated_temperature(grid, geoL.projection[II].angle, geoL.projection[II].point1, geoL.projection[II].point2, TL, II, periodic_x, periodic_y)
            dTL = normal_gradient(geoL.projection[II].d1, geoL.projection[II].d2, T_1, T_2, θ_d)
        else
            T_1 = interpolated_temperature(grid, geoL.projection[II].angle, geoL.projection[II].point1, TL, II, periodic_x, periodic_y)
            dTL = normal_gradient(geoL.projection[II].d1, T_1, θ_d)
        end
        # V[II] = dTL #+ dTS

        print("\n grad",II,"val",dTL,"val",geoL.projection[II].flag,"val",geoL.projection[II].angle,"val", geoL.projection[II].point1,"val", geoL.projection[II].point2,"val", T_1,"val",T_2,"val",θ_d)
        
        xp1=geoL.projection[II].point1.x/num.plot_xscale
        yp1=geoL.projection[II].point1.y/num.plot_xscale

        xp2=geoL.projection[II].point2.x/num.plot_xscale
        yp2=geoL.projection[II].point2.y/num.plot_xscale

        print("\n test",xp1,geoL.projection[II].point1.x)


        print("\n plotscale",num.plot_xscale)
        # print("\n test",xp1)
        # print("\n test",yp1)
        # print("\n test",xp2)
        # print("\n test",yp2)

        #TODO cannot use ' but "

        plt.scatter(xp1,yp1,c="b")
        plt.scatter(xp2,yp2,c="r")
   
    end

    plt.savefig(num.plot_prefix*"jump.pdf")
    plt.close(fig1)


    return nothing
end


function plot_python_bc(num,x_array,y_array,field,figname,prefix,grid)
    # plot_levelset,isocontour,levels,range,cmap,,gp,cbarlabel)

    # PyPlot.rc("text", usetex=true)
    # rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    # rcParams["text.latex.preamble"] = raw"\usepackage{siunitx}"

    # print("range",range)
    fig1, ax2 = plt.subplots(layout="constrained")
    # CS = ax2.contourf(x_array,y_array,max.((phL.trans_scal[:,:,1] .-c0_H2)./c0_H2,0.0), 10, cmap=cmap)

    ms=0.5
    mcolor="r"

    # plt.plot(y_array,field)
    # plt.scatter(y_array,field,s=ms,c=mcolor)

    lw=0.1
    lcolor="k"

    plt.axhline((-num.shifted_y+num.R)/num.plot_xscale,lw=lw,c=lcolor)

    plt.axhline((-num.shifted_y-num.R)/num.plot_xscale,lw=lw,c=lcolor)

    plt.plot(field,y_array,lw=lw,c=mcolor)
    plt.scatter(field,y_array,s=ms,c=mcolor)




    # plt.plot(y_array,grid.LS[1].u[:,1])


    
    # if levels==0
    #     CS = ax2.contourf(x_array,y_array,field, 
    #     levels=range, #10, 
    #     cmap=cmap)
    # else
    #     CS = ax2.contourf(x_array,y_array,field, 
    #     levels=levels,
    #     cmap=cmap)
    #     # print("levels " ,levels)
    # end

    # # ax2.set_title("Title")
    # # ax2.set_xlabel(L"$x (\mu m)$")
    # # ax2.set_ylabel(L"$y (\mu m)$")

    # ax2.set_xlabel(raw"$x ( \unit{\um})$")
    # ax2.set_ylabel(L"$y (\mu m)$")

    # # Make a colorbar for the ContourSet returned by the contourf call.
    # cbar = fig1.colorbar(CS)
    # cbar.ax.set_ylabel(cbarlabel)
    # # Add the contour line levels to the colorbar
    # if isocontour
    #     CS2 = ax2.contour(CS, 
    #     # levels=CS.levels[::2], 
    #     # levels=
    #     colors="r")
    #     cbar.add_lines(CS2)
    # end
    # if plot_levelset
    #     # CSlvl = ax2.contourf(x_array,y_array,(phL.trans_scal[:,:,1] .-c0_H2)./c0_H2, levels=0.0, cmap=cmap)
    #     # CS2 = ax2.contour(CSlvl, 
    #     # # levels=CS.levels[::2], 
    #     # # levels=
    #     # colors="r")
    #     # cbar.add_lines(CS2)
    #     CSlvl = ax2.contour(x_array,y_array, gp.LS[1].u, [0.0],colors="r")
    # end


    # plt.axis("equal")

    plt.savefig(prefix*figname*".pdf")
    plt.close(fig1)

end

function plot_last_iter_python_pdf(field,figname,prefix,plot_levelset,isocontour,levels,range,cmap,x_array,y_array,gp,cbarlabel)

    # PyPlot.rc("text", usetex=true)
    # rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    # rcParams["text.latex.preamble"] = raw"\usepackage{siunitx}"

    # print("range",range)
    fig1, ax2 = plt.subplots(layout="constrained")
    # CS = ax2.contourf(x_array,y_array,max.((phL.trans_scal[:,:,1] .-c0_H2)./c0_H2,0.0), 10, cmap=cmap)
    
    if levels==0
        CS = ax2.contourf(x_array,y_array,field, 
        levels=range, #10, 
        cmap=cmap)
    else
        CS = ax2.contourf(x_array,y_array,field, 
        levels=levels,
        cmap=cmap)
        # print("levels " ,levels)
    end

    # ax2.set_title("Title")
    # ax2.set_xlabel(L"$x (\mu m)$")
    # ax2.set_ylabel(L"$y (\mu m)$")

    ax2.set_xlabel(raw"$x ( \unit{\um})$")
    ax2.set_ylabel(L"$y (\mu m)$")

    # Make a colorbar for the ContourSet returned by the contourf call.
    cbar = fig1.colorbar(CS)
    cbar.ax.set_ylabel(cbarlabel)
    # Add the contour line levels to the colorbar
    if isocontour
        CS2 = ax2.contour(CS, 
        # levels=CS.levels[::2], 
        # levels=
        colors="r")
        cbar.add_lines(CS2)
    end
    if plot_levelset
        # CSlvl = ax2.contourf(x_array,y_array,(phL.trans_scal[:,:,1] .-c0_H2)./c0_H2, levels=0.0, cmap=cmap)
        # CS2 = ax2.contour(CSlvl, 
        # # levels=CS.levels[::2], 
        # # levels=
        # colors="r")
        # cbar.add_lines(CS2)
        CSlvl = ax2.contour(x_array,y_array, gp.LS[1].u, [0.0],colors="r")
    end

    # if plot_levelset
    #     gp.LS[1].u .= sqrt.((gp.x .+ xcoord).^ 2 + (gp.y .+ ycoord) .^ 2) - (radius) * ones(gp);

    #     CSlvl = ax2.contour(x_array,y_array, gp.LS[1].u, [0.0],colors="r")
    # end

    plt.axis("equal")

    plt.savefig(prefix*figname*".pdf")
    plt.close(fig1)

end




# function plot_python_pdf(field,name,plot_levelset,range)
#     # print("range",range)
#     fig1, ax2 = plt.subplots(layout="constrained")
#     # CS = ax2.contourf(x_array,y_array,max.((phL.trans_scal[:,:,1] .-c0_H2)./c0_H2,0.0), 10, cmap=cmap)
#     CS = ax2.contourf(x_array,y_array,field, 
#     levels=range, #10, 
#     cmap=cmap)

#     # ax2.set_title("Title")
#     ax2.set_xlabel(L"$x (\mu m)$")
#     ax2.set_ylabel(L"$y (\mu m)$")

#     # Make a colorbar for the ContourSet returned by the contourf call.
#     cbar = fig1.colorbar(CS)
#     cbar.ax.set_ylabel(cbarlabel)
#     # Add the contour line levels to the colorbar
#     if isocontour
#         CS2 = ax2.contour(CS, 
#         # levels=CS.levels[::2], 
#         # levels=
#         colors="r")
#         cbar.add_lines(CS2)
#     end
#     if plot_levelset
#         # CSlvl = ax2.contourf(x_array,y_array,(phL.trans_scal[:,:,1] .-c0_H2)./c0_H2, levels=0.0, cmap=cmap)
#         # CS2 = ax2.contour(CSlvl, 
#         # # levels=CS.levels[::2], 
#         # # levels=
#         # colors="r")
#         # cbar.add_lines(CS2)
#         CSlvl = ax2.contour(x_array,y_array, gp.LS[1].u, [0.0],colors="r")
#     end

#     # if plot_levelset
#     #     gp.LS[1].u .= sqrt.((gp.x .+ xcoord).^ 2 + (gp.y .+ ycoord) .^ 2) - (radius) * ones(gp);

#     #     CSlvl = ax2.contour(x_array,y_array, gp.LS[1].u, [0.0],colors="r")
#     # end

#     plt.axis("equal")

#     plt.savefig(prefix*name*".pdf")
#     plt.close(fig1)

# end

######################################################################################################


# function plot_python_several_pdf(field,name,plot_levelset,size_frame)

#     for isnap in 1:size_frame

#         fig1, ax2 = plt.subplots(layout="constrained")
#         # CS = ax2.contourf(x_array,y_array,max.((phL.trans_scal[:,:,1] .-c0_H2)./c0_H2,0.0), 10, cmap=cmap)
#         # CS = ax2.contourf(x_array,y_array,(fwd.trans_scal[isnap,:,:,1] .-c0_H2)./c0_H2, 10, cmap=cmap)
#         CS = ax2.contourf(x_array,y_array, field[isnap,:,:], 
#         10, 
#         cmap=cmap)
#         # print("nx ny ", nx,ny)
#         nplot=5
#         nplotx=nx ÷ nplot 
#         nploty=ny ÷ nplot
#         fontsize=5
#         startmod=1
#         ms=2
#         mcolor="w"

#         # strtitle = @sprintf "t %.2e radius %.2e" fwd.t[isnap] fwd.radius[isnap]
#         # strtitle = @sprintf "t %.2e (ms) radius %.2e (mm)" fwd.t[isnap]*1e3 fwd.radius[isnap]*1e6
#         strtitle=strtitlefunc(isnap)

#         plt.title(strtitle)

#         for i in 1:gp.nx
#             for j in 1:gp.ny
#                 if (i%nplotx==startmod) && (j%nploty==startmod)

#                     # print("\nplot i ",i," j ",j)

#                     str=@sprintf "%.2e" fwd.trans_scal[isnap,j,i,1]
#                     # str=@sprintf "%.2e %.2e" x_array[i] y_array[j]
#                     # str=@sprintf "%.2e" x_array[i]

#                     # str=@sprintf "%.5i %.5i" i j 


#                     ax2.annotate(str,(x_array[i],y_array[j]),fontsize=fontsize)
#                     plt.scatter(x_array[i],y_array[j],s=ms,c=mcolor)

#                 end
#             end
#         end
#         # ax2.set_title("Title")
#         ax2.set_xlabel(L"$x (\mu m)$")
#         ax2.set_ylabel(L"$y (\mu m)$")

#         # Make a colorbar for the ContourSet returned by the contourf call.
#         cbar = fig1.colorbar(CS)
#         cbar.ax.set_ylabel(cbarlabel)
#         # Add the contour line levels to the colorbar
#         if isocontour
#             CS2 = ax2.contour(CS, 
#             # levels=CS.levels[::2], 
#             # levels=
#             colors="r")
#             cbar.add_lines(CS2)
#         end
#         if plot_levelset
#             # CSlvl = ax2.contourf(x_array,y_array,(phL.trans_scal[:,:,1] .-c0_H2)./c0_H2, levels=0.0, cmap=cmap)
#             # CS2 = ax2.contour(CSlvl, 
#             # # levels=CS.levels[::2], 
#             # # levels=
#             # colors="r")
#             # cbar.add_lines(CS2)
#             CSlvl = ax2.contour(x_array,y_array, fwd.u[1,isnap,:,:], [0.0],colors="r")
#         end

#         # if plot_levelset
#         #     gp.LS[1].u .= sqrt.((gp.x .+ xcoord).^ 2 + (gp.y .+ ycoord) .^ 2) - (radius) * ones(gp);

#         #     CSlvl = ax2.contour(x_array,y_array, gp.LS[1].u, [0.0],colors="r")
#         # end

#         plt.axis("equal")

#         plt.savefig(prefix*name*"_"*string(isnap)*".pdf")
#         plt.close(fig1)
#     end

# end

# function python_movie(field,name,plot_levelset,size_frame)

#     fig1, ax2 = plt.subplots(layout="constrained")

#     CS = ax2.contourf(x_array,y_array,field[1,:,:], 
#     # levels=10, 
#     levels=range(0,1400,step=200),
#     cmap=cmap)

#     # Make a colorbar for the ContourSet returned by the contourf call.
#     cbar = fig1.colorbar(CS)
#     cbar.ax.set_ylabel(cbarlabel)

#     function make_frame(i)
#         # ax1.clear()
#         ax2.clear()
#         # ax1.imshow(A[:,:,i+1, 1])

#         # CS = ax2.contourf(x_array,y_array,max.((phL.trans_scal[:,:,1] .-c0_H2)./c0_H2,0.0), 10, cmap=cmap)

#         CS = ax2.contourf(x_array,y_array,field[i+1,:,:], 
#         # levels=10, 
#         levels=range(0,1400,step=200),
#         cmap=cmap)


#         strtitle=strtitlefunc(i+1)
#         plt.title(strtitle)

#         # CS = ax2.contourf(x_array,y_array,(fwd.trans_scal[i+1,:,:,1] .-c0_H2)./c0_H2, 
#         # # levels=10, 
#         # levels=range(0,1400,step=200),
#         # cmap=cmap)

#         plt.axis("equal")


#         # ax2.imshow(A[:,:,i+1])

#         # ax2.set_title("Title")
#         ax2.set_xlabel(L"$x (\mu m)$")
#         ax2.set_ylabel(L"$y (\mu m)$")

#         # Make a colorbar for the ContourSet returned by the contourf call.
#         #  cbar = fig1.colorbar(CS)
#         #  cbar.ax.set_ylabel(cbarlabel)

#         #https://stackoverflow.com/questions/5180518/duplicated-colorbars-when-creating-an-animation

#         if (i==0) #(i+1==1) python starts at 0
#             # # Make a colorbar for the ContourSet returned by the contourf call.
#             # cbar = fig1.colorbar(CS)
#             # cbar.ax.set_ylabel(cbarlabel)
            
#             # Add the contour line levels to the colorbar
#             if isocontour
#                 CS2 = ax2.contour(CS, 
#                 # levels=CS.levels[::2], 
#                 # levels=
#                 colors="r")
#                 cbar.add_lines(CS2)
#             end
#         end
    
#         if plot_levelset
#             CSlvl = ax2.contour(x_array,y_array, fwd.u[1,i+1,:,:], [0.0],colors="r")
#         end


#         # plt.savefig(prefix*"H2.pdf")
#         # plt.close(fig1)
#     end

#     myanim = anim.FuncAnimation(fig1, make_frame, frames=size_frame, interval=size_frame, blit=false)

#     # myanim.save(prefix*"test.gif")
#     myanim.save(prefix*name*".mp4")

#     plt.close("all")

# end

# function python_movie_zoom(field,name,plot_levelset,size_frame,i0,i1,j0,j1,lmin::Integer,lmax::Integer,step::Integer)
#     # if step!=0
#     #     levels=range(lmin,lmax,step)
#     # else
#     #     levels=10
#     # end
#     nlevels=10

#     fig1, ax2 = plt.subplots(layout="constrained")

#     x_arr=x_array[i0:i1]
#     y_arr=y_array[j0:j1]

#     if step!=0
#         CS = ax2.contourf(x_array,y_array,field[1,:,:], 
#         levels=range(trunc(lmin),trunc(lmax),step),
#         cmap=cmap)
#     else
#         CS = ax2.contourf(x_array,y_array,field[1,:,:], 
#         levels=nlevels,
#         cmap=cmap)
#     end
#     # CS = ax2.contourf(x_array,y_array,field[1,:,:], 
#     # # levels=10, 
#     # levels=levels,
#     # cmap=cmap)

#     # Make a colorbar for the ContourSet returned by the contourf call.
#     cbar = fig1.colorbar(CS)
#     cbar.ax.set_ylabel(cbarlabel)

#     function make_frame(i)
#         # ax1.clear()
#         ax2.clear()
#         # ax1.imshow(A[:,:,i+1, 1])

#         # CS = ax2.contourf(x_array,y_array,max.((phL.trans_scal[:,:,1] .-c0_H2)./c0_H2,0.0), 10, cmap=cmap)

#         if step!=0
#             CS = ax2.contourf(x_array,y_array,field[i+1,:,:], 
#             levels=range(trunc(lmin),trunc(lmax),step),
#             cmap=cmap)
#         else
#             CS = ax2.contourf(x_array,y_array,field[i+1,:,:], 
#             levels=nlevels,
#             cmap=cmap)
#         end

#         # CS = ax2.contourf(x_array,y_array,field[i+1,:,:], 
#         # # levels=10, 
#         # levels=levels,
#         # cmap=cmap)


#         strtitle=strtitlefunc(i+1)
#         plt.title(strtitle)

#         # CS = ax2.contourf(x_array,y_array,(fwd.trans_scal[i+1,:,:,1] .-c0_H2)./c0_H2, 
#         # # levels=10, 
#         # levels=range(0,1400,step=200),
#         # cmap=cmap)

#         plt.axis("equal")


#         # ax2.imshow(A[:,:,i+1])

#         # ax2.set_title("Title")
#         ax2.set_xlabel(L"$x (\mu m)$")
#         ax2.set_ylabel(L"$y (\mu m)$")

#         # Make a colorbar for the ContourSet returned by the contourf call.
#         #  cbar = fig1.colorbar(CS)
#         #  cbar.ax.set_ylabel(cbarlabel)

#         #https://stackoverflow.com/questions/5180518/duplicated-colorbars-when-creating-an-animation

#         if (i==0) #(i+1==1) python starts at 0
#             # # Make a colorbar for the ContourSet returned by the contourf call.
#             # cbar = fig1.colorbar(CS)
#             # cbar.ax.set_ylabel(cbarlabel)
            
#             # Add the contour line levels to the colorbar
#             if isocontour
#                 CS2 = ax2.contour(CS, 
#                 # levels=CS.levels[::2], 
#                 # levels=
#                 colors="r")
#                 cbar.add_lines(CS2)
#             end
#         end
    
#         if plot_levelset
#             CSlvl = ax2.contour(x_array,y_array, fwd.u[1,i+1,:,:], [0.0],colors="r")
#         end


#         # plt.savefig(prefix*"H2.pdf")
#         # plt.close(fig1)
#     end

#     myanim = anim.FuncAnimation(fig1, make_frame, frames=size_frame, interval=size_frame, blit=false)

#     # myanim.save(prefix*"test.gif")
#     myanim.save(prefix*name*".mp4")

#     plt.close("all")

# end
# ###############################################################################################

# end