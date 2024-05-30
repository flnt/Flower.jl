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

function plot_python_pdf(field,figname,prefix,plot_levelset,isocontour,levels,range,cmap,x_array,y_array,gp,cbarlabel)

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



function python_movie_zoom(field0,figname,prefix,plot_levelset,isocontour,levels,range,cmap,x_array,y_array,gp,cbarlabel,size_frame,i0,i1,j0,j1,fwd)

    # if step!=0
    #     levels=range(lmin,lmax,step)
    # else
    #     levels=10
    # end
    nlevels=10

    fig1, ax2 = plt.subplots(layout="constrained")

    x_arr=x_array[i0:i1]
    y_arr=y_array[j0:j1]

    field = field0[:,i0:i1,j0:j1]

    if levels==0
        CS = ax2.contourf(x_arr,y_arr,field[1,:,:], 
        levels=range,
        cmap=cmap)
    else
        CS = ax2.contourf(x_arr,y_arr,field[1,:,:], 
        levels=levels,
        cmap=cmap)
    end
    # CS = ax2.contourf(x_arr,y_arr,field[1,:,:], 
    # # levels=10, 
    # levels=levels,
    # cmap=cmap)

    # Make a colorbar for the ContourSet returned by the contourf call.
    cbar = fig1.colorbar(CS)
    cbar.ax.set_ylabel(cbarlabel)

    function make_frame(i)
        # ax1.clear()
        ax2.clear()
        # ax1.imshow(A[:,:,i+1, 1])

        # CS = ax2.contourf(x_arr,y_arr,max.((phL.trans_scal[:,:,1] .-c0_H2)./c0_H2,0.0), 10, cmap=cmap)

        if levels==0
            CS = ax2.contourf(x_arr,y_arr,field[i+1,:,:], 
            levels=range,
            cmap=cmap)
        else
            CS = ax2.contourf(x_arr,y_arr,field[i+1,:,:], 
            levels=levels,
            cmap=cmap)
        end

        # CS = ax2.contourf(x_arr,y_arr,field[i+1,:,:], 
        # # levels=10, 
        # levels=levels,
        # cmap=cmap)


        strtitle=strtitlefunc(i+1,fwd)
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
            CSlvl = ax2.contour(x_arr,y_arr, fwd.u[1,i+1,i0:i1,j0:j1], [0.0],colors="r")
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