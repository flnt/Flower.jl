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
    strtitle = @sprintf "t %.2e (ms) radius %.2e (mm)" fwd.t[isnap]*1e3 fwd.radius[isnap]*1e6
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



function python_movie_zoom(field,figname,prefix,plot_levelset,isocontour,levels,range,cmap,x_array,y_array,gp,cbarlabel,size_frame,i0,i1,j0,j1,fwd)

    # if step!=0
    #     levels=range(lmin,lmax,step)
    # else
    #     levels=10
    # end
    nlevels=10

    fig1, ax2 = plt.subplots(layout="constrained")

    x_arr=x_array[i0:i1]
    y_arr=y_array[j0:j1]

    if levels==0
        CS = ax2.contourf(x_array,y_array,field[1,:,:], 
        levels=range,
        cmap=cmap)
    else
        CS = ax2.contourf(x_array,y_array,field[1,:,:], 
        levels=levels,
        cmap=cmap)
    end
    # CS = ax2.contourf(x_array,y_array,field[1,:,:], 
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

        # CS = ax2.contourf(x_array,y_array,max.((phL.trans_scal[:,:,1] .-c0_H2)./c0_H2,0.0), 10, cmap=cmap)

        if step!=0
            CS = ax2.contourf(x_array,y_array,field[i+1,:,:], 
            levels=range,
            cmap=cmap)
        else
            CS = ax2.contourf(x_array,y_array,field[i+1,:,:], 
            levels=levels,
            cmap=cmap)
        end

        # CS = ax2.contourf(x_array,y_array,field[i+1,:,:], 
        # # levels=10, 
        # levels=levels,
        # cmap=cmap)


        strtitle=strtitlefunc(i+1,fwd)
        plt.title(strtitle)

        # CS = ax2.contourf(x_array,y_array,(fwd.trans_scal[i+1,:,:,1] .-c0_H2)./c0_H2, 
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
            CSlvl = ax2.contour(x_array,y_array, fwd.u[1,i+1,:,:], [0.0],colors="r")
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