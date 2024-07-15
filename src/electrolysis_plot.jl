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


function plot_current_wall()

    fig, ax = plt.subplots(layout="constrained")

    # fig.subplots_adjust(right=0.75)

    # varx = [0, 1, 2]
    varx = gp.y[:,1]/xscale
    # vecb_L(phL.trans_scalD[:,3], gp)
    label1 = raw"$c\left(H_2O\right)$"
    label2 = raw"$-\eta ~\text{(-overpotential)}$ "
    label3 = "Current"
    alpha = 0.5
    # ls= (0, (5, 10)) #":" #"--" #"--"
    # # ls2="loosely dashed" #"--"
    # # ls2 = [0, [5, 10]]
    # ls2 = (5, (10, 3)) #"dashdot"
    # ls3 = (5, (10, 3)) #"dotted"

    # ls  = (0, (5, 10)) 
    # ls2 = (5, (5, 10)) #"dashdot"
    # ls3 = (10, (5, 10)) #"dotted"

    # ls  = (0, (5, 10)) 
    # ls2 = (10, (5, 10)) #"dashdot"
    # ls3 = (10, (5, 10)) #"dotted"

    ls  = (0, (3, 6)) 
    ls2 = (3, (3, 6)) #"dashdot"
    ls3 = (6, (3, 6)) #"dotted"

    # print("current", phL.i_current_mag[:,1])

    # print(colors)

    twin1 = ax.twinx()
    twin2 = ax.twinx()

    # Offset the right spine of twin2.  The ticks and label have already been
    # placed on the right by twinx above.
    twin2.spines.right.set_position(("axes", 1.2))
    #colors "C0", "C1", "C2"
    p1, = ax.plot(varx, phL.trans_scal[:,1,3],colors[1], label=label1,ls=ls)
    p2, = twin1.plot(varx, phL.phi_ele[:,1] .- phi_ele1, colors[2], label=label2,ls=ls2)
    p3, = twin2.plot(varx, phL.i_current_mag[:,1], colors[3], label=label3,ls=ls3)

    # p1, = ax.plot(varx, ones(gp.ny),colors[1], label=label1,ls=ls)
    # p2, = twin1.plot(varx, ones(gp.ny), colors[2], label=label2,ls=ls2)
    # p3, = twin2.plot(varx, ones(gp.ny), colors[3], label=label3,ls=ls3)

    ax.set(
        # xlim=(0, 2),
        # ylim=(0, 2),
        xlabel=raw"$y ( \unit{\um})$",
        ylabel=label1)
    twin1.set(
        # ylim=(0, 4), 
    ylabel=label2)
    twin2.set(
        # ylim=(1, 65), 
    ylabel=label3)

    ax.yaxis.label.set_color(p1.get_color())
    twin1.yaxis.label.set_color(p2.get_color())
    twin2.yaxis.label.set_color(p3.get_color())

    ax.tick_params(axis="y", colors=p1.get_color())
    twin1.tick_params(axis="y", colors=p2.get_color())
    twin2.tick_params(axis="y", colors=p3.get_color())

    # ax.legend(handles=[p1, p2, p3],
    # # loc = "center left",
    # loc = "outside upper left",
    # )

    fig.legend(handles=[p1, p2, p3],
    # loc = "center left",
    loc = "outside upper left",
    )

    plt.savefig(prefix*"electrode_c_eta_i.pdf")
    plt.close(fig)

end

function plot_bc(iter_list,vec,grid,xscale,figname,prefix,fwd)

    fig, ax = plt.subplots(layout="constrained")

    varx = grid.x[1,:]/xscale
    # print("\n varx ",varx)

    print(iter_list)
    for isnap in iter_list
        strtime = @sprintf "%.2e" fwd.t[isnap]*1e3 
        str = "t "*strtime*raw"$( \unit{\ms})$"
        # str=@sprintf "%.5i" i
        # print("\n vec ",vecb_T(vec[i,:],grid))
        plt.plot(varx,vecb_T(vec[isnap,:],grid),label = str)
    end
    eps=1e-4
    ylim0=0.16-eps
    ylim1=0.16+eps

    ax.set_ylim(ylim0,ylim1)
    plt.legend()
    plt.savefig(prefix*figname*".pdf")

    plt.close(fig)

end

function plot_bc2(iter_list,vec,grid,xscale,figname,prefix,fwd)

    fig, ax = plt.subplots(layout="constrained")

    varx = grid.x[1,:]/xscale
    # print("\n varx ",varx)

    print(iter_list)
    for isnap in iter_list
        strtime = @sprintf "%.2e" fwd.t[isnap]*1e3 
        str = "t "*strtime*raw"$( \unit{\ms})$"
        # str=@sprintf "%.5i" i
        # print("\n vec ",vecb_T(vec[i,:],grid))
        plt.plot(varx,vecb_B(vec[isnap,:],grid),label = str)
    end
    eps=1e-4
    ylim0=0.16-eps
    ylim1=0.16+eps

    ax.set_ylim(ylim0,ylim1)
    plt.legend()
    plt.savefig(prefix*figname*".pdf")

    plt.close(fig)

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

    ax2.spines["right"].set_visible(false)
    ax2.spines["top"].set_visible(false)
   
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
    # ax2.set_xlabel(raw"$x ( \unit{\um})$")

    # ax2.set_ylabel(raw"$y ( \unit{\um})$")


    ax2.set_xlabel(raw"$x ( \unit{\um})$")
    ax2.set_ylabel(raw"$y ( \unit{\um})$")


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
     
        if typeof(grid) == GridCC #isCC(grid)
            CSlvl = ax2.contour(x_arr,y_arr, fwd.u[1,indLS,j0:j1,i0:i1], [0.0],colors="r")
        elseif typeof(grid) == GridFCx #isFCx(grid)
            CSlvl = ax2.contour(x_arr,y_arr, fwd.ux[1,indLS,j0:j1,i0:i1], [0.0],colors="r")
        elseif typeof(grid) == GridFCy #isFCy(grid)
            CSlvl = ax2.contour(x_arr,y_arr, fwd.uy[1,indLS,j0:j1,i0:i1], [0.0],colors="r")
        end

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

function plot_python_pdf_full2(itmp,field0,field0D,figname,prefix,plot_levelset,isocontour,plot_grid,plot_mode,levels,range,cmap,x_array,y_array,gp,cbarlabel,ii0,ii1,jj0,jj1,fwd,fwdL,xscale,fontsize,printmode,plotcase,num,plotbc)

    # PyPlot.rc("text", usetex=true)
    # rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    # rcParams["text.latex.preamble"] = raw"\usepackage{siunitx}"

    i0=ii0
    i1=ii1
    j0=jj0
    j1=jj1
    i0tmp = i0 + 1
    j0tmp = j0 + 1
    i1tmp = i1 + 1
    j1tmp = j1 + 1

    i0tmp2 = i0 + 1
    j0tmp2 = j0 + 1
    i1tmp2 = i1 + 1
    j1tmp2 = j1 + 1

    if itmp<0
        i=-itmp
        field = field0
        field1 = field0
     else
        i=itmp
        field1 = field0[i,:,:]
        # print("\n i0 ",i0," i1 ",i1," j0 ",j0," j1 ",j1)
    end




    x_arr=x_array[i0:i1]
    y_arr=y_array[j0:j1]

    vecb_l = false
    vecb_r = false
    vecb_b = false
    vecb_t = false

    fieldtmp = zeros(gp.ny + 2, gp.nx + 2)

    #TODO distinguish u, v, w grids

    if plotbc
        if ii0 == 1
            vecb_l=true
            i1+=1
            i0tmp2-=1 

            pushfirst!(x_arr,x_array[i0]-0.5*gp.dx[1,1]/xscale)
        end

        if ii1 == gp.nx
            i1+=1
            vecb_r=true
            push!(x_arr,x_array[end]+0.5*gp.dx[1,end]/xscale)
            i1tmp2+=1 

        end

        if jj0 == 1
            vecb_b=true
            j1+=1
            j0tmp2-=1
            pushfirst!(y_arr,y_array[j0]-0.5*gp.dy[1,1]/xscale)
        end

        if jj1 == gp.ny
            vecb_t=true
            j1+=1
            push!(y_arr,y_array[end]+0.5*gp.dy[end,1]/xscale)
            j1tmp2+=1
        end

        if vecb_l 
            # @views fieldtmp[2:gp.ny+1,1] = vecb_L(field0D[i,:],gp)
            # @views vecb_L(field0D[i,:],gp)[128-68+1]=68
            # @views vecb_L(field0D[i,:],gp)[128-69+1]=69

            fieldtmp[2:gp.ny+1,1] = vecb_L(field0D[i,:],gp) #reverse(vecb_L(field0D[i,:],gp),dims=1)#[1:gp.ny]
            # @views fieldtmp[gp.ny+1,1]=1
            # @views fieldtmp[gp.ny,1]=2
            # # @views fieldtmp[68,1]=68
            # print("\n vectest ", reverse(vecb_L(field0D[i,:],gp))[65:70])

            # print("\n vectest 2 ", vecb_L(field0D[i,:],gp)[128-70+1:128-65+1])

            # print("\n vectest 2 ", fieldtmp[68:69,1])
            # vec = [1,2,3]
            # print("\n vectest ", reverse(vec), "no rev ", vec)


            # print("\n vecb", reverse(vecb_L(field0D[i,:],gp)))
        end

        if vecb_r
            @views fieldtmp[2:gp.ny+1,end] = vecb_R(field0D[i,:],gp) #right direction
        end

        if vecb_b
            @views fieldtmp[1,2:gp.nx+1] = vecb_B(field0D[i,:],gp)  #right direction
        end

        if vecb_t
            @views fieldtmp[end,2:gp.nx+1] = vecb_T(field0D[i,:],gp) #reverse(vecb_T(field0D[i,:],gp),dims=1)#[1:gp.nx]
        end
    end

    if vecb_l || vecb_r || vecb_b || vecb_t

        # printstyled(color=:green, @sprintf "\n indices : %.5i %.5i %.5i %.5i %.5i %.5i %.5i %.5i \n" ii0 ii1 jj0 jj1 i0 i1 j0 j1 )
        # printstyled(color=:green, @sprintf "\n indices 2: %.5i %.5i %.5i %.5i \n" i0tmp i1tmp j0tmp j1tmp )

        # @views fieldtmp[2:jj1+1,2:ii1+1] = field1[jj0:jj1,ii0:ii1]
        # field = fieldtmp[j0:j1,i0:i1]

        fieldtmp[j0tmp:j1tmp,i0tmp:i1tmp] = field1[jj0:jj1,ii0:ii1]
        # field = fieldtmp[j0:j1,i0:i1]
        field = fieldtmp[j0tmp2:j1tmp2,i0tmp2:i1tmp2]


        # print("\n test array ", fieldtmp[68,1], " ", fieldtmp[69,1],)
        # print("\n test array ", field[:,1])
        # print("\n test array ", fieldtmp[65:70,1])



    else
        # field = field1[j0:j1,i0:i1]
        # field1 = field0
        @views field = field0[i,j0:j1,i0:i1]

        # fieldLS = fwd.u[1,indLS,j0:j1,i0:i1]
        # i1=i1-i0+1
        # i0=1
        # j1=j1-j0+1
        # j0=1

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
        if plot_mode == "contourf"
            CS = ax2.contourf(x_arr,y_arr,field, 
            levels=levels, #10, 
            cmap=cmap)
        else

            mpl_levels = mpl_tickers.MaxNLocator(nbins=levels).tick_values(minimum(field), maximum(field))
            norm = mpl_colors.BoundaryNorm(mpl_levels, ncolors=cmap.N, clip=true)
            CS = ax2.pcolormesh(x_arr,y_arr,field, cmap=cmap, norm=norm)
        end
        # fig.colorbar(im, ax=ax0)
    end

    
    lcolor= "k" #"w" #"k"
    lw=0.5
    ms=0.5
    

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
                igrid=igrid0-i0+1 #TODO
                jgrid = jgrid0-j0+1

                if ((igrid0 == i0 && vecb_l) || (igrid0 == i1 && vecb_r) || (jgrid0 == j0 && vecb_b) || (jgrid0 == j1 && vecb_t) )
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
                if printmode == "val"
                    str=@sprintf "%.4e" field[jgrid,igrid]
                    # str=@sprintf "%.4e %.3i %.3i %.3i" field[jgrid,igrid] igrid0 jgrid0 jgrid

                elseif printmode == "ij"
                    str=@sprintf "%.3i %.3i" igrid0 jgrid0
                else 
                    str=@sprintf "%.2e %.3i %.3i" field[jgrid,igrid] igrid0 jgrid0
                end

                ax2.annotate(str,(x_arr[igrid],y_arr[jgrid]),fontsize=fontsize,c=lcolor,ha="center",va=va)

            end
        end
    end


    # ax2.set_title("Title")
    # ax2.set_xlabel(raw"$x ( \unit{\um})$")

    # ax2.set_ylabel(raw"$y ( \unit{\um})$")


    ax2.set_xlabel(raw"$x ( \unit{\um})$")
    ax2.set_ylabel(raw"$y ( \unit{\um})$")


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
     


        if plotcase=="circle"
            # plt.plot( , colors="g")
            theta1=0 #-90
            theta2=90

            radius = fwd.radius[i]/num.plot_xscale
            arc = matplotlib.patches.Arc((num.xcoord/num.plot_xscale, num.ycoord/num.plot_xscale), radius*2, radius*2, color="g", theta1=theta1, theta2=theta2,ls="--")

            ax2.add_patch(arc)

        end

        # printstyled(color=:green, @sprintf "\n j: %5i %5i %5i %5i\n" i0 i1 j0 j1)

        # CSlvl = ax2.contour(x_arr,y_arr, fwd.u[1,indLS,j0:j1,i0:i1], [0.0],colors="r")
        CSlvl = ax2.contour(x_array[ii0:ii1],y_array[jj0:jj1], fwd.u[1,indLS,jj0:jj1,ii0:ii1], [0.0],colors="r")


    end

    ax2.spines["right"].set_visible(false)
    ax2.spines["top"].set_visible(false)

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

    ax2.set_aspect("equal")

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
    ax2.spines["right"].set_visible(false)
    ax2.spines["top"].set_visible(false)

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
        ax2.set_xlabel(raw"$x ( \unit{\um})$")

        ax2.set_ylabel(raw"$y ( \unit{\um})$")


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
    ax2.set_xlabel(raw"$x ( \unit{\um})$")

    ax2.set_ylabel(raw"$y ( \unit{\um})$")


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
    # # ax2.set_xlabel(raw"$x ( \unit{\um})$")

    # # ax2.set_ylabel(raw"$y ( \unit{\um})$")


    # ax2.set_xlabel(raw"$x ( \unit{\um})$")
    # ax2.set_ylabel(raw"$y ( \unit{\um})$")


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
    # ax2.set_xlabel(raw"$x ( \unit{\um})$")

    # ax2.set_ylabel(raw"$y ( \unit{\um})$")


    ax2.set_xlabel(raw"$x ( \unit{\um})$")
    ax2.set_ylabel(raw"$y ( \unit{\um})$")


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


function debug_border_top(iplot,jplot,A,B,rhs,geo,inside,ind,grid,bc,iscal,num,op,ni,nb,ph,Bx,By)

    printstyled(color=:magenta, @sprintf "\n debug border \n")

    ny = grid.ny
    @views scalar_debug_border(geo.dcap, ny,bc[iscal], inside, ind ,num,grid,iplot,jplot,op)

    testb = iplot
    testn = 2*grid.ny+ 2* grid.nx - testb + 1

    tmp = grid.nx - testb + 1

    print("\n test",testn," testb ",testb)
    # print("\n test -1 ",ph.trans_scalD[end-nb+testn-1,iscal])
    printstyled(color=:green, @sprintf "\n jtmp : %.5i j : %.5i chi_b %.2e  chi_b adim %.2e border %.2e\n" testn testb op.χ_b[end-nb+testn,end-nb+testn] op.χ_b[end-nb+testn,end-nb+testn]/grid.dy[1,1] ph.trans_scalD[end-nb+testn,iscal])
    printstyled(color=:cyan, @sprintf "\n BC %.5e rhs %.5e rhs %.5e \n" bc[iscal].top.val bc[iscal].top.val*op.χ_b[end-nb+testn,end-nb+testn] rhs[end-nb+testn])

    # print("\n rhs ", rhs, "\n")
    # print("\n B ", maximum(B[testb,:])," \n ")

    # Border BCs
    print("\n A[end-nb+testn,1:ni]", A[end-nb+testn,1:ni], "\n")
    print("\n A[end-nb+testn,ni+1:2*ni]", A[end-nb+testn,ni+1:2*ni], "\n")
    print("\n A[end-nb+testn,end-nb+1:end]", A[end-nb+testn,end-nb+1:end], "\n")

    # pII = lexicographic(II, grid.ny)

    II = CartesianIndex(jplot, iplot) #(id_y, id_x)
    pII = lexicographic(II, grid.ny)

    print("\n pII ", pII)

    print("\n A[pII,1:ni]", A[pII,1:ni], "\n")
    print("\n A[pII,ni+1:2*ni]", A[pII,ni+1:2*ni], "\n")
    print("\n A[pII,end-nb+1:end]", A[pII,end-nb+1:end], "\n")

    print("\n A[ni+pII,1:ni]", A[ni+pII,1:ni], "\n")
    print("\n A[ni+pII,ni+1:2*ni]", A[ni+pII,ni+1:2*ni], "\n")
    print("\n A[ni+pII,end-nb+1:end]", A[ni+pII,end-nb+1:end], "\n")

    # print("\n b_b ", b_b[testb,testb]) 1 
    # print("\n a1_b ", a1_b[testb,testb]) 0

    # print("\n HxT_b", op.HxT_b[testb,:])
    # print("\n iMx_b'", op.iMx_b'[testb,:])
    # print("\n Hx_b", op.Hx_b[testb,:])

    # print("\n HyT_b", op.HyT_b[testb,:])
    # print("\n iMy_b'", op.iMy_b'[testb,:])
    # print("\n Hy_b", op.Hy_b[testb,:])


    # print("\n mult", op.HxT_b[testb,testb]*op.iMx_b'[testb,testb]*op.Hx_b[testb,testb])

    # print("\n mult", op.HyT_b[testb,testb]*op.iMy_b'[testb,testb]*op.Hy_b[testb,testb])

end

function debug_border_left(iplot,jplot,A,B,rhs,geo,inside,ind,grid,bc,iscal,num,op,ni,nb,ph,Bx,By)

    printstyled(color=:magenta, @sprintf "\n debug border \n")

    ny = grid.ny
    @views scalar_debug_border(geo.dcap, ny,bc[iscal], inside, ind ,num,grid,iplot,jplot,op)

    testb = jplot
    testn = ny-testb+1
    print("\n test",testn," testb ",testb)
    printstyled(color=:green, @sprintf "\n jtmp : %.5i j : %.5i chi_b %.2e  chi_b adim %.2e border %.2e\n" testn testb op.χ_b[end-nb+testn,end-nb+testn] op.χ_b[end-nb+testn,end-nb+testn]/grid.dy[1,1] vecb_L(ph.trans_scalD[:,iscal], grid)[testn])
    printstyled(color=:cyan, @sprintf "\n BC %.5e rhs %.5e rhs %.5e \n" bc[iscal].left.val[testn] bc[iscal].left.val[testn]*op.χ_b[end-nb+testn,end-nb+testn] vecb_L(rhs, grid)[testn])

    print("\n B ", maximum(B[testb,:])," \n ")

    print("\n A[end-nb+testn,1:ni]", A[end-nb+testn,1:ni], "\n")
    print("\n A[end-nb+testn,ni+1:2*ni]", A[end-nb+testn,ni+1:2*ni], "\n")
    print("\n A[end-nb+testn,end-nb+1:end]", A[end-nb+testn,end-nb+1:end], "\n")

    print("\n A[jplot,1:ni]", A[jplot,1:ni], "\n")
    print("\n A[jplot,ni+1:2*ni]", A[jplot,ni+1:2*ni], "\n")
    print("\n A[jplot,end-nb+1:end]", A[jplot,end-nb+1:end], "\n")

    print("\n A[ni+jplot,1:ni]", A[ni+jplot,1:ni], "\n")
    print("\n A[ni+jplot,ni+1:2*ni]", A[ni+jplot,ni+1:2*ni], "\n")
    print("\n A[ni+jplot,end-nb+1:end]", A[ni+jplot,end-nb+1:end], "\n")

    # print("\n b_b ", b_b[testb,testb]) 1 
    # print("\n a1_b ", a1_b[testb,testb]) 0

    print("\n HxT_b", op.HxT_b[testb,:])
    print("\n iMx_b'", op.iMx_b'[testb,:])
    print("\n Hx_b", op.Hx_b[testb,:])

    print("\n HyT_b", op.HyT_b[testb,:])
    print("\n iMy_b'", op.iMy_b'[testb,:])
    print("\n Hy_b", op.Hy_b[testb,:])


    print("\n mult", op.HxT_b[testb,testb]*op.iMx_b'[testb,testb]*op.Hx_b[testb,testb])

    print("\n mult", op.HyT_b[testb,testb]*op.iMy_b'[testb,testb]*op.Hy_b[testb,testb])

end

function scalar_debug_border(cap, n, BC, inside, ind,num,grid,iplot,jplot,op)

    @unpack nx, ny = grid
    @unpack b_left, b_bottom, b_right, b_top = ind

    printstyled(color=:green, @sprintf "\n scalar debug \n")

    II = CartesianIndex(jplot, iplot) #(id_y, id_x)
    pII = lexicographic(II, grid.ny)


    prefix = num.plot_prefix
    figname = "scalar_border"
    xscale = num.plot_xscale
    alpha = 0.5

    idx = ny-jplot+1

    A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, b_left[1][idx])

    printstyled(color=:cyan, @sprintf "\n c: %.2e %.2e %.2e %.2e %.2e %.2e \n" A1/xscale A2/xscale A3/xscale A4/xscale B1/xscale B2/xscale)

    # print("\n vol, ", cap[II,8:11],cap[II,8:11]./(xscale^2))

    fig1, ax2 = plt.subplots(layout="constrained")


    width = grid.dx[jplot,iplot] / xscale
    height = grid.dy[jplot,iplot] / xscale
    xc = (grid.x[jplot,iplot] -0.5*grid.dx[1,1])/xscale
    yc = grid.y[jplot,iplot]/xscale



    xy = (xc-width/2,yc-height/2)
    square = matplotlib.patches.Rectangle(xy, width, height, color="k", ls="--",alpha=alpha)
    ax2.add_patch(square)

    x1 = xc - width/2 
    x2 = xc + width/2 
    y1 = yc + width/2 - A1 / xscale
    y2 = yc + width/2 - A3 / xscale


    # print("\n x1 ",x1)
    # print("\n y1 ",y1)
    # print("\n x2 ",x2)
    # print("\n y2 ",y2)
    # print("\n point ",xc, " y " ,yc)
    # print("\n point ",xc - width/2 , " y " ,yc - width/2)


    #interface segment
    ax2.plot([x1, x2], [y1, y2], marker = "o",lw=1)

    #x
    x1 = xc - width/2 
    x2 = x1 
    y1 = yc + width/2 - A1 / xscale
    y2 = yc + width/2 
    ax2.plot([x1, x2], [y1, y2], marker = "o",lw=1)
    #x
    x1 = xc + width/2 
    x2 = x1 
    y1 = yc + width/2 - A3 / xscale
    y2 = yc + width/2 
    ax2.plot([x1, x2], [y1, y2], marker = "o",lw=1)

    #y
    x1 = xc - width/2 
    x2 = xc - width/2 + A2 /xscale 
    y1 = yc - width/2
    y2 = y1
    ax2.plot([x1, x2], [y1, y2], marker = "o",lw=1)
    #y
    x1 = xc - width/2 
    x2 = xc - width/2 + A4 /xscale 
    y1 = yc + width/2
    y2 = y1
    ax2.plot([x1, x2], [y1, y2], marker = "o",lw=1)


    #TODO plot at centroid, not xy yc
    #B
    x1 = xc 
    x2 = x1 
    y1 = yc + width/2 - B1 / xscale
    y2 = yc + width/2 
    ax2.plot([x1, x2], [y1, y2], marker = "o",lw=1)

    #B
    x1 = xc - width/2 
    x2 = xc - width/2 + B2 /xscale 
    y1 = yc 
    y2 = y1
    ax2.plot([x1, x2], [y1, y2], marker = "o",lw=1)


    ax2.scatter(grid.x[jplot,iplot]/xscale,grid.y[jplot,iplot]/xscale)


    widthvol = sqrt(cap[b_left[1][idx],8]./(xscale^2))
    heightvol = widthvol
    # xy=(xc-width/2-widthvol/2,yc)
    xy=(xc-width/2,yc)
    square = matplotlib.patches.Rectangle(xy, widthvol, heightvol, color="g", ls="--",alpha=alpha)
    ax2.add_patch(square)


    theta1=75
    theta2=95

    radius = num.R/num.plot_xscale
    arc = matplotlib.patches.Arc((num.xcoord/num.plot_xscale, num.ycoord/num.plot_xscale), radius*2, radius*2, color="g", theta1=theta1, theta2=theta2,ls="--")

    ax2.add_patch(arc)


    plt.axis("equal")

    ax2.set_aspect("equal")

    str_it = @sprintf "_%.5i" pII
    plt.savefig(prefix*figname*str_it*".pdf")
    plt.close(fig1)


end

function scalar_debug!(::Dirichlet, O, B, u, v, Dx, Dy, Du, Dv, cap, n, BC, inside, b_left, b_bottom, b_right, b_top,num,grid,iplot,jplot)

    printstyled(color=:green, @sprintf "\n scalar debug \n")

    II = CartesianIndex(jplot, iplot) #(id_y, id_x)
    pII = lexicographic(II, grid.ny)


    prefix = num.plot_prefix
    figname = "scalar_conv"
    xscale = num.plot_xscale
    alpha = 0.5
    A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II)
    u1, v2, u3, v4 = u[II], v[II], u[δx⁺(II)], v[δy⁺(II)]

    printstyled(color=:cyan, @sprintf "\n c: %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e\n" u1 v2 u3 v4 Du[II] Dv[II] A1 A2 A3 A4 B1 B2 Dx[II] Dy[II])

    printstyled(color=:cyan, @sprintf "\n c: %.2e %.2e %.2e %.2e %.2e %.2e \n" A1/xscale A2/xscale A3/xscale A4/xscale B1/xscale B2/xscale)



    print("\n vol, ", cap[II,8:11],cap[II,8:11]./(xscale^2))

    fig1, ax2 = plt.subplots(layout="constrained")


    width = grid.dx[jplot,iplot] / xscale
    height = grid.dy[jplot,iplot] / xscale
    xc = grid.x[jplot,iplot]/xscale
    yc = grid.y[jplot,iplot]/xscale
    xy = (xc-width/2,yc-height/2)
    square = matplotlib.patches.Rectangle(xy, width, height, color="k", ls="--",alpha=alpha)
    ax2.add_patch(square)

    x1 = xc - width/2 
    x2 = xc + width/2 
    y1 = yc + width/2 - A1 / xscale
    y2 = yc + width/2 - A3 / xscale


    # print("\n x1 ",x1)
    # print("\n y1 ",y1)
    # print("\n x2 ",x2)
    # print("\n y2 ",y2)
    # print("\n point ",xc, " y " ,yc)
    # print("\n point ",xc - width/2 , " y " ,yc - width/2)


    #interface segment
    ax2.plot([x1, x2], [y1, y2], marker = "o",lw=1)

    #x
    x1 = xc - width/2 
    x2 = x1 
    y1 = yc + width/2 - A1 / xscale
    y2 = yc + width/2 
    ax2.plot([x1, x2], [y1, y2], marker = "o",lw=1)
    #x
    x1 = xc + width/2 
    x2 = x1 
    y1 = yc + width/2 - A3 / xscale
    y2 = yc + width/2 
    ax2.plot([x1, x2], [y1, y2], marker = "o",lw=1)

    #y
    x1 = xc - width/2 
    x2 = xc - width/2 + A2 /xscale 
    y1 = yc - width/2
    y2 = y1
    ax2.plot([x1, x2], [y1, y2], marker = "o",lw=1)
    #y
    x1 = xc - width/2 
    x2 = xc - width/2 + A4 /xscale 
    y1 = yc + width/2
    y2 = y1
    ax2.plot([x1, x2], [y1, y2], marker = "o",lw=1)


    #TODO plot at centroid, not xy yc
    #B
    x1 = xc 
    x2 = x1 
    y1 = yc + width/2 - B1 / xscale
    y2 = yc + width/2 
    ax2.plot([x1, x2], [y1, y2], marker = "o",lw=1)

    #B
    x1 = xc - width/2 
    x2 = xc - width/2 + B2 /xscale 
    y1 = yc 
    y2 = y1
    ax2.plot([x1, x2], [y1, y2], marker = "o",lw=1)


    ax2.scatter(grid.x[jplot,iplot]/xscale,grid.y[jplot,iplot]/xscale)


    #Square with same vol as cap 8
    # widthvol = sqrt(cap[II,8]./(xscale^2))
    # heightvol = widthvol
    # xy=(xc-width/2-widthvol/2,yc)
    # square = matplotlib.patches.Rectangle(xy, widthvol, heightvol, color="g", ls="--",alpha=alpha)
    # ax2.add_patch(square)

    #plot arc
    # theta1=75
    # theta2=90
    # radius = num.R/num.plot_xscale
    # arc = matplotlib.patches.Arc((num.xcoord/num.plot_xscale, num.ycoord/num.plot_xscale), radius*2, radius*2, color="g", theta1=theta1, theta2=theta2,ls="--")
    # ax2.add_patch(arc)


    plt.axis("equal")

    ax2.set_aspect("equal")

    str_it = @sprintf "_%.5i" pII
    plt.savefig(prefix*figname*str_it*".pdf")
    plt.close(fig1)


    # B .= 0.0
    # O .= 0.0
    # @inbounds for II in inside
    #     pII = lexicographic(II, n)
    #     A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II)
    #     u1, v2, u3, v4 = u[II], v[II], u[δx⁺(II)], v[δy⁺(II)]

    #     @inbounds O[pII,pII] = 0.5 * (A3 * u3 - A1 * u1 + A4 * v4 - A2 * v2)
    #     @inbounds O[pII,pII+n] = 0.5 * A3 * u3
    #     @inbounds O[pII,pII-n] = -0.5 * A1 * u1
    #     @inbounds O[pII,pII+1] = 0.5 * A4 * v4
    #     @inbounds O[pII,pII-1] = -0.5 * A2 * v2

    #     @inbounds O[pII,pII] += -0.5 * ((A3 - B1) * Du[δx⁺(II)] + (B1 - A1) * Du[II])
    #     @inbounds O[pII,pII] += -0.5 * ((A4 - B2) * Dv[δy⁺(II)] + (B2 - A2) * Dv[II])

    #     @inbounds B[pII] += -0.5 * Dx[II] * ((A3 - B1) * Du[δx⁺(II)] + (B1 - A1) * Du[II])
    #     @inbounds B[pII] += -0.5 * Dy[II] * ((A4 - B2) * Dv[δy⁺(II)] + (B2 - A2) * Dv[II])
    # end

    # @inbounds for II in vcat(b_left, b_bottom[2:end-1], b_right, b_top[2:end-1])
    #     pII = lexicographic(II, n)
    #     A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II)
    #     u1, v2, u3, v4 = u[II], v[II], u[δx⁺(II)], v[δy⁺(II)]

    #     @inbounds O[pII,pII] = 0.5 * (A3 * u3 - A1 * u1 + A4 * v4 - A2 * v2)

    #     @inbounds O[pII,pII] += -0.5 * ((A3 - B1) * Du[δx⁺(II)] + (B1 - A1) * Du[II])
    #     @inbounds O[pII,pII] += -0.5 * ((A4 - B2) * Dv[δy⁺(II)] + (B2 - A2) * Dv[II])

    #     @inbounds B[pII] += -0.5 * Dx[II] * ((A3 - B1) * Du[δx⁺(II)] + (B1 - A1) * Du[II])
    #     @inbounds B[pII] += -0.5 * Dy[II] * ((A4 - B2) * Dv[δy⁺(II)] + (B2 - A2) * Dv[II])
    # end
    # @inbounds for II in vcat(b_left, b_bottom[2:end-1], b_top[2:end-1])
    #     pII = lexicographic(II, n)
    #     A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II)

    #     @inbounds O[pII,pII+n] = 0.5 * A3 * u[δx⁺(II)]
    # end
    # @inbounds for II in vcat(b_bottom[2:end-1], b_right, b_top[2:end-1])
    #     pII = lexicographic(II, n)
    #     A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II)

    #     @inbounds O[pII,pII-n] = -0.5 * A1 * u[II]
    # end
    # @inbounds for II in vcat(b_left[2:end-1], b_bottom, b_right[2:end-1])
    #     pII = lexicographic(II, n)
    #     A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II)

    #     @inbounds O[pII,pII+1] = 0.5 * A4 * v[δy⁺(II)]
    # end
    # @inbounds for II in vcat(b_left[2:end-1], b_right[2:end-1], b_top)
    #     pII = lexicographic(II, n)
    #     A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II)

    #     @inbounds O[pII,pII-1] = -0.5 * A2 * v[II]
    # end

    # if is_periodic(BC.left) && is_periodic(BC.right)
    #     @inbounds for (II,JJ) in zip(b_right, b_left)
    #         pII = lexicographic(II, n)
    #         pJJ = lexicographic(JJ, n)
    #         A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II)
    
    #         @inbounds O[pII,pJJ] = 0.5 * A3 * u[δx⁺(II)]
    #     end
    #     @inbounds for (II,JJ) in zip(b_left, b_right)
    #         pII = lexicographic(II, n)
    #         pJJ = lexicographic(JJ, n)
    #         A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II)
    
    #         @inbounds O[pII,pJJ] = -0.5 * A1 * u[II]
    #     end
    # end
    # if is_periodic(BC.bottom) && is_periodic(BC.top)
    #     @inbounds for (II,JJ) in zip(b_top, b_bottom)
    #         pII = lexicographic(II, n)
    #         pJJ = lexicographic(JJ, n)
    #         A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II)
    
    #         @inbounds O[pII,pJJ] = 0.5 * A4 * v[δy⁺(II)]
    #     end
    #     @inbounds for (II,JJ) in zip(b_bottom, b_top)
    #         pII = lexicographic(II, n)
    #         pJJ = lexicographic(JJ, n)
    #         A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II)
    
    #         @inbounds O[pII,pJJ] = -0.5 * A2 * v[II]
    #     end
    # end

    # @inbounds _A1 = @view cap[:,:,1]
    # @inbounds _A2 = @view cap[:,:,2]
    # @inbounds _A3 = @view cap[:,:,3]
    # @inbounds _A4 = @view cap[:,:,4]
    # @inbounds _B1 = @view cap[:,:,6]
    # @inbounds _B2 = @view cap[:,:,7]

    # set_sca_conv_bnd!(dir, BC.left, O, δx⁺, _A1, _A3, _B1, Du, n, b_left, b_right)
    # set_sca_conv_bnd!(dir, BC.bottom, O, δy⁺, _A2, _A4, _B2, Dv, n, b_bottom, b_top)
    # set_sca_conv_bnd!(dir, BC.right, O, δx⁺, _A1, _A3, _B1, Du, n, b_right, b_left)
    # set_sca_conv_bnd!(dir, BC.top, O, δy⁺, _A2, _A4, _B2, Dv, n, b_top, b_bottom)

    # return nothing
end


function plot_radial_vel()

    us,vs = interpolate_grid_liquid(gp,gu,gv,phL.u,phL.v)
   
    # us .= 0.0
    # vs .= 0.0
    # for ju in 1:gp.ny
    #     for iu in 1:gp.nx
    #         xcell = gp.x[ju,iu]
    #         ycell = gp.y[ju,iu]

    #         vec0 = [xcoord, ycoord]
    #         vec1 = [xcell, ycell]

    #         vecr = vec1-vec0
    #         normr = norm(vecr)
    #         if normr>radius
    #             vecr .*= 1.0/normr
    #             factor = 1.0/normr
    #             #factor = 1.0 
    #             us[ju,iu] = factor * vecr[1]
    #             vs[ju,iu] = factor * vecr[2]

    #             print("\n vecr",vecr," vec0 ", vec0, " vec1 ", vec1)
    #             #print("\n i ",iu," j ",ju," vec",vecr," y ",ycell," ycoord ",ycoord," v ",vs)

    #         end
    #     end
    # end

    fig, ax = plt.subplots()
    q = ax.quiver(x_array,y_array,us,vs,
    #color = "red",
    )
    # # ax.quiverkey(q, X=0.3, Y=1.1, U=10,
    # #              label='Quiver key, length = 10', labelpos='E')

    # # plt.show()
    plt.axis("equal")

    plt.savefig(prefix*"vector0.pdf")

end

function plot_vector()

    us,vs = interpolate_grid_liquid(gp,gu,gv,phL.u,phL.v)
   
    fig, ax = plt.subplots()
    q = ax.quiver(x_array,y_array,us,vs,
    #color = "red",
    )
    # # ax.quiverkey(q, X=0.3, Y=1.1, U=10,
    # #              label='Quiver key, length = 10', labelpos='E')

    # # plt.show()
    plt.axis("equal")

    plt.savefig(prefix*"vector0.pdf")

    # us,vs = interpolate_grid_liquid(gp,gu,gv,phL.u,phL.v)

    # fig, ax = plt.subplots()
    # q = ax.quiver(x_array,y_array,us,vs)
    # # ax.quiverkey(q, X=0.3, Y=1.1, U=10,
    # #              label='Quiver key, length = 10', labelpos='E')
    # plt.axis("equal")

    # # plt.show()

    # plt.savefig(prefix*"vector.pdf")
    # plt.close(fig)


    # print("\n test u ", vecb_L(phL.uD, gu))
    # print("\n test u ", phL.u[1,:])
    # print("\n test u ", phL.u[:,1])

    # vecb_L(phL.uD, gu) .=0.0
    # vecb_R(phL.uD, gu) .=0.0
    # vecb_T(phL.uD, gu) .=0.0
    # vecb_B(phL.uD, gu) .=0.0
    # vec1(phL.uD, gu) .=0.0


    # phL.u .= reshape(vec1(phL.uD,gu), gu)
    # phL.v .= reshape(vec1(phL.vD,gv), gv)

    # us,vs = interpolate_grid_liquid(gp,gu,gv,phL.u,phL.v)

    # fig, ax = plt.subplots()
    # q = ax.quiver(x_array,y_array,us,vs)
    # # ax.quiverkey(q, X=0.3, Y=1.1, U=10,
    # #              label='Quiver key, length = 10', labelpos='E')

    # plt.axis("equal")

    # plt.show()

    # plt.savefig(prefix*"vector_test1.pdf")

    # # print("\n test u ", vecb_L(phL.uD, gu))
end


#u0
# ####################################################################################################
# # u_array = phL.u 
# u_array = phL.u ./velscale

# fig1, ax2 = plt.subplots(layout="constrained")
# CS = ax2.contourf(xu,yu,u_array, 10, cmap=cmap)

# # Note that in the following, we explicitly pass in a subset of the contour
# # levels used for the filled contours.  Alternatively, we could pass in
# # additional levels to provide extra resolution, or leave out the *levels*
# # keyword argument to use all of the original levels.

# # CS2 = ax2.contour(CS, 
# # # levels=CS.levels[::2], 
# # # levels=
# # colors="r")

# # ax2.set_title("Title")
# ax2.set_xlabel(L"$x (\mu m)$")
# ax2.set_ylabel(L"$y (\mu m)$")

# # Make a colorbar for the ContourSet returned by the contourf call.
# cbar = fig1.colorbar(CS)
# cbar.ax.set_ylabel("u")
# # Add the contour line levels to the colorbar
# # cbar.add_lines(CS2)

# # cbar.formatter.set_powerlimits((0, 0))
# # # to get 10^3 instead of 1e3
# # # cbar.formatter.set_useMathText(True)
# # cbar.formatter.set_useMathText(1)

# cbar.ax.set_title(L"$10^{-4}$")

# plt.axis("equal")


# plt.savefig(prefix*"u0.pdf")
# plt.close(fig1)

# ####################################################################################################

# v_array = phL.v ./velscale

# fig1, ax2 = plt.subplots(layout="constrained")
# CS = ax2.contourf(xv,yv,phL.v ./velscale, 10, cmap=cmap)


# # ax2.set_title("Title")
# ax2.set_xlabel(L"$x (\mu m)$")
# ax2.set_ylabel(L"$y (\mu m)$")

# # Make a colorbar for the ContourSet returned by the contourf call.
# cbar = fig1.colorbar(CS)
# cbar.ax.set_ylabel("v")

# cbar.ax.set_title(L"$10^{-4}$")
# plt.axis("equal")

# plt.savefig(prefix*"v0.pdf")

# plt.close(fig1)
# ####################################################################################################





# function plot_python_pdf(field,name,plot_levelset,range)
#     # print("range",range)
#     fig1, ax2 = plt.subplots(layout="constrained")
#     # CS = ax2.contourf(x_array,y_array,max.((phL.trans_scal[:,:,1] .-c0_H2)./c0_H2,0.0), 10, cmap=cmap)
#     CS = ax2.contourf(x_array,y_array,field, 
#     levels=range, #10, 
#     cmap=cmap)

#     # ax2.set_title("Title")
#     ax2.set_xlabel(raw"$x ( \unit{\um})$")

#     ax2.set_ylabel(raw"$y ( \unit{\um})$")


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
#         ax2.set_xlabel(raw"$x ( \unit{\um})$")

#         ax2.set_ylabel(raw"$y ( \unit{\um})$")


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
#         ax2.set_xlabel(raw"$x ( \unit{\um})$")

#         ax2.set_ylabel(raw"$y ( \unit{\um})$")


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
#         ax2.set_xlabel(raw"$x ( \unit{\um})$")

#         ax2.set_ylabel(raw"$y ( \unit{\um})$")


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