import h5py
import matplotlib.pyplot as plt
import numpy as np
import os
import yaml
import matplotlib.ticker as mpl_tickers
import matplotlib.colors as mpl_colors
import pandas as pd
import sys
from termcolor import colored


# from termcolor import colored

# from matplotlib import rc_params

# Run the script
# https://www.geeksforgeeks.org/python-import-module-from-different-directory/
# export PYTHONPATH="/local/home/pr277828/flower/Flower.jl/examples"
# python3 -c "import plot_flower; plot_flower.plot_all_h5()"


# Colors

# OkabeIto=["#E69F00", #0 orange clair 230, 159, 0
# "#56B4E9", #1 bleu clair 86, 180, 233
# "#009E73", #2 vert 0, 158, 115
# "#F0E442", #3 jaune 240, 228, 66
# "#0072B2", #4 bleu 0, 114, 178
# "#D55E00", #5 orange 213, 94, 0
# "#CC79A7", #6 rose 204, 121, 171
# "#000000"] #7 noir 0 0 0

# colors=OkabeIto

# colors=["#000000" for color in OkabeIto]

# colors[1]="#000000"
# colors[2]=OkabeIto[5] #bleu
# colors[3]=OkabeIto[6] #orange


# Latex

plt.rc("text", usetex=True)
# rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
# rc_params["text.latex.preamble"] = [r"\usepackage{siunitx}"]

plt.rc('text.latex', preamble=r"\usepackage{siunitx}")
#  matplotlib.verbose.level = 'debug-annoying'


#  fig, axs = plt.subplots(nvary, 2,
#                               # sharey=True,
#                               # sharex=True,
#                               sharey=sharey,
#                               sharex='col',
#          figsize=set_size(latexlinewidth, fraction=figfraction,ratio=2,nvary=nvary,ratio2=ratio2),
#          layout="constrained",
#          squeeze=False,
#          )


def set_size(width,fraction=1,ratio=1,nvary=1,ratio2=4.8/6.4,height=None):
    """Set figure dimensions to avoid scaling in LaTeX.

    Parameters
    ----------
    width: float
            Document textwidth or columnwidth in pts
    fraction: float, optional
            Fraction of the width which you wish the figure to occupy

    Returns
    -------
    fig_dim: tuple
            Dimensions of figure in inches
    """
    #  https://jwalton.info/Embed-Publication-Matplotlib-Latex/#comment-4904502735

    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio to set aesthetic figure height
    # https://disq.us/p/2940ij3
    golden_ratio = (5**.5 - 1) / 2

    if height == None:
        # Width of figure (in pts)
        fig_width_pt = width * fraction

        # Figure width in inches
        fig_width_in = fig_width_pt * inches_per_pt
        # Figure height in inches
        fig_height_in = fig_width_in * golden_ratio

        # fig_dim = (fig_width_in, fig_height_in)

        #  print('figdim')
        #  print(fig_dim,fig_width_in/fig_height_in)
        #  print(6.4, 4.8,6.4/4.8)

        fig_dim = (fig_width_in, fig_width_in*ratio2/ratio*nvary)
    else:

        fig_height_in = height * fraction

        fig_width_in = fig_height_in / golden_ratio

        fig_dim = (fig_width_in, fig_height_in*ratio2/ratio*nvary)

    return fig_dim


def plot_all_fig():
    """
    Plot all fig in YAML file
    """
    # List all files in the current directory
    all_files = os.listdir(".")

    # Filter for .h5 files
    h5_files = [file for file in all_files if file.endswith(".h5")]

    # print(sys.argv)

    try:
        yamlfile = sys.argv[1]
        if ".yml" not in yamlfile:
            yamlfile += ".yml"
    except Exception as error:
        print(error)
        print(colored("error", "red"))

    with open(yamlfile, "r") as file:
        yml = yaml.safe_load(file)

    # print(yml)

    mesh = yml["flower"]["mesh"]
    plotpar = yml["plot"]

    plotpar["scale_time"] = float(plotpar["scale_time"])

    mesh["nx"] = int(mesh["nx"])
    mesh["ny"] = int(mesh["ny"])

    mesh["xmax"] = float(mesh["xmax"])
    mesh["xmin"] = float(mesh["xmin"])

    mesh["ymax"] = float(mesh["ymax"])
    mesh["ymin"] = float(mesh["ymin"])

    mesh["dx"] = (mesh["xmax"] - mesh["xmin"]) / mesh["nx"]
    mesh["dy"] = (mesh["ymax"] - mesh["ymin"]) / mesh["ny"]

    xp = np.linspace(float(mesh["xmin"]), float(mesh["xmax"]), int(mesh["nx"]))
    yp = np.linspace(float(mesh["ymin"]), float(mesh["ymax"]), int(mesh["ny"]))

    dx = (float(mesh["xmax"]) - float(mesh["xmin"])) / int(mesh["nx"])
    dy = (float(mesh["ymax"]) - float(mesh["ymin"])) / int(mesh["ny"])

    xu = xp + dx / 2
    xu = np.insert(xu, 0, xp[0] - dx / 2)
    # yu = yp # xv = xp
    yv = yp + dy / 2
    yv = np.insert(yv, 0, yp[0] - dy / 2)

    scale_x = float(plotpar["scale_x"])
    scale_y = float(plotpar["scale_y"])

    xp /= scale_x
    yp /= scale_y
    xu /= scale_x
    yv /= scale_y

    # len_data = len(data)
    # len_border = 2*nx + 2*ny
    # # TODO use yaml data selection
    # n_saved_bulk_and_interfaces = (len_data-len_border)//(nx*ny)
    # print(len_data,len_data-len_border,n_saved_bulk_and_interfaces)

    # for ifield in range(n_saved_bulk_and_interfaces):

    #     # Bulk value: ifield =1
    #     field=veci(data,nx,ny,ifield)

    #     file_name=file_name_1 + "_" + str(ifield) #"_bulk"

    for file_name in h5_files:

        print(file_name)
        # Load the HDF5 file
        with h5py.File(file_name, "r") as file:
            print(file.keys())

            # data = file['data'][:]

            time = file["time"][()]
            nstep = file["nstep"][()]

            print("time", time, "nstep", nstep)

            # Figures defined in YAML file
            for figpar in plotpar["figures"]:

                print(colored(figpar["var"] + " " + figpar["file"], "cyan"))

                if (
                    figpar["var"] == "velocity_x"
                ):  # plot vector with velocity interpolated on scalar grid

                    plot_vector(file, xp, yp, time, nstep, plotpar, figpar)

                elif (
                    figpar["var"] == "i_current_x"
                ):  # plot vector with velocity interpolated on scalar grid

                    plot_current_lines(file, xp, yp, mesh, time, nstep, plotpar, figpar)

                elif figpar["var"] in file.keys():
                    key = figpar["var"]

                    print(key)

                    if "_1D" in key:

                        if "zoom" in figpar.keys():
                            plot_python_pdf_full2(
                                file,
                                key,
                                xp,
                                yp,
                                xu,
                                yv,
                                yml,
                                mesh,
                                time,
                                nstep,
                                plotpar,
                                figpar,
                            )

                        else:
                            plot_file(
                                file,
                                key,
                                xp,
                                yp,
                                xu,
                                yv,
                                yml,
                                mesh,
                                time,
                                nstep,
                                plotpar,
                                figpar,
                            )

                    else:
                        if key in plotpar["no_2D_plot"]:
                            continue  # no plot

                        plot_file(
                            file,
                            key,
                            xp,
                            yp,
                            xu,
                            yv,
                            yml,
                            mesh,
                            time,
                            nstep,
                            plotpar,
                            figpar,
                        )


def plot_file(
    file,
    key,
    xp,
    yp,
    xu,
    yv,
    yml,
    mesh,
    time,
    nstep,
    plotpar,
    figpar=None,
):
    """Plot one figure for field"""

    nx = mesh["nx"]
    ny = mesh["ny"]
    if key=="u_1D":
        nx=nx+1
        x_1D = xu 
        y_1D = yp
        key_LS = "levelset_u"

    elif key=="v_1D":
        ny=ny+1
        x_1D = xp
        y_1D = yv 
        key_LS = "levelset_v"

    else:
        x_1D = xp
        y_1D = yp
        key_LS = "levelset_p"


    # if field.shape[0] == mesh["ny"]+1:
    #     print("v mesh")
    #     x_1D = xp
    #     y_1D = yv
    # elif field.shape[1] == mesh["nx"]+1:
    #     print("u mesh")
    #     x_1D = xu
    #     y_1D = yp
    # else:
    #     print("p mesh")
    #     x_1D = xp
    #     y_1D = yp

    data = file[key][:]
    ifield = 1 # bulk value
    if data.ndim ==1:
        data = veci(data,nx,ny,ifield)
        field=data
    else:
        field = data.transpose()


    # file_name_1 = key
    file_name = figpar['file']

    print(key,"max ",np.max(data))

    # print(field.shape,nx,ny,len(x_1D),len(y_1D))


    # print(field.shape[0])
    # print(field.shape[1])
    if figpar == None:
        figpar = plotpar

    print("plotting ", file_name)

    # fig1, ax2 = plt.subplots(layout="constrained")
    fig1, ax2 = plt.subplots(figsize=set_size(plotpar["latex_frame_width"], fraction=plotpar["fig_fraction"],
                                              ratio=1,nvary=1,ratio2=1,height=plotpar["latex_frame_height"]),
                                              layout="constrained")
    ax2.spines["right"].set_visible(False)
    ax2.spines["top"].set_visible(False)

    if 'plot_mode' not in figpar.keys():
        figpar['plot_mode'] = plotpar['plot_mode']

    scale_time = float(plotpar["scale_time"])
    scale_x = float(plotpar["scale_x"])
    cmap = plt.get_cmap(plotpar["cmap"])

    cbarlabel = plotpar["cbarlabel"]
    isocontour = plotpar["isocontour"]

    prefix = plotpar["prefix"]

    time /= scale_time 
    # radius /= scale_x

    shading="nearest"

    print(field.shape,nx,ny,len(x_1D),len(y_1D))


    if figpar['plot_mode'] == "contourf":
        if figpar['levels']==0:
            CS = ax2.contourf(x_1D,y_1D,field, 
            levels=figpar['range'], #10, 
            cmap=plotpar['cmap'],)
        else:
            CS = ax2.contourf(x_1D,y_1D,field, 
            levels=figpar['levels'],
            cmap=plotpar['cmap'])           

    elif figpar["plot_mode"] == "pcolormesh":
        if figpar["levels"] == 0:
            norm = mpl_colors.BoundaryNorm(range, ncolors=cmap.N, clip=True)
            CS = ax2.pcolormesh(
                x_1D, y_1D, field, cmap=plotpar["cmap"], norm=norm, shading=shading
            )
        else:
            levels = mpl_tickers.MaxNLocator(nbins=figpar["levels"]).tick_values(
                np.min(field), np.max(field)
            )
            norm = mpl_colors.BoundaryNorm(levels, ncolors=cmap.N, clip=True)
            CS = ax2.pcolormesh(
                x_1D, y_1D, field, cmap=plotpar["cmap"], norm=norm, shading=shading
            )

    # CS = ax2.pcolormesh(field, cmap=plotpar['cmap'])

    # Make a colorbar for the ContourSet returned by the contourf call.
    cbar = fig1.colorbar(CS)
    cbar.ax.set_ylabel(cbarlabel)
    # Add the contour line levels to the colorbar
    if isocontour:
        CS2 = ax2.contour(CS, 
        # levels=CS.levels[::2], 
        # levels=
        colors="r")
        cbar.add_lines(CS2)

    # indLS = max(i-1,1)

    if figpar['plot_levelset']:
        # CSlvl = ax2.contourf(xp,yp,(phL.trans_scal[:,:,1] -c0_H2)./c0_H2, levels=0.0, cmap=plotpar['cmap'])
        # CS2 = ax2.contour(CSlvl,
        # # levels=CS.levels[::2],
        # # levels=
        # colors="r")
        # cbar.add_lines(CS2)
        # CSlvl = ax2.contour(xp,yp, gp.LS[1].u, [0.0],colors="r")
        # CSlvl = ax2.contour(xp,yp, fwd.u[1,i,i0:i1,j0:j1], [0.0],colors="r")

        LSdat = file[key_LS][:]
        LSdat = LSdat.transpose()
        CSlvl = ax2.contour(x_1D, y_1D, LSdat, [0.0],colors="r")

    str_time = '{:.2e}'.format(time*plotpar['scale_time'])
    # strrad = '{:.2e}'.format(radius)
    # str_iter = "{:05}".format(i)

    plt.title("t "+str_time +r"$(\unit{s})$")

    ax2.spines["right"].set_visible(False)
    ax2.spines["top"].set_visible(False)

    ax2.set_xlabel(r"$x ( \unit{\um})$")
    ax2.set_ylabel(r"$y ( \unit{\um})$")

    # plt.savefig(prefix+figname+str_iter+img_format)

    # Save the plot as a PNG file with the same name as the .h5 file
    # plt.savefig(f"{file_name.split('.')[0]}"+img_format, dpi=300)

    # plt.savefig(f"{file_name.split('.')[0]}"+img_format)

    plt.axis("equal")
    # print("saving ", file_name)
    # plt.show()
    str_nstep = str(nstep)
    plt.savefig(file_name+'_'+str_nstep+ "." + plotpar["img_format"])
    # Close the plot to avoid overlapping issues
    # plt.close()
    plt.close(fig1)


# TODO ticks
def plot_vector(file,xp,yp,time,nstep,plotpar,figpar):
    """plot vectors on scalar grid
    """

    us = file["velocity_x"][:].transpose()
    vs = file["velocity_y"][:].transpose()
    # file_name = "vectors"
    file_name = figpar['file']


    fig1, ax2 = plt.subplots(layout="constrained")
    ax2.spines["right"].set_visible(False)
    ax2.spines["top"].set_visible(False)

    scale_units=plotpar["quiver_scale_unit"]
    scale_units = None if scale_units == 'None' else scale_units

    skip_every = int(figpar['skip_every'])
    skip = (slice(None, None, skip_every), slice(None, None, skip_every))
    skip1D = slice(None, None, skip_every)

    q = ax2.quiver(xp[skip1D],yp[skip1D],us[skip],vs[skip],
    scale=float(plotpar["quiver_scale"]),
    scale_units=scale_units,
    #color = "red",
    )

    # q = ax2.quiver(xp,yp,us,vs,
    # scale=float(plotpar["quiver_scale"]),
    # scale_units=scale_units,
    # #color = "red",
    # )
    # # ax.quiverkeys(q, X=0.3, Y=1.1, U=10,
    # #              label='Quiver key, length = 10', labelpos='E')

    ax2.set_xlabel(r"$x ( \unit{\um})$")
    ax2.set_ylabel(r"$y ( \unit{\um})$")

    # plt.axis("equal")

    # print('figpar',figpar['xlim'],figpar['ylim'])
    ax2.set_xlim([float(x0) for x0 in figpar['xlim']])
    ax2.set_ylim([float(x0) for x0 in figpar['ylim']])
    # print(xp[0],xp[-1],yp[0],yp[-1])
    # plt.axis("equal")
    # print(xp)
    # print(yp)

    # ax2.axis('square')

    ax2.set_xlim([float(x0) for x0 in figpar['xlim']])
    ax2.set_ylim([float(x0) for x0 in figpar['ylim']])
    ax2.set_aspect('equal', 'box')

    # str_nstep = str1="{:05}".format(nstep)
    str_nstep = str(nstep)
    str_time = '{:.2e}'.format(time*plotpar['scale_time'])
    plt.title("t "+str_time +r"$(\unit{s})$")

    plt.savefig(file_name+'_'+str_nstep+ "." + plotpar["img_format"])
    plt.close(fig1)


def plot_current_lines(file,
    xp,
    yp,
    mesh,
    time,
    nstep,
    plotpar,
    figpar,
):
    """plot current lines
    args:

    """
    phi_array = file["i_current_x"][:].transpose()
    Eus = file["i_current_x"][:].transpose()
    Evs = file["i_current_y"][:].transpose()

    nx = mesh["nx"]
    ny = mesh["ny"]
    ifield = 1 #bulk
    data = file["phi_ele_1D"][:]

    field=veci(data,nx,ny,ifield)

    # file_name = "current_lines"
    file_name = figpar['file']



    # https://matplotlib.org/stable/gallery/images_contours_and_fields/contourf_demo.html

    fig1, ax2 = plt.subplots(layout="constrained")
    ax2.spines["right"].set_visible(False)
    ax2.spines["top"].set_visible(False)

    CS = ax2.contourf(xp, yp, phi_array, 10, cmap=plotpar["cmap"])

    CS2 = ax2.contour(
        CS,
        # levels=CS.levels[::2],
        # levels=
        colors="r",
    )

    # ax2.set_title("Title")
    ax2.set_xlabel(r"$x ( \unit{\um})$")
    ax2.set_ylabel(r"$y ( \unit{\um})$")

    # Make a colorbar for the ContourSet returned by the contourf call.
    cbar = fig1.colorbar(CS)
    cbar.ax.set_ylabel("Electrical potential")
    # Add the contour line levels to the colorbar
    cbar.add_lines(CS2)

    plt.streamplot(xp, yp, -Eus, -Evs, color="w")
    # Eus[last_it,:,:], Evs[last_it,:,:])#, color=(.75,.90,.93)) #do no transpose, python row major
    
    ax2.set_xlim([float(x0) for x0 in figpar['xlim']])
    ax2.set_ylim([float(x0) for x0 in figpar['ylim']])
    ax2.set_aspect('equal', 'box')

    str_nstep = str(nstep)
    plt.savefig(file_name + "_" + str_nstep + "." + plotpar["img_format"])
    plt.close(fig1)


def plot_radius(time_list,radius_list):
    # fwd.t[1:size_frame]
    # fwd.radius[1:size_frame].
    df = pd.DataFrame([time_list, radius_list*1.e6],columns=["t","r"])

    print(df)

    df.to_pickle("/radius.pkl")

    fig1, ax2 = plt.subplots(layout="constrained")

    # print("t",fwd.t)
    # print("radius",fwd.radius.*1.e6)
    # print("current_i", current_i)
    # print("radius ",fwd.radius[1:current_i+1])
    # print("\nradius ",fwd.radius)

    plt.plot(time_list,radius_list*1.e6)

    # print()

    # ilog0 = size_frame ÷ 2
    # ilog1 =
    # xlog =

    # ax2.set_title("Title")
    ax2.set_xlabel(r"$t (s)$")
    # ax2.set_ylabel(L"$R (m)$")
    ax2.set_ylabel(r"$R (\mu m)$")

    ax2.set_xscale("log")
    ax2.set_yscale("log")

    plt.axis("equal")

    plt.savefig("R.pdf")
    plt.close(fig1)


def plot_python_pdf_full2(
    file,
    key,
    xp,
    yp,
    xu,
    yv,
    yml,
    mesh,
    time,
    nstep,
    plotpar,
    figpar,
):
    """Plot one figure for field, with BC
    args:
    """

    nx = mesh["nx"]
    ny = mesh["ny"]
    if key == "u_1D":
        nx = nx + 1
        key_LS = "levelset_u"
        x_1D = xu
        y_1D = yp
    elif key == "v_1D":
        ny = ny + 1
        key_LS = "levelset_v"
        x_1D = xp
        y_1D = yv
    else:
        key_LS = "levelset_p"
        x_1D = xp
        y_1D = yp


    # nx = mesh["nx"]
    # ny = mesh["ny"]

    # if figpar["var"] == "u_1D":
    #     print("u mesh")
    #     nx = nx + 1
    #     x_1D = xu
    #     y_1D = yp
    # elif figpar["var"] == "v_1D":
    #     print("v mesh")
    #     ny = ny + 1
    #     x_1D = xp
    #     y_1D = yv
    # else:
    #     print("p mesh")
    #     x_1D = xp
    #     y_1D = yp

    data_1D = file[key][:]
    # file_name_1 = key
    # file_name = file_name_1
    # file_name = key
    file_name = figpar['file']

    print(key,"max ",np.max(data_1D))

    print("plotting ", file_name)

    cmap = plt.get_cmap(plotpar["cmap"])

    ii0, ii1 = figpar["zoom"][0]
    jj0, jj1 = figpar["zoom"][1]
    i0 = ii0
    i1 = ii1
    j0 = jj0
    j1 = jj1
    i0tmp = i0 + 1
    j0tmp = j0 + 1
    i1tmp = i1 + 1
    j1tmp = j1 + 1
    i0tmp2 = i0 + 1
    j0tmp2 = j0 + 1
    i1tmp2 = i1 + 1
    j1tmp2 = j1 + 1



    ifield = 1 # bulk value
    field0 = veci(data_1D,nx,ny,ifield)

    field = field0
    field1 = field0

    x_arr=x_1D[i0:i1]
    y_arr=y_1D[j0:j1]

    vecb_l = False
    vecb_r = False
    vecb_b = False
    vecb_t = False

    fieldtmp = np.zeros((ny + 2, nx + 2))

    # TODO distinguish u, v, w grids even though it is a dummy position to plot BC

    if figpar['plot_bc'] == 'True':
        if ii0 == 1:
            vecb_l=True
            i1+=1
            i0tmp2-=1 

            x_arr.insert(0,x_1D[i0]-0.5*mesh['dx']/plotpar['scale_x'])

        if ii1 == gp.nx:
            i1+=1
            vecb_r=True
            x_arr.append(x_1D[end]+0.5*mesh['dx']/plotpar['scale_x'])
            i1tmp2+=1 

        if jj0 == 1:
            vecb_b=True
            j1+=1
            j0tmp2-=1
            y_arr.insert(0,y_1D[j0]-0.5*mesh['dy']/plotpar['scale_x'])

        if jj1 == gp.ny:
            vecb_t=True
            j1+=1
            y_arr.append(y_1D[end]+0.5*mesh['dy']/plotpar['scale_x'])
            j1tmp2+=1

        if vecb_l: 
            fieldtmp[2:gp.ny+1,1] = vecb_L(data_1D,nx,ny) 
        if vecb_r:
            fieldtmp[2:gp.ny+1,end] = vecb_R(data_1D,nx,ny) 
        if vecb_b:
            fieldtmp[1,2:gp.nx+1] = vecb_B(data_1D,nx,ny)
        if vecb_t:
            fieldtmp[end,2:gp.nx+1] = vecb_T(data_1D,nx,ny)

    if vecb_l or vecb_r or vecb_b or vecb_t:
        fieldtmp[j0tmp:j1tmp,i0tmp:i1tmp] = field1[jj0:jj1,ii0:ii1]
        field = fieldtmp[j0tmp2:j1tmp2,i0tmp2:i1tmp2]
    else:
        field = field0[j0:j1,i0:i1]

    fig1, ax2 = plt.subplots(layout="constrained")
    ax2.spines["right"].set_visible(False)
    ax2.spines["top"].set_visible(False)

    # print(field.shape,nx,ny,len(x_1D),len(y_1D))


    if figpar['levels']==0: 
        CS = ax2.contourf(x_arr,y_arr,field, 
        levels=figpar['range'], 
        cmap=plotpar['cmap'])
    else:

        if figpar['plot_mode'] == "contourf":
            CS = ax2.contourf(x_arr,y_arr,field, 
            levels=figpar['levels'], #10, 
            cmap=plotpar['cmap'])
        else:

            mpl_levels = mpl_tickers.MaxNLocator(nbins=figpar['levels']).tick_values(np.min(field), np.max(field))
            norm = mpl_colors.BoundaryNorm(mpl_levels, ncolors=cmap.N, clip=True)
            CS = ax2.pcolormesh(x_arr,y_arr,field, cmap=plotpar['cmap'], norm=norm)

        # fig.colorbar(im, ax=ax0)

    lcolor= "k" #"w" #"k"
    lw=0.5
    ms=0.5

    if figpar['plot_grid']=='True':
        # for igrid in i0:i1
        #     ax2.axvline(x_1D[igrid],c=lcolor,lw=lw)
        # end
        # for igrid in j0:j1
        #     ax2.axhline(y_1D[igrid],c=lcolor,lw=lw)
        # end

        for igrid0 in range(i0,i1+1):        
            for jgrid0 in range(j0,j1+1): 
                # print("\n",x_1D[igrid],y_1D[jgrid])
                # print("\n ",igrid," ",jgrid)
                igrid=igrid0-i0+1 #TODO
                jgrid = jgrid0-j0+1

                if ((igrid0 == i0 and vecb_l) or (igrid0 == i1 and vecb_r) or (jgrid0 == j0 and vecb_b) or (jgrid0 == j1 and vecb_t) ):
                    lcolor= "k"
                else:
                    lcolor= "w"

                if igrid%2 == 0:
                    va="top"
                else:
                    va="bottom"

                ax2.scatter(x_arr[igrid],y_arr[jgrid],
                c=lcolor,
                s=ms,
                )

                # str1=@sprintf "%.2e" fieldtmp[jgrid0,igrid0]
                if figpar['printmode'] == "val":
                    str1='{:.2e}'.format(field[jgrid,igrid])
                    # str1=@sprintf "%.4e %.3i %.3i %.3i" field[jgrid,igrid] igrid0 jgrid0 jgrid

                elif figpar['printmode'] == "ij":
                    str1="{:03} {:03}".format(igrid0,jgrid0)
                else: 
                    str1="{:.2e} {:03} {:03}".format(field[jgrid,igrid],igrid0,jgrid0)                         

                ax2.annotate(str1,(x_arr[igrid],y_arr[jgrid]),fontsize=plotpar['fontsize'],c=lcolor,ha="center",va=va)

    # ax2.set_title("Title")

    str_time = '{:.2e}'.format(time*plotpar['scale_time'])
    plt.title("t "+str_time +r"$(\unit{s})$")


    ax2.set_xlabel(r"$x ( \unit{\um})$")
    ax2.set_ylabel(r"$y ( \unit{\um})$")

    # Make a colorbar for the ContourSet returned by the contourf call.
    cbar = fig1.colorbar(CS)
    cbar.ax.set_ylabel(figpar['cbarlabel'])
    # Add the contour line levels to the colorbar
    if figpar['isocontour'] == 'True':
        CS2 = ax2.contour(CS, 
        # levels=CS.levels[::2], 
        # levels=
        colors="r")
        cbar.add_lines(CS2)

    if figpar['plot_levelset']:
        # CSlvl = ax2.contourf(x_arr,y_arr,(phL.trans_scal[:,:,1] -c0_H2)./c0_H2, levels=0.0, cmap=plotpar['cmap'])
        # CS2 = ax2.contour(CSlvl,
        # # levels=CS.levels[::2],
        # # levels=
        # colors="r")
        # cbar.add_lines(CS2)
        # CSlvl = ax2.contour(x_arr,y_arr, gp.LS[1].u, [0.0],colors="r")
        # CSlvl = ax2.contour(x_arr,y_arr, fwd.u[1,i,i0:i1,j0:j1], [0.0],colors="r")

        if 'plot_case' in figpar.keys():
            if figpar['plot_case']=="circle":
                theta1=figpar['theta1']
                theta2=figpar['theta2']

                radius = fwd.radius[i]/plotpar['scale_x']
                arc = matplotlib.patches.Arc(
                    (
                        yml["flower"]["physics"]["intfc_x"] / plotpar["scale_x"],
                        yml["flower"]["physics"]["intfc_y"] / plotpar["scale_x"],
                    ),
                    radius * 2,
                    radius * 2,
                    color="g",
                    theta1=theta1,
                    theta2=theta2,
                    ls="--",
                )

                ax2.add_patch(arc)

        LSdat = file[key_LS][:]
        LSdat = LSdat.transpose()
        CSlvl = ax2.contour(x_1D[ii0:ii1],y_1D[jj0:jj1], LSdat[jj0:jj1,ii0:ii1], [0.0],colors="r")

    ax2.spines["right"].set_visible(False)
    ax2.spines["top"].set_visible(False)

    # strtitle=strtitlefunc(indLS,fwd)
    # isnap = indLS
    # strtime = @sprintf "%.2e" fwd.t[isnap]*1e3
    # strrad = @sprintf "%.2e" fwd.radius[isnap]*1e6

    str_time = '{:.2e}'.format(time*plotpar['scale_time'])
    # strrad = '{:.2e}'.format(radius)
    # str_iter = "{:05}".format(i)

    # plt.title("t "*strtime*r"$( \unit{\ms})$"*"radius "*strrad*r"$( \unit{\um})$")

    ax2.set_xlim([float(x0) for x0 in figpar['xlim']])
    ax2.set_ylim([float(x0) for x0 in figpar['ylim']])
    ax2.set_aspect('equal', 'box')
    str_nstep = str(nstep)
    plt.savefig(file_name+'_'+str_nstep+ "." + plotpar["img_format"])

    plt.close(fig1)


def python_movie_zoom(
    field0,
    figname,
    prefix,
    plot_levelset,
    isocontour,
    levels,
    range,
    cmap,
    xp,yp,xu,yv,
    gp,
    cbarlabel,
    size_frame,
    i0,
    i1,
    j0,
    j1,
    fwd,
):

    fig1, ax2 = plt.subplots(layout="constrained")
    ax2.spines["right"].set_visible(False)
    ax2.spines["top"].set_visible(False)

    x_arr = x_1D[i0:i1]
    y_arr = y_1D[j0:j1]

    field = field0[:, j0:j1, i0:i1]

    if figpar['plot_mode'] == "contourf":
        contourf = True
        pcolormesh = False
        if figpar['levels']==0:
            CS = ax2.contourf(x_arr,y_arr,field[1,:,:], 
            levels=figpar['range'], #10, 
            cmap=plotpar['cmap'])
        else:
            CS = ax2.contourf(x_arr,y_arr,field[1,:,:], 
            levels=figpar['levels'],
            cmap=plotpar['cmap'])           

    elif figpar['plot_mode'] == "pcolormesh":
        contourf = False
        pcolormesh = True
        if figpar['levels']==0:         
            norm = mpl_colors.BoundaryNorm(range, ncolors=cmap.N, clip=True)
            CS = ax2.pcolormesh(x_arr,y_arr,field[1,:,:], cmap=plotpar['cmap'], norm=norm)
        else:          
            mpl_levels = mpl_tickers.MaxNLocator(nbins=figpar['levels']).tick_values(np.min(field[1,:,:]), np.max(field[1,:,:]))
            norm = mpl_colors.BoundaryNorm(mpl_levels, ncolors=cmap.N, clip=True)
            CS = ax2.pcolormesh(x_arr,y_arr,field[1,:,:], cmap=plotpar['cmap'], norm=norm)

    # Make a colorbar for the ContourSet returned by the contourf call.
    cbar = fig1.colorbar(CS)
    cbar.ax.set_ylabel(cbarlabel)

    def make_frame(i):
        ax2.clear()

        if figpar['plot_mode'] == "contourf":
            if figpar['levels']==0:
                CS = ax2.contourf(x_arr,y_arr,field[i+1,:,:], 
                levels=figpar['range'], #10, 
                cmap=plotpar['cmap'])
            else:
                CS = ax2.contourf(x_arr,y_arr,field[i+1,:,:], 
                levels=figpar['levels'],
                cmap=plotpar['cmap'])           

        elif figpar['plot_mode'] == "pcolormesh":
            if figpar['levels']==0:           
                norm = mpl_colors.BoundaryNorm(range, ncolors=cmap.N, clip=True)
                CS = ax2.pcolormesh(x_arr,y_arr,field[i+1,:,:], cmap=plotpar['cmap'], norm=norm)
            else:         
                mpl_levels = mpl_tickers.MaxNLocator(nbins=figpar['levels']).tick_values(np.min(field[i+1,:,:]), np.max(field[i+1,:,:]))
                norm = mpl_colors.BoundaryNorm(mpl_levels, ncolors=cmap.N, clip=True)
                CS = ax2.pcolormesh(x_arr,y_arr,field[i+1,:,:], cmap=plotpar['cmap'], norm=norm)

        indLS = i 
        if i ==0: 
            indLS = 1

        strtitle=strtitlefunc(indLS,fwd)
        plt.title(strtitle)

        plt.axis("equal")

        ax2.set_xlabel(r"$x ( \unit{\um})$")
        ax2.set_ylabel(r"$y ( \unit{\um})$")

        # Make a colorbar for the ContourSet returned by the contourf call.
        #  cbar = fig1.colorbar(CS)
        #  cbar.ax.set_ylabel(cbarlabel)

        if step!=0:
            fig1.colorbar(CS,cax=cbar.ax)

        # https://stackoverflow.com/questions/5180518/duplicated-colorbars-when-creating-an-animation

        if (i==0): #(i+1==1) python starts at 0
            # # Make a colorbar for the ContourSet returned by the contourf call.
            # cbar = fig1.colorbar(CS)
            # cbar.ax.set_ylabel(cbarlabel)

            # Add the contour line levels to the colorbar
            if isocontour:
                CS2 = ax2.contour(CS, 
                # levels=CS.levels[::2], 
                # levels=
                colors="r")
                cbar.add_lines(CS2)

        if plot_levelset:
            # CSlvl = ax2.contour(x_arr,y_arr, fwd.u[1,i+1,i0:i1,j0:j1], [0.0],colors="r")
            CSlvl = ax2.contour(x_arr,y_arr, fwd.u[1,indLS,j0:j1,i0:i1], [0.0],colors="r")
            # print("check levelset")

        # plt.savefig(prefix*"H2.pdf")
        # plt.close(fig1)

    myanim = anim.FuncAnimation(fig1, make_frame, frames=size_frame, interval=size_frame, blit=False)

    # myanim.save(prefix*"test.gif")
    myanim.save(prefix*figname*".mp4")

    plt.close("all")


def plot_current_wall():

    fig, ax = plt.subplots(layout="constrained")

    # fig.subplots_adjust(right=0.75)

    # varx = [0, 1, 2]
    varx = gp.y[:,1]/plotpar['scale_x']
    # vecb_L(phL.trans_scalD[:,3], gp)
    label1 = r"$c\left(H_2O\right)$"
    label2 = r"$-\eta ~\text{(-overpotential)}$ "
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
    p2, = twin1.plot(varx, phL.phi_ele[:,1] - phi_ele1, colors[2], label=label2,ls=ls2)
    p3, = twin2.plot(varx, phL.i_current_mag[:,1], colors[3], label=label3,ls=ls3)

    # p1, = ax.plot(varx, ones(gp.ny),colors[1], label=label1,ls=ls)
    # p2, = twin1.plot(varx, ones(gp.ny), colors[2], label=label2,ls=ls2)
    # p3, = twin2.plot(varx, ones(gp.ny), colors[3], label=label3,ls=ls3)

    ax.set(
        # xlim=(0, 2),
        # ylim=(0, 2),
        xlabel=r"$y ( \unit{\um})$",
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


def plot_bc(iter_list,vec,grid,plotpar,figname,prefix,time):

    fig, ax = plt.subplots(layout="constrained")

    varx = grid.x[1,:]/plotpar['scale_x']
    # print("\n varx ",varx)

    print(iter_list)
    for isnap in iter_list:

        str_time = '{:.2e}'.format(time*plotpar['scale_time'])  #TODO scale time
        str1 = "t "+str_time+r"$( \unit{\ms})$"
        # print("\n vec ",vecb_T(vec[i,:],nx,ny))
        plt.plot(varx,vecb_T(vec[isnap,:],nx,ny),label = str1)
    
    eps=1e-4
    ylim0=0.16-eps
    ylim1=0.16+eps

    ax.set_ylim(ylim0,ylim1)
    plt.legend()
    plt.savefig(prefix*figname*".pdf")

    plt.close(fig)


def plot_bc2(iter_list,vec,grid,plotpar,figname,prefix,time):

    fig, ax = plt.subplots(layout="constrained")

    varx = grid.x[1,:]/plotpar['scale_x']
    # print("\n varx ",varx)

    print(iter_list)
    for isnap in iter_list:
        str_time = '{:.2e}'.format(time*plotpar['scale_time'])  #TODO scale time
        str1 = "t "+str_time+r"$( \unit{\ms})$"
        plt.plot(varx,vecb_B(vec[isnap,:],nx,ny),label = str1)
    
    eps=1e-4
    ylim0=0.16-eps
    ylim1=0.16+eps

    ax.set_ylim(ylim0,ylim1)
    plt.legend()
    plt.savefig(prefix*figname*".pdf")

    plt.close(fig)


# def example_docstring_function(a: str, b:bool, c:int):
#     ''' the function doc string.
#     Here is a bit more.

#     args:
#         a (str): a random string goes here
#         b (bool): lets describe something binary
#         c (int): we have a whole number

#     return:
#         gives something back


#     '''

#     a = a + ' world'
#     b = 5 * b
#     c = 10 + c


#     return c


def veci(data,nx,ny,ifield):
    """Returns ith field stored in the 1D vector like in Flower.jl code
    
    args:
        ifield (int): index of bulk or interface data
    return: 
        bulk (i=1) or i-th interface field
    """
    field = data[ifield*ny*nx:(ifield+1)*ny*nx] 
    field = np.reshape(field, (nx, ny))            
    field = field.transpose()
    return field


def vecb(data, nx, ny):
    """BC at border"""
    # where {G<:Grid} = @view a[end-2*g.ny-2*g.nx+1:end]
    # return data[-1-2*ny-2*nx+1:-1]
    return data[-2 * ny - 2 * nx : -1]


def vecb_L(data, nx, ny):
    """BC at left border"""
    data = vecb(data, nx, ny)
    return data[0 : ny - 1]


def vecb_B(data, nx, ny):
    """BC at bottom border"""
    data = vecb(data, nx, ny)
    return data[ny : nx + ny - 1]


def vecb_R(data, nx, ny):
    """BC at right border"""
    data = vecb(data, nx, ny)
    return data[nx + ny : 2 * ny + nx - 1]


def vecb_T(data, nx, ny):
    """# BC at top border"""
    data = vecb(data, nx, ny)
    return data[2 * ny + nx : 2 * nx + 2 * ny - 1]


# vecb_L(a,g::G) where {G<:Grid} = @view vecb(a, g)[1:g.ny]
# vecb_B(a,g::G) where {G<:Grid} = @view vecb(a, g)[g.ny+1:g.ny+g.nx]
# vecb_R(a,g::G) where {G<:Grid} = @view vecb(a, g)[g.ny+g.nx+1:2*g.ny+g.nx]
# vecb_T(a,g::G) where {G<:Grid} = @view vecb(a, g)[2*g.ny+g.nx+1:2*g.ny+2*g.nx]
