import h5py
import matplotlib.pyplot as plt
import numpy as np
import os
import yaml
import matplotlib.ticker as mpl_tickers
import matplotlib.colors as mpl_colors

import sys


# from termcolor import colored

# from matplotlib import rc_params

# Run the script
#https://www.geeksforgeeks.org/python-import-module-from-different-directory/
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



def plot_all_h5():
    """
    Plot all h5 in folder
    """
    # List all files in the current directory
    all_files = os.listdir('.')

    # Filter for .h5 files
    h5_files = [file for file in all_files if file.endswith('.h5')]

    # print(sys.argv)

    try:
        yamlfile = sys.argv[1]
        if '.yml' not in yamlfile:
            yamlfile+='.yml'    
    except Exception as error:
        print(error)   
        print(colored("error", "red"))

    # with open('flower.yml', 'r') as file:
    with open(yamlfile, 'r') as file:
        yml = yaml.safe_load(file)


    # print(yml)

    mesh    = yml["flower"]["mesh"]
    plotpar = yml["plot"]

    mesh["nx"] = int(mesh["nx"])
    mesh["ny"] = int(mesh["ny"])

    xp =np.linspace(float(mesh["xmin"]),float(mesh["xmax"]),int(mesh["nx"]))
    yp =np.linspace(float(mesh["ymin"]),float(mesh["ymax"]),int(mesh["ny"]))

    dx = (float(mesh["xmax"])- float(mesh["xmin"]))/int(mesh["nx"])
    dy = (float(mesh["ymax"])- float(mesh["ymin"]))/int(mesh["ny"])

    xu = xp + dx/2
    xu = np.insert(xu,0,xp[0] - dx/2)
    # yu = yp # xv = xp 
    yv = yp + dy/2
    yv = np.insert(yv,0,yp[0] - dy/2)

    scale_x = float(plotpar["scale_x"])
    scale_y = float(plotpar["scale_y"])

    xp /= scale_x
    yp /= scale_y
    xu /= scale_x
    yv /= scale_y


    

    for file_name in h5_files:

        print(file_name)
        # Load the HDF5 file
        with h5py.File(file_name, 'r') as file:            
            print(file.keys())

            # data = file['data'][:]

            time = file['time'][()]
            print("time",time)            

            for key in file.keys():
                print(key)
                if key =="time":continue #no plot

                data = file[key][:]
                file_name_1 = key
                print("max ",np.max(data))

                if "_1D" in key:
                    nx = mesh["nx"]
                    ny = mesh["ny"]
                    if key=="u_1D":nx=nx+1
                    if key=="v_1D":ny=ny+1

                    len_data = len(data)
                    len_border = 2*nx + 2*ny
                    #TODO use yaml data selection
                    n_saved_bulk_and_interfaces = (len_data-len_border)//(nx*ny) 
                    print(len_data,len_data-len_border,n_saved_bulk_and_interfaces)

                    for ifield in range(n_saved_bulk_and_interfaces):                    

                        #Bulk value
                        field = data[ifield*ny*nx:(ifield+1)*ny*nx] 
                        #stop so +1 in python,like range 

                        field = np.reshape(field, (nx, ny))            
                        field = field.transpose()

                        file_name=file_name_1 + "_" + str(ifield) #"_bulk"
                        plot_file(field,file_name,xp,yp,xu,yv,yml,mesh,plotpar,time) 
                    
                elif key == "velocity_x": #plot vector with velocity interpolated on scalar grid
                   
                    us = file["velocity_x"][:].transpose()
                    vs = file["velocity_y"][:].transpose()
                    file_name = "vectors"

                    #  print(plotpar["figures"])
                    for figparam in plotpar["figures"]:
                        if figparam['file'] == 'velocity_vectors':
                            plot_vector(us,vs,file_name,xp,yp,plotpar,figparam)

                    
                    # plot_vector(us,vs,file_name,xp,yp,plotpar)
                   

                elif key == "i_current_x": #plot vector with velocity interpolated on scalar grid

                    Eus = file["i_current_x"][:].transpose()
                    Evs = file["i_current_y"][:].transpose()

                    nx = mesh["nx"]
                    ny = mesh["ny"]
                    ifield = 1 #bulk
                    data = file["phi_ele_1D"][:]
                    field = data[ifield*ny*nx:(ifield+1)*ny*nx] 

                    field = np.reshape(field, (nx, ny))            
                    field = field.transpose()

                    file_name = "current_lines"

                    # plot_vector(Eus,Evs,file_name,xp,yp,plotpar)   

                else:
                    if key in plotpar['no_2D_plot']:continue #no plot
                    field = data.transpose()
                   
                    plot_file(field,file_name_1,xp,yp,xu,yv,yml,mesh,plotpar,time) 
                # try:
                #     data = file[key][:]
                #     file_name = key
                #     plot_file(data,file_name,xp,yp,xu,yv,yml,mesh)     
                # except Exception as error:
                #     print(error)   
                #     print(colored("error", "red"))


def plot_file(field,file_name,xp,yp,xu,yv,yml,mesh,plotpar,time):
    """ Plot one figure for field
    """
    # Create the plot
  
    # print(field.shape[0])
    # print(field.shape[1])

    print("plotting ", file_name)

    if field.shape[0] == mesh["ny"]+1:
        print("v mesh")
        x = xp
        y = yv 
    elif field.shape[1] == mesh["nx"]+1:
        print("u mesh")
        x = xu 
        y = yp
    else: 
        print("p mesh")
        x = xp
        y = yp

    fig1, ax2 = plt.subplots(layout="constrained")
    

    plot_mode = plotpar["plot_mode"]
    levels = 10  
    # levels = 0
    # range = np.linspace(0,1e-4,10)

    scale_time = float(plotpar["scale_time"])
    scale_x = float(plotpar["scale_x"])
    cmap = plt.get_cmap(plotpar["cmap"])
    
    cbarlabel = plotpar["cbarlabel"]
    isocontour = plotpar["isocontour"]

    prefix = plotpar["prefix"]

    time /= scale_time 
    # radius /= scale_x
    i = 0

    shading="nearest"
    
    if plot_mode == "contourf":
        if levels==0:
            CS = ax2.contourf(x,y,field, 
            levels=range, #10, 
            cmap=cmap)
        else:
            CS = ax2.contourf(x,y,field, 
            levels=levels,
            cmap=cmap)           
        
    elif plot_mode == "pcolormesh":
        if levels==0:            
            norm = mpl_colors.BoundaryNorm(range, ncolors=cmap.N, clip=True)
            CS = ax2.pcolormesh(x,y,field, cmap=cmap, norm=norm,shading=shading)
        else:          
            levels = mpl_tickers.MaxNLocator(nbins=levels).tick_values(np.min(field), np.max(field))
            norm = mpl_colors.BoundaryNorm(levels, ncolors=cmap.N, clip=True)
            CS = ax2.pcolormesh(x,y,field, cmap=cmap, norm=norm,shading=shading)

    # CS = ax2.pcolormesh(field, cmap=cmap)


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
   

    # if plot_levelset:
    #     # CSlvl = ax2.contourf(xp,yp,(phL.trans_scal[:,:,1] .-c0_H2)./c0_H2, levels=0.0, cmap=cmap)
    #     # CS2 = ax2.contour(CSlvl, 
    #     # # levels=CS.levels[::2], 
    #     # # levels=
    #     # colors="r")
    #     # cbar.add_lines(CS2)
    #     # CSlvl = ax2.contour(xp,yp, gp.LS[1].u, [0.0],colors="r")
    #     # CSlvl = ax2.contour(xp,yp, fwd.u[1,i,i0:i1,j0:j1], [0.0],colors="r")
     
    #     if typeof(grid) == GridCC: #isCC(grid)
    #         CSlvl = ax2.contour(xp,yp, fwd.u[1,indLS,j0:j1,i0:i1], [0.0],colors="r")
    #     elif typeof(grid) == GridFCx: #isFCx(grid)
    #         CSlvl = ax2.contour(xp,yp, fwd.ux[1,indLS,j0:j1,i0:i1], [0.0],colors="r")
    #     elif typeof(grid) == GridFCy: #isFCy(grid)
    #         CSlvl = ax2.contour(xp,yp, fwd.uy[1,indLS,j0:j1,i0:i1], [0.0],colors="r")
        


    str_time = '{:.2e}'.format(time)
    # strrad = '{:.2e}'.format(radius)
    # striter = "{:05}".format(i)


    plt.title("t "+str_time +r"$(\unit{s})$")


    ax2.spines["right"].set_visible(False)
    ax2.spines["top"].set_visible(False)


    ax2.set_xlabel(r"$x ( \unit{\um})$")
    ax2.set_ylabel(r"$y ( \unit{\um})$")
    
    # plt.savefig(prefix+figname+str_it+img_format)

    # Save the plot as a PNG file with the same name as the .h5 file
    # plt.savefig(f"{file_name.split('.')[0]}"+img_format, dpi=300)

    # plt.savefig(f"{file_name.split('.')[0]}"+img_format)

    plt.axis("equal")
    # print("saving ", file_name)
    # plt.show()
    plt.savefig(file_name+ "." + plotpar["img_format"])


    # Close the plot to avoid overlapping issues
    # plt.close()
    plt.close(fig1)


#TODO ticks
def plot_vector(us,vs,file_name,xp,yp,plotpar,figparam):
    fig1, ax2 = plt.subplots(layout="constrained")
    scale_units=plotpar["quiver_scale_unit"]
    scale_units = None if scale_units == 'None' else scale_units

    skip_every = int(figparam['skip_every'])
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
    # # ax.quiverkey(q, X=0.3, Y=1.1, U=10,
    # #              label='Quiver key, length = 10', labelpos='E')


    ax2.set_xlabel(r"$x ( \unit{\um})$")
    ax2.set_ylabel(r"$y ( \unit{\um})$")

    # plt.axis("equal")

    # print('figparam',figparam['xlim'],figparam['ylim'])
    ax2.set_xlim([float(x0) for x0 in figparam['xlim']])
    ax2.set_ylim([float(x0) for x0 in figparam['ylim']])
    # print(xp[0],xp[-1],yp[0],yp[-1])
    # plt.axis("equal")
    # print(xp)
    # print(yp)

    # ax2.axis('square')

    ax2.set_aspect('equal', 'box')

    # # plt.show()

    plt.savefig(file_name+ "." + plotpar["img_format"])
    plt.close(fig1)

def plot_vector1(us,vs,file_name,xp,yp,plotpar):
    fig1, ax2 = plt.subplots(layout="constrained")
    scale_units=plotpar["quiver_scale_unit"]
    scale_units = None if scale_units == 'None' else scale_units

    skip_every = 4
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
    # # ax.quiverkey(q, X=0.3, Y=1.1, U=10,
    # #              label='Quiver key, length = 10', labelpos='E')


    ax2.set_xlabel(r"$x ( \unit{\um})$")
    ax2.set_ylabel(r"$y ( \unit{\um})$")


    ax2.set_xlim(xp[0],xp[-1])
    ax2.set_ylim(yp[0],yp[-1])
    print(xp[0],xp[-1],yp[0],yp[-1])

    # # plt.show()
    plt.axis("equal")

    plt.savefig(file_name+ "." + plotpar["img_format"])
    plt.close(fig1)


def plot_current_lines(phi_array,Eus,Evs,file_name,xp,yp,plotpar):
    """
    plot current lines
    """
    #https://matplotlib.org/stable/gallery/images_contours_and_fields/contourf_demo.html

    fig1, ax2 = plt.subplots(layout="constrained")
    CS = ax2.contourf(xp,yp,phi_array, 10, cmap=cmap)

    CS2 = ax2.contour(CS, 
    # levels=CS.levels[::2], 
    # levels=
    colors="r")

    # ax2.set_title("Title")
    ax2.set_xlabel(r"$x ( \unit{\um})$")
    ax2.set_ylabel(r"$y ( \unit{\um})$")

    # Make a colorbar for the ContourSet returned by the contourf call.
    cbar = fig1.colorbar(CS)
    cbar.ax.set_ylabel("Electrical potential")
    # Add the contour line levels to the colorbar
    cbar.add_lines(CS2)

    plt.streamplot(xp,yp, -Eus,-Evs,color="w")
    # Eus[last_it,:,:], Evs[last_it,:,:])#, color=(.75,.90,.93)) #do no transpose, python row major
    plt.axis("equal")

    plt.savefig(file_name+ "." + plotpar["img_format"])
    plt.close(fig1)

# def plot_radius()
#     df = pd.DataFrame([fwd.t[1:size_frame] fwd.radius[1:size_frame].*1.e6],columns=["t","r"])

#     print(df)

#     df.to_pickle(prefix*"/radius.pkl")

#     fig1, ax2 = plt.subplots(layout="constrained")

#     # print("t",fwd.t)
#     # print("radius",fwd.radius.*1.e6)
#     # print("current_i", current_i)
#     # print("radius ",fwd.radius[1:current_i+1])
#     # print("\nradius ",fwd.radius)

#     plt.plot(fwd.t[1:size_frame],fwd.radius[1:size_frame].*1.e6)

#     print()

#     ilog0 = size_frame ÷ 2 
#     ilog1 = 
#     xlog = 
#     print("\n")

#     # ax2.set_title("Title")
#     ax2.set_xlabel(L"$t (s)$")
#     # ax2.set_ylabel(L"$R (m)$")
#     ax2.set_ylabel(L"$R (\mu m)$")

#     ax2.set_xscale("log")
#     ax2.set_yscale("log")

#     plt.axis("equal")

#     plt.savefig(prefix*"R.pdf")
#     plt.close(fig1)