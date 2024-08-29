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
import matplotlib.animation as animation
import functools

import matplotlib.font_manager as fm

import matplotlib.ticker as mticker


#Find font path

# #!/usr/bin/env python3

# import matplotlib.font_manager
# from PIL import ImageFont

# # Iterate over all font files known to matplotlib
# for filename in matplotlib.font_manager.findSystemFonts(): 
#     # Avoid these two trouble makers - don't know why they are problematic
#     if "Emoji" not in filename and "18030" not in filename:
#         # Look up what PIL knows about the font
#         font = ImageFont.FreeTypeFont(filename)
#         name, weight = font.getname()
#         print(f'File: {filename}, fontname: {name}, weight: {weight}')

##########################################################################
plt.rcParams['text.usetex'] = True
plt.rcParams["font.size"] = "11"

plt.rc('text.latex', preamble="\n".join([ # plots will use this preamble
        r'\usepackage{amsmath}',
        r'\usepackage{booktabs}',
        r"\usepackage{siunitx}",
        r"\setlength{\abovedisplayskip}{0pt}",
        r" \setlength{\belowdisplayskip}{0pt}",
        r"\setlength{\belowdisplayshortskip}{0pt}",
        r"\setlength{\abovedisplayshortskip}{0pt}",
        r"\addtolength{\jot}{-4pt}",
        r"\usepackage{mhchem}",
       ])
)

 #   r"\listfiles"
      #   r"\sisetup{"
      #   r"table-alignment-mode = format,"
      #   r"table-number-alignment = center,"
      #   r"table-auto-round"
      #   r"}"


# "pgf.preamble": "\n".join([ # plots will use this preamble
#         r"\usepackage[utf8]{inputenc}",
#         r"\usepackage[T1]{fontenc}",
#         r"\usepackage[detect-all,locale=DE]{siunitx}"


# import os
# print(os.environ['PATH'])


# fontpath='/mnt/c/Program Files/MiKTeX/fonts/opentype/public/tex-gyre/texgyrepagella-regular.otf'


# /gpfs/workdir/regnaultp/latex/texmf-dist/fonts/opentype/public/tex-gyre-math

def apply_font(fontpath1):
    fontpath2 = 'public/tex-gyre/texgyrepagella-regular.otf'
    fontpath = fontpath1 + fontpath2

    fe = fm.FontEntry( 
    fname=fontpath,
    name='TeX Gyre Pagella Math'
    )
    fm.fontManager.ttflist.insert(0, fe) # or append is fine
    plt.rcParams['font.family'] = fe.name # = 'your custom ttf font name'

    prop = fm.FontProperties(fname=fontpath)

try:
    apply_font('/usr/share/fonts/opentype/')
    
except:
    apply_font('/gpfs/workdir/regnaultp/latex/texmf-dist/fonts/opentype/')



# fontpath2 = 'public/tex-gyre/texgyrepagella-regular.otf'
# fontpath = fontpath1 + fontpath2

# fe = fm.FontEntry( 
# fname=fontpath,
# name='TeX Gyre Pagella Math'
# )
# fm.fontManager.ttflist.insert(0, fe) # or append is fine
# plt.rcParams['font.family'] = fe.name # = 'your custom ttf font name'

# prop = fm.FontProperties(fname=fontpath)


# fe = fm.FontEntry( 
# fname=fontpath,
# name='TeX Gyre Pagella Math'
# )
# fm.fontManager.ttflist.insert(0, fe) # or append is fine
# plt.rcParams['font.family'] = fe.name # = 'your custom ttf font name'

# prop = fm.FontProperties(fname=fontpath)





##########################################################################

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


def compute_slope(ax,xls,yls,x,y,slopes,R2,param_line,colors,alpha):
   #https://math.stackexchange.com/questions/3500898/understanding-the-least-squares-regression-formula
   #https://towardsdatascience.com/linear-regression-using-least-squares-a4c3456e8570
   # print('least-squares',len(xls),len(yls))

   # print(xls)
   # print(yls)
   logslope=True
   if logslope:
      X_mean1 = np.mean(xls)
      Y_mean1 = np.mean(yls)
      xls=np.log10(xls)
      yls=np.log10(yls)

   X_mean = np.mean(xls)
   Y_mean = np.mean(yls)
   X_min = np.min(xls)
   Y_min = np.min(yls)
   X_max = np.max(xls)
   Y_max = np.max(yls)

   num = 0
   den = 0
   for i in range(len(xls)):
      num += (xls[i] - X_mean)*(yls[i] - Y_mean)
      den += (xls[i] - X_mean)**2
   m = num / den
   c = Y_mean - m*X_mean

   Y_pred = m*xls + c

   rms=np.linalg.norm(yls-Y_pred, ord=2)

   rmsx=np.linalg.norm(xls-X_mean, ord=2)
   rmsy=np.linalg.norm(yls-Y_mean, ord=2)

   corr=num/(rmsx*rmsy)
   # print('corr',corr)

   if logslope:
      # ax.plot([10**(X_mean), 10**(max(xls))], [10**(X_mean*m+c), 10**(max(Y_pred))], color='black',alpha=0.5)
      line1,=ax.plot([10**(min(xls)), 10**(max(xls))], [10**(min(Y_pred)), 10**(max(Y_pred))],
                                 # color='black',
                                 color=colors,
                                 alpha=alpha,
                                 # label=str(m)
                                 )
      xy=(X_mean1,Y_mean1)
      text='Slope={:.2g}\nR²={:.2g}'.format(float(m),float(rms))
      ax.annotate(text=text,xy=xy,ha='left',va='top')

      
      # param_line.append(line1)
    #   ax.annotate('Slope '+"{:.1f}".format(m),xy=(X_mean1,Y_mean1))

    

      # print ('least-squares',x,y,m,c,10**(X_mean),10**(Y_mean))

      # ax.scatter(x=10**(X_mean),y=10**(Y_mean))
      # ax.scatter(x=10**(X_min),y=10**(Y_min))
      # ax.scatter(x=10**(X_max),y=10**(Y_max))

      # print(xls)
      # print(yls)


      # plt.text(0.8,0.9,
      # 'Slope={:.2g}\nR²{:.2g}'.format(float(m),float(rms)),
      # transform =ax.transAxes,fontsize=16,ha='center')
      slopes=m
      R2=corr#rms

      # ax.annotate('Slope={:.2g}\nR²={:.2g}'.format(float(m),float(rms)),
      # xy=(0.8,0.1),xycoords='axes fraction',ha='center')


      # plt.legend(('data', 'line-regression r={}'.format(r_value)), 'best')
      # legendslope=plt.legend(handles=[],labels=('Slope={:.2g}\nR²{:.2g}'.format(float(m),float(rms))), loc='best')
      # ax.add_artist(legendslope)

   else:   
      if (min(Y_pred)<0):
         ax.plot([X_mean, max(xls)], [X_mean*m+c, max(Y_pred)], color='black',
               #   alpha=0.5,
                  lw=lw) # predicted
      else:
         ax.plot([min(xls), max(xls)], [min(Y_pred), max(Y_pred)], color='black',
               #   alpha=0.5,
                  lw=lw) # predicted
   
      # ax.annotate("{:.1f}".format(m)+'*x'+"{:.3f}".format(c),xy=(X_mean,Y_mean))
      ax.annotate('Slope '+"{:.1f}".format(m),xy=(X_mean,Y_mean))

   print ('least-squares',x,y,m,c,corr)#,10**(X_mean),10**(Y_mean))
   
   # print(X_mean,Y_mean,
   #       # len(xls),len(yls)
   #       ,min(xls),max(xls),min(Y_pred),max(Y_pred))

   print(min(Y_pred),max(Y_pred),Y_min,Y_max,10**(min(xls)),10**(max(xls)))   
   return(ax)


def plot_radius_from_h5():
    """
    Plot radius from h5 files, with slope
    """

    # print('arg', len(sys.argv),sys.argv)
    if len(sys.argv) == 2:
        # List all files in the current directory
        all_files = os.listdir(".")
        h5_files = [file for file in all_files if file.endswith(".h5")]
    else:
        h5_files = sys.argv[2::]

    print(h5_files)
    h5_files = sorted(h5_files)
    print('\n sorted \n')
    print(h5_files)

    nsteps = len(h5_files)
    time_list =[]
    radius_list=[]
    for file_name in h5_files:
        with h5py.File(file_name, "r") as file:
            print(file.keys())
            time = file["time"][()]
            radius = file["radius"][()]
            print(time,radius)
            time_list.append(time)
            radius_list.append(radius)

    df = pd.DataFrame({'t': time_list, 'r': radius_list})

    # new_row = {'t': time_list, 'r': radius_list}
    # df._append(new_row,ignore_index=True)

    # df = pd.DataFrame(columns=['t','r'])

    print(df)

    plot_radius_from_pandas(df)

        
def plot_radius_from_pkl():
    """
    Plot radius from pkl file, with slope
    """

    df = pd.read_pickle("./radius.pkl")

    print(df)

    plot_radius_from_pandas(df)


def plot_radius_from_pandas(df):
    """
    Plot radius pandas DF, with slope
    """

    fig1, ax2 = plt.subplots(layout="constrained")

    # print("t",fwd.t)
    # print("radius",fwd.radius.*1.e6)
    # print("current_i", current_i)
    # print("radius ",fwd.radius[1:current_i+1])
    # print("\nradius ",fwd.radius)

    color="#4d5156"

    plt.plot(df["t"],df["r"],color=color)


    xls=df["t"].to_numpy()
    yls=df["r"].to_numpy()

    # xls=xls[1:]
    # yls=yls[1:]

    tstart=2e-4
    tstop=1e-3
    istart = 0
    istop = 0

    print("first")
    for i in range(len(xls)):
        print(i,xls[i],tstart,tstop)
        if xls[i]>tstart:
            istart=i
            break

    print("second")

    for i in range(len(xls)):
        print(i,xls[i],tstart,tstop)
        if xls[i]>tstop:
            istop=i
            break
        istop=i
    print(istart,istop)

    print(xls)

    x="t"
    y="r"

    param_line=[]
    slopes=0
    R2=0
    colors="k"
    colors="#fa8b2b"
    alpha=1
    # \definecolor{mdlsorange}{HTML}{fa8b2b}

    # % \definecolor{mdlsgray}{HTML}{8e98a4}
    # \definecolor{mdlsgray}{HTML}{4d5156}



    xls=xls[istart:istop]
    yls=yls[istart:istop]

    compute_slope(ax2,xls,yls,x,y,slopes,R2,param_line,colors,alpha)


    # ilog0 = size_frame ÷ 2 
    # ilog1 = 
    # xlog = 
    # print("\n")

    # ax2.set_title("Title")
    ax2.set_xlabel(r"$t (s)$")
    # ax2.set_ylabel(L"$R (m)$")
    ax2.set_ylabel(r"$R (\mu m)$")

    ax2.set_xscale("log")
    ax2.set_yscale("log")



    plt.axis("equal")

    prefix="./"

    plt.savefig(prefix+"R.pdf")
    plt.close(fig1)



def parse_is_true(val):
    return bool(val)

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

        fig_height_pt = height * fraction

        # Figure height in inches
        fig_height_in = fig_height_pt * inches_per_pt

        fig_width_in = fig_height_in / golden_ratio

        fig_dim = (fig_width_in, fig_height_in*ratio2/ratio*nvary)

        # print('height mode',fig_dim,ratio,ratio2,nvary,width,height,fig_height_in/inches_per_pt,fraction)

    return fig_dim


def plot_all_fig():
    """
    Plot all fig in YAML file
    """

    # print('arg', len(sys.argv),sys.argv)
    if len(sys.argv) == 2:
        # List all files in the current directory
        all_files = os.listdir(".")
        h5_files = [file for file in all_files if file.endswith(".h5")]
    else:
        h5_files = sys.argv[2::]
          

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

    plotpar["scale_x"] = float(plotpar["scale_x"])
    plotpar["scale_y"] = float(plotpar["scale_y"])

    xp /= scale_x
    yp /= scale_y
    xu /= scale_x
    yv /= scale_y

    # len_data = len(data)
    # len_border = 2*nx + 2*ny
    # # TODO use yaml data selection
    # n_saved_bulk_and_interfaces = (len_data-len_border)//(nx*ny)
    # print(len_data,len_data-len_border,n_saved_bulk_and_interfaces)

    # for field_index in range(n_saved_bulk_and_interfaces):

    #     # Bulk value: field_index =1
    #     field=veci(data,nx,ny,field_index)

    #     file_name=file_name_1 + "_" + str(field_index) #"_bulk"

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

                    plot_vector(file, xp, yp, time, nstep, yml, plotpar, figpar)

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


def plot_all_films():
    """
    Plot all films in YAML file
    """

    # print('arg', len(sys.argv),sys.argv)
    if len(sys.argv) == 2:
        # List all files in the current directory
        all_files = os.listdir(".")
        h5_files = [file for file in all_files if file.endswith(".h5")]
    else:
        h5_files = sys.argv[2::]
          
    
    # print(h5_files)
    h5_files = sorted(h5_files)
    print(h5_files)

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

    for figpar in plotpar["films"]:

        key = figpar['var']
        python_movie_zoom(
        h5_files,
        key,
        xp,
        yp,
        xu,
        yv,
        yml,
        mesh,
        plotpar,
        figpar,
        )



def plot_all_films_func():
    """
    Plot all films in YAML file
    """

    # print('arg', len(sys.argv),sys.argv)
    if len(sys.argv) == 2:
        # List all files in the current directory
        all_files = os.listdir(".")
        h5_files = [file for file in all_files if file.endswith(".h5")]
    else:
        h5_files = sys.argv[2::]
          
    
    # print(h5_files)
    h5_files = sorted(h5_files)
    print(h5_files)

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

    # func = getattr('plot_file')

    func = globals()['plot_file']

    for figpar in plotpar["films"]:

        key = figpar['var']
        python_movie_zoom_func(
        h5_files,
        key,
        xp,
        yp,
        xu,
        yv,
        yml,
        mesh,
        func,
        plotpar,
        figpar,
        )


# def plot_segments(intfc_vtx_x,intfc_vtx_y,intfc_vtx_field,intfc_vtx_connectivities,intfc_vtx_num):

#     #plot scatter

#     ax2.scatter(x=intfc_vtx_x,y=intfc_vtx_y,c=intfc_vtx_field)

#     # order segments





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
    mode='close',
    fig1=None,
    ax2=None,
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

    data = file[key][:]
    
    if 'field_index' in figpar.keys():
        field_index = figpar['field_index']
    else:
        field_index = 1 # bulk value

    print(key,nstep,time,"max ",np.max(data),'min',np.min(data))

    # print(data.shape)

    if data.ndim ==1:
        data = veci(data,nx,ny,field_index)
        field=data
    elif data.shape[1] == 1:
        data = veci(data[:,0],nx,ny,field_index)
        field=data
    elif data.shape[0] == 1:

        # print(key,"max ",np.max(data),'min',np.min(data))

        # field_index = 2
        # np.set_printoptions(threshold=sys.maxsize)
        # # print(data)
        # print(data[0,:])

        data = veci(data[0,:],nx,ny,field_index)

        field=data
        # np.set_printoptions(threshold=sys.maxsize)
        # print(data)
        
    else:
        field = data.transpose()

    # file_name_1 = key
    file_name = figpar['file']

    # print(key,"max ",np.max(data),'min',np.min(data))

    if figpar == None:
        figpar = plotpar

    if mode == 'film' or mode == 'first':
        ax2.clear()
    else:
        # fig1, ax2 = plt.subplots(layout="constrained")
        fig1, ax2 = plt.subplots(figsize=set_size(plotpar["latex_frame_width"], fraction=float(plotpar["fig_fraction"]),
                                                ratio=1,nvary=1,ratio2=1,height=float(plotpar["latex_frame_height"])),
                                                layout="constrained")


    if 'plot_mode' not in figpar.keys():
        figpar['plot_mode'] = plotpar['plot_mode']

    scale_time = float(plotpar["scale_time"])
    scale_x = float(plotpar["scale_x"])
    cmap = plt.get_cmap(plotpar["cmap"])

    cbarlabel = plotpar["cbarlabel"]
    isocontour = plotpar["isocontour"]

    time /= scale_time 
    # radius /= scale_x

    shading="nearest"

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

    # Make a colorbar for the ContourSet returned by the contourf call.
    if mode !='film':
        cbar = fig1.colorbar(CS)
        cbar.ax.set_ylabel(r""+cbarlabel)
    # Add the contour line levels to the colorbar
    if isocontour:
        CS2 = ax2.contour(CS, 
        # levels=CS.levels[::2], 
        # levels=
        colors="r")
        cbar.add_lines(CS2)

    if figpar['plot_levelset']:
        LSdat = file[key_LS][:]
        LSdat = LSdat.transpose()
        CSlvl = ax2.contour(x_1D, y_1D, LSdat, [0.0],colors="r",linewidths=figpar['linewidth'],linestyles=figpar['linestyle'])

    if figpar['plot_levelset_segments']:
        intfc_vtx_num = file["intfc_vtx_num"][()]
        print("intfc_vtx_num ",intfc_vtx_num)

        intfc_vtx_x = file["intfc_vtx_x"][:]
        intfc_vtx_y = file["intfc_vtx_y"][:]
        intfc_vtx_field = file["intfc_vtx_field"][:]
        intfc_vtx_connectivities = file["intfc_vtx_connectivities"][:]

        intfc_vtx_x /= plotpar['scale_x']
        intfc_vtx_y /= plotpar['scale_y']

        print(intfc_vtx_x)
        print(intfc_vtx_y)
        print(intfc_vtx_field)

        ms = 0.2 #0.5

        ax2.scatter(x=intfc_vtx_x,y=intfc_vtx_y,
                    # c=intfc_vtx_field,
                    s=ms,
                    )

        # plot_segments(intfc_vtx_x,intfc_vtx_y,intfc_vtx_field,intfc_vtx_connectivities,intfc_vtx_num):

#     #plot scatter

#     ax2.scatter(x=intfc_vtx_x,y=intfc_vtx_y,c=intfc_vtx_field)

#     # order segments





    str_time = '{:.2e}'.format(time/plotpar['scale_time'])
    # strrad = '{:.2e}'.format(radius)
    # str_iter = "{:05}".format(i)

    # plt.title("t "+str_time +r"$(\unit{s})$")

    plt.title('Time '+r"$\SI{{{0:.2e}}}".format(time/plotpar['scale_time'])+'{'+plotpar['unit_time']+'}$')

    if mode =='first' or mode =='close':
        ax2.spines["right"].set_visible(False)
        ax2.spines["top"].set_visible(False)

        ax2.set_xlabel(r"$x ( \unit{\um})$")
        ax2.set_ylabel(r"$y ( \unit{\um})$")

        ax2.set_xlim([float(x0) for x0 in figpar['xlim']])
        ax2.set_ylim([float(x0) for x0 in figpar['ylim']])
        ax2.set_aspect('equal', 'box')

    if mode == 'close':
        str_nstep = str(nstep)
        plt.savefig(file_name+'_'+str_nstep+ "." + plotpar["img_format"],dpi=plotpar['dpi'])
        plt.close(fig1)
        return

    return(fig1,ax2,CS)


# TODO ticks
def plot_vector(file,x_1D,y_1D,time,nstep,yml,plotpar,figpar):
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

    q = ax2.quiver(x_1D[skip1D],y_1D[skip1D],us[skip],vs[skip],
    scale=float(plotpar["quiver_scale"]),
    scale_units=scale_units,
    angles=scale_units,
    #color = "red",
    )

    # txt1 = "My name is {fname}, I'm {age}".format(fname = "John", age = 36)
    # txt2 = "My name is {0}, I'm {1}".format("John",36)
    # txt3 = "My name is {}, I'm {}".format("John",36) 
    # print(txt1)
    # print(txt2)
    # print(txt3)


    if parse_is_true(figpar['quiverkey']):
        v_inlet = float(yml['flower']['physics']['v_inlet'])

        # textquiver = "$\SI{{{0:.2e}}}".format(v_inlet)
        # textquiver += '{m/s}'
        # textquiver += '$'
        # textquiver = r''+textquiver
        # print(textquiver)

        qk = ax2.quiverkey(
            q,
            figpar['quiver_x'],
            figpar['quiver_y'],
            v_inlet,
            # r"$\SI{{{0:.2e}}}{{}}$".format(v_inlet,'m/s'),
            # textquiver,
            r"$\SI{{{0:.2e}}}".format(v_inlet)+'{'+figpar['quiver_unit']+'}$',
            labelpos="E",
            coordinates="figure",
        )

    if figpar['plot_levelset']:
        key_LS = 'levelset_p'
        LSdat = file[key_LS][:]
        LSdat = LSdat.transpose()
        CSlvl = ax2.contour(x_1D, y_1D, LSdat, [0.0],colors="r",linewidths=figpar['linewidth'],linestyles=figpar['linestyle'])

    ax2.set_xlabel(r"$x ( \unit{\um})$")
    ax2.set_ylabel(r"$y ( \unit{\um})$")

    ax2.set_xlim([float(x0) for x0 in figpar['xlim']])
    ax2.set_ylim([float(x0) for x0 in figpar['ylim']])
    ax2.set_aspect('equal', 'box')

    # str_nstep = str1="{:05}".format(nstep)
    str_nstep = str(nstep)
    # str_time = '{:.2e}'.format(time/plotpar['scale_time'])
    # plt.title("t "+str_time +r"$(\unit{s})$")

    plt.title('Time '+r"$\SI{{{0:.2e}}}".format(time/plotpar['scale_time'])+'{'+plotpar['unit_time']+'}$')


    plt.savefig(file_name+'_'+str_nstep+ "." + plotpar["img_format"],dpi=plotpar['dpi'])
    plt.close(fig1)


def plot_current_lines(file,
    xp,
    yp,
    mesh,
    time,
    nstep,
    plotpar,
    figpar,
    mode='close',
    fig1=None,
    ax2=None,
):
    """plot current lines
    args:

    """
    phi_array = file["i_current_x"][:].transpose()
    Eus = file["i_current_x"][:].transpose()
    Evs = file["i_current_y"][:].transpose()

    nx = mesh["nx"]
    ny = mesh["ny"]
    field_index = 1 #bulk
    data = file["phi_ele_1D"][:]

    field=veci(data,nx,ny,field_index)

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
    plt.savefig(file_name + "_" + str_nstep + "." + plotpar["img_format"],dpi=plotpar['dpi'])
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

    if 'field_index' in figpar.keys():
        field_index = figpar['field_index']
    else:
        field_index = 1 # bulk value


    if data_1D.ndim ==1:
        field0 = veci(data_1D,nx,ny,field_index)
    elif data_1D.shape[1] == 1:
        print('data_1D.shape[1] == 1')
        data_1D = veci(data_1D[:,0],nx,ny,field_index)
        field = data_1D
    elif data_1D.shape[0] == 1:
        # print('data_1D.shape[0] == 1')
        # print(data_1D)
        # print(data_1D[0]) 
        #[0,:]
        data_1D = data_1D[0]
        field0 = veci(data_1D,nx,ny,field_index)

    #TODO slice vector trans_scal

    # field0 = veci(data_1D,nx,ny,field_index)

    field = field0
    field1 = field0

    x_arr=x_1D[i0:i1+1]
    y_arr=y_1D[j0:j1+1]

    vecb_l = False
    vecb_r = False
    vecb_b = False
    vecb_t = False

    fieldtmp = np.zeros((ny + 2, nx + 2))

    # TODO distinguish u, v, w grids even though it is a dummy position to plot BC

    if parse_is_true(figpar['plot_bc']):
        if ii0 == 0:
            vecb_l=True
            i1+=1
            i0tmp2-=1 

            x_arr = np.insert(x_arr,0,x_1D[i0]-0.5*mesh['dx']/plotpar['scale_x'])

        if ii1 == nx-1:
            i1+=1
            vecb_r=True
            # x_arr.append(x_1D[end]+0.5*mesh['dx']/plotpar['scale_x'])
            np.append(x_arr,x_1D[-1]+0.5*mesh['dx']/plotpar['scale_x'])
            i1tmp2+=1 

        if jj0 == 0:
            vecb_b=True
            j1+=1
            j0tmp2-=1
            y_arr= np.insert(y_arr,0,y_1D[j0]-0.5*mesh['dy']/plotpar['scale_x'])

        if jj1 == ny-1:
            vecb_t=True
            j1+=1
            # y_arr.append(y_1D[end]+0.5*mesh['dy']/plotpar['scale_x'])
            np.append(y_arr, y_1D[-1]+0.5*mesh['dy']/plotpar['scale_x'])
            j1tmp2+=1

        # print(np.size(data_1D),nx,ny,vecb_l,vecb_r,vecb_b,vecb_t)
        # if vecb_l: 
        #     # print('len',len(fieldtmp[2:ny+1,0]), len(vecb_L(data_1D,nx,ny)))
        #     fieldtmp[2:ny+1,0] = vecb_L(data_1D,nx,ny) 
        # if vecb_r:
        #     fieldtmp[2:ny+1,-1] = vecb_R(data_1D,nx,ny) 
        # if vecb_b:
        #     fieldtmp[0,2:nx+1] = vecb_B(data_1D,nx,ny)
        # if vecb_t:
        #     fieldtmp[-1,2:nx+1] = vecb_T(data_1D,nx,ny)

        # print(vecbprint(data_1D, nx, ny))
        # testprint = vecbprint(data_1D, nx, ny)
        # print(len(testprint),2*nx+2*ny)
        # for jtest,valtest in enumerate(testprint):
        #     print(jtest,valtest)
        # print('test',ny,testprint[0],testprint[ny-1],testprint[ny],testprint[ny+1],testprint[ny+2])

        # print(vecb_L(data_1D,nx,ny))
        # print(vecb_R(data_1D,nx,ny))
        # print(vecb_B(data_1D,nx,ny))
        # print(vecb_T(data_1D,nx,ny))


        if vecb_l: 
            # print('len',len(fieldtmp[1:ny+1,0]), len(vecb_L(data_1D,nx,ny)),ny)
            # test = np.array(range(0,ny+2))
            # print(test)
            # print(test[1:ny])
            fieldtmp[1:ny+1,0] = vecb_L(data_1D,nx,ny) 
        if vecb_r:
            fieldtmp[1:ny+1,-1] = vecb_R(data_1D,nx,ny) 
        if vecb_b:
            fieldtmp[0,1:nx+1] = vecb_B(data_1D,nx,ny)
        if vecb_t:
            fieldtmp[-1,1:nx+1] = vecb_T(data_1D,nx,ny)

        #TODO end excluded recheck


    if vecb_l or vecb_r or vecb_b or vecb_t:
        fieldtmp[j0tmp:j1tmp+1,i0tmp:i1tmp+1] = field1[jj0:jj1+1,ii0:ii1+1]
        field = fieldtmp[j0tmp2:j1tmp2+1,i0tmp2:i1tmp2+1]
    else:
        field = field0[j0:j1+1,i0:i1+1]

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

    if parse_is_true(figpar['plot_grid']):
        # for igrid in i0:i1
        #     ax2.axvline(x_1D[igrid],c=lcolor,lw=lw)
        # end
        # for igrid in j0:j1
        #     ax2.axhline(y_1D[igrid],c=lcolor,lw=lw)
        # end


        if 'color_annot_bc' in figpar.keys():
            color_annot_bc = figpar['color_annot_bc']
        else:
            color_annot_bc = 'k'

        if 'color_annot_bulk' in figpar.keys():
            color_annot_bulk = figpar['color_annot_bulk']
        else:
            color_annot_bulk = 'w'



        for igrid0 in range(i0,i1+1):        
            for jgrid0 in range(j0,j1+1): 
                # print("\n",x_1D[igrid],y_1D[jgrid])
                # print("\n ",igrid," ",jgrid)
                igrid=igrid0-i0 #TODO
                jgrid = jgrid0-j0

                if ((igrid0 == i0 and vecb_l) or (igrid0 == i1 and vecb_r) or (jgrid0 == j0 and vecb_b) or (jgrid0 == j1 and vecb_t) ):
                    lcolor= color_annot_bc
                else:
                    lcolor= color_annot_bulk

                if igrid%2 == 0:
                    va="top"
                else:
                    va="bottom"

                ax2.scatter(x_arr[igrid],y_arr[jgrid],
                c=lcolor,
                s=ms,
                )

                # str1=@sprintf "%.2e" fieldtmp[jgrid0,igrid0]
                if figpar['print_mode'] == "val":
                    str1='{:.2e}'.format(field[jgrid,igrid])
                    # str1=@sprintf "%.4e %.3i %.3i %.3i" field[jgrid,igrid] igrid0 jgrid0 jgrid

                elif figpar['print_mode'] == "ij":
                    str1="{:03} {:03}".format(igrid0,jgrid0)
                else: 
                    str1="{:.2e} {:03} {:03}".format(field[jgrid,igrid],igrid0,jgrid0)                         

                ax2.annotate(str1,(x_arr[igrid],y_arr[jgrid]),fontsize=figpar['fontsize'],c=lcolor,ha="center",va=va)

    # ax2.set_title("Title")

    # str_time = '{:.2e}'.format(time/plotpar['scale_time'])
    # plt.title("t "+str_time +r"$(\unit{s})$")
    # str_time = '{:.2e}'.format(time/plotpar['scale_time'])

    plt.title('Time '+r"$\SI{{{0:.2e}}}".format(time/plotpar['scale_time'])+'{'+plotpar['unit_time']+'}$')


    ax2.set_xlabel(r"$x ( \unit{\um})$")
    ax2.set_ylabel(r"$y ( \unit{\um})$")

    # Make a colorbar for the ContourSet returned by the contourf call.
    cbar = fig1.colorbar(CS)
    cbar.ax.set_ylabel(figpar['cbarlabel'])
    # Add the contour line levels to the colorbar
    if str(figpar['isocontour']) == 'True':
        CS2 = ax2.contour(CS, 
        # levels=CS.levels[::2], 
        # levels=
        colors="r")
        cbar.add_lines(CS2)

    if figpar["plot_levelset"]:
        if "plot_case" in figpar.keys():
            if figpar["plot_case"] == "circle":
                theta1 = figpar["theta1"]
                theta2 = figpar["theta2"]

                radius = fwd.radius[i] / plotpar["scale_x"]
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
        CSlvl = ax2.contour(
            x_1D[ii0:ii1+1], y_1D[jj0:jj1+1], LSdat[jj0:jj1+1, ii0:ii1+1], [0.0], colors="r",linewidths=figpar['linewidth'],linestyles=figpar['linestyle']
        )



    # ax2.set_aspect('equal', 'box')


    if vecb_l:
        x = x_arr[i0]
        ticks_loc = ax2.get_xticks().tolist()        
        labels = [w.get_text() for w in ax2.get_xticklabels()]
        labels+=[r'$BC$']
        ticks_loc+=[x]
        ax2.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
        ax2.set_xticklabels(labels)

    if vecb_r:
        x = x_arr[i1]
        ticks_loc = ax2.get_xticks().tolist()        
        labels = [w.get_text() for w in ax2.get_xticklabels()]
        labels+=[r'$BC$']
        ticks_loc+=[x]
        ax2.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
        ax2.set_xticklabels(labels)
    if vecb_b:

        x = y_arr[j0]

        ticks_loc = ax2.get_yticks().tolist()        
        labels = [w.get_text() for w in ax2.get_yticklabels()]
        labels+=[r'$BC$']
        ticks_loc+=[x]
        ax2.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
        ax2.set_yticklabels(labels)

    if vecb_t:
        x = y_arr[j1]

        ticks_loc = ax2.get_yticks().tolist()        
        labels = [w.get_text() for w in ax2.get_yticklabels()]
        labels+=[r'$BC$']
        ticks_loc+=[x]
        ax2.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
        ax2.set_yticklabels(labels)

    # isnap = indLS
    # strtime = @sprintf "%.2e" fwd.t[isnap]*1e3
    # strrad = @sprintf "%.2e" fwd.radius[isnap]*1e6

    # strrad = '{:.2e}'.format(radius)
    # str_iter = "{:05}".format(i)

    # plt.title("t "*strtime*r"$( \unit{\ms})$"*"radius "*strrad*r"$( \unit{\um})$")

    # ax2.set_xlim([float(x0) for x0 in figpar['xlim']]) #zoom
    # ax2.set_ylim([float(x0) for x0 in figpar['ylim']])
    str_nstep = str(nstep)
    plt.savefig(file_name+'_'+str_nstep+ "." + plotpar["img_format"])

    plt.close(fig1)

def python_movie_zoom(
    h5_files,
    key,
    xp,
    yp,
    xu,
    yv,
    yml,
    mesh,
    plotpar,
    figpar,
):
    # # print(h5_files)
    # h5_files = sorted(h5_files)
    # print(h5_files)

    size_frame = len(h5_files)

    print('size_frame',size_frame)
    
    file_name = h5_files[0]

    fig1, ax2 = plt.subplots(figsize=set_size(plotpar["latex_frame_width"], fraction=float(plotpar["fig_fraction"]),
                                                ratio=1,nvary=1,ratio2=1,height=float(plotpar["latex_frame_height"])),
                                                layout="constrained")

    with h5py.File(file_name, "r") as file:

        time = file["time"][()]
        nstep = file["nstep"][()]

        fig1,ax2,CS = plot_file(
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
        figpar=figpar,
        mode='first',
        fig1=fig1,
        ax2=ax2,
        )


        # Make a colorbar for the ContourSet returned by the contourf call.
        # cbar = fig1.colorbar(CS)
        # cbar.ax.set_ylabel(r""+cbarlabel)

    def animate(i,fig1,ax2):
        # ax2.clear()

        file_name = h5_files[i]

        print(file_name)



        with h5py.File(file_name, "r") as file:

            # print(file.keys())

            time = file["time"][()]
            nstep = file["nstep"][()]           

            fig1,ax2,CS = plot_file(
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
            figpar=figpar,
            mode='film',
            fig1=fig1,
            ax2=ax2,
            )

            # if step!=0:
            #     fig1.colorbar(CS,cax=cbar.ax)

            # https://stackoverflow.com/questions/5180518/duplicated-colorbars-when-creating-an-animation

            # if (i==0): 
            #     # # Make a colorbar for the ContourSet returned by the contourf call.
            #     # cbar = fig1.colorbar(CS)
            #     # cbar.ax.set_ylabel(r""+cbarlabel)

            #     # Add the contour line levels to the colorbar
            #     if isocontour:
            #         CS2 = ax2.contour(CS, 
            #         # levels=CS.levels[::2], 
            #         # levels=
            #         colors="r")
            #         cbar.add_lines(CS2)

            # if plot_levelset:
            #     # CSlvl = ax2.contour(x_arr,y_arr, fwd.u[1,i+1,i0:i1,j0:j1], [0.0],colors="r")
            #     CSlvl = ax2.contour(x_arr,y_arr, fwd.u[1,indLS,j0:j1,i0:i1], [0.0],colors="r")
            #     # print("check levelset")
    anim = functools.partial(animate,fig1=fig1,ax2=ax2)
    ani = animation.FuncAnimation(fig1, anim, frames=size_frame, interval=size_frame, blit=False)

    ani.save(key + "." + figpar["img_format"],dpi=plotpar['dpi'])  # mp4, gif

    # with open(key+'_jshtml'+'.html', "w") as f:
    #     print(ani.to_jshtml(), file=f)

    # with open(key+'_html5'+'.html', "w") as f:
    #     print(ani.to_html5_video(), file=f)


    plt.close("all")


def python_movie_zoom_func(
    h5_files,
    key,
    xp,
    yp,
    xu,
    yv,
    yml,
    mesh,
    func,
    plotpar,
    figpar,
):
    # # print(h5_files)
    # h5_files = sorted(h5_files)
    # print(h5_files)

    size_frame = len(h5_files)

    print('size_frame',size_frame)
    
    file_name = h5_files[0]

    fig1, ax2 = plt.subplots(figsize=set_size(plotpar["latex_frame_width"], fraction=float(plotpar["fig_fraction"]),
                                                ratio=1,nvary=1,ratio2=1,height=float(plotpar["latex_frame_height"])),
                                                layout="constrained")

    with h5py.File(file_name, "r") as file:

        time = file["time"][()]
        nstep = file["nstep"][()]

        fig1,ax2,CS = func(
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
        figpar=figpar,
        mode='first',
        fig1=fig1,
        ax2=ax2,
        )


        # Make a colorbar for the ContourSet returned by the contourf call.
        # cbar = fig1.colorbar(CS)
        # cbar.ax.set_ylabel(r""+cbarlabel)

    def animate(i,fig1,ax2):
        # ax2.clear()

        file_name = h5_files[i]

        print(file_name)



        with h5py.File(file_name, "r") as file:

            # print(file.keys())

            time = file["time"][()]
            nstep = file["nstep"][()]           

            fig1,ax2,CS = func(
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
            figpar=figpar,
            mode='film',
            fig1=fig1,
            ax2=ax2,
            )

            # if step!=0:
            #     fig1.colorbar(CS,cax=cbar.ax)

            # https://stackoverflow.com/questions/5180518/duplicated-colorbars-when-creating-an-animation

            # if (i==0): 
            #     # # Make a colorbar for the ContourSet returned by the contourf call.
            #     # cbar = fig1.colorbar(CS)
            #     # cbar.ax.set_ylabel(r""+cbarlabel)

            #     # Add the contour line levels to the colorbar
            #     if isocontour:
            #         CS2 = ax2.contour(CS, 
            #         # levels=CS.levels[::2], 
            #         # levels=
            #         colors="r")
            #         cbar.add_lines(CS2)

            # if plot_levelset:
            #     # CSlvl = ax2.contour(x_arr,y_arr, fwd.u[1,i+1,i0:i1,j0:j1], [0.0],colors="r")
            #     CSlvl = ax2.contour(x_arr,y_arr, fwd.u[1,indLS,j0:j1,i0:i1], [0.0],colors="r")
            #     # print("check levelset")
    anim = functools.partial(animate,fig1=fig1,ax2=ax2)
    ani = animation.FuncAnimation(fig1, anim, frames=size_frame, interval=size_frame, blit=False)

    ani.save(key + "." + figpar["img_format"],dpi=plotpar['dpi'])  # mp4, gif

    # with open(key+'_jshtml'+'.html', "w") as f:
    #     print(ani.to_jshtml(), file=f)

    # with open(key+'_html5'+'.html', "w") as f:
    #     print(ani.to_html5_video(), file=f)


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

    # p1, = ax.plot(varx, ones(ny),colors[1], label=label1,ls=ls)
    # p2, = twin1.plot(varx, ones(ny), colors[2], label=label2,ls=ls2)
    # p3, = twin2.plot(varx, ones(ny), colors[3], label=label3,ls=ls3)

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

        str_time = '{:.2e}'.format(time/plotpar['scale_time'])  #TODO scale time
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
        str_time = '{:.2e}'.format(time/plotpar['scale_time'])  #TODO scale time
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


def veci(data,nx,ny,field_index):
    """Returns ith field stored in the 1D vector like in Flower.jl code
    
    args:
        field_index (int): index of bulk or interface data
    return: 
        bulk (i=1) or i-th interface field
    """
    # print(data.shape)
    field = data[(field_index-1)*ny*nx:field_index*ny*nx] 
    # print(field.shape)

    field = np.reshape(field, (nx, ny))            
    field = field.transpose()
    return field


def vecb(data, nx, ny):
    """BC at border"""
    extract = data[-2 * ny - 2 * nx:]
    # print('vecb',extract,'len',len(extract))
    return extract

def vecbprint(data, nx, ny):
    """BC at border"""
    extract = data[-2 * ny - 2 * nx:]
    # print('vecb',extract,'len',len(extract))
    return extract

def vecb_L(data, nx, ny):
    """BC at left border"""
    data = vecb(data, nx, ny)
    extract = data[0 : ny]
    # print('vecb_L',extract,'len',len(extract))
    return extract


def vecb_B(data, nx, ny):
    """BC at bottom border"""
    data = vecb(data, nx, ny)
    extract = data[ny : nx + ny]
    # print('vecb_B',extract,'len',len(extract))
    return extract


def vecb_R(data, nx, ny):
    """BC at right border"""
    data = vecb(data, nx, ny)
    extract = data[nx + ny : 2 * ny + nx]
    # print('vecb_R',extract,'len',len(extract))
    return extract



def vecb_T(data, nx, ny):
    """# BC at top border"""
    data = vecb(data, nx, ny)
    extract = data[2 * ny + nx : 2 * nx + 2 * ny]
    # print('vecb_T',extract,'len',len(extract))
    return extract



# vecb_L(a,g::G) where {G<:Grid} = @view vecb(a, g)[1:g.ny]
# vecb_B(a,g::G) where {G<:Grid} = @view vecb(a, g)[g.ny+1:g.ny+g.nx]
# vecb_R(a,g::G) where {G<:Grid} = @view vecb(a, g)[g.ny+g.nx+1:2*g.ny+g.nx]
# vecb_T(a,g::G) where {G<:Grid} = @view vecb(a, g)[2*g.ny+g.nx+1:2*g.ny+2*g.nx]
