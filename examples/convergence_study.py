import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import yaml
import sys
import h5py
import math
from termcolor import colored

import matplotlib.patches as patches

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.transforms as transforms
from matplotlib.ticker import MaxNLocator


# import numpy as np
from scipy.stats import pearsonr

# module to plot files from Flower.jl
# from plot_flower import * 
from plot_flower import set_size, init_fig, compute_slope, roundlog, \
   logticks,reshape_data,veci,vecb_L,reshape_data_veci,plot_current_lines,\
   plot_python_pdf_full2,plot_file,plot_schematics,plot_schematics_full,add_schematics,compute_slope_lin_or_log,plot_vector 

plt.rcParams["text.parse_math"] = False #necessary for mhchem


plt.rcParams['text.usetex'] = True
plt.rcParams["font.size"] = "11"
# plt.rcParams["font.size"] = "8"

# plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

plt.rc('text.latex', preamble="\n".join([ # plots will use this preamble
        r'\usepackage{amsmath}',
        r'\usepackage{booktabs}',
        r"\usepackage{siunitx}",
        r"\setlength{\abovedisplayskip}{0pt}",
        r" \setlength{\belowdisplayskip}{0pt}",
        r"\setlength{\belowdisplayshortskip}{0pt}",
        r"\setlength{\abovedisplayshortskip}{0pt}"
        r"\addtolength{\jot}{-4pt}"
        r"\usepackage{mhchem}",
       ])
)

# Colors

OkabeIto=["#E69F00", #0 orange clair 230, 159, 0
"#56B4E9", #1 bleu clair 86, 180, 233
"#009E73", #2 vert 0, 158, 115
"#F0E442", #3 jaune 240, 228, 66
"#0072B2", #4 bleu 0, 114, 178
"#D55E00", #5 orange 213, 94, 0
"#CC79A7", #6 rose 204, 121, 171
"#000000"] #7 noir 0 0 0

# colors=OkabeIto

# colors=["#000000" for color in OkabeIto]

# colors[1]="#000000"
# colors[1]=OkabeIto[5] #bleu
# colors[2]=OkabeIto[6] #orange



# "#000000", for ref

colors= [
   "#0072B2", #bleu
   "#D55E00", #orange
   "#009E73", #vert,
   "#CC79A7",
]

orange_Okabe = "#D55E00"
blue_Okabe = "#0072B2"




def create_camembert_markers(markers,num_slices,startangle0 = np.pi/2):

   size=60

   numfill=50
   # numfill = 100

   startangle=np.pi/2

   # camembert_linewidth=0.025
   camembert_linewidth=0


   # linewidth = 0.0

   linewidth=0.5

   linewidth2=linewidth/2

   sizes = []

   for i in range(num_slices):
      # first define the ratios
      # r1 = i/(num_slices-1)
      r1 = 1/(num_slices)

      # startangle= startangle0[i] -2 * np.pi*r1
      startangle= startangle0 -2 * np.pi * i /(num_slices)

      # calculate the points of the first pie marker
      # these are just the origin (0, 0) + some (cos, sin) points on a circle
      # x1 = np.cos(startangle-2 * np.pi * np.linspace(0, r1,num=numfill))
      # y1 = np.sin(startangle-2 * np.pi * np.linspace(0, r1,num=numfill))

      x1 = np.cos(startangle+2 * np.pi * np.linspace(0, r1,num=numfill))
      y1 = np.sin(startangle+2 * np.pi * np.linspace(0, r1,num=numfill))

      xy1 = np.row_stack([[0, 0], np.column_stack([x1, y1])])
      # s1 = np.abs(xy1).max()
      markers[i]=xy1

      s1 = np.abs(xy1).max()
      sizetemp=s1**2 * size
      sizetemp=np.sqrt(sizetemp)
      sizetemp=sizetemp-linewidth2-camembert_linewidth
      sizetemp=sizetemp**2
      sizes.append(sizetemp)

   
   markers=markers[::-1]
   sizes=sizes[::-1]

   return markers,sizes


def plot_errors(domain_length,nx_list,l1_rel_error,l2_rel_error,linfty_rel_error):

   fig, ax = plt.subplots(layout="constrained")

   x = nx_list

   print(x)
   print(l1_rel_error)
   print(l2_rel_error)
   print(linfty_rel_error)


   plt.plot(x, l1_rel_error,label = r'$l_1$')
   plt.plot(x, l2_rel_error,label = r'$l_1$')
   plt.plot(x, linfty_rel_error,label = r'$l_1$')

   ax.scatter(x, l1_rel_error,label = r'$l_1$')
   ax.scatter(x, l2_rel_error,label = r'$l_1$')
   ax.scatter(x, linfty_rel_error,label = r'$l_1$')

   ax.set_xscale("log", base=10)
   ax.set_yscale("log", base=10)

   plt.legend()

   plt.savefig("errors.pdf",transparent=True)
   print("plot python")

   plt.close(fig)



def plot_errors_from_pandas(df,figpar,plotpar,colors,filename):
   """
   Plot radius pandas DF, with slope
   """

   fig1,ax2 = init_fig(plotpar,figpar)

   # print("t",fwd.t)
   # print("radius",fwd.radius.*1.e6)
   # print("current_i", current_i)
   # print("radius ",fwd.radius[1:current_i+1])
   # print("\nradius ",fwd.radius)

   color="#4d5156"

   df['1/n'] = 1/df["nx_list"]
   # plt.plot(df["1/n"],df["l1_rel_error"],color=color)

   alpha=0.75

   markers = ['<','>','^']

   markers, sizes  = create_camembert_markers(markers,3,np.pi/2)



   xls=df["1/n"].to_numpy()



   print(df)

   plot_slope = True
   # plot_slope = False

   if plot_slope:
 
      xlsindex=df["nx_list"].to_numpy()

      tstart=float(figpar['slope_start'])
      tstop=float(figpar['slope_stop'])


      # tstop=1/float(figpar['slope_start'])
      # tstart=1/float(figpar['slope_stop'])

      print('tstart',tstart,tstop)
      istart = 0
      istop = 0

      print("first")
      for i in range(len(xlsindex)):
         print(i,xlsindex[i],tstart,tstop)
         # if xlsindex[i]<=tstart:
         #    break
         # elif xlsindex[i]>=tstart:
         #    # print('larger',xlsindex[i],istart)
         #    istart=i
         #    break
         # istart=i
         if xlsindex[i]>=tstart:
            # print('larger',xlsindex[i],istart)
            istart=i
            break
      # print('larger',istart)

      print("second")

      for i in range(len(xlsindex)):
         print(i,xlsindex[i],tstart,tstop)
         if xlsindex[i]>tstop:
            istop=i
            break
         istop=i
      print('istart',istart,istop)

      print(xls)

      x="t"
      y="r"

      param_line=[]
      slope_and_correlation=0
      R2=0
      color_line="k"
      # colors="#fa8b2b"
      alpha= 1 
      # \definecolor{mdlsorange}{HTML}{fa8b2b}
      # % \definecolor{mdlsgray}{HTML}{8e98a4}
      # \definecolor{mdlsgray}{HTML}{4d5156}
      print('xls',xls)
      xls=xls[istart:istop+1]

      print('xls',xls)

      for i,err in enumerate(figpar['var']):
         yls=df[err].to_numpy()
         yls=yls[istart:istop+1]

         print('yls',yls)

         slope_and_correlation=[0,0]
         print('test slope_and_correlation ',slope_and_correlation)
         compute_slope(ax2,xls,yls,
                     #   x,y,
                       slope_and_correlation,R2,param_line,color_line,alpha,plot_text=False)
         print(slope_and_correlation)
         if 'l1' in err:
            label = r'$l_1$'
         elif 'l2' in err:
            label = r'$l_2$'
         elif 'linfty' in err:
            label = r'$l_\infty$'

         # text_slope=', slope={:.2g}, R²={:.2g}'.format(float(slope_and_correlation[0]),float(slope_and_correlation[1]))

         # text_slope = 'Time '+r"$\SI[retain-zero-exponent=true]{{{0:.2e}}}".format(time/plotpar['scale_time'])+'{'+plotpar['unit_time']+'}$'
         # text_slope = ', slope '+r"$\SI[retain-zero-exponent=true]{{{0:.2e}}}, R^2 {1:.2e}".format(slope_and_correlation[0],slope_and_correlation[1]) + '$'
         # text_slope = ', slope '+r"${0:.2f}, R^2 {1:.2f}".format(slope_and_correlation[0],slope_and_correlation[1]) + '$'
         # text_slope = r"$\mathrm{{, slope}} {0:.2f}, R^2 {1:.2f}".format(slope_and_correlation[0],slope_and_correlation[1]) + '$'
         # text_slope = r"$\mathrm{{{}}} {1:.2f}, R^2 {2:.2f}".format(', slope',slope_and_correlation[0],slope_and_correlation[1]) + '$'
         if math.isnan(slope_and_correlation[1]): 
            text_slope = r"$\mathrm{{{text}}}: {slope:.2f}".format(text=', slope',slope=slope_and_correlation[0]) + '$'
         else:
            text_slope = r"$\mathrm{{{text}}}: {slope:.2f}, R^2: {R2:.2f}".format(text=', slope',slope=slope_and_correlation[0],R2=slope_and_correlation[1]) + '$'
         # text_slope = r"$ text {slope:.2f}, R^2 {R2:.2f}".format(slope=slope_and_correlation[0],R2=slope_and_correlation[1]) + '$'

         print('text_slope',text_slope)

         ax2.scatter(df["1/n"],df[err],
                     color=colors[i],
                     alpha=alpha,
                     marker=markers[i],
                     s = sizes[i],

                     linewidth=0,
                     # linewidth=camembert_linewidth, #edge pie 
                     # edgecolors='white',
                     clip_on = True,

                     label = label+text_slope,

                     zorder=5,
                     )
         

         # print(colors)
   

   
   # ax2.scatter(df["1/n"],df["l1_rel_error"],color=color,alpha=alpha,marker='<',label = r'$l_1$')

   # ax2.scatter(df["1/n"],df["l2_rel_error"],color=colors[2],alpha=alpha,marker='>',label = r'$l_2$')

   # ax2.scatter(df["1/n"],df["linfty_rel_error"],color=colors[3],alpha=alpha,marker='^',label = r'$l_\infty$')



   # ilog0 = size_frame ÷ 2 
   # ilog1 = 
   # xlog = 
   # print("\n")

   # ax2.set_title("Title")
   ax2.set_xlabel('1/Number of points per direction')
   # ax2.set_ylabel(L"$R (m)$")
   # ax2.set_ylabel(r"$l_1 relative error l_1=\frac{\sum{S_i \lvert p_i- p_{i}^e \rvert } }{\sum{S_i\lvert p_{i}^e \rvert}}$")
   # ax2.set_ylabel(r"$l_1 relative error$")
   # ax2.set_ylabel(r"$l_1=\frac{\sum{S_i \lvert p_i- p_{i}^e \rvert } }{\sum{S_i\lvert p_{i}^e \rvert}}$")
   # ax2.set_ylabel(r"$\sum{S_i \lvert p_i- p_{i}^e \rvert}/\sum{S_i\lvert p_{i}^e \rvert}$")
   ax2.set_ylabel('Relative error')




   ax2.set_xscale("log")
   ax2.set_yscale("log")




   plt.legend()
   # plt.axis("equal")

   # prefix="./"

   prefix = filename.replace(".h5", "") +'_'



   if 'xlim' in figpar:
      # ax2.set_xlim(figpar['xlim'])
      xlim=np.zeros(2)

      xlim[0] = figpar['xlim'][0]
      xlim[1] = figpar['xlim'][1]
      print('xlim',xlim)
      # roundlog(xlim)
      ax2.set_xlim(xlim)
      # ax2.set_xlim([(10**-1,10**0) )])

      # plt.xlim( [10**-3,10**-1] )
      
      # plt.xlim( [10**-5,10**-1] )


      # xticks,xlabels,xminorticks,xminorticklabels=logticks(xlim)
      
   
      # ax2.set_clip_on(False)
      # xminorticklabels='' #deactivate

      # try:
      #    # plt.savefig(filename+'bug'+str(bug)+'.'+ext)
      #    # if plotscatter:
      #    minorlabels=ax2.set_xticks(ticks=xminorticks,
      #                labels=xminorticklabels,
      #                minor=True,
      #                #   rotation=axrotation,
      #                #   rotation_mode="anchor",
      #                #   ha=axpos,
      #                )
      # except:
      #    print('Error ticks')

   # plt.xlim( [10**-5,10**-1] )


   plt.savefig(prefix+figpar['file']+".pdf",transparent=True)
   plt.savefig(prefix+figpar['file']+".svg",transparent=True)

   plt.close(fig1)

   print(df.to_latex(index=False,

                  formatters={"name": str.upper},

                  float_format="{:.2e}".format,))
   
   print(df.to_html(index=False,

               formatters={"name": str.upper},

               float_format="{:.2e}".format,))

def plot_errors_from_h5():
   """
   Plot errors from h5 files
   """

   # print('arg', len(sys.argv),sys.argv)
   if len(sys.argv) == 2:
      # List all files in the current directory
      all_files = os.listdir(".")
      h5_files = [file for file in all_files if file.endswith(".h5")]
   else:
      h5_files = sys.argv[2::]

   try:
      yamlfile = sys.argv[1]
      if ".yml" not in yamlfile:
         yamlfile += ".yml"
   except Exception as error:
      print(error)
      print(colored("error", "red"))

   with open(yamlfile, "r") as file:
      yml = yaml.safe_load(file)

      plotpar = yml["plot"]

      # print(h5_files)
      h5_files = sorted(h5_files)
      print('\n sorted \n')
      print(h5_files)
      
      nsteps = len(h5_files)

      for figpar in plotpar['curves']:

         print(colored(figpar['file'],'cyan'))
         # print(figpar)

         if 'radius' in figpar['var']: #we do not plot the figures with radius
            continue
         
         if 'func' not in figpar.keys():
            continue
         if 'plot_errors_from_h5' != figpar['func']:
            continue

         time_list =[]
         radius_list=[]
         for file_name in h5_files:
               with h5py.File(file_name, "r") as file:
                  print(file.keys())
                  nx_list = file["nx_list"][()]
                  l1_rel_error = file["l1_rel_error"][()]
                  l2_rel_error = file["l2_rel_error"][()]
                  linfty_rel_error = file["linfty_rel_error"][()]
                  
            

                  print(nx_list,l1_rel_error)
                  # time_list.append(time)
                  # radius_list.append(radius)

                  # df = pd.DataFrame({'nx_list': nx_list, 'l1_rel_error': l1_rel_error})
                  # df['l2_rel_error'] = l2_rel_error
                  # df['linfty_rel_error'] = linfty_rel_error

                  df = pd.DataFrame({'nx_list': nx_list})


                  for i,err in enumerate(figpar['var']):
                     print(colored(figpar['var'],'red'))
                     df[err] = file[err][()]

                  print(df)

                  plot_errors_from_pandas(df,figpar,plotpar,colors,file_name)



def plot_schematics_func():
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

   plotpar["scale_x"] = float(plotpar["scale_x"])
   plotpar["scale_y"] = float(plotpar["scale_y"])

   scale_x = float(plotpar["scale_x"])
   scale_y = float(plotpar["scale_y"])

   xp /= scale_x
   yp /= scale_y
   xu /= scale_x
   yv /= scale_y



    # Schematics describing BC
   for figpar in plotpar['schematics']:


      if 'func' in figpar.keys():
         func = globals()[figpar['func']] #'plot_current_lines'
         func(figpar,plotpar)
      else:
         plot_schematics(figpar,plotpar)



def plot_convergence_study_func():
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

   plotpar["scale_x"] = float(plotpar["scale_x"])
   plotpar["scale_y"] = float(plotpar["scale_y"])

   scale_x = float(plotpar["scale_x"])
   scale_y = float(plotpar["scale_y"])

   xp /= scale_x
   yp /= scale_y
   xu /= scale_x
   yv /= scale_y



  





   for figpar in plotpar["curves"]:
      
      try:
         # print(figpar)
         if 'func' not in figpar.keys():
            continue
            
         # print(figpar)

         if 'func' in figpar.keys():
            func = globals()[figpar['func']] #'plot_current_lines'
         else:
            func = globals()['plot_file']

         if (not isinstance(figpar['var'], str)) and len(figpar['var'])>0:
            key = figpar['var'][0]
         else:
            key = figpar['var']

         # print('key',key)

         # print(figpar)

         print(colored(figpar['file'], "cyan"))


         
         plot_convergence_func(
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
      except:
         print(colored('Failed '+figpar['file'], "red"))   
         raise # was: pass


   for figpar in plotpar['figures']:
      
      try:
         # print(figpar)
         if 'func' not in figpar.keys():
            continue
            
         # print(figpar)

         if 'func' in figpar.keys():
            func = globals()[figpar['func']] #'plot_current_lines'
         else:
            func = globals()['plot_file']

         if (not isinstance(figpar['var'], str)) and len(figpar['var'])>0:
            key = figpar['var'][0]
         else:
            key = figpar['var']

         # print('key',key)

         # print(figpar)

         print(colored(figpar['file'], "cyan"))


         
         plot_convergence_func_new_ax(
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
      except:
         print(colored('Failed '+figpar['file'], "red"))
         raise # was: pass

 


def plot_convergence_func(
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

   print('size_frame',size_frame,key)

   # file_name = h5_files[0]

   fig1,ax2 = init_fig(plotpar,figpar)

   if 'logplot' in figpar.keys():
      figpar['error_list'] = []
      figpar['x_list'] = []
      # print()

   if 'streamplot_color' in figpar.keys():
      # cbar=[None,None]
      cbar=[]

   else:
      cbar=None

   cbar=[]

   # Create dataframe for latex and html
   columns = ['Iteration', 
            #   '$\phi$', 
              'Residual',
              'Variation']
   df = pd.DataFrame(columns=columns)
   df['Iteration'].astype("Int64")
   

   for i,file_name in enumerate(h5_files):
      yml['study']['iter'] = i
      with h5py.File(file_name, "r") as file:

         try:
            time = file["time"][()]
            nstep = file["nstep"][()]
         except:
            time = 0
            nstep = 0 
            print("time not available")

         print(key,file_name)


         try:
            # Fill dataframe for latex and html
            k = file['poisson_iter'][()]

            # phi_wall = file['variation_electrical_potential'][()]

            residual = file['residual_electrical_potential'][()]

            variation = file['variation_electrical_potential'][()]

            # Create the list of values
            data_list = [k, 
                        #  phi_wall,
                        residual, 
                        variation]

            # Append the list to the DataFrame
            df.loc[len(df)] = data_list
         except:
            print('no poisson_iter')


         nx = file['nx'][()]
         ny = nx 

         # if key=="u_1D":
         #    nx=nx+1
         #    # x_1D = xu 
         #    # y_1D = yp
         #    # key_LS = "levelset_u"

         # elif key=="v_1D":
         #    ny=ny+1
            # x_1D = xp
            # y_1D = yv 
            # key_LS = "levelset_v"

         # else:
         #    x_1D = xp
         #    y_1D = yp
         #    key_LS = "levelset_p"

         # print('nx',nx)

         mesh["nx"] = nx
         mesh["ny"] = ny

         print(colored('mesh '+str(nx)+" "+str(mesh["nx"]),'red'))

         mesh["xmax"] = float(mesh["xmax"])
         mesh["xmin"] = float(mesh["xmin"])

         mesh["ymax"] = float(mesh["ymax"])
         mesh["ymin"] = float(mesh["ymin"])

         mesh["dx"] = (mesh["xmax"] - mesh["xmin"]) / mesh["nx"]
         mesh["dy"] = (mesh["ymax"] - mesh["ymin"]) / mesh["ny"]

         # import numpy as np

         xp = np.linspace(float(mesh["xmin"]), float(mesh["xmax"]), int(mesh["nx"]))
         yp = np.linspace(float(mesh["ymin"]), float(mesh["ymax"]), int(mesh["ny"]))

         print('TODO mesh check')

         dx = (float(mesh["xmax"]) - float(mesh["xmin"])) / int(mesh["nx"])
         dy = (float(mesh["ymax"]) - float(mesh["ymin"])) / int(mesh["ny"])


         # print(xp)
         # print(yp)

         print('len(xp)',len(xp),len(yp))

         xp[0] = dx/2
         yp[0] = dy/2

         for i in range(1,len(xp)):
            xp[i] = xp[i-1]+dx

         for i in range(1,len(yp)):
            yp[i] = yp[i-1]+dx

         # print('xp',xp)
         # print('yp',yp)
         
         print('len(xp)',len(xp),len(yp))

         xu = xp + dx / 2
         xu = np.insert(xu, 0, xp[0] - dx / 2)
         # yu = yp # xv = xp
         yv = yp + dy / 2
         yv = np.insert(yv, 0, yp[0] - dy / 2)

         plotpar["scale_x"] = float(plotpar["scale_x"])
         plotpar["scale_y"] = float(plotpar["scale_y"])

         scale_x = float(plotpar["scale_x"])
         scale_y = float(plotpar["scale_y"])

         xp /= scale_x
         yp /= scale_y
         xu /= scale_x
         yv /= scale_y

         if i == 0:
            mode = 'first'
         else:
            mode = 'next'

         fig1,ax2,cbar = func(
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
         mode=mode,
         fig1=fig1,
         ax2=ax2,
         cbar=cbar,
         )
       
   file_name = figpar['file']

   # plt.savefig(file_name+ "." + plotpar["img_format"],dpi=plotpar['dpi']) #also for film for latex display
   # plt.savefig(file_name+ ".svg",dpi=plotpar['dpi']) #also for film for latex display

   if 'fontsize' in figpar.keys():
      font_size = figpar['fontsize']
   else:
      font_size = plotpar['fontsize']

   if 'add_schematics' in figpar.keys():
      ax2 = add_schematics_curve(fig1,ax2,font_size,figpar)


   if 'logplot' in figpar.keys():
     
      xls = figpar['x_list']
      yls = figpar['error_list']

      q_values = []
      n = len(yls)

      for k in range(2, n-1):
     
         # Calculate the numerator and denominator
         numerator = (yls[k+1] - yls[k]) / (yls[k] - yls[k-1])
         denominator = (yls[k] - yls[k-1]) / (yls[k-1] - yls[k-2])

         # Compute q
         q = np.log(abs(numerator)) / np.log(abs(denominator))
         print(k,'order',q)
         q_values.append(q)
      
      print(colored(q_values,'red'))

      xls = xls[figpar['slope_start']:figpar['slope_stop']]
      yls = yls[figpar['slope_start']:figpar['slope_stop']]

      print(xls)
      print(yls)

      for j,yelem in enumerate(yls):
         print(yelem**2)

      # if plot_slope:
      #    slope_and_correlation = [0,0]
      #    R2 = 0
      #    color_line = 'k'  
      #    alpha = 1
      #    param_line=[]

      #    compute_slope_lin_or_log(ax2,xls,yls,
      #                #   x,y,
      #                slope_and_correlation,R2,param_line,color_line,alpha,plot_text=False)
      
     


      
      ax2.xaxis.set_major_locator(MaxNLocator(integer=True))
      
      # plt.legend()

      if 'legend' in figpar.keys():
         if not figpar['legend']:
            # print()
         # else:
            ax2.get_legend().remove()
      # print()


   if 'macro_file_name' in figpar.keys():
      # print(figpar['macro_file_name'])
      # plt.savefig(eval(figpar['macro_file_name']),dpi=plotpar['dpi'])

      for macro in figpar['macro_file_name']:
         # print(macro)
         plt.savefig(eval(macro),dpi=plotpar['dpi'],transparent=True)
   
   else:
      plt.savefig(file_name+ "." + plotpar["img_format"],dpi=plotpar['dpi'],transparent=True) #also for film for latex display
         



   plt.close("all")

   print_latex =True

   try:
      if print_latex:

         print()
         print(colored('Tables','cyan'))

         # Convert the DataFrame to a LaTeX table with scientific notation and three significant digits
         latex_table = df.to_latex(index=False,column_format='ccc',formatters={
         'Iteration': lambda x: f'{x:.0f}',
         'Residual': lambda x: f'{x:.3e}',
         'Variation': lambda x: f'{x:.3e}'
         })

         # Print the LaTeX table
         print(latex_table)


         print(df.to_html(index=False,
            formatters={
            'Iteration': lambda x: f'{x:.0f}',
            'Residual': lambda x: f'{x:.3e}',
            'Variation': lambda x: f'{x:.3e}'
            }))
         
         print()
   except:
      print('cannot print table poisson_iter')


def add_schematics_curve(fig,ax2,fontsize,figpar):

   # x1, x2, y1, y2 = 0,1, 55,57  # subregion of the original image
   x1, x2, y1, y2 = figpar['add_schematics_coords']


   #  inset_ax = ax2.inset_axes(
   #  [-1, 0.125, figpar['schematics_width'], figpar['schematics_height']],

   #  # [0.5, 0.5, 0.47, 0.47],
   #  xlim=(x1, x2), ylim=(y1, y2), xticklabels=[], yticklabels=[])


   inset_ax = inset_axes(ax2,width=figpar['schematics_width'], height=figpar['schematics_height'], loc=figpar['schematics_loc'])

   # axins.imshow(Z2, extent=extent, origin="lower")

   lw_inset = 0.5
   color_inset = orange_Okabe


   #  rect, lines = ax2.indicate_inset_zoom(inset_ax, edgecolor=color_inset,
   #                                        lw=lw_inset,
   #                                        alpha=1)

   #  # rect.set_edgecolor('none')
   #  # lines.set_width(0.125)
   #  plt.setp(lines, linewidth=lw_inset,color=color_inset)



   xmin = -0.1

   inset_ax.set_xlim(xmin, 1)
   inset_ax.set_ylim(0, 1)
   # inset_ax.axis('off')

   inset_ax.get_xaxis().set_visible(False)
   inset_ax.get_yaxis().set_visible(False)

   for key, spine in inset_ax.spines.items():
      spine.set_edgecolor(orange_Okabe)
      spine.set_linewidth(lw_inset)


   # Draw the gradient in the inset
   # inset_ax.fill_between([0, 1], [0, 0], [1, 0.5], color='cyan', alpha=0.3)


   gray = 0.3
   # Draw the electrode
   electrode = patches.Rectangle((xmin, 0), -xmin, 1.0, 
                              #   edgecolor='black', 
                              facecolor= (gray, gray, gray) #'gray'
                              )
   inset_ax.add_patch(electrode)

   inset_ax.text(xmin/2, 0.5, 'Electrode', color='w',fontsize=fontsize,rotation=90, va='center',ha='center')



   liq_height=0.5

   liq_height_2 = 0.75

   liq_height = 1.0
   liq_height_2 =1.0
   # inset_ax.fill_between([0, 1-liq_height], [1, 1], [1, 1-liq_height_2], color='cyan', alpha=0.3)

   #region Plot liquid
   inset_ax.fill_between([0, 1],[1,1] ,[1-liq_height, 1-liq_height_2])

   #  inset_ax.plot([0, 1], [1-liq_height, 1-liq_height_2], ls='-',color='r')
   #endregion


   if 'macro_plot_BC' in figpar.keys():
      exec(figpar['macro_plot_BC'],
         )
   else:
      #phi
      inset_ax.text(             
      0,# xmin/2
      0.65, r'$\frac{\partial \phi }{\partial n} = \frac{i}{\kappa}$', fontsize=fontsize,color='w',ha='left',va='center')

      inset_ax.text(0.65, 0.0, r'$\frac{\partial \phi }{\partial n} = 0$', fontsize=fontsize,va='bottom',ha='center',color='w')
      inset_ax.text(1.0, 0.5, r'$\phi = 0$', fontsize=fontsize,va='bottom',ha='right',color='w')


 


   # inset_ax.text(0.1, 0.35, r'$c_{H_2} = c_{H_2, 0}$', fontsize=fontsize)
   # inset_ax.text(xmin/2, 0.5, 'Electrode', color='w',fontsize=12,rotation=90, va='center',ha='center')
   # inset_ax.text(xmin/2, 0.75, r'$\frac{\partial c_{H_2}}{\partial n} = -\frac{i}{2FD}$', fontsize=fontsize,color='w')

   if 'macro_show_slice' in figpar.keys():
      exec(figpar['macro_show_slice'],
         )
   # print(figpar['macro_show_slice'])

   return ax2

def plot_convergence_func_new_ax(
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

   print('size_frame',size_frame,key)

   # file_name = h5_files[0]

   fig1,ax2 = init_fig(plotpar,figpar)

   if 'streamplot_color' in figpar.keys():
      # cbar=[None,None]
      cbar=[]

   else:
      cbar=None

   cbar=[]

   for i,file_name in enumerate(h5_files):
      yml['study']['iter'] = i
      with h5py.File(file_name, "r") as file:

         try:
            time = file["time"][()]
            nstep = file["nstep"][()]
         except:
            time = 0
            nstep = 0 
            print("time not available")



         print(colored(key+' '+file_name,'cyan'))


         nx = file['nx'][()]
         ny = nx 

         # if key=="u_1D":
         #    nx=nx+1
         #    # x_1D = xu 
         #    # y_1D = yp
         #    # key_LS = "levelset_u"

         # elif key=="v_1D":
         #    ny=ny+1
            # x_1D = xp
            # y_1D = yv 
            # key_LS = "levelset_v"

         # else:
         #    x_1D = xp
         #    y_1D = yp
         #    key_LS = "levelset_p"

         # print('nx',nx)

         mesh["nx"] = nx
         mesh["ny"] = ny

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


         print(xp)
         print(yp)

         print(len(xp),len(yp))

         xp[0] = dx/2
         yp[0] = dy/2

         for i in range(1,len(xp)):
            xp[i] = xp[i-1]+dx

         for i in range(1,len(yp)):
            yp[i] = yp[i-1]+dx

         print(xp)
         print(yp)
         
         print(len(xp),len(yp))



         print('TODO mesh check')

         xu = xp + dx / 2
         xu = np.insert(xu, 0, xp[0] - dx / 2)
         # yu = yp # xv = xp
         yv = yp + dy / 2
         yv = np.insert(yv, 0, yp[0] - dy / 2)

         plotpar["scale_x"] = float(plotpar["scale_x"])
         plotpar["scale_y"] = float(plotpar["scale_y"])

         scale_x = float(plotpar["scale_x"])
         scale_y = float(plotpar["scale_y"])

         xp /= scale_x
         yp /= scale_y
         xu /= scale_x
         yv /= scale_y

         # if i == 0:
         #    mode = 'first'
         # else:
         #    mode = 'next'
         
         # mode = ' '

         # fig1,ax2,cbar = func(
         # file,
         # key,
         # xp,
         # yp,
         # xu,
         # yv,
         # yml,
         # mesh,
         # time,
         # nstep,
         # plotpar,
         # figpar=figpar,
         # mode=mode,
         # fig1=fig1,
         # ax2=ax2,
         # cbar=cbar,
         # )

         func(
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
            mode='close',
            fig1=None,
            ax2=None,
            cbar=None,
            )

         # plt.close("all")
       
   # file_name = figpar['file']

   # plt.savefig(file_name+ "." + plotpar["img_format"],dpi=plotpar['dpi']) #also for film for latex display

   # plt.close("all")


def plot_1D(
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
    cbar=None,
):
   """Plot one figure for field"""

   file_name = figpar['file']

   nx = file['nx'][()]
   ny = nx 
   from matplotlib.colors import LinearSegmentedColormap, ListedColormap
   cmap = ListedColormap(colors)

   # print('nx',nx,ny)

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

   print(file.keys())

   field_index = 1 #bulk
   file_name = figpar['file']

   if 'field_index' in figpar.keys():
      field_index = figpar['field_index']
   else:
      field_index = 1 # bulk value
      
   # print(key,nstep,time,"max ",np.max(data),'min',np.min(data))

   #  print('key',key)



   if figpar == None:
      figpar = plotpar

   if mode == 'first':
      ax20 = ax2
      # print(ax20)
      # ax20.cla()
      if len(figpar['var'])>1:
         twin1 = ax20.twinx()
         twin2 = ax20.twinx()
   else:
      ax20 = ax2

   # elif mode == 'film':
   #    # print(ax2)
   #    ax20,twin1,twin2 = ax2
   #    # ax20.cla()     
   #    if len(figpar['var'])>1:   
   #       twin1.cla()
   #       twin2.cla()
   #    # twin1 = ax20.twinx()
   #    # twin2 = ax20.twinx()

   # else:
   #    fig1,ax20 = init_fig(plotpar,figpar)
   #    if len(figpar['var'])>1:
   #       twin1 = ax20.twinx()
   #       twin2 = ax20.twinx()



   ###########################################



   if 'plot_mode' not in figpar.keys():
      figpar['plot_mode'] = plotpar['plot_mode']

   scale_time = float(plotpar["scale_time"])
   scale_x = float(plotpar["scale_x"])
   cmap = plt.get_cmap(plotpar["cmap"])

   # cbarlabel = plotpar["cbarlabel"]
   isocontour = plotpar["isocontour"]

   time /= scale_time 
   # radius /= scale_x

   # fig.subplots_adjust(right=0.75)

  

   for i, varxy in enumerate(figpar['var']):

      varx = varxy[0]

      if varx == 'x_1D':
         varx = x_1D
         print(x_1D)
      elif varx == 'y_1D':
         varx = y_1D
      elif varx == 'poisson_iter':
         varx = file['poisson_iter'][()] 
      else:
         varx = x_1D



      vary = varxy[1]

      label_i = r""+figpar['labels'][i][1]

      label1 = str(nx)

      ls  = eval(figpar['linestyles'][i])

      lw = figpar['linewidth']

      labelx = figpar['labels'][i][0]
  
      if len(figpar['var'])>1:
         # Offset the right spine of twin2.  The ticks and label have already been
         # placed on the right by twinx above.
         twin2.spines.right.set_position(("axes", figpar['axis_offset']))

      data = file[vary][:]

      # data= reshape_data_veci(data,nx,ny,field_index)
      # slice_1D = data[0,:]

      slice_1D = eval(figpar['macro_slice'])

      print('data',slice_1D)
      try:
         print('len data',len(slice_1D))
      except:
         print('one point')



      try:
         tick0 = list(eval(figpar['ticks'][0]))
      except:
         print('no ticks')

      if 'logplot_x' in figpar.keys():
         if figpar['logplot_x']:
            ax2.set_xscale("log")
      
      if 'logplot_y' in figpar.keys():
         if figpar['logplot_y']:
            ax2.set_yscale("log")
         # print('plot log')



      # print('varx',varx,len(varx))
      # print('slice_1D',slice_1D,len(slice_1D))

      print('mesh number',yml['study']['iter'])

   
      if 'plot_ref' in figpar.keys() and yml['study']['iter'] == 0:
         print('plotting ref')
         ref = eval(figpar['plot_ref'])
         # print('ref',ref)
         # print(4* yml["flower"]["physics"]["v_inlet"]*x_1D*scale_x/(mesh["xmax"]-mesh["xmin"])*(1-x_1D*scale_x/(mesh["xmax"]-mesh["xmin"])))
         # print(yml["flower"]["physics"]["v_inlet"])
         # print((mesh["xmax"]-mesh["xmin"]))
         print('i test',i,figpar['linestyles'][i+1],figpar['linestyles'][i])
         # ls  = eval(figpar['linestyles'][i+1])
         ax20.plot(varx, ref, 
         'k',
         label='Reference solution',
         ls=eval(figpar['linestyles'][i+1]),
         lw=lw)
         print('ref',ref)
         print('len ref',len(ref))

         print('y_1D',y_1D)
         print('y_1D',y_1D*scale_x)


      if 'macro' in figpar.keys():
         # X = varx
         # print(X)
         # local_context = {}
         exec(figpar['macro'],
            #   ,globals(),
            # globals(),
            # # locals(), 
            # local_context
            )
         
         # label1 = local_context['label1']
         label1 = label2

   
      #    print(label1)
      # print(label2)
      # print(label1)
      # print(figpar['macro'])
      # print(figpar)
      if 'logplot' in figpar.keys():

         ax20.scatter(x=varx, y=slice_1D, 
                     #  s=10,
                     #  marker='+',
         #  colors[i+1], #color wrt variable
         color=colors[(yml['study']['iter'])%len(colors)],
         # cmap=cmap,
         label=label1,
         # ls=ls, #creates bug
         lw=lw)

         figpar['error_list'].append(slice_1D)
         figpar['x_list'].append(varx)

      

      else:
         p1, = ax20.plot(varx, slice_1D, 
         #  colors[i+1], #color wrt variable
         colors[(yml['study']['iter'])%len(colors)],
         # cmap=cmap,
         label=label1,ls=ls,lw=lw)

      ax20.set(
      # xlim=(0, 2),
      # ylim=(0, 2),
      xlabel=labelx, #r""+plotpar['xlabel'],
      ylabel=label_i)
         
      if 'legend_pos' in figpar.keys():
         plt.legend(loc=figpar['legend_pos'])
      else:
         plt.legend()

   # tick0 = list(eval(figpar['ticks'][0]))
   # ax20.yaxis.set_major_locator(mticker.FixedLocator(tick0))
   # # ax2.yaxis.set_minor_locator(mticker.FixedLocator(tick0))
   # ax20.yaxis.set_ticks(tick0)


   # twin1.yaxis.set_major_locator(mticker.FixedLocator(eval(figpar['ticks'][1])))
   # twin1.yaxis.set_ticks(eval(figpar['ticks'][1]))

   # twin2.yaxis.set_major_locator(mticker.FixedLocator(eval(figpar['ticks'][2])))

   # twin2.yaxis.set_ticks(eval(figpar['ticks'][2]))

   # ax20.set_title('Time '+r"$\SI[retain-zero-exponent=true]{{{0:.2e}}}".format(time/plotpar['scale_time'])+'{'+plotpar['unit_time']+'}$')

   # ax20.yaxis.label.set_color(p1.get_color())
   # twin1.yaxis.label.set_color(p2.get_color())
   # twin2.yaxis.label.set_color(p3.get_color())

   # twin1.spines["right"].set_color(p2.get_color())
   # twin2.spines["right"].set_color(p3.get_color())

   # ax20.set(
   # # xlim=(0, 2),
   # # ylim=(0, 2),
   # xlabel=r""+plotpar['ylabel'],
   # ylabel=label1)
   # twin1.set(
   #    # ylim=(0, 4), 
   # ylabel=label2)
   # twin2.set(
   #    # ylim=(1, 65), 
   # ylabel=label3)

   # ax20.tick_params(axis="y", right = False, colors=p1.get_color())
   # twin1.tick_params(axis="y", right = True, labelright = True, left = False, labelleft = False, colors=p2.get_color())
   # twin2.tick_params(axis="y", right = True, labelright = True, left = False, labelleft = False, colors=p3.get_color())

   # twin1.yaxis.set_label_position("right")
   # twin2.yaxis.set_label_position("right")


   if 'plot_legend' in figpar.keys():
      if parse_is_true(figpar['plot_legend']):
         fig1.legend(handles=[p1, p2, p3],
         # loc = "center left",
         loc = "outside upper left",
         )


   ###########################################


   # return(fig1,ax2,cbar)


   # if mode =='first' or mode =='close':

   
   #    str_nstep = str(nstep)
   #    plt.savefig(file_name+'_'+str_nstep+ "." + plotpar["img_format"],dpi=plotpar['dpi']) #also for film for latex display

   # if mode == 'close':
   #    # str_nstep = str(nstep)
   #    # plt.savefig(file_name+'_'+str_nstep+ "." + plotpar["img_format"],dpi=plotpar['dpi'])
   #    plt.close(fig1)
   #    return

   if mode == 'first': 
      if len(figpar['var'])>1:
         ax2list = [ax2,twin1,twin2]
      else:
         ax2list = ax2
      # print(mode)
      # print(ax2list)
      return(fig1,ax2list,cbar)    
   else:
      # print(ax2)
      # print([ax20,twin1,twin2])
      return(fig1,ax2,cbar)
   


      



