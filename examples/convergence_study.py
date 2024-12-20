import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import yaml
import sys
import h5py
import math

# module to plot files from Flower.jl
# from plot_flower import * 
from plot_flower import set_size, init_fig, compute_slope, roundlog, logticks 



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

colors= [
   "#0072B2", #bleu
   "#D55E00", #orange
   "#009E73", #vert
]

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

   plt.savefig("errors.pdf")
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
         compute_slope(ax2,xls,yls,x,y,slope_and_correlation,R2,param_line,color_line,alpha,plot_text=False)
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

         print(figpar['file'])
         print(figpar)

         if 'radius' in figpar['var']: #we do not plot the figures with radius
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
                     df[err] = file[err][()]

                  print(df)

                  plot_errors_from_pandas(df,figpar,plotpar,colors,file_name)