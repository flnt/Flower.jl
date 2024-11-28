import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm


# Latex
font_size = 14
font_size = 18

font_size = 24

plt.rcParams['text.usetex'] = True
plt.rcParams["font.size"] = str(font_size)

# plt.rc("text", usetex=True)
# plt.rc('text.latex', preamble=r"\usepackage{siunitx}")
#  matplotlib.verbose.level = 'debug-annoying'


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

def apply_font(fontpath1):

    fontpath2 = 'public/tex-gyre/texgyrepagella-regular.otf'
    fontpath = fontpath1 + fontpath2

    print('looking for font at',fontpath)

    fe = fm.FontEntry( 
    fname=fontpath,
    name='TeX Gyre Pagella Math'
    )
    fm.fontManager.ttflist.insert(0, fe) # or append is fine
    plt.rcParams['font.family'] = fe.name # = 'your custom ttf font name'

    prop = fm.FontProperties(fname=fontpath)

fontpath2 = 'public/tex-gyre/texgyrepagella-regular.otf'

try :
    fontpath1 = '/usr/share/fonts/opentype/'
    fontpath = fontpath1 + fontpath2
    os.path.isfile(fontpath)
except:
    fontpath1 = '/gpfs/workdir/regnaultp/latex/texmf-dist/fonts/opentype/'
    fontpath = fontpath1 + fontpath2
    os.path.isfile(fontpath)


# apply_font(fontpath)

OkabeIto=["#E69F00", #0 orange clair 230, 159, 0
"#56B4E9", #1 bleu clair 86, 180, 233
"#009E73", #2 vert 0, 158, 115
"#F0E442", #3 jaune 240, 228, 66
"#0072B2", #4 bleu 0, 114, 178
"#D55E00", #5 orange 213, 94, 0
"#CC79A7", #6 rose 204, 121, 171
"#000000"] #7 noir 0 0 0

colors=OkabeIto

colors=["#000000" for color in OkabeIto]

colors[1]="#000000"
colors[2]=OkabeIto[5] #bleu
colors[3]=OkabeIto[6] #orange

lightorange = "#E69F00"
lightblue = "#56B4E9"

blue = "#0072B2"
orange = '#D55E00'

def demo_con_style(ax, connectionstyle,points,color='k',mutation_scale=5):

    # x1, y1 = 0.3, 0.2
    # x2, y2 = 0.8, 0.6
    x1, y1 = points[0]
    x2, y2 = points[1]

    ax.plot([x1, x2], [y1, y2], ".")
    ax.annotate("",
                xy=(x1, y1), xycoords='data',
                xytext=(x2, y2), textcoords='data',
               
                arrowprops=dict(mutation_scale=mutation_scale,arrowstyle="->",linewidth=4, color=color,
                                shrinkA=5, shrinkB=5,
                                patchA=None, patchB=None,
                                connectionstyle=connectionstyle,
                                ),
                )
    
    

    # ax.text(.05, .95, connectionstyle.replace(",", ",\n"),
    #         transform=ax.transAxes, ha="left", va="top")

def plot_fig(Vmode,Wx,plot_wall_annotation=False):

    if Vmode != '':
        V = True
    else:
        V = False

    figname0=''
    if V:
        figname0+='V'
        figname0+=Vmode
    
    if Wx:
        figname0+='Wx'

    fig, ax = plt.subplots(figsize=(16,4))

    # ax.set_facecolor(lightblue) #background color

    nx    = int(4)                       # number of mesh cells in x
    ny    = int(4)                       # number of mesh cells in y
    lx    = 4                    # domain length in x
    ly    = 4                       # domain length in y
    dx    = lx/(nx)                  # mesh spacing in x
    dy    = ly/(ny)                  # mesh spacing in y
    xp     = np.arange(0,lx,dx)
    yp     = np.arange(0,ly,dy)


    startx = 0
    startx = -1

    sy = ny-1

    eps = 0.1


    lx1    = 2
    ly1    = 3

    ly0 = 1

    ly00 = 3

    lx00 = -1

    

    plt.axis('equal')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)


    plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off

    plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left=False,      # ticks along the bottom edge are off
    right=False,         # ticks along the top edge are off
    labelleft=False) # labels along the bottom edge are off


    # p-grid
    for i in range(startx,nx+1):
        plt.plot([(i)*dx, (i)*dx],[ly00, ly],'-k',lw=2.0)
    for j in range(sy,ny+1):
        plt.plot([lx00, lx],[(j)*dy, (j)*dy],'-k',lw=2.0)


    #wall
    plt.plot([0.-eps,0.-eps], [ly00,ly], '-r', lw=4)




    x2=0 #-0.1
    ly2 = ly+0.5
    ly2=ly1
    ly2=ly
    # ly2=ly+1

    yarr = [ly00,ly]


    if plot_wall_annotation:
        ax.annotate(r'Wall at $0^-$',
                    (x2,ly2),
                    fontsize=font_size,
                    c='r',
                    ha="right",
                    va='bottom')

    # ax.annotate(r"$c_0^\gamma$",
    #             # (1,2),
    #             (0-eps,3.5),
    #             fontsize=font_size,
    #             c='r',
    #             ha="left",
    #             va='bottom',)



   



    if V:
        ax.scatter(0.5,3.5,
            marker='x',
            color = lightorange, 
            )

        if Vmode == 'c':
            text = r"\begin{tabular}{c}" + \
                        r" $\mathcal{V}_{1}$ \\"+ \
                        r" $x^\omega$ \\"+ \
                        r"\end{tabular}"

            ax.annotate(
                        text,
                        # 'Centroid i=1',
                        (0.5,3.5),
                        fontsize=font_size,
                        c=lightorange,
                        ha="left",
                        va='center')
        # else:
            # text = r"\begin{tabular}{c}" + \
            #             r" $\mathcal{V}_{1}$ \\"+ \
            #             r" $h_xh_y$ \\"+ \
            #             r"\end{tabular}"

            # ax.annotate(
            #             text,
            #             # 'Centroid i=1',
            #             (0.5,3),
            #             fontsize=font_size,
            #             c=lightorange,
            #             ha="left",
            #             va='top')


        ax.fill_betweenx(yarr, 0, 1,
                        alpha = 0,
                        hatch='/',
                        )



    plt.savefig(figname0+'_'+'mesh_square_bubble_wall_1'+'.svg')

    # plt.show()


    if Wx:
        ax.fill_betweenx(yarr, 0-eps, 0.5,
                    alpha = 0,
                    hatch='-',
                    #  hatch='.',
                    )
        ax.scatter(0.25,3.5,
            marker='x',
            color = lightblue, 
            )

        # text = r"\begin{tabular}{c}" + \
        #             r" $\mathcal{W}_{x,0}$ \\"+ \
        #             r"\end{tabular}"
        #             # r" $x^\omega$ \\"+ \

        # ax.annotate(
        #             text,
        #             # 'Centroid i=1',
        #             (0.25,3.5),
        #             fontsize=font_size,
        #             c=lightblue,
        #             ha="left",
        #             va='center')
        # ax.annotate(
        #             r'$\mathcal{W}_{x,i=0}$',
        #             # 'Centroid i=0',
        #             (0,3.5),
        #             fontsize=font_size,
        #             c=blue,
        #             ha="left",
        #             va='top')

    plt.plot([0.5,0.5], [ly00,ly], '-',color=blue, lw=4)
    
    # text = r"\begin{tabular}{c}" + \
    #                 r" \toprule" + \
    #                 r" $\mathcal{B}_{x,1}$ \\"+ \
    #                 r" $h_y$ \\"+ \
    #                 r"\end{tabular}"
    # ax.annotate(
    #             text,
    #             # 'Centroid i=0',
    #             (0.5,3),
    #             fontsize=font_size,
    #             c=blue,
    #             ha="center",
    #             va='top')

    plt.savefig(figname0+'_'+'mesh_square_bubble_wall_2'+'.svg')


   
    if Wx:
        ax.fill_betweenx(yarr, 0.5,1.5,
                    alpha = 0,
                    #  hatch='-',
                    hatch='|',
                    )
       
        ax.scatter(1,3.5,
            marker='x',
            color = lightblue, 
            )

        # text = r"\begin{tabular}{c}" + \
        #             r" $\mathcal{W}_{x,1}$ \\"+ \
        #             r"\end{tabular}"
        #             # r" $x^\omega$ \\"+ \

        # ax.annotate(
        #             text,
        #             # 'Centroid i=1',
        #             (1,3.5),
        #             fontsize=font_size,
        #             c=lightblue,
        #             ha="left",
        #             va='center')
        
    if V:
        ax.fill_betweenx(yarr, -1, 0,
                        alpha = 0,
                        hatch='\\',
                        )
        
        ax.scatter(-0.5,3.5,
            marker='x',
            color = lightorange, 
            )


        if Vmode == 'c':
            text = r"\begin{tabular}{c}" + \
                    r" $\mathcal{V}_{0}$ \\"+ \
                    r" $0$ \\"+ \
                    r" $x^\omega$ \\"+ \
                    r"\end{tabular}"

            ax.annotate(
                    text,
                    # 'Centroid i=1',
                    (-0.5,3.5),
                    fontsize=font_size,
                    c=lightorange,
                    ha="left",
                    va='center')
        else:
            # text = r"\begin{tabular}{c}" + \
            #         r" $\mathcal{V}_{0}$ \\"+ \
            #         r" $0$ \\"+ \
            #         r"\end{tabular}"

            text = r"\begin{tabular}{cc}" + \
                    r" \toprule" + \
                    r" $\mathcal{B}_{x,0}$ & $\mathcal{V}_{x,0}$ \\"+ \
                    r" $0$ & $0$ \\"+ \
                    r"\end{tabular}"

            ax.annotate(
                    text,
                    # 'Centroid i=1',
                    (-0.5,3),
                    fontsize=font_size,
                    c='k',
                    ha="center",
                    va='top')
            
          


       

        ax.fill_betweenx(yarr, 1, 2,
                        alpha = 0,
                        hatch='\\',
                        )
        
        ax.scatter(1.5,3.5,
            marker='x',
            color = lightorange, 
            )
        
        if Vmode == 'c':
            text = r"\begin{tabular}{c}" + \
                    r" $\mathcal{V}_{2}$ \\"+ \
                    r" $h_xh_y$ \\"+ \
                    r" $x^\omega$ \\"+ \
                    r"\end{tabular}"

            ax.annotate(
                    text,
                    # 'Centroid i=1',
                    (1.5,3.5),
                    fontsize=font_size,
                    c=lightorange,
                    ha="left",
                    va='center')

        # else:
            # text = r"\begin{tabular}{c}" + \
            #         r" $\mathcal{V}_{2}$ \\"+ \
            #         r" $h_xh_y$ \\"+ \
            #         r"\end{tabular}"

            # ax.annotate(
            #         text,
            #         # 'Centroid i=1',
            #         (1.5,3),
            #         fontsize=font_size,
            #         c=lightorange,
            #         ha="left",
            #         va='top')

       
    # text = r"\begin{tabular}{c}" + \
    #                 r" \toprule" + \
    #                 r" $\mathcal{A}_{x,1}$ \\"+ \
    #                 r" $h_y$ \\"+ \
    #                 r"\end{tabular}"
    # ax.annotate(
    #             text,
    #             # 'Centroid i=0',
    #             (0,3),
    #             fontsize=font_size,
    #             c=orange,
    #             ha="center",
    #             va='top')

    text = r"\begin{tabular}{cc}" + \
            r" $\mathcal{A}_{x,1}$ & $\mathcal{W}_{x,1}$\\"+ \
            r" $h_y$ & $\frac{h_xh_y}{2}$ \\"+ \
            r"\end{tabular}"
    ax.annotate(
                text,
                # 'Centroid i=0',
                (0,4),
                fontsize=font_size,
                c='k',
                ha="center",
                va='bottom')
    

    plt.plot([-0.5,-0.5], [ly00,ly], '-',color=blue, lw=4)

    plt.plot([0,0], [ly00,ly], '-',color=orange, lw=4)
    
    text = r"\begin{tabular}{cc}" + \
                    r" \toprule" + \
                    r" $\mathcal{B}_{x,1}$ & $\mathcal{V}_{x,1}$ \\"+ \
                    r" $h_y$ & $h_xh_y$ \\"+ \
                    r"\end{tabular}"

                    # r"\bottomrule" + \

    ax.annotate(
                text,
                # 'Centroid i=0',
                (0.5,3),
                fontsize=font_size,
                c='k',
                ha="center",
                va='top')

    plt.savefig(figname0+'_'+'mesh_square_bubble_wall_3'+'.svg')


 
    # 0  

    plt.plot([-1,-1], [ly00,ly], '-',color=orange, lw=4)

    text = r"\begin{tabular}{cc}" + \
            r" $\mathcal{A}_{x,0}$ & $\mathcal{W}_{x,0}$\\"+ \
            r" $0$ & $0$ \\"+ \
            r"\end{tabular}"
    ax.annotate(
                text,
                # 'Centroid i=0',
                (-1,4),
                fontsize=font_size,
                c='k',
                ha="center",
                va='bottom')

    # 2
    plt.plot([1,1], [ly00,ly], '-',color=orange, lw=4)

    text = r"\begin{tabular}{cc}" + \
            r" $\mathcal{A}_{x,2}$ & $\mathcal{W}_{x,2}$\\"+ \
            r" $h_y$ & $h_xh_y$ \\"+ \
            r"\end{tabular}"
    ax.annotate(
                text,
                # 'Centroid i=0',
                (1,4),
                fontsize=font_size,
                c='k',
                ha="center",
                va='bottom')


    plt.savefig(figname0+'_'+'mesh_square_bubble_wall_3'+'.svg')

    # ax.annotate(
    #             r'$i=1$',
    #             # 'Centroid i=0',
    #             (1,4),
    #             fontsize=font_size,
    #             c=blue,
    #             ha="center",
    #             va='bottom')

    # text = r"\begin{alignat*}{4}" + \
    #                               r"\mathrm{{Slopes}}\,& \num[round-mode = places,round-precision=2]{%s}\,\, \num[round-mode = places,round-precision=2]{%s}\,\, \num[round-mode = places,round-precision=2]{%s}\\"% tuple(slopes[:]) + \
    #                               r"R^2\,& \num[round-mode = places,round-precision=2]{%s}\,\, \num[round-mode = places,round-precision=2]{%s}\,\, \num[round-mode = places,round-precision=2]{%s}"% tuple(R2[:]) + \
    #                               r"\end{alignat*}"

    text = r"\begin{alignat*}{1}" + \
        r"&\mathcal{A}_{x,0}\\" + \
        r"&0" + \
        r"\end{alignat*}"

    text = r"\begin{tabular}{c}" + \
                r" \toprule" + \
                r" $\mathcal{A}_{x,0}$ \\"+ \
                r" $0$ \\"+ \
                r"\end{tabular}"


    ax.annotate(
                text,
                # 'Centroid i=0',
                (0,3),
                fontsize=font_size,
                c=blue,
                ha="center",
                va='top')

    # ax.annotate(
    #             r'\noindent$\mathcal{A}_{x,0}\\0$',
    #             # 'Centroid i=0',
    #             (0,3),
    #             fontsize=font_size,
    #             c=blue,
    #             ha="center",
    #             va='top')


    text = r"\begin{alignat*}{1}" + \
        r"&\mathcal{A}_{x,1}\\" + \
        r"&h_y" + \
        r"\end{alignat*}"

    text = r"\begin{tabular}{c}" + \
                r" \toprule" + \
                r" $\mathcal{A}_{x,1}$ \\"+ \
                r" $h_y$ \\"+ \
                r"\end{tabular}"
    ax.annotate(
                text,
                # 'Centroid i=0',
                (1,3),
                fontsize=font_size,
                c=blue,
                ha="center",
                va='top')
    
    
    plt.plot([0,0.5], [3,3], '-g', lw=4,
            #  zorder = 2,
             )

    plt.plot([0.5,1], [3,3], '-g', lw=4,
                        #   zorder = 2,
                )


    ax.scatter(0.25,3,
        marker='x',
        color = 'k', 
        zorder=10,
        )
    
    ax.scatter(0.75,3,
        marker='x',
        color = 'k', 
        zorder =10,
        )


    mutation_scale = 30
    demo_con_style(ax, "angle3,angleA=90,angleB=0",points=[[0.25,3],[0.25,3.5]],color='k',mutation_scale=mutation_scale)

    demo_con_style(ax, "angle3,angleA=90,angleB=0",points=[[0.75,3],[1,3.5]],color='k',mutation_scale=mutation_scale)

    
    plt.savefig(figname0+'_'+'mesh_square_bubble_wall_4'+'.svg')


    # ax.annotate(
    #             r'$\mathcal{A}_{x,i=1}$',
    #             # 'Centroid i=0',
    #             (1,3),
    #             fontsize=font_size,
    #             c=blue,
    #             ha="center",
    #             va='top')

    # text = r"\begin{alignat*}{1}" + \
    #        r"&\mathcal{B}_{x,1}\\" + \
    #        r"&h_y" + \
    #        r"\end{alignat*}"

    # slopetext = r"\begin{tabular}{cccc}" + \
    #             r" \toprule" + \
    #             r" Slopes & \num[round-mode = places,round-precision=2]{%s} & \num[round-mode = places,round-precision=2]{%s} & \num[round-mode = places,round-precision=2]{%s} \\"% tuple(slopes[:]) + \
    #             r" \midrule" + \
    #             r"$R^2$ & \num[round-mode = places,round-precision=2]{%s} & \num[round-mode = places,round-precision=2]{%s} & \num[round-mode = places,round-precision=2]{%s} \\"% tuple(R2[:]) + \
    #             r"\bottomrule" + \
    #             r"\end{tabular}"
    # text = r"\begin{tabular}{c}" + \
    #             r" \toprule" + \
    #             r" $\mathcal{B}_{x,1}$ \\"+ \
    #             r" $h_y$ \\"+ \
    #             r"\end{tabular}"
    # ax.annotate(
    #             text,
    #             # 'Centroid i=0',
    #             (0.5,3),
    #             fontsize=font_size,
    #             c=blue,
    #             ha="center",
    #             va='top')


    # i =0

    # ax.scatter(-0.5,3.5,
    #            marker='x'
    #            )

    # ax.annotate(
    #             r'$\mathcal{B}_{x,i=1}$',
    #             # 'Centroid i=0',
    #             (0.5,3),
    #             fontsize=font_size,
    #             c=blue,
    #             ha="left",
    #             va='top')

    # {'/', '\', '|', '-', '+', 'x', 'o', 'O', '.', '*'}


    # plt.show()

    filename = 'mesh_square_bubble_wall'


    plt.savefig(figname0+'_'+filename+'.svg')

    # plt.show()

plot_fig(Vmode='c',Wx=False)

plot_fig(Vmode='',Wx=True)

plot_fig(Vmode='t',Wx=False)
