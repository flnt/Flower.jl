import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

from matplotlib.patches import Rectangle

import matplotlib.transforms as transforms

from matplotlib.patches import Arc



plt.rcParams['text.usetex'] = True
plt.rcParams["font.size"] = "11"

font_size = 14
font_size = 18

font_size = 24

lw=4


plt.rcParams["font.size"] = str(font_size)


plt.rc('text.latex', preamble="\n".join([ # plots will use this preamble
        r'\usepackage{amsmath}',
        r'\usepackage{booktabs}',
        r"\usepackage{siunitx}",
        r"\setlength{\abovedisplayskip}{0pt}",
        r"\setlength{\belowdisplayskip}{0pt}",
        r"\setlength{\belowdisplayshortskip}{0pt}",
        r"\setlength{\abovedisplayshortskip}{0pt}",
        r"\addtolength{\jot}{-4pt}",
        # r"\usepackage{mhchem}",
        r"\usepackage[version=]{mhchem}",
        # r"\usepackage[version=4,arrows=pgf-filled,textfontname=sffamily,mathfontname=mathsf]{mhchem}",
       ])
)

# Latex

plt.rc("text", usetex=True)
plt.rc('text.latex', preamble=r"\usepackage{siunitx}")
#  matplotlib.verbose.level = 'debug-annoying'



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


lightorange = "#E69F00"
lightblue = "#56B4E9"

blue = "#0072B2"
orange = '#D55E00'

plot_interface = False
plot_wall = 'detailed'
plot_wall = 'nodetails'

case = 'Poisson_square'


# ax.set_facecolor(lightblue) #background color


nx    = int(4)                       # number of mesh cells in x
ny    = int(4)                       # number of mesh cells in y
lx    = 4                    # domain length in x
ly    = 4                       # domain length in y
dx    = lx/(nx)                  # mesh spacing in x
dy    = ly/(ny)                  # mesh spacing in y
xp     = np.arange(0,lx,dx)
yp     = np.arange(0,ly,dy)

zorder=1

eps = 0.1

lx1    = 2
ly1    = 3

ly0 = 1


def plot_figure(plot_interface,plot_wall,case):

    if 'BC' in case:
        fig, ax = plt.subplots(layout='constrained')
    else:
        fig, ax = plt.subplots(layout='constrained',figsize=(4,4))

    if '_BC' in case:
        radius = 0.5
    if '_radius' in case:
        radius = 4*1/2.5
    else:
        radius = 1




    plt.plot([0.,lx], [0.,0.], '-k', lw=lw,zorder=zorder)
    plt.plot([0.,0.], [0.,ly], '-k', lw=lw,zorder=zorder)
    plt.plot([0.,lx], [0.,0.], '-k', lw=lw,zorder=zorder)
    plt.plot([0.,lx], [ly,ly], '-k', lw=lw,zorder=zorder)
    plt.plot([lx,lx], [0.,ly], '-k', lw=lw,zorder=zorder)

    # p-grid
    # for i in range(0,nx+1):
    #     plt.plot([(i)*dx, (i)*dx],[0., ly],'-k',lw=2.0)
    # for j in range(0,ny+1):
    #     plt.plot([0., lx],[(j)*dy, (j)*dy],'-k',lw=2.0)
    # # u-grid
    # for i in range(0,nx+2):
    #     plt.plot([(i-0.5)*dx, (i-0.5)*dx],[0., ly],'--r',lw=1.0)
    # for j in range(0,ny+1):
    #     plt.plot([-dx/2., lx+dx/2.],[(j-0.0)*dy, (j-0.0)*dy],'--r',lw=1.0)
    # # v-grid
    # for i in range(0,nx+1):
    #     plt.plot([(i-0.0)*dx, (i-0.0)*dx],[-dy/2., ly+dy/2.],'--g',lw=1.0)
    # for j in range(0,ny+2):
    #     plt.plot([0., lx],[(j-0.5)*dy, (j-0.5)*dy],'--g',lw=1.0)



    if plot_interface:
        # plt.plot([0.,lx], [0.,0.], '-g', lw=lw)
        # plt.plot([0.,0.], [ly0,ly], '-g', lw=lw)
        plt.plot([0.,lx1], [ly0-eps,ly0-eps], '-g', lw=lw)
        plt.plot([0.,lx1], [ly1+eps,ly1+eps], '-g', lw=lw)
        plt.plot([lx1+eps,lx1+eps], [ly0,ly1], '-g', lw=lw)

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


    if plot_interface:
        ax.annotate('Bubble',
                    # (1,2),
                    (1,2),
                    fontsize=font_size,
                    c='g',
                    ha="left",
                    va='bottom')

        #wall
        plt.plot([0.-eps,0.-eps], [0,ly], '-r', lw=lw)

    if plot_interface:

        ax.annotate('Interface shifted',
                    # (1,2),
                    (1,3+eps),
                    fontsize=font_size,
                    c='g',
                    ha="left",
                    va='bottom')

        ax.annotate(r"$c_1^\gamma$",
                    # (1,2),
                    (0.5,3+eps),
                    fontsize=font_size,
                    c='g',
                    ha="left",
                    va='bottom')



    x2=0 #-0.1
    ly2 = ly+0.5
    ly2=ly1
    ly2=ly
    # ly2=ly+1


    if case == 'tuto':
        text = r"\begin{tabular}{c}" + \
                r" Wall at $0^-$\\"+ \
                r" $c_0^\gamma$ \\"+ \
                r"\end{tabular}"
        ax.annotate(
            text,
            (x2,ly2),
            fontsize=font_size,
            c='r',
            ha="right",
            va='center')
    elif '_BC' in case:

        print(case)

        xtmp = 0
        ytmp = ly/2

        offset = lw

        # # shift the object over 2 points, and down 2 points
        # shift = 10
        # dx = -shift/72.
        # dy = 0
        # offset = transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
        # shadow_transform = ax.transData + offset

        

        # now plot the same data with our offset transform;
        # use the zorder to make sure we are below the line
        # ax.plot(x, y, lw=3, color='gray',
        #         transform=shadow_transform,
        #         zorder=0.5*line.get_zorder())


        text = r"\begin{tabular}{c}" + \
                r" Neumann \\"+ \
                r" $\frac{\partial c}{\partial n}\neq 0$ \\"+ \
                r"\end{tabular}"
        ax.annotate(
            text,
            (xtmp,ytmp),
            xytext=(-offset, 0), textcoords='offset points',
            # transform=shadow_transform,
            fontsize=font_size,
            c='k',
            ha="right",
            va='center')    
        
        xtmp = lx/2
        ytmp = 0
        # dx = 0
        # dy = -shift/72.
        # offset = transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
        # shadow_transform = ax.transData + offset

        text = r"\begin{tabular}{c}" + \
                r" Neumann \\"+ \
                r" $\frac{\partial c}{\partial n} = 0$ \\"+ \
                r"\end{tabular}"
        ax.annotate(
            text,
            (xtmp,ytmp),
            xytext=(0,-offset), textcoords='offset points',
            # transform=shadow_transform,
            fontsize=font_size,
            c='k',
            ha="center",
            va='top')    
        
        xtmp = lx/2
        ytmp = ly

        # dx = shift/72.
        # dy = 0
        # offset = transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
        # shadow_transform = ax.transData + offset

        text = r"\begin{tabular}{c}" + \
                r" Neumann \\"+ \
                r" $\frac{\partial c}{\partial n} = 0$ \\"+ \
                r"\end{tabular}"
        ax.annotate(
            text,
            (xtmp,ytmp),
            xytext=(0,+offset), textcoords='offset points',
            # transform=shadow_transform,
            fontsize=font_size,
            c='k',
            ha="center",
            va='bottom')  

        xtmp = lx
        ytmp = ly/2

        # dx = 0
        # dy = shift/72.
        # offset = transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
        # shadow_transform = ax.transData + offset

        text = r"\begin{tabular}{c}" + \
                r" Dirichlet \\"+ \
                r" $c=0$ \\"+ \
                r"\end{tabular}"
        ax.annotate(
            text,
            (xtmp,ytmp),
            xytext=(+offset,0), textcoords='offset points',
            # transform=shadow_transform,
            fontsize=font_size,
            c='k',
            ha="left",
            va='center')    
        

        if 'arc_BC_dir' in case:
            
            xtmp = radius
            ytmp = ly/2

            # dx = 0
            # dy = shift/72.
            # offset = transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
            # shadow_transform = ax.transData + offset

            text = r"\begin{tabular}{c}" + \
                    r" Dirichlet \\"+ \
                    r" $c\neq 0$ \\"+ \
                    r"\end{tabular}"
            ax.annotate(
                text,
                (xtmp,ytmp),
                xytext=(+offset,0), textcoords='offset points',
                # transform=shadow_transform,
                fontsize=font_size,
                c='k',
                ha="left",
                va='center') 
            
        elif '_BC_dir' in case:
            
            if '_radius' in case:
                xtmp = lx/2
                ytmp = ly/2
                ha="center"
                va='center'
            else:
                xtmp = lx/2
                ytmp = ly/2-radius
                ha="center"
                va='top'

            # dx = 0
            # dy = shift/72.
            # offset = transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
            # shadow_transform = ax.transData + offset

            text = r"\begin{tabular}{c}" + \
                    r" Dirichlet \\"+ \
                    r" $c\neq 0$ \\"+ \
                    r"\end{tabular}"
            ax.annotate(
                text,
                (xtmp,ytmp),
                xytext=(+offset,0), textcoords='offset points',
                # transform=shadow_transform,
                fontsize=font_size,
                c='k',
                ha=ha,
                va=va,
            )    

        if 'arc_BC_neu' in case:
            
            xtmp = radius
            ytmp = ly/2

            # dx = 0
            # dy = shift/72.
            # offset = transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
            # shadow_transform = ax.transData + offset

            text = r"\begin{tabular}{c}" + \
                    r" Neumann \\"+ \
                    r" $\frac{\partial c}{\partial n}\neq 0$ \\"+ \
                    r"\end{tabular}"
            ax.annotate(
                text,
                (xtmp,ytmp),
                xytext=(+offset,0), textcoords='offset points',
                # transform=shadow_transform,
                fontsize=font_size,
                c='k',
                ha="left",
                va='center')    
            
        elif '_BC_neu' in case:
            
            if '_radius' in case:
                xtmp = lx/2
                ytmp = ly/2
                ha="center"
                va='center'
            else:
                xtmp = lx/2
                ytmp = ly/2-radius
                ha="center"
                va='top'

            # dx = 0
            # dy = shift/72.
            # offset = transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
            # shadow_transform = ax.transData + offset

            text = r"\begin{tabular}{c}" + \
                    r" Neumann \\"+ \
                    r" $\frac{\partial c}{\partial n}\neq 0$ \\"+ \
                    r"\end{tabular}"
            ax.annotate(
                text,
                (xtmp,ytmp),
                xytext=(+offset,0), textcoords='offset points',
                # transform=shadow_transform,
                fontsize=font_size,
                c='k',
                ha=ha,
                va=va,
                )    


        



    # ax.annotate(r'Wall at $0^-$',
    #             (x2,ly2),
    #             fontsize=font_size,
    #             c='r',
    #             ha="right",
    #             va='bottom')

    # ax.annotate(r"$c_0^\gamma$",
    #             # (1,2),
    #             (0-eps,0),
    #             fontsize=font_size,
    #             c='r',
    #             ha="center",
    #             va='top')

    # ax.annotate(r"$c_0^\gamma$",
    #             # (1,2),
    #             (0-eps,3.5),
    #             fontsize=font_size,
    #             c='r',
    #             ha="right",
    #             va='bottom')

    # ax.annotate(r'Wall at $0^- \\ c_0^\gamma$',
    #             (x2,ly2),
    #             fontsize=font_size,
    #             c='r',
    #             ha="right",
    #             va='bottom')

    plt.axis('equal')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    if 'BC' in case:
        ax.set_xlim([-1.5,4])


    # plt.show()
    # plt.savefig('mesh_square_bubble.pdf') #Documenter does not look for .pdf
    # plt.savefig('mesh_square_bubble.png')
    # plt.savefig('mesh_square_bubble.svg')

    if plot_wall == 'detailed':
        ax.add_patch(Rectangle((-1, 3), 1-eps, 1, fill=False, 
                            #    hatch=h,
                            color='r',
                            lw=lw,
                            zorder=10,
                            ))

    if plot_interface:
        ax.add_patch(Rectangle((0, 3+eps), 1, 1-eps, fill=False, 
                            #    hatch=h,
                            color='g',
                            lw=lw,
                            zorder=10,
                            ))

    if 'BC' in case:
        ax.set_xlim([-1.5,4])


 
    if '_circle' in case and ( not ('arc' in case)):

        xc = 2
        yc = 2
        

        circle1 = plt.Circle((xc, yc), radius, fill=False, color='k',hatch="\\",lw=lw,)

        ax.add_patch(circle1)

        if 'arrow' in case:
            arrow_length = radius/2
            vec1 = [arrow_length,arrow_length]
            thetas = np.linspace(0, 2*np.pi, 10)
            for theta in thetas:
               
                # plt.arrow(xc, yc, xc+vec1[0], yc+vec1[1], head_width=0.15, color='k', length_includes_head=True)
                plt.arrow(xc+radius*np.cos(theta), 
                          yc+radius*np.sin(theta), 
                          vec1[0]*np.cos(theta), 
                          vec1[1]*np.sin(theta), 
            
                          head_width=0.15, color='k', length_includes_head=True)
                
    elif '_arc' in case:
        xc = 0
        yc = 2
        thetalim = [-90,90]
        narrows = 5

        arc = Arc((xc, yc), radius*2, radius*2, color='k', theta1=thetalim[0], theta2=thetalim[1],hatch="\\",lw=lw)
        ax.add_patch(arc)

        if 'arrow' in case:
            arrow_length = radius/2
            vec1 = [arrow_length,arrow_length]
            thetas = np.linspace(thetalim[0]*np.pi/180, thetalim[1]*np.pi/180, narrows)

            for theta in thetas:
                
                # plt.arrow(xc, yc, xc+vec1[0], yc+vec1[1], head_width=0.15, color='k', length_includes_head=True)
                plt.arrow(xc+radius*np.cos(theta), 
                            yc+radius*np.sin(theta), 
                            vec1[0]*np.cos(theta), 
                            vec1[1]*np.sin(theta), 
            
                            head_width=0.15, color='k', length_includes_head=True)


    plt.savefig(case+'.svg', transparent=True)
    plt.savefig(case+'.pdf', transparent=True)

    plt.close('all')


plot_figure(plot_interface,plot_wall,'Poisson_square_BC')

plot_figure(plot_interface,plot_wall,'Poisson_square_circle_BC')

plot_figure(plot_interface,plot_wall,'Poisson_square_circle_BC_dir')

plot_figure(plot_interface,plot_wall,'Poisson_square_circle_BC_neu')

plot_figure(plot_interface,plot_wall,'Poisson_square_circle_arc_BC_dir')

plot_figure(plot_interface,plot_wall,'Poisson_square_circle_arc_BC_neu')

#For same proportions as radius = 1 L = 2.5
plot_figure(plot_interface,plot_wall,'Poisson_square_circle_BC_radius')

plot_figure(plot_interface,plot_wall,'Poisson_square_circle_BC_dir_radius')

plot_figure(plot_interface,plot_wall,'Poisson_square_circle_BC_neu_radius')

plot_figure(plot_interface,plot_wall,'Poisson_square_circle_arc_BC_dir_radius')

plot_figure(plot_interface,plot_wall,'Poisson_square_circle_arc_BC_neu_radius')


plot_figure(plot_interface,plot_wall,'Poisson_square')

plot_figure(plot_interface,plot_wall,'Poisson_square_circle')

plot_figure(plot_interface,plot_wall,'Poisson_square_circle_arrow')

plot_figure(plot_interface,plot_wall,'Poisson_square_circle_arc')

plot_figure(plot_interface,plot_wall,'Poisson_square_circle_arc_arrow')
