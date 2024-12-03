import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

from matplotlib.patches import Rectangle

plt.rcParams['text.usetex'] = True
plt.rcParams["font.size"] = "11"

font_size = 14
font_size = 18

font_size = 24


plt.rcParams["font.size"] = str(font_size)


plt.rc('text.latex', preamble="\n".join([ # plots will use this preamble
        r'\usepackage{amsmath}',
        r'\usepackage{booktabs}',
        r"\usepackage{siunitx}",
        r"\setlength{\abovedisplayskip}{0pt}",
        r" \setlength{\belowdisplayskip}{0pt}",
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

bubblecolor = "#009E73" 

wallcolor = "#CC79A7" #pink
wallcolor = orange 


fig, ax = plt.subplots(layout='constrained')

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

plt.plot([0.,lx], [0.,0.], '-k', lw=4,zorder=zorder)
plt.plot([0.,0.], [0.,ly], '-k', lw=4,zorder=zorder)
plt.plot([0.,lx], [0.,0.], '-k', lw=4,zorder=zorder)
plt.plot([0.,lx], [ly,ly], '-k', lw=4,zorder=zorder)
plt.plot([lx,lx], [0.,ly], '-k', lw=4,zorder=zorder)

# p-grid
for i in range(0,nx+1):
    plt.plot([(i)*dx, (i)*dx],[0., ly],'-k',lw=2.0)
for j in range(0,ny+1):
    plt.plot([0., lx],[(j)*dy, (j)*dy],'-k',lw=2.0)
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

eps = 0.1


lx1    = 2
ly1    = 3

ly0 = 1

# plt.plot([0.,lx], [0.,0.], '-',  color = bubblecolor , lw=4)
# plt.plot([0.,0.], [ly0,ly], '-',  color = bubblecolor , lw=4)
plt.plot([0.,lx1], [ly0-eps,ly0-eps], '-',  color = bubblecolor , lw=4)
plt.plot([0.,lx1], [ly1+eps,ly1+eps], '-',  color = bubblecolor , lw=4)
plt.plot([lx1+eps,lx1+eps], [ly0,ly1], '-',  color = bubblecolor , lw=4)

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


ax.annotate('Bubble',
            # (1,2),
            (1,2),
            fontsize=font_size,
            c=bubblecolor ,
            ha="left",
            va='bottom')

#wall
plt.plot([0.-eps,0.-eps], [0,ly], '-', color = wallcolor, lw=4)

ax.annotate('Interface shifted',
            # (1,2),
            (1,3+eps),
            fontsize=font_size,
            c=bubblecolor ,
            ha="left",
            va='bottom')

ax.annotate(r"$c_1^\gamma$",
            # (1,2),
            (0.5,3+eps),
            fontsize=font_size,
            c=bubblecolor ,
            ha="left",
            va='bottom')



x2=0 #-0.1
ly2 = ly+0.5
ly2=ly1
ly2=ly
# ly2=ly+1

text = r"\begin{tabular}{c}" + \
        r" Wall at $0^-$\\"+ \
        r" $c_0^\gamma$ \\"+ \
        r"\end{tabular}"
ax.annotate(
    text,
    (x2,ly2),
    fontsize=font_size,
    c=wallcolor,
    ha="right",
    va='center')

# ax.annotate(r'Wall at $0^-$',
#             (x2,ly2),
#             fontsize=font_size,
#             c=wallcolor,
#             ha="right",
#             va='bottom')

# ax.annotate(r"$c_0^\gamma$",
#             # (1,2),
#             (0-eps,0),
#             fontsize=font_size,
#             c=wallcolor,
#             ha="center",
#             va='top')

# ax.annotate(r"$c_0^\gamma$",
#             # (1,2),
#             (0-eps,3.5),
#             fontsize=font_size,
#             c=wallcolor,
#             ha="right",
#             va='bottom')

# ax.annotate(r'Wall at $0^- \\ c_0^\gamma$',
#             (x2,ly2),
#             fontsize=font_size,
#             c=wallcolor,
#             ha="right",
#             va='bottom')

plt.axis('equal')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)

ax.set_xlim([-1.5,4])


# plt.show()
# plt.savefig('mesh_square_bubble.pdf') #Documenter does not look for .pdf
# plt.savefig('mesh_square_bubble.png')
plt.savefig('mesh_square_bubble.svg')

ax.add_patch(Rectangle((-1, 3), 1-eps, 1, fill=False, 
                    #    hatch=h,
                    color=wallcolor,
                    lw=4,
                    zorder=10,
                       ))

ax.add_patch(Rectangle((0, 3+eps), 1, 1-eps, fill=False, 
                    #    hatch=h,
                    color=bubblecolor ,
                    lw=4,
                    zorder=10,
                       ))

ax.set_xlim([-1.5,4])

plt.savefig('mesh_square_bubble_2.svg')

