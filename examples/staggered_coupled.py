import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches
import matplotlib.transforms as transforms
import matplotlib.gridspec as gridspec
import os

#region Latex font

import matplotlib.font_manager as fm

##########################################################################
plt.rcParams['text.usetex'] = True
plt.rcParams["font.size"] = "14" #"11"

plt.rcParams["text.parse_math"] = False #necessary for mhchem


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
        # r"\usepackage[version=4]{mhchem}",
        # r"\usepackage[version=4,arrows=pgf-filled,textfontname=sffamily,mathfontname=mathsf]{mhchem}",
       ])
)

plt.rc("text", usetex=True)

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

apply_font(fontpath)

#endregion Latex font


orange_Okabe = "#D55E00"
blue_Okabe = "#0072B2"

p_color = 'k'
u_color = orange_Okabe
v_color = blue_Okabe

radius = 11
marker_size =radius**2
# marker_size = 100
# marker_size = 200

fontsize = 14

# Create a 3x3 meshgrid
x = np.array([0.5,1.5,2.5])

y = x
X, Y = np.meshgrid(x, y)

xu = [1,2] #x -0.5
Xu, Yu = np.meshgrid(xu, y)

yv = xu
Xv, Yv = np.meshgrid(x, yv)

subplots_mode = True
# subplots_mode = False
if subplots_mode:
        # Create a figure and two subplots
        fig, axs = plt.subplots(1, 3, figsize=(12, 6),width_ratios=[1, 0.1,1],layout="constrained")

        ax_left = axs[0]
        ax_right = axs[2]
        legend_ax = axs[1]
else:
        # Create a figure
        fig = plt.figure(figsize=(15, 5))

        # Create a GridSpec with specific width ratios for the subplots
        gs = gridspec.GridSpec(1, 3, width_ratios=[3, 1, 3],wspace=0)

        # Create subplots using the GridSpec
        ax_left = plt.subplot(gs[0])
        legend_ax = plt.subplot(gs[1])
        ax_right = plt.subplot(gs[2])


xu = [0.25,1,2,2.75]
# Plot on the first subplot
ax_left.scatter(X, Y, color=p_color,s=marker_size)
ax_left.scatter(Xu, Yu, color=u_color, marker='>',s=marker_size)
ax_left.scatter(Xv, Yv, color=v_color, marker='^',s=marker_size)


#
ax_left.scatter(x, [0]*len(x), color='None', marker='o', edgecolors=p_color,s=marker_size)
ax_left.scatter(x, [3]*len(x), color='None', marker='o', edgecolors=p_color,s=marker_size)

ax_left.scatter([0]*len(x),x, color='None', marker='o', edgecolors=p_color,s=marker_size)
ax_left.scatter([3]*len(x),x, color='None', marker='o', edgecolors=p_color,s=marker_size)
#


ax_left.scatter(xu, [0]*len(xu), color='None', marker='>', edgecolors=u_color,s=marker_size)
ax_left.scatter(xu, [3]*len(xu), color='None', marker='>', edgecolors=u_color,s=marker_size)
ax_left.scatter([0]*len(xu), xu, color='None', marker='^', edgecolors=v_color,s=marker_size)

x2 = np.array([0.5,1.5,2.5])
ax_left.scatter([0]*len(x2),x2, color='None', marker='>', edgecolors=u_color,s=marker_size)
ax_left.scatter([3]*len(x2),x2, color='None', marker='>', edgecolors=u_color,s=marker_size)
ax_left.scatter(x2,[0]*len(x2), color='None', marker='^', edgecolors=v_color,s=marker_size)
ax_left.scatter(x2,[3]*len(x2), color='None', marker='^', edgecolors=v_color,s=marker_size)

#region first half cells
ax_left.scatter([0.25]*len(x2),x2, color=u_color, marker='>',s=marker_size, 
        #        edgecolors=u_color,
               )
ax_left.scatter([2.75]*len(x2),x2, color=u_color, marker='>',s=marker_size, 
        #        edgecolors=u_color,
               )

ax_left.scatter(x2,[0.25]*len(x2), color=v_color, marker='^',s=marker_size, 
        #        edgecolors=u_color,
               )
ax_left.scatter(x2,[2.75]*len(x2), color=v_color, marker='^',s=marker_size, 
        #        edgecolors=u_color,
               )
#endregion first half cells





ax_left.scatter([3]*len(xu), xu, color='None', marker='^', edgecolors=v_color,s=marker_size)

xtext, ytext = 3, 0
linewidth_points = 1 #lw_inset
# ax = ax_right
dx, dy = 0.0, -linewidth_points/72.
offset = transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
shadow_transform = ax_left.transData + offset

# ax_right.annotate(r'$\Gamma$', xy=(xtext, ytext),
#                 xytext=(xtext, ytext),
#                 ha='center', va='top',
#                 fontsize=fontsize, color='black',
#                 transform=shadow_transform)

ax_left.text(xtext, ytext, r'$\Gamma$', fontsize=fontsize, ha='center', va='top',
        transform = shadow_transform)


# Set labels and title for the first subplot
ax_left.set_xlabel('X-axis')
ax_left.set_ylabel('Y-axis')
# ax_left.set_title('3x3 Mesh with Scatter Plot')

# Turn off the axes for the first subplot
ax_left.axis('off')

# Add square from (0,0) to (2,2) using a patch for the first subplot
square = patches.Rectangle((0, 0), 3, 3, linewidth=1, edgecolor='k', facecolor='none')
ax_left.add_patch(square)
ax_left.axis('equal')


xu = [1,2] #x -0.5


# Plot on the second subplot
ax_right.scatter(X, Y, color=p_color,s=marker_size)
ax_right.scatter(Xu, Yu, color=u_color, marker='>',s=marker_size)
ax_right.scatter(Xv, Yv, color=v_color, marker='^',s=marker_size)
ax_right.scatter(xu, [0]*len(xu), color='None', marker='>', edgecolors=u_color,s=marker_size)
ax_right.scatter(xu, [3]*len(xu), color='None', marker='>', edgecolors=u_color,s=marker_size)
ax_right.scatter([0]*len(xu), xu, color='None', marker='^', edgecolors=v_color,s=marker_size)
ax_right.scatter([3]*len(xu), xu, color='None', marker='^', edgecolors=v_color,s=marker_size)


x2 = np.array([0.5,1.5,2.5])
ax_right.scatter([0]*len(x2),x2, color='None', marker='>', edgecolors=u_color,s=marker_size)
ax_right.scatter([3]*len(x2),x2, color='None', marker='>', edgecolors=u_color,s=marker_size)
ax_right.scatter(x2,[0]*len(x2), color='None', marker='^', edgecolors=v_color,s=marker_size)
ax_right.scatter(x2,[3]*len(x2), color='None', marker='^', edgecolors=v_color,s=marker_size)

xtext, ytext = 3, 0
linewidth_points = 1 #lw_inset
# ax = ax_right
dx, dy = 0.0, -linewidth_points/72.
offset = transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
shadow_transform = ax_right.transData + offset

# ax_right.annotate(r'$\Gamma$', xy=(xtext, ytext),
#                 xytext=(xtext, ytext),
#                 ha='center', va='top',
#                 fontsize=fontsize, color='black',
#                 transform=shadow_transform)

ax_right.text(xtext, ytext, r'$\Gamma$', fontsize=fontsize, ha='center', va='top',
        transform = shadow_transform)

# Set labels and title for the second subplot
ax_right.set_xlabel('X-axis')
ax_right.set_ylabel('Y-axis')
# ax_right.set_title('3x3 Mesh with Scatter Plot')

# Turn off the axes for the second subplot
ax_right.axis('off')

# Add square from (0,0) to (2,2) using a patch for the second subplot
square = patches.Rectangle((0, 0), 3, 3, linewidth=1, edgecolor='k', facecolor='none')
ax_right.add_patch(square)
ax_right.axis('equal')

# # Draw an arrow between the two subplots
# # Coordinates are in figure coordinates (0,0 is bottom-left, 1,1 is top-right)
# fig.annotate('', xy=(0.6, 0.5), xytext=(0.4, 0.5),
#              xycoords='figure fraction', textcoords='figure fraction',
#              arrowprops=dict(facecolor='black', shrink=0.05))

# Create legend handles
interior_handle = ax_left.scatter([], [], color=p_color, label='bulk pressure',s=marker_size)
u_handle = ax_left.scatter([], [], color=u_color, marker='>', label='bulk u',s=marker_size)
v_handle = ax_left.scatter([], [], color=v_color, marker='^', label='bulk v',s=marker_size)
bc_handle_p = ax_left.scatter([], [], color='none',edgecolors=p_color, label='pressure BC',s=marker_size)
bc_handle = ax_left.scatter([], [], color='None', marker='>', edgecolors=u_color, label='u BC',s=marker_size)
bc_handle_v = ax_left.scatter([], [], color='None', marker='^', edgecolors=v_color, label='v BC',s=marker_size)

handles = [interior_handle, u_handle, v_handle, bc_handle_p,bc_handle,bc_handle_v]

# Add legend to both subplots
# ax_left.legend(handles=[interior_handle, u_handle, v_handle, bc_handle,bc_handle_v], loc='upper right')
# ax_right.legend(handles=[interior_handle, u_handle, v_handle, bc_handle], loc='upper right')


# legend_ax = axs[1]

# Create a separate Axes for the legend
# legend_ax = fig.add_axes([0.5, 0.1, 0.1, 0.8])  # Position and size of the legend box
# legend_ax.axis('off')  # Turn off the axes for the legend box

legend_ax.axis('off')



# Add legend to the separate Axes
legend_ax.legend(handles=handles, 
                #  loc='lower center',
                # loc='center',
                loc='upper center', bbox_to_anchor=(0.5, 0.4)
                 )

# legend_ax.annotate('', xy=(0.75, 0.5), xytext=(0.25, 0.5),
#                   xycoords='axes fraction', textcoords='axes fraction',
#                   arrowprops=dict(facecolor='black', arrowstyle='->'))

legend_ax.annotate('', xy=(1, 0.5), xytext=(0., 0.5),
                  xycoords='axes fraction', textcoords='axes fraction',
                  arrowprops=dict(facecolor='black', arrowstyle='->',linewidth=2))

# plt.tight_layout()
# Save the figure
plt.savefig('staggered_coupled_comparison.pdf', 
            transparent=True,
            )

# Show plot
plt.show()


