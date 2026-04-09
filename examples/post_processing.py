import numpy as np
# from pdi import access, Error, event, finalize, init, reclaim, release, version, share, OUT, IN, INOUT, NONE
import pdi

def share_test(name, data, access):
    print('[share_test]',name,data,access)
    if (isinstance(data, np.ndarray)):
        print("isinstance")
        data_np_array = data
    elif (access == OUT or access == NONE):
        # data is not numpy array
        print('access == OUT or access == NONE')
        try:
            data_np_array = np.array(data)
        except:
            raise Error("`" + name + "' share: Type is not supported by PDI, cannot insert it into numpy array")
    else:
        raise Error("`" + name + "' share: IN and INOUT can be only done with numpy array data type")
    print('data_np_array',data_np_array,type(data_np_array))
    try:
        share(name, data_np_array, access)
    except Exception as e:
        print(e)
        print('alignment',data_np_array.ctypes.data % 8)
    print('[end share_test]')

def multi_expose_test(event_name, expose_list):
    # from pdi import access, Error, event, finalize, init, reclaim, release, version, OUT, IN, INOUT, NONE
    # from pdi import share
    # import pdi._pdi
    import numpy as np
    import inspect
    exposed = []
    print('[multi_expose_test]')
    print('expose_list',expose_list)
    try:
        for (name, data, access) in expose_list:
            print('name',name)
            share_test(name, data, access) #share
            exposed.append(name)
        event(event_name)
        print('event ',event_name)
    except Exception as e:
        print('error multi_expose_test')
        print('Error code ',e)
        pass
    final_error = ()
    for name in exposed:
        try:
            reclaim(name)
        except Exception as e:
            final_error += (e)
    if (final_error != ()):
        raise final_error

"""
TODO constant mesh spacing
"""
def calculate_centroid(x, y, volume_cell):
    # x and y are scalar node coord

    # print(y)

    # print(volume_cell)

    integral_x = np.sum(x * volume_cell)
    integral_y = np.sum(y * volume_cell)
    integral_1 = np.sum(volume_cell)
    if integral_1 > 0.0:
        x_c = integral_x / integral_1
        y_c = integral_y / integral_1
    else:
        x_c = 0.0
        y_c = 0.0
        from termcolor import colored
        print(colored('error calculate_centroid','red'))

    return (x_c,y_c)

def calculate_circularity(perimeter_bubble, area):
    # area = pi r 2
    # r = sqrt(area/pi)
    # perim= 2 sqrt(area*pi)
    # perim = 2pi r
    # area = volume_fraction * dx dy ou dcap
    if perimeter_bubble > 0.0:
        perimeter_circle = 2 * np.sqrt(np.pi * area)
        circularity = perimeter_circle / perimeter_bubble
    else:
        circularity = 0.0
        from termcolor import colored
        print(colored('error calculate_circularity','red'))
    return circularity

def calculate_rise_velocity(v, volume_cell):
    # v velocity
    integral_u = np.sum(v * volume_cell)
    integral_1 = np.sum( volume_cell)
    if integral_1 > 0.0:
        U_c = integral_u / integral_1
    else:
        U_c = 0.0
        from termcolor import colored
        print(colored('error calculate_rise_velocity','red'))
    return U_c

def calculate_rise_velocity_plot(v, volume_cell,nx=None,ny=None):
    import matplotlib.pyplot as plt
    # v velocity
    integral_u = np.sum(v * volume_cell)
    integral_1 = np.sum( volume_cell)

    # np.set_printoptions(suppress=False, precision=2)
    # np.set_printoptions(formatter={'float': lambda x: format(x, '.2e')})

    # np.set_printoptions(edgeitems=30, linewidth=100000,
    #       formatter={'float': lambda x: format(x, '.2e')})

    # # print('integral_1',integral_1,3.14*0.25**2)
    # print('multiply',v )

    # print('multiply',v * volume_cell)

    nx = 40
    ny = 80
    xp = np.linspace(0,1,nx)
    yp = np.linspace(0,2,ny)
    # xu= np.linspace(0,1,nx+1)

    # yv= np.linspace(0,2,ny+1)

    fig,ax = plt.subplots()
    v = v.transpose()
    CS = ax.contourf(xp,yp,v)
    cbar = fig.colorbar(CS)

    plt.savefig('calculate_rise_velocity_plot.svg')

    fig,ax = plt.subplots()

    volume_cell = volume_cell.transpose()
    CS = ax.contourf(xp,yp,v * volume_cell)
    cbar = fig.colorbar(CS)

    plt.savefig('calculate_rise_velocity_plot_2.svg')

    # plt.show()


    # yml = None
    # mesh = None 
    # time = 0
    # nstep = 0 
    # plotpar = 0
    # figpar = 0
    
    # plot_zoom(
    # v,
    # xp,
    # yp,
    # xu,
    # yv,
    # yml,
    # time,
    # nstep,
    # plotpar,
    # figpar,
    # mode='close',
    # fig1=None,
    # ax2=None,
    # cbar=None,
    # )


    U_c = integral_u / integral_1

    return U_c

def find_one_minimum(slice,x_1D,eps): #TODO change order after
    # print('len',len(slice),slice)
    min_dist = np.min(abs(slice))
    min_dist_tmp = np.max(abs(slice))
    for i in range(len(slice)):
        current_abs = abs(slice[i])
        if current_abs <min_dist_tmp:
            min_dist_tmp = current_abs
            found_minimum = abs(current_abs - min_dist)<eps*min_dist
            # print('i',i,abs(slice[i]),min_dist,found_minimum)
            i1 = i
            if found_minimum:
                # print('min is at',i,slice[i],x_1D[i])
                # print(slice[i-1:i+1])
                # print(i-1,i+1)

                if abs(slice[i-1]) < abs(slice[i+1]):
                    i2 = i-1
                    # print('i-1',i-1,slice[i-1],x_1D[i-1])
                else:
                    i2=i+1
                    # print('i+1')                    
                    # print('i2',i2,slice[i2],x_1D[i2])
                break
    
    return(i1,i2)

def find_sign_changes(slice,x_1D,eps):
    # print('len',len(slice),slice)
    # min_dist = np.min(abs(slice))
    # min_dist_tmp = np.max(abs(slice))
    for i in range(len(slice)):
        if (slice[i] * slice[i+1]) < 0:
            i1 = i
            i2 = i+1
            break

    return(i1,i2)


def compute_radius_from_levelset_slice(slice,x_1D,eps):      
   
    dx = x_1D[1]-x_1D[0]

    # # print(colored('first','red'))
    # i1,i2 = find_one_minimum(slice,x_1D,eps)
    # print('i1 i2',i1,i2)

    # # print(colored('second','red'))
    # itmp = max(i1,i2)
    # # print('itmp',itmp)
    # slice2 = slice[itmp+1:] 
    # i3,i4 = find_one_minimum(slice2,x_1D,eps)
    # i3+= itmp+1
    # i4+= itmp+1
    # print('i3 i4',i3,i4)

    i1,i2 = find_sign_changes(slice,x_1D,eps)
    # print('i1 i2',i1,i2)
    itmp = max(i1,i2)
    slice2 = slice[itmp+1:]
    i3,i4 = find_sign_changes(slice2,x_1D,eps)
    i3+= itmp+1
    i4+= itmp+1
    # print('i3 i4',i3,i4)


    a = (slice[i1]-slice[i2])/((x_1D[i1]-x_1D[i2]))
    interp1 = x_1D[i1]-slice[i1]/a
    # print('x1',x_1D[i1],x_1D[i2],interp1)

    a = (slice[i3]-slice[i4])/((x_1D[i3]-x_1D[i4]))
    interp2 = x_1D[i3]-slice[i3]/a
    # print('x1',x_1D[i3],x_1D[i4],interp2)

    radius = abs(interp2-interp1)/2

    return radius

     

def plot_zoom(
    # file,
    data,
    # key,
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
    mode='close',
    fig1=None,
    ax2=None,
    cbar=None,
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
        key_LS_wall = "levelset_p_wall"
        x_1D = xp
        y_1D = yp

    if 'mesh_macro' in figpar.keys():
   
        exec(figpar['mesh_macro'])
        x_1D = x_1D_2
        y_1D = y_1D_2




    # try:
    #     data_1D = file[key][:]
    # except:
    #     print(colored('Failed to open '+key+' in '+figpar['file'],'red'))
    #     print(file.keys())
        

    file_name = figpar['file']

    # print(key,"max ",np.max(data_1D))

    cmap = plt.get_cmap(plotpar["cmap"])

    if figpar['zoom_mode'] == 'coord':

        i0=0
        i1=0
        j0=0
        j1=0

        if figpar["zoom"][0]>figpar["zoom"][1]:
            print('error zoom')


        for i,x in enumerate(x_1D):
            if figpar["zoom"][0][0]<x:
                break
            i0=i

        for i,x in enumerate(x_1D):
            # print('test zoom',figpar["zoom"][0][1],x,(figpar["zoom"][0][1]<x))
            i1=i
            if figpar["zoom"][0][1]<x:
                break

        for j,y in enumerate(y_1D):
            if figpar["zoom"][1][0]<y:
                break
            j0=j

        for j,y in enumerate(y_1D):
            if figpar["zoom"][1][1]<y:
                j1=j
                break

        print('zoom',figpar["zoom"],i0,i1,j0,j1,x_1D[i0],x_1D[i1],y_1D[j0],y_1D[j1])
        x_arr=x_1D[i0:i1+1]
        y_arr=y_1D[j0:j1+1]

        ii0, ii1 = i0,i1
        jj0, jj1 = j0,j1

    else:
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

    # reshape_data_veci()
    plot_bc_possible_based_on_dim = True

    # print("dim",data_1D.ndim)

    data_1D = data

    if data_1D.ndim ==1:
        # print('data_1D.ndim == 1')
        field0 = veci(data_1D,nx,ny,field_index)
    elif data_1D.shape[1] == 1:
        # print('data_1D.shape[1] == 1')
        data_1D = veci(data_1D[:,0],nx,ny,field_index)
        field = data_1D
    elif data_1D.shape[0] == 1:
        # print('data_1D.shape[0] == 1')
        # print(data_1D)
        # print(data_1D[0]) 
        #[0,:]
        # print("data_1D.shape[0] == 1")
        data_1D = data_1D[0]
        field0 = veci(data_1D,nx,ny,field_index)
    elif data_1D.ndim ==2:
        # print("data_1D.ndim ==2")
        field0 = data_1D.transpose()
        plot_bc_possible_based_on_dim = False

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

    if parse_is_true(get_value_from_dicts('plot_bc',figpar,plotpar)) and plot_bc_possible_based_on_dim:
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

        if vecb_l: 
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

    if mode == 'film' or mode == 'first':
        ax2.clear()
    else:
        fig1,ax2 = init_fig(plotpar,figpar)
   

    if 'range' in figpar.keys():
        CS = ax2.contourf(x_arr,y_arr,field, 
        # levels=figpar['range'], 
        levels=eval(figpar['range']),
        cmap=plotpar['cmap'],extend=plotpar['extend'],)
    else:

        if get_value_from_dicts('plot_mode',figpar,plotpar) == "contourf":
          

            CS = ax2.contourf(x_arr,y_arr,field, 
            levels=get_value_from_dicts('levels',figpar,plotpar), #10, 
            cmap=plotpar['cmap'],extend=plotpar['extend'],)
        else:

            mpl_levels = mticker.MaxNLocator(nbins=get_value_from_dicts('levels',figpar,plotpar)).tick_values(np.min(field), np.max(field))
            norm = mpl_colors.BoundaryNorm(mpl_levels, ncolors=cmap.N, clip=True)
            CS = ax2.pcolormesh(x_arr,y_arr,field, cmap=plotpar['cmap'], norm=norm)


    lcolor= "k" #"w" #"k"
    lw=0.5
    ms=0.5

    if parse_is_true(get_value_from_dicts('plot_grid',figpar,plotpar)):
       

        color_annot_bc = get_value_from_dicts('color_annot_bc',figpar,plotpar)
        color_annot_bulk = get_value_from_dicts('color_annot_bulk',figpar,plotpar)
    
        for igrid0 in range(i0,i1+1):        
            for jgrid0 in range(j0,j1+1): 
                
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

                if 'print_mode' in figpar.keys():
                    if figpar['print_mode'] == "val":
                        str1='{:.2e}'.format(field[jgrid,igrid])
                    elif figpar['print_mode'] == "valres":
                        str1=figpar['print_res'].format(field[jgrid,igrid])
                    elif figpar['print_mode'] == "val10":
                        str1='{:.10e}'.format(field[jgrid,igrid])
                    elif figpar['print_mode'] == "ij":
                        str1="{:03} {:03}".format(igrid0,jgrid0)
                    elif figpar['print_mode'] == "ijval": 
                        str1="{:.2e} {:03} {:03}".format(field[jgrid,igrid],igrid0,jgrid0)      
                    elif figpar['print_mode'] == "ijcoord": 
                        str1="{:.2e} {:.2e} {:03} {:03}".format(x_arr[igrid],y_arr[jgrid],igrid0,jgrid0)        
                    elif figpar['print_mode'] == "ijx": 
                        str1="{:.2e} {:03} {:03}".format(x_arr[igrid],igrid0,jgrid0)         
                    elif figpar['print_mode'] == "ijy": 
                        str1="{:.2e} {:03} {:03}".format(y_arr[jgrid],igrid0,jgrid0)  
                else:
                    str1='{:.2e}'.format(field[jgrid,igrid])

                if 'fontsize' in figpar.keys():
                    fontsize = get_value_from_dicts('fontsize',figpar,plotpar)
                else:
                    fontsize = plotpar['fontsize']

                ax2.annotate(str1,(x_arr[igrid],y_arr[jgrid]),fontsize=fontsize,c=lcolor,ha="center",va=va)

        if 'plot_capacities' in figpar.keys():
            for igrid0 in range(i0,i1+1):        
                for jgrid0 in range(j0,j1+1): 
                   
                    igrid =igrid0-i0 #TODO
                    jgrid = jgrid0-j0

                    try:
                        xc = (x_arr[igrid]+x_arr[igrid-1])/2
                        yc = (y_arr[jgrid]+y_arr[jgrid-1])/2

                        xc2 = (x_arr[igrid]+x_arr[igrid+1])/2
                        yc2 = (y_arr[jgrid]+y_arr[jgrid+1])/2
                    
                        dcap_1 = file['dcap_1'][:].transpose()
                        dcap_2 = file['dcap_2'][:].transpose()
                        dcap_3 = file['dcap_3'][:].transpose()
                        dcap_4 = file['dcap_4'][:].transpose()

                        if [igrid,jgrid] in figpar['plot_capacities_ijlist']:
                            print(igrid,jgrid,field[jgrid,igrid]/plotpar["scale_y"],dcap_1[jgrid,igrid]/plotpar["scale_y"],x_arr[igrid],y_arr[jgrid])


                            
                            ax2.plot([xc,xc],[yc,yc+dcap_1[jgrid,igrid]/plotpar["scale_y"]],color="pink",lw=lw)

                            ax2.plot([xc2,xc2],[yc,yc+dcap_3[jgrid,igrid]/plotpar["scale_y"]],color="green",lw=lw)

                            ax2.plot([xc,xc+dcap_2[jgrid,igrid]/plotpar["scale_x"]],[yc,yc],color='black',lw=lw)

                            ax2.plot([xc,xc+dcap_4[jgrid,igrid]/plotpar["scale_x"]],[yc2,yc2],color='blue',lw=lw)

                    except:
                        print('not plotted')
   

    ax2.set_title('Time '+r"$\SI[retain-zero-exponent=true]{{{0:.2e}}}".format(time/plotpar['scale_time'])+'{'+plotpar['unit_time']+'}$',color=plotpar['text_color'])


   
    # # Make a colorbar for the ContourSet returned by the contourf call.
    # cbar = fig1.colorbar(CS)
    # cbar.ax.set_ylabel(r""+figpar['cbarlabel'])

    if mode !='film':
        cbar = fig1.colorbar(CS)
        cbar.ax.set_ylabel(r""+figpar['cbarlabel'],color=plotpar['text_color'])
    # Add the contour line levels to the colorbar

    else:
        cbar = plt.colorbar(CS,cax=cbar.ax)
        cbar.ax.set_ylabel(r""+figpar['cbarlabel'],color=plotpar['text_color'])
        # if 'ticks_format' in figpar:
        if get_value_from_dicts('ticks_format',figpar,plotpar)!=None:
            cbar.ax.yaxis.set_major_formatter(mticker.FormatStrFormatter(get_value_from_dicts('ticks_format',figpar,plotpar)))

    # Add the contour line levels to the colorbar
    if str(get_value_from_dicts('isocontour',figpar,plotpar)) == 'True':
        CS2 = ax2.contour(CS, 
        # levels=CS.levels[::2], 
        # levels=
        colors="r")
        cbar.add_lines(CS2)

    if get_value_from_dicts('plot_levelset',figpar,plotpar):
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

        # print("test ii0 ",ii0,ii1+1,jj0,jj1+1)
        try:
            linewidths = get_value_from_dicts('linewidth',figpar,plotpar)
            linestyles=get_value_from_dicts('linestyle',figpar,plotpar)
        except:
            linewidths = plotpar['linewidth']
            linestyles=plotpar['linestyle']

        CSlvl = ax2.contour(
            x_1D[ii0:ii1+1], y_1D[jj0:jj1+1], LSdat[jj0:jj1+1, ii0:ii1+1], [0.0], colors="r",
            linewidths=linewidths,linestyles=linestyles,
        )


    if 'plot_wall' in figpar.keys():
        if figpar['plot_wall']:

            wallii0 = ii0
            wallii1 = ii1
            walljj0 = jj0
            walljj1 = jj1


            try:
                LSdat = file[key_LS_wall][:]

                LSdat = LSdat.transpose()

                CSlvl = ax2.contour(
                x_1D[wallii0:wallii1+1], y_1D[walljj0:walljj1+1], LSdat[walljj0:walljj1+1, wallii0:wallii1+1], [0.0], 
                colors="orange",linewidths=get_value_from_dicts('linewidth',figpar,plotpar),linestyles=get_value_from_dicts('linestyle',figpar,plotpar),
                clip_on=True,
                )

                CSlvlwall = ax2.contourf(x_1D[wallii0:wallii1+1], y_1D[walljj0:walljj1+1], LSdat[walljj0:walljj1+1, wallii0:wallii1+1], 
                            levels=0,
                            # clip_on=True,
                            ) #not very precise

                # # print(CSlvlwall.levels)

                # # A = plt.Polygon(np.array([(0, 0), (50, 100), (100, 0)]), color='w', ec='k')
                # # B = plt.Polygon(np.array([(120, 0), (170, 100), (220, 0)]), color='w', ec='k')
                # # C = plt.Polygon(np.array([(240, 0), (290, 100), (340, 0)]), color='w', ec='k')
            
                # A = plt.Polygon(np.array([(figpar["zoom"][0][0], figpar["zoom"][1][0]), 
                #                           (figpar["zoom"][0][1], figpar["zoom"][1][0]),
                #                           (figpar["zoom"][0][1], figpar["zoom"][1][1]), 
                #                           (figpar["zoom"][0][0], figpar["zoom"][1][1])
                #                         ]), color='w', ec='k')


                # # fig, ax = plt.subplots()
                # # all_polys = [A, B, C]
                # all_polys= [A]
                # [ax2.add_patch(i) for i in all_polys]
                # vertices = np.concatenate([i.get_path().vertices for i in all_polys])
                # codes = np.concatenate([i.get_path().codes for i in all_polys])

                # # dots = ax2.scatter(points[:, 0], points[:, 1], zorder=3)
                # # CSlvlwall.set_clip_path(PathPatch(Path(vertices, codes), transform=ax2.transData))
                # # plt.show()


                for level, collection in zip(CSlvlwall.levels[1:], CSlvlwall.collections):
                    # print(f"Levelset contours",level)
                    if level>0:
                        collection.remove()
                    else:
                        collection.set_facecolor(plotpar['color_wall'])  

                    # collection.set_clip_path(clip) 
            
            except:
                file_wall_name = "flower_00000000.h5"
                # print('loading file to plot wall: ',file_wall_name)
                with h5py.File(file_wall_name, "r") as file_wall:
                    LSdat = file_wall[key_LS_wall][:]

                    LSdat = LSdat.transpose()

                    CSlvl = ax2.contour(
                    x_1D[wallii0:wallii1+1], y_1D[walljj0:walljj1+1], LSdat[walljj0:walljj1+1, wallii0:wallii1+1], [0.0], colors="orange",linewidths=get_value_from_dicts('linewidth',figpar,plotpar),linestyles=get_value_from_dicts('linestyle',figpar,plotpar)
                    )
                    
                    # CSlvlwall = ax2.contourf(x_1D, y_1D, LSdat, levels=1,colors='gray') #not very precise
                    CSlvlwall = ax2.contourf(x_1D[wallii0:wallii1+1], y_1D[walljj0:walljj1+1], LSdat[walljj0:walljj1+1, wallii0:wallii1+1], levels=0,
                                            #  zorder=3,
                                             ) #not very precise

                    # print(CSlvlwall.levels)

                    if "clip" in figpar.keys():
                        dxpatch = mesh["dx"]/2/float(plotpar["scale_x"])
                        dypatch = mesh["dy"]/2/float(plotpar["scale_x"])

                        print("patch ", dxpatch,dypatch)


                        A = plt.Polygon(np.array([(figpar["zoom"][0][0]-dxpatch, figpar["zoom"][1][0]-dypatch), 
                                                    (figpar["zoom"][0][1]+dxpatch, figpar["zoom"][1][0]-dypatch),
                                                    (figpar["zoom"][0][1]+dxpatch, figpar["zoom"][1][1]+dypatch), 
                                                    (figpar["zoom"][0][0]-dxpatch, figpar["zoom"][1][1]+dypatch)
                                                ]), color='w', ec='k',alpha=0)


                        # fig, ax = plt.subplots()
                        # all_polys = [A, B, C]
                        all_polys= [A]
                        [ax2.add_patch(i) for i in all_polys]
                        vertices = np.concatenate([i.get_path().vertices for i in all_polys])
                        codes = np.concatenate([i.get_path().codes for i in all_polys])

                    # dots = ax2.scatter(points[:, 0], points[:, 1], zorder=3)
                    # CSlvlwall.set_clip_path(PathPatch(Path(vertices, codes), transform=ax2.transData))
                    # plt.show()

                    for level, collection in zip(CSlvlwall.levels[1:], CSlvlwall.collections):
                        # print(f"Levelset contours",level)
                        if level>0:
                            collection.remove()
                        else:
                            collection.set_facecolor(plotpar['color_wall'])  
                        
                        if "clip" in figpar.keys():
                            collection.set_clip_path(PathPatch(Path(vertices, codes), transform=ax2.transData))

    if get_value_from_dicts('plot_levelset_segments',figpar,plotpar):
        ax2 = plot_segments(file,plotpar,figpar,ax2)


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

    if 'aspect_ratio' in figpar.keys():
        ax2.set_aspect(aspect=get_value_from_dicts('aspect_ratio',figpar,plotpar),adjustable=get_value_from_dicts('aspect_box',figpar,plotpar))

    #debug subplots with  plot_children(fig1)

    if mode =='first' or mode =='close':
        ax2.spines["right"].set_visible(False)
        ax2.spines["top"].set_visible(False)

        ax2.set_xlabel(r""+plotpar['xlabel'],color=plotpar['text_color'])
        ax2.set_ylabel(r""+plotpar['ylabel'],color=plotpar['text_color'])

        # ax2.set_xlim([float(x0) for x0 in get_value_from_dicts('xlim',figpar,plotpar)])
        # ax2.set_ylim([float(x0) for x0 in get_value_from_dicts('ylim',figpar,plotpar)])
        # ax2.set_aspect('equal', 'box')

        str_nstep = str(nstep)
        # plt.savefig(file_name+'_'+str_nstep+ "." + plotpar["img_format"],dpi=plotpar['dpi']) #also for film for latex display


        for macro in get_value_from_dicts('macro_file_name',figpar,plotpar):
            # print(macro)
            # print(colored(file_name+"_"+(mesh["nx"])+"_"+("VOF" if yml["flower"]["simulation"]["surface_tension"] == 0 else "LS")+"_"+plotpar["theme"]+ ".pdf","red"))
            # print(colored(file_name+"_"+str(mesh["nx"])+"_"+("VOF" if yml["flower"]["simulation"]["surface_tension"] == 0 else "LS")+"_"+plotpar["theme"]+ ".pdf","red"))
            
            plt.savefig(eval(macro),dpi=plotpar['dpi'],transparent=True)
            
            if 'svg' in macro:
                gen_name = eval(macro).split('.')[0]
                #print(gen_name)
                call_inkscape(gen_name)


    if mode == 'close':
        # str_nstep = str(nstep)
        # plt.savefig(file_name+'_'+str_nstep+ "." + plotpar["img_format"],dpi=plotpar['dpi'])
        plt.close(fig1)
        return
    
    return(fig1,ax2,cbar)
