from termcolor import colored

#Debug with PDI + python instead of recompiling code in Julia (for quick tests)


# def print_bc(var):
#     # vec = vecb_T(vec[isnap,:],nx,ny)
#     print(var)



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
    print(colored('[Debug]', "red")) 
    print(extract)

    print('bottom',vecb_B(data, nx, ny))


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
