import numpy as np

def read_grid_data(fname,ngx,ngy,ngz,Lx,Ly,Lz):
    """
    Reads grid.data files
    
    :param fname: Name of fortran binary file to read (path + grid.data)
    :param ng[x,y,z]: Number of grid points in [x,y,z]
    :param L[x,y,z]: Length of domain in [x,y,z]. Lx and Ly will be determined
                     based on the grid so they can be set to Lx=None and Ly=None.

    :return: A tuple containing:

            - XX, YY:  Mesh-grid objects in x and y. 
            - xx, yy, zz: 1D arrays of grid in x,y, and z
    """
    # Read in X grid
    XX = np.fromfile(fname,dtype=float,count=ngx*ngy,offset=0)
    XX = np.transpose(np.reshape(XX,(ngy,ngx))) 
    xx = XX[:,0]
    Lx = xx[-1]-xx[0]

    # Read in Y grid
    YY = np.fromfile(fname,dtype=float,count=ngx*ngy,offset=ngx*ngy*8)
    YY = np.transpose(np.reshape(YY,(ngy,ngx)))
    yy = YY[0,:]
    Ly = yy[-1]-yy[0]

    # uniform z grid 
    zz = np.linspace(0,Lz,ngz)

    return(XX,YY,xx,yy,zz)

def write_grid_data(xx, yy, fname='grid.data'):
    """
    Write grid data to a binary file.

    Args:
        xx (ndarray): 1D array of x-coordinates.
        yy (ndarray): 1D array of y-coordinates.
        fname (str): Name of the file to write the data to.

    Returns:
        None

    This function creates mesh grids from the 1D arrays of x-coordinates and y-coordinates.
    It then flattens the mesh grids and concatenates them. The concatenated data is saved to a binary file.

    """
    # Create mesh grids from the 1D arrays
    XX,YY = np.meshgrid(xx,yy)

    # Flatten the mesh grids and concatenate them
    fout = open(fname, 'wb')
    fout.write(bytes(XX.flatten()))
    fout.write(bytes(YY.flatten()))
    fout.close()

def get_grid_wHalos(rr):
    """
    Gets the grid locations and returns the face centered grid with two halos on either ends

    Parameters
    ----------
    rr : 1D array of floats
        Grid points

    Returns
    -------
    1D array of floats
        Collocated grid with two halos on either ends
    """
    return np.concatenate([[2 * rr[0] - rr[1]], rr, [2 * rr[-1] - rr[-2]]])

def get_grid_center_wHalos(rr):
    """
    Gets the grid locations and returns the center grid with two halos on either ends

    Parameters
    ----------
    rr : 1D array of floats
        Grid points

    Returns
    -------
    1D array of floats
        Center grid with two halos on either ends
    """
    r =  0.5 * (rr[1:] + rr[:-1])
    return np.concatenate([[2 * r[0] - r[1]], r, [2 * r[-1] - r[-2]]])

def get_staggered_grid(xx, yy, zz, option, flag_interior=False):
    if option == 1:  # for u
      x = get_grid_wHalos(xx)
      y = get_grid_center_wHalos(yy)
      z = get_grid_center_wHalos(zz)
    elif option == 2:  # for v
      x = get_grid_center_wHalos(xx)
      y = get_grid_wHalos(yy)
      z = get_grid_center_wHalos(zz)
    elif option == 3:  # for w
      x = get_grid_center_wHalos(xx)
      y = get_grid_center_wHalos(yy)
      z = get_grid_wHalos(zz)
    elif option == 4:  # for T
      x = get_grid_center_wHalos(xx)
      y = get_grid_center_wHalos(yy)
      z = get_grid_center_wHalos(zz)
    if flag_interior:
        x = x[1:-1]
        y = y[1:-1]
        z = z[1:-1]
    return x, y, z

def get_observation_grid(su, sv, sw, start,end,skip):
    """
    Extracts a subset of the observation grid data for the given range and step.

    Parameters
    ----------
    su : 1D array-like
        Observation data along the u direction.
    sv : 1D array-like
        Observation data along the v direction.
    sw : 1D array-like
        Observation data along the w direction.
    start : int
        The starting index of the subset.
    end : int
        The ending index of the subset.
    skip : int
        The step size between each index in the subset.

    Returns
    -------
    tuple of 1D array-like
        Subsets of the observation data (su, sv, sw) from start to end with the specified skip.
    """

    return su[start:end+1:skip], sv[start:end+1:skip], sw[start:end+1:skip]

def get_surface_areas(xx,yy,zz):
    """
    Calculate the surface areas for a 3D grid based on the provided x, y, and z coordinates.

    :param xx: 1D array of x coordinates.
    :param yy: 1D array of y coordinates.
    :param zz: 1D array of z coordinates.

    :return: A tuple containing:
        - sf_u: 3D array representing the surface area of grid faces perpendicular to the x-axis.
        - sf_v: 3D array representing the surface area of grid faces perpendicular to the y-axis.
        - sf_w: 3D array representing the surface area of grid faces perpendicular to the z-axis.
    """

    ngx = len(xx); ngy = len(yy); ngz = len(zz)
    dz   = zz[1]-zz[0]
    dx   = xx[1]-xx[0]
    dy   = np.diff(get_grid_wHalos(yy))
    sf_u = np.tile(dy*dz, (ngz+1,ngx+2,1))
    sf_v = np.tile(dx*dz, (ngz+1,ngx+1,ngy+2))
    sf_w = np.tile(dy*dx, (ngz+2,ngx+1,1))

    return sf_u, sf_v, sf_w

# ^ Function to generate a cosine stretched grid
def get_cosine_stretched_grid(Ny, r, Ly):
    """
    Generate a cosine stretched grid

    Parameters
    ----------
    Ny : int
        The number of grid points in the y direction
    r : float
        The ratio of the maximum to minimum grid spacing
    Ly : float
        The total height of the domain

    Returns
    -------
    yy : 1D array, float
        The y coordinates of the grid points

    Notes
    -----
    The grid is generated by mapping the uniform grid in the x direction
    to the y direction using the cosine function. The grid spacing is
    proportional to the cosine of the x coordinate. The grid is then
    stretched to cover the entire domain of height Ly.
    """
    pi      = np.pi
    dy_tmp  = 0.01                      #^ Initial guess for the minimum grid spacing
    coeff_B = (r + 1) / (r - 1)         #^ Coefficient B in the cosine function
    coeff_A = dy_tmp / (coeff_B - 1)    #^ Coefficient A in the cosine function

        # --> Initialize the arrays
    arg    = np.zeros((Ny - 1))  #^ Argument of the cosine function
    dy_tmp = np.zeros_like(arg)  #^ Temporary array for the grid spacing
    for j in range(Ny - 1):
            # -> Calculate the argument of the cosine function
        arg[j] = pi * (j / (Ny - 1))
            # -> Calculate the grid spacing using the cosine function
        dy_tmp[j] = coeff_A * (coeff_B - np.cos(arg[j]))

        #--> Initialize the y-coordinate array
    yy = np.zeros((Ny))
    for j in range(1, Ny):
            #--> Calculate the y-coordinate by summing the grid spacings
        yy[j] = yy[j - 1] + dy_tmp[j - 1]

        #--> Scale the grid to cover the entire domain
    yy = yy * Ly / yy[-1]

    return yy

# ^ Function to calculate the local stretching
def get_local_stretching(yy):
    """
    Calculate the local stretching of a grid.

    Parameters
    ----------
    yy : 1D array, float
        The y coordinates of the grid points.

    Returns
    -------
    local_stretching : 1D array, float
        The local stretching factors of the grid points.

    Notes
    -----
    The local stretching is calculated as the ratio of the grid spacing
    to the left and right of a grid point. The first and last grid points
    have a stretching factor of 1.0.
    """
    Ny               = len(yy)  # Number of grid points
    local_stretching = np.zeros((Ny))  # Initialize the local stretching array

    for j in range(Ny):
        if j == 0 or j == Ny - 1:
            # Set the stretching factor for the first and last points
            local_stretching[j] = 1.0
        else:
            # Calculate local stretching for interior points
            local_stretching[j] = (yy[j + 1] - yy[j]) / (yy[j] - yy[j - 1])

    return local_stretching