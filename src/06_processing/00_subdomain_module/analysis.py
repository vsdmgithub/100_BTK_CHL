# Functions for analysis
import numpy as np
from numba import jit, prange
from grid_operations import get_staggered_grid
from scipy.interpolate import interpn
from numpy.linalg import norm as l2

def get_rms_error(f1,f2):
    """
    Compute the root mean square error between two arrays.

    Parameters:
        f1 (ndarray): The first array.
        f2 (ndarray): The second array.

    Returns:
        float: The root mean square error between the two arrays.
    """
    
    dim = f1.size
    if f1.shape != f2.shape:
        raise ValueError("Arrays must have the same shape")
        exit
    return l2(f1-f2)/np.sqrt(dim) 

def get_collocated_velocities(uf, vf, wf, xx, yy, zz, interp_method = 'linear'):
    """
    Interpolates the given velocity fields `uf`, `vf`, and `wf` onto the cell-centered grid defined by `xx`, `yy`, and `zz`.

    Parameters:
        uf (ndarray): The velocity field in the x-direction.
        vf (ndarray): The velocity field in the y-direction.
        wf (ndarray): The velocity field in the z-direction.
        xx (ndarray): The x-coordinates of the grid points.
        yy (ndarray): The y-coordinates of the grid points.
        zz (ndarray): The z-coordinates of the grid points.

    Returns:
        tuple: A tuple containing the interpolated velocity fields `uc`, `vc`, and `wc` on the cell-centered grid.
    """

    # * Obtain grids for each component
    xu, yu, zu = get_staggered_grid(xx, yy, zz, 1)
    xv, yv, zv = get_staggered_grid(xx, yy, zz, 2)
    xw, yw, zw = get_staggered_grid(xx, yy, zz, 3)

    Z, X, Y =np.meshgrid(zz, xx, yy, indexing='ij')

    # * Flatten the new grid points
    pts = np.array([Z.flatten(), X.flatten(), Y.flatten()]).T

    # * Interpolate the values at the new grid points
    pts_u=(zu[:],xu[:],yu[:])
    pts_v=(zv[:],xv[:],yv[:])
    pts_w=(zw[:],xw[:],yw[:])

    uc = interpn(pts_u, uf , pts, method=interp_method, bounds_error=False, fill_value=None)
    vc = interpn(pts_v, vf , pts, method=interp_method, bounds_error=False, fill_value=None)
    wc = interpn(pts_w, wf , pts, method=interp_method, bounds_error=False, fill_value=None)

    uc = uc.reshape(X.shape)
    vc = vc.reshape(X.shape)
    wc = wc.reshape(X.shape)

    return uc, vc, wc

def get_collocated_pressure(pf, xx, yy, zz, interp_method = 'linear'):

    # * Obtain grids for each component
    xp, yp, zp = get_staggered_grid(xx, yy, zz, 4, flag_interior=True) # * Get the pressure inside the domain

    Z, X, Y = np.meshgrid(zz, xx, yy, indexing='ij')

    # * Flatten the new grid points
    pts = np.array([Z.flatten(), X.flatten(), Y.flatten()]).T

    # * Interpolate the values at the new grid points
    pts_p = (zp[:], xp[:], yp[:])

    pf = pf[1:-1, 1:-1, 1:-1] # * Remove the ghost cells
    pc = interpn(pts_p, pf, pts, method=interp_method, bounds_error=False, fill_value=None)

    pc = pc.reshape(X.shape)
    
    return pc

def get_collocated_velocities_and_pressure(uf, vf, wf, pf, xx, yy, zz, interp_method = 'linear'):
    """
    Interpolates the given velocity fields `uf`, `vf`, and `wf` onto the cell-centered grid defined by `xx`, `yy`, and `zz`.

    Parameters:
        uf (ndarray): The velocity field in the x-direction.
        vf (ndarray): The velocity field in the y-direction.
        wf (ndarray): The velocity field in the z-direction.
        xx (ndarray): The x-coordinates of the grid points.
        yy (ndarray): The y-coordinates of the grid points.
        zz (ndarray): The z-coordinates of the grid points.
        pf (ndarray): The pressure field.
        interp_method (str): The interpolation method to use.

    Returns:
        tuple: A tuple containing the interpolated velocity fields `uc`, `vc`, and `wc` and the interpolated pressure field `pc` on the cell-centered grid.
    """

    # * Obtain grids for each component
    xu, yu, zu = get_staggered_grid(xx, yy, zz, 1)
    xv, yv, zv = get_staggered_grid(xx, yy, zz, 2)
    xw, yw, zw = get_staggered_grid(xx, yy, zz, 3)
    xp, yp, zp = get_staggered_grid(xx, yy, zz, 4, flag_interior=True) # * Get the pressure inside the domain

    Z, X, Y = np.meshgrid(zz, xx, yy, indexing='ij')

    # * Flatten the new grid points
    pts = np.array([Z.flatten(), X.flatten(), Y.flatten()]).T

    # * Interpolate the values at the new grid points
    pts_u = (zu[:], xu[:], yu[:])
    pts_v = (zv[:], xv[:], yv[:])
    pts_w = (zw[:], xw[:], yw[:])
    pts_p = (zp[:], xp[:], yp[:])

    uc = interpn(pts_u, uf, pts, method=interp_method, bounds_error=False, fill_value=None)
    vc = interpn(pts_v, vf, pts, method=interp_method, bounds_error=False, fill_value=None)
    wc = interpn(pts_w, wf, pts, method=interp_method, bounds_error=False, fill_value=None)
    pf = pf[1:-1, 1:-1, 1:-1] # * Remove the ghost cells
    pc = interpn(pts_p, pf, pts, method=interp_method, bounds_error=False, fill_value=None)

    uc = uc.reshape(X.shape)
    vc = vc.reshape(X.shape)
    wc = wc.reshape(X.shape)
    pc = pc.reshape(X.shape)
    
    return uc, vc, wc, pc

def stagger_the_grid(rr, flag):
    from grid_operations import get_grid_wHalos, get_grid_center_wHalos
    if flag == 0:
        ru = get_grid_wHalos(rr)
        rv = get_grid_center_wHalos(rr)
        rw = get_grid_center_wHalos(rr)
    elif flag == 1:
        ru = get_grid_center_wHalos(rr)
        rv = get_grid_wHalos(rr)
        rw = get_grid_center_wHalos(rr)
    elif flag == 2:
        ru = get_grid_center_wHalos(rr)
        rv = get_grid_center_wHalos(rr)
        rw = get_grid_wHalos(rr)
    else:
        ru = get_grid_center_wHalos(rr)
        rv = get_grid_center_wHalos(rr)
        rw = get_grid_center_wHalos(rr)

    return ru, rv, rw

def get_hor_ver_grid(xx, yy, zz, flag):
    if flag == 1 or flag == 2:
        horGrid = xx; verGrid = zz
        stag_flag = [0,2]
    elif flag == 3 or flag == 4:
        horGrid = zz; verGrid = yy
        stag_flag = [2,1]
    elif flag == 5 or flag == 6:
        horGrid = xx; verGrid = yy
        stag_flag = [0,1]
    else:
        print('Invalid flag')
        exit()
    return horGrid, verGrid, stag_flag

def get_collocated_velocities_at_boundary(u, v, w, xx, yy, zz, bcflag=-1, interp_method = 'linear'):
    # * CREATE THE MESHGRID FOR THE FACE
    """
    Interpolates the given velocity fields `u`, `v`, and `w` onto collocation points on a face of a 3D grid.

    Parameters:
        u (ndarray): The velocity field in the x-direction.
        v (ndarray): The velocity field in the y-direction.
        w (ndarray): The velocity field in the z-direction.
        xx (ndarray): The x-coordinates of the grid points.
        yy (ndarray): The y-coordinates of the grid points.
        zz (ndarray): The z-coordinates of the grid points.
        flag (int): The face for which to perform the interpolation. The faces are numbered as:
            1: Positive x-face
            2: Negative x-face
            3: Positive y-face
            4: Negative y-face
            5: Positive z-face
            6: Negative z-face
        interp_method (str): The interpolation method. Default is 'linear'.

    Returns:
        list: A list containing the interpolated velocity fields `uc`, `vc`, and `wc` on the collocated grid.
    """
    horGrid, verGrid, stag_flag = get_hor_ver_grid(xx, yy, zz, bcflag)
    hor = stagger_the_grid(horGrid, stag_flag[0])
    ver = stagger_the_grid(verGrid, stag_flag[1])
    
        # * Create data points and query points in a more concise manner
    data_pts = [(np.array(ver[i]), np.array(hor[i])) for i in range(3)]
    verMesh, horMesh = np.meshgrid(verGrid, horGrid, indexing='ij')
    query_pts = np.array([verMesh.flatten(), horMesh.flatten()]).T

        # * Interpolation to collocation points using a loop
    interp_results = [
        interpn((data_pts[i][0][:], data_pts[i][1][:]), [u, v, w][i], query_pts,
                method=interp_method, bounds_error=False, fill_value=None).reshape(len(verGrid), len(horGrid)) 
        for i in range(3) ]
    return interp_results

def get_wallStress(u, w, yy,Re_c):
   """
   Computes the wall shear stress from the velocity fields `u` and `w`

   Parameters:
       u (ndarray): The velocity field in the x-direction.
       w (ndarray): The velocity field in the z-direction.
       xx (ndarray): The x-coordinates of the grid points.
       yy (ndarray): The y-coordinates of the grid points.
       zz (ndarray): The z-coordinates of the grid points.
       visc (float): The kinematic viscosity of the fluid.

   Returns:
       tuple: A tuple containing the wall shear stress in the x- and z-directions, `Sxy` and `Szy`.

   """
   Sxy = np.sqrt(np.abs(u[:,:,1]-u[:,:,0])*Re_c/yy[1])
   Szy = np.sqrt(np.abs(w[:,:,1]-w[:,:,0])*Re_c/yy[1])

   return Sxy[1:,1:-1], Szy[1:-1,1:]

def get_wallData(fileroot,ngx,ngy,ngz,yy,num,Re_c,clip_halos=True):
    """
    Reads U,W velocity data and pressure data from binary fortran files 
    and computes the friction Reynolds number, the wall-normal gradient of 
    U and W, and the wall pressure.

    :param fileroot:  Location of data files
    :param ng[x,y,z]: number of grid points in each dimension
    :param yy:        Wall-normal coordinates
    :param num:       Time step to extract data
    :param Re_c:      Centerline Reynolds number

    :return: The friction Reynolds number, the wall-normal gradient of 
             U and W, and the wall pressure as numpy arrays 
    """
    fname    = '01uc00.dble.'
    filename = fileroot + fname + "{:07d}".format(num)
    Ufield   = np.fromfile(filename,dtype=float,count=(ngx+2)*(ngz+1)*(ngy+1),offset=0)
    Ufield   = np.transpose(np.reshape(Ufield,(ngy+1,ngx+2,ngz+1))) 

    fname    = '03wc00.dble.'
    filename = fileroot + fname + "{:07d}".format(num)
    Wfield   = np.fromfile(filename,dtype=float,count=(ngx+1)*(ngz+2)*(ngy+1),offset=0)
    Wfield   = np.transpose(np.reshape(Wfield,(ngy+1,ngx+1,ngz+2))) 

    # fname    = 'pres_p.dble.'
    # filename = fileroot + fname + "{:07d}".format(num)
    # Pfield   = np.fromfile(filename,dtype=float,count=(ngx+1)*(ngz+1)*(ngy+1),offset=0)
    # Pfield   = np.transpose(np.reshape(Pfield,(ngy+1,ngx+1,ngz+1))) 

    if clip_halos:
        #^ Clip the halos off
        Ufield = Ufield[1:-1,1:-1,:]
        Wfield = Wfield[1:-1,1:-1,:]

    du_dy    = (Ufield[:,:,1]-Ufield[:,:,0])/yy[1]
    dw_dy    = (Wfield[:,:,1]-Wfield[:,:,0])/yy[1]
    # Re_tau   = np.sqrt(np.abs(du_dy)*Re_c)
    # pWall    = (Pfield[1:-1,1:-1,1] + Pfield[1:-1,1:-1,0])/2.0

    # return (Re_tau, du_dy/Re_c, dw_dy/Re_c, pWall) 
    return du_dy/Re_c, dw_dy/Re_c

def decompose_pressure(Pn, xp, ngx, ngy, ngz, gradx_Pi, remove_gradx_Pi=False):
    """
    Decompose the pressure field into fluctuating and mean components.

    Parameters
    ----------
    Pn : numpy.ndarray
        The input pressure field.
    xp : numpy.ndarray
        The x-coordinates of the points in the domain.
    ngx : int
        Number of grid points in the x-direction.
    ngy : int
        Number of grid points in the y-direction.
    ngz : int
        Number of grid points in the z-direction.
    gradx_Pi : float
        The gradient of the mean pressure field.

    Returns
    -------
    Pn_fluc : numpy.ndarray
        The fluctuating component of the pressure field.
    net_mean_gradx_Ptot : float
        The net mean x-gradient of the total pressure field.
    """

    Gg                      = np.mean(Pn)
    Pn                      = Pn - Gg
    mean_gradx_Pn           = np.mean(calculate_x_gradient(Pn, xp[1] - xp[0]))
    net_mean_gradx_Ptot     = -gradx_Pi + mean_gradx_Pn
    if remove_gradx_Pi:
        Pn_fluc                 = Pn - get_linear_field(-gradx_Pi,xp,ngx,ngy,ngz) 
    else:
        Pn_fluc                 = Pn - get_linear_field(mean_gradx_Pn,xp,ngx,ngy,ngz) 

    return Pn_fluc, mean_gradx_Pn 

def get_linear_field(grad,X_ar,nx,ny,nz):
    """
    Create a linear field of the form lin_field[:,:,i] = grad * X_ar[i].
    
    Parameters
    ----------
    grad : float
        The gradient of the linear field.
    X_ar : array
        A 1D array of size nx, containing the x-coordinates of the points.
    nx : int
        The number of points in the x-direction.
    ny : int
        The number of points in the y-direction.
    nz : int
        The number of points in the z-direction.
    
    Returns
    -------
    lin_field : array
        A 3D array of size (nz,ny,nx) containing the linear field.
    """
    lin_field = np.zeros((nz,nx,ny))
    if (X_ar.shape[0] != nx):
        print("Dimension mismatch! size(X_ar) != nx")
    for i in range(0,nx): 
        lin_field[:,i,:] = grad * X_ar[i] 
    mean      = np.mean(lin_field)
    lin_field = lin_field - mean # Remove the mean

    return lin_field

def calculate_x_gradient(array, dx):
    """
    Calculate the x-gradient of a 3D array while keeping the same shape.
    
    Parameters:
        array (numpy.ndarray): The input 3D array with shape (nz, nx, ny).
        dx (float): The grid spacing in the x-direction.
    
    Returns:
        numpy.ndarray: The x-gradient of the input array with the same shape.
    """
    grad_x = np.zeros_like(array)

    # Central difference for interior points
    grad_x[:, 1:-1, :] = (array[:, 2:, :] - array[:, :-2, :]) / (2 * dx)

    # Forward difference for the first point
    grad_x[:, 0, :] = (array[:,  1, :] - array[:, 0, :]) / dx

    # Backward difference for the last point
    grad_x[:, -1, :] = (array[:, -1, :] - array[:, -2, :]) / dx

    return grad_x

def extend_grid(arr, N_ext_start, N_ext_end):
    """
    Extend an array by symmetrically mirroring the ends.

    Parameters:
    arr (numpy.ndarray): The input array to be extended.
    N_ext_start (int): The number of points to extend at the start of the array.
    N_ext_end (int): The number of points to extend at the end of the array.

    Returns:
    numpy.ndarray: The extended array.
    """
    if N_ext_start == 0 and N_ext_end == 0:
        return arr
    while N_ext_start != 0:
        arr = np.insert(arr, 0, 2*arr[0] - arr[1])
        N_ext_start = N_ext_start - 1
    while N_ext_end != 0: 
        arr = np.insert(arr, np.size(arr), 2*arr[-1] - arr[-2])
        N_ext_end = N_ext_end - 1
    return arr
     
def extend_array_2d(arr, N_ext_bot, N_ext_top, N_ext_lft, N_ext_rht):
    """
    Extend a 2D array by padding specified numbers of rows and columns.
    Edges are padded by copying boundary values. Corners are filled with 
    the mean of adjacent corner values.

    Parameters:
    arr (numpy.ndarray): 2D array to be extended.
    N_ext_top (int): Number of rows to add at the top.
    N_ext_bottom (int): Number of rows to add at the bottom.
    N_ext_left (int): Number of columns to add to the left.
    N_ext_right (int): Number of columns to add to the right.

    Returns:
    numpy.ndarray: The extended 2D array.
    """
    ny, nx = arr.shape

    # Create the extended array filled with zeros first
    ext_arr = np.zeros((ny + N_ext_top + N_ext_bot, 
                        nx + N_ext_lft + N_ext_rht), dtype=arr.dtype)

    # Insert the original array
    ext_arr[N_ext_bot:N_ext_bot+ny, N_ext_lft:N_ext_lft+nx] = arr

    # Top and bot edge rows
    if N_ext_bot > 0:
        ext_arr[:N_ext_bot, N_ext_lft:N_ext_lft+nx] = arr[0:1, :]
    if N_ext_top > 0:
        ext_arr[N_ext_bot+ny:, N_ext_lft:N_ext_lft+nx] = arr[-1:, :]

    # lft and rht edge columns
    if N_ext_lft > 0:
        ext_arr[:, :N_ext_lft] = ext_arr[:, N_ext_lft:N_ext_lft+1]
    if N_ext_rht > 0:
        ext_arr[:, N_ext_lft+nx:] = ext_arr[:, N_ext_lft+nx-1:N_ext_lft+nx]

    # Corners: use the actual values of two corners from the original array
    if N_ext_bot > 0 and N_ext_lft > 0:
        ext_arr[:N_ext_bot, :N_ext_lft] = arr[0,0]
    if N_ext_bot > 0 and N_ext_rht > 0:
        ext_arr[:N_ext_bot, -N_ext_rht:] = arr[0,-1]
    if N_ext_top > 0 and N_ext_lft > 0:
        ext_arr[-N_ext_top:, :N_ext_lft] = arr[-1,0]
    if N_ext_top > 0 and N_ext_rht > 0:
        ext_arr[-N_ext_top:, -N_ext_rht:] = arr[-1,-1]

    return ext_arr

def extend_array_2d_taper_to_zero(arr, N_ext, taper_value=0.0):

    """
    Extend a 2D array by adding tapered rows and columns that gradually reduce to zero.

    This function extends a given 2D array by adding `N_ext` rows and columns on each side.
    The values in the added rows and columns are tapered to gradually reduce to a specified
    `taper_value`, creating a smooth transition from the original array to the extended edges.

    Parameters:
    arr (numpy.ndarray): 2D array to be extended.
    N_ext (int): Number of rows and columns to add on each side of the array.
    taper_value (float): The value to which the extended regions taper. Default is 0.0.

    Returns:
    numpy.ndarray: The extended 2D array with tapered edges.
    """

    ny, nx = arr.shape
    ext_shape = (ny + 2 * N_ext, nx + 2 * N_ext)
    ext_arr = np.zeros(ext_shape, dtype=arr.dtype)

    # Insert original array
    ext_arr[N_ext:N_ext+ny, N_ext:N_ext+nx] = arr

    # Taper top and bot edge rows
    for i in range(0, N_ext):
        weight = np.cos(0.5 * np.pi * (N_ext - 1 - i) / (N_ext - 1))
        weight_C = 1.0 - weight
        ext_arr[i, N_ext:N_ext+nx] = arr[0:1, :] * weight + taper_value * weight_C
        ext_arr[N_ext+ny+N_ext-i-1, N_ext:N_ext+nx] = arr[-1:, :] * weight + taper_value * weight_C

    # Taper lft and rht edge columns
    for i in range(0, N_ext):
        weight = np.cos(0.5 * np.pi * (N_ext - 1 - i) / (N_ext - 1))
        weight_C = 1.0 - weight
        ext_arr[N_ext:N_ext+ny, i] = ext_arr[N_ext:N_ext+ny, N_ext] * weight + taper_value * weight_C
        ext_arr[N_ext:N_ext+ny, N_ext+nx+N_ext-i-1] = ext_arr[N_ext:N_ext+ny, N_ext+nx-1] * weight + taper_value * weight_C
    
    for i in range(0,N_ext):
        for j in range(0,N_ext):
            weight_i = np.cos(0.5 * np.pi * (N_ext - 1 - i) / (N_ext - 1))
            weight_j = np.cos(0.5 * np.pi * (N_ext - 1 - j) / (N_ext - 1))
            weight = weight_i * weight_j
            weight_C = 1.0 - weight
            ext_arr[i, j] = ext_arr[N_ext, N_ext] * weight + taper_value * weight_C
            ext_arr[i, N_ext+nx+N_ext-j-1] = ext_arr[N_ext, N_ext+nx-1] * weight + taper_value * weight_C
            ext_arr[N_ext+ny+N_ext-i-1, j] = ext_arr[N_ext+ny-1, N_ext] * weight + taper_value * weight_C
            ext_arr[N_ext+ny+N_ext-i-1, N_ext+nx+N_ext-j-1] = ext_arr[N_ext+ny-1, N_ext+nx-1] * weight + taper_value * weight_C

    return ext_arr

def get_vel_statistics(u,v,w,num,interior_flag=True):
    """
    Calculate velocity statistics.

    Parameters
    ----------
    u : 3D numpy array
        x-component of velocity
    v : 3D numpy array
        y-component of velocity
    w : 3D numpy array
        z-component of velocity
    num : int
        Number of samples
    interior_flag : bool, optional
        Whether to remove the ghost cells

    Returns
    -------
    num : int
        Number of samples
    u_mean, v_mean, w_mean : 1D numpy arrays
        Mean of x, y, and z components of velocity
    uu, vv, ww : 1D numpy arrays
        Variance of x, y, and z components of velocity
    uv : 1D numpy array
        Covariance of x and y components of velocity
    """
    u_cc = 0.5 *(u[1:-1,1:-2,1:-1] + u[1:-1,2:-1,1:-1])
    v_cc = 0.5 *(v[1:-1,1:-1,1:-2] + v[1:-1,1:-1,2:-1])
    if interior_flag:
        u = u[1:-1,2:-2,1:-1]
        v = v[1:-1,1:-1,2:-2]
        w = w[2:-2,1:-1,1:-1]
    axis   = (0,1)
    u_mean = np.mean(u,axis=axis)
    v_mean = np.mean(v,axis=axis)
    w_mean = np.mean(w,axis=axis)
    uu     = np.mean(u**2,axis=axis)
    vv     = np.mean(v**2,axis=axis)
    ww     = np.mean(w**2,axis=axis)
    uv     = np.mean(u_cc*v_cc,axis=axis)
    return num,u_mean,v_mean,w_mean,uu,vv,ww,uv

def get_p_statistics(P, num, interior_flag=True):
    """
    Calculate pressure statistics

    Parameters
    ----------
    P : 3D numpy array
        Pressure field
    num : int
        Number of samples
    interior_flag : bool, optional
        Whether to remove the ghost cells

    Returns
    -------
    num : int
        Number of samples
    pp : 1D numpy array
        Variance of pressure field
    """
    if interior_flag:
        P = P[1:-1,1:-1,1:-1] # | Remove the ghost cells
    P     = P - np.mean(P)
    pp    = np.mean(P**2, axis=(0,1))

    return num, pp

def calculate_pdf(array, bin_edges):
    # Calculate histogram
    """
    Calculate probability density function (PDF) of an array using a histogram.

    Parameters
    ----------
    array : 1D numpy array
        Input array
    bin_edges : 1D numpy array
        Bin edges for the histogram

    Returns
    -------
    pdf : 2D numpy array
        PDF of the input array
    """
   
    hist, bin_edges = np.histogram(array, bins=bin_edges, density=True)
    bin_sizes = bin_edges[1:] - bin_edges[:-1]
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    pdf_bin = hist

    # Checking the normalization
    # norm = np.sum(pdf_bin*bin_sizes)
    # print("Norm: ", norm)

    # PDF is the normalized histogram
    pdf = np.column_stack((bin_centers, pdf_bin)) # bin_centers, hist

    return pdf

def get_vel_and_p_statistics(u,v,w,P,num,interior_flag=True):
    """
    Calculate velocity and pressure statistics.

    Parameters
    ----------
    u : 3D numpy array
        x-component of velocity
    v : 3D numpy array
        y-component of velocity
    w : 3D numpy array
        z-component of velocity
    P : 3D numpy array
        Pressure
    num : int
        Number of samples
    interior_flag : bool, optional
        Whether to remove the ghost cells

    Returns
    -------
    num : int
        Number of samples
    u_mean, v_mean, w_mean : 1D numpy arrays
        Mean of x, y, and z components of velocity
    uu, vv, ww : 1D numpy arrays
        Variance of x, y, and z components of velocity
    uv : 1D numpy array
        Covariance of x and y components of velocity
    pp : 1D numpy array
        Variance of pressure
    """
    u_cc = 0.5 *(u[1:-1,1:-2,1:-1] + u[1:-1,2:-1,1:-1])
    v_cc = 0.5 *(v[1:-1,1:-1,1:-2] + v[1:-1,1:-1,2:-1])
    if interior_flag:
        u = u[1:-1,2:-2,1:-1]
        v = v[1:-1,1:-1,2:-2]
        w = w[2:-2,1:-1,1:-1]
        P = P[1:-1,1:-1,1:-1] # | Remove the ghost cells
    # P     = P - np.mean(P) # | Remove the mean
    axis   = (0,1) # | Horizontal axes
    u_mean = np.mean(u,axis=axis)
    v_mean = np.mean(v,axis=axis)
    w_mean = np.mean(w,axis=axis)
    p_mean = np.mean(P,axis=(0,1))
    uu     = np.mean(u**2,axis=axis)
    vv     = np.mean(v**2,axis=axis)
    ww     = np.mean(w**2,axis=axis)
    uv     = np.mean(u_cc*v_cc,axis=axis)
    pp     = np.mean(P**2, axis=(0,1))

    return num,u_mean,v_mean,w_mean,p_mean,uu,vv,ww,uv,pp

def finite_diff_weights(xi, x0, order=1):
    """
    Compute finite difference weights for derivative of given order at x0
    using the stencil points in xi.
    """
    n = len(xi)
    A = np.zeros((n, n))
    dx = xi - x0
    for i in range(n):
        A[i] = dx**i / np.math.factorial(i)
    b = np.zeros(n)
    b[order] = 1.0
    return np.linalg.solve(A, b)

def custom_gradient(f, x):
    """
    Compute gradient of f(x) with:
    - 4th-order accuracy for interior points
    - 2nd-order accuracy for boundaries
    """
    f = np.asarray(f, dtype=float)
    x = np.asarray(x, dtype=float)
    n = len(f)
    df = np.zeros_like(f)

    # 2nd-order forward difference at i = 0
    xi = x[:3]
    fi = f[:3]
    w = finite_diff_weights(xi, x[0], order=1)
    df[0] = np.dot(w, fi)

    # 2nd-order backward difference at i = n-1
    xi = x[-3:]
    fi = f[-3:]
    w = finite_diff_weights(xi, x[-1], order=1)
    df[-1] = np.dot(w, fi)

    # 2nd-order centered at i = 1
    xi = x[0:3]
    fi = f[0:3]
    w = finite_diff_weights(xi, x[1], order=1)
    df[1] = np.dot(w, fi)

    # 2nd-order centered at i = n-2
    xi = x[-3:]
    fi = f[-3:]
    w = finite_diff_weights(xi, x[-2], order=1)
    df[-2] = np.dot(w, fi)

    # 4th-order interior points
    for i in range(2, n - 2):
        xi = x[i - 2:i + 3]
        fi = f[i - 2:i + 3]
        w = finite_diff_weights(xi, x[i], order=1)
        df[i] = np.dot(w, fi)

    return df
    
# -============================================================================
# - NUMBA-OPTIMIZED OBSERVATION KERNEL APPLICATION
# -============================================================================
@jit(nopython=True, parallel=True, cache=True)
def apply_observation_kernels(u, v, w, u_obs, v_obs, w_obs,
                               OKU, OKV, OKW,
                               ju_obs, jv_obs, jw_obs,
                               istart, iend, iskip,
                               kstart, kend, kskip,
                               noy, del_ii, del_jj, del_kk):
    """
    Apply observation kernels to velocity fields using Numba JIT compilation.
    This function is optimized for parallel execution on CPU.

    Parameters:
    -----------
    u, v, w : velocity fields
    u_obs, v_obs, w_obs : observation arrays (modified in-place)
    OKU, OKV, OKW : observation kernels
    ju_obs, jv_obs, jw_obs : y-indices for observations
    istart, iend, iskip : x-direction grid parameters
    kstart, kend, kskip : z-direction grid parameters
    noy : number of observations in y-direction
    del_ii, del_jj, del_kk : kernel stencil half-widths
    """
    # Loop over y-direction observations (parallel)
    for jj in prange(noy):
        ctr = jj % 3

        # Compute observation locations
        j_obs_u = ju_obs[jj] + 1
        j_obs_v = jv_obs[jj] + 1
        j_obs_w = jw_obs[jj] + 1

        # Loop over x-direction
        for i_obs in range(istart+1, iend+2, iskip):
            # Loop over z-direction
            for k_obs in range(kstart+1, kend+2, kskip):

                # Initialize accumulator for this observation point
                sum_u = 0.0
                sum_v = 0.0
                sum_w = 0.0

                # Loop over kernel stencil
                for ker_jj in range(-del_jj, del_jj+2):
                    ind_ju = j_obs_u + ker_jj
                    ind_jv = j_obs_v + ker_jj
                    ind_jw = j_obs_w + ker_jj

                    for ker_kk in range(-del_kk, del_kk+2):
                        ind_ku = k_obs + ker_kk
                        ind_kv = k_obs + ker_kk
                        ind_kw = k_obs + ker_kk

                        for ker_ii in range(-del_ii, del_ii+2):
                            ind_iu = i_obs + ker_ii
                            ind_iv = i_obs + ker_ii
                            ind_iw = i_obs + ker_ii

                            # Apply kernels
                            sum_u += OKU[ctr, ker_kk+del_kk, ker_ii+del_ii, ind_ju] * u[ind_ku, ind_iu, ind_ju]
                            sum_v += OKV[ctr, ker_kk+del_kk, ker_ii+del_ii, ind_jv] * v[ind_kv, ind_iv, ind_jv]
                            sum_w += OKW[ctr, ker_kk+del_kk, ker_ii+del_ii, ind_jw] * w[ind_kw, ind_iw, ind_jw]

                # Store results
                u_obs[k_obs, i_obs, j_obs_u] = sum_u
                v_obs[k_obs, i_obs, j_obs_v] = sum_v
                w_obs[k_obs, i_obs, j_obs_w] = sum_w

# Run a test or example if executed directly
if __name__ == "__main__":
    pass

def compute_error_with_measurements(u_ob, v_ob, w_ob,
                                    u_mn, v_mn, w_mn, 
                                    nox, noy, noz, 
                                    istart, iend, iskip, 
                                    kstart, kend, kskip, 
                                    ju_obs, jv_obs, jw_obs):

    """
    Compute the error between the observation values and the kernel computed values.

    Parameters
    ----------
    u_ob : numpy.ndarray
        The kernel computed u values
    v_ob : numpy.ndarray
        The kernel computed v values
    w_ob : numpy.ndarray
        The kernel computed w values
    u_mn : numpy.ndarray
        The experimental u values
    v_mn : numpy.ndarray
        The experimental v values
    w_mn : numpy.ndarray
        The experimental w values
    nox : int
        The number of points in the x direction
    noy : int
        The number of points in the y direction
    noz : int
        The number of points in the z direction
    istart : int
        The starting index in the x direction
    iend : int
        The ending index in the x direction
    iskip : int
        The number of points to skip between each index in the x direction
    kstart : int
        The starting index in the z direction
    kend : int
        The ending index in the z direction
    kskip : int
        The number of points to skip between each index in the z direction
    ju_obs : numpy.ndarray
        The observation locations in the y direction for u
    jv_obs : numpy.ndarray
        The observation locations in the y direction for v
    jw_obs : numpy.ndarray
        The observation locations in the y direction for w

    Returns
    -------
    err_u : numpy.ndarray
        The error between the observation values and the kernel computed u values
    err_v : numpy.ndarray
        The error between the observation values and the kernel computed v values
    err_w : numpy.ndarray
        The error between the observation values and the kernel computed w values
    """
    err_u = np.zeros((noz, nox, noy))
    err_v = np.zeros((noz, nox, noy))
    err_w = np.zeros((noz, nox, noy))

    for c1 in range(noz):
        for c2 in range(nox):
            for c3 in range(noy):
                    #--> Loop over the observation points

                u_avg = u_ob[kstart+ c1*kskip+1, istart+c2*iskip+1, ju_obs[c3]+1] + \
                        u_ob[kstart+ c1*kskip+1, istart+c2*iskip+2, ju_obs[c3]+1] + \
                        u_ob[kstart+ c1*kskip+1, istart+c2*iskip+1, ju_obs[c3]+2] + \
                        u_ob[kstart+ c1*kskip+1, istart+c2*iskip+2, ju_obs[c3]+2]

                v_avg = v_ob[kstart+ c1*kskip+1, istart+c2*iskip+1, jv_obs[c3]+1] + \
                        v_ob[kstart+ c1*kskip+1, istart+c2*iskip+1, jv_obs[c3]+2]

                w_avg = w_ob[kstart+ c1*kskip+1, istart+c2*iskip+1, jw_obs[c3]+1] + \
                        w_ob[kstart+ c1*kskip+1, istart+c2*iskip+1, jw_obs[c3]+2] + \
                        w_ob[kstart+ c1*kskip+2, istart+c2*iskip+1, jw_obs[c3]+1] + \
                        w_ob[kstart+ c1*kskip+2, istart+c2*iskip+1, jw_obs[c3]+2]

                err_u[c1, c2, c3] =  (u_avg - u_mn[kstart+ c1*kskip+1, istart+c2*iskip+1, ju_obs[c3]+1])
                err_v[c1, c2, c3] =  (v_avg - v_mn[kstart+ c1*kskip+1, istart+c2*iskip+1, jv_obs[c3]+1])
                err_w[c1, c2, c3] =  (w_avg - w_mn[kstart+ c1*kskip+1, istart+c2*iskip+1, jw_obs[c3]+1])

    return err_u, err_v, err_w