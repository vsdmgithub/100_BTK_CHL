# Pre-processing functions
import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import interpn
from scipy.interpolate import interp1d

def get_advection_velocity_profile(loc,y):
    """
    Interpolates the advection velocity profile to the given y locations.

    Parameters
    ----------
    loc : str
        File location of the advection velocity profile data.
    y : array_like
        y locations to interpolate the advection velocity profile to.

    Returns
    -------
    q : array_like
        Advection velocity values at the given y locations.
    """
    data = np.loadtxt(loc)
    y_loc = data[:,0]
    q_loc = data[:,1]
    spline = UnivariateSpline(y_loc,q_loc,k=3,s=0)
    return spline(y)

def get_obs_weights(tC, nOBS, weights):
    """
    Return neighboring observation indices and weights for weighted interpolation.

    Parameters:
    tC (int): Current observation index.
    nOBS (int): Total number of observations.
    nbh (int): Number of neighbors for interpolation.

    Returns:
    tuple: Neighbor indices and their weights.
    """
    # weights = [100, 15, 12, 10]  # Predefined weights
    nbh = len(weights)-1
    obs_nbh = np.arange(-nbh, nbh + 1) + tC
    obs_nbh = obs_nbh[(obs_nbh >= 0) & (obs_nbh < nOBS)]  # Clip to valid indices
    w_nbh = np.array([weights[abs(i - tC)] for i in obs_nbh])
    return obs_nbh, w_nbh

def add_halo_Neumann_BC(arr, width=1):
    """
    Add a 1D halo with Neumann boundary conditions to the input array.

    Parameters:
    arr (numpy.ndarray): 1D input array

    Returns:
    numpy.ndarray: 1D array with added halo
    """
    if width == 0:
        return arr

    # Create a new array with an additional point on each side
    halo_arr = np.zeros(len(arr) + 2)
    
    # Copy the original array into the new array
    halo_arr[1:-1] = arr
    
    # Apply Neumann boundary conditions (i.e., copy the edge values)
    halo_arr[0] = arr[1]
    halo_arr[-1] = arr[-2]
    
    return halo_arr

def add_grid_halo(x, width=1):
    """
    Adds two extra halo points on the left and right of a 1D uniform grid.

    Parameters:
    x (numpy.ndarray): Input 1D array of locations (assumed uniform spacing).

    Returns:
    numpy.ndarray: Extended array with one extra point on each side.
    """
    if width == 0:
        return x
    dx      = x[1] - x[0]  # Compute uniform spacing
    x_left  = x[0] - dx  # New leftmost point
    x_right = x[-1] + dx  # New rightmost point
    x_halo  = np.hstack(([x_left], x, [x_right])) # Concatenate new points with original array

    return x_halo

def weighted_interpolation(q_obs, xq_ob, xq, advDistance, w_nbh, obs_nbh_rel, smoothness=0.01, halo_width=0, degree=3): 
        # * Add a halo point
    xH = add_grid_halo(xq_ob, halo_width)

        # * Line interpolation
    weights = np.array([])
    x_flat  = np.array([])
    q_flat  = np.array([])
    for n, nbh in enumerate(obs_nbh_rel):
        qH = add_halo_Neumann_BC(q_obs[n, :], halo_width)
        x_flat      = np.concatenate((x_flat,xH - nbh*advDistance))
        q_flat      = np.concatenate((q_flat,qH))
        weights_nbh = w_nbh[n] * np.ones_like(xH)
        weights     = np.concatenate((weights,weights_nbh))
    sort_idx     = np.argsort(x_flat)
    spline       = UnivariateSpline(x_flat[sort_idx], q_flat[sort_idx], k=degree, s=smoothness,w=weights[sort_idx])

    return spline(xq)

def interp_2D_DNSOBS_to_DNSDNS(xq,sq_ob,sq,q,interp_method='linear',type=1):
    if type == 2:
        # Direct 2D interpolation
        Xq, Sq    = np.meshgrid(xq,sq)
        query_pts = np.array([Sq.flatten(), Xq.flatten()]).T
        data_pts  = (sq_ob[:],xq[:])
        q_interp  = interpn(data_pts, q, query_pts, method=interp_method, bounds_error=False, fill_value=None)
        q_interp    = q_interp.reshape(Xq.shape)
    
    if type == 1:
        # Line interpolation
        q_interp = np.zeros((len(sq),len(xq)))
        for i in range(len(xq)):
            spline = UnivariateSpline(sq_ob,q[:,i],k=3,s=0)
            q_interp[:,i] = spline(sq)
    
    return q_interp

def get_inflow_data(sq,xq_ob_nearInflow,c_adv,dt_obs,dt_dns,frames_obs,frames,q_data,smoothness=0.0,degree=3):
    q_inflow   = np.zeros((frames.size,sq.size))
    x0         = xq_ob_nearInflow[0]
    for ind in range(sq.size):
        tq_m = []
        tq   = []
        q_m  = []
        for tC in range(frames_obs.size):
            for i,x in enumerate(xq_ob_nearInflow): # * Taking all grid points crossed in the observation window
                adv_dist = x - c_adv[ind] * dt_obs * tC - x0
                tq_m.append(adv_dist)
                q_m.append(q_data[tC,ind,i]) # * Stacking the inflow data
        tq_m       = np.array(tq_m)
        q_m        = np.array(q_m)
        sort_idx_m = np.argsort(tq_m) # * Sort in ascending order
        for tF in range(frames.size):
            adv_dist = - c_adv[ind] * dt_dns * tF
            tq.append(adv_dist)# * Stacking the inflow DNS positions
        tq       = np.array(tq)
        sort_idx = np.argsort(tq) # * Sort in ascending order

        spline = UnivariateSpline(tq_m[sort_idx_m], q_m[sort_idx_m], k=degree, s=smoothness)
        # spline = interp1d(xq_m[sort_idx],q_m[sort_idx], kind='cubic', fill_value="extrapolate") 
        q_inflow[:,ind] = np.flip(spline(tq[sort_idx])) # * Interpolated inflow data 

    return q_inflow

def solve_advSoln_obsForced(phi,phi_in,dt,dx,rhs,c_adv,CN_flag=False,sig=0.0,weight=0.0):
    # * rhs --> The forcing term, in combined form.
    
    # ~ Explicit upwind scheme for advection term
    """
    Solve the advection equation with a forcing term using the explicit upwind scheme.

    Parameters:
    phi (np.array): The solution at the current time step
    phi_in (float): The inflow boundary condition
    dt (float): Time step
    dx (float): Grid spacing
    rhs (np.array): The forcing term
    c_adv (float): Advection velocity

    Returns:
    np.array: The solution at the next time step
    """
    Nx      = len(phi)
    cfl     = c_adv * dt/dx 
    phi_new = np.copy(phi)
    tol     = 1e-16

    # for i in range(1, Nx - 1):  # ~ Upwind for positive ubar
    #     phi_new[i] = phi[i] - cfl * (phi[i] - phi[i - 1])
    # !!! WARNING: This approach produces numerical diffusion errors

    # ~ Explicit upwind scheme for convection term
        # & Second-order scheme for i=1 and i=Nx-2
    phi_new[1] = phi[1] - cfl * 0.5 * (phi[2] - phi[0]) + (cfl**2 / 2) * (phi[2] - 2 * phi[1] + phi[0])
    phi_new[Nx-2] = phi[Nx-2] - cfl * 0.5 * (phi[Nx-1] - phi[Nx-3]) + (cfl**2 / 2) * (phi[Nx-1] - 2 * phi[Nx-2] + phi[Nx-3])
    
        # & Third-order scheme for rest
    for i in range(2, Nx - 2):
        phi_new[i] = (phi[i] - (cfl / 6) * (2 * phi[i+1] + 3 * phi[i] - 6 * phi[i-1] + phi[i-2]) 
                      + (cfl**2 / 2) * (phi[i+1] - 2 * phi[i] + phi[i-1])
                      - (cfl**3 / 6) * (phi[i+2] - 3 * phi[i+1] + 3 * phi[i] - phi[i-1])) 
    
    if CN_flag:
        # ~ Crank-Nicolson scheme for forcing term 
        phi_new[1:] = (2*phi_new[1:] + dt * (rhs[1:] + weight[1:]) )/( 2.0 + dt * sig)
    else:
        # ~ Explicit update for the forcing term (Forward Euler)
        phi_new[1:] = phi_new[1:] + dt * rhs[1:]
    
    # ~ Apply boundary conditions
    phi_new[0]  = phi_in                             # & Inflow boundary condition
    phi_new[-1] = phi[-1]- cfl * (phi[-1] - phi[-2]) # * First order
    # phi_new[-1] = phi[-1] - cfl * (1.5 * phi[-1] - 2.0 * phi[-2] + 0.5 * phi[-3])  # * Second-order
    # phi_new[-1] = phi[-1] - cfl * (11/6 * phi[-1] - 3 * phi[-2] + 1.5 * phi[-3] - 1/3 * phi[-4])  # * Third-order

    
    return phi_new

def gaussian_filter(t, t0,norm_flag=0):
    if norm_flag==1: # * Normalized
        return np.exp(-((t - 0) ** 2) / (2 * t0 ** 2)) / np.sqrt(2 * np.pi * t0 ** 2) 
    else: # * Not normalized
        return np.exp(-((t - 0) ** 2) / (2 * t0 ** 2))

def get_WI_advected(q_WI, t_dns, t_obs, obsNt, obs_nbh, x, c_adv):
    """
    Compute the advected weighted interpolation (WI) data for neighboring observations.

    Parameters:
    ----------
    q_WI : numpy.ndarray
        Observation data for weighted interpolation (shape: N_nbh x nx).
    t_dns : numpy.ndarray
        DNS time frames (fine resolution).
    t_obs : numpy.ndarray
        Observation time frames (coarse resolution).
    obsNt : int
        Number of DNS frames per observation frame.
    obs_nbh : list
        Indices of neighboring observations.
    x : numpy.ndarray
        Spatial grid points.
    c_adv : float
        Advection velocity.

    Returns:
    -------
    q_WI_adv : numpy.ndarray
        Advected weighted interpolation data (shape: N_nbh x (N_nbh * obsNt + 1) x nx).
    """
    # Initialize variables
    dx = x[1] - x[0]  # Grid spacing
    nx = len(x)  # Number of spatial grid points
    N_nbh = len(obs_nbh)  # Number of neighboring observations
    dt_obs = t_obs[1] - t_obs[0]  # Time step between observation frames
    tol = 1e-10  # Tolerance for extrapolation warnings

    # Initialize the output array for advected data
    q_WI_adv = np.zeros((N_nbh, N_nbh * obsNt + 1, nx))

    # Loop over neighboring observations
    for tC, frame_C in enumerate(obs_nbh):
        # Initialize current data and position for the neighbor
        cur_data = q_WI[tC, :]
        cur_pos = x.copy()

        # Pad data from past observations
        for s in range(tC - 1, -1, -1):
            cur_data = np.append(cur_data, q_WI[s, -1])  # Append the last value of past observation
            cur_pos = np.append(cur_pos, cur_pos[-1] + c_adv * dt_obs)  # Append advected position
            if not is_increasing(cur_pos):  # Ensure positions are strictly increasing
                raise ValueError(f"Error: Non-increasing positions detected (tC={tC}, s={s}, cur_pos={cur_pos[-3:]})")

        # Pad data from future observations
        for s in range(tC + 1, N_nbh):
            cur_data = np.insert(cur_data, 0, q_WI[s, 0])  # Prepend the first value of future observation
            cur_pos = np.insert(cur_pos, 0, cur_pos[0] - c_adv * dt_obs)  # Prepend advected position
            if not is_increasing(cur_pos):  # Ensure positions are strictly increasing
                raise ValueError(f"Error: Non-increasing positions detected (tC={tC}, s={s}, cur_pos={cur_pos[:3]})")

        # Compute relative time for the neighbor over fine frames
        t_loc = t_dns[obs_nbh[0] * obsNt : obs_nbh[-1] * obsNt + 1] - t_obs[frame_C]

        # Loop over fine frames to compute advected data
        for tF, frame in enumerate(t_loc):
            # Compute advected positions for the current fine frame
            x_adv = cur_pos + c_adv * frame

            # Interpolate the advected data using cubic interpolation
            spline = interp1d(x_adv, cur_data, kind="cubic", fill_value="extrapolate")

            # Check for extrapolation warnings
            if x_adv[0] > x[0] + tol:
                print(f"WARNING: Extrapolation at start (tC={tC}, tF={tF}, delta={(x_adv[0] - x[0]) / dx:.2f})")
            if x_adv[-1] < x[-1] - tol:
                print(f"WARNING: Extrapolation at end (tC={tC}, tF={tF}, delta={(x[-1] - x_adv[-1]) / dx:.2f})")

            # Interpolate the advected data onto the original grid
            q_WI_adv[tC, tF, :] = spline(x)

    return q_WI_adv

def is_increasing(arr):
    """Check if the array is strictly increasing."""
    return np.all(np.diff(arr) > 0)

def Gaussian_weight(t, t0,norm_flag=0):
    if norm_flag==1: # * Normalized
        return np.exp(-((t - 0) ** 2) / (2 * t0 ** 2)) / np.sqrt(2 * np.pi * t0 ** 2) 
    else: # * Not normalized
        return np.exp(-((t - 0) ** 2) / (2 * t0 ** 2))

def solve_advSoln_forced(q0,q_WI,q_inflow,dt_dns,frames,frames_obs,obsNt,xqI,c_adv,weights,sig0,tau,CN_flag=False):
    tC_max    = -1
    nOBS      = len(frames_obs)
    t_dns     = dt_dns*frames
    t_obs     = dt_dns*frames_obs
    dx        = xqI[1] - xqI[0]
    q         = np.zeros((len(frames),len(xqI)))
    q[0,:]    = q0 # * Initialize the solution at t=0
    for num, frame in enumerate(frames[:-1]):
        tC = int(num/obsNt) # * Denotes the closest observation frame less than the current DNS frame
        if (tC > tC_max):  
            obs_nbh, weights_nbh = get_obs_weights(tC,nOBS,weights)
                # * Get the moving target for each of the neighboring observations.
            q_WI_adv             = get_WI_advected(q_WI[obs_nbh,:],t_dns,t_obs,obsNt,obs_nbh,xqI,c_adv)
            tC_max               = tC
        frc     = np.zeros_like(xqI)
        num_rel = num - (obs_nbh[0])*obsNt # * Relative time index of moving target
        for n,nbh in enumerate(obs_nbh):
            discrep = q[num,:] - q_WI_adv[n,num_rel,:]
            frc     = frc - sig0 * Gaussian_weight(t_dns[num]-t_obs[nbh],tau) * discrep
        if CN_flag: 
            q_fac = np.zeros_like(xqI)
            sig   = 0.0
            for n,nbh in enumerate(obs_nbh):
                sig   += sig0 * Gaussian_weight(t_dns[num]-t_obs[nbh],tau) 
                q_fac += sig0 * Gaussian_weight(t_dns[num]-t_obs[nbh],tau) * q_WI_adv[n,num_rel+1,:] 
             # * Get the solution at next time step
            q[num+1,:] = march_advSoln_forced(q[num,:],q_inflow[num+1],dt_dns,dx,frc,c_adv,CN_flag,sig,q_fac)
        else:
            q[num+1,:] = march_advSoln_forced(q[num,:],q_inflow[num+1],dt_dns,dx,frc,c_adv,CN_flag)
    
    return q

def march_advSoln_forced(phi,phi_in,dt,dx,rhs,c_adv,CN_flag=False,sig=0.0,q_fac=0.0):
    # * rhs --> The forcing term, in combined form.
    
    # ~ Explicit upwind scheme for advection term
    """
    Solve the advection equation with a forcing term using the explicit upwind scheme.

    Parameters:
    phi (np.array): The solution at the current time step
    phi_in (float): The inflow boundary condition
    dt (float): Time step
    dx (float): Grid spacing
    rhs (np.array): The forcing term
    c_adv (float): Advection velocity

    Returns:
    np.array: The solution at the next time step
    """
    Nx      = len(phi)
    cfl     = c_adv * dt/dx 
    phi_new = np.copy(phi)
    tol     = 1e-16

    # for i in range(1, Nx - 1):  # ~ Upwind for positive ubar
    #     phi_new[i] = phi[i] - cfl * (phi[i] - phi[i - 1])
    # !!! WARNING: This approach produces numerical diffusion errors

    # ~ Explicit upwind scheme for convection term
        # & Second-order scheme for i=1 and i=Nx-2
    phi_new[1] = phi[1] - cfl * 0.5 * (phi[2] - phi[0]) + (cfl**2 / 2) * (phi[2] - 2 * phi[1] + phi[0])
    phi_new[Nx-2] = phi[Nx-2] - cfl * 0.5 * (phi[Nx-1] - phi[Nx-3]) + (cfl**2 / 2) * (phi[Nx-1] - 2 * phi[Nx-2] + phi[Nx-3])
    
        # & Third-order scheme for rest
    for i in range(2, Nx - 2):
        phi_new[i] = (phi[i] - (cfl / 6) * (2 * phi[i+1] + 3 * phi[i] - 6 * phi[i-1] + phi[i-2]) 
                      + (cfl**2 / 2) * (phi[i+1] - 2 * phi[i] + phi[i-1])
                      - (cfl**3 / 6) * (phi[i+2] - 3 * phi[i+1] + 3 * phi[i] - phi[i-1])) 
    
    if CN_flag:
        # ~ Crank-Nicolson scheme for forcing term 
        phi_new[1:] = (2*phi_new[1:] + dt * (rhs[1:] + q_fac[1:]) )/( 2.0 + dt * sig)
    else:
        # ~ Explicit update for the forcing term (Forward Euler)
        phi_new[1:] = phi_new[1:] + dt * rhs[1:]
    
    # ~ Apply boundary conditions
    phi_new[0]  = phi_in                             # & Inflow boundary condition
    # phi_new[-1] = phi[-1]- cfl * (phi[-1] - phi[-2]) # * First order
    # phi_new[-1] = phi[-1] - cfl * (1.5 * phi[-1] - 2.0 * phi[-2] + 0.5 * phi[-3])  # * Second-order
    phi_new[-1] = phi[-1] - cfl * (11/6 * phi[-1] - 3 * phi[-2] + 1.5 * phi[-3] - 1/3 * phi[-4])  # * Third-order

    return phi_new

def calculate_Re_tau_del(Re_del):
    """
    Calculate Re_tau based on the given relation between Reynolds numbers.

    Parameters:
    Re_del (float): Reynolds number

    Returns:
    float: Re_tau
    """
    k3=0.005;b1=0.9;b2=1.2;k4=np.power(k3,b1-0.5)
    return k4 * np.power(Re_del, b1) * np.power((1.0 + np.power(k3 * Re_del, -b2)), (b1 - 0.5) / b2)

def wall_model(y_plus, kappa=0.41):
    """
    Calculate the velocity u in the near-wall region of turbulent flow.
    
    Parameters:
    y_plus (float or np.array): Dimensionless wall distance (y+)
    kappa (float, optional): von Kármán constant, default is 0.41
    
    Returns:
    float or np.array: Velocity u at the given y+
    """
    term1 = (1 / kappa) * np.log(1 + kappa * y_plus)
    term2 = 7.8 * (1 - np.exp(-y_plus / 11) - ( y_plus / 11 ) * np.exp(-y_plus / 3))
    u = (term1 + term2)
    
    return u

def get_wall_modelled_profile(y_del, u_del, w_del, y_grid, delta):
    """
    Returns the wall-modelled profile for a given y-location, y-del, and velocity components, u-del and w-del.

    Parameters
    ----------
    y_del : float
        y-location of the wall
    u_del : float
        u-velocity component at y-del
    w_del : float
        w-velocity component at y-del
    y_grid : array_like
        array of y-locations to generate the wall-modelled profile
    delta : float
        viscous layer thickness scale

    Returns
    -------
    u : array_like
        u-velocity component of the wall-modelled profile
    w : array_like
        w-velocity component of the wall-modelled profile
    """
    y_del_plus = y_del/delta
    horVel     = np.sqrt(u_del**2 + w_del**2)
    sin        = w_del/horVel
    cos        = u_del/horVel
    Re_tau_del = calculate_Re_tau_del(horVel * y_del_plus)
    u_tau      = Re_tau_del/y_del_plus
    u, w       = np.zeros(len(y_grid)), np.zeros(len(y_grid))
    for j in range(len(y_grid)):
        y_plus   = y_grid[j] * Re_tau_del / y_del
        horVel_y = u_tau * wall_model(y_plus)
        u[j]     = horVel_y * cos
        w[j]     = horVel_y * sin
    
    return u, w

def get_wall_modelled_profile_with_vQuadratic(y_del, u_del, v_del, w_del, y_grid, delta):
    """
    Returns the wall-modelled profile for a given y-location, y-del, and velocity components, u-del and w-del.
    Returns a quadratic profile for the vertical velocity

    Parameters
    ----------
    y_del : float
        y-location of the wall
    u_del : float
        u-velocity component at y-del
    v_del : float
        v-velocity component at y-del
    w_del : float
        w-velocity component at y-del
    y_grid : array_like
        array of y-locations to generate the wall-modelled profile
    delta : float
        viscous layer thickness scale

    Returns
    -------
    u : array_like
        u-velocity component of the wall-modelled profile
    v : array_like
        v-velocity component of the wall-modelled profile
    w : array_like
        w-velocity component of the wall-modelled profile
    """
    y_del_plus = y_del/delta
    horVel     = np.sqrt(u_del**2 + w_del**2)
    sin        = w_del/horVel
    cos        = u_del/horVel
    Re_tau_del = calculate_Re_tau_del(horVel * y_del_plus)
    u_tau      = Re_tau_del/y_del_plus
    u, v, w    = np.zeros(len(y_grid)), np.zeros(len(y_grid)), np.zeros(len(y_grid))
    for j in range(len(y_grid)):
        y_plus   = y_grid[j] * Re_tau_del / y_del
        horVel_y = u_tau * wall_model(y_plus)
        u[j]     = horVel_y * cos
        w[j]     = horVel_y * sin
        v[j]     = (v_del / (y_del**2)) * (y_grid[j]**2) 
    
    return u, v, w