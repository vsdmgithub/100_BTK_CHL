# %%
import numpy as np
import os

def get_filename(number):
    """
    Generate a filename based on the given number.

    Parameters:
        number (int): The number to be used in the filename.

    Returns:
        str: The generated filename in the format "Flow_{number:05d}.dat", where {number:05d} is a zero-padded 5-digit integer representation of the input number.
    """
    # return f"Flow_{number:05d}.dat" #^ Yuhui's data (Re_tau ~ 3000)
    return f"Vel_{number:05d}.dat"    #^ Sukbong's data (Re_tau ~ 1500)

def get_measurement_location(file_addr, nx, ny, nz, debug=False):
    """
    Reads measurement locations from a specified file and returns the 1D arrays
    of the measurement grid points.

    Parameters:
    - file_addr (str): Path to the file containing measurement data.
    - nx (int): Number of grid points along the x-axis.
    - ny (int): Number of grid points along the y-axis.
    - nz (int): Number of grid points along the z-axis.
    - debug (bool): If True, print debug information. Default is False.

    Returns:
    - xx (ndarray): 1D array of x-coordinates of the measurement grid points.
    - yy (ndarray): 1D array of y-coordinates of the measurement grid points.
    - zz (ndarray): 1D array of z-coordinates of the measurement grid points.
    """
    # Check if the specified file exists
    if not os.path.isfile(file_addr):
        raise FileNotFoundError(f"The file at {file_addr} does not exist.")
    else:
        if debug:
            print(f"Reading measurement locations from {file_addr}")

    # Load the data from the file, skipping the first 3 rows
    data = np.loadtxt(file_addr, skiprows=3)

    # Extract each column into separate variables
    X_flat = data[:, 0]
    Y_flat = data[:, 1]
    Z_flat = data[:, 2]

    # Reshape the flat arrays back into their 3D grid representations
    X = X_flat.reshape((nz, ny, nx))
    Y = Y_flat.reshape((nz, ny, nx))
    Z = Z_flat.reshape((nz, ny, nx))

    # Extract the 1D arrays xx, yy, zz from the 3D grids
    xx = X[0, 0, :]
    yy = Y[0, :, 0]
    zz = Z[:, 0, 0]

    if debug:
        print("Done reading measurement locations")

    return xx, yy, zz

def get_measurement_velocities(file_addr, nx, ny, nz, debug=False, skiprows=3):
    # Check if the specified file exists
    """
    Reads velocity components from a specified file and returns the 3D arrays
    of the velocity components at the measurement grid points.

    Parameters:
    - file_addr (str): Path to the file containing measurement data.
    - nx (int): Number of grid points along the x-axis.
    - ny (int): Number of grid points along the y-axis.
    - nz (int): Number of grid points along the z-axis.
    - debug (bool): If True, print debug information. Default is False.

    Returns:
    - u (ndarray): 3D array of velocity components in the x-direction at the measurement grid points.
    - v (ndarray): 3D array of velocity components in the y-direction at the measurement grid points.
    - w (ndarray): 3D array of velocity components in the z-direction at the measurement grid points.
    """
    if not os.path.isfile(file_addr):
        raise FileNotFoundError(f"The file at {file_addr} does not exist.")
    else:
        if debug:
            print(f"Reading measurement locations from {file_addr}")

        # * Load the data from the file
    data = np.loadtxt(file_addr, skiprows=skiprows)

        # * Extract the velocity components from the reshuffled data array
    U = data[:, 3]  # Velocity in the x-direction
    V = data[:, 4]  # Velocity in the y-direction
    W = data[:, 5]  # Velocity in the z-direction

        # * Reshape the velocity arrays to the grid  
    u = U.reshape((nz, ny, nx))
    v = V.reshape((nz, ny, nx)) 
    w = W.reshape((nz, ny, nx)) 

    return u, v, w

def get_measurement_velocities_and_pressure(file_addr, nx, ny, nz, debug=False):
    if not os.path.isfile(file_addr):
        raise FileNotFoundError(f"The file at {file_addr} does not exist.")
    else:
        if debug:
            print(f"Reading measurement locations from {file_addr}")

        # * Load the data from the file
    data = np.loadtxt(file_addr, skiprows=3)

        # * Extract the velocity components from the reshuffled data array
    U = data[:, 3]  # Velocity in the x-direction
    V = data[:, 4]  # Velocity in the y-direction
    W = data[:, 5]  # Velocity in the z-direction
    P = data[:, 6]  # Pressure

        # * Reshape the velocity arrays to the grid  
    u = U.reshape((nz, ny, nx))
    v = V.reshape((nz, ny, nx)) 
    w = W.reshape((nz, ny, nx)) 
    P = P.reshape((nz, ny, nx))

    return u, v, w, P

def get_measurement_pressure(file_addr, nx, ny, nz, debug=False):
    """
    Reads pressure data from a specified file and returns the 3D array
    of pressure values at the measurement grid points.

    Parameters:
    - file_addr (str): Path to the file containing measurement data.
    - nx (int): Number of grid points along the x-axis.
    - ny (int): Number of grid points along the y-axis.
    - nz (int): Number of grid points along the z-axis.
    - debug (bool): If True, print debug information. Default is False.

    Returns:
    - P (ndarray): 3D array of pressure values at the measurement grid points.
    """
    if not os.path.isfile(file_addr):
        raise FileNotFoundError(f"The file at {file_addr} does not exist.")
    else:
        if debug:
            print(f"Reading pressure data from {file_addr}")

    # Load the data from the file
    data = np.loadtxt(file_addr, skiprows=3)

    # Extract the pressure data from the reshuffled data array
    P_flat = data[:, 6]  # Pressure

    # Reshape the pressure array to the grid
    P = P_flat.reshape((nz, ny, nx))

    return P

def process_expData(ind, frame, nex, ney, nez, root_data, skiprows=3):
    """
    Process the experimental data at a given frame.

    Parameters:
    - ind (int): Index of the frame.
    - frame (int): Frame number.
    - root_data (str): Root directory of the experimental data.

    Returns:
    - ind (int): Index of the frame.
    - u (ndarray): 3D array of velocity components in the x-direction.
    - v (ndarray): 3D array of velocity components in the y-direction.
    - w (ndarray): 3D array of velocity components in the z-direction.
    """
    file_name = get_filename(frame) 
    file_addr = root_data + file_name
    u,v,w     = get_measurement_velocities(file_addr, nex, ney, nez, skiprows=skiprows)

    return ind, u, v, w

def get_measurement_statistics(num,frame,nex,ney,nez,root_data):
    file_addr = root_data + get_filename(frame)
    u_exp, v_exp, w_exp, P_exp = get_measurement_velocities_and_pressure(file_addr,nex,ney,nez)
    P_exp      = P_exp - np.mean(P_exp)
    u_mean_exp = np.mean(u_exp,axis=(0,2))
    v_mean_exp = np.mean(v_exp,axis=(0,2))
    w_mean_exp = np.mean(w_exp,axis=(0,2))
    pp_exp     = np.mean(P_exp**2,axis=(0,2)) 
    uu_fluc    = np.mean(u_exp**2,axis=(0,2))    - u_mean_exp**2
    vv_fluc    = np.mean(v_exp**2,axis=(0,2))    - v_mean_exp**2
    ww_fluc    = np.mean(w_exp**2,axis=(0,2))    - w_mean_exp**2
    uv_fluc    = np.mean(u_exp*v_exp,axis=(0,2)) - u_mean_exp*v_mean_exp

    uu_rms = np.sqrt(uu_fluc)
    vv_rms = np.sqrt(vv_fluc)
    ww_rms = np.sqrt(ww_fluc)
    uv_rms = np.sqrt(-uv_fluc)

    return num, uu_rms, vv_rms, ww_rms, uv_rms

def get_measurement_statistics_onlyVel(num,frame,nex,ney,nez,root_data):
    file_addr = root_data + get_filename(frame)
    u_exp, v_exp, w_exp = get_measurement_velocities(file_addr,nex,ney,nez)
    u_mean_exp = np.mean(u_exp,axis=(0,2))
    v_mean_exp = np.mean(v_exp,axis=(0,2))
    w_mean_exp = np.mean(w_exp,axis=(0,2))
    uu_fluc    = np.mean(u_exp**2,axis=(0,2))    - u_mean_exp**2
    vv_fluc    = np.mean(v_exp**2,axis=(0,2))    - v_mean_exp**2
    ww_fluc    = np.mean(w_exp**2,axis=(0,2))    - w_mean_exp**2
    uv_fluc    = np.mean(u_exp*v_exp,axis=(0,2)) - u_mean_exp*v_mean_exp

    uu_rms = np.sqrt(uu_fluc)
    vv_rms = np.sqrt(vv_fluc)
    ww_rms = np.sqrt(ww_fluc)
    uv_rms = np.sqrt(np.abs(uv_fluc))

    return num, uu_rms, vv_rms, ww_rms, uv_rms