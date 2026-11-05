"""
Numba-accelerated IO operations for CFD data processing.

This module provides optimized functions for reading/writing binary Fortran files
used in subdomain CFD simulations. Numba JIT compilation is used to accelerate
array operations (transpose, reshape, slicing) while file I/O remains in Python.

Usage:
------
    from io_operations_Numba import (
        read_ucvcwc_file, write_ucvcwc_file,
        read_boundary_file, write_boundary_file,
        get_boundary_data_parallel,
        # ... other functions
    )

    # Read velocity fields
    num, U, V, W = read_ucvcwc_file(fileroot, ngx, ngy, ngz, timestep)

    # Write velocity fields
    write_ucvcwc_file(U, V, W, num=timestep, addr='./output/')

    # Read boundary data in parallel (fast for many timesteps)
    u, v, w = get_boundary_data_parallel(loc, ngx, ngy, ngz, frames, flag)

Notes:
------
- First call to JIT-compiled functions will be slower due to compilation
- Subsequent calls benefit from compiled code cache
- For parallel functions, ensure arrays are contiguous (use np.ascontiguousarray if needed)
"""

import numpy as np
from numba import njit, prange
import concurrent.futures
from tqdm import tqdm

# =============================================================================
# NUMBA-ACCELERATED ARRAY OPERATIONS
# =============================================================================

@njit(cache=True)
def _transpose_3d(arr):
    """
    Transpose a 3D array (swap axes 0 and 2).
    Equivalent to np.transpose(arr) for 3D arrays.
    """
    nz, nx, ny = arr.shape
    out = np.empty((ny, nx, nz), dtype=arr.dtype)
    for i in range(nz):
        for j in range(nx):
            for k in range(ny):
                out[k, j, i] = arr[i, j, k]
    return out

@njit(cache=True)
def _transpose_2d(arr):
    """
    Transpose a 2D array.
    """
    n0, n1 = arr.shape
    out = np.empty((n1, n0), dtype=arr.dtype)
    for i in range(n0):
        for j in range(n1):
            out[j, i] = arr[i, j]
    return out

@njit(cache=True)
def _reshape_and_transpose_fortran(data, shape):
    """
    Reshape flat array to 3D (Fortran order) and transpose.
    Combines reshape + transpose in one pass for efficiency.

    Parameters:
    -----------
    data : 1D array
        Flat input data
    shape : tuple (ny, nx, nz)
        Target shape before transpose (Fortran ordering from file)

    Returns:
    --------
    3D array with shape (nz, nx, ny)
    """
    ny, nx, nz = shape
    out = np.empty((nz, nx, ny), dtype=data.dtype)
    idx = 0
    for k in range(nz):
        for j in range(nx):
            for i in range(ny):
                out[k, j, i] = data[idx]
                idx += 1
    return out

@njit(cache=True, parallel=True)
def _reshape_and_transpose_fortran_parallel(data, shape):
    """
    Parallel version of reshape + transpose for large arrays.
    """
    ny, nx, nz = shape
    out = np.empty((nz, nx, ny), dtype=data.dtype)
    for k in prange(nz):
        for j in range(nx):
            for i in range(ny):
                idx = i + j * ny + k * ny * nx
                out[k, j, i] = data[idx]
    return out

@njit(cache=True)
def _prepare_for_write(arr):
    """
    Prepare 3D array for writing to Fortran binary file.
    Returns transposed array as contiguous bytes-compatible array.
    """
    nz, nx, ny = arr.shape
    out = np.empty((ny, nx, nz), dtype=arr.dtype)
    for i in range(nz):
        for j in range(nx):
            for k in range(ny):
                out[k, j, i] = arr[i, j, k]
    return out

@njit(cache=True)
def _prepare_2d_for_write(arr):
    """
    Prepare 2D array for writing to Fortran binary file.
    """
    n0, n1 = arr.shape
    out = np.empty((n1, n0), dtype=arr.dtype)
    for i in range(n0):
        for j in range(n1):
            out[j, i] = arr[i, j]
    return out

@njit(cache=True, parallel=True)
def _extract_boundary_bot_top(U, V, W, flag):
    """
    Extract bottom (flag=1) or top (flag=2) boundary from 3D fields.
    Returns 2D slices for U, V, W.
    """
    nz, nx, ny_u = U.shape
    _, _, ny_v = V.shape
    nz_w, _, _ = W.shape

    if flag == 1:  # bottom
        Ubc = np.empty((nz, nx), dtype=U.dtype)
        Vbc = np.empty((nz, nx), dtype=V.dtype)
        Wbc = np.empty((nz_w, nx), dtype=W.dtype)
        for i in prange(nz):
            for j in range(nx):
                Ubc[i, j] = U[i, j, 0]
        for i in prange(nz):
            for j in range(nx):
                Vbc[i, j] = V[i, j, 1]
        for i in prange(nz_w):
            for j in range(nx):
                Wbc[i, j] = W[i, j, 0]
    else:  # top (flag == 2)
        Ubc = np.empty((nz, nx), dtype=U.dtype)
        Vbc = np.empty((nz, nx), dtype=V.dtype)
        Wbc = np.empty((nz_w, nx), dtype=W.dtype)
        for i in prange(nz):
            for j in range(nx):
                Ubc[i, j] = U[i, j, ny_u - 1]
        for i in prange(nz):
            for j in range(nx):
                Vbc[i, j] = V[i, j, ny_v - 2]
        for i in prange(nz_w):
            for j in range(nx):
                Wbc[i, j] = W[i, j, ny_u - 1]

    return Ubc, Vbc, Wbc

@njit(cache=True, parallel=True)
def _extract_boundary_lft_rht(U, V, W, flag):
    """
    Extract left (flag=3) or right (flag=4) boundary from 3D fields.
    """
    nz, nx_u, ny = U.shape
    _, nx_v, ny_v = V.shape
    nz_w, _, _ = W.shape

    if flag == 3:  # left
        Ubc = np.empty((nz, ny), dtype=U.dtype)
        Vbc = np.empty((nz, ny_v), dtype=V.dtype)
        Wbc = np.empty((nz_w, ny), dtype=W.dtype)
        for i in prange(nz):
            for k in range(ny):
                Ubc[i, k] = U[i, 1, k]
        for i in prange(nz):
            for k in range(ny_v):
                Vbc[i, k] = V[i, 0, k]
        for i in prange(nz_w):
            for k in range(ny):
                Wbc[i, k] = W[i, 0, k]
    else:  # right (flag == 4)
        Ubc = np.empty((nz, ny), dtype=U.dtype)
        Vbc = np.empty((nz, ny_v), dtype=V.dtype)
        Wbc = np.empty((nz_w, ny), dtype=W.dtype)
        for i in prange(nz):
            for k in range(ny):
                Ubc[i, k] = U[i, nx_u - 2, k]
        for i in prange(nz):
            for k in range(ny_v):
                Vbc[i, k] = V[i, nx_v - 1, k]
        for i in prange(nz_w):
            for k in range(ny):
                Wbc[i, k] = W[i, nx_v - 1, k]

    return Ubc, Vbc, Wbc

@njit(cache=True, parallel=True)
def _extract_boundary_bck_frt(U, V, W, flag):
    """
    Extract back (flag=5) or front (flag=6) boundary from 3D fields.
    """
    nz, nx_u, ny = U.shape
    _, nx_v, ny_v = V.shape
    nz_w, _, _ = W.shape

    if flag == 5:  # back
        Ubc = np.empty((nx_u, ny), dtype=U.dtype)
        Vbc = np.empty((nx_v, ny_v), dtype=V.dtype)
        Wbc = np.empty((nx_v, ny), dtype=W.dtype)
        for j in prange(nx_u):
            for k in range(ny):
                Ubc[j, k] = U[0, j, k]
        for j in prange(nx_v):
            for k in range(ny_v):
                Vbc[j, k] = V[0, j, k]
        for j in prange(nx_v):
            for k in range(ny):
                Wbc[j, k] = W[1, j, k]
    else:  # front (flag == 6)
        Ubc = np.empty((nx_u, ny), dtype=U.dtype)
        Vbc = np.empty((nx_v, ny_v), dtype=V.dtype)
        Wbc = np.empty((nx_v, ny), dtype=W.dtype)
        for j in prange(nx_u):
            for k in range(ny):
                Ubc[j, k] = U[nz - 1, j, k]
        for j in prange(nx_v):
            for k in range(ny_v):
                Vbc[j, k] = V[nz - 1, j, k]
        for j in prange(nx_v):
            for k in range(ny):
                Wbc[j, k] = W[nz_w - 2, j, k]

    return Ubc, Vbc, Wbc

# =============================================================================
# WRITE FUNCTIONS
# =============================================================================

def write_ucvcwc_file(U, V, W, num=0, addr='./'):
    """
    Writes ucvcwc restart file at a given timestep.

    Parameters:
    -----------
    addr : str
        Name of the folder to write to
    U, V, W : ndarray
        Velocity components
    num : int
        Number of the file to write (default 0)
    """
    # Use Numba-accelerated transpose
    U_out = _prepare_for_write(np.ascontiguousarray(U))
    V_out = _prepare_for_write(np.ascontiguousarray(V))
    W_out = _prepare_for_write(np.ascontiguousarray(W))

    fname = addr + '01uc00.dble.' + "{:07d}".format(num)
    with open(fname, 'wb') as fout:
        fout.write(U_out.tobytes())

    fname = addr + '02vc00.dble.' + "{:07d}".format(num)
    with open(fname, 'wb') as fout:
        fout.write(V_out.tobytes())

    fname = addr + '03wc00.dble.' + "{:07d}".format(num)
    with open(fname, 'wb') as fout:
        fout.write(W_out.tobytes())

def write_boundary_file(U, V, W, num=0, addr='./', flag=-1):
    """
    Writes boundary files for subdomain calculation.

    Parameters:
    -----------
    U, V, W : ndarray
        Velocity components
    num : int
        Timestep
    addr : str
        Name of the folder to write to
    flag : int
        Flag indicating which boundary to write:
        -1: write all boundaries
        1: bot boundary
        2: top boundary
        3: left boundary
        4: right boundary
        5: back boundary
        6: front boundary
    """
    U = np.ascontiguousarray(U)
    V = np.ascontiguousarray(V)
    W = np.ascontiguousarray(W)

    if flag == -1 or flag == 1:
        Ubc, Vbc, Wbc = _extract_boundary_bot_top(U, V, W, 1)
        fname = addr + 'bc_bot.dble.' + "{:07d}".format(num)
        with open(fname, 'wb') as fout:
            fout.write(_prepare_2d_for_write(Ubc).tobytes())
            fout.write(_prepare_2d_for_write(Vbc).tobytes())
            fout.write(_prepare_2d_for_write(Wbc).tobytes())

    if flag == -1 or flag == 2:
        Ubc, Vbc, Wbc = _extract_boundary_bot_top(U, V, W, 2)
        fname = addr + 'bc_top.dble.' + "{:07d}".format(num)
        with open(fname, 'wb') as fout:
            fout.write(_prepare_2d_for_write(Ubc).tobytes())
            fout.write(_prepare_2d_for_write(Vbc).tobytes())
            fout.write(_prepare_2d_for_write(Wbc).tobytes())

    if flag == -1 or flag == 3:
        Ubc, Vbc, Wbc = _extract_boundary_lft_rht(U, V, W, 3)
        fname = addr + 'bc_lft.dble.' + "{:07d}".format(num)
        with open(fname, 'wb') as fout:
            fout.write(_prepare_2d_for_write(Ubc).tobytes())
            fout.write(_prepare_2d_for_write(Vbc).tobytes())
            fout.write(_prepare_2d_for_write(Wbc).tobytes())

    if flag == -1 or flag == 4:
        Ubc, Vbc, Wbc = _extract_boundary_lft_rht(U, V, W, 4)
        fname = addr + 'bc_rht.dble.' + "{:07d}".format(num)
        with open(fname, 'wb') as fout:
            fout.write(_prepare_2d_for_write(Ubc).tobytes())
            fout.write(_prepare_2d_for_write(Vbc).tobytes())
            fout.write(_prepare_2d_for_write(Wbc).tobytes())

    if flag == -1 or flag == 5:
        Ubc, Vbc, Wbc = _extract_boundary_bck_frt(U, V, W, 5)
        fname = addr + 'bc_bck.dble.' + "{:07d}".format(num)
        with open(fname, 'wb') as fout:
            fout.write(_prepare_2d_for_write(Ubc).tobytes())
            fout.write(_prepare_2d_for_write(Vbc).tobytes())
            fout.write(_prepare_2d_for_write(Wbc).tobytes())

    if flag == -1 or flag == 6:
        Ubc, Vbc, Wbc = _extract_boundary_bck_frt(U, V, W, 6)
        fname = addr + 'bc_frt.dble.' + "{:07d}".format(num)
        with open(fname, 'wb') as fout:
            fout.write(_prepare_2d_for_write(Ubc).tobytes())
            fout.write(_prepare_2d_for_write(Vbc).tobytes())
            fout.write(_prepare_2d_for_write(Wbc).tobytes())

def write_boundary_file_from_six_boundaries(Ulft, Vlft, Wlft,
                                            Urht, Vrht, Wrht,
                                            Utop, Vtop, Wtop,
                                            Ubot, Vbot, Wbot,
                                            Ubck, Vbck, Wbck,
                                            Ufrt, Vfrt, Wfrt,
                                            num=0, bc_addr='./'):
    """
    Writes boundary files for the given boundary conditions.

    Parameters:
    -----------
    Ulft, Vlft, Wlft : ndarray
        Left boundary velocity components
    Urht, Vrht, Wrht : ndarray
        Right boundary velocity components
    Ubot, Vbot, Wbot : ndarray
        Bottom boundary velocity components
    Utop, Vtop, Wtop : ndarray
        Top boundary velocity components
    Ubck, Vbck, Wbck : ndarray
        Back boundary velocity components
    Ufrt, Vfrt, Wfrt : ndarray
        Front boundary velocity components
    num : int
        The timestep number (default is 0)
    bc_addr : str
        The address of the boundary file directory (default is './')
    """
    # Left boundary
    fname = bc_addr + 'bc_lft.dble.' + "{:07d}".format(num)
    with open(fname, 'wb') as fout:
        fout.write(_prepare_2d_for_write(np.ascontiguousarray(Ulft[:, 0, :])).tobytes())
        fout.write(_prepare_2d_for_write(np.ascontiguousarray(Vlft[:, 0, :])).tobytes())
        fout.write(_prepare_2d_for_write(np.ascontiguousarray(Wlft[:, 0, :])).tobytes())

    # Right boundary
    fname = bc_addr + 'bc_rht.dble.' + "{:07d}".format(num)
    with open(fname, 'wb') as fout:
        fout.write(_prepare_2d_for_write(np.ascontiguousarray(Urht[:, 0, :])).tobytes())
        fout.write(_prepare_2d_for_write(np.ascontiguousarray(Vrht[:, 0, :])).tobytes())
        fout.write(_prepare_2d_for_write(np.ascontiguousarray(Wrht[:, 0, :])).tobytes())

    # Bottom boundary
    fname = bc_addr + 'bc_bot.dble.' + "{:07d}".format(num)
    with open(fname, 'wb') as fout:
        fout.write(_prepare_2d_for_write(np.ascontiguousarray(Ubot[:, :, 0])).tobytes())
        fout.write(_prepare_2d_for_write(np.ascontiguousarray(Vbot[:, :, 0])).tobytes())
        fout.write(_prepare_2d_for_write(np.ascontiguousarray(Wbot[:, :, 0])).tobytes())

    # Top boundary
    fname = bc_addr + 'bc_top.dble.' + "{:07d}".format(num)
    with open(fname, 'wb') as fout:
        fout.write(_prepare_2d_for_write(np.ascontiguousarray(Utop[:, :, 0])).tobytes())
        fout.write(_prepare_2d_for_write(np.ascontiguousarray(Vtop[:, :, 0])).tobytes())
        fout.write(_prepare_2d_for_write(np.ascontiguousarray(Wtop[:, :, 0])).tobytes())

    # Back boundary
    fname = bc_addr + 'bc_bck.dble.' + "{:07d}".format(num)
    with open(fname, 'wb') as fout:
        fout.write(_prepare_2d_for_write(np.ascontiguousarray(Ubck[0, :, :])).tobytes())
        fout.write(_prepare_2d_for_write(np.ascontiguousarray(Vbck[0, :, :])).tobytes())
        fout.write(_prepare_2d_for_write(np.ascontiguousarray(Wbck[0, :, :])).tobytes())

    # Front boundary
    fname = bc_addr + 'bc_frt.dble.' + "{:07d}".format(num)
    with open(fname, 'wb') as fout:
        fout.write(_prepare_2d_for_write(np.ascontiguousarray(Ufrt[0, :, :])).tobytes())
        fout.write(_prepare_2d_for_write(np.ascontiguousarray(Vfrt[0, :, :])).tobytes())
        fout.write(_prepare_2d_for_write(np.ascontiguousarray(Wfrt[0, :, :])).tobytes())

def write_boundary_file_from_boundary(Ubc, Vbc, Wbc, num=0, bc_addr='./', flag=-1, adjFlag=False):
    """
    Writes boundary files for the given boundary conditions.

    Parameters:
    -----------
    Ubc, Vbc, Wbc : ndarray
        Velocity components of the boundary
    num : int
        The timestep number (default is 0)
    bc_addr : str
        The address of the boundary file directory (default is './')
    flag : int
        The flag indicating which boundary to write:
        1: bot boundary
        2: top boundary
        3: left boundary
        4: right boundary
        5: back boundary
        6: front boundary
    adjFlag : bool
        If True, write adjoint boundary files
    """
    boundary_names = {
        1: ('bc_bot', 'adjbot'),
        2: ('bc_top', 'adjtop'),
        3: ('bc_lft', 'adjlft'),
        4: ('bc_rht', 'adjrht'),
        5: ('bc_bck', 'adjbck'),
        6: ('bc_frt', 'adjfrt')
    }

    if flag not in boundary_names:
        print("Invalid flag value", flush=True)
        return

    prefix = boundary_names[flag][1] if adjFlag else boundary_names[flag][0]
    fname = bc_addr + prefix + '.dble.' + "{:07d}".format(num)

    Ubc_out = _prepare_2d_for_write(np.ascontiguousarray(Ubc[:, :]))
    Vbc_out = _prepare_2d_for_write(np.ascontiguousarray(Vbc[:, :]))
    Wbc_out = _prepare_2d_for_write(np.ascontiguousarray(Wbc[:, :]))

    with open(fname, 'wb') as fout:
        fout.write(Ubc_out.tobytes())
        fout.write(Vbc_out.tobytes())
        fout.write(Wbc_out.tobytes())

# =============================================================================
# READ FUNCTIONS
# =============================================================================

def read_ucvcwc_file(fileroot, ngx, ngy, ngz, num, adjointFlag=0, component=-1):
    """
    Reads U,V,W velocity data from binary Fortran files.

    Parameters:
    -----------
    fileroot : str
        Location of data files
    ngx, ngy, ngz : int
        Number of grid points in each dimension
    num : int
        Time step to extract data
    adjointFlag : int
        [1] read adjoint field [else] read forward field
    component : int
        [-1] all, [0] U, [1] V, [2] W

    Returns:
    --------
    tuple
        (num, U, V, W) or (num, single_component) depending on component flag
    """
    if component != 1 and component != 2:
        fname = '01ucad.dble.' if adjointFlag == 1 else '01uc00.dble.'
        filename = fileroot + fname + "{:07d}".format(num)
        Udata = np.fromfile(filename, dtype=np.float64)
        Ufield = _reshape_and_transpose_fortran_parallel(
            Udata, (ngy + 1, ngx + 2, ngz + 1)
        )
        if component == 0:
            return num, Ufield

    if component != 0 and component != 2:
        fname = '02vcad.dble.' if adjointFlag == 1 else '02vc00.dble.'
        filename = fileroot + fname + "{:07d}".format(num)
        Vdata = np.fromfile(filename, dtype=np.float64)
        Vfield = _reshape_and_transpose_fortran_parallel(
            Vdata, (ngy + 2, ngx + 1, ngz + 1)
        )
        if component == 1:
            return num, Vfield

    if component != 0 and component != 1:
        fname = '03wcad.dble.' if adjointFlag == 1 else '03wc00.dble.'
        filename = fileroot + fname + "{:07d}".format(num)
        Wdata = np.fromfile(filename, dtype=np.float64)
        Wfield = _reshape_and_transpose_fortran_parallel(
            Wdata, (ngy + 1, ngx + 1, ngz + 2)
        )
        if component == 2:
            return num, Wfield

    return (num, Ufield, Vfield, Wfield)

def read_icvec_file(fileroot, ngx, ngy, ngz, num):
    """
    Reads U,V,W velocity data from new Fortran ICvec binary file.

    Parameters:
    -----------
    fileroot : str
        Path to directory containing ICvec files
    ngx, ngy, ngz : int
        Number of grid points in each direction
    num : int
        Time step / iteration number

    Returns:
    --------
    tuple
        (num, U, V, W) numpy arrays
    """
    filename = f"{fileroot}ICvec.{num:03d}"

    size_u = (ngz - 1) * (ngx - 2) * (ngy - 1)
    size_v = (ngz - 1) * (ngx - 1) * (ngy - 2)
    size_w = (ngz - 2) * (ngx - 1) * (ngy - 1)

    data = np.fromfile(filename, dtype=np.float64)

    if data.size != size_u + size_v + size_w:
        raise ValueError(
            f"File {filename} has {data.size} entries, "
            f"expected {size_u + size_v + size_w}"
        )

    u_flat = data[0:size_u]
    v_flat = data[size_u:size_u + size_v]
    w_flat = data[size_u + size_v:]

    U = np.reshape(u_flat, (ngz - 1, ngx - 2, ngy - 1), order='F')
    V = np.reshape(v_flat, (ngz - 1, ngx - 1, ngy - 2), order='F')
    W = np.reshape(w_flat, (ngz - 2, ngx - 1, ngy - 1), order='F')

    return num, U, V, W

def read_ufvfwf_file(fileroot, ngx, ngy, ngz, num):
    """
    Reads U,V,W velocity flux data from binary Fortran files.

    Parameters:
    -----------
    fileroot : str
        Location of data files
    ngx, ngy, ngz : int
        Number of grid points in each dimension
    num : int
        Time step to extract data

    Returns:
    --------
    tuple
        (num, Ufield, Vfield, Wfield)
    """
    fname = '04uu00.dble.'
    filename = fileroot + fname + "{:07d}".format(num)
    Udata = np.fromfile(filename, dtype=np.float64)
    Ufield = _reshape_and_transpose_fortran_parallel(
        Udata, (ngy + 1, ngx + 2, ngz + 1)
    )

    fname = '05vv00.dble.'
    filename = fileroot + fname + "{:07d}".format(num)
    Vdata = np.fromfile(filename, dtype=np.float64)
    Vfield = _reshape_and_transpose_fortran_parallel(
        Vdata, (ngy + 2, ngx + 1, ngz + 1)
    )

    fname = '06ww00.dble.'
    filename = fileroot + fname + "{:07d}".format(num)
    Wdata = np.fromfile(filename, dtype=np.float64)
    Wfield = _reshape_and_transpose_fortran_parallel(
        Wdata, (ngy + 1, ngx + 1, ngz + 2)
    )

    return (num, Ufield, Vfield, Wfield)

def read_ucvcwcPres_file(fileroot, ngx, ngy, ngz, num, adjointFlag=0, component=-1):
    """
    Reads U,V,W velocity and pressure data from binary Fortran files.

    Parameters:
    -----------
    fileroot : str
        Location of data files
    ngx, ngy, ngz : int
        Number of grid points in each dimension
    num : int
        Time step to extract data
    adjointFlag : int
        [1] read adjoint field [else] read forward field
    component : int
        [-1] all, [0] U, [1] V, [2] W, [3] Pressure

    Returns:
    --------
    tuple
        Velocity and pressure field components
    """
    fname = 'phiadj.dble.' if adjointFlag == 1 else 'pres_p.dble.'
    filename = fileroot + fname + "{:07d}".format(num)
    Pdata = np.fromfile(filename, dtype=np.float64)
    Pfield = _reshape_and_transpose_fortran_parallel(
        Pdata, (ngy + 1, ngx + 1, ngz + 1)
    )
    if component == 3:
        return num, Pfield

    fname = '01ucad.dble.' if adjointFlag == 1 else '01uc00.dble.'
    filename = fileroot + fname + "{:07d}".format(num)
    Udata = np.fromfile(filename, dtype=np.float64)
    Ufield = _reshape_and_transpose_fortran_parallel(
        Udata, (ngy + 1, ngx + 2, ngz + 1)
    )
    if component == 0:
        return num, Ufield

    fname = '02vcad.dble.' if adjointFlag == 1 else '02vc00.dble.'
    filename = fileroot + fname + "{:07d}".format(num)
    Vdata = np.fromfile(filename, dtype=np.float64)
    Vfield = _reshape_and_transpose_fortran_parallel(
        Vdata, (ngy + 2, ngx + 1, ngz + 1)
    )
    if component == 1:
        return num, Vfield

    fname = '03wcad.dble.' if adjointFlag == 1 else '03wc00.dble.'
    filename = fileroot + fname + "{:07d}".format(num)
    Wdata = np.fromfile(filename, dtype=np.float64)
    Wfield = _reshape_and_transpose_fortran_parallel(
        Wdata, (ngy + 1, ngx + 1, ngz + 2)
    )
    if component == 2:
        return num, Wfield

    return (num, Ufield, Vfield, Wfield, Pfield)

def read_obs_file(fileroot, ngx, ngy, ngz, num):
    """
    Read the observed velocity fields from binary files for a given time step.

    Parameters:
    -----------
    fileroot : str
        The root location of the data files
    ngx, ngy, ngz : int
        The number of grid points in each direction
    num : int
        The time step to extract data for

    Returns:
    --------
    tuple
        (Ufield, Vfield, Wfield) observed velocity fields
    """
    fname = '04uobs.dble.'
    filename = fileroot + fname + "{:07d}".format(num)
    Udata = np.fromfile(filename, dtype=np.float64, count=ngx * ngy * ngz)
    Ufield = _reshape_and_transpose_fortran_parallel(Udata, (ngy, ngx, ngz))

    fname = '05vobs.dble.'
    filename = fileroot + fname + "{:07d}".format(num)
    Vdata = np.fromfile(filename, dtype=np.float64, count=ngx * ngy * ngz)
    Vfield = _reshape_and_transpose_fortran_parallel(Vdata, (ngy, ngx, ngz))

    fname = '06wobs.dble.'
    filename = fileroot + fname + "{:07d}".format(num)
    Wdata = np.fromfile(filename, dtype=np.float64, count=ngx * ngy * ngz)
    Wfield = _reshape_and_transpose_fortran_parallel(Wdata, (ngy, ngx, ngz))

    return (Ufield, Vfield, Wfield)

def return_gradientBC_prefix(case):
    """Return gradient BC prefix based on case number."""
    prefixes = {0: 'adj', 1: 'ap1', 2: 'ap2', 3: 'ap3', 4: 'ap4', 5: 'ap5'}
    return prefixes.get(case, 'adj')

def read_boundary_file(fileroot, ngx, ngy, ngz, num, flag=-1, adjFlag=False, adjOpt=0):
    """
    Read boundary data files for/from subdomain Fortran code.

    Parameters:
    -----------
    fileroot : str
        Location of data files
    ngx, ngy, ngz : int
        Number of grid points in each dimension
    num : int
        Time step to extract data
    flag : int
        Boundary flag (1-6 or -1 for all)
    adjFlag : bool
        If True, read adjoint boundary files
    adjOpt : int
        Adjoint option for gradient BC prefix

    Returns:
    --------
    tuple
        Boundary velocity data
    """
    results = {}

    if flag == -1 or flag == 1:
        fname = return_gradientBC_prefix(adjOpt) + 'bot.dble.' if adjFlag else 'bc_bot.dble.'
        filename = fileroot + fname + "{:07d}".format(num)

        offset = 0
        Ubot = np.fromfile(filename, dtype=np.float64, count=(ngx + 2) * (ngz + 1), offset=offset)
        Ubot = _reshape_and_transpose_fortran(Ubot, (1, ngx + 2, ngz + 1))

        offset = (ngz + 1) * (ngx + 2) * 8
        Vbot = np.fromfile(filename, dtype=np.float64, count=(ngx + 1) * (ngz + 1), offset=offset)
        Vbot = _reshape_and_transpose_fortran(Vbot, (1, ngx + 1, ngz + 1))

        offset = (ngz + 1) * (ngx + 2) * 8 + (ngz + 1) * (ngx + 1) * 8
        Wbot = np.fromfile(filename, dtype=np.float64, count=(ngx + 1) * (ngz + 2), offset=offset)
        Wbot = _reshape_and_transpose_fortran(Wbot, (1, ngx + 1, ngz + 2))

        if flag == 1:
            return (num, Ubot, Vbot, Wbot)
        results['bot'] = (Ubot, Vbot, Wbot)

    if flag == -1 or flag == 2:
        fname = return_gradientBC_prefix(adjOpt) + 'top.dble.' if adjFlag else 'bc_top.dble.'
        filename = fileroot + fname + "{:07d}".format(num)

        offset = 0
        Utop = np.fromfile(filename, dtype=np.float64, count=(ngx + 2) * (ngz + 1), offset=offset)
        Utop = _reshape_and_transpose_fortran(Utop, (1, ngx + 2, ngz + 1))

        offset = (ngz + 1) * (ngx + 2) * 8
        Vtop = np.fromfile(filename, dtype=np.float64, count=(ngx + 1) * (ngz + 1), offset=offset)
        Vtop = _reshape_and_transpose_fortran(Vtop, (1, ngx + 1, ngz + 1))

        offset = (ngz + 1) * (ngx + 2) * 8 + (ngz + 1) * (ngx + 1) * 8
        Wtop = np.fromfile(filename, dtype=np.float64, count=(ngx + 1) * (ngz + 2), offset=offset)
        Wtop = _reshape_and_transpose_fortran(Wtop, (1, ngx + 1, ngz + 2))

        if flag == 2:
            return (num, Utop, Vtop, Wtop)
        results['top'] = (Utop, Vtop, Wtop)

    if flag == -1 or flag == 3:
        fname = return_gradientBC_prefix(adjOpt) + 'lft.dble.' if adjFlag else 'bc_lft.dble.'
        filename = fileroot + fname + "{:07d}".format(num)

        offset = 0
        Ulft = np.fromfile(filename, dtype=np.float64, count=(ngz + 1) * (ngy + 1), offset=offset)
        Ulft = _reshape_and_transpose_fortran(Ulft, (ngy + 1, 1, ngz + 1))

        offset = (ngz + 1) * (ngy + 1) * 8
        Vlft = np.fromfile(filename, dtype=np.float64, count=(ngz + 1) * (ngy + 2), offset=offset)
        Vlft = _reshape_and_transpose_fortran(Vlft, (ngy + 2, 1, ngz + 1))

        offset = ((ngz + 1) * (ngy + 1) + (ngz + 1) * (ngy + 2)) * 8
        Wlft = np.fromfile(filename, dtype=np.float64, count=(ngz + 2) * (ngy + 1), offset=offset)
        Wlft = _reshape_and_transpose_fortran(Wlft, (ngy + 1, 1, ngz + 2))

        if flag == 3:
            return (num, Ulft, Vlft, Wlft)
        results['lft'] = (Ulft, Vlft, Wlft)

    if flag == -1 or flag == 4:
        fname = return_gradientBC_prefix(adjOpt) + 'rht.dble.' if adjFlag else 'bc_rht.dble.'
        filename = fileroot + fname + "{:07d}".format(num)

        offset = 0
        Urht = np.fromfile(filename, dtype=np.float64, count=(ngz + 1) * (ngy + 1), offset=offset)
        Urht = _reshape_and_transpose_fortran(Urht, (ngy + 1, 1, ngz + 1))

        offset = (ngz + 1) * (ngy + 1) * 8
        Vrht = np.fromfile(filename, dtype=np.float64, count=(ngz + 1) * (ngy + 2), offset=offset)
        Vrht = _reshape_and_transpose_fortran(Vrht, (ngy + 2, 1, ngz + 1))

        offset = (ngz + 1) * (ngy + 1) * 8 + (ngz + 1) * (ngy + 2) * 8
        Wrht = np.fromfile(filename, dtype=np.float64, count=(ngz + 2) * (ngy + 1), offset=offset)
        Wrht = _reshape_and_transpose_fortran(Wrht, (ngy + 1, 1, ngz + 2))

        if flag == 4:
            return (num, Urht, Vrht, Wrht)
        results['rht'] = (Urht, Vrht, Wrht)

    if flag == -1 or flag == 5:
        fname = return_gradientBC_prefix(adjOpt) + 'bck.dble.' if adjFlag else 'bc_bck.dble.'
        filename = fileroot + fname + "{:07d}".format(num)

        offset = 0
        Ubck = np.fromfile(filename, dtype=np.float64, count=(ngx + 2) * (ngy + 1), offset=offset)
        Ubck = _reshape_and_transpose_fortran(Ubck, (ngy + 1, ngx + 2, 1))

        offset = (ngy + 1) * (ngx + 2) * 8
        Vbck = np.fromfile(filename, dtype=np.float64, count=(ngx + 1) * (ngy + 2), offset=offset)
        Vbck = _reshape_and_transpose_fortran(Vbck, (ngy + 2, ngx + 1, 1))

        offset = (ngy + 1) * (ngx + 2) * 8 + (ngx + 1) * (ngy + 2) * 8
        Wbck = np.fromfile(filename, dtype=np.float64, count=(ngx + 1) * (ngy + 1), offset=offset)
        Wbck = _reshape_and_transpose_fortran(Wbck, (ngy + 1, ngx + 1, 1))

        if flag == 5:
            return (num, Ubck, Vbck, Wbck)
        results['bck'] = (Ubck, Vbck, Wbck)

    if flag == -1 or flag == 6:
        fname = return_gradientBC_prefix(adjOpt) + 'frt.dble.' if adjFlag else 'bc_frt.dble.'
        filename = fileroot + fname + "{:07d}".format(num)

        offset = 0
        Ufrt = np.fromfile(filename, dtype=np.float64, count=(ngx + 2) * (ngy + 1), offset=offset)
        Ufrt = _reshape_and_transpose_fortran(Ufrt, (ngy + 1, ngx + 2, 1))

        offset = (ngy + 1) * (ngx + 2) * 8
        Vfrt = np.fromfile(filename, dtype=np.float64, count=(ngx + 1) * (ngy + 2), offset=offset)
        Vfrt = _reshape_and_transpose_fortran(Vfrt, (ngy + 2, ngx + 1, 1))

        offset = (ngy + 1) * (ngx + 2) * 8 + (ngx + 1) * (ngy + 2) * 8
        Wfrt = np.fromfile(filename, dtype=np.float64, count=(ngx + 1) * (ngy + 1), offset=offset)
        Wfrt = _reshape_and_transpose_fortran(Wfrt, (ngy + 1, ngx + 1, 1))

        if flag == 6:
            return (num, Ufrt, Vfrt, Wfrt)
        results['frt'] = (Ufrt, Vfrt, Wfrt)

    # Return all boundaries
    Ulft, Vlft, Wlft = results['lft']
    Urht, Vrht, Wrht = results['rht']
    Utop, Vtop, Wtop = results['top']
    Ubot, Vbot, Wbot = results['bot']
    Ubck, Vbck, Wbck = results['bck']
    Ufrt, Vfrt, Wfrt = results['frt']

    return (num, Ulft, Vlft, Wlft, Urht, Vrht, Wrht, Utop, Vtop, Wtop,
            Ubot, Vbot, Wbot, Ubck, Vbck, Wbck, Ufrt, Vfrt, Wfrt)

# =============================================================================
# OBSERVATION FUNCTIONS
# =============================================================================

def get_observation_at_boundary(loc, nox, noy, noz, frames_obs, sf_u, sf_v, sf_w, flag=-1):
    """
    Reads the observation data for a given boundary and frame number.

    Parameters:
    -----------
    loc : str
        The root location of the observation data files
    nox, noy, noz : int
        The number of grid points in the x, y, and z directions
    frames_obs : list
        The list of time frames for which to read the observations
    sf_u, sf_v, sf_w : ndarray
        The flux-to-cartesian conversion factors for u, v, and w
    flag : int
        The flag specifying which boundary to read from (1-6)

    Returns:
    --------
    tuple
        (u_obs, v_obs, w_obs) observed velocity fields
    """
    nframes = len(frames_obs)

    if flag == 1 or flag == 2:
        u_obs = np.zeros((nframes, noz, nox))
        v_obs = np.zeros_like(u_obs)
        w_obs = np.zeros_like(u_obs)

        for num, frame in enumerate(frames_obs):
            u_obs_f, v_obs_f, w_obs_f = read_obs_file(loc, nox, noy, noz, frame)
            if flag == 1:
                u_obs[num, :, :] = u_obs_f[:, :, 0] / sf_u[:, :, 0]
                v_obs[num, :, :] = v_obs_f[:, :, 0] / sf_v[:, :, 0]
                w_obs[num, :, :] = w_obs_f[:, :, 0] / sf_w[:, :, 0]
            else:
                u_obs[num, :, :] = u_obs_f[:, :, -1] / sf_u[:, :, -1]
                v_obs[num, :, :] = v_obs_f[:, :, -1] / sf_v[:, :, -1]
                w_obs[num, :, :] = w_obs_f[:, :, -1] / sf_w[:, :, -1]

    elif flag == 3 or flag == 4:
        u_obs = np.zeros((nframes, noz, noy))
        v_obs = np.zeros_like(u_obs)
        w_obs = np.zeros_like(u_obs)

        for num, frame in enumerate(frames_obs):
            u_obs_f, v_obs_f, w_obs_f = read_obs_file(loc, nox, noy, noz, frame)
            if flag == 3:
                u_obs[num, :, :] = u_obs_f[:, 0, :] / sf_u[:, 0, :]
                v_obs[num, :, :] = v_obs_f[:, 0, :] / sf_v[:, 0, :]
                w_obs[num, :, :] = w_obs_f[:, 0, :] / sf_w[:, 0, :]
            else:
                u_obs[num, :, :] = u_obs_f[:, -1, :] / sf_u[:, -1, :]
                v_obs[num, :, :] = v_obs_f[:, -1, :] / sf_v[:, -1, :]
                w_obs[num, :, :] = w_obs_f[:, -1, :] / sf_w[:, -1, :]

    elif flag == 5 or flag == 6:
        u_obs = np.zeros((nframes, nox, noy))
        v_obs = np.zeros_like(u_obs)
        w_obs = np.zeros_like(u_obs)

        for num, frame in enumerate(frames_obs):
            u_obs_f, v_obs_f, w_obs_f = read_obs_file(loc, nox, noy, noz, frame)
            if flag == 5:
                u_obs[num, :, :] = u_obs_f[0, :, :] / sf_u[0, :, :]
                v_obs[num, :, :] = v_obs_f[0, :, :] / sf_v[0, :, :]
                w_obs[num, :, :] = w_obs_f[0, :, :] / sf_w[0, :, :]
            else:
                u_obs[num, :, :] = u_obs_f[-1, :, :] / sf_u[-1, :, :]
                v_obs[num, :, :] = v_obs_f[-1, :, :] / sf_v[-1, :, :]
                w_obs[num, :, :] = w_obs_f[-1, :, :] / sf_w[-1, :, :]
    else:
        raise Exception("Invalid observation boundary flag value")

    return u_obs, v_obs, w_obs

def get_boundary_data(loc, ngx, ngy, ngz, frames, flag, adjFlag=False, adjOpt=0):
    """
    Read boundary data for multiple frames using parallel processing.

    Parameters:
    -----------
    loc : str
        Location of data files
    ngx, ngy, ngz : int
        Grid dimensions
    frames : list
        List of frame numbers to read
    flag : int
        Boundary flag (1-6)
    adjFlag : bool
        If True, read adjoint files
    adjOpt : int
        Adjoint option

    Returns:
    --------
    tuple
        (u, v, w) boundary velocity arrays with shape (nframes, ...)
    """
    total_frames = len(frames)
    skip = frames[1] - frames[0] if len(frames) > 1 else 1

    if flag == 1 or flag == 2:
        u = np.zeros((total_frames, ngz + 1, ngx + 2))
        v = np.zeros((total_frames, ngz + 1, ngx + 1))
        w = np.zeros((total_frames, ngz + 2, ngx + 1))

        with concurrent.futures.ProcessPoolExecutor() as executor:
            future_to_frame = {
                executor.submit(read_boundary_file, loc, ngx, ngy, ngz, frame, flag, adjFlag, adjOpt): frame
                for frame in frames
            }

            for future in tqdm(concurrent.futures.as_completed(future_to_frame),
                               total=total_frames,
                               desc="Reading boundary data (flag 1/2)"):
                num, U, V, W = future.result()
                idx = (num - frames[0]) // skip
                u[idx, :, :] = U[:, :, 0]
                v[idx, :, :] = V[:, :, 0]
                w[idx, :, :] = W[:, :, 0]

        return u, v, w

    if flag == 3 or flag == 4:
        u = np.zeros((total_frames, ngz + 1, ngy + 1))
        v = np.zeros((total_frames, ngz + 1, ngy + 2))
        w = np.zeros((total_frames, ngz + 2, ngy + 1))

        with concurrent.futures.ProcessPoolExecutor() as executor:
            future_to_frame = {
                executor.submit(read_boundary_file, loc, ngx, ngy, ngz, frame, flag, adjFlag, adjOpt): frame
                for frame in frames
            }

            for future in tqdm(concurrent.futures.as_completed(future_to_frame),
                               total=total_frames,
                               desc="Reading boundary data (flag 3/4)"):
                num, U, V, W = future.result()
                idx = (num - frames[0]) // skip
                u[idx, :, :] = U[:, 0, :]
                v[idx, :, :] = V[:, 0, :]
                w[idx, :, :] = W[:, 0, :]

        return u, v, w

    if flag == 5 or flag == 6:
        u = np.zeros((total_frames, ngx + 2, ngy + 1))
        v = np.zeros((total_frames, ngx + 1, ngy + 2))
        w = np.zeros((total_frames, ngx + 1, ngy + 1))

        with concurrent.futures.ProcessPoolExecutor() as executor:
            future_to_frame = {
                executor.submit(read_boundary_file, loc, ngx, ngy, ngz, frame, flag, adjFlag, adjOpt): frame
                for frame in frames
            }

            for future in tqdm(concurrent.futures.as_completed(future_to_frame),
                               total=total_frames,
                               desc="Reading boundary data (flag 5/6)"):
                num, U, V, W = future.result()
                idx = (num - frames[0]) // skip
                u[idx, :, :] = U[0, :, :]
                v[idx, :, :] = V[0, :, :]
                w[idx, :, :] = W[0, :, :]

        return u, v, w

# =============================================================================
# TECPLOT WRITERS (with Numba-accelerated array preparation)
# =============================================================================

@njit(cache=True)
def _format_tecplot_data_I_plane(Z, Y, u, v, w, p):
    """
    Prepare I-plane data for Tecplot output.
    Returns flattened arrays in correct order.
    """
    nz, ny = Z.shape
    n_points = nz * ny
    data = np.empty((n_points, 6), dtype=np.float64)

    idx = 0
    for iz in range(nz):
        for iy in range(ny):
            data[idx, 0] = Z[iz, iy]
            data[idx, 1] = Y[iz, iy]
            data[idx, 2] = u[iz, iy]
            data[idx, 3] = v[iz, iy]
            data[idx, 4] = w[iz, iy]
            data[idx, 5] = p[iz, iy]
            idx += 1

    return data

@njit(cache=True)
def _format_tecplot_data_JK_plane(X, coord2, u, v, w, p):
    """
    Prepare J or K plane data for Tecplot output.
    """
    nx, n2 = X.shape
    n_points = nx * n2
    data = np.empty((n_points, 6), dtype=np.float64)

    idx = 0
    for ix in range(nx):
        for i2 in range(n2):
            data[idx, 0] = X[ix, i2]
            data[idx, 1] = coord2[ix, i2]
            data[idx, 2] = u[ix, i2]
            data[idx, 3] = v[ix, i2]
            data[idx, 4] = w[ix, i2]
            data[idx, 5] = p[ix, i2]
            idx += 1

    return data

def write_I_plane_data_ASCII_tecplot(file_loc, u, v, w, p, z, y, i_ind=0, t_ind=0, time=0.0):
    """
    Write a Z-Y plane slice (I-plane) of 3D flow data to ASCII Tecplot file.

    Parameters:
    -----------
    file_loc : str
        Directory to save the output file (end with '/')
    u, v, w : 2D ndarray
        Velocity components in the Z-Y plane (shape: [nz, ny])
    p : 2D ndarray
        Pressure values in the Z-Y plane (shape: [nz, ny])
    z, y : 1D ndarray
        Grid coordinates along z and y directions
    i_ind : int
        Index of the I-plane
    t_ind : int
        Time step index
    time : float
        Physical solution time
    """
    filename = file_loc + f"I_{i_ind:03d}_T_{t_ind:03d}.dat"
    Z, Y = np.meshgrid(z, y, indexing='ij')

    # Use Numba-accelerated data preparation
    data = _format_tecplot_data_I_plane(
        np.ascontiguousarray(Z), np.ascontiguousarray(Y),
        np.ascontiguousarray(u), np.ascontiguousarray(v),
        np.ascontiguousarray(w), np.ascontiguousarray(p)
    )

    with open(filename, "w") as f:
        f.write(f"TITLE = \"Data at i = {i_ind} \"\n")
        f.write("VARIABLES = \"Z\", \"Y\", \"U\", \"V\", \"W\", \"P\"\n")
        f.write(
            f"ZONE I={len(y)}, J={len(z)}, F=POINT, "
            f"STRANDID={i_ind}, SOLUTIONTIME={time:.2f}, T=\"t_{t_ind}\"\n"
        )
        for row in data:
            f.write(f"{row[0]:.6f} {row[1]:.6f} {row[2]:.6e} {row[3]:.6e} {row[4]:.6e} {row[5]:.6e}\n")

def write_J_plane_data_ASCII_tecplot(file_loc, u, v, w, p, x, z, j_ind=0, t_ind=0, time=0.0):
    """
    Write an X-Z plane slice (J-plane) of 3D flow data to ASCII Tecplot file.

    Parameters:
    -----------
    file_loc : str
        Directory to save the output file (end with '/')
    u, v, w : 2D ndarray
        Velocity components in the X-Z plane (shape: [nx, nz])
    p : 2D ndarray
        Pressure values in the X-Z plane (shape: [nx, nz])
    x, z : 1D ndarray
        Grid coordinates along x and z directions
    j_ind : int
        Index of the J-plane
    t_ind : int
        Time step index
    time : float
        Physical solution time
    """
    filename = file_loc + f"J_{j_ind:03d}_T_{t_ind:03d}.dat"
    X, Z = np.meshgrid(x, z, indexing='ij')

    data = _format_tecplot_data_JK_plane(
        np.ascontiguousarray(X), np.ascontiguousarray(Z),
        np.ascontiguousarray(u), np.ascontiguousarray(v),
        np.ascontiguousarray(w), np.ascontiguousarray(p)
    )

    with open(filename, "w") as f:
        f.write(f"TITLE = \"Data at j = {j_ind} \"\n")
        f.write("VARIABLES = \"X\", \"Z\", \"U\", \"V\", \"W\", \"P\"\n")
        f.write(
            f"ZONE I={len(z)}, J={len(x)}, F=POINT, "
            f"STRANDID={j_ind}, SOLUTIONTIME={time:.6f}, T=\"t_{t_ind}\"\n"
        )
        for row in data:
            f.write(f"{row[0]:.6f} {row[1]:.6f} {row[2]:.6e} {row[3]:.6e} {row[4]:.6e} {row[5]:.6e}\n")

def write_K_plane_data_ASCII_tecplot(file_loc, u, v, w, p, x, y, k_ind=0, t_ind=0, time=0.0):
    """
    Write an X-Y plane slice (K-plane) of 3D flow data to ASCII Tecplot file.

    Parameters:
    -----------
    file_loc : str
        Directory to save the output file (end with '/')
    u, v, w : 2D ndarray
        Velocity components in the X-Y plane (shape: [nx, ny])
    p : 2D ndarray
        Pressure values in the X-Y plane (shape: [nx, ny])
    x, y : 1D ndarray
        Grid coordinates along x and y directions
    k_ind : int
        Index of the K-plane
    t_ind : int
        Time step index
    time : float
        Physical solution time
    """
    filename = file_loc + f"K_{k_ind:03d}_T_{t_ind:03d}.dat"
    X, Y = np.meshgrid(x, y, indexing='ij')

    data = _format_tecplot_data_JK_plane(
        np.ascontiguousarray(X), np.ascontiguousarray(Y),
        np.ascontiguousarray(u), np.ascontiguousarray(v),
        np.ascontiguousarray(w), np.ascontiguousarray(p)
    )

    with open(filename, "w") as f:
        f.write(f"TITLE = \"Data at k = {k_ind} \"\n")
        f.write("VARIABLES = \"X\", \"Y\", \"U\", \"V\", \"W\", \"P\"\n")
        f.write(
            f"ZONE I={len(y)}, J={len(x)}, F=POINT, "
            f"STRANDID={k_ind}, SOLUTIONTIME={time:.6f}, T=\"t_{t_ind}\"\n"
        )
        for row in data:
            f.write(f"{row[0]:.6f} {row[1]:.6f} {row[2]:.6e} {row[3]:.6e} {row[4]:.6e} {row[5]:.6e}\n")

@njit(cache=True)
def _format_wall_data(X, Z, Sxy, Szy, Pw):
    """Prepare wall data for Tecplot output."""
    nx, nz = X.shape
    n_points = nx * nz
    data = np.empty((n_points, 5), dtype=np.float64)

    idx = 0
    for ix in range(nx):
        for iz in range(nz):
            data[idx, 0] = X[ix, iz]
            data[idx, 1] = Z[ix, iz]
            data[idx, 2] = Sxy[ix, iz]
            data[idx, 3] = Szy[ix, iz]
            data[idx, 4] = Pw[ix, iz]
            idx += 1

    return data

def write_wall_data_ASCII_tecplot(file_loc, Sxy, Szy, Pw, x, z, t_ind=0, time=0.0):
    """
    Write wall data to ASCII Tecplot file.

    Parameters:
    -----------
    file_loc : str
        Directory to save the output file
    Sxy, Szy : 2D ndarray
        Wall shear stress components
    Pw : 2D ndarray
        Wall pressure
    x, z : 1D ndarray
        Grid coordinates
    t_ind : int
        Time step index
    time : float
        Physical solution time
    """
    filename = file_loc + f"wall_T_{t_ind:03d}.dat"
    X, Z = np.meshgrid(x, z, indexing='ij')

    data = _format_wall_data(
        np.ascontiguousarray(X), np.ascontiguousarray(Z),
        np.ascontiguousarray(Sxy), np.ascontiguousarray(Szy),
        np.ascontiguousarray(Pw)
    )

    with open(filename, "w") as f:
        f.write(f"TITLE = \"Data at wall \"\n")
        f.write("VARIABLES = \"X\", \"Z\", \"Sxy\", \"Szy\",\"Pw\"\n")
        f.write(
            f"ZONE I={len(z)}, J={len(x)}, F=POINT, "
            f"STRANDID={1}, SOLUTIONTIME={time:.6f}, T=\"t_{t_ind}\"\n"
        )
        for row in data:
            f.write(f"{row[0]:.6f} {row[1]:.6f} {row[2]:.6e} {row[3]:.6e} {row[4]:.6e}\n")

# =============================================================================
# OBSERVATION WRITE FUNCTIONS
# =============================================================================

def write_obs_file(U, V, W, num=0, dir='./'):
    """
    Writes observation file for subdomain.

    Parameters:
    -----------
    U, V, W : ndarray
        Velocity components (should be in flux form)
    num : int
        Timestep number
    dir : str
        Output directory
    """
    U_out = _prepare_for_write(np.ascontiguousarray(U))
    V_out = _prepare_for_write(np.ascontiguousarray(V))
    W_out = _prepare_for_write(np.ascontiguousarray(W))

    fname = '04uobs.dble.' + "{:07d}".format(num)
    with open(dir + fname, 'wb') as fout:
        fout.write(U_out.tobytes())

    fname = '05vobs.dble.' + "{:07d}".format(num)
    with open(dir + fname, 'wb') as fout:
        fout.write(V_out.tobytes())

    fname = '06wobs.dble.' + "{:07d}".format(num)
    with open(dir + fname, 'wb') as fout:
        fout.write(W_out.tobytes())

def write_obs_fullsized_file(U, V, W, num=0, addr='./'):
    """
    Writes ucvcwc restart file at a given timestep.

    Parameters:
    -----------
    addr : str
        Name of the folder to write to
    U, V, W : ndarray
        Velocity components
    num : int
        Number of the file to write
    """
    U_out = _prepare_for_write(np.ascontiguousarray(U))
    V_out = _prepare_for_write(np.ascontiguousarray(V))
    W_out = _prepare_for_write(np.ascontiguousarray(W))

    fname = addr + '01ucmn.dble.' + "{:07d}".format(num)
    with open(fname, 'wb') as fout:
        fout.write(U_out.tobytes())

    fname = addr + '02vcmn.dble.' + "{:07d}".format(num)
    with open(fname, 'wb') as fout:
        fout.write(V_out.tobytes())

    fname = addr + '03wcmn.dble.' + "{:07d}".format(num)
    with open(fname, 'wb') as fout:
        fout.write(W_out.tobytes())

def read_obs_fullsized_file(fileroot, ngx, ngy, ngz, num=0, flag=1):
    """
    Read full-sized observation files.

    Parameters:
    -----------
    fileroot : str
        Location of data files
    ngx, ngy, ngz : int
        Grid dimensions
    num : int
        Timestep number
    flag : int
        1 for observation files, else for mean files

    Returns:
    --------
    tuple
        (Ufield, Vfield, Wfield)
    """
    fname = '01ucob.dble.' if flag == 1 else '01ucmn.dble.'
    filename = fileroot + fname + "{:07d}".format(num)
    Udata = np.fromfile(filename, dtype=np.float64)
    Ufield = _reshape_and_transpose_fortran_parallel(
        Udata, (ngy + 1, ngx + 2, ngz + 1)
    )

    fname = '02vcob.dble.' if flag == 1 else '02vcmn.dble.'
    filename = fileroot + fname + "{:07d}".format(num)
    Vdata = np.fromfile(filename, dtype=np.float64)
    Vfield = _reshape_and_transpose_fortran_parallel(
        Vdata, (ngy + 2, ngx + 1, ngz + 1)
    )

    fname = '03wcob.dble.' if flag == 1 else '03wcmn.dble.'
    filename = fileroot + fname + "{:07d}".format(num)
    Wdata = np.fromfile(filename, dtype=np.float64)
    Wfield = _reshape_and_transpose_fortran_parallel(
        Wdata, (ngy + 1, ngx + 1, ngz + 2)
    )

    return Ufield, Vfield, Wfield

# =============================================================================
# OBSERVATION KERNEL FUNCTIONS
# =============================================================================

def write_obsKern(Ku, Kv, Kw, addr='./'):
    """
    Write observation kernel files.

    Parameters:
    -----------
    Ku, Kv, Kw : ndarray
        Observation kernels
    addr : str
        Output directory
    """
    num = 0

    Ku_out = _prepare_for_write(np.ascontiguousarray(Ku))
    Kv_out = _prepare_for_write(np.ascontiguousarray(Kv))
    Kw_out = _prepare_for_write(np.ascontiguousarray(Kw))

    fname = addr + 'obKeru.dble.' + "{:07d}".format(num)
    with open(fname, 'wb') as fout:
        fout.write(Ku_out.tobytes())

    fname = addr + 'obKerv.dble.' + "{:07d}".format(num)
    with open(fname, 'wb') as fout:
        fout.write(Kv_out.tobytes())

    fname = addr + 'obKerw.dble.' + "{:07d}".format(num)
    with open(fname, 'wb') as fout:
        fout.write(Kw_out.tobytes())

    print(f'Observation kernel written to {addr}obKer*.dble.{num}')

def write_flags(Fu, Fv, Fw, addr='./'):
    """
    Write observation flag files.

    Parameters:
    -----------
    Fu, Fv, Fw : ndarray (int32)
        Observation flags
    addr : str
        Output directory
    """
    num = 0

    # Flags are integers, need different handling
    Fu_out = np.transpose(np.ascontiguousarray(Fu.astype(np.int32)))
    Fv_out = np.transpose(np.ascontiguousarray(Fv.astype(np.int32)))
    Fw_out = np.transpose(np.ascontiguousarray(Fw.astype(np.int32)))

    fname = addr + 'obFlgu.intg.' + "{:07d}".format(num)
    with open(fname, 'wb') as fout:
        fout.write(Fu_out.tobytes())

    fname = addr + 'obFlgv.intg.' + "{:07d}".format(num)
    with open(fname, 'wb') as fout:
        fout.write(Fv_out.tobytes())

    fname = addr + 'obFlgw.intg.' + "{:07d}".format(num)
    with open(fname, 'wb') as fout:
        fout.write(Fw_out.tobytes())

    print(f'Observation flags written to {addr}obFlg*.intg.{num}')

def read_obsKern(addr, ngx, ngy, ngz, num=0):
    """
    Read observation kernel files.

    Parameters:
    -----------
    addr : str
        Directory containing files
    ngx, ngy, ngz : int
        Grid dimensions
    num : int
        File number

    Returns:
    --------
    tuple
        (Ku, Kv, Kw)
    """
    fname = addr + 'obKeru.dble.' + "{:07d}".format(num)
    Ku_data = np.fromfile(fname, dtype=np.float64)
    Ku = _reshape_and_transpose_fortran_parallel(
        Ku_data, (ngy + 1, ngx + 2, ngz + 1)
    )

    fname = addr + 'obKerv.dble.' + "{:07d}".format(num)
    Kv_data = np.fromfile(fname, dtype=np.float64)
    Kv = _reshape_and_transpose_fortran_parallel(
        Kv_data, (ngy + 2, ngx + 1, ngz + 1)
    )

    fname = addr + 'obKerw.dble.' + "{:07d}".format(num)
    Kw_data = np.fromfile(fname, dtype=np.float64)
    Kw = _reshape_and_transpose_fortran_parallel(
        Kw_data, (ngy + 1, ngx + 1, ngz + 2)
    )

    print(f'Observation kernel read from {addr}obKer*.dble.{num}')
    return Ku, Kv, Kw

def read_flags(addr, ngx, ngy, ngz, num=0):
    """
    Read observation flag files.

    Parameters:
    -----------
    addr : str
        Directory containing files
    ngx, ngy, ngz : int
        Grid dimensions
    num : int
        File number

    Returns:
    --------
    tuple
        (uflag, vflag, wflag)
    """
    fname = addr + 'obFlgu.intg.' + "{:07d}".format(num)
    uflag_data = np.fromfile(fname, dtype=np.int32)
    uflag = np.transpose(np.reshape(uflag_data, (ngy + 1, ngx + 2, ngz + 1)))

    fname = addr + 'obFlgv.intg.' + "{:07d}".format(num)
    vflag_data = np.fromfile(fname, dtype=np.int32)
    vflag = np.transpose(np.reshape(vflag_data, (ngy + 2, ngx + 1, ngz + 1)))

    fname = addr + 'obFlgw.intg.' + "{:07d}".format(num)
    wflag_data = np.fromfile(fname, dtype=np.int32)
    wflag = np.transpose(np.reshape(wflag_data, (ngy + 1, ngx + 1, ngz + 2)))

    print(f'Observation flags read from {addr}obKer*.flag.{num}')
    return uflag, vflag, wflag

# =============================================================================
# STATISTICS FUNCTIONS
# =============================================================================

def read_DA_statistics(file_addr):
    """
    Read data assimilation statistics from text files.

    Parameters:
    -----------
    file_addr : str
        Directory containing statistics files

    Returns:
    --------
    tuple
        Statistics arrays (u_mean, v_mean, w_mean, p_mean, uu, vv, ww, uv, pp)
    """
    u_mean = np.loadtxt(file_addr + 'u_mean.txt')
    v_mean = np.loadtxt(file_addr + 'v_mean.txt')
    w_mean = np.loadtxt(file_addr + 'w_mean.txt')
    p_mean = np.loadtxt(file_addr + 'w_mean.txt')
    uu = np.loadtxt(file_addr + 'uu.txt')
    vv = np.loadtxt(file_addr + 'vv.txt')
    ww = np.loadtxt(file_addr + 'ww.txt')
    uv = np.loadtxt(file_addr + 'uv.txt')
    pp = np.loadtxt(file_addr + 'pp.txt')
    return u_mean, v_mean, w_mean, p_mean, uu, vv, ww, uv, pp

def read_exp_statistics(file_addr):
    """
    Read experimental statistics from files.

    Parameters:
    -----------
    file_addr : str
        Path to the directory containing the experimental statistics files

    Returns:
    --------
    tuple
        (u_mean_exp, v_mean_exp, w_mean_exp, uu_exp, vv_exp, ww_exp, uv_exp)
    """
    u_mean_exp = np.loadtxt(file_addr + 'u_mean.txt')
    v_mean_exp = np.loadtxt(file_addr + 'v_mean.txt')
    w_mean_exp = np.loadtxt(file_addr + 'w_mean.txt')
    uu_exp = np.loadtxt(file_addr + 'uu.txt')
    vv_exp = np.loadtxt(file_addr + 'vv.txt')
    ww_exp = np.loadtxt(file_addr + 'ww.txt')
    uv_exp = np.loadtxt(file_addr + 'uv.txt')

    return u_mean_exp, v_mean_exp, w_mean_exp, uu_exp, vv_exp, ww_exp, uv_exp

# =============================================================================
# PRE-COMPILATION FUNCTION
# =============================================================================

def precompile_numba_functions():
    """
    Pre-compile all Numba JIT functions with small dummy arrays.
    Call this once at module import or application startup to avoid
    compilation delays during actual use.

    Usage:
    ------
        from io_operations_Numba import precompile_numba_functions
        precompile_numba_functions()  # Call once at startup
    """
    # Small test arrays
    arr_3d = np.zeros((4, 4, 4), dtype=np.float64)
    arr_2d = np.zeros((4, 4), dtype=np.float64)
    arr_1d = np.zeros(64, dtype=np.float64)

    # Trigger compilation
    _ = _transpose_3d(arr_3d)
    _ = _transpose_2d(arr_2d)
    _ = _reshape_and_transpose_fortran(arr_1d, (4, 4, 4))
    _ = _reshape_and_transpose_fortran_parallel(arr_1d, (4, 4, 4))
    _ = _prepare_for_write(arr_3d)
    _ = _prepare_2d_for_write(arr_2d)

    # Boundary extraction functions
    U = np.zeros((5, 6, 5), dtype=np.float64)
    V = np.zeros((5, 5, 6), dtype=np.float64)
    W = np.zeros((6, 5, 5), dtype=np.float64)
    _ = _extract_boundary_bot_top(U, V, W, 1)
    _ = _extract_boundary_lft_rht(U, V, W, 3)
    _ = _extract_boundary_bck_frt(U, V, W, 5)

    # Tecplot formatting functions
    Z, Y = np.meshgrid(np.arange(4), np.arange(4), indexing='ij')
    _ = _format_tecplot_data_I_plane(Z.astype(np.float64), Y.astype(np.float64),
                                      arr_2d, arr_2d, arr_2d, arr_2d)
    _ = _format_tecplot_data_JK_plane(Z.astype(np.float64), Y.astype(np.float64),
                                       arr_2d, arr_2d, arr_2d, arr_2d)
    _ = _format_wall_data(Z.astype(np.float64), Y.astype(np.float64),
                          arr_2d, arr_2d, arr_2d)

    print("Numba functions pre-compiled successfully.")
