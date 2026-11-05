# %%
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ MODULES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import os
import numpy as np
import yaml
from scipy.interpolate import interpn
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ FUNCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def get_grid_wHalos(rr):
    return np.concatenate([[2 * rr[0] - rr[1]], rr, [2 * rr[-1] - rr[-2]]])

def get_grid_center_wHalos(rr):
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

def read_grid_data(fname,ngx,ngy,ngz,Lz):
    XX = np.fromfile(fname,dtype=float,count=ngx*ngy,offset=0)
    XX = np.transpose(np.reshape(XX,(ngy,ngx)))
    xx = XX[:,0]

    YY = np.fromfile(fname,dtype=float,count=ngx*ngy,offset=ngx*ngy*8)
    YY = np.transpose(np.reshape(YY,(ngy,ngx)))
    yy = YY[0,:]

    zz = np.linspace(0,Lz,ngz)   # uniform periodic z grid
    return(XX,YY,xx,yy,zz)

def read_ucvcwc_file(fileroot,ngx,ngy,ngz,num,adjointFlag=0,component=-1):
    fname =  '01uc00.dble.'
    if (adjointFlag == 1):
        fname =  '01ucad.dble.'

    filename = fileroot + fname + "{:07d}".format(num)
    Ufield =  np.fromfile(filename,dtype=float,count=(ngx+2)*(ngz+1)*(ngy+1),offset=0)
    Ufield =  np.transpose(np.reshape(Ufield,(ngy+1,ngx+2,ngz+1)))
    if component == 0:
        return num,Ufield

    fname = '02vc00.dble.'
    if (adjointFlag == 1):
        fname = '02vcad.dble.'

    filename = fileroot + fname + "{:07d}".format(num)
    Vfield =  np.fromfile(filename,dtype=float,count=(ngx+1)*(ngz+1)*(ngy+2),offset=0)
    Vfield =  np.transpose(np.reshape(Vfield,(ngy+2,ngx+1,ngz+1)))
    if component == 1:
        return num,Vfield

    fname = '03wc00.dble.'
    if (adjointFlag == 1):
        fname = '03wcad.dble.'

    filename = fileroot + fname + "{:07d}".format(num)
    Wfield =  np.fromfile(filename,dtype=float,count=(ngx+1)*(ngz+2)*(ngy+1),offset=0)
    Wfield =  np.transpose(np.reshape(Wfield,(ngy+1,ngx+1,ngz+2)))
    if component == 2:
        return num,Wfield

    return (num,Ufield,Vfield,Wfield)

def write_ucvcwc_file(U,V,W,num=0,addr='./'):
    """
    Writes ucvcwc restart file at a given timestep.
    Fields are [z,x,y]; the main code expects y,x,z (z fastest) in the file.
    """
    for prefix, F in (('01uc00.dble.', U), ('02vc00.dble.', V), ('03wc00.dble.', W)):
        np.ascontiguousarray(np.transpose(F)).tofile(addr + prefix + "{:07d}".format(num))

def interpolate_field_3d(field_src, grid_src, grid_tgt, method="linear"):
    """
    Interpolate a [z,x,y] field from the source staggered grid to the target one.
    grid_* = (x, y, z) arrays as returned by get_staggered_grid.
    Target halo points outside the source range are linearly extrapolated.
    """
    xs, ys, zs = grid_src
    xt, yt, zt = grid_tgt
    ZT, XT, YT = np.meshgrid(zt, xt, yt, indexing="ij")
    pts = np.column_stack((ZT.ravel(), XT.ravel(), YT.ravel()))
    out = interpn((zs, xs, ys), field_src, pts, method=method, bounds_error=False, fill_value=None)
    return out.reshape(len(zt), len(xt), len(yt))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ LOGISTICS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
here     = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.normpath(os.path.join(here, '../../data'))
dir_prod = os.path.join(data_dir, 'sim_package_production') + '/'   # source (fine)
dir_test = os.path.join(data_dir, 'sim_package_testing') + '/'      # target (coarse)

def load_params(folder):
    with open(folder + 'grid_parameters.yaml') as f:
        return yaml.safe_load(f)

p_fine = load_params(dir_prod)
p_crse = load_params(dir_test)

ngx_fine, ngy_fine, ngz_fine = p_fine['ngx'], p_fine['ngy'], p_fine['ngz']
ngx_crse, ngy_crse, ngz_crse = p_crse['ngx'], p_crse['ngy'], p_crse['ngz']
Lz_fine, Lz_crse = p_fine['Lz'], p_crse['Lz']

assert np.isclose(p_fine['Lx'], p_crse['Lx']) and np.isclose(Lz_fine, Lz_crse) and np.isclose(p_fine['Ly'], p_crse['Ly']), \
    "production and testing domains must be identical"

print(f"fine  (production): {ngx_fine} x {ngy_fine} x {ngz_fine}")
print(f"coarse (testing)  : {ngx_crse} x {ngy_crse} x {ngz_crse}")

#^ Read grids
XX_fine, YY_fine, xx_fine, yy_fine, zz_fine = read_grid_data(dir_prod + 'grid.data', ngx_fine, ngy_fine, ngz_fine, Lz_fine)
XX_crse, YY_crse, xx_crse, yy_crse, zz_crse = read_grid_data(dir_test + 'grid.data', ngx_crse, ngy_crse, ngz_crse, Lz_crse)

grid_u_fine = get_staggered_grid(xx_fine, yy_fine, zz_fine, 1)
grid_v_fine = get_staggered_grid(xx_fine, yy_fine, zz_fine, 2)
grid_w_fine = get_staggered_grid(xx_fine, yy_fine, zz_fine, 3)

grid_u_crse = get_staggered_grid(xx_crse, yy_crse, zz_crse, 1)
grid_v_crse = get_staggered_grid(xx_crse, yy_crse, zz_crse, 2)
grid_w_crse = get_staggered_grid(xx_crse, yy_crse, zz_crse, 3)

#^ Read production velocity
num = 0
_, u_fine, v_fine, w_fine = read_ucvcwc_file(dir_prod, ngx_fine, ngy_fine, ngz_fine, num, adjointFlag=0, component=-1)

#^ Interpolate fine (production) to coarse (testing) grid
interp_method = 'linear'
u_crse = interpolate_field_3d(u_fine, grid_u_fine, grid_u_crse, method=interp_method)
v_crse = interpolate_field_3d(v_fine, grid_v_fine, grid_v_crse, method=interp_method)
w_crse = interpolate_field_3d(w_fine, grid_w_fine, grid_w_crse, method=interp_method)

#^ Write testing velocity
write_ucvcwc_file(u_crse, v_crse, w_crse, num=num, addr=dir_test)
print("written 01uc00 / 02vc00 / 03wc00 to", dir_test)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ PLOTTING
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Wall units assume Re = Re_bulk (nu = 1/Re_bulk, h = 1) and the Re_tau in grid_parameters.yaml
Re_tau = p_fine['Re_tau']
delta_nu = 1.0 / Re_tau
u_tau = Re_tau / p_fine['Re_bulk']

k_mid_fine = ngz_fine // 2
k_mid_crse = ngz_crse // 2

def fluc(u):
    return u - u.mean(axis=(0, 1))[None, None, :]

def plot_uv_w_planes(u_plane, v_plane, w_plane, g_u, g_v, g_w, out_path):
    fig, axes = plt.subplots(3, 1, figsize=(6, 10), constrained_layout=True)
    n_lev = 201
    items = ((u_plane, g_u, r"$u'^+$", 4), (v_plane, g_v, r"$v^+$", 2), (w_plane, g_w, r"$w^+$", 2))
    for ax, (plane, (x, y, _), title, lim) in zip(axes, items):
        X, Y = np.meshgrid(x / delta_nu, y / delta_nu, indexing="ij")
        im = ax.contourf(X, Y, plane / u_tau, levels=np.linspace(-lim, lim, n_lev), cmap="RdBu", vmin=-lim, vmax=lim, extend="both")
        ax.set_ylabel(r"$y^+$")
        ax.set_title(title)
        fig.colorbar(im, ax=ax, orientation="vertical")
    axes[-1].set_xlabel(r"$x^+$")
    fig.savefig(out_path, dpi=200)
    plt.close(fig)

plot_uv_w_planes(fluc(u_fine)[k_mid_fine], v_fine[k_mid_fine], w_fine[k_mid_fine],
                 grid_u_fine, grid_v_fine, grid_w_fine, os.path.join(here, "uvw_xy_kmid_fine.png"))
plot_uv_w_planes(fluc(u_crse)[k_mid_crse], v_crse[k_mid_crse], w_crse[k_mid_crse],
                 grid_u_crse, grid_v_crse, grid_w_crse, os.path.join(here, "uvw_xy_kmid_coarse.png"))
