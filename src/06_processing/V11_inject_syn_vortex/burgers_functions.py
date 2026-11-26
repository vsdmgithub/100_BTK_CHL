# %% [markdown]
'''
> __Project__    : Injecting artificial vortices into channel flow
> __Author__     : Sugan Velazhagan
> __Date__       : 25/11/2026
> __Description__: Functions for the Burgers vortex packet (used by main.py).
>
>                  Local frame (a, s1, s2), right-handed, a = vortex axis:
>                      axis = x -> (x, y, z),  axis = y -> (y, z, x),  axis = z -> (z, x, y)
>                  with r^2 = s1^2 + s2^2, x_a = coordinate along the axis (relative to the center).
>
>                  Burgers vortex packet (t = 0):
>                      g(x_a)   = exp(-x_a^2 / (2 L_0^2))
>                      u_theta  = g(x_a) Gamma / (2 pi r) [1 - exp(-r^2 / delta^2)]
>                      omega_a  = g(x_a) Gamma / (pi delta^2) exp(-r^2 / delta^2)
>                      omega_r  = (x_a / L_0^2) g(x_a) Gamma / (2 pi r) [1 - exp(-r^2 / delta^2)]
>                  Vector potential (u = curl A):
>                      A_a = g(x_a) psi(r),  psi = -Gamma / (4 pi) [ln(xi) + E1(xi)],  xi = r^2 / delta^2
>
>                  Arrays are stored as [z, x, y] (with halos), as in the solver files.
'''
# %%
from framework import *
from scipy.special import exp1

    # --> Vortex axis -> global indices (x=0, y=1, z=2) of the local frame (a, s1, s2)
AXIS_MAP = {'x': (0, 1, 2), 'y': (1, 2, 0), 'z': (2, 0, 1)}

# %%
# ~ HELPER FUNCTIONS
# ~ ~~~~~~~~~~~~~~~~~~~~~
def get_local_coordinates(X, Y, Z, center, axis):
    """
    Global coordinates -> local coordinates (x_a, s1, s2) relative to the vortex center.
    """
    m = AXIS_MAP[axis]
    d = (X - center[0], Y - center[1], Z - center[2])
    return d[m[0]], d[m[1]], d[m[2]]

def get_global_components(F_loc, axis):
    """
    Local vector components (F_a, F_1, F_2) -> global components (F_x, F_y, F_z).
    """
    m = AXIS_MAP[axis]
    F = [None, None, None]
    for i in range(3):
        F[m[i]] = F_loc[i]
    return F

def get_swirl_factor(rr2, Gamma, spread):
    """
    q(r) = u_theta / r = Gamma / (2 pi r^2) [1 - exp(-r^2/delta^2)]   (regular at r = 0).
    """
    rr2_safe = np.maximum(rr2, 1e-30)
    q = Gamma / (2.0 * np.pi * rr2_safe) * (-np.expm1(-rr2_safe / spread**2))
    return np.where(rr2 > 1e-30, q, Gamma / (2.0 * np.pi * spread**2))

def get_stream_function(rr2, Gamma, spread):
    """
    psi(r) = -Gamma / (4 pi) [ln(xi) + E1(xi)],  xi = r^2 / delta^2   (d psi / dr = -u_theta).
    """
    xi = np.maximum(rr2 / spread**2, 1e-14)   #^ ln(xi) + E1(xi) -> -gamma_Euler as xi -> 0
    return -Gamma / (4.0 * np.pi) * (np.log(xi) + exp1(xi))

def get_window(rr2, win_factor, spread):
    """
    Radial window exp(-(r / W)^4), W = win_factor * spread.
    """
    return np.exp(-(rr2 / (win_factor * spread)**2)**2)

def get_wall_damping(Y, wall_damp):
    """
    f(y) = [1 - exp(-(y_b/delta_w)^2)] [1 - exp(-(y_t/delta_w)^2)],  f = df/dy = 0 at both walls.
    """
    return (-np.expm1(-((Y - y_wall_bot) / wall_damp)**2)) * (-np.expm1(-((y_wall_top - Y) / wall_damp)**2))

def get_vector_potential(X, Y, Z, center, Gamma, spread, length, axis, win_factor, wall_damp, y_wall_near,
                         flag_localize=True, flag_image=False):
    """
    Vector potential (Ax, Ay, Az) of the Burgers vortex packet, including:
        - the image vortex across the near wall (only if flag_image):  A_img(x) = -P A(P x),  P: y -> 2 y_w - y
        - the radial window and the wall damping (only if flag_localize)
    """
    def get_A_global(Yq):
        xa, s1, s2 = get_local_coordinates(X, Yq, Z, center, axis)
        g   = np.exp(-xa**2 / (2.0 * length**2))
        A_a = g * get_stream_function(s1**2 + s2**2, Gamma, spread)
        return get_global_components((A_a, 0.0, 0.0), axis)       #^ Scalar zeros: no memory for empty components

    A     = get_A_global(Y)
    if flag_image:
        A_img = get_A_global(2.0 * y_wall_near - Y)
        A     = [A[0] - A_img[0], A[1] + A_img[1], A[2] - A_img[2]]   #^ Pseudo-vector mirror: tangential A = 0 on the wall
        del A_img

    if flag_localize:
        _, s1, s2 = get_local_coordinates(X, Y, Z, center, axis)
        fac = get_window(s1**2 + s2**2, win_factor, spread) * get_wall_damping(Y, wall_damp)
        A   = [fac * A[i] for i in range(3)]
    shape = np.broadcast_shapes(np.shape(X), np.shape(Y), np.shape(Z))
    return [np.broadcast_to(A[i], shape) for i in range(3)]          #^ Zero components: read-only views, no memory

def get_curl_staggered(Ax, Ay, Az):
    """
    Discrete curl on the staggered grid (div of the result = 0 to round-off):
        Ax at (x_cc, y_fc, z_fc), Ay at (x_fc, y_cc, z_fc), Az at (x_fc, y_fc, z_cc)
        -> u at (x_fc, y_cc, z_cc), v at (x_cc, y_fc, z_cc), w at (x_cc, y_cc, z_fc)
    """
    dx = np.diff(xx_fc)[None, :, None]
    dy = np.diff(yy_fc)[None, None, :]
    dz = np.diff(zz_fc)[:, None, None]
    u = np.diff(Az, axis=2) / dy - np.diff(Ay, axis=0) / dz
    v = np.diff(Ax, axis=0) / dz - np.diff(Az, axis=1) / dx
    w = np.diff(Ay, axis=1) / dx - np.diff(Ax, axis=2) / dy
    return u, v, w

# %%
# ~ 1. VORTICITY OF A VELOCITY FIELD
# ~ ~~~~~~~~~~~~~~~~~~~~~
def get_interp_axis(F, r_src, r_tgt, axis):
    """
    1D linear interpolation of F along 'axis' from the points r_src to r_tgt
    (non-uniform grids; linear extrapolation outside r_src, e.g. into the halos).
    """
    i = np.clip(np.searchsorted(r_src, r_tgt) - 1, 0, len(r_src) - 2)
    t = (r_tgt - r_src[i]) / (r_src[i + 1] - r_src[i])
    shape = [1, 1, 1]; shape[axis] = -1
    t = t.reshape(shape)
    return np.take(F, i, axis=axis) * (1.0 - t) + np.take(F, i + 1, axis=axis) * t

def get_interp_3d(F, src, tgt):
    """
    Trilinear interpolation of F[z, x, y] from the grid src = (x, y, z) to tgt = (x, y, z).
    """
    F = get_interp_axis(F, src[2], tgt[2], 0)
    F = get_interp_axis(F, src[0], tgt[0], 1)
    return get_interp_axis(F, src[1], tgt[1], 2)

def get_vorticity(u, v, w):
    """
    Vorticity of a staggered velocity field [z, x, y].

    Natural (2nd-order, compact) finite-difference locations:
        omega_x = dw/dy - dv/dz  at (x_cc, yy, zz)   -> shape (ngz,   ngx+1, ngy)
        omega_y = du/dz - dw/dx  at (xx, y_cc, zz)   -> shape (ngz,   ngx,   ngy+1)
        omega_z = dv/dx - du/dy  at (xx, yy, z_cc)   -> shape (ngz+1, ngx,   ngy)
    then interpolated to the staggered velocity points:
        omega_x -> u-grid (xu, yu, zu),  omega_y -> v-grid (xv, yv, zv),  omega_z -> w-grid (xw, yw, zw)

    Args:
        u, v, w : (ngz+1, ngx+2, ngy+1), (ngz+1, ngx+1, ngy+2), (ngz+2, ngx+1, ngy+1)

    Returns:
        omx, omy, omz : same shapes as u, v, w
    """
    dx_c = np.diff(xx_cc)[None, :, None]      #^ Spacing between cell centers -> derivative on the nodes
    dy_c = np.diff(yy_cc)[None, None, :]
    dz_c = np.diff(zz_cc)[:, None, None]

        # --> Natural points (u: x_fc, y_cc, z_cc | v: x_cc, y_fc, z_cc | w: x_cc, y_cc, z_fc; x_fc[1:-1] = xx)
    omx_n = np.diff(w[1:-1], axis=2) / dy_c - np.diff(v[:, :, 1:-1], axis=0) / dz_c      #^ (x_cc, yy, zz)
    omy_n = np.diff(u[:, 1:-1], axis=0) / dz_c - np.diff(w[1:-1], axis=1) / dx_c          #^ (xx, y_cc, zz)
    omz_n = np.diff(v[:, :, 1:-1], axis=1) / dx_c - np.diff(u[:, 1:-1], axis=2) / dy_c    #^ (xx, yy, z_cc)

        # --> Interpolate to the staggered velocity points
    omx = get_interp_3d(omx_n, (xx_cc, yy, zz), (xu, yu, zu))
    omy = get_interp_3d(omy_n, (xx, yy_cc, zz), (xv, yv, zv))
    omz = get_interp_3d(omz_n, (xx, yy, zz_cc), (xw, yw, zw))
    return omx, omy, omz

# %%
# ~ 2. BURGERS VORTEX PACKET: VORTICITY
# ~ ~~~~~~~~~~~~~~~~~~~~~
def get_burgers_vorticity(center=(x_c, y_c, z_c), Gamma=Gamma, spread=spread, length=length, axis=axis,
                          flag_image=False, y_wall_near=y_wall_near):
    """
    Analytical vorticity of the Burgers vortex packet on the staggered velocity points
    (omega_x on the u-grid, omega_y on the v-grid, omega_z on the w-grid; same as get_vorticity).
        omega_a = g Gamma / (pi delta^2) exp(-r^2/delta^2)
        omega_r = (x_a / L_0^2) g u_theta^0(r)   ->  omega_1 = omega_r s1 / r,  omega_2 = omega_r s2 / r

    Args:
        center     : (x_c, y_c, z_c) of the vortex core
        flag_image : also add the image vortex across the near wall, omega_img(x) = -P omega(P x)
                     (no window / wall damping: these are only applied on the injected velocity)

    Returns:
        omx, omy, omz : shapes (ngz+1, ngx+2, ngy+1), (ngz+1, ngx+1, ngy+2), (ngz+2, ngx+1, ngy+1)
    """
    def get_omega_global(X, Y, Z):
        xa, s1, s2 = get_local_coordinates(X, Y, Z, center, axis)
        rr2  = s1**2 + s2**2
        g    = np.exp(-xa**2 / (2.0 * length**2))
        om_a = g * Gamma / (np.pi * spread**2) * np.exp(-rr2 / spread**2)
        om_s = xa / length**2 * g * get_swirl_factor(rr2, Gamma, spread)    #^ omega_r / r
        return get_global_components((om_a, om_s * s1, om_s * s2), axis)

    omega = []
    for n, (x_g, y_g, z_g) in enumerate([(xu, yu, zu), (xv, yv, zv), (xw, yw, zw)]):
        Z, X, Y = np.meshgrid(z_g, x_g, y_g, indexing='ij', sparse=True)   #^ Broadcasting grids (low memory)
        om_n = np.zeros((len(z_g), len(x_g), len(y_g)))
        om_n += get_omega_global(X, Y, Z)[n]
        if flag_image:
            sign = 1.0 if n == 1 else -1.0      #^ Pseudo-vector mirror: -P omega(P x)
            om_n += sign * get_omega_global(X, 2.0 * y_wall_near - Y, Z)[n]
        omega.append(om_n)
    return tuple(omega)

# %%
# ~ 3. BURGERS VORTEX PACKET: VELOCITY
# ~ ~~~~~~~~~~~~~~~~~~~~~
def get_burgers_velocity(center=(x_c, y_c, z_c), Gamma=Gamma, spread=spread, length=length, axis=axis,
                         win_factor=win_factor, wall_damp=wall_damp, y_wall_near=y_wall_near, flag_localize=True,
                         flag_image=flag_image):
    """
    Velocity of the Burgers vortex packet on the staggered simulation grid, u' = curl(A),
    with radial window, wall damping and optional image vortex (discretely divergence-free, u' = 0 at walls).
    flag_localize = False switches off the window and the wall damping, e.g. to compare curl(u')
    with get_burgers_vorticity(flag_image=flag_image).

    Returns:
        u, v, w : staggered velocities [z, x, y]:
                  (ngz+1, ngx+2, ngy+1), (ngz+1, ngx+1, ngy+2), (ngz+2, ngx+1, ngy+1)
    """
    prm = (center, Gamma, spread, length, axis, win_factor, wall_damp, y_wall_near, flag_localize, flag_image)

        # --> Each component of A on its own cell edge
        #^ Sparse (broadcasting) grids keep the memory to a few full-size arrays
    Z, X, Y = np.meshgrid(zz_fc, xx_cc, yy_fc, indexing='ij', sparse=True)
    Ax = get_vector_potential(X, Y, Z, *prm)[0]
    Z, X, Y = np.meshgrid(zz_fc, xx_fc, yy_cc, indexing='ij', sparse=True)
    Ay = get_vector_potential(X, Y, Z, *prm)[1]
    Z, X, Y = np.meshgrid(zz_cc, xx_fc, yy_fc, indexing='ij', sparse=True)
    Az = get_vector_potential(X, Y, Z, *prm)[2]

    u = np.zeros((ngz + 1, ngx + 2, ngy + 1))
    v = np.zeros((ngz + 1, ngx + 1, ngy + 2))
    w = np.zeros((ngz + 2, ngx + 1, ngy + 1))
    du, dv, dw = get_curl_staggered(Ax, Ay, Az)
    del Ax, Ay, Az
    u += du; v += dv; w += dw
    return u, v, w

# %%
# ~ 3b. MEAN PROFILE (log law)
# ~ ~~~~~~~~~~~~~~~~~~~~~
def get_mean_profile(y, kappa=0.41, B=5.2):
    """
    Mean streamwise velocity U(y) (outer units) from Reichardt's law (viscous sublayer -> buffer -> log law):
        U+ = ln(1 + kappa y+)/kappa + (B - ln(kappa)/kappa) [1 - exp(-y+/11) - (y+/11) exp(-y+/3)]
    y+ is the distance to the nearest wall; odd extension into the ghost cells (U = 0 midway, at the wall).
    """
    d  = np.minimum(y - y_wall_bot, y_wall_top - y) / delta_nu
    dp = np.abs(d)
    Up = np.log1p(kappa * dp) / kappa + (B - np.log(kappa) / kappa) * (1.0 - np.exp(-dp / 11.0) - dp / 11.0 * np.exp(-dp / 3.0))
    return np.sign(d) * Up * u_tau

# %%
# ~ 4. DIAGNOSTICS
# ~ ~~~~~~~~~~~~~~~~~~~~~
def get_divergence(u, v, w):
    """
    Discrete divergence of a staggered velocity field, at the cell centers (with halos).
    """
    dx = np.diff(xx_fc)[None, :, None]
    dy = np.diff(yy_fc)[None, None, :]
    dz = np.diff(zz_fc)[:, None, None]
    return np.diff(u, axis=1) / dx + np.diff(v, axis=2) / dy + np.diff(w, axis=0) / dz

# %%
# ~ 5. 3D VISUALIZATION (PyVista)
# ~ ~~~~~~~~~~~~~~~~~~~~~
    # --> Default camera view (applied in this order: preset -> view_up -> azimuth/elevation/roll -> zoom)
VIEW_DEFAULT = dict(
    preset    = 'xy',        #^ Starting view: 'xy' (look down z), 'xz' (look down y), 'yz' (look down x), 'iso' (oblique)
    up        = (0, 1, 0),   #^ Direction pointing up in the image: (0, 1, 0) -> y+ up (wall at the bottom)
    azimuth   = 30.0,        #^ Rotate about the up direction [deg]; > 0 swings the camera to the right
    elevation = 20.0,        #^ Tilt up (> 0, look from above) or down (< 0) [deg]
    roll      = 0.0,         #^ Spin the image in its own plane [deg]
    zoom      = 1.0,         #^ > 1 zoom in, < 1 zoom out
    exact     = None)        #^ Exact camera [(eye), (focus), (up)] in wall units; overrides all of the above

def get_cell_centered(omx, omy, omz):
    """
    Staggered vorticity (omega_x on u-grid, omega_y on v-grid, omega_z on w-grid) -> cell centers (with halos).
    """
    return (0.5 * (omx[:, 1:, :] + omx[:, :-1, :]),
            0.5 * (omy[:, :, 1:] + omy[:, :, :-1]),
            0.5 * (omz[1:] + omz[:-1]))

def get_crop_indices(r_cc, r_c, half):
    """
    Indices of the cell centers within r_c +/- half.
    """
    return np.where(np.abs(r_cc - r_c) <= half)[0]

def plot_vorticity_iso(om, center, title, fig_name=None, crop_fac=3.0, iso_frac=0.3, y_excl_plus=0.0, view=VIEW_DEFAULT,
                       half=None):
    """
    Live PyVista view of the |omega| iso-surface colored by omega_z, cropped around the vortex center
    (coordinates in wall units). The last view (as left in the window) is exported as PNG + PDF.

    Args:
        om       : staggered vorticity (omega_x on u-grid, omega_y on v-grid, omega_z on w-grid)
        center   : (x_c, y_c, z_c) of the vortex core
        title    : text shown on the top-left
        fig_name : path without extension for the PNG/PDF export (None: no export)
        crop_fac : crop half-width = crop_fac * spread (transverse), crop_fac * length (axial)
        iso_frac : initial iso-level, as a fraction of max|omega| (slider)
        y_excl_plus : exclude the region within y+ < y_excl_plus of both walls from the crop
                      (the wall shear layers otherwise set max|omega|)
        view     : camera view (see VIEW_DEFAULT)
        half     : crop half-widths (hx, hy, hz); overrides crop_fac (e.g. (Lx, Ly, Lz) -> entire domain)
    """
    import pyvista as pv

        # --> Vorticity at the cell centers, cropped around the vortex center (clipped at the walls)
    om_c = get_cell_centered(*om)
    if half is None:
        half = [crop_fac * spread] * 3
        half[AXIS_MAP[axis][0]] = crop_fac * length
    i_cr = get_crop_indices(xx_cc, center[0], half[0])
    j_cr = get_crop_indices(yy_cc, center[1], half[1])
    y_ex = y_excl_plus * delta_nu
    j_cr = j_cr[(yy_cc[j_cr] >= y_wall_bot + y_ex) & (yy_cc[j_cr] <= y_wall_top - y_ex)]
    k_cr = get_crop_indices(zz_cc, center[2], half[2])
    sub  = np.ix_(k_cr, i_cr, j_cr)
    om_mag = np.sqrt(om_c[0][sub]**2 + om_c[1][sub]**2 + om_c[2][sub]**2)
    om_z   = om_c[2][sub]

        # --> Rectilinear grid in wall units (x+, y+ from the bottom wall, z+), fields as point data (x fastest)
    grid = pv.RectilinearGrid(xx_cc[i_cr] / delta_nu, (yy_cc[j_cr] - y_wall_bot) / delta_nu, zz_cc[k_cr] / delta_nu)
    grid.point_data['omega_mag'] = np.transpose(om_mag, (1, 2, 0)).ravel(order='F')
    grid.point_data['omega_z']   = np.transpose(om_z,   (1, 2, 0)).ravel(order='F')
    log.info(f"Crop {len(i_cr)} x {len(j_cr)} x {len(k_cr)} points, max|omega| = {om_mag.max():.4e}")

    om_max  = float(om_mag.max())
    omz_lim = float(np.max(np.abs(om_z))) + 1e-30

    pl = pv.Plotter(window_size=(1200, 800))
    pl.add_mesh(grid.outline(), color='k')

    def update_iso(frac):
        """
        Redraw the |omega| iso-surface at frac * max|omega|, colored by omega_z.
        """
        iso = grid.contour([frac * om_max], scalars='omega_mag')
        if iso.n_points == 0:
            pl.remove_actor('iso')
            return
        pl.add_mesh(iso, scalars='omega_z', cmap='RdBu_r', clim=(-omz_lim, omz_lim), smooth_shading=True,
                    name='iso', scalar_bar_args={'title': r'$\omega_z$', 'vertical': True, 'width': 0.05, 'height': 0.35,
                                                 'position_x': 0.88, 'position_y': 0.30, 'fmt': '%.2f',
                                                 'title_font_size': 30, 'label_font_size': 21})

    pl.add_slider_widget(update_iso, rng=(0.05, 0.95), value=iso_frac,
                         title=r'$\Vert\omega\Vert \,/\, \Vert\omega\Vert_{max}$',   #^ No '|' in VTK math text
                         pointa=(0.6, 0.92), pointb=(0.95, 0.92), style='modern', title_height=0.045)
    pl.add_text(title, font_size=15)
    pl.show_grid(xtitle='x+', ytitle='y+', ztitle='z+', fmt='%.0f', font_size=27)   #^ Integer ticks in wall units

        # --> Camera
    if view['exact'] is not None:
        pl.camera_position = view['exact']
    else:
        pl.camera_position = view['preset']
        pl.set_viewup(view['up'])
        pl.camera.azimuth   = view['azimuth']
        pl.camera.elevation = view['elevation']
        pl.camera.roll      = pl.camera.roll + view['roll']
        pl.camera.zoom(view['zoom'])
    pl.show(auto_close=False)          #^ Keep the scene after the window is closed, to export it

        # --> Export the last view: high-resolution PNG + vector PDF
    if fig_name is not None:
        pl.screenshot(fig_name + '.png', scale=2)
        pl.save_graphic(fig_name + '.pdf')
        log.info(f"Figures saved: {fig_name}.png / .pdf")
    pl.close()

# %%
