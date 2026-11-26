# %% [markdown]
'''
> __Project__    : Injecting artificial vortices into channel flow
> __Author__     : Sugan Velazhagan
> __Date__       : 26/11/2026
> __Description__: Functions for the hairpin vortex with a Burgers core (used by main.py, viz_hairpin.py).
>                  Successor of V11_inject_syn_vortex/burgers_functions.py (doc/Analytical_vortices_p2.md, Sec. 1-4).
>
>                  Centreline (Gaussian hairpin, head at X_h = (x_c, y_c, z_c), s_w = +1 bottom / -1 top wall):
>                      X(q)  = X_h - (H sin(alpha), s_w H cos(alpha), 0) [1 - exp(-q^2/2)] + (0, 0, B q)
>                      X'(q) = (-H sin(alpha) q e^{-q^2/2}, -s_w H cos(alpha) q e^{-q^2/2}, B)
>                      kappa_0 = H / B^2 (head curvature), T = X' / |X'|  (T = +z at the head)
>                  Burgers core along the centreline, circulation tapered on the legs:
>                      omega = Gamma g(q) / (pi delta^2) exp(-r^2 / delta^2) T,   g(q) = exp(-(q / q_m)^6)
>                  realized as a line of 3D Gaussian vortex blobs (-> the above for a slender tube), each screened by a
>                  broad Gaussian of radius lambda = W = win_factor * delta (opposite sign, zero net circulation far away):
>                      omega(x) = Gamma int g(q) X'(q) [G_delta(rho) - G_lambda(rho)] dq,  rho = |x - X(q)|
>                      G_s(rho) = exp(-rho^2 / s^2) / (pi^{3/2} s^3)
>                  Vector potential (screened Biot-Savart; local, continuous, decays like erfc(rho/lambda)/rho):
>                      A(x) = Gamma / (4 pi) int g(q) X'(q) [erf(rho / delta) - erf(rho / lambda)] / rho dq
>                  u = curl A on the staggered grid (discretely div-free), A damped at the walls (u' = 0 there).
>                  curl(u) = omega + grad(div A): the tapered tube is projected onto a solenoidal field (leg ends).
>
>                  Arrays are stored as [z, x, y] (with halos), as in the solver files.
'''
# %%
from framework import *
from scipy.special import erf

# %%
# ~ HELPER FUNCTIONS
# ~ ~~~~~~~~~~~~~~~~~~~~~
def get_centreline(q, center=(x_c, y_c, z_c), height=height, span=span, incline=incline):
    """
    Hairpin centreline X(q) and its derivative X'(q), each of shape (3, len(q)).
    """
    al = np.deg2rad(incline)
    e  = np.exp(-0.5 * q**2)
    X  = np.array([center[0] - height * np.sin(al) * (1.0 - e),
                   center[1] - s_wall * height * np.cos(al) * (1.0 - e),
                   center[2] + span * q])
    dX = np.array([-height * np.sin(al) * q * e,
                   -s_wall * height * np.cos(al) * q * e,
                   span * np.ones_like(q)])
    return X, dX

def get_taper(q, q_max=q_max):
    """
    Circulation taper along the legs, g(q) = exp(-(q / q_m)^6)  (g = 1 on the head, -> 0 beyond q_m).
    """
    return np.exp(-(q / q_max)**6)

def get_samples(center, Gamma, spread, height, span, incline, q_max, flag_image, y_wall_near):
    """
    Quadrature of the centreline (midpoint rule in q, arc-length step <= delta / 4):
        P : sample points (n, 3);  Wq : weights Gamma g(q) X'(q) dq (n, 3);  P_real : points of the real hairpin.
    Image hairpin (flag_image): points P y -> 2 y_w - y, weights (-W_x, W_y, -W_z)  (A_img(x) = -P A(P x)).
    """
    q_end = 1.6 * q_max                                                  #^ g(q_end) ~ 4e-8
    speed = np.sqrt(span**2 + 0.37 * height**2)                          #^ max |X'(q)| (q e^{-q^2/2} <= 0.61)
    nq    = int(np.ceil(2.0 * q_end * speed / (0.25 * spread)))
    dq    = 2.0 * q_end / nq
    q     = -q_end + (np.arange(nq) + 0.5) * dq
    X, dX = get_centreline(q, center, height, span, incline)
    P, Wq = X.T, (Gamma * get_taper(q, q_max) * dq * dX).T
    P_real = P[get_taper(q, q_max) > 1e-3]
    if flag_image:
        P_img  = P * np.array([1.0, -1.0, 1.0]) + np.array([0.0, 2.0 * y_wall_near, 0.0])
        W_img  = Wq * np.array([-1.0, 1.0, -1.0])
        P, Wq  = np.vstack([P, P_img]), np.vstack([Wq, W_img])
    return P, Wq, P_real

def get_kernel_potential(rr2, spread, shield):
    """
    Potential of a unit screened blob, [erf(rho / delta) - erf(rho / lambda)] / rho  (shield = lambda; None: no screening)
    (-> 2 / sqrt(pi) (1 / delta - 1 / lambda) as rho -> 0).
    """
    rho = np.maximum(np.sqrt(rr2), 1e-12)
    if shield is None:
        return np.where(rho > 1e-12, erf(rho / spread) / rho, 2.0 / (np.sqrt(np.pi) * spread))
    return np.where(rho > 1e-12, (erf(rho / spread) - erf(rho / shield)) / rho,
                    2.0 / np.sqrt(np.pi) * (1.0 / spread - 1.0 / shield))

def get_kernel_vorticity(rr2, spread, shield):
    """
    Unit screened blob G_delta(rho) - G_lambda(rho),  G_s = exp(-rho^2 / s^2) / (pi^{3/2} s^3)  (shield None: G_delta).
    """
    G = np.exp(-rr2 / spread**2) / (np.pi**1.5 * spread**3)
    return G if shield is None else G - np.exp(-rr2 / shield**2) / (np.pi**1.5 * shield**3)

def get_distance2(X, Y, Z, P_real):
    """
    Squared distance from the grid points to the centreline (minimum over its samples).
    """
    d2 = np.full(np.broadcast_shapes(np.shape(X), np.shape(Y), np.shape(Z)), np.inf)
    for p in P_real:
        d2 = np.minimum(d2, (X - p[0])**2 + (Y - p[1])**2 + (Z - p[2])**2)
    return d2

def get_wall_damping(Y, wall_damp):
    """
    f(y) = [1 - exp(-(y_b/delta_w)^2)] [1 - exp(-(y_t/delta_w)^2)],  f = df/dy = 0 at both walls.
    """
    return (-np.expm1(-((Y - y_wall_bot) / wall_damp)**2)) * (-np.expm1(-((y_wall_top - Y) / wall_damp)**2))

def get_box_indices(r, lo, hi):
    """
    Index slice of the (sorted) coordinates r within [lo, hi].
    """
    i = np.where((r >= lo) & (r <= hi))[0]
    return slice(i[0], i[-1] + 1) if len(i) else slice(0, 0)

def get_hairpin_extent(margin=0.0, center=(x_c, y_c, z_c), spread=spread, height=height, span=span, incline=incline,
                       q_max=q_max):
    """
    Bounding box of the hairpin centreline (|q| <= 1.6 q_m) plus a margin: returns (center, half-widths).
    """
    _, _, P_real = get_samples(center, 1.0, spread, height, span, incline, q_max, False, y_wall_near)
    lo, hi = P_real.min(axis=0) - margin, P_real.max(axis=0) + margin
    return tuple(0.5 * (lo + hi)), tuple(0.5 * (hi - lo))

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
# ~ 2. HAIRPIN: VORTICITY (line of Gaussian blobs = Burgers core along the centreline)
# ~ ~~~~~~~~~~~~~~~~~~~~~
def get_hairpin_vorticity(center=(x_c, y_c, z_c), Gamma=Gamma, spread=spread, height=height, span=span,
                          incline=incline, q_max=q_max, win_factor=win_factor, flag_image=False,
                          y_wall_near=y_wall_near, flag_localize=True):
    """
    Analytical vorticity of the hairpin on the staggered simulation grid (no wall damping):
        omega(x) = sum_k Gamma g(q_k) X'(q_k) dq [G_delta(rho_k) - G_lambda(rho_k)],  lambda = win_factor * delta
        -> Gamma g / (pi delta^2) exp(-r^2 / delta^2) T (+ shield)  for a slender tube (kappa delta << 1).
    flag_image    : also add the image hairpin across the near wall
    flag_localize : False -> no shield (unscreened blobs)

    Returns:
        omx, omy, omz : omega_x on the u-grid, omega_y on the v-grid, omega_z on the w-grid [z, x, y]
    """
    P, Wq, _ = get_samples(center, Gamma, spread, height, span, incline, q_max, flag_image, y_wall_near)
    shield = win_factor * spread if flag_localize else None
    om = []
    for n, (xg, yg, zg) in enumerate([(xu, yu, zu), (xv, yv, zv), (xw, yw, zw)]):
        Z, X, Y = np.meshgrid(zg, xg, yg, indexing='ij', sparse=True)
        F = np.zeros((len(zg), len(xg), len(yg)))
        for p, wq in zip(P, Wq):
            F += wq[n] * get_kernel_vorticity((X - p[0])**2 + (Y - p[1])**2 + (Z - p[2])**2, spread, shield)
        om.append(F)
    return tuple(om)

# %%
# ~ 3. HAIRPIN: VELOCITY
# ~ ~~~~~~~~~~~~~~~~~~~~~
def get_vector_potential(n, grid, P, Wq, P_real, spread, shield, wall_damp, flag_localize=True):
    """
    Component n of the vector potential on the 1D grid = (x, y, z), shape [z, x, y]:
        A_n(x) = 1 / (4 pi) sum_k W_k,n [erf(rho_k / delta) - erf(rho_k / lambda)] / rho_k,   times the wall damping.
    With flag_localize, A is computed only in the box around the hairpin (A ~ erfc(rho / lambda) / rho outside);
    without it: unscreened potential over the whole domain, no wall damping.
    """
    xg, yg, zg = grid
    A = np.zeros((len(zg), len(xg), len(yg)))
    if flag_localize:
        W  = 2.5 * shield                                                #^ erfc(2.5) ~ 4e-4
        lo, hi = P_real.min(axis=0) - W, P_real.max(axis=0) + W
        sl = (get_box_indices(zg, lo[2], hi[2]), get_box_indices(xg, lo[0], hi[0]), get_box_indices(yg, lo[1], hi[1]))
    else:
        sl, shield = (slice(None), slice(None), slice(None)), None
    Z, X, Y = np.meshgrid(zg[sl[0]], xg[sl[1]], yg[sl[2]], indexing='ij', sparse=True)
    F = np.zeros(np.broadcast_shapes(X.shape, Y.shape, Z.shape))
    for p, wq in zip(P, Wq):
        F += wq[n] * get_kernel_potential((X - p[0])**2 + (Y - p[1])**2 + (Z - p[2])**2, spread, shield)
    F /= 4.0 * np.pi
    if flag_localize:
        F *= get_wall_damping(Y, wall_damp)
    A[sl] = F
    return A

def get_hairpin_velocity(center=(x_c, y_c, z_c), Gamma=Gamma, spread=spread, height=height, span=span,
                         incline=incline, q_max=q_max, win_factor=win_factor, wall_damp=wall_damp,
                         y_wall_near=y_wall_near, flag_localize=True, flag_image=flag_image):
    """
    Velocity of the hairpin on the staggered simulation grid, u' = curl(A),
    with screening (shield radius win_factor * spread), wall damping and optional image hairpin
    (discretely divergence-free, u' = 0 at walls).
    flag_localize = False switches off the screening and the wall damping (A over the whole domain), e.g. to
    compare curl(u') with get_hairpin_vorticity(flag_localize=False).

    Returns:
        u, v, w : staggered velocities [z, x, y]:
                  (ngz+1, ngx+2, ngy+1), (ngz+1, ngx+1, ngy+2), (ngz+2, ngx+1, ngy+1)
    """
    P, Wq, P_real = get_samples(center, Gamma, spread, height, span, incline, q_max, flag_image, y_wall_near)
    prm = (P, Wq, P_real, spread, win_factor * spread, wall_damp, flag_localize)

        # --> Each component of A on its own cell edge (only the box around the hairpin, if flag_localize)
    Ax = get_vector_potential(0, (xx_cc, yy_fc, zz_fc), *prm)
    Ay = get_vector_potential(1, (xx_fc, yy_cc, zz_fc), *prm)
    Az = get_vector_potential(2, (xx_fc, yy_fc, zz_cc), *prm)

    u = np.zeros((ngz + 1, ngx + 2, ngy + 1))
    v = np.zeros((ngz + 1, ngx + 1, ngy + 2))
    w = np.zeros((ngz + 2, ngx + 1, ngy + 1))
    du, dv, dw = get_curl_staggered(Ax, Ay, Az)
    del Ax, Ay, Az
    u += du; v += dv; w += dw
    return u, v, w

# %%
# ~ 4. MEAN PROFILE (log law)
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
# ~ 5. DIAGNOSTICS
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
# ~ 6. 3D VISUALIZATION (PyVista)
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

def plot_vorticity_iso(om, title, fig_name=None, crop_fac=1.0, iso_frac=0.3, y_excl_plus=0.0, view=VIEW_DEFAULT,
                       center=None, half=None):
    """
    Live PyVista view of the |omega| iso-surface colored by omega_z, cropped around the hairpin
    (coordinates in wall units). The last view (as left in the window) is exported as PNG + PDF.

    Args:
        om       : staggered vorticity (omega_x on u-grid, omega_y on v-grid, omega_z on w-grid)
        title    : text shown on the top-left
        fig_name : path without extension for the PNG/PDF export (None: no export)
        crop_fac : margin around the hairpin bounding box, in units of spread
        iso_frac : initial iso-level, as a fraction of max|omega| (slider)
        y_excl_plus : exclude the region within y+ < y_excl_plus of both walls from the crop
                      (the wall shear layers otherwise set max|omega|)
        view     : camera view (see VIEW_DEFAULT)
        center, half : crop box center and half-widths (hx, hy, hz); None -> hairpin bounding box + margin
    """
    import pyvista as pv

        # --> Vorticity at the cell centers, cropped around the hairpin (clipped at the walls)
    om_c = get_cell_centered(*om)
    if half is None:
        center, half = get_hairpin_extent(margin=crop_fac * spread)
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
