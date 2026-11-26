# %% [markdown]
'''
> __Project__    : Injecting artificial vortices into channel flow
> __Author__     : Sugan Velazhagan
> __Date__       : 06/11/2026
> __Description__: Catalog of analytical vortices (see doc/Analytical_vortices.md), defined as python
>                  functions, with an interactive viewer (sliders) to inspect them near the wall.
>
>                  Every vortex is defined through its vector potential A, so that u' = curl(A) is
>                  divergence-free by construction. The image vortex (across the near wall), the
>                  localization window and the wall damping are all applied on A (not on u'), hence
>                  they never break incompressibility.
>
>                  Units are outer units (h = 1, U_b = 1). The base level of strength is u_tau:
>                  strength = peak perturbation velocity / u_tau.
'''
# %%
# ~  LOAD THE FRAMEWORK
#~ ~~~~~~~~~~~~~~~~~~~~~
from framework import *
from scipy.special import exp1
ts0 = tt.time()

# %%
# ~ CONSTANTS
# ~ ~~~~~~~~~~~~~~~~~~~~~
    # --> Peak of (1 - exp(-x^2)) / x, i.e. u_theta_max = C_LO * Gamma / (2 pi r_c) for the Lamb-Oseen profile
_x_LO = np.linspace(0.5, 2.0, 200001)
C_LO  = np.max((1.0 - np.exp(-_x_LO**2)) / _x_LO)       #^ ~ 0.6382 at r = 1.1209 r_c

    # --> Vortex axis -> global indices (x=0, y=1, z=2) of the local right-handed frame (a, s1, s2)
AXIS_MAP = {'x': (0, 1, 2),   #^ a = x, s1 = y, s2 = z  (quasi-streamwise vortex)
            'y': (1, 2, 0),   #^ a = y, s1 = z, s2 = x  (wall-normal vortex)
            'z': (2, 0, 1)}   #^ a = z, s1 = x, s2 = y  (spanwise vortex, e.g. hairpin head)

# %%
# ~ HELPER FUNCTIONS
# ~ ~~~~~~~~~~~~~~~~~~~~~
def get_lamb_oseen_psi(rr2, Gamma, r_c):
    """
    Stream function of the Lamb-Oseen swirl, psi(r), such that u_theta = -dpsi/dr.

    Args:
        rr2  : squared distance from the vortex axis
        Gamma: circulation
        r_c  : core radius

    Returns:
        psi (same shape as rr2)
    """
    xi = np.maximum(rr2 / r_c**2, 1e-14)  #^ ln(xi) + E1(xi) -> -gamma_Euler as xi -> 0
    return -Gamma / (4.0 * np.pi) * (np.log(xi) + exp1(xi))

def get_lamb_oseen_q(rr2, Gamma, r_c):
    """
    Swirl factor q(r) = u_theta / r = Gamma / (2 pi r^2) * (1 - exp(-r^2/r_c^2)); regular at r = 0.
    """
    rr2_safe = np.maximum(rr2, 1e-30)
    q = Gamma / (2.0 * np.pi * rr2_safe) * (-np.expm1(-rr2_safe / r_c**2))
    return np.where(rr2 > 1e-30, q, Gamma / (2.0 * np.pi * r_c**2))

# %%
# ~ ANALYTICAL VORTICES (local frame: a = axis, s1, s2 = transverse)
# ~ ~~~~~~~~~~~~~~~~~~~~~
#@ Each function returns the vector potential (A_a, A_1, A_2), or the analytical velocity
#@ (u_a, u_1, u_2) if flag_vel = True. Coordinates are relative to the vortex center.

def vortex_lamb_oseen(xa, s1, s2, prm, flag_vel=False):
    """
    1. Lamb-Oseen vortex tube (axis a, uniform along a)
        u_theta = Gamma / (2 pi r) [1 - exp(-r^2/r_c^2)],   r_c = spread
        Gamma is set such that max(u_theta) = U_ref.
    """
    r_c   = prm['spread']
    Gamma = 2.0 * np.pi * r_c * prm['U_ref'] / C_LO
    rr2   = s1**2 + s2**2
    zero  = np.zeros_like(rr2 + xa)
    if flag_vel:
        q = get_lamb_oseen_q(rr2, Gamma, r_c)
        return zero, -q * s2 + zero, q * s1 + zero
    return get_lamb_oseen_psi(rr2, Gamma, r_c) + zero, zero, zero

def vortex_burgers(xa, s1, s2, prm, flag_vel=False):
    """
    2. Burgers vortex: Lamb-Oseen swirl with the core set by the strain/diffusion balance
        delta_B = sqrt(4 nu / a) = spread   ->   a = 4 nu / spread^2
       Optional axial strain (flag_strain): u_a = a x_a, u_s = -a s / 2 (windowed when injected).
    """
    A_a, A_1, A_2 = vortex_lamb_oseen(xa, s1, s2, prm, flag_vel)
    if prm['flag_strain']:
        a_strain = 4.0 * prm['nu'] / prm['spread']**2
        if flag_vel:
            return A_a + a_strain * xa, A_1 - 0.5 * a_strain * s1, A_2 - 0.5 * a_strain * s2
        #^ A_strain = (a x_a / 2) (0, -s2, s1)
        return A_a, A_1 - 0.5 * a_strain * xa * s2, A_2 + 0.5 * a_strain * xa * s1
    return A_a, A_1, A_2

def vortex_burgers_packet(xa, s1, s2, prm, flag_vel=False):
    """
    3. Burgers vortex packet (axially modulated Burgers vortex, t = 0)
        omega_a = A(x_a) exp(-r^2/delta_B^2),   A(x_a) = A_0 exp(-x_a^2 / (2 L_0^2)),   L_0 = length
       The radial vorticity omega_r = -d(u_theta)/dx_a appears automatically through curl(A).
    """
    g = np.exp(-xa**2 / (2.0 * prm['length']**2))
    u_a, u_1, u_2 = vortex_lamb_oseen(xa, s1, s2, prm, flag_vel)
    return g * u_a, g * u_1, g * u_2

def vortex_hill(xa, s1, s2, prm, flag_vel=False):
    """
    4. Hill's spherical vortex (lab frame, fluid at rest at infinity), radius a = spread
        Stokes stream function Psi = s^2 Phi(R):
            Phi = U (5/4 - 3 R^2 / (4 a^2))   (R < a)
            Phi = (U / 2) (a / R)^3            (R > a)
       It translates along +a with speed U. U is set such that max|u| = 5U/2 = U_ref.
    """
    a_H = prm['spread']
    U   = prm['U_ref'] / 2.5
    ss2 = s1**2 + s2**2
    RR2 = np.maximum(xa**2 + ss2, 1e-30)
    inside = RR2 < a_H**2
    if flag_vel:
        u_a_in  = U * (2.5 - 3.0 * ss2 / a_H**2 - 1.5 * xa**2 / a_H**2)
        u_s_in  = 1.5 * U * xa / a_H**2                              #^ u_s / s
        u_a_out = 0.5 * U * a_H**3 * (2.0 / RR2**1.5 - 3.0 * ss2 / RR2**2.5)
        u_s_out = 1.5 * U * a_H**3 * xa / RR2**2.5                   #^ u_s / s
        u_a = np.where(inside, u_a_in, u_a_out)
        u_s = np.where(inside, u_s_in, u_s_out)
        return u_a, u_s * s1, u_s * s2
    Phi = np.where(inside, U * (1.25 - 0.75 * RR2 / a_H**2), 0.5 * U * (a_H**2 / RR2)**1.5)
    #^ A = (Psi / s) e_phi = Phi (0, -s2, s1)
    return np.zeros_like(Phi), -Phi * s2, Phi * s1

def vortex_taylor_green(xa, s1, s2, prm, flag_vel=False):
    """
    5. Taylor-Green vortex array (axis a), cell size = spread, k = pi / spread
        u_1 = -U cos(k s1) sin(k s2) cos(k x_a)
        u_2 =  U sin(k s1) cos(k s2) cos(k x_a),   u_a = 0
       Centered such that a vortex core sits at the origin. Periodic -> localized by the window.
    """
    k = np.pi / prm['spread']
    U = prm['U_ref']
    if flag_vel:
        u_1 = -U * np.cos(k * s1) * np.sin(k * s2) * np.cos(k * xa)
        u_2 =  U * np.sin(k * s1) * np.cos(k * s2) * np.cos(k * xa)
        return np.zeros_like(u_1), u_1, u_2
    A_a = U / k * np.cos(k * s1) * np.cos(k * s2) * np.cos(k * xa)
    return A_a, np.zeros_like(A_a), np.zeros_like(A_a)

def vortex_abc(xa, s1, s2, prm, flag_vel=False):
    """
    6. Arnold-Beltrami-Childress flow (A = B = C = U), k = pi / spread
        u_a = U (sin k s2 + cos k s1),  u_1 = U (sin k x_a + cos k s2),  u_2 = U (sin k s1 + cos k x_a)
       Beltrami: curl(u) = k u  ->  vector potential A = u / k. max|u_i| = 2U = U_ref.
       Periodic -> localized by the window.
    """
    k = np.pi / prm['spread']
    U = 0.5 * prm['U_ref']
    u_a = U * (np.sin(k * s2) + np.cos(k * s1))
    u_1 = U * (np.sin(k * xa) + np.cos(k * s2))
    u_2 = U * (np.sin(k * s1) + np.cos(k * xa))
    if flag_vel:
        return u_a, u_1, u_2
    return u_a / k, u_1 / k, u_2 / k

    # --> Catalog of the vortex choices
    #^ flag_axial_window: the vortex is not (well) localized along its axis -> window also applied along the axis
    #^ flag_image       : add the image vortex across the near wall (off for Taylor-Green: a periodic lattice whose
    #^                    image can coincide with itself and cancel; it has its own no-penetration planes)
VORTEX_CHOICES = {
    'lamb_oseen'    : {'func': vortex_lamb_oseen,     'label': 'Lamb-Oseen tube',        'flag_axial_window': False, 'flag_image': True},
    'burgers'       : {'func': vortex_burgers,        'label': 'Burgers vortex',         'flag_axial_window': False, 'flag_image': True},
    'burgers_packet': {'func': vortex_burgers_packet, 'label': 'Burgers vortex packet',  'flag_axial_window': False, 'flag_image': True},
    'hill'          : {'func': vortex_hill,           'label': "Hill's spherical vortex",'flag_axial_window': True, 'flag_image': True},
    'taylor_green'  : {'func': vortex_taylor_green,   'label': 'Taylor-Green array',     'flag_axial_window': True, 'flag_image': False},
    'abc'           : {'func': vortex_abc,            'label': 'ABC flow',               'flag_axial_window': True, 'flag_image': True},
}

# %%
# ~ PARAMETERS, WINDOW AND WALL DAMPING
# ~ ~~~~~~~~~~~~~~~~~~~~~
def get_vortex_parameters(vortex_choice, strength, sense, spread, length, y_core, wall, axis,
                          x_core, z_core, win_factor, wall_damp, flag_strain):
    """
    Collect all vortex parameters in a dictionary (outer units).
    x_core / z_core = None -> middle of the domain in x / z.
    """
    if vortex_choice not in VORTEX_CHOICES:
        raise ValueError(f"Unknown vortex_choice '{vortex_choice}'. Options: {list(VORTEX_CHOICES)}")
    if wall not in ('bottom', 'top'):
        raise ValueError(f"Unknown wall '{wall}'. Options: bottom, top")
    if axis not in AXIS_MAP:
        raise ValueError(f"Unknown axis '{axis}'. Options: {list(AXIS_MAP)}")

    x_c = 0.5 * (xx[0] + xx[-1]) if x_core is None else x_core
    z_c = 0.5 * (zz[0] + zz[-1]) if z_core is None else z_core
    y_c = y_wall_bot + y_core if wall == 'bottom' else y_wall_top - y_core

    return dict(vortex_choice = vortex_choice,
                U_ref         = sense * strength * u_tau,   #^ Base level of strength: u_tau
                spread        = spread,
                length        = length,
                nu            = nu,
                center        = (x_c, y_c, z_c),
                axis          = axis,
                win_factor    = win_factor,
                wall_damp     = wall_damp,
                flag_strain   = flag_strain,
                y_wall_bot    = y_wall_bot,
                y_wall_top    = y_wall_top,
                y_wall_near   = y_wall_bot if wall == 'bottom' else y_wall_top)

def get_window(d_loc, prm):
    """
    Smooth window exp(-(r / W)^4), W = win_factor * spread, r = distance from the vortex axis.
    The axial direction is windowed only for the choices not localized along their axis.
    """
    W   = prm['win_factor'] * prm['spread']
    arg = ((d_loc[1]**2 + d_loc[2]**2) / W**2)**2
    if VORTEX_CHOICES[prm['vortex_choice']]['flag_axial_window'] or \
       (prm['vortex_choice'] == 'burgers' and prm['flag_strain']):   #^ The axial strain grows linearly along the axis
        arg = arg + (d_loc[0] / W)**4
    return np.exp(-arg)

def get_wall_damping(Y, prm):
    """
    f(y) = [1 - exp(-(y_b/delta_w)^2)] [1 - exp(-(y_t/delta_w)^2)],  y_b, y_t = distance to the walls.
    f = df/dy = 0 at both walls -> A = 0 and dA/dy = 0 -> u' = v' = w' = 0 (no-slip, no-penetration).
    """
    d_w = prm['wall_damp']
    return (-np.expm1(-((Y - prm['y_wall_bot']) / d_w)**2)) * (-np.expm1(-((prm['y_wall_top'] - Y) / d_w)**2))

def get_vector_potential(X, Y, Z, prm, flag_window=True, flag_wall=True, flag_vel=False):
    """
    Evaluate the vector potential (Ax, Ay, Az) of the chosen vortex in global coordinates.

    Args:
        X, Y, Z    : global coordinates (any broadcastable shape)
        prm        : vortex parameters (get_vortex_parameters)
        flag_window: apply the localization window
        flag_wall  : add the image vortex across the near wall and apply the wall damping
        flag_vel   : return the raw analytical velocity instead (no image, window, damping)

    Returns:
        (Ax, Ay, Az) or (u, v, w)
    """
    m = AXIS_MAP[prm['axis']]
    def get_field_global(Yq):
        d     = (X - prm['center'][0], Yq - prm['center'][1], Z - prm['center'][2])
        d_loc = (d[m[0]], d[m[1]], d[m[2]])
        F_loc = VORTEX_CHOICES[prm['vortex_choice']]['func'](*d_loc, prm, flag_vel=flag_vel)
        F = [None, None, None]
        for i in range(3):
            F[m[i]] = F_loc[i] + np.zeros_like(X + Yq + Z, dtype=float)
        return F, d_loc

    F, d_loc = get_field_global(Y)
    if flag_vel:
        return tuple(F)

    if flag_wall and VORTEX_CHOICES[prm['vortex_choice']]['flag_image']:
            # --> Image vortex across the near wall: A_img(x) = -P A(P x), P = mirror in y (A is a pseudo-vector)
            #^ Tangential A (Ax, Az) vanishes on the wall -> v' = 0 there, and the far field decays as a dipole
        F_img, _ = get_field_global(2.0 * prm['y_wall_near'] - Y)
        F = [F[0] - F_img[0], F[1] + F_img[1], F[2] - F_img[2]]

    fac = np.ones_like(F[0])
    if flag_window:
        fac = fac * get_window(d_loc, prm)
    if flag_wall:
        fac = fac * get_wall_damping(Y, prm)
    return tuple(fac * F[i] for i in range(3))

# %%
# ~ CURL / DIVERGENCE OPERATORS
# ~ ~~~~~~~~~~~~~~~~~~~~~
def get_curl_collocated(Ax, Ay, Az, x, y, z):
    """
    Curl of a collocated field stored as [z, x, y] (np.gradient, 2nd order). Used for plotting.
    """
    dAx_dz, _, dAx_dy = np.gradient(Ax, z, x, y)
    dAy_dz, dAy_dx, _ = np.gradient(Ay, z, x, y)
    _, dAz_dx, dAz_dy = np.gradient(Az, z, x, y)
    return dAz_dy - dAy_dz, dAx_dz - dAz_dx, dAy_dx - dAx_dy

def get_curl_staggered(Ax, Ay, Az, x_fc, y_fc, z_fc):
    """
    Discrete curl on the staggered (MAC) grid, all arrays stored as [z, x, y] (with halos).
        Ax at (x_cc, y_fc, z_fc), Ay at (x_fc, y_cc, z_fc), Az at (x_fc, y_fc, z_cc)
        -> u at (x_fc, y_cc, z_cc), v at (x_cc, y_fc, z_cc), w at (x_cc, y_cc, z_fc)
    The discrete divergence of the result is zero to round-off (difference operators commute).
    """
    dx = np.diff(x_fc)[None, :, None]
    dy = np.diff(y_fc)[None, None, :]
    dz = np.diff(z_fc)[:, None, None]
    u = np.diff(Az, axis=2) / dy - np.diff(Ay, axis=0) / dz
    v = np.diff(Ax, axis=0) / dz - np.diff(Az, axis=1) / dx
    w = np.diff(Ay, axis=1) / dx - np.diff(Ax, axis=2) / dy
    return u, v, w

def get_divergence_staggered(u, v, w, x_fc, y_fc, z_fc):
    """
    Discrete divergence at the cell centers (with halos) of a staggered field [z, x, y].
    """
    dx = np.diff(x_fc)[None, :, None]
    dy = np.diff(y_fc)[None, None, :]
    dz = np.diff(z_fc)[:, None, None]
    return np.diff(u, axis=1) / dx + np.diff(v, axis=2) / dy + np.diff(w, axis=0) / dz

def get_syn_vortex_staggered(prm, x_fc, y_fc, z_fc, x_cc, y_cc, z_cc):
    """
    Synthetic vortex velocity (u', v', w') on the staggered grid of the solver.
    Each component of A is evaluated on its own cell edge, then u' = curl(A) discretely.
    """
    Z, X, Y = np.meshgrid(z_fc, x_cc, y_fc, indexing='ij')
    Ax = get_vector_potential(X, Y, Z, prm)[0]
    Z, X, Y = np.meshgrid(z_fc, x_fc, y_cc, indexing='ij')
    Ay = get_vector_potential(X, Y, Z, prm)[1]
    Z, X, Y = np.meshgrid(z_cc, x_fc, y_fc, indexing='ij')
    Az = get_vector_potential(X, Y, Z, prm)[2]
    return get_curl_staggered(Ax, Ay, Az, x_fc, y_fc, z_fc)

# %%
# ~ VERIFICATION: curl(A) vs analytical velocity
# ~ ~~~~~~~~~~~~~~~~~~~~~
def check_vortex_choices(prm_base, n_pts=41):
    """
    For every vortex choice, compare curl(A) (no window, no damping) with the analytical velocity
    in a box of size 3 * spread around the core. Logs the relative L2 error.
    """
    x_c, y_c, z_c = prm_base['center']
    h  = 1.5 * prm_base['spread']
    x1 = x_c + np.linspace(-h, h, n_pts)
    y1 = y_c + np.linspace(-h, h, n_pts)
    z1 = z_c + np.linspace(-h, h, n_pts)
    Z, X, Y = np.meshgrid(z1, x1, y1, indexing='ij')
    sl = (slice(2, -2),) * 3
    for choice in VORTEX_CHOICES:
        prm = dict(prm_base, vortex_choice=choice)
        A   = get_vector_potential(X, Y, Z, prm, flag_window=False, flag_wall=False)
        u_n = get_curl_collocated(*A, x1, y1, z1)
        u_a = get_vector_potential(X, Y, Z, prm, flag_vel=True)
        err = np.sqrt(sum(np.sum((u_n[i][sl] - u_a[i][sl])**2) for i in range(3)) /
                      sum(np.sum(u_a[i][sl]**2) for i in range(3)))
        log.info(f"  {choice:<15s}: rel. L2 error |curl(A) - u_analytical| = {err:.2e}")

# %%
if __name__ == '__main__':
    # ~ INTERACTIVE VIEWER
    # ~ ~~~~~~~~~~~~~~~~~~~~~
    from matplotlib.widgets import Slider, RadioButtons

    prm = get_vortex_parameters(vortex_choice, strength, sense, spread, length, y_core, 'bottom', axis,
                                x_core, z_core, win_factor, wall_damp, flag_strain)
    log.info(f"Verification of the vortex choices (axis = {axis}):")
    check_vortex_choices(prm)

        # --> Plotting planes (bottom wall, through the core)
    x_c, y_c, z_c = prm['center']
    h_slab = 1e-3
    y_pl   = y_wall_bot + np.linspace(0.0, 1.0, 161)
    z_pl   = z_c + np.linspace(-0.6, 0.6, 161)
    x_pl   = x_c + np.linspace(-1.5, 1.5, 241)
    x_slab = x_c + h_slab * np.arange(-2.0, 3.0)
    z_slab = z_c + h_slab * np.arange(-2.0, 3.0)

    def get_planes(prm):
        """
        Velocity u' = curl(A) and vorticity w' = curl(u') (with image, window and wall damping)
        on the y-z plane at x_c and on the x-y plane at z_c (5-point slabs across each plane).
        """
        y_c = prm['center'][1]
            # --> y-z plane (thin slab in x)
        Z, X, Y = np.meshgrid(z_pl, x_slab, y_pl, indexing='ij')
        uvw = get_curl_collocated(*get_vector_potential(X, Y, Z, prm), x_slab, y_pl, z_pl)
        omg = get_curl_collocated(*uvw, x_slab, y_pl, z_pl)
        u_yz, v_yz, w_yz = (c[:, 2, :] for c in uvw)                               #^ [z, y]
        omx_yz = omg[0][:, 2, :]
            # --> x-y plane (thin slab in z)
        Z, X, Y = np.meshgrid(z_slab, x_pl, y_pl, indexing='ij')
        uvw = get_curl_collocated(*get_vector_potential(X, Y, Z, prm), x_pl, y_pl, z_slab)
        omg = get_curl_collocated(*uvw, x_pl, y_pl, z_slab)
        u_xy, v_xy, w_xy = (c[2] for c in uvw)                                    #^ [x, y]
        omm_xy = np.sqrt(sum(c[2]**2 for c in omg))
        return dict(u_yz=u_yz, v_yz=v_yz, w_yz=w_yz, omx_yz=omx_yz,
                    u_xy=u_xy, v_xy=v_xy, w_xy=w_xy, omm_xy=omm_xy,
                    j_core=np.argmin(np.abs(y_pl - y_c)))

    def plot_quiver(ax, x, y, U, V, sk=8):
        """
        Quiver with a fixed scale: the largest arrow is ~ 1/12 of the axis width.
        """
        Us, Vs = U[::sk, ::sk].T, V[::sk, ::sk].T
        mag = np.max(np.hypot(Us, Vs))
        if mag > 0:
            ax.quiver(x[::sk], y[::sk], Us, Vs, color='k', scale=12 * mag, width=0.003)

        # --> Figure layout
    fig = plt.figure(figsize=(15, 8.5))
    gs  = GridSpec(2, 3, figure=fig, left=0.05, right=0.97, top=0.93, bottom=0.36, wspace=0.35, hspace=0.45)
    ax_yz = fig.add_subplot(gs[:, 0])
    ax_xy = fig.add_subplot(gs[0, 1:])
    ax_pz = fig.add_subplot(gs[1, 1])
    ax_py = fig.add_subplot(gs[1, 2])
    cax_yz = make_axes_locatable(ax_yz).append_axes('right', size='4%', pad=0.05)
    cax_xy = make_axes_locatable(ax_xy).append_axes('right', size='2%', pad=0.05)
    info   = fig.text(0.70, 0.04, '', fontsize=9, family='monospace', va='bottom')

    ax_radio_v = fig.add_axes([0.02, 0.03, 0.15, 0.24])
    ax_radio_a = fig.add_axes([0.18, 0.03, 0.05, 0.24])
    radio_v = RadioButtons(ax_radio_v, list(VORTEX_CHOICES), active=list(VORTEX_CHOICES).index(vortex_choice))
    radio_a = RadioButtons(ax_radio_a, list(AXIS_MAP), active=list(AXIS_MAP).index(axis))
    ax_radio_v.set_title('vortex', fontsize=9)
    ax_radio_a.set_title('axis', fontsize=9)

    slider_specs = [  # name, label, min, max, init
        ('strength',   r'strength [$u_\tau$]', 0.0,  20.0, strength),
        ('spread',     'spread [h]',           0.01, 0.30, spread),
        ('y_core',     'y_core [h]',           0.02, 0.90, y_core),
        ('length',     'length [h]',           0.05, 2.00, length),
        ('win_factor', 'win_factor',           1.5,  10.0, win_factor),
        ('wall_damp',  'wall_damp [h]',        0.005, 0.15, wall_damp),
    ]
    sliders = {}
    for i, (name, label, vmin, vmax, vinit) in enumerate(slider_specs):
        ax_s = fig.add_axes([0.31, 0.25 - i * 0.04, 0.30, 0.025])
        sliders[name] = Slider(ax_s, label, vmin, vmax, valinit=vinit)

    def update(_=None):
        """
        Rebuild the vortex from the widgets and redraw all panels.
        """
        prm = get_vortex_parameters(radio_v.value_selected, sliders['strength'].val, sense,
                                    sliders['spread'].val, sliders['length'].val, sliders['y_core'].val,
                                    'bottom', radio_a.value_selected, x_core, z_core,
                                    sliders['win_factor'].val, sliders['wall_damp'].val, flag_strain)
        pl  = get_planes(prm)
        y_c = prm['center'][1]
        yy_pl = y_pl - y_wall_bot

            # --> y-z plane: omega_x + (w, v)
        ax_yz.clear(); cax_yz.clear()
        om_max = np.max(np.abs(pl['omx_yz'])) + 1e-30
        cf = ax_yz.contourf(z_pl - z_c, yy_pl, pl['omx_yz'].T, levels=np.linspace(-om_max, om_max, 41), cmap='RdBu_r')
        plot_quiver(ax_yz, z_pl - z_c, yy_pl, pl['w_yz'], pl['v_yz'])
        fig.colorbar(cf, cax=cax_yz, label=r"$\omega_x'$")
        ax_yz.set_xlabel(r'$z - z_c$'); ax_yz.set_ylabel(r'$y$ (from wall)')
        ax_yz.set_title(r"y-z plane at $x_c$: $\omega_x'$, $(w', v')$")
        ax_yz.set_aspect('equal')

            # --> x-y plane: |omega| + (u, v)
        ax_xy.clear(); cax_xy.clear()
        om_max = np.max(pl['omm_xy']) + 1e-30
        cf = ax_xy.contourf(x_pl - x_c, yy_pl, pl['omm_xy'].T, levels=np.linspace(0.0, om_max, 41), cmap='viridis')
        plot_quiver(ax_xy, x_pl - x_c, yy_pl, pl['u_xy'], pl['v_xy'])
        fig.colorbar(cf, cax=cax_xy, label=r"$|\omega'|$")
        ax_xy.set_xlabel(r'$x - x_c$'); ax_xy.set_ylabel(r'$y$ (from wall)')
        ax_xy.set_title(r"x-y plane at $z_c$: $|\omega'|$, $(u', v')$")
        ax_xy.set_aspect('equal')

            # --> Profiles through the core (normalized by u_tau)
        j = pl['j_core']
        ax_pz.clear()
        ax_pz.plot(z_pl - z_c, pl['u_yz'][:, j] / u_tau, label="u'")
        ax_pz.plot(z_pl - z_c, pl['v_yz'][:, j] / u_tau, label="v'")
        ax_pz.plot(z_pl - z_c, pl['w_yz'][:, j] / u_tau, label="w'")
        ax_pz.set_xlabel(r'$z - z_c$'); ax_pz.set_ylabel(r"$u_i'/u_\tau$"); ax_pz.set_title('along z at core')
        ax_pz.legend(fontsize=8); ax_pz.grid(alpha=0.3)

        k = np.argmin(np.abs(z_pl - z_c))
        ax_py.clear()
        ax_py.plot(pl['u_yz'][k, :] / u_tau, yy_pl, label="u'")
        ax_py.plot(pl['v_yz'][k, :] / u_tau, yy_pl, label="v'")
        ax_py.plot(pl['w_yz'][k, :] / u_tau, yy_pl, label="w'")
        ax_py.axhline(y_c - y_wall_bot, color='gray', ls='--', lw=0.8)
        ax_py.set_xlabel(r"$u_i'/u_\tau$"); ax_py.set_ylabel(r'$y$ (from wall)'); ax_py.set_title('along y at core')
        ax_py.legend(fontsize=8); ax_py.grid(alpha=0.3)

            # --> Derived quantities
        uvw_max = max(np.max(np.abs(pl[c])) for c in ('u_yz', 'v_yz', 'w_yz', 'u_xy', 'v_xy', 'w_xy'))
        sp = prm['spread']
        info.set_text(
            f"{VORTEX_CHOICES[prm['vortex_choice']]['label']}\n"
            f"U_ref        = {prm['U_ref']:.4f} = {prm['U_ref']/u_tau:.2f} u_tau\n"
            f"max|u'| (planes)= {uvw_max/u_tau:.2f} u_tau\n"
            f"Gamma (LO)   = {2*np.pi*sp*prm['U_ref']/C_LO:.4e}\n"
            f"spread+      = {sp/delta_nu:.1f}\n"
            f"y_core+      = {(y_c-y_wall_bot)/delta_nu:.1f}\n"
            f"length+      = {prm['length']/delta_nu:.1f}\n"
            f"Burgers a    = {4*nu/sp**2:.4f}\n"
            f"window W+    = {prm['win_factor']*sp/delta_nu:.1f}")
        fig.suptitle(f"{VORTEX_CHOICES[prm['vortex_choice']]['label']}  (axis = {prm['axis']}, "
                     f"Re_tau = {Re_tau:.0f}, outer units)")
        fig.canvas.draw_idle()

    for s in sliders.values():
        s.on_changed(update)
    radio_v.on_clicked(update)
    radio_a.on_clicked(update)
    update()

    log.info(f"Total time taken is {tt.time()-ts0:.2f} seconds")
    plt.show()

# %%
