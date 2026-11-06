# %% [markdown]
'''
> __Project__    : Injecting artificial vortices into channel flow
> __Author__     : Sugan Velazhagan
> __Date__       : 06/11/2026
> __Description__: Read a turbulent channel flow state (initial condition), add an analytical vortex
>                  (vortex_functions.py) near the bottom/top wall, plot contours, and save the
>                  result as a new initial condition in data/syn_vortex/type_<vortex_choice>/.
>
>                  Usage: python add_syn_vortex.py [vortex_choice]   (overrides settings.yaml)
'''
# %%
# ~  LOAD THE FRAMEWORK
#~ ~~~~~~~~~~~~~~~~~~~~~
from framework import *
from vortex_functions import (VORTEX_CHOICES, get_vortex_parameters, get_syn_vortex_staggered,
                                 get_divergence_staggered)
ts0 = tt.time()

# %%
# ~ VORTEX PARAMETERS
# ~ ~~~~~~~~~~~~~~~~~~~~~
prm = get_vortex_parameters(vortex_choice, strength, sense, spread, length, y_core, wall, axis,
                            x_core, z_core, win_factor, wall_damp, flag_strain)
x_c, y_c, z_c = prm['center']

log.info(f"Vortex           : {VORTEX_CHOICES[vortex_choice]['label']} (axis = {axis}, wall = {wall})")
log.info(f"U_ref            = {prm['U_ref']:.4e}  ({sense * strength:.2f} u_tau)")
log.info(f"spread           = {spread:.4f}  ({spread / delta_nu:.1f} wall units)")
log.info(f"y_core           = {y_core:.4f}  ({y_core / delta_nu:.1f} wall units)")
log.info(f"center (x, y, z) = ({x_c:.4f}, {y_c:.4f}, {z_c:.4f})")

    # --> Resolution check: grid points per core radius in each direction
j_c = np.argmin(np.abs(yy - y_c))
n_core = {'x': spread / (xx[1] - xx[0]),
          'y': spread / (yy[j_c + 1] - yy[j_c]),
          'z': spread / (zz[1] - zz[0])}
for key, val in n_core.items():
    log.info(f"Points per spread along {key}: {val:.2f}")
    if val < 2.0 and key != axis:
        log.warning(f"The vortex is under-resolved along {key} ({val:.2f} points per spread < 2). "
                    f"Consider a larger spread or a finer grid.")

# %%
# ~ READ THE INITIAL CONDITION
# ~ ~~~~~~~~~~~~~~~~~~~~~
from io_operations import read_ucvcwc_file, write_ucvcwc_file

_, u0, v0, w0 = read_ucvcwc_file(field_folder, ngx, ngy, ngz, frame_IC)
log.info(f"Read initial condition (frame {frame_IC}): u {u0.shape}, v {v0.shape}, w {w0.shape}")
log.info(f"Bulk velocity (approx.) = {np.mean(u0[1:-1, 1:-1, 1:-1]):.4f}")

# %%
# ~ BUILD AND ADD THE SYNTHETIC VORTEX
# ~ ~~~~~~~~~~~~~~~~~~~~~
ts = tt.time()
    # --> u' = curl(A) on the staggered grid (discretely divergence-free)
u_p, v_p, w_p = get_syn_vortex_staggered(prm, xx_fc, yy_fc, zz_fc, xx_cc, yy_cc, zz_cc)
assert u_p.shape == u0.shape and v_p.shape == v0.shape and w_p.shape == w0.shape, "Staggered shapes mismatch"
log.info(f"Synthetic vortex built in {tt.time() - ts:.2f} seconds")

u1 = u0 + u_p
v1 = v0 + v_p
w1 = w0 + w_p

    # --> Diagnostics
sl_in   = (slice(1, -1), slice(1, -1), slice(1, -1))        #^ Interior cells (exclude halos)
div_p   = get_divergence_staggered(u_p, v_p, w_p, xx_fc, yy_fc, zz_fc)[sl_in]
div_0   = get_divergence_staggered(u0, v0, w0, xx_fc, yy_fc, zz_fc)[sl_in]
div_1   = get_divergence_staggered(u1, v1, w1, xx_fc, yy_fc, zz_fc)[sl_in]
div_ref = np.abs(prm['U_ref']) / spread
log.info(f"max|div(u')| / (U_ref/spread) = {np.max(np.abs(div_p)) / div_ref:.3e}")
log.info(f"max|div(u)|  before / after   = {np.max(np.abs(div_0)):.3e} / {np.max(np.abs(div_1)):.3e}")

j_w = 0 if wall == 'bottom' else -1                         #^ Wall-adjacent u, w (first interior cell)
log.info(f"max|u'|, |v'|, |w'| / u_tau   = {np.max(np.abs(u_p))/u_tau:.3f}, "
         f"{np.max(np.abs(v_p))/u_tau:.3f}, {np.max(np.abs(w_p))/u_tau:.3f}")
log.info(f"max|v'| on the walls / u_tau  = {max(np.max(np.abs(v_p[:, :, 1])), np.max(np.abs(v_p[:, :, -2])))/u_tau:.3e}")
log.info(f"max|u'|, |w'| first cell / u_tau = {np.max(np.abs(u_p[:, :, 1 if j_w == 0 else -2]))/u_tau:.3e}, "
         f"{np.max(np.abs(w_p[:, :, 1 if j_w == 0 else -2]))/u_tau:.3e}")
log.info(f"max|u'_i| on x/z domain edges / u_tau = "
         f"{max(np.max(np.abs(u_p[:, :2, :])), np.max(np.abs(u_p[:, -2:, :])), np.max(np.abs(w_p[:2])), np.max(np.abs(w_p[-2:])))/u_tau:.3e}")

    # --> Kinetic energy of the perturbation vs the turbulent fluctuations (cell-centered, interior)
dV   = (np.diff(zz_fc)[1:-1, None, None] * np.diff(xx_fc)[None, 1:-1, None] * np.diff(yy_fc)[None, None, 1:-1])
def get_cc(u, v, w):
    """
    Interpolate staggered velocities to the cell centers (with halos) -> arrays [z, x, y] of (ngz+1, ngx+1, ngy+1).
    """
    return 0.5 * (u[:, 1:, :] + u[:, :-1, :]), 0.5 * (v[:, :, 1:] + v[:, :, :-1]), 0.5 * (w[1:] + w[:-1])
uc_p, vc_p, wc_p = get_cc(u_p, v_p, w_p)
E_p = 0.5 * np.sum((uc_p**2 + vc_p**2 + wc_p**2)[sl_in] * dV)
log.info(f"Kinetic energy of the vortex  = {E_p:.4e} (domain volume {np.sum(dV):.3f})")

# %%
# ~ SAVE THE NEW INITIAL CONDITION
# ~ ~~~~~~~~~~~~~~~~~~~~~
if flag_output:
    if os.path.exists(save_loc):
        log.warning(f"Save location exists, files will be overwritten: {save_loc}")
    create_directory(save_loc)

        # --> Copy the simulation package (grid, inputs, archive, ...) except the velocity fields
    vel_files = [pre + "{:07d}".format(frame_IC) for pre in ('01uc00.dble.', '02vc00.dble.', '03wc00.dble.')]
    for fname in sorted(os.listdir(field_folder)):
        if fname in vel_files or not os.path.isfile(field_folder + fname):
            continue
        shutil.copy2(field_folder + fname, save_loc + fname)
        log.debug(f"Copied {fname}")

        # --> Write the velocity fields with the vortex
    write_ucvcwc_file(u1, v1, w1, num=frame_IC, addr=save_loc)

        # --> Verify by reading back
    _, u_chk, v_chk, w_chk = read_ucvcwc_file(save_loc, ngx, ngy, ngz, frame_IC)
    assert np.array_equal(u_chk, u1) and np.array_equal(v_chk, v1) and np.array_equal(w_chk, w1), "Read-back mismatch"
    log.info("Velocity fields written and verified (read-back identical)")

        # --> Save the vortex parameters for reproducibility
    vortex_info = {key: settings[key] for key in ('sim_name', 'frame_IC', 'vortex_choice', 'strength', 'sense',
                   'spread', 'length', 'y_core', 'wall', 'axis', 'win_factor', 'wall_damp', 'flag_strain')}
    vortex_info.update(dict(vortex_choice=vortex_choice,
                            x_core=float(x_c), y_center=float(y_c), z_core=float(z_c),
                            U_ref=float(prm['U_ref']), u_tau=float(u_tau), nu=float(nu),
                            Re_tau=float(Re_tau), Re_bulk=float(Re_bulk),
                            max_div_perturbation=float(np.max(np.abs(div_p))),
                            kinetic_energy=float(E_p)))
    with open(save_loc + 'vortex_parameters.yaml', 'w') as file:
        yaml.safe_dump(vortex_info, file, sort_keys=False)
    log.info(f"Saved new initial condition to: {save_loc}")

# %%
# ~ CONTOUR PLOTS
# ~ ~~~~~~~~~~~~~~~~~~~~~
if flag_plot:
    create_directory(fig_loc)

    def get_vorticity_cc(uc, vc, wc):
        """
        Vorticity (omega_x, omega_y, omega_z) of a cell-centered field [z, x, y] (np.gradient, 2nd order).
        """
        omx = np.gradient(wc, yy_cc, axis=2) - np.gradient(vc, zz_cc, axis=0)
        omy = np.gradient(uc, zz_cc, axis=0) - np.gradient(wc, xx_cc, axis=1)
        omz = np.gradient(vc, xx_cc, axis=1) - np.gradient(uc, yy_cc, axis=2)
        return omx, omy, omz

        # --> Total (before / after) and perturbation vorticity at the cell centers
    om_0 = get_vorticity_cc(*get_cc(u0, v0, w0))
    om_1 = get_vorticity_cc(*get_cc(u1, v1, w1))
    om_p = get_vorticity_cc(uc_p, vc_p, wc_p)
    uvw_p = (uc_p, vc_p, wc_p)

        # --> Plane indices through the core (cell centers) and distance from the chosen wall
    i_c   = np.argmin(np.abs(xx_cc - x_c))
    j_c   = np.argmin(np.abs(yy_cc - y_c))
    k_c   = np.argmin(np.abs(zz_cc - z_c))
    y_dw  = (yy_cc - y_wall_bot) if wall == 'bottom' else (y_wall_top - yy_cc)
    j_pl  = np.where((y_dw >= 0.0) & (y_dw <= 1.0))[0]       #^ Half channel next to the chosen wall
    i_win = np.where(np.abs(xx_cc - x_c) <= 1.5)[0]
    k_win = np.where(np.abs(zz_cc - z_c) <= 0.75)[0]
    i_win = i_win if len(i_win) > 2 else np.arange(len(xx_cc))
    k_win = k_win if len(k_win) > 2 else np.arange(len(zz_cc))

    def plot_plane(ax, X, Y, F, title, xlabel, ylabel, F_lim=None, flag_sym=True, cmap=None):
        """
        Filled contour of F[x, y] on the axis ax.
        Color limits: F_lim if given, else the 99.5th percentile of |F| (so that thin wall layers
        do not saturate the scale); values beyond the limits are shown with extended colors.
        """
        F_lim = F_lim or (np.percentile(np.abs(F), 99.5) + 1e-30)
        if flag_sym:
            levels, cmap = np.linspace(-F_lim, F_lim, 41), cmap or 'RdBu_r'
        else:
            levels, cmap = np.linspace(0.0, F_lim, 41), cmap or 'viridis'
        cf = ax.contourf(X, Y, F.T, levels=levels, cmap=cmap, extend='both' if flag_sym else 'max')
        cax = make_axes_locatable(ax).append_axes('right', size='3%', pad=0.05)
        plt.colorbar(cf, cax=cax)
        ax.set_title(title, fontsize=11)
        ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
        ax.set_aspect('equal')
        return cf

    def plot_quiver(ax, X, Y, U, V, sk=1):
        """
        Quiver with a fixed scale: the largest arrow is ~ 1/12 of the axis width.
        """
        Us, Vs = U[::sk, ::sk].T, V[::sk, ::sk].T
        mag = np.max(np.hypot(Us, Vs))
        if mag > 0:
            ax.quiver(X[::sk], Y[::sk], Us, Vs, color='k', scale=12 * mag, width=0.003)

    def plot_plane_set(axs, X, Y, get_slice, n_om, n_q, labels, file_name, title):
        """
        Four panels for one plane: total vorticity normal to the plane (before / after, same color
        scale), perturbation vorticity with in-plane perturbation velocity, in-plane |u'| / u_tau.

        Args:
            get_slice: function F[z, x, y] -> F[X, Y] on the plane
            n_om     : index of the vorticity component normal to the plane (0: x, 1: y, 2: z)
            n_q      : indices (horizontal, vertical) of the in-plane velocity components
            labels   : (component name, horizontal label, vertical label, in-plane names)
        """
        om_b, om_a, om_d = get_slice(om_0[n_om]), get_slice(om_1[n_om]), get_slice(om_p[n_om])
        q_h, q_v = get_slice(uvw_p[n_q[0]]), get_slice(uvw_p[n_q[1]])
        lim = max(np.percentile(np.abs(om_b), 99.5), np.percentile(np.abs(om_a), 99.5))
        c, xl, yl, q = labels
        plot_plane(axs[0], X, Y, om_b, rf"$\omega_{c}$ (before)", xl, yl, F_lim=lim)
        plot_plane(axs[1], X, Y, om_a, rf"$\omega_{c}$ (after)",  xl, yl, F_lim=lim)
        plot_plane(axs[2], X, Y, om_d, rf"$\omega_{c}'$ and $({q})'$", xl, yl)
        plot_quiver(axs[2], X, Y, q_h, q_v)
        plot_plane(axs[3], X, Y, np.hypot(q_h, q_v) / u_tau, rf"$|({q})'|/u_\tau$", xl, yl, flag_sym=False)
        fig = axs[0].figure
        fig.suptitle(f"{title}\n{tag}", fontsize=11)
        fig.tight_layout()
        fig.savefig(fig_loc + file_name, dpi=200, bbox_inches='tight')

    y_lab = r'$y$ (from wall)'
    tag   = f"{VORTEX_CHOICES[vortex_choice]['label']} (axis = {axis}, {wall} wall, " \
            f"strength = {sense * strength:.1f} $u_\\tau$, spread = {spread}, $y_{{core}}$ = {y_core})"
    zk, yj, xi = zz_cc[k_win] - z_c, y_dw[j_pl], xx_cc[i_win] - x_c

        # --> y-z plane at x_c: omega_x, (w', v')
    _, axs = plt.subplots(2, 2, figsize=(12, 7.5))
    plot_plane_set(axs.ravel(), zk, yj, lambda F: F[:, i_c, :][np.ix_(k_win, j_pl)], 0, (2, 1),
                   ('x', r'$z - z_c$', y_lab, 'w, v'), 'contour_yz_plane.png', f"y-z plane at $x_c$ = {xx_cc[i_c]:.3f}")

        # --> x-y plane at z_c: omega_z, (u', v')
    _, axs = plt.subplots(4, 1, figsize=(12, 10))
    plot_plane_set(axs, xi, yj, lambda F: F[k_c][np.ix_(i_win, j_pl)], 2, (0, 1),
                   ('z', r'$x - x_c$', y_lab, 'u, v'), 'contour_xy_plane.png', f"x-y plane at $z_c$ = {zz_cc[k_c]:.3f}")

        # --> x-z plane at y_core: omega_y, (u', w')
    _, axs = plt.subplots(4, 1, figsize=(12, 9))
    plot_plane_set(axs, xi, zk, lambda F: F[:, :, j_c][np.ix_(k_win, i_win)].T, 1, (0, 2),
                   ('y', r'$x - x_c$', r'$z - z_c$', 'u, w'), 'contour_xz_plane.png',
                   f"x-z plane at $y$ = {yy_cc[j_c]:.3f} (core height)")

    log.info(f"Contour plots saved to: {fig_loc}")
    if mpl.get_backend().lower() != 'agg':
        plt.show()

log.info(f"Total time taken is {tt.time() - ts0:.2f} seconds")
# %%
