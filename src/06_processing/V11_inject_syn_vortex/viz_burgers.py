# %% [markdown]
'''
> __Project__    : Injecting artificial vortices into channel flow
> __Author__     : Sugan Velazhagan
> __Date__       : 25/11/2026
> __Description__: Verify the Burgers vortex packet on the staggered simulation grid:
>                  analytical vorticity vs curl of the injected velocity (relative L2 error).
'''
# %%
# ~  LOAD THE FRAMEWORK
# ~ ~~~~~~~~~~~~~~~~~~~~~
from burgers_functions import *
ts = tt.time()

# %%
#* ==========================
#* INTERIOR SLICES
#* ==========================
    # --> Leave 1 ghost/halo cell on either side in all three directions
sl_in = (slice(1, -1), slice(1, -1), slice(1, -1))

# %%
#* ==========================
#* VORTEX PARAMETERS
#* ==========================
y_core_viz = 1.0                                     #^ Height of the core from the chosen wall
om_choice  = 'om_a'                                  #^ Vorticity to visualize: 'om_a' (analytical) or 'om_n' (curl of u')
crop_fac   = 3.0                                     #^ Crop half-width: crop_fac * spread (transverse), crop_fac * length (axial)

    # --> Camera view (fixed and reproducible; applied in this order: preset -> view_up -> azimuth/elevation/roll -> zoom)
view_preset    = 'xy'         #^ Starting view: 'xy' (look down z, x-y plane), 'xz' (look down y), 'yz' (look down x), 'iso' (oblique)
view_up        = (0, 1, 0)    #^ Direction pointing up in the image: (0, 1, 0) -> y+ up (wall at the bottom), (0, 0, 1) -> z+ up
view_azimuth   = 30.0         #^ Rotate the camera about the up direction [deg]; > 0 swings it to the right, 0 = keep the preset
view_elevation = 20.0         #^ Tilt the camera up (> 0, look from above) or down (< 0) [deg]
view_roll      = 0.0          #^ Spin the image in its own plane [deg]
view_zoom      = 1.0          #^ > 1 zoom in, < 1 zoom out
view_exact     = None         #^ Exact camera, overrides all of the above: [(eye_x, eye_y, eye_z), (focus_x, focus_y, focus_z), (up_x, up_y, up_z)]
                              #^ (wall units). To copy an interactive view: rotate it, run print(pl.camera_position), paste here.
center_viz = (0.5 * Lx, y_wall_near + (y_core_viz if wall == 'bottom' else -y_core_viz), 0.5 * Lz + z0)
prm_viz    = dict(center=center_viz, Gamma=Gamma, spread=spread, length=length, axis=axis)

log.info(f"Center = ({center_viz[0]:.4f}, {center_viz[1]:.4f}, {center_viz[2]:.4f}), axis = {axis}")

# %%
# ~ 1. ZERO FIELDS (staggered, [z, x, y] with halos, as read by read_ucvcwc_file)
# ~ ~~~~~~~~~~~~~~~~~~~~~
u   = np.zeros((ngz + 1, ngx + 2, ngy + 1))
v   = np.zeros((ngz + 1, ngx + 1, ngy + 2))
w   = np.zeros((ngz + 2, ngx + 1, ngy + 1))
omx = np.zeros_like(u)          #^ omega_x on the u-grid
omy = np.zeros_like(v)          #^ omega_y on the v-grid
omz = np.zeros_like(w)          #^ omega_z on the w-grid

# %%
# ~ 2. VORTICITY: ANALYTICAL
# ~ ~~~~~~~~~~~~~~~~~~~~~
    #^ Image vortex as in the injected velocity (flag_image in settings.yaml)
om_a = get_burgers_vorticity(**prm_viz, flag_image=flag_image, y_wall_near=y_wall_near)
omx += om_a[0]; omy += om_a[1]; omz += om_a[2]
om_a = (omx, omy, omz)

# %%
# ~ 3. VELOCITY -> CURL -> VORTICITY
# ~ ~~~~~~~~~~~~~~~~~~~~~
    #^ Window and wall damping switched off: same field as the analytical vorticity
du, dv, dw = get_burgers_velocity(**prm_viz, win_factor=win_factor, wall_damp=wall_damp, y_wall_near=y_wall_near,
                                  flag_localize=False, flag_image=flag_image)
u += du; v += dv; w += dw
del du, dv, dw
om_n = get_vorticity(u, v, w)

# %%
# ~ COMPARISON (relative L2 error, interior points)
# ~ ~~~~~~~~~~~~~~~~~~~~~
    # --> Core region: within 2 spreads of the vortex axis
grids = [(xu, yu, zu), (xv, yv, zv), (xw, yw, zw)]
m     = AXIS_MAP[axis]

log.info('=' * 70)
log.info(f"{'comp.':<8}{'max|omega_a|':>14}{'max|omega_n|':>14}{'L2 err (all)':>14}{'L2 err (core)':>15}")
for n, (om_ref, om_num, comp) in enumerate(zip((omx, omy, omz), om_n, ('omega_x', 'omega_y', 'omega_z'))):
    Z, X, Y = np.meshgrid(grids[n][2], grids[n][0], grids[n][1], indexing='ij', sparse=True)
    d       = (X - center_viz[0], Y - center_viz[1], Z - center_viz[2])
    core    = np.broadcast_to(np.hypot(d[m[1]], d[m[2]]) < 2.0 * spread, om_ref.shape)[sl_in]
    a, b    = om_ref[sl_in], om_num[sl_in]
    norm    = np.linalg.norm(a)
    err_all = np.linalg.norm(b - a) / norm if norm > 0 else np.nan
    norm_c  = np.linalg.norm(a[core])
    err_c   = np.linalg.norm(b[core] - a[core]) / norm_c if norm_c > 0 else np.nan
    log.info(f"{comp:<8}{np.max(np.abs(a)):>14.4e}{np.max(np.abs(b)):>14.4e}{err_all:>14.4f}{err_c:>15.4f}")
log.info('=' * 70)

# %%
# ~ 3D VISUALIZATION (PyVista): |omega| iso-surfaces colored by omega_z
# ~ ~~~~~~~~~~~~~~~~~~~~~
import pyvista as pv

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

def get_rectilinear_grid(F_dict, x, y, z):
    """
    PyVista rectilinear grid (non-uniform x, y, z) with the fields F[z, x, y] attached as point data.
    """
    grid = pv.RectilinearGrid(x, y, z)
    for name, F in F_dict.items():
        grid.point_data[name] = np.transpose(F, (1, 2, 0)).ravel(order='F')   #^ [z, x, y] -> x fastest
    return grid

    # --> Vorticity at the cell centers, cropped around the vortex center
om_c   = get_cell_centered(*(om_a if om_choice == 'om_a' else om_n))
half   = [crop_fac * spread] * 3
half[AXIS_MAP[axis][0]] = crop_fac * length
i_cr   = get_crop_indices(xx_cc, center_viz[0], half[0])
j_cr   = get_crop_indices(yy_cc, center_viz[1], half[1])
j_cr   = j_cr[(yy_cc[j_cr] >= y_wall_bot) & (yy_cc[j_cr] <= y_wall_top)]
k_cr   = get_crop_indices(zz_cc, center_viz[2], half[2])
sub    = np.ix_(k_cr, i_cr, j_cr)
om_mag = np.sqrt(om_c[0][sub]**2 + om_c[1][sub]**2 + om_c[2][sub]**2)

    #^ Coordinates in wall units: x+ = x / delta_nu, y+ = (distance from the bottom wall) / delta_nu, z+ = z / delta_nu
grid = get_rectilinear_grid({'omega_mag': om_mag, 'omega_z': om_c[2][sub]},
                            xx_cc[i_cr] / delta_nu, (yy_cc[j_cr] - y_wall_bot) / delta_nu, zz_cc[k_cr] / delta_nu)
log.info(f"Visualizing {om_choice}: crop {len(i_cr)} x {len(j_cr)} x {len(k_cr)} points, max|omega| = {om_mag.max():.4e}")

    # --> Live view: iso-surface of |omega| (slider: fraction of max|omega|), colored by omega_z
om_max  = float(om_mag.max())
omz_lim = float(np.max(np.abs(om_c[2][sub]))) + 1e-30

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

pl.add_slider_widget(update_iso, rng=(0.05, 0.95), value=0.3, title=r'$\Vert\omega\Vert \,/\, \Vert\omega\Vert_{max}$',  #^ No '|' in VTK math text (cell separator)
                     pointa=(0.6, 0.92), pointb=(0.95, 0.92), style='modern', title_height=0.045)
pl.add_text(rf"{om_choice}: iso-surface of $\Vert\omega\Vert$, colored by $\omega_z$ (axis = {axis})", font_size=15)
pl.show_grid(xtitle='x+', ytitle='y+', ztitle='z+', fmt='%.0f', font_size=27)     #^ Integer ticks in wall units
    # --> Camera (see 'Camera view' settings at the top)
if view_exact is not None:
    pl.camera_position = view_exact
else:
    pl.camera_position = view_preset
    pl.set_viewup(view_up)
    pl.camera.azimuth   = view_azimuth
    pl.camera.elevation = view_elevation
    pl.camera.roll      = pl.camera.roll + view_roll
    pl.camera.zoom(view_zoom)
pl.show(auto_close=False)          #^ Keep the scene after the window is closed, to export it

    # --> Export the last view (as left in the window): high-resolution PNG + vector PDF
create_directory(root_out)
fig_name = root_out + f'burgers_iso_{om_choice}_axis{axis}'
pl.screenshot(fig_name + '.png', scale=2)
pl.save_graphic(fig_name + '.pdf')
pl.close()
log.info(f"Figures saved: {fig_name}.png / .pdf")

log.info(f"Total time taken: {tt.time() - ts:.2f} seconds")
# %%
