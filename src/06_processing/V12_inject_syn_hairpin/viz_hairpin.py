# %% [markdown]
'''
> __Project__    : Injecting artificial vortices into channel flow
> __Author__     : Sugan Velazhagan
> __Date__       : 26/11/2026
> __Description__: Verify the hairpin vortex (Burgers core) on the staggered simulation grid:
>                  analytical vorticity vs curl of the injected velocity (relative L2 error), and 3D view.
'''
# %%
# ~  LOAD THE FRAMEWORK
# ~ ~~~~~~~~~~~~~~~~~~~~~
from hairpin_functions import *
ts = tt.time()

# %%
#* ==========================
#* INTERIOR SLICES
#* ==========================
    # --> Leave 1 ghost/halo cell on either side in all three directions
sl_in = (slice(1, -1), slice(1, -1), slice(1, -1))

# %%
#* ==========================
#* VISUALIZATION PARAMETERS
#* ==========================
om_choice  = 'om_a'                                  #^ Vorticity to visualize: 'om_a' (analytical) or 'om_n' (curl of u')
crop_fac   = 3.0                                     #^ Crop: hairpin bounding box + crop_fac * spread
iso_frac   = 0.3                                     #^ Initial iso-level, fraction of max|omega| in the crop (slider)

    # --> Camera view (fixed and reproducible; applied in this order: preset -> view_up -> azimuth/elevation/roll -> zoom)
view = dict(VIEW_DEFAULT)
view['preset']    = 'xy'         #^ Starting view: 'xy' (look down z, x-y plane), 'xz' (look down y), 'yz' (look down x), 'iso' (oblique)
view['up']        = (0, 1, 0)    #^ Direction pointing up in the image: (0, 1, 0) -> y+ up (wall at the bottom), (0, 0, 1) -> z+ up
view['azimuth']   = -40.0        #^ Rotate the camera about the up direction [deg]; < 0 swings it upstream (head in front, legs behind)
view['elevation'] = 25.0         #^ Tilt the camera up (> 0, look from above) or down (< 0) [deg]
view['roll']      = 0.0          #^ Spin the image in its own plane [deg]
view['zoom']      = 1.0          #^ > 1 zoom in, < 1 zoom out
view['exact']     = None         #^ Exact camera, overrides all of the above: [(eye), (focus), (up)] in wall units
                                 #^ To copy an interactive view: rotate it, print(pl.camera_position), paste here.

center = (x_c, y_c, z_c)
log.info(f"Head = ({center[0]:.4f}, {center[1]:.4f}, {center[2]:.4f})")

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
    #^ Screened core (shield radius win_factor * spread) and image hairpin as in the injected velocity
om_a = get_hairpin_vorticity(center=center, flag_image=flag_image, y_wall_near=y_wall_near)
omx += om_a[0]; omy += om_a[1]; omz += om_a[2]
om_a = (omx, omy, omz)

# %%
# ~ 3. VELOCITY -> CURL -> VORTICITY
# ~ ~~~~~~~~~~~~~~~~~~~~~
    #^ Injected field (screening and wall damping on); flag_localize=False in both calls -> unscreened tube
du, dv, dw = get_hairpin_velocity(center=center, flag_image=flag_image)
u += du; v += dv; w += dw
del du, dv, dw
om_n = get_vorticity(u, v, w)
log.info(f"max|div(u')| / (U_ref/spread) = {np.max(np.abs(get_divergence(u, v, w)[sl_in])) / (abs(U_ref) / spread):.3e}")

# %%
# ~ COMPARISON (relative L2 error, interior points)
# ~ ~~~~~~~~~~~~~~~~~~~~~
    # --> Core region: within 2 spreads of the untapered part of the centreline (|q| <= q_m / 2, g >= 0.98)
    #^ Near the tapered leg ends curl(u') = omega_a + grad(div A): the injected vorticity is the solenoidal
    #^ projection of the tapered tube (it spreads out at the leg ends), so omega_n != omega_a there by construction
X_core, _ = get_centreline(np.linspace(-0.5 * q_max, 0.5 * q_max, 400), center)
grids = [(xu, yu, zu), (xv, yv, zv), (xw, yw, zw)]

log.info('=' * 70)
log.info(f"{'comp.':<8}{'max|omega_a|':>14}{'max|omega_n|':>14}{'L2 err (all)':>14}{'L2 err (core)':>15}")
for n, (om_ref, om_num, comp) in enumerate(zip((omx, omy, omz), om_n, ('omega_x', 'omega_y', 'omega_z'))):
    Z, X, Y = np.meshgrid(grids[n][2], grids[n][0], grids[n][1], indexing='ij', sparse=True)
    core    = (get_distance2(X, Y, Z, X_core.T) < (2.0 * spread)**2)[sl_in]
    a, b    = om_ref[sl_in], om_num[sl_in]
    norm    = np.linalg.norm(a)
    err_all = np.linalg.norm(b - a) / norm if norm > 0 else np.nan
    norm_c  = np.linalg.norm(a[core])
    err_c   = np.linalg.norm(b[core] - a[core]) / norm_c if norm_c > 0 else np.nan
    log.info(f"{comp:<8}{np.max(np.abs(a)):>14.4e}{np.max(np.abs(b)):>14.4e}{err_all:>14.4f}{err_c:>15.4f}")
log.info('=' * 70)

# %%
# ~ 3D VISUALIZATION (PyVista): |omega| iso-surface colored by omega_z
# ~ ~~~~~~~~~~~~~~~~~~~~~
create_directory(root_out)
plot_vorticity_iso(om_a if om_choice == 'om_a' else om_n,
                   title=rf"{om_choice}: iso-surface of $\Vert\omega\Vert$, colored by $\omega_z$ (hairpin)",
                   fig_name=root_out + f'hairpin_iso_{om_choice}',
                   crop_fac=crop_fac, iso_frac=iso_frac, view=view)

log.info(f"Total time taken: {tt.time() - ts:.2f} seconds")
# %%
