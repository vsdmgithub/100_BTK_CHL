# %% [markdown]
'''
> __Project__    : Injecting artificial vortices into channel flow
> __Author__     : Sugan Velazhagan
> __Date__       : 25/11/2026
> __Description__: Insert a Burgers vortex packet into a turbulent channel flow state:
>                  read the IC from the template package, add the vortex velocity, save the new IC in
>                  root_sim (same format), and visualize the |omega| iso-surface around the vortex.
>
>                  Usage: python main.py [sim_name]   (overrides settings.yaml)
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
#* VISUALIZATION PARAMETERS
#* ==========================
om_choice = 'om_fluct'     #^ Vorticity to visualize: 'om_fluct' (IC with vortex, x-z plane mean removed),
                           #^ 'om_total' (IC with vortex), 'om_vortex' (vortex only)
crop_fac  = 3.0            #^ Crop half-width: crop_fac * spread (transverse), crop_fac * length (axial)
iso_frac  = 0.3            #^ Initial iso-level, fraction of max|omega| in the crop (slider)
y_excl_plus = 30.0         #^ Exclude y+ < y_excl_plus from the crop (near-wall shear layers dominate max|omega|)
view      = dict(VIEW_DEFAULT)   #^ Camera view; change e.g. view['azimuth'] = 0.0 (see VIEW_DEFAULT in burgers_functions)

# %%
# ~ READ THE INITIAL CONDITION (template package, t = 0)
# ~ ~~~~~~~~~~~~~~~~~~~~~
from io_operations import read_ucvcwc_file, write_ucvcwc_file

_, *uvw = read_ucvcwc_file(root_package, ngx, ngy, ngz, frame_IC)
u, v, w = uvw
log.info(f"Read IC (frame {frame_IC}) from {root_package}: u {u.shape}, v {v.shape}, w {w.shape}")
log.info(f"Bulk velocity (approx.) = {np.mean(u[sl_in]):.4f}")

    # --> Replace the turbulent state by the mean profile: log law for u, v = w = 0 (shapes kept)
if flag_mean_padded:
    u = np.broadcast_to(get_mean_profile(yu)[None, None, :], u.shape).copy()
    v = np.zeros_like(v)
    w = np.zeros_like(w)
    log.info(f"Mean-padded IC (log law): U_cl = {u.max():.4f}, bulk velocity (approx.) = {np.mean(u[sl_in]):.4f}")

# %%
# ~ ADD THE BURGERS VORTEX PACKET
# ~ ~~~~~~~~~~~~~~~~~~~~~
    # --> Vortex velocity u' = curl(A) on the staggered grid (image vortex, window, wall damping)
center = (x_c, y_c, z_c)
du, dv, dw = get_burgers_velocity(center=center)
u += du; v += dv; w += dw

    # --> Diagnostics
div_p = get_divergence(du, dv, dw)[sl_in]
div_1 = get_divergence(u, v, w)[sl_in]
log.info(f"max|div(u')| / (U_ref/spread) = {np.max(np.abs(div_p)) / (abs(U_ref) / spread):.3e}")
log.info(f"max|div(u)| (IC with vortex)  = {np.max(np.abs(div_1)):.3e}")
log.info(f"max|u'|, |v'|, |w'| / u_tau   = {np.max(np.abs(du))/u_tau:.3f}, "
         f"{np.max(np.abs(dv))/u_tau:.3f}, {np.max(np.abs(dw))/u_tau:.3f}")

# %%
# ~ SAVE THE NEW INITIAL CONDITION (root_sim, same format)
# ~ ~~~~~~~~~~~~~~~~~~~~~
if flag_output:
        # --> Copy the template package (grid, inputs, run scripts, ...) except the velocity fields
    vel_files = [pre + "{:07d}".format(frame_IC) for pre in ('01uc00.dble.', '02vc00.dble.', '03wc00.dble.')]
    for fname in sorted(os.listdir(root_package)):
        if fname not in vel_files and os.path.isfile(root_package + fname):
            shutil.copy2(root_package + fname, root_sim + fname)

        # --> Write the velocity fields with the vortex, and verify by reading back
    write_ucvcwc_file(u, v, w, num=frame_IC, addr=root_sim)
    _, *uvw_chk = read_ucvcwc_file(root_sim, ngx, ngy, ngz, frame_IC)
    assert all(np.array_equal(a, b) for a, b in zip(uvw_chk, (u, v, w))), "Read-back mismatch"

        # --> Save the vortex parameters for reproducibility
    vortex_info = {key: settings[key] for key in ('grid_name', 'package_name', 'frame_IC', 'strength', 'sense', 'spread',
                                                  'length', 'y_core', 'wall', 'axis', 'win_factor', 'wall_damp', 'flag_image', 'flag_mean_padded')}
    vortex_info.update(center=[float(c) for c in center], U_ref=float(U_ref), Gamma=float(Gamma),
                       omega_0=float(omega_0), u_tau=float(u_tau), nu=float(nu))
    with open(root_sim + 'burgers_parameters.yaml', 'w') as file:
        yaml.safe_dump(vortex_info, file, sort_keys=False)
    log.info(f"New IC saved and verified in: {root_sim}")

# %%
# ~ 3D VISUALIZATION: |omega| iso-surface around the vortex, colored by omega_z
# ~ ~~~~~~~~~~~~~~~~~~~~~
if flag_plot:
    if om_choice == 'om_vortex':
        om = get_vorticity(du, dv, dw)
    else:
        om = get_vorticity(u, v, w)
        if om_choice == 'om_fluct':
            om = tuple(c - np.mean(c[1:-1, 1:-1, :], axis=(0, 1), keepdims=True) for c in om)   #^ Remove the mean (x-z average) profile
    create_directory(root_out)
    plot_vorticity_iso(om, center,
                       title=rf"{sim_name} ({om_choice}): iso-surface of $\Vert\omega\Vert$, colored by $\omega_z$",
                       fig_name=root_out + f'burgers_iso_{om_choice}_axis{axis}',
                       crop_fac=crop_fac, iso_frac=iso_frac, y_excl_plus=y_excl_plus, view=view)

log.info(f"Total time taken: {tt.time() - ts:.2f} seconds")
# %%
