# %% [markdown]
'''
> __Project__    : Injecting artificial vortices into channel flow
> __Author__     : Sugan Velazhagan
> __Date__       : 25/11/2026
> __Description__: Insert multiple Burgers vortex packets into a channel flow state:
>                  random x-z centers with a minimum distance, random (Gaussian) height, strength, spread and length.
>                  Read the IC from the template package (or use the mean profile only), add the vortices,
>                  save the new IC in root_sim (same format), and visualize |omega| over the entire domain.
'''
# %%
# ~  LOAD THE FRAMEWORK
# ~ ~~~~~~~~~~~~~~~~~~~~~
from framework import *              #^ Adds func_path to sys.path
from burgers_functions import *      #^ Shared with V11_inject_syn_vortex; uses this folder's framework.py
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
om_choice   = 'om_fluct'   #^ Vorticity to visualize: 'om_fluct' (IC with vortices, x-z plane mean removed),
                           #^ 'om_total' (IC with vortices), 'om_vortex' (vortices only)
iso_frac    = 0.3          #^ Initial iso-level, fraction of max|omega| in the domain (slider)
y_excl_plus = 30.0 if not flag_mean_padded else 0.0   #^ Exclude y+ < y_excl_plus (turbulent near-wall layers dominate max|omega|)
view        = dict(VIEW_DEFAULT)   #^ Camera view; change e.g. view['preset'] = 'iso' (see VIEW_DEFAULT in burgers_functions)
view.update(azimuth=-35.0, elevation=35.0)   #^ Oblique view of the whole channel

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
# ~ ADD THE BURGERS VORTEX PACKETS
# ~ ~~~~~~~~~~~~~~~~~~~~~
    # --> Sum of the vortex velocities u' = curl(A) (each: window, wall damping, optional image)
du = np.zeros_like(u); dv = np.zeros_like(v); dw = np.zeros_like(w)
for n in range(n_vortex):
    du_n, dv_n, dw_n = get_burgers_velocity(center=(x_cs[n], y_cs[n], z_cs[n]), Gamma=Gammas[n],
                                            spread=spreads[n], length=lengths[n])
    du += du_n; dv += dv_n; dw += dw_n
    log.info(f"Vortex {n + 1:>3}/{n_vortex} added")
del du_n, dv_n, dw_n
u += du; v += dv; w += dw

    # --> Diagnostics
div_p = get_divergence(du, dv, dw)[sl_in]
div_1 = get_divergence(u, v, w)[sl_in]
log.info(f"max|div(u')|                  = {np.max(np.abs(div_p)):.3e}")
log.info(f"max|div(u)| (IC with vortices) = {np.max(np.abs(div_1)):.3e}")
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

        # --> Write the velocity fields with the vortices, and verify by reading back
    write_ucvcwc_file(u, v, w, num=frame_IC, addr=root_sim)
    _, *uvw_chk = read_ucvcwc_file(root_sim, ngx, ngy, ngz, frame_IC)
    assert all(np.array_equal(a, b) for a, b in zip(uvw_chk, (u, v, w))), "Read-back mismatch"

        # --> Save the settings and the drawn parameters of every vortex, for reproducibility
    vortex_info = {key: settings[key] for key in ('grid_name', 'package_name', 'frame_IC', 'n_vortex', 'min_dist',
                                                  'rng_seed', 'y_core_stat', 'strength_stat', 'spread_stat',
                                                  'length_stat', 'sense', 'wall', 'axis', 'win_factor', 'wall_damp',
                                                  'flag_image', 'flag_mean_padded')}
    vortex_info.update(u_tau=float(u_tau), nu=float(nu))
    vortex_info['vortices'] = [dict(center=[float(x_cs[n]), float(y_cs[n]), float(z_cs[n])],
                                    y_core=float(y_cores[n]), strength=float(strengths[n]),
                                    spread=float(spreads[n]), length=float(lengths[n]),
                                    U_ref=float(U_refs[n]), Gamma=float(Gammas[n]), omega_0=float(omega_0s[n]))
                               for n in range(n_vortex)]
    with open(root_sim + 'burgers_parameters.yaml', 'w') as file:
        yaml.safe_dump(vortex_info, file, sort_keys=False)
    log.info(f"New IC saved and verified in: {root_sim}")

# %%
# ~ 3D VISUALIZATION: |omega| iso-surface over the entire domain, colored by omega_z
# ~ ~~~~~~~~~~~~~~~~~~~~~
if flag_plot:
    if om_choice == 'om_vortex':
        om = get_vorticity(du, dv, dw)
    else:
        om = get_vorticity(u, v, w)
        if om_choice == 'om_fluct':
            om = tuple(c - np.mean(c[1:-1, 1:-1, :], axis=(0, 1), keepdims=True) for c in om)   #^ Remove the mean (x-z average) profile
    create_directory(root_out)
    center_dom = (0.5 * (xx[0] + xx[-1]), 0.5 * (y_wall_bot + y_wall_top), 0.5 * (zz[0] + zz[-1]))
    plot_vorticity_iso(om, center_dom,
                       title=rf"{sim_name} ({om_choice}, {n_vortex} vortices): iso-surface of $\Vert\omega\Vert$, colored by $\omega_z$",
                       fig_name=root_out + f'burgers_iso_{om_choice}_n{n_vortex}_axis{axis}',
                       iso_frac=iso_frac, y_excl_plus=y_excl_plus, view=view,
                       half=(0.5 * Lx, 0.5 * Ly, 0.5 * Lz))    #^ Entire domain (interior cell centers)

log.info(f"Total time taken: {tt.time() - ts:.2f} seconds")
# %%
