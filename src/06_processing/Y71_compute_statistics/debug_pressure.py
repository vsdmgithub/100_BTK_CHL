# %% [markdown]
'''
> __Project__    : DNS with a channel subdomain
> __Author__     : Sugan Velazhagan                                                                        \
> __Date__       : 08/07/2026                                                                               \
> __Description__: Diagnostic script for pressure decomposition at full-domain and subdomain levels.
                    Computes pressure mean and second moment of the decomposed subdomain pressure.
                    flag_truth selects the domain (FD or extended subdomain).
'''
# %%
from framework import *
from io_operations import read_ucvcwcPres_file
from analysis import decompose_pressure
from grid_operations import get_grid_center_wHalos

ts0 = tt.time()

# %%
# ~ DOMAIN-SPECIFIC SETUP BASED ON flag_truth
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if flag_truth:
    # Full DNS channel domain
    Ngx, Ngy, Ngz    = ngx_FD, ngy_FD, ngz_FD
    field_folder_use  = field_folder_truth
    xx_cc_full        = get_grid_center_wHalos(xx_FD)  # cell centers + halos, FD grid
    ip_sub = slice(istart_meas+1, iend_meas+2, 1)
    jp_sub = slice(0, jend+2, 1)
    kp_sub = slice(kstart_meas+1, kend_meas+2, 1)
else:
    # Extended subdomain
    Ngx, Ngy, Ngz    = ngx, ngy, ngz
    field_folder_use  = field_folder
    xx_cc_full        = xx_cc                           # cell centers + halos, subdomain grid
    ip_sub = slice(istart+1, iend+2, 1)
    jp_sub = slice(0, jend+2, 1)
    kp_sub = slice(kstart+1, kend+2, 1)

# x-coords for subdomain decompose_pressure: same slice as the pressure subdomain i-index
xx_cc_sub = xx_cc_full[ip_sub]

# y-dimension for pressure statistics
ngy_SD  = jend + 2          # ngy_SD - 1 = jend interior y-points after ghost-cell stripping
ngy_p   = ngy_SD - 2

log.info(f"flag_truth = {flag_truth}")
log.info(f"Domain: Ngx={Ngx}, Ngy={Ngy}, Ngz={Ngz}")
log.info(f"Pressure stats y-points: {ngy_p}")

# %%
# ~ INITIALIZE ACCUMULATOR ARRAYS
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p_mean = np.zeros(ngy_p)   # ensemble-accumulated pressure mean (decomposed)
pp     = np.zeros(ngy_p)   # ensemble-accumulated pressure second moment (decomposed)

# Temporal arrays: decomposition scalars at every time step
gauge_full_t = np.zeros(N_frames)   # gauge (spatial mean) of full-domain pressure
gradx_full_t = np.zeros(N_frames)   # mean x-pressure-gradient of full-domain
gauge_sub_t  = np.zeros(N_frames)   # gauge of subdomain pressure (before chop decomp)
gradx_sub_t  = np.zeros(N_frames)   # mean x-pressure-gradient of subdomain

# Temporal arrays: frame-by-frame pressure statistics of decomposed subdomain pressure
p_mean_t = np.zeros((N_frames, ngy_p))
pp_t     = np.zeros((N_frames, ngy_p))


# %%
# ~ FRAME PROCESSING FUNCTION
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~
def process_pressure_frame(ind, frame):
    """
    Process a single frame: decompose pressure at full-domain then subdomain level.

    Steps:
      1. Read full pressure (velocities read but discarded).
      2. Scale by 1/dt.
      3. Decompose full pressure → gauge_full, gradx_full.
      4. Chop to subdomain.
      5. Decompose subdomain pressure → gauge_sub, gradx_sub.
      6. Compute p_mean and pp from interior of decomposed subdomain pressure.

    Returns:
        (ind, gauge_full, gradx_full, gauge_sub, gradx_sub, p_mean_frame, pp_frame)
    """
    # Read velocities and pressure (velocities are intentionally ignored)
    _, _u, _v, _w, P_full = read_ucvcwcPres_file(field_folder_use, Ngx, Ngy, Ngz, frame)
    P_full = P_full / dt

    # --- Full-domain decomposition ---
    gauge_full = float(np.mean(P_full))
    _, gradx_full = decompose_pressure(
        P_full, xx_cc_full, Ngx+1, Ngy+1, Ngz+1, gradx_Pi=0.0)
    gradx_full = float(gradx_full)

    # --- Chop to subdomain ---
    P_sub = P_full[kp_sub, ip_sub, jp_sub]

    # --- Subdomain decomposition ---
    ngz_sub_p, ngx_sub_p, ngy_sub_p = P_sub.shape
    gauge_sub = float(np.mean(P_sub))
    p_sub_decomp, gradx_sub = decompose_pressure(
        P_sub, xx_cc_sub, ngx_sub_p, ngy_sub_p, ngz_sub_p, gradx_Pi=0.0)
    gradx_sub = float(gradx_sub)

    # --- Pressure statistics on interior of decomposed subdomain pressure ---
    p_int        = p_sub_decomp[1:-1, 1:-1, 1:-1]       # strip ghost cells: (k, i, j)
    p_mean_frame = np.mean(p_int, axis=(0, 1))           # average over x-z → (ngy_p,)
    pp_frame     = np.mean(p_int**2, axis=(0, 1))

    return ind, gauge_full, gradx_full, gauge_sub, gradx_sub, p_mean_frame, pp_frame


# %%
# ~ MAIN PROCESSING LOOP
# ~ ~~~~~~~~~~~~~~~~~~~~~~
ts = tt.time()
log.info(f"Processing {N_frames} frames with {num_workers} workers...")

results = Parallel(n_jobs=num_workers, verbose=10)(
    delayed(process_pressure_frame)(ind, frame)
    for ind, frame in enumerate(tqdm(frames_DNS, desc="Pressure frames"))
)

# Aggregate results
for res in results:
    ind, gf, df, gs, ds, pm, pp_fr = res
    gauge_full_t[ind]       = gf
    gradx_full_t[ind]       = df
    gauge_sub_t[ind]        = gs
    gradx_sub_t[ind]        = ds
    p_mean_t[ind, :len(pm)] = pm
    pp_t[ind, :len(pp_fr)]  = pp_fr
    p_mean                 += pm
    pp                     += pp_fr

te = tt.time()
log.info(f"Processing time: {te-ts:.2f} s")


# %%
# ~ NORMALIZE
# ~ ~~~~~~~~~~
N = N_frames

p_mean /= (N * norm_pres)
pp     /= (N * norm_pres**2)

p_mean_t /= norm_pres
pp_t     /= norm_pres**2

gauge_full_t /= norm_pres
gradx_full_t /= norm_pres
gauge_sub_t  /= norm_pres
gradx_sub_t  /= norm_pres

# y-coordinates in wall units (same convention as main.py)
yuIp = yu[1:ngy_SD] / delta_nu


# %%
# ~ SAVE
# ~ ~~~~~
if flag_truth:
    stat_folder = root_out_mount + truth_dir + '/' + save_dir + '/'
else:
    stat_folder = save_loc

create_directory(stat_folder)

np.savez(stat_folder + 'pressure_stats_diag.npz',
         p_mean       = p_mean,          # ensemble pressure mean  [norm_pres]
         pp           = pp,              # ensemble pressure second moment  [norm_pres^2]
         p_mean_t     = p_mean_t,        # per-frame pressure mean  (N_frames, ngy_p)  [norm_pres]
         pp_t         = pp_t,            # per-frame pressure second moment  (N_frames, ngy_p)  [norm_pres^2]
         gauge_full_t = gauge_full_t,    # gauge of full-domain pressure per frame  [norm_pres]
         gradx_full_t = gradx_full_t,    # mean x-gradient of full-domain pressure per frame  [norm_pres]
         gauge_sub_t  = gauge_sub_t,     # gauge of subdomain pressure per frame  [norm_pres]
         gradx_sub_t  = gradx_sub_t,     # mean x-gradient of subdomain pressure per frame  [norm_pres]
         yuIp         = yuIp,            # y-coordinates in wall units
         frames       = frames_DNS,
         N            = N)

log.info(f"Saved pressure_stats_diag.npz to {stat_folder}")

te = tt.time()
log.info(f"Total time: {te-ts0:.2f} s")
# %%
