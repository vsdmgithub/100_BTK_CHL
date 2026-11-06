# %% [markdown]
'''
> __Project__    : DNS with a channel subdomain
> __Author__     : Sugan Velazhagan                                                                        \
> __Date__       : 01/07/2025                                                                               \
> __Description__: Compute statistics in the subdomain
'''
# %% 
# ~  LOAD THE FRAMEWORK
#~ ~~~~~~~~~~~~~~~~~~~~~
from framework import * 
ts0=tt.time()
# %%
#* ==========================
#* FOR EXTENDED DOMAIN
#* ==========================
    # --> Subdomain velocity data extraction slices from extended grid
iu_sub = slice(istart+1, iend+3, 1)
iv_sub = slice(istart+1, iend+2, 1)
iw_sub = slice(istart+1, iend+2, 1)
ip_sub = slice(istart+1, iend+2, 1)
ju_sub = slice(0, jend + 2, 1)
jv_sub = slice(0, jend + 3, 1)
jw_sub = slice(0, jend + 2, 1)
jp_sub = slice(0, jend + 2, 1)
ku_sub = slice(kstart+1, kend+2, 1)
kv_sub = slice(kstart+1, kend+2, 1)
kw_sub = slice(kstart+1, kend+3, 1)
kp_sub = slice(kstart+1, kend+2, 1)
#* ==========================
# * FOR FULL DOMAIN
#* ==========================
     # --> Subdomain velocity data extraction slices from full channel grid
iu_sub_FD = slice(istart_meas+1, iend_meas+3, 1)
iv_sub_FD = slice(istart_meas+1, iend_meas+2, 1)
iw_sub_FD = slice(istart_meas+1, iend_meas+2, 1)
ip_sub_FD = slice(istart_meas+1, iend_meas+2, 1)
ju_sub_FD = slice(0, jend+2, 1)
jv_sub_FD = slice(0, jend+3, 1)
jw_sub_FD = slice(0, jend+2, 1)
jp_sub_FD = slice(0, jend+2, 1)
ku_sub_FD = slice(kstart_meas+1, kend_meas+2, 1)
kv_sub_FD = slice(kstart_meas+1, kend_meas+2, 1)
kw_sub_FD = slice(kstart_meas+1, kend_meas+3, 1)
kp_sub_FD = slice(kstart_meas+1, kend_meas+2, 1)

# ~ INITIALIZE ARRAYS
# ~ ~~~~~~~~~~~~~~~~~~~~~
#@ ngy_SD + 1 = jend + 2 => ngy_SD = jend + 1
ngy_SD = jend + 1
u_mean = np.zeros((ngy_SD-1))
v_mean = np.zeros((ngy_SD-2))
w_mean = np.zeros((ngy_SD-1))
p_mean = np.zeros((ngy_SD-1))
uu     = np.zeros((ngy_SD-1))
vv     = np.zeros((ngy_SD-2))
ww     = np.zeros((ngy_SD-1))
uv     = np.zeros((ngy_SD-1))
pp     = np.zeros((ngy_SD-1))
tau_xy_list = []
tau_zy_list = []

# %%
# ~ PROCESS STATISTICS
# ~ ~~~~~~~~~~~~~~~~~~~~~
from io_operations import read_ucvcwc_file, read_ucvcwcPres_file
from analysis import get_vel_and_p_statistics

def process_statistics(ind, frame, flag_truth=False):
    """
    Process velocity and pressure statistics for a single frame.
    
    Args:
        ind: Index of the current frame
        frame: Frame number to process
        
    Returns:
        tuple: (index, statistics_list, tau_xy, tau_zy)
    """
    if flag_truth:
        # Read velocity and pressure fields
        log.info(f"Processing truth frame: {frame}")
        _, *uvwP = read_ucvcwcPres_file(field_folder_truth, ngx_FD, ngy_FD, ngz_FD, frame)
        u = uvwP[0][ku_sub_FD, iu_sub_FD, ju_sub_FD]
        v = uvwP[1][kv_sub_FD, iv_sub_FD, jv_sub_FD]
        w = uvwP[2][kw_sub_FD, iw_sub_FD, jw_sub_FD]
        p = uvwP[3][ku_sub_FD, iw_sub_FD, ju_sub_FD]
    else:
            #--> Read assimilated velocity field
        _, *uvwP = read_ucvcwcPres_file(field_folder, ngx, ngy, ngz, frame)
        u = uvwP[0][ku_sub, iu_sub, ju_sub]
        v = uvwP[1][kv_sub, iv_sub, jv_sub]
        w = uvwP[2][kw_sub, iw_sub, jw_sub]
        p = uvwP[3][ku_sub, iw_sub, ju_sub]

    p = p / dt
    
    # Get statistics for current frame
    stats_list = get_vel_and_p_statistics(u, v, w, p, ind, interior_flag=True)

    # Compute shear stresses
    tau_xy = u[1:-1, 2:-2, 1].flatten()
    tau_zy = w[2:-2, 1:-1, 1].flatten()
    tau_p  = p[1:-1, 1:-1, 1].flatten()

    tau_p = tau_p - np.mean(tau_p)

    return ind, stats_list, tau_xy, tau_zy, tau_p

# =============================================================================
# MAIN PROCESSING LOOP
# =============================================================================
# Initialize frame-by-frame storage (3D array: [variable, frame, y-position])
n_vars = 9  # u_mean, v_mean, w_mean, p_mean, uu, vv, ww, uv, pp
frame_stats = np.zeros((n_vars, len(frames_DNS), ngy_SD-1))

tau_xy_list = []
tau_zy_list = []
tau_p_list  = []

ts=tt.time()
# Log the flag_truth value before parallel processing
log.info(f"flag_truth value before parallel processing: {flag_truth}")

# Process frames in parallel
results = Parallel(n_jobs=num_workers, verbose=10)(
    delayed(process_statistics)(ind, frame, flag_truth=flag_truth)
    for ind, frame in enumerate(tqdm(frames_DNS, desc="Processing frames"))
)

# Aggregate results
for res in results:
    ind, stats_list, tau_xy, tau_zy, tau_p  = res
    
    # Accumulate ensemble averages
    u_mean += stats_list[1]
    v_mean += stats_list[2]
    w_mean += stats_list[3] 
    p_mean += stats_list[4]
    uu += stats_list[5] 
    vv += stats_list[6] 
    ww += stats_list[7] 
    uv += stats_list[8] 
    pp += stats_list[9]
    
    # Store frame-by-frame statistics
    frame_stats[0, ind, :len(stats_list[1])] = stats_list[1]  # u_mean
    frame_stats[1, ind, :len(stats_list[2])] = stats_list[2]  # v_mean
    frame_stats[2, ind, :len(stats_list[3])] = stats_list[3]  # w_mean
    frame_stats[3, ind, :len(stats_list[4])] = stats_list[4]  # p_mean
    frame_stats[4, ind, :len(stats_list[5])] = stats_list[5]  # uu
    frame_stats[5, ind, :len(stats_list[6])] = stats_list[6]  # vv
    frame_stats[6, ind, :len(stats_list[7])] = stats_list[7]  # ww
    frame_stats[7, ind, :len(stats_list[8])] = stats_list[8]  # uv
    frame_stats[8, ind, :len(stats_list[9])] = stats_list[9]  # pp

    tau_xy_list.append(tau_xy)
    tau_zy_list.append(tau_zy)
    tau_p_list.append(tau_p)

te = tt.time()
log.info(f"Statistics processing time: {te-ts:.2f} seconds")

# =============================================================================
# POST-PROCESSING AND SAVING
# =============================================================================

# Normalize ensemble averages
N = N_frames
u_mean /= (N * u_tau)
v_mean /= (N * u_tau)
w_mean /= (N * u_tau)
uu /= (N * u_tau**2)
vv /= (N * u_tau**2)
ww /= (N * u_tau**2)
uv /= (N * u_tau**2)
pp /= (N * u_tau**4)
p_mean /= (N * u_tau**2)

# Normalize frame-by-frame statistics
frame_stats[0:3, :, :] /= u_tau       # Mean velocities
frame_stats[3:7, :, :] /= u_tau**2    # Second moments
frame_stats[7, :, :] /= u_tau**4      # Pressure variance
frame_stats[8, :, :] /= u_tau**2      # Pressure variance

# Prepare y-coordinates
yuIp = yu[1:ngy_SD]/delta_nu
yvIp = yv[2:ngy_SD]/delta_nu
ywIp = yw[1:ngy_SD]/delta_nu

# =============================================================================
# SAVE RESULTS
# =============================================================================
if flag_truth:
    stat_folder = root_out_mount + truth_dir + '/' + save_dir + '/'
else:
    stat_folder = save_loc

create_directory(stat_folder)

# Save ensemble averages as .npz file
np.savez(stat_folder + 'stats_ensemble.npz',
         u_mean=u_mean,
         v_mean=v_mean,
         w_mean=w_mean,
         uu=uu,
         vv=vv,
         ww=ww,
         uv=uv,
         pp=pp,
         p_mean=p_mean,
         yuIp=yuIp,
         yvIp=yvIp,
         ywIp=ywIp,
         N=N)

# Save frame-by-frame statistics
np.savez(stat_folder + 'stats_everyFrame.npz',
         u_mean=frame_stats[0],
         v_mean=frame_stats[1],
         w_mean=frame_stats[2],
         uu=frame_stats[3],
         vv=frame_stats[4],
         ww=frame_stats[5],
         uv=frame_stats[6],
         pp=frame_stats[7],
         p_mean=frame_stats[8],
         yuIp=yuIp,
         yvIp=yvIp,
         ywIp=ywIp,
         frames=frames)

log.info("Ensemble averages and frame-by-frame statistics saved")
# %%
# ~ PDF OF WALL STRESS
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~
from analysis import calculate_pdf

# Flatten all arrays after the loop
tau_xy_all = np.concatenate(tau_xy_list)
tau_zy_all = np.concatenate(tau_zy_list)
tau_p_all  = np.concatenate(tau_p_list)

# Remove NaNs
tau_xy_valid = tau_xy_all[~np.isnan(tau_xy_all)]
tau_zy_valid = tau_zy_all[~np.isnan(tau_zy_all)]
tau_p_valid  = tau_p_all[~np.isnan(tau_p_all)]

tau_norm = u_tau**2
tau_xy_valid = (tau_xy_valid / (Re_cen * yu[1]))/ tau_norm
tau_zy_valid = (tau_zy_valid / (Re_cen * yw[1]))/ tau_norm
tau_p_valid = (tau_p_valid / tau_norm)

pdf_bins = 50
bins_tau_xy = np.linspace(-1.0,5.0,pdf_bins+1)
bins_tau_zy = np.linspace(-3.0,3.0,pdf_bins+1)
bins_tau_p  = np.linspace(-12,12,pdf_bins+1)
tau_xy_pdf = calculate_pdf(tau_xy_valid, bins_tau_xy)
tau_zy_pdf = calculate_pdf(tau_zy_valid, bins_tau_zy)
tau_p_pdf  = calculate_pdf(tau_p_valid, bins_tau_p)

# Save tau PDFs as .npz file
np.savez(stat_folder + 'pdf_tau.npz',
         tau_xy_bins=tau_xy_pdf[:,0],
         tau_xy_pdf=tau_xy_pdf[:,1],
         tau_zy_bins=tau_zy_pdf[:,0],
         tau_zy_pdf=tau_zy_pdf[:,1],
         tau_p_bins=tau_p_pdf[:,0],
         tau_p_pdf=tau_p_pdf[:,1])
log.info("PDF of wall stress saved")

# plt.plot(tau_xy_pdf[:,0], tau_xy_pdf[:,1])
# plt.plot(tau_zy_pdf[:,0], tau_zy_pdf[:,1])
# plt.plot(tau_p_pdf[:,0], tau_p_pdf[:,1])
# plt.savefig(stat_folder+'pdf_tau.png',dpi=75,transparent=False,bbox_inches='tight')

te=tt.time()
log.info(f"Total time taken is {te-ts0} seconds")
# %%
