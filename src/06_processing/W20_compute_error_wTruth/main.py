# %% [markdown]
'''
> __Project__    : Data assimilation with synthetic measurements of turbulent channel flow
> __Author__     : Sugan Velazhagan
> __Description__: Compute error and correlation with respect to the DNS truth data.
>                  Errors and correlations are stored as RAW components for flexible post-processing.
'''
# %%
# ~  LOAD THE FRAMEWORK
# ~ ~~~~~~~~~~~~~~~~~~~~~
from framework import *
ts = tt.time()

# %%
#* ==========================
#* FOR EXTENDED DOMAIN
#* ==========================
    # --> Subdomain velocity data extraction slices from extended grid
iu_sub = slice(istart+1, iend+2, 1)
iv_sub = slice(istart+1, iend+2, 1)
iw_sub = slice(istart+1, iend+2, 1)
ip_sub = slice(istart+1, iend+2, 1)
ju_sub = slice(1, jend + 2, 1)
jv_sub = slice(1, jend + 2, 1)
jw_sub = slice(1, jend + 2, 1)
jp_sub = slice(1, jend + 2, 1)
ku_sub = slice(kstart+1, kend+2, 1)
kv_sub = slice(kstart+1, kend+2, 1)
kw_sub = slice(kstart+1, kend+2, 1)
kp_sub = slice(kstart+1, kend+2, 1)
#* ==========================
# * FOR FULL DOMAIN
#* ==========================
     # --> Subdomain velocity data extraction slices from full channel grid
iu_sub_FD = slice(istart_meas+1, iend_meas+2, 1)
iv_sub_FD = slice(istart_meas+1, iend_meas+2, 1)
iw_sub_FD = slice(istart_meas+1, iend_meas+2, 1)
ip_sub_FD = slice(istart_meas+1, iend_meas+2, 1)
ju_sub_FD = slice(1, jend+2, 1)
jv_sub_FD = slice(1, jend+2, 1)
jw_sub_FD = slice(1, jend+2, 1)
jp_sub_FD = slice(1, jend+2, 1)
ku_sub_FD = slice(kstart_meas+1, kend_meas+2, 1)
kv_sub_FD = slice(kstart_meas+1, kend_meas+2, 1)
kw_sub_FD = slice(kstart_meas+1, kend_meas+2, 1)
kp_sub_FD = slice(kstart_meas+1, kend_meas+2, 1)

# %%
# ~ COMPARE WITH TRUTH
# ~ ~~~~~~~~~~~~~~~~~~~~~
from io_operations import read_ucvcwcPres_file

horAxes = (0, 1)
spline = UnivariateSpline(yy, U_yy, k=3, s=0.0)
U_yu_sub = spline(yu[ju_sub])

def process_error(ind, frame):
    """
    Compute RAW error and correlation components between assimilated and truth fields.

    Returns RAW sums (before sqrt and normalization) for correct temporal averaging.
    """
    #--> Read assimilated velocity field
    _, *uvwP = read_ucvcwcPres_file(field_folder, ngx, ngy, ngz, frame)
    u = uvwP[0][ku_sub, iu_sub, ju_sub] - U_yu_sub[None, None, :] #^ Fluctuations
    v = uvwP[1][kv_sub, iv_sub, jv_sub]
    w = uvwP[2][kw_sub, iw_sub, jw_sub]
    G = np.mean(uvwP[3][kp_sub, ip_sub, jp_sub]) #^ Gauge pressure
    p = uvwP[3][kp_sub, ip_sub, jp_sub] - G

        #--> Read truth velocity field
    _, *uvwP_r = read_ucvcwcPres_file(field_folder_truth, ngx_FD, ngy_FD, ngz_FD, frame)
    u_r = uvwP_r[0][ku_sub_FD, iu_sub_FD, ju_sub_FD] - U_yu_sub[None, None, :]
    v_r = uvwP_r[1][kv_sub_FD, iv_sub_FD, jv_sub_FD]
    w_r = uvwP_r[2][kw_sub_FD, iw_sub_FD, jw_sub_FD]
    G   = np.mean(uvwP_r[3])
    p_r = uvwP_r[3][kp_sub_FD, ip_sub_FD, jp_sub_FD] - G

        #--> Interior slices (exclude boundary cells)
    # sl_u = (slice(1,-1), slice(1,-1), slice(1,-1))
    # sl_v = (slice(1,-1), slice(1,-1), slice(2,-1))
    # sl_w = (slice(1,-1), slice(1,-1), slice(1,-1))
    # sl_p = (slice(1,-1), slice(1,-1), slice(1,-1))

    # components = [
    #     (u[sl_u], u_r[sl_u]),
    #     (v[sl_v], v_r[sl_v]),
    #     (w[sl_w], w_r[sl_w]),
    #     (p[sl_p], p_r[sl_p]),
    # ]
    components = [ (u, u_r), (v, v_r), (w, w_r), (p, p_r) ]

        #--> Compute RAW error and correlation components for each velocity/pressure component
    n_y = jend + 1
    error_sum   = np.zeros((n_y, 4))
    cross_prod  = np.zeros((n_y, 4))
    sq_model    = np.zeros((n_y, 4))
    sq_truth    = np.zeros((n_y, 4))

    n_zEdge  = ( 2 - 1 ) * kskip + 1
    n_Inflow = ( 8 - 1 ) * iskip + 1
    for i, (vel, vel_r) in enumerate(components):

            #--> Exclude boundary cells
        if flag_cut_zEdge:
            vel   =   vel[n_zEdge:-n_zEdge, :, :]
            vel_r = vel_r[n_zEdge:-n_zEdge, :, :]
        if flag_cut_Inflow:
            vel   =   vel[:, n_Inflow:, :]
            vel_r = vel_r[:, n_Inflow:, :]

        error_sum[:, i]  = np.sum((vel - vel_r)**2, axis=horAxes)
        cross_prod[:, i] = np.sum(vel * vel_r,      axis=horAxes)
        sq_model[:, i]   = np.sum(vel**2,            axis=horAxes)
        sq_truth[:, i]   = np.sum(vel_r**2,          axis=horAxes)

    return ind, error_sum, cross_prod, sq_model, sq_truth

# %%
#~ MAIN PROCESSING LOOP
#~ ~~~~~~~~~~~~~~~~~~~~~
    #--> Process frames in parallel
ts_processing = tt.time()
results = Parallel(n_jobs=num_workers)(
    delayed(process_error)(ind, frame)
    for ind, frame in enumerate(tqdm(frames_DNS, desc="Processing frames"))
)

    #--> Preallocate arrays
n_vertical_points = jend + 1
n_components = 4
error_sum_all  = np.zeros((N_frames, n_vertical_points, n_components))
cross_prod_all = np.zeros((N_frames, n_vertical_points, n_components))
sq_model_all   = np.zeros((N_frames, n_vertical_points, n_components))
sq_truth_all   = np.zeros((N_frames, n_vertical_points, n_components))

    #--> Aggregate results
for ind, err_sum, cross, sq_m, sq_t in results:
    error_sum_all[ind]  = err_sum
    cross_prod_all[ind] = cross
    sq_model_all[ind]   = sq_m
    sq_truth_all[ind]   = sq_t

te_processing = tt.time()
log.info(f"Error processing time: {te_processing - ts_processing:.2f} seconds")

# %%
#~ SAVE THE RESULTS
#~ ~~~~~~~~~~~~~~~~~~~~~
yp_subdomain = yp[:jend+1]
output_data = {
    # RAW error components (sum of squared differences per y-level)
    'error_sum': error_sum_all,        # (N_frames, ngy-1, 4) - u, v, w, p

    # Cross-correlation products
    'cross_prod': cross_prod_all,      # (N_frames, ngy-1, 4)

    # Model squared terms
    'sq_model': sq_model_all,          # (N_frames, ngy-1, 4)

    # Truth squared terms
    'sq_truth': sq_truth_all,          # (N_frames, ngy-1, 4)

    # Metadata
    'yp_sub': yp_subdomain,
    'frames': frames_DNS,
    'n_frames': N_frames,
    'n_meas': (iend - istart + 1) * (kend - kstart + 1),
}
suffix = ''
if flag_cut_zEdge:
    suffix += '_cutEdge'
if flag_cut_Inflow:
    suffix += '_cutInflow'
np.savez(export_folder + f'error_wTruth_raw{suffix}.npz', **output_data)

te_total = tt.time()
log.info(" >>> RAW errors and correlations with truth saved ")
log.info(f"Total time taken: {te_total - ts:.2f} seconds")
# %%
