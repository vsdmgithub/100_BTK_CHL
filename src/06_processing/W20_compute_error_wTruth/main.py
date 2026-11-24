# %% [markdown]
'''
> __Project__    : Data assimilation with synthetic measurements of turbulent channel flow
> __Author__     : Sugan Velazhagan
> __Description__: Compute error and correlation with respect to the DNS truth data.
>                  Instantaneous y-profiles and y-averaged values are stored at every frame.
'''
# %%
# ~  LOAD THE FRAMEWORK
# ~ ~~~~~~~~~~~~~~~~~~~~~
from framework import *
ts = tt.time()

# %%
#* ==========================
#* INTERIOR SLICES
#* ==========================
    # --> Leave 1 ghost/halo cell on either side in all three directions
sl_in = (slice(1, -1), slice(1, -1), slice(1, -1))
yu_in = yu[1:-1]

    # --> y-locations of each component (v is wall-normal staggered -> 1 extra point)
y_comp   = [yu[1:-1], yv[1:-1], yw[1:-1], yu[1:-1]]  #^ u, v, w, p (p at cell centers)
n_y_comp = [len(y_c) for y_c in y_comp]

# %%
#* ==========================
#* SYNTHETIC MEAN PROFILE
#* ==========================
from preprocessing_operations import wall_model

    # --> Reichardt profile for Re_tau = 180 (wall distance from nearest wall)
y_wall    = np.minimum(yu_in, 2.0 - yu_in)
U_yu      = wall_model(y_wall / delta_nu)  #^ U+ (wall units)

# %%
# ~ COMPARE WITH TRUTH
# ~ ~~~~~~~~~~~~~~~~~~~~~
from io_operations import read_ucvcwcPres_file

horAxes = (0, 1)
allAxes = (0, 1, 2)

def process_error(ind, frame):
    """
    Compute error and correlation between assimilated and truth fields.

    Returns instantaneous y-profiles and y-averaged values for u, v, w, p.
    """
    #--> Read assimilated velocity field
    _, *uvwP = read_ucvcwcPres_file(root_sim, ngx, ngy, ngz, frame)
    u = uvwP[0][sl_in] / u_tau - U_yu[None, None, :] #^ Fluctuations (wall units)
    v = uvwP[1][sl_in] / u_tau
    w = uvwP[2][sl_in] / u_tau
    p = ( uvwP[3][sl_in] - np.mean(uvwP[3][sl_in]) ) / norm_pres #^ Gauge pressure removed (wall units)

        #--> Read truth velocity field
    _, *uvwP_r = read_ucvcwcPres_file(root_truth, ngx, ngy, ngz, frame+frame0_truth)
    u_r = uvwP_r[0][sl_in] / u_tau - U_yu[None, None, :]
    v_r = uvwP_r[1][sl_in] / u_tau
    w_r = uvwP_r[2][sl_in] / u_tau
    p_r = ( uvwP_r[3][sl_in] - np.mean(uvwP_r[3][sl_in]) ) / norm_pres

    components = [ (u, u_r), (v, v_r), (w, w_r), (p, p_r) ]

        #--> Compute error and correlation for each velocity/pressure component
    error_y = [np.zeros(n_y) for n_y in n_y_comp]  #^ one profile per component (different y-sizes)
    corr_y  = [np.zeros(n_y) for n_y in n_y_comp]
    error_m = np.zeros(4)
    corr_m  = np.zeros(4)

    for i, (vel, vel_r) in enumerate(components):

            #--> Instantaneous y-profile (absolute RMS error, wall units)
        error_y[i] = np.sqrt( np.mean((vel - vel_r)**2, axis=horAxes) )
        num_y      = np.sum(vel * vel_r, axis=horAxes)
        den_y      = np.sqrt( np.sum(vel**2, axis=horAxes) * np.sum(vel_r**2, axis=horAxes) )
        corr_y[i]  = np.divide(num_y, den_y, out=np.zeros_like(num_y), where=den_y > 0) #^ 0 where field vanishes (e.g. v at wall)

            #--> y-averaged
        error_m[i] = np.sqrt( np.mean((vel - vel_r)**2, axis=allAxes) )
        num_m      = np.sum(vel * vel_r, axis=allAxes)
        den_m      = np.sqrt( np.sum(vel**2, axis=allAxes) * np.sum(vel_r**2, axis=allAxes) )
        corr_m[i]  = num_m / den_m if den_m > 0 else 0.0 #^ 0 if field vanishes

    return ind, error_y, corr_y, error_m, corr_m

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
n_components = 4
comp_names   = ['u', 'v', 'w', 'p']
error_y_all  = [np.zeros((N_frames, n_y)) for n_y in n_y_comp]
corr_y_all   = [np.zeros((N_frames, n_y)) for n_y in n_y_comp]
error_m_all  = np.zeros((N_frames, n_components))
corr_m_all   = np.zeros((N_frames, n_components))

    #--> Aggregate results
for ind, err_y, cor_y, err_m, cor_m in results:
    for i in range(n_components):
        error_y_all[i][ind] = err_y[i]
        corr_y_all[i][ind]  = cor_y[i]
    error_m_all[ind] = err_m
    corr_m_all[ind]  = cor_m

te_processing = tt.time()
log.info(f"Error processing time: {te_processing - ts_processing:.2f} seconds")

# %%
#~ SAVE THE RESULTS
#~ ~~~~~~~~~~~~~~~~~~~~~
output_data = {
    # y-averaged values
    'error_m': error_m_all,        # (N_frames, 4) - u, v, w, p
    'corr_m' : corr_m_all,         # (N_frames, 4)

    # Metadata
    'U_mean'  : U_yu,
    'frames'  : frames_DNS,
    'time_p'  : time_p,
    'n_frames': N_frames,
}
    #--> Instantaneous y-profiles, per component (different y-sizes)
for i, name in enumerate(comp_names):
    output_data[f'error_y_{name}'] = error_y_all[i]       # (N_frames, n_y_{name})
    output_data[f'corr_y_{name}']  = corr_y_all[i]        # (N_frames, n_y_{name})
    output_data[f'yp_{name}']      = y_comp[i] / delta_nu
if flag_output:
    np.savez(root_out + 'error_wTruth.npz', **output_data)

te_total = tt.time()
log.info(" >>> Errors and correlations with truth saved ")
log.info(f"Total time taken: {te_total - ts:.2f} seconds")
# %%
