# %% [markdown]
'''
> __Project__    : Data assimilation with experimental measurements of turbulent channel flow \
> __Author__     : Sugan Velazhagan                                                              \
> __Date__       : 06/22/2026                                                                     \
> __Description__: Export publication-ready statistics in wall units as HDF5.
>                  Exports: DNS reference, truth DNS, and DA-assimilated statistics.
>                  Processes second moments into fluctuating components:
>                    uu_fluc = uu - U^2,  vv_fluc = vv - V^2,  ww_fluc = ww - W^2
>                    uv_fluc = uv[:-1] - 0.5*(U[1:]+U[:-1])*V
>                    p_rms   = sqrt(pp - p_mean^2)
>                  Toggle flag_dns and flag_truth to skip already-exported reference data.
>                  DA output file is named by sim_name.
'''

# %%
# ~ LOAD THE FRAMEWORK
# ~ ~~~~~~~~~~~~~~~~~~~~~
from framework import *
import h5py
from datetime import datetime
ts = tt.time()

# %%
# ~ TOGGLES
# ~ ~~~~~~~~~~~~~~~~~~~~~
flag_dns   = False    #@ True: export DNS reference statistics
flag_truth = False

# %%
# ~ EXPORT DIRECTORY
# ~ ~~~~~~~~~~~~~~~~~~~~~
pub_dir = base_data + 'publication_data/'
create_directory(pub_dir)
export_date = datetime.now().strftime('%Y-%m-%d')

log.info(f"Simulation     : {sim_name}")
log.info(f"flag_dns       : {flag_dns}")
log.info(f"flag_truth     : {flag_truth}")
log.info(f"Export dir     : {pub_dir}")

# %%
# ~ HELPER — compute fluctuating Reynolds stresses from second moments
# ~ ~~~~~~~~~~~~~~~~~~~~~
def get_fluctuations(sd):
    """
    From a stats_ensemble dict (keys: u_mean, v_mean, w_mean, uu, vv, ww, uv, pp, p_mean),
    return (U, uu_f, vv_f, ww_f, uv_f, p_rms) all in wall units.
    uv_f is on the v-grid (one element shorter than u-grid).
    """
    U      = sd['u_mean']
    V      = sd['v_mean']
    W      = sd['w_mean']
    uu_f   = sd['uu']    - U**2
    vv_f   = sd['vv']    - V**2
    ww_f   = sd['ww']    - W**2
    uv_f   = sd['uv'][:-1] - 0.5 * (U[1:] + U[:-1]) * V
    pp_f   = sd['pp']    - sd['p_mean']**2
    p_rms  = np.sqrt(np.clip(pp_f, 0.0, None))
    return U, uu_f, vv_f, ww_f, uv_f, p_rms

def load_stats(npz_file):
    data = np.load(npz_file, allow_pickle=True)
    sd   = {key: data[key] for key in data.files}
    data.close()
    return sd

# %%
# ~ READ DNS REFERENCE STATISTICS
# ~ ~~~~~~~~~~~~~~~~~~~~~
if flag_dns:
    log.info("Reading DNS reference statistics...")
    ngy_SD       = 193
    mean_profile = np.loadtxt(root_out + 'statistics/avg_yp.dat')
    yp_dns       = mean_profile[:ngy_SD, 0]
    u_mean_dns   = mean_profile[:ngy_SD, 1]
    var_profile  = np.loadtxt(root_out + 'statistics/var_yp.dat')
    uu_dns       = var_profile[:ngy_SD, 1]
    vv_dns       = var_profile[:ngy_SD, 2]
    ww_dns       = var_profile[:ngy_SD, 3]
    uv_dns       = var_profile[:ngy_SD, 4]
    pp_dns       = var_profile[:ngy_SD, 5]
    p_rms_dns    = np.sqrt(np.clip(pp_dns, 0.0, None))
    log.info(f"DNS: {ngy_SD} y-points, y+ = [{yp_dns[0]:.3f}, {yp_dns[-1]:.1f}]")

# %%
# ~ READ TRUTH STATISTICS
# ~ ~~~~~~~~~~~~~~~~~~~~~
if flag_truth:
    log.info("Reading truth statistics...")
    truth_file = root_out_mount + truth_dir + '/' + save_dir + '/stats_ensemble.npz'
    if not os.path.isfile(truth_file):
        log.error(f"Truth statistics not found: {truth_file}")
        sys.exit(1)
    sd_TR          = load_stats(truth_file)
    yuIp_TR        = sd_TR['yuIp']
    yvIp_TR        = sd_TR['yvIp']
    U_TR, uu_TR, vv_TR, ww_TR, uv_TR, p_rms_TR = get_fluctuations(sd_TR)
    log.info(f"Truth: {len(yuIp_TR)} u-grid points, y+ = [{yuIp_TR[0]:.3f}, {yuIp_TR[-1]:.1f}]")

# %%
# ~ READ DA STATISTICS
# ~ ~~~~~~~~~~~~~~~~~~~~~
log.info("Reading DA statistics...")
da_file = export_loc + 'stats_ensemble.npz'
if not os.path.isfile(da_file):
    log.error(f"DA statistics not found: {da_file}")
    sys.exit(1)
sd_DA          = load_stats(da_file)
yuIp_DA        = sd_DA['yuIp']
yvIp_DA        = sd_DA['yvIp']
U_DA, uu_DA, vv_DA, ww_DA, uv_DA, p_rms_DA = get_fluctuations(sd_DA)
N_frames_DA    = int(sd_DA.get('N', 0))
log.info(f"DA ({sim_name}): {len(yuIp_DA)} u-grid points, y+ = [{yuIp_DA[0]:.3f}, {yuIp_DA[-1]:.1f}]")

# %%
# ~ SHARED METADATA HELPER
# ~ ~~~~~~~~~~~~~~~~~~~~~
def _write_meta(hf, description):
    hf.attrs['project']     = '32_DA_CHL'
    hf.attrs['Re_tau']      = float(Re_tau)
    hf.attrs['Re_cen']      = float(Re_cen)
    hf.attrs['u_tau']       = float(u_tau)
    hf.attrs['delta_nu']    = float(delta_nu)
    hf.attrs['export_date'] = export_date
    hf.attrs['units']       = ('Wall units: y+=y/delta_nu, U+=U/u_tau, '
                               'R_ij+=R_ij/u_tau^2, p_rms+=p_rms/u_tau^2')
    hf.attrs['description'] = description

def _verify(h5_path):
    with h5py.File(h5_path, 'r') as hf:
        for ds_name in hf:
            item = hf[ds_name]
            if hasattr(item, 'keys'):   # group
                shapes = {k: item[k].shape[0] for k in item.keys()}
                log.info(f"  /{ds_name:<10}: {shapes}")
            else:                        # top-level dataset
                log.info(f"  {ds_name:<12}: shape={item.shape}")

# %%
# ~ EXPORT 1 — DNS REFERENCE
# ~ ~~~~~~~~~~~~~~~~~~~~~
if flag_dns:
    h5_dns = pub_dir + 'statistics_DNS.h5'
    log.info(f"Writing DNS HDF5: {h5_dns}")
    with h5py.File(h5_dns, 'w') as hf:
        _write_meta(hf, 'Global DNS reference statistics in wall units. '
                        'p_rms = sqrt(pp), where pp is the pressure variance from the .dat file.')
        hf.create_dataset('yp',    data=yp_dns,    compression='gzip')
        hf.create_dataset('U',     data=u_mean_dns, compression='gzip')
        hf.create_dataset('uu',    data=uu_dns,    compression='gzip')
        hf.create_dataset('vv',    data=vv_dns,    compression='gzip')
        hf.create_dataset('ww',    data=ww_dns,    compression='gzip')
        hf.create_dataset('uv',    data=uv_dns,    compression='gzip')
        hf.create_dataset('p_rms', data=p_rms_dns, compression='gzip')
    log.info(f"DNS written: {h5_dns}  ({len(yp_dns)} points)")
    _verify(h5_dns)

# %%
# ~ EXPORT 2 — TRUTH (DNS GROUND-TRUTH)
# ~ ~~~~~~~~~~~~~~~~~~~~~
if flag_truth:
    h5_tr = pub_dir + 'statistics_Truth.h5'
    log.info(f"Writing Truth HDF5: {h5_tr}")
    with h5py.File(h5_tr, 'w') as hf:
        _write_meta(hf, 'Truth DNS statistics in wall units. '
                        'uu/vv/ww/p_rms on u-grid (yuIp); uv on v-grid (yvIp). '
                        'Second moments converted to fluctuating components.')
        hf.attrs['source'] = truth_file
        hf.create_dataset('yuIp',  data=yuIp_TR,  compression='gzip')
        hf.create_dataset('yvIp',  data=yvIp_TR,  compression='gzip')
        hf.create_dataset('U',     data=U_TR,     compression='gzip')
        hf.create_dataset('uu',    data=uu_TR,    compression='gzip')
        hf.create_dataset('vv',    data=vv_TR,    compression='gzip')
        hf.create_dataset('ww',    data=ww_TR,    compression='gzip')
        hf.create_dataset('uv',    data=uv_TR,    compression='gzip')
        hf.create_dataset('p_rms', data=p_rms_TR, compression='gzip')
    log.info(f"Truth written: {h5_tr}  ({len(yuIp_TR)} u-grid, {len(yvIp_TR)} v-grid points)")
    _verify(h5_tr)

# %%
# ~ EXPORT 3 — DA ASSIMILATED
# ~ ~~~~~~~~~~~~~~~~~~~~~
h5_da = pub_dir + f'statistics_{sim_name}.h5'
log.info(f"Writing DA HDF5: {h5_da}")
with h5py.File(h5_da, 'w') as hf:
    _write_meta(hf, '4DVar assimilated statistics in wall units. '
                    'uu/vv/ww/p_rms on u-grid (yuIp); uv on v-grid (yvIp). '
                    'Second moments converted to fluctuating components.')
    hf.attrs['sim_name'] = sim_name
    hf.attrs['N_frames'] = N_frames_DA
    hf.attrs['source']   = da_file
    hf.create_dataset('yuIp',  data=yuIp_DA,  compression='gzip')
    hf.create_dataset('yvIp',  data=yvIp_DA,  compression='gzip')
    hf.create_dataset('U',     data=U_DA,     compression='gzip')
    hf.create_dataset('uu',    data=uu_DA,    compression='gzip')
    hf.create_dataset('vv',    data=vv_DA,    compression='gzip')
    hf.create_dataset('ww',    data=ww_DA,    compression='gzip')
    hf.create_dataset('uv',    data=uv_DA,    compression='gzip')
    hf.create_dataset('p_rms', data=p_rms_DA, compression='gzip')
log.info(f"DA written: {h5_da}  ({len(yuIp_DA)} u-grid, {len(yvIp_DA)} v-grid points)")
_verify(h5_da)

te = tt.time()
log.info('=' * 70)
log.info('EXPORT SUMMARY')
log.info('=' * 70)
if flag_dns:
    log.info(f"  DNS    : {h5_dns}")
if flag_truth:
    log.info(f"  Truth  : {h5_tr}")
log.info(f"  DA     : {h5_da}")
log.info(f"Total time : {te - ts:.2f} s")
log.info('=' * 70)

# %%
