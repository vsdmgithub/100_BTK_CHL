# %% [markdown]
'''
> __Project__    : Data assimilation with experimental measurements of turbulent channel flow \
> __Author__     : Sugan Velazhagan                                                              \
> __Date__       : 06/22/2026                                                                     \
> __Description__: Export publication-ready wall PDFs (tau_xy+, tau_zy+, p_wall+) as HDF5.
>                  Four separate files:
>                    wallPDF_convFD.h5          — converged full-domain DNS (flag_conv_FD)
>                    wallPDF_Truth.h5           — truth single-window DNS   (flag_truth)
>                    wallPDF_{sim_name}.h5      — DA assimilated (always; rename via settings)
>                  Source PDFs from ref_pdf.py (P71) and channel_wallpdf.py (Y71).
>                  Keys per file: tau_xy_bins/pdf, tau_zy_bins/pdf, tau_p_bins/pdf
>                  Toggle flag_conv_FD and flag_truth to skip already-exported reference data.
>                  Re-run with different sim_suffix in settings.yaml for each DA simulation.
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
flag_conv_FD = False    #@ True: export converged full-domain DNS reference PDF
flag_truth   = False    #@ True: export truth single-window DNS PDF

# %%
# ~ EXPORT DIRECTORY
# ~ ~~~~~~~~~~~~~~~~~~~~~
pub_dir = base_data + 'publication_data/'
create_directory(pub_dir)
export_date = datetime.now().strftime('%Y-%m-%d')

log.info(f"Simulation     : {sim_name}")
log.info(f"flag_conv_FD   : {flag_conv_FD}")
log.info(f"flag_truth     : {flag_truth}")
log.info(f"Export dir     : {pub_dir}")

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
    hf.attrs['units']       = ('Wall units: tau+ = tau/u_tau^2, p+ = p/u_tau^2. '
                               'bins are bin-centre values; pdf integrates to ≈ 1.')
    hf.attrs['description'] = description

def _write_pdf_datasets(hf, bins_xy, pdf_xy, bins_zy, pdf_zy, bins_p, pdf_p):
    hf.create_dataset('tau_xy_bins', data=bins_xy, compression='gzip')
    hf.create_dataset('tau_xy_pdf',  data=pdf_xy,  compression='gzip')
    hf.create_dataset('tau_zy_bins', data=bins_zy, compression='gzip')
    hf.create_dataset('tau_zy_pdf',  data=pdf_zy,  compression='gzip')
    hf.create_dataset('tau_p_bins',  data=bins_p,  compression='gzip')
    hf.create_dataset('tau_p_pdf',   data=pdf_p,   compression='gzip')

def _verify(h5_path):
    with h5py.File(h5_path, 'r') as hf:
        for ds_name in hf:
            item = hf[ds_name]
            if isinstance(item, h5py.Dataset):
                log.info(f"  {ds_name:<18}: shape={item.shape}")

# %%
# ~ READ + EXPORT 1 — CONVERGED FULL-DOMAIN DNS REFERENCE
# ~ ~~~~~~~~~~~~~~~~~~~~~
if flag_conv_FD:
    log.info("Reading converged full-domain reference PDFs...")
    pdf_file_conv = root_out_mount + truth_dir + '/statistics/pdf_wallData_converged_FD.npz'

    if not os.path.isfile(pdf_file_conv):
        log.error(f"Converged FD PDF not found: {pdf_file_conv}")
        log.error("Run channel_wallpdf.py (Y71) first.")
        sys.exit(1)

    conv = np.load(pdf_file_conv)
    log.info(f"  N_frames={int(conv['N_frames'])},  pdf_bins={int(conv['pdf_bins'])},  "
             f"u_tau={float(conv['u_tau']):.5e}")

    h5_conv = pub_dir + 'wallPDF_convFD.h5'
    log.info(f"Writing: {h5_conv}")
    with h5py.File(h5_conv, 'w') as hf:
        _write_meta(hf, 'Converged full-domain DNS wall PDFs in wall units. '
                        '161 frames, full channel, time-averaged until statistical convergence.')
        hf.attrs['source']    = pdf_file_conv
        hf.attrs['N_frames']  = int(conv['N_frames'])
        hf.attrs['pdf_bins']  = int(conv['pdf_bins'])
        hf.attrs['u_tau_conv'] = float(conv['u_tau'])
        _write_pdf_datasets(hf,
                            conv['tau_xy_bins'], conv['tau_xy_pdf'],
                            conv['tau_zy_bins'], conv['tau_zy_pdf'],
                            conv['tau_p_bins'],  conv['tau_p_pdf'])
    log.info(f"convFD written: {h5_conv}")
    _verify(h5_conv)

# %%
# ~ READ + EXPORT 2 — TRUTH (SINGLE-WINDOW DNS)
# ~ ~~~~~~~~~~~~~~~~~~~~~
if flag_truth:
    log.info("Reading truth wall PDFs...")
    pdf_file_TR = root_out_mount + truth_dir + '/' + save_dir + '/pdf_tau.npz'

    if not os.path.isfile(pdf_file_TR):
        log.error(f"Truth wall PDF not found: {pdf_file_TR}")
        log.error("Run ref_pdf.py (P71) with flag_truth=true first.")
        sys.exit(1)

    tr = np.load(pdf_file_TR)

    h5_tr = pub_dir + 'wallPDF_Truth.h5'
    log.info(f"Writing: {h5_tr}")
    with h5py.File(h5_tr, 'w') as hf:
        _write_meta(hf, 'Truth single-window DNS wall PDFs in wall units. '
                        'Computed from the same window length as the DA simulation.')
        hf.attrs['source'] = pdf_file_TR
        _write_pdf_datasets(hf,
                            tr['tau_xy_bins'], tr['tau_xy_pdf'],
                            tr['tau_zy_bins'], tr['tau_zy_pdf'],
                            tr['tau_p_bins'],  tr['tau_p_pdf'])
    log.info(f"Truth written: {h5_tr}")
    _verify(h5_tr)

# %%
# ~ READ + EXPORT 3 — DA ASSIMILATED  (named by sim_name)
# ~ ~~~~~~~~~~~~~~~~~~~~~
log.info("Reading DA wall PDFs...")
pdf_file_DA = export_loc + 'pdf_tau.npz'

if not os.path.isfile(pdf_file_DA):
    log.error(f"DA wall PDF not found: {pdf_file_DA}")
    log.error("Run ref_pdf.py (P71) first.")
    sys.exit(1)

da = np.load(pdf_file_DA)

h5_da = pub_dir + f'wallPDF_{sim_name}.h5'
log.info(f"Writing: {h5_da}")
with h5py.File(h5_da, 'w') as hf:
    _write_meta(hf, '4DVar assimilated wall PDFs in wall units.')
    hf.attrs['sim_name'] = sim_name
    hf.attrs['source']   = pdf_file_DA
    _write_pdf_datasets(hf,
                        da['tau_xy_bins'], da['tau_xy_pdf'],
                        da['tau_zy_bins'], da['tau_zy_pdf'],
                        da['tau_p_bins'],  da['tau_p_pdf'])
log.info(f"DA written: {h5_da}")
_verify(h5_da)

# %%
# ~ SUMMARY
# ~ ~~~~~~~~~~~~~~~~~~~~~
te = tt.time()
log.info('=' * 70)
log.info('EXPORT SUMMARY')
log.info('=' * 70)
if flag_conv_FD:
    log.info(f"  convFD : {h5_conv}")
if flag_truth:
    log.info(f"  Truth  : {h5_tr}")
log.info(f"  DA     : {h5_da}")
log.info(f"Total time : {te - ts:.2f} s")
log.info('=' * 70)

# %%
