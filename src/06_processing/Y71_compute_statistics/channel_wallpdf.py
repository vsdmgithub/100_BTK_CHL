# %% [markdown]
'''
> __Project__    : Data assimilation of turbulent channel flow DNS                  \
> __Author__     : Sugan Velazhagan                                                  \
> __Date__       : 06/22/2026                                                        \
> __Description__: Compute converged PDFs of wall quantities (tau_xy+, tau_zy+, p+)
>                  over the FULL channel domain.
>                  Averages over DNS frames 0..6400 in steps of 40 (161 frames).
>                  Histograms are accumulated per frame (memory-efficient) and
>                  normalised to PDFs after all frames are processed.
>                  Saves 'pdf_wallData_converged_FD.npz' under ground_truth/statistics/.
'''
# %%
# ~ LOAD THE FRAMEWORK
# ~ ~~~~~~~~~~~~~~~~~~~~~
from framework import *
from io_operations import read_ucvcwcPres_file
ts0 = tt.time()

# %%
# ~ FRAME RANGE FOR CONVERGED STATISTICS
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
frames_DNS_ch = np.arange(0, 6401, 40)    #^ 0, 40, 80, ..., 6400  →  161 frames
N_frames_ch   = len(frames_DNS_ch)

log.info(f"channel_wallpdf: full-domain converged wall PDF reference")
log.info(f"Frames: {frames_DNS_ch[0]}..{frames_DNS_ch[-1]}, step {frames_DNS_ch[1]-frames_DNS_ch[0]}, N={N_frames_ch}")

# %%
# ~ NORMALISATION AND BIN SETUP
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tau_norm = u_tau**2    #^ wall-unit reference  τ+ = τ / u_τ²

#* 200 bins per quantity — finer resolution than main.py (100 bins)
pdf_bins    = 200
bins_tau_xy = np.linspace(-1.0,  5.0, pdf_bins + 1) 
bins_tau_zy = np.linspace(-3.0, 3.0, pdf_bins + 1)  
bins_tau_p  = np.linspace(-12.0, 12.0, pdf_bins + 1) 

log.info(f"pdf_bins={pdf_bins}  |  tau_norm={tau_norm:.5e}")
log.info(f"bins_tau_xy: [{bins_tau_xy[0]:.2f}, {bins_tau_xy[-1]:.2f}]")
log.info(f"bins_tau_zy: [{bins_tau_zy[0]:.2f}, {bins_tau_zy[-1]:.2f}]")
log.info(f"bins_tau_p:  [{bins_tau_p[0]:.2f},  {bins_tau_p[-1]:.2f}]")

# %%
# ~ PER-FRAME HISTOGRAM COMPUTATION
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def process_frame(frame):
    """
    Read full-channel truth field at 'frame', extract wall-adjacent layer,
    normalise to wall units, and bin into histograms.

    Wall quantities (same convention as main.py):
        tau_xy+ = u[j=1] / (Re_cen * yu[1]) / u_τ²   (streamwise wall stress)
        tau_zy+ = w[j=1] / (Re_cen * yw[1]) / u_τ²   (spanwise   wall stress)
        p+      = p[j=1] / u_τ²  − ⟨p[j=1]⟩/u_τ²    (wall pressure fluctuation)

    Halo stripping mirrors main.py:
        U[1:-1, 2:-2, 1]  →  removes 1 z-ghost, 2 x-ghosts on each side, picks y=1
        W[2:-2, 1:-1, 1]  →  removes 2 z-ghosts, 1 x-ghost on each side, picks y=1
        P[1:-1, 1:-1, 1]  →  removes 1 ghost on each side in x and z, picks y=1

    Returns
    -------
    counts_xy, counts_zy, counts_p : ndarray (pdf_bins,)
        Histogram counts in each bin.
    """
    _, U, V, W, P = read_ucvcwcPres_file(field_folder_truth, ngx_FD, ngy_FD, ngz_FD, frame)
    P = P / dt

    #* Strip halos, take wall-adjacent y-level (j=1) — full (z, x) extent
    tau_xy_raw = U[1:-1, 2:-2, 1]    #^ (ngz_FD, ngx_FD)  staggered u near wall
    tau_zy_raw = W[2:-2, 1:-1, 1]    #^ (ngz_FD, ngx_FD)  staggered w near wall
    tau_p_raw  = P[1:-1, 1:-1, 1]    #^ (ngz_FD, ngx_FD)  pressure near wall

    #* Normalise to wall units  (yu[1], yw[1] from subdomain = same as FD, y-grid shared)
    tau_xy = (tau_xy_raw.ravel() / (Re_cen * yu[1])) / tau_norm
    tau_zy = (tau_zy_raw.ravel() / (Re_cen * yw[1])) / tau_norm
    tau_p  = tau_p_raw.ravel() / tau_norm
    tau_p -= np.mean(tau_p)

    #* Remove NaNs before binning
    tau_xy = tau_xy[~np.isnan(tau_xy)]
    tau_zy = tau_zy[~np.isnan(tau_zy)]
    tau_p  = tau_p[ ~np.isnan(tau_p )]

    #* Accumulate histogram counts (memory-efficient: avoids storing raw values)
    counts_xy, _ = np.histogram(tau_xy, bins=bins_tau_xy)
    counts_zy, _ = np.histogram(tau_zy, bins=bins_tau_zy)
    counts_p,  _ = np.histogram(tau_p,  bins=bins_tau_p )

    return counts_xy, counts_zy, counts_p

# %%
# ~ PARALLEL PROCESSING OVER FRAMES
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
log.info(f"Starting parallel histogram computation ...")
ts = tt.time()
results = Parallel(n_jobs=num_workers)(
    delayed(process_frame)(frame)
    for frame in tqdm(frames_DNS_ch, desc="Processing frames")
)
te = tt.time()
log.info(f"Processing time: {te - ts:.2f} s")

# %%
# ~ SUM HISTOGRAMS ACROSS FRAMES AND NORMALISE TO PDF
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
counts_xy = sum(r[0] for r in results)
counts_zy = sum(r[1] for r in results)
counts_p  = sum(r[2] for r in results)

def counts_to_pdf(counts, bins):
    """Normalise histogram counts to a PDF  (∫ PDF(x) dx = 1)."""
    bin_width = bins[1] - bins[0]
    centers   = 0.5 * (bins[:-1] + bins[1:])
    pdf       = counts / (counts.sum() * bin_width)
    return centers, pdf

centers_xy, pdf_xy = counts_to_pdf(counts_xy, bins_tau_xy)
centers_zy, pdf_zy = counts_to_pdf(counts_zy, bins_tau_zy)
centers_p,  pdf_p  = counts_to_pdf(counts_p,  bins_tau_p )

#* Parseval-style check: ∫ PDF dx ≈ 1
bw_xy = bins_tau_xy[1] - bins_tau_xy[0]
bw_zy = bins_tau_zy[1] - bins_tau_zy[0]
bw_p  = bins_tau_p[1]  - bins_tau_p[0]
log.info(f"[CHECK] ∫ PDF(tau_xy+) = {(pdf_xy * bw_xy).sum():.5f}  (should be 1)")
log.info(f"[CHECK] ∫ PDF(tau_zy+) = {(pdf_zy * bw_zy).sum():.5f}  (should be 1)")
log.info(f"[CHECK] ∫ PDF(p+)      = {(pdf_p  * bw_p ).sum():.5f}  (should be 1)")

# %%
# ~ SAVE .NPZ
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
save_loc_ch = root_out_mount + truth_dir + '/statistics/'
create_directory(save_loc_ch)
output_file = save_loc_ch + 'pdf_wallData_converged_FD.npz'

np.savez(output_file,
    #* ── PDF arrays: (pdf_bins,) — bin centres and normalised density ─────────
    tau_xy_bins = centers_xy,     #^ bin centres  τ_xy+
    tau_xy_pdf  = pdf_xy,         #^ PDF(τ_xy+)
    tau_zy_bins = centers_zy,     #^ bin centres  τ_zy+
    tau_zy_pdf  = pdf_zy,         #^ PDF(τ_zy+)
    tau_p_bins  = centers_p,      #^ bin centres  p+
    tau_p_pdf   = pdf_p,          #^ PDF(p+)
    #* ── Bin edges (for re-computing PDFs with different normalisation) ────────
    bins_tau_xy = bins_tau_xy,
    bins_tau_zy = bins_tau_zy,
    bins_tau_p  = bins_tau_p,
    #* ── Raw histogram counts (additive — useful for combining windows) ────────
    counts_xy   = counts_xy,
    counts_zy   = counts_zy,
    counts_p    = counts_p,
    #* ── Scaling ──────────────────────────────────────────────────────────────
    u_tau    = u_tau,
    delta_nu = delta_nu,
    tau_norm = tau_norm,           #^ u_τ²  — reference for normalisation
    #* ── Ensemble metadata ────────────────────────────────────────────────────
    pdf_bins   = pdf_bins,
    N_frames   = N_frames_ch,
    frames_DNS = frames_DNS_ch,
)
log.info(f"Wall PDF saved to: {output_file}")
log.info(f"Total time: {tt.time() - ts0:.2f} s")

# %%
