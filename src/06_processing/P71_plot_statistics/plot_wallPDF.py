# %% [markdown]
'''
> __Project__    : Data assimilation with experimental measurements of turbulent channel flow \
> __Author__     : Sugan Velazhagan                                                              \
> __Date__       : 06/17/2026                                                                     \
> __Description__: Log-scale wall PDFs (tau_xy+, tau_zy+, p_wall+) — 1×3 figure.
>                  Col 1: tau_xy+   Col 2: tau_zy+   Col 3: p_wall+
>                  Lines: 4DVar (solid), truth DNS (dashed, toggle via flag_truth),
>                         JHTDB DNS (dashed grey, toggle via flag_JHTDB).
>                  PDFs computed by ref_pdf.py and stored as pdf_tau.npz.
>                  Keys: tau_xy_bins, tau_xy_pdf, tau_zy_bins, tau_zy_pdf,
>                        tau_p_bins, tau_p_pdf
'''

# %%
# ~  LOAD THE FRAMEWORK
# ~ ~~~~~~~~~~~~~~~~~~~~~
from framework import *
from scipy.integrate import trapezoid
ts = tt.time()

# %%
# ~ TOGGLES
# ~ ~~~~~~~~~~~~~~~~~~~~~
flag_truth  = True       #@ True: load and overlay truth (DNS) PDF
flag_JHTDB  = True       #@ True: load and overlay JHTDB DNS reference PDF
grad_method = 'second'   #@ JHTDB gradient method: 'first' or 'second'

# %%
# ~ LOGISTICS
# ~ ~~~~~~~~~~~~~~~~~~~~~
plot_dir = save_loc
create_directory(plot_dir)

#^ JHTDB PDF files — written by main.py / main_pressure.py in S45_spectral_analysis_JHTDB
jhtdb_dir = root_out + 'JHTDB_Retau_1000/'

log.info(f"Simulation     : {sim_name}")
log.info(f"flag_truth     : {flag_truth}")
log.info(f"flag_JHTDB     : {flag_JHTDB}")
log.info(f"Plot dir       : {plot_dir}")
log.info(f"JHTDB data dir : {jhtdb_dir}")

# %%
# ~ READ DA WALL PDF DATA
# ~ ~~~~~~~~~~~~~~~~~~~~~
log.info("Reading DA wall PDFs...")
pdf_file_DA = export_loc + 'pdf_tau.npz'

if not os.path.isfile(pdf_file_DA):
    log.error(f"DA wall PDF file not found: {pdf_file_DA}")
    log.error("Run ref_pdf.py (P71) first.")
    sys.exit(1)

da = np.load(pdf_file_DA)

bins_tau_xy_DA = da['tau_xy_bins']
pdf_tau_xy_DA  = da['tau_xy_pdf']
bins_tau_zy_DA = da['tau_zy_bins']
pdf_tau_zy_DA  = da['tau_zy_pdf']
bins_tau_p_DA  = da['tau_p_bins']
pdf_tau_p_DA   = da['tau_p_pdf']

log.info(f"DA PDFs loaded from: {pdf_file_DA}")

# %%
# ~ READ TRUTH PDF DATA  (optional, controlled by flag_truth)
# ~ ~~~~~~~~~~~~~~~~~~~~~
has_truth = False
bins_tau_xy_TR = bins_tau_zy_TR = bins_tau_p_TR = None
pdf_tau_xy_TR  = pdf_tau_zy_TR  = pdf_tau_p_TR  = None

if flag_truth:
    log.info("Reading truth wall PDFs...")
    pdf_file_TR = root_out_mount + truth_dir + '/' + save_dir + '/pdf_tau.npz'

    if not os.path.isfile(pdf_file_TR):
        log.warning(f"Truth wall PDF not found: {pdf_file_TR}")
        log.warning("Run ref_pdf.py with flag_truth=true first.")
    else:
        tr = np.load(pdf_file_TR)
        bins_tau_xy_TR = tr['tau_xy_bins']
        pdf_tau_xy_TR  = tr['tau_xy_pdf']
        bins_tau_zy_TR = tr['tau_zy_bins']
        pdf_tau_zy_TR  = tr['tau_zy_pdf']
        bins_tau_p_TR  = tr['tau_p_bins']
        pdf_tau_p_TR   = tr['tau_p_pdf']
        has_truth = True
        log.info(f"Truth PDFs loaded from: {pdf_file_TR}")

# %%
# ~ READ JHTDB REFERENCE PDFs  (stress + pressure)
# ~ ~~~~~~~~~~~~~~~~~~~~~
has_JHTDB = False
bins_tau_xy_JHTDB = bins_tau_zy_JHTDB = bins_p_JHTDB = None
pdf_tau_xy_JHTDB  = pdf_tau_zy_JHTDB  = pdf_p_JHTDB  = None

if flag_JHTDB:
    log.info("Reading JHTDB reference PDFs...")
    npz_jhtdb_xy = jhtdb_dir + f'pdf_tau_xy_{grad_method}order.npz'
    npz_jhtdb_zy = jhtdb_dir + f'pdf_tau_zy_{grad_method}order.npz'
    npz_jhtdb_p  = jhtdb_dir + 'pdf_p_wall.npz'

    missing = [f for f in [npz_jhtdb_xy, npz_jhtdb_zy, npz_jhtdb_p]
               if not os.path.isfile(f)]
    if missing:
        for f in missing:
            log.warning(f"JHTDB PDF not found: {f}")
        log.warning("Run main.py / main_pressure.py in S45_spectral_analysis_JHTDB.")
    else:
        jhtdb_xy = np.load(npz_jhtdb_xy)
        jhtdb_zy = np.load(npz_jhtdb_zy)
        jhtdb_p  = np.load(npz_jhtdb_p)

        bins_tau_xy_JHTDB = jhtdb_xy['bins']
        pdf_tau_xy_JHTDB  = jhtdb_xy['pdf']
        bins_tau_zy_JHTDB = jhtdb_zy['bins']
        pdf_tau_zy_JHTDB  = jhtdb_zy['pdf']
        bins_p_JHTDB      = jhtdb_p['bins']
        pdf_p_JHTDB       = jhtdb_p['pdf']

        has_JHTDB = True
        log.info(f"JHTDB tau_xy PDF: {npz_jhtdb_xy}  ({len(bins_tau_xy_JHTDB)} pts)")
        log.info(f"JHTDB tau_zy PDF: {npz_jhtdb_zy}  ({len(bins_tau_zy_JHTDB)} pts)")
        log.info(f"JHTDB p_wall PDF: {npz_jhtdb_p}   ({len(bins_p_JHTDB)} pts)")
        log.info(f"  Re_tau = {float(jhtdb_xy['Re_tau']):.1f},  "
                 f"u_tau = {float(jhtdb_xy['u_tau']):.5e},  "
                 f"gradient method = '{str(jhtdb_xy['gradient_method'])}'")

# %%
# ~ NORMALIZATION CHECKS
# ~ ~~~~~~~~~~~~~~~~~~~~~
log.info("=" * 60)
log.info("PDF NORMALIZATION CHECKS  (integral should be ≈ 1.0)")
log.info("=" * 60)

norm_pairs = [
    ("DA   tau_xy", bins_tau_xy_DA, pdf_tau_xy_DA),
    ("DA   tau_zy", bins_tau_zy_DA, pdf_tau_zy_DA),
    ("DA   p_wall", bins_tau_p_DA,  pdf_tau_p_DA),
]
if has_truth:
    norm_pairs += [
        ("TR   tau_xy", bins_tau_xy_TR, pdf_tau_xy_TR),
        ("TR   tau_zy", bins_tau_zy_TR, pdf_tau_zy_TR),
        ("TR   p_wall", bins_tau_p_TR,  pdf_tau_p_TR),
    ]
if has_JHTDB:
    norm_pairs += [
        ("JHTDB tau_xy", bins_tau_xy_JHTDB, pdf_tau_xy_JHTDB),
        ("JHTDB tau_zy", bins_tau_zy_JHTDB, pdf_tau_zy_JHTDB),
        ("JHTDB p_wall", bins_p_JHTDB,      pdf_p_JHTDB),
    ]

for label, b, p in norm_pairs:
    intg = trapezoid(p, b)
    log.info(f"  {label}: {intg:.6f}  (deviation {abs(1.0 - intg)*100:.3f}%)")

log.info("=" * 60)

# %%
# ~ PLOT SETTINGS
# ~ ~~~~~~~~~~~~~~~~~~~~~
from plot_functions import plot_settings_JFM, plot_settings_woLatex
if username == 'smuruga4':
    plot_settings_woLatex()
    mpl.rcParams["text.usetex"] = False
else:
    plot_settings_JFM()

#+ Colors
colB1   = '#03045e'   #+ deep navy  — tau_xy
colB2   = '#0077b6'   #+ ocean blue — tau_zy
colB3   = '#023e8a'   #+ dark royal blue — p_wall
col_TR  = '#d62828'   #+ red — ground truth
col_ref = '#404040'   #+ dark grey — JHTDB reference

#- Line weights and styles
dpi    = 300
lw_DA  = 3.0
lw_TR  = 4.5
lw_REF = 4.0
ls_DA  = '-'
ls_TR  = '--'
ls_REF = '--'

ls_label = 32
ts_tick  = 28

# %%
# ~ HELPER — shared axis styling (linear x, log y)
# ~ ~~~~~~~~~~~~~~~~~~~~~
def _style_ax_log(ax, x_major, x_minor):
    """Linear x-axis, log y-axis; in-pointing ticks on all four sides."""
    ax.xaxis.set_major_locator(ticker.MultipleLocator(x_major))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(x_minor))
    ax.tick_params(axis='x', which='major',
                   direction='in', colors='k', labelsize=ts_tick,
                   width=1.5, length=6, pad=6,
                   bottom=True, top=True)
    ax.tick_params(axis='x', which='minor',
                   direction='in', colors='k',
                   width=1.0, length=3, pad=3,
                   bottom=True, top=True)
    ax.tick_params(axis='y', which='major',
                   direction='in', colors='k', labelsize=ts_tick,
                   width=1.5, length=6, pad=6,
                   left=True, right=True)
    ax.tick_params(axis='y', which='minor',
                   direction='in', colors='k',
                   width=1.0, length=3, pad=3,
                   left=True, right=True)
    ax.grid(True, alpha=0.20, which='major')
    ax.grid(True, alpha=0.08, which='minor')

# %%
# ~ FIGURE: 1 × 3 LOG-SCALE SUBPLOTS  (tau_xy | tau_zy | p_wall)
# ~ ~~~~~~~~~~~~~~~~~~~~~
fig, (ax_xy, ax_zy, ax_p) = plt.subplots(1, 3, figsize=(24, 6))
fig.subplots_adjust(wspace=0.40)

# =============================================================================
# PANEL 1 — tau_xy+
# =============================================================================

if has_JHTDB:
    ax_xy.plot(bins_tau_xy_JHTDB, pdf_tau_xy_JHTDB,
               ls_REF, color=col_ref, lw=lw_REF,
               zorder=1, alpha=0.6, label='JHTDB DNS')

if has_truth:
    ax_xy.plot(bins_tau_xy_TR, pdf_tau_xy_TR,
               ls_TR, color=col_TR, lw=lw_TR, zorder=2, label='Truth (DNS)')

ax_xy.plot(bins_tau_xy_DA, pdf_tau_xy_DA,
           ls_DA, color=colB1, lw=lw_DA, zorder=3, label='4DVar')

ax_xy.set_yscale('log')
ax_xy.set_xlabel(r'$\tau_{xy}^+$',               fontsize=ls_label, labelpad=10)
ax_xy.set_ylabel(r'$\mathsf{PDF}(\tau_{xy}^+)$', fontsize=ls_label, labelpad=15)
ax_xy.set_xlim(-1.0, 5.0)
ax_xy.set_ylim(2e-4, 1e1)
ax_xy.set_yticks([1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1])

_style_ax_log(ax_xy, x_major=1.0, x_minor=0.2)

# =============================================================================
# PANEL 2 — tau_zy+
# =============================================================================

if has_JHTDB:
    ax_zy.plot(bins_tau_zy_JHTDB, pdf_tau_zy_JHTDB,
               ls_REF, color=col_ref, lw=lw_REF,
               zorder=1, alpha=0.6, label='JHTDB DNS')

if has_truth:
    ax_zy.plot(bins_tau_zy_TR, pdf_tau_zy_TR,
               ls_TR, color=col_TR, lw=lw_TR, zorder=2, label='Truth (DNS)')

ax_zy.plot(bins_tau_zy_DA, pdf_tau_zy_DA,
           ls_DA, color=colB2, lw=lw_DA, zorder=3, label='4DVar')

ax_zy.set_yscale('log')
ax_zy.set_xlabel(r'$\tau_{zy}^+$',               fontsize=ls_label, labelpad=10)
ax_zy.set_ylabel(r'$\mathsf{PDF}(\tau_{zy}^+)$', fontsize=ls_label, labelpad=15)
ax_zy.set_xlim(-3.0, 3.0)
ax_zy.set_ylim(2e-4, 1e1)
ax_zy.set_yticks([1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1])

_style_ax_log(ax_zy, x_major=1.0, x_minor=0.2)

# =============================================================================
# PANEL 3 — p_wall+  (DA + truth + JHTDB; pressure has no wall model)
# =============================================================================

if has_JHTDB:
    ax_p.plot(bins_p_JHTDB, pdf_p_JHTDB,
              ls_REF, color=col_ref, lw=lw_REF,
              zorder=1, alpha=0.6, label='JHTDB DNS')

if has_truth:
    ax_p.plot(bins_tau_p_TR, pdf_tau_p_TR,
              ls_TR, color=col_TR, lw=lw_TR, zorder=2, label='Truth (DNS)')

ax_p.plot(bins_tau_p_DA, pdf_tau_p_DA,
          ls_DA, color=colB3, lw=lw_DA, zorder=3, label='4DVar')

ax_p.set_yscale('log')
ax_p.set_xlabel(r'$p_w^+$',               fontsize=ls_label, labelpad=10)
ax_p.set_ylabel(r'$\mathsf{PDF}(p_w^+)$', fontsize=ls_label, labelpad=15)
ax_p.set_xlim(-12.0, 12.0)
ax_p.set_ylim(2e-4, 1e0)
ax_p.set_yticks([1e-4, 1e-3, 1e-2, 1e-1, 1e0])

_style_ax_log(ax_p, x_major=4.0, x_minor=1.0)

plt.tight_layout()

save_name = f'wallPDF_{sim_name}.png'
plt.savefig(plot_dir + save_name,
            dpi=dpi, bbox_inches='tight',
            transparent=False, facecolor='w', edgecolor='k')
plt.show()
log.info(f"Saved: {plot_dir}{save_name}")
plt.close()

# %%
# ~ SUMMARY
# ~ ~~~~~~~~~~~~~~~~~~~~~
te = tt.time()
log.info("=" * 60)
log.info(f"Total time        : {te - ts:.2f} s")
log.info(f"Plot saved to     : {plot_dir}{save_name}")
log.info(f"JHTDB pdf dir     : {jhtdb_dir}")
log.info(f"flag_truth        : {flag_truth}  (has_truth={has_truth})")
log.info(f"flag_JHTDB        : {flag_JHTDB}  (has_JHTDB={has_JHTDB})")
log.info(f"grad_method       : '{grad_method}'")
log.info("=" * 60)

# %%
