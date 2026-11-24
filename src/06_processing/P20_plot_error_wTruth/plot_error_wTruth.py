# %% [markdown]
"""
Plot publication-ready error and correlation against DNS truth (JFM style).

Figure 1 — temporal evolution (single row of 3 subplots: u, v, w), twin y-axes.
  Left  y-axis  — y-averaged RMS error vs truth  E_comp^{+}  (solid red)
  Right y-axis  — y-averaged correlation         C_comp      (dashed dark blue)

Figure 2 — y-profile of the error (single row of 3 subplots: u, v, w).
  Several time instants: early times thin & transparent, later times thick & opaque.

Output: PDF + PNG (1200 dpi).
"""

# %%
# ~ LOAD THE FRAMEWORK
# ~ ═══════════════════════════════════════════════════════════════════════════
from framework import *

# %%
# ~ LOAD DATA (from W20)
# ~ ═══════════════════════════════════════════════════════════════════════════
data    = np.load(data_loc + 'error_wTruth.npz')
time_p  = data['time_p']
frames  = data['frames']
error_m = data['error_m']    #+ (N_frames, 4) - u, v, w, p
corr_m  = data['corr_m']     #+ (N_frames, 4)

component_names = ['u', 'v', 'w']
yp_comp      = [data[f'yp_{c}']      for c in component_names]   #+ (n_y_c,)
error_y_comp = [data[f'error_y_{c}'] for c in component_names]   #+ (N_frames, n_y_c)

N_frames = len(frames)
log.info(f"Loaded {N_frames} frames, t+ = {time_p[0]:.1f} to {time_p[-1]:.1f}")

# %%
# ~ MATPLOTLIB / LATEX CONFIGURATION  (JFM — mirrors ref_settings_plot.py)
# ~ ═══════════════════════════════════════════════════════════════════════════
mpl.rcParams.update({
    "text.usetex":              True,
    "text.latex.preamble":      r"\usepackage{mathrsfs}",
    "font.family":              "serif",
    "font.serif":               ["Computer Modern Roman"],
    "mathtext.fontset":         "cm",

    "font.size":                12,
    "axes.labelsize":           16,
    "xtick.labelsize":          12,
    "ytick.labelsize":          12,
    "legend.fontsize":          10,

    "xtick.direction":          "in",
    "ytick.direction":          "in",
    "xtick.major.size":         4,
    "xtick.minor.size":         2,
    "ytick.major.size":         4,
    "ytick.minor.size":         2,
    "xtick.major.width":        0.6,
    "xtick.minor.width":        0.4,
    "ytick.major.width":        0.6,
    "ytick.minor.width":        0.4,
    "xtick.top":                True,
    "ytick.right":              False,   #+ twinx owns the right side
    "xtick.minor.visible":      True,
    "ytick.minor.visible":      True,

    "axes.linewidth":           0.6,

    "figure.dpi":               1200,
    "savefig.dpi":              1200,
    "savefig.bbox":             "tight",
    "savefig.pad_inches":       0.02,
})

# ~ COLOUR & STYLE CONSTANTS
# ~ ═══════════════════════════════════════════════════════════════════════════
col_err  = '#B5332E'    #+ matte red  — all error lines (solid)
col_corr = '#1F3A68'    #+ dark blue  — all correlation lines (dashed)

lw_main  = 1.4          #+ temporal plot lines

#* y-profile time instants: thin/transparent (early) -> thick/opaque (late)
lw_min,  lw_max  = 0.5, 2.0
alp_min, alp_max = 0.2, 1.0

lab_size = 12.5
tck_size = 12
txt_size = 12
dpi      = 1200

#* Correlation axis limits
ylim_corr_ax = [0, 1.05]

# ~ FIGURE LAYOUT
# ~ ═══════════════════════════════════════════════════════════════════════════
fig_w = 11.0
fig_h = 3.0

panel_labels_t = [r'$(a)$', r'$(b)$', r'$(c)$']
panel_labels_y = [r'$(d)$', r'$(e)$', r'$(f)$']

def format_axis(ax, xlabel):
    """Common left-axis formatting (ticks inward, minor ticks, x-label)."""
    ax.tick_params(axis='both', which='major',
                   labelsize=tck_size, labelcolor='k',
                   direction='in', width=0.6, length=4)
    ax.tick_params(axis='both', which='minor',
                   direction='in', width=0.4, length=2)
    ax.set_xlabel(xlabel, fontsize=lab_size)

# %%
# ~ FIGURE 1 — TEMPORAL ERROR & CORRELATION
# ~ ═══════════════════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(fig_w, fig_h))
gs  = fig.add_gridspec(1, 3, wspace=0.35)
axes = [fig.add_subplot(gs[0, i]) for i in range(3)]

for i, (ax, comp, plabel) in enumerate(zip(axes, component_names, panel_labels_t)):

    ax2 = ax.twinx()

    # ─── error (left y-axis, solid red) ──────────────────────────────────
    ax.plot(time_p, error_m[:, i],
            color=col_err, lw=lw_main, ls='-', zorder=3)

    # ─── correlation (right y-axis, dashed dark blue) ────────────────────
    ax2.plot(time_p, corr_m[:, i],
             color=col_corr, lw=lw_main, ls='--', zorder=3)

    # ─── left axis formatting ────────────────────────────────────────────
    ax.set_xlim([time_p[0], time_p[-1]])
    ax.set_ylim(bottom=0)
    format_axis(ax, r'$t^+$')
    ax.tick_params(axis='y', which='major', labelcolor=col_err)

    # ─── right axis formatting ───────────────────────────────────────────
    ax2.set_ylim(ylim_corr_ax)
    ax2.tick_params(axis='y', which='major',
                    labelsize=tck_size, labelcolor=col_corr,
                    direction='in', width=0.6, length=4)
    ax2.tick_params(axis='y', which='minor',
                    direction='in', width=0.4, length=2)
    #* suppress x-ticks on twin (primary already draws them)
    ax2.tick_params(axis='x', which='both', bottom=False, top=False)

    # ─── axis titles (above each y-axis, component-specific) ─────────────
    ax.text(-0.04, 1.02,
            rf'$\mathscr{{E}}_{{{comp}}}^{{+}}$',
            fontsize=lab_size, color=col_err,
            transform=ax.transAxes, ha='center', va='bottom')
    ax2.text(1.04, 1.02,
             rf'$\mathscr{{C}}_{{{comp}}}$',
             fontsize=lab_size, color=col_corr,
             transform=ax.transAxes, ha='center', va='bottom')

    # ─── panel label ─────────────────────────────────────────────────────
    ax.text(-0.12, 1.02, plabel, transform=ax.transAxes,
            fontsize=txt_size, va='top', ha='right')

# ~ SAVE
# ~ ═══════════════════════════════════════════════════════════════════════════
if flag_output:
    out_stem = save_loc + 'fig_error_wTruth_time'
    fig.savefig(out_stem + '.pdf')
    fig.savefig(out_stem + '.png', dpi=dpi)
    log.info(f"Saved: {out_stem}.pdf / .png")

# %%
# ~ FIGURE 2 — ERROR y-PROFILES AT SEVERAL TIMES
# ~ ═══════════════════════════════════════════════════════════════════════════
    #--> Time instants (evenly spread, always including the final frame)
t_inds = np.unique(np.linspace(0, N_frames - 1, n_profiles).astype(int))
n_t    = len(t_inds)
w_t    = np.linspace(0, 1, n_t) if n_t > 1 else np.ones(1)  #^ 0 = earliest, 1 = latest

fig = plt.figure(figsize=(fig_w, fig_h))
gs  = fig.add_gridspec(1, 3, wspace=0.30)
axes = [fig.add_subplot(gs[0, i]) for i in range(3)]

for i, (ax, comp, plabel) in enumerate(zip(axes, component_names, panel_labels_y)):

    # ─── error profiles (solid red, thin/transparent -> thick/opaque) ────
    for w, it in zip(w_t, t_inds):
        ax.plot(yp_comp[i], error_y_comp[i][it, :],
                color=col_err, ls='-', zorder=3,
                lw    = lw_min  + w * (lw_max  - lw_min),
                alpha = alp_min + w * (alp_max - alp_min))

    # ─── axis formatting ─────────────────────────────────────────────────
    ax.set_xlim([yp_comp[i][0], yp_comp[i][-1]])
    ax.set_ylim(bottom=0)
    format_axis(ax, r'$y^+$')

    # ─── axis title (component-specific) ─────────────────────────────────
    ax.text(-0.04, 1.02,
            rf'$\mathscr{{E}}_{{{comp}}}^{{+}}$',
            fontsize=lab_size, color=col_err,
            transform=ax.transAxes, ha='center', va='bottom')

    # ─── panel label ─────────────────────────────────────────────────────
    ax.text(-0.12, 1.02, plabel, transform=ax.transAxes,
            fontsize=txt_size, va='top', ha='right')

#* Note the time range shown (on last panel)
axes[-1].text(0.97, 0.95,
              rf'$t^+ = {time_p[t_inds[0]]:.0f} \rightarrow {time_p[t_inds[-1]]:.0f}$',
              transform=axes[-1].transAxes, fontsize=txt_size-2,
              ha='right', va='top')

# ~ SAVE
# ~ ═══════════════════════════════════════════════════════════════════════════
if flag_output:
    out_stem = save_loc + 'fig_error_wTruth_yProfile'
    fig.savefig(out_stem + '.pdf')
    fig.savefig(out_stem + '.png', dpi=dpi)
    log.info(f"Saved: {out_stem}.pdf / .png")

# %%
