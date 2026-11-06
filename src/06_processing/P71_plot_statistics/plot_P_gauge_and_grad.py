# %% [markdown]
'''
> __Project__    : DNS with a channel subdomain
> __Author__     : Sugan Velazhagan                                                                        \
> __Date__       : 08/07/2026                                                                               \
> __Description__: Plot gauge pressure G(t) and mean pressure gradient g(t) from
                    pressure_stats_diag.npz (written by debug_pressure.py).
                    Full-domain quantities plotted with thin lines;
                    subdomain quantities with thick lines.
'''
# %%
from framework import *
from plot_functions import plot_settings, plot_settings_woLatex

ts = tt.time()

# %%
# ~ SELECT DOMAIN  (mirrors flag_truth logic in debug_pressure.py)
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
flag_truth = False  #! Toggle: True → read from ground_truth folder, False → read from save_loc

if flag_truth:
    stat_folder = root_out_mount + truth_dir + '/' + save_dir + '/'
    save_loc  = root_out + truth_dir + '/' + save_dir + '/'
    create_directory(save_loc)
else:
    stat_folder = export_loc 

diag_file = stat_folder + 'pressure_stats_diag.npz'
log.info(f"Reading: {diag_file}")

data = np.load(diag_file)

gauge_full_t = data['gauge_full_t']   # (N_frames,)  full-domain gauge  [norm_pres]
gradx_full_t = data['gradx_full_t']   # (N_frames,)  full-domain mean x-gradient  [norm_pres]
gauge_sub_t  = data['gauge_sub_t']    # (N_frames,)  subdomain gauge  [norm_pres]
gradx_sub_t  = data['gradx_sub_t']    # (N_frames,)  subdomain mean x-gradient  [norm_pres]

# %%
# ~ TIME AXIS
# ~ ~~~~~~~~~~
t_plus = frames_DNS * dt / t_delta    # viscous time units t^+

# %%
# ~ PLOT PARAMETERS
# ~ ~~~~~~~~~~~~~~~~~
colB1    = '#03045e'
dpi      = 300
ls       = 32   # label font size
ts_font  = 28   # tick font size

colG = 'k'              # G(t) panel colour: black
colg = '#0077b6'        # g(t) panel colour: blue

lw_full  = 1.5   # thin  → full domain
lw_sub   = 4.5   # thick → subdomain

alp_line = 1.0
alp_sub  = 0.7

# %%
# ~~ PLOT
# ~ ~~~~~~
if username == 'smuruga4':
    plot_settings_woLatex()
    mpl.rcParams["text.usetex"] = False
else:
    plot_settings()

fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

# ── Top panel: gauge pressure G^+(t) ──────────────────────────────────────────
axes[0].plot(t_plus, gauge_full_t,
             '-', color=colG, linewidth=lw_full, alpha=alp_line,
             label=r'$G^+_{\mathrm{full}}(t)$')
axes[0].plot(t_plus, gauge_sub_t,
             '-', color=colG, linewidth=lw_sub, alpha=alp_sub,
             label=r'$G^+_{\mathrm{sub}}(t)$')

axes[0].set_ylabel(r'$G^+(t)$', fontsize=ls, rotation=0, labelpad=40)
axes[0].set_ylim(-1.5, 1.5)
# axes[0].legend(fontsize=ts_font - 4, framealpha=0.6)
axes[0].tick_params(axis='both', which='major', direction='in',
                    colors='k', labelsize=ts_font, width=1.5, length=6, pad=6,
                    bottom=True, top=True, left=True, right=True)
axes[0].tick_params(axis='both', which='minor', direction='in',
                    colors='k', labelsize=ts_font, width=1.0, length=3, pad=3,
                    bottom=True, top=True, left=True, right=True)
axes[0].grid(True, alpha=0.2, which='major')
axes[0].grid(True, alpha=0.08, which='minor')

# ── Bottom panel: mean pressure gradient g^+(t) ───────────────────────────────
axes[1].plot(t_plus, gradx_full_t,
             '-', color=colg, linewidth=lw_full, alpha=alp_line,
             label=r'$g^+_{\mathrm{full}}(t)$')
axes[1].plot(t_plus, gradx_sub_t,
             '-', color=colg, linewidth=lw_sub, alpha=alp_sub,
             label=r'$g^+_{\mathrm{sub}}(t)$')

axes[1].set_ylabel(r'$g^+(t)$', fontsize=ls, rotation=0, labelpad=40)
axes[1].set_xlabel(r'$t^+$', fontsize=ls, rotation=0, labelpad=10)
axes[1].set_ylim(-2.5, 2.5)
# axes[1].legend(fontsize=ts_font - 4, framealpha=0.6)
axes[1].tick_params(axis='both', which='major', direction='in',
                    colors='k', labelsize=ts_font, width=1.5, length=6, pad=6,
                    bottom=True, top=True, left=True, right=True)
axes[1].tick_params(axis='both', which='minor', direction='in',
                    colors='k', labelsize=ts_font, width=1.0, length=3, pad=3,
                    bottom=True, top=True, left=True, right=True)
axes[1].grid(True, alpha=0.2, which='major')
axes[1].grid(True, alpha=0.08, which='minor')

plt.tight_layout()
plt.savefig(save_loc + 'P_gauge_and_grad.png',
            bbox_inches='tight', dpi=dpi, transparent=False,
            facecolor='w', edgecolor='k')
plt.show()
log.info(f"Saved: {save_loc}P_gauge_and_grad.png")

# %%
te = tt.time()
log.info(f"Total time: {te - ts:.2f} seconds")
# %%
