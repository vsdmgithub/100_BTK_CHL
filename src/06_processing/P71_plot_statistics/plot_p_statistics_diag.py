# %% [markdown]
'''
> __Project__    : DNS with a channel subdomain
> __Author__     : Sugan Velazhagan                                                                        \
> __Date__       : 08/07/2026                                                                               \
> __Description__: Plot pressure variance from pressure_stats_diag.npz (written by debug_pressure.py).
                    Reads both DA and truth datasets and plots them together.
                    Mirrors the pressure variance panel style in plot_statistics.py.
'''
# %%
from framework import *
# from plot_functions import plot_settings, plot_settings_woLatex

ts = tt.time()

# %%
# ~ READ PRESSURE STATISTICS: DA and TRUTH
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
file_DA    = export_loc + 'pressure_stats_diag.npz'
file_TR    = root_out_mount + truth_dir + '/' + save_dir + '/pressure_stats_diag.npz'

log.info(f"Reading DA  : {file_DA}")
log.info(f"Reading Truth: {file_TR}")

data_DA = np.load(file_DA)
data_TR = np.load(file_TR)

# yuIp = data_TR['yuIp']
yuIp = yu[1:jend+1] / delta_nu

pp_fluc_DA = data_DA['pp'] - data_DA['p_mean']**2
pp_fluc_TR = data_TR['pp'] - data_TR['p_mean']**2

# %%
# ~ READ DNS REFERENCE DATA
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~
ngy_SD       = 193
mean_profile = np.loadtxt(root_out + 'statistics/avg_yp.dat')
yp_dns       = mean_profile[:ngy_SD, 0]
var_profile  = np.loadtxt(root_out + 'statistics/var_yp.dat')
pp_dns       = var_profile[:ngy_SD, 5]

# %%
# ~ PLOT PARAMETERS
# ~ ~~~~~~~~~~~~~~~~~
colB1 = '#03045e'
dpi   = 300
ls    = 36   # label font size
ts_sz = 32   # tick font size

lw_line  = 2.0
alp_line = 1.0
lw_dns   = 1.25
alp_dns  = 0.35
alp_TR   = 0.85

ms            = 6
markevery     = (0.0, 0.02)
markers       = ['o', 's', 'D', '^', 'v', '<', '>', 'p', '*', 'h', 'X', 'P', '8', 'H', 'd']
dashdot_loose = (0, (4, 8, 2, 8))

# %%
# ~~ PLOT: pressure variance
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~
def plot_settings():
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "sans-serif",  # or 'sans-serif'
        "font.serif": ["Computer Modern Roman"],  # default LaTeX font
        "text.latex.preamble": r"\usepackage{newtxmath}"  # enable cmbright
    })
# if username == 'smuruga4':
#     plot_settings_woLatex()
#     mpl.rcParams["text.usetex"] = False
# else:
#     plot_settings()
plot_settings()

fig = plt.figure(figsize=(8, 6))

# DA: solid line
plt.plot(yuIp, pp_fluc_DA,
         linestyle='-', color=colB1, linewidth=lw_line, alpha=alp_line, zorder=1)

# Truth: thin line with markers
plt.plot(yuIp, pp_fluc_TR,
         '-', color=colB1,
         markerfacecolor=colB1, markeredgecolor=colB1,
         linewidth=0.3*lw_line, marker=markers[9],
         markersize=ms, markevery=markevery,
         markeredgewidth=0.0, alpha=alp_TR, zorder=2)

# DNS reference: dash-dot
plt.plot(yp_dns, pp_dns,
         linestyle=dashdot_loose, color='k', linewidth=lw_dns, alpha=alp_dns)

plt.ylabel(r'$pp^+$', fontsize=ls, rotation=90, labelpad=26)
plt.xlabel(r'$y^+$', rotation=0, labelpad=10, fontsize=ls)

plt.xlim(0.40, 6E2)
plt.ylim(0, 9)
plt.xscale('log')

plt.tick_params(axis='both', which='major', direction='in', colors='k', labelsize=ts_sz,
                width=1.5, length=6, labelcolor='k', pad=6,
                bottom=True, top=True, left=True, right=True)
plt.tick_params(axis='both', which='minor', direction='in', colors='k', labelsize=ts_sz,
                width=1.0, length=3, labelcolor='k', pad=3,
                bottom=True, top=True, left=True, right=True)

plt.gca().yaxis.set_major_locator(ticker.MultipleLocator(3.0))
plt.gca().yaxis.set_minor_locator(ticker.MultipleLocator(0.6))
plt.gca().xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=10))
plt.gca().xaxis.set_minor_locator(LogLocator(base=10.0, subs=range(2, 10), numticks=10))
plt.grid(True, alpha=0.2, which='major')
plt.grid(True, alpha=0.08, which='minor')

plt.tight_layout()
plt.savefig(save_loc + 'p_var_diag.png',
            bbox_inches='tight', dpi=dpi, transparent=False,
            facecolor='w', edgecolor='k', orientation='landscape')
plt.show()
log.info(f"Saved: {save_loc}p_var_diag.png")

# %%
te = tt.time()
log.info(f"Total time: {te - ts:.2f} seconds")
# %%
