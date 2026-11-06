# %% [markdown]
'''
> __Project__    : DNS with a channel subdomain
> __Author__     : Sugan Velazhagan                                                                        \
> __Date__       : 01/07/2025                                                                               \
> __Description__: This code is for plotting the velocity and pressure statistics
  
'''
# %% 
# ~  LOAD THE FRAMEWORK
#~ ~~~~~~~~~~~~~~~~~~~~~
from framework import * 
ts=tt.time()

# ~ READ GLOBAL STATISTICS
# ~ ~~~~~~~~~~~~~~~~~~~~~
ngy_SD = 193
mean_profile = np.loadtxt(root_out + 'statistics/avg_yp.dat')
u_mean_dns   = mean_profile[:ngy_SD,1]
yp_dns       = mean_profile[:ngy_SD,0]
var_profile =  np.loadtxt(root_out + 'statistics/var_yp.dat')
uu_dns = var_profile[:ngy_SD,1]
vv_dns = var_profile[:ngy_SD,2]
ww_dns = var_profile[:ngy_SD,3]
uv_dns = var_profile[:ngy_SD,4]
pp_dns = var_profile[:ngy_SD,5]

# ~ FUNCTION TO READ STATISTICS
def load_stats(npz_file):
    data = np.load(npz_file, allow_pickle=True)
    stats_dict = {key: data[key] for key in data.files}
    data.close()
    
    return stats_dict

def get_correlations(stats_data):
    u_mean_AF = (stats_data['u_mean'])[:]  # Shape: (n_frames, n_y_points)
    v_mean_AF = (stats_data['v_mean'])[:]  
    w_mean_AF = (stats_data['w_mean'])[:]  
    uu_AF = (stats_data['uu'])[:]  
    vv_AF = (stats_data['vv'])[:]  
    ww_AF = (stats_data['ww'])[:]  
    uv_AF = (stats_data['uv'])[:]  
    pp_AF = (stats_data['pp'])[:]  
    uu_fluc = uu_AF - u_mean_AF**2
    vv_fluc = vv_AF - v_mean_AF**2
    ww_fluc = ww_AF - w_mean_AF**2
    uv_fluc = uv_AF[:-1]-0.5*(u_mean_AF[1:]+u_mean_AF[:-1])*v_mean_AF[:]

    return uu_fluc, vv_fluc, ww_fluc, uv_fluc

   # --> Estimated
file_addr = export_folder +  'stats_ensemble.npz'
stats_data = load_stats(file_addr)
    # --> Initial guess
if flag_plotIG:
    file_addr = root_out_mount + sim_name_IG + '/' + save_dir + '/stats_ensemble.npz'
    stats_data_IG = load_stats(file_addr)
    # --> Truth statistics
file_addr = root_out_mount + truth_dir +'/' +  save_dir + '/stats_ensemble.npz'
stats_data_truth = load_stats(file_addr)    

# for key in stats_data.keys():
#     print(key, stats_data[key].shape)

# %% [markdown]
# ~ PLOT PARAMETERS
# ~ ~~~~~~~~~~~~~~~~~~~~~
    # --> Colors
colB1 = '#03045e'
# colB2 = '#023e8a'
colB2 = '#0077b6'
colB3 = '#0096c7'
colO = '#d95f02'  # Orange
color = colB1

dpi=300

    # --> Define 15 distinct markers
markers = ['o', 's', 'D', '^', 'v', '<', '>', 'p', '*', 'h', 'X', 'P', '8', 'H', 'd']

    # -->  Define 15 shades of blue and scramble them
colB_arr = plt.cm.Blues(np.linspace(0.3, 0.6, 15))
np.random.shuffle(colB_arr)

def plot_settings():
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "sans-serif",  # or 'sans-serif'
        "font.serif": ["Computer Modern Roman"],  # default LaTeX font
        "text.latex.preamble": r"\usepackage{newtxmath}"  # enable cmbright
    })

# Set very large font sizes
ls = 36  # Label size
ts = 32  # Tick size

lw_line  = 2.0
alp_line = 1.0
lw_dns   = 1.25
alp_dns  = 0.35

ms = 6  # Marker size
alp_mark = 1.0
alp_TR = 0.85
mew_exp  = 1.2 # Experimental marker edge width

# Create a custom dash-dot linestyle with loose spacing
# Dash-dot: (dash, space, dot, space)
dashdot_loose = (0, (4, 8, 2, 8))  # Increase numbers for looser spacing
dash_tight = (0, (5, 2))  # Tight dash pattern
# dash_tight = '-'
# Example: (8, 8, 2, 8) for even looser dash-dot

# ~~ PLOTTING
# ~~~~~~~~~~~~~~~~~~~~~
# ^ |> AVERAGE STREAMWISE VELOCITY
plot_settings()
fig = plt.figure(figsize=(8,6))

   # --> Declare data
yuIp = stats_data['yuIp']
u_mean_DA = stats_data['u_mean']
if flag_plotIG:
    u_mean_IG = stats_data_IG['u_mean']
u_mean_TR = stats_data_truth['u_mean']

plt.plot(yuIp,u_mean_DA, color = colB1, linewidth=lw_line,zorder=1)
if flag_plotIG:
    plt.plot(yuIp,u_mean_IG, linestyle=dash_tight, color = colB1, linewidth=lw_line)
plt.plot(yuIp,u_mean_TR,
                '-', color=color,
                markerfacecolor=color, markeredgecolor=color,
                linewidth=0.3*lw_line, marker='o',
                markersize=ms, markevery=0.03,
                markeredgewidth=mew_exp, alpha=alp_TR, zorder=2
        )

#-- Plot DNS data
lw2=1.25; alp2=0.75
plt.plot(yp_dns, u_mean_dns,
         linestyle=dashdot_loose, color = 'k', linewidth=lw2,alpha=alp_dns) 

plt.ylabel(r'$U^+$',fontsize=ls,rotation=0,labelpad=26)
plt.xlabel(r'$y^+$',rotation=0,labelpad=10,fontsize=ls)

plt.ylim(0,25)
plt.xlim(0.40,6E2)

plt.xscale('log')

plt.tick_params(axis='both',which='major',direction='in',colors='k',labelsize=ts,width=1.5,length=6,labelcolor='k',pad=6,\
                bottom=True,top=True,left=True,right=True)
plt.tick_params(axis='both',which='minor',direction='in',colors='k',labelsize=ts,width=1.0,length=3,labelcolor='k',pad=3,\
                bottom=True,top=True,left=True,right=True)

plt.gca().yaxis.set_major_locator(ticker.MultipleLocator(5.0))
plt.gca().yaxis.set_minor_locator(ticker.MultipleLocator(1.0))
plt.gca().xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=10))
# plt.gca().xaxis.set_minor_locator(ticker.LogLocator(base=10.0, numticks=10))
plt.gca().xaxis.set_minor_locator(LogLocator(base=10.0, subs=range(2, 10), numticks=10))
plt.grid(True,alpha=0.2,which='major')
plt.grid(True,alpha=0.08,which='minor')

plt.tight_layout()
plt.savefig( save_loc+ 'mean_velocity.png', bbox_inches='tight',dpi=dpi,transparent=False,\
            facecolor='w',edgecolor='k',orientation='landscape')
plt.show()

# %%
# ^ |> CORRELATIONS
plt.close()
plot_settings()

# === Figure & Subplots ===
fig, axes = plt.subplots(2, 2, figsize=(12, 10), sharex=False)
axes = axes.flatten()   # [0]=uu, [1]=vv, [2]=ww, [3]=uv

# === Helper function for consistent plotting ===
def plot_with_styles(ax, y_DA, R_DA, y_IG, R_IG, y_TR, R_TR, y_dns, R_dns,
                     color, marker, label):
    mew_exp = 0.0; ms=6
    freq_meas = (0.0, 0.030)
    # Assimilated
    ax.plot(y_DA,R_DA,
            linestyle='-', color=color, linewidth=lw_line, alpha=alp_line
            # ,label=f'Model {label}'
            )
    # Initial Guess
    if flag_plotIG:
        ax.plot(y_IG,R_IG,
                linestyle=dash_tight, color=color, linewidth=lw_line, alpha=alp_line
                # ,label=f'IG {label}'
                )
    # DNS reference
    ax.plot(y_dns, R_dns,
            linestyle=dashdot_loose, color='k', linewidth=lw_dns, alpha=alp_dns
            #, label=f'DNS {label}'
            )
    # Truth
    ax.plot(y_TR,R_TR,
            '-', color=color,
            markerfacecolor=color, markeredgecolor=color,
            linewidth=0.3*lw_line, marker=marker,
            markersize=ms, markevery=freq_meas,
            markeredgewidth=mew_exp, alpha=alp_mark
    )

    ax.set_xscale('log')
    ax.set_xlabel(r'$y^+$', fontsize=ls)
    ax.set_ylabel(label, fontsize=ls)

    ax.set_xlim(0.40,6E2)

    ax.tick_params(axis='both', which='major', direction='in',
                   colors='k', labelsize=ts, width=1.0, length=4, pad=6)
    ax.tick_params(axis='both', which='minor', direction='in',
                   colors='k', labelsize=ts, width=0.5, length=2, pad=3)
    ax.tick_params(which='both', top=True, right=True)

    # Major ticks at 10^0, 10^1, 10^2, ...
    # ax.xaxis.set_major_locator(LogLocator(base=10.0, subs=None, numticks=10))
    # Minor ticks at 2,3,...9 × 10^n
    ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=range(2, 10), numticks=10))

    ax.grid(True, alpha=0.3, which='major')
    ax.grid(True, alpha=0.1, which='minor')

   # --> Get data 
yuIp = stats_data_truth['yuIp']
yvIp = stats_data_truth['yvIp']
uu_fluc_TR, vv_fluc_TR, ww_fluc_TR, uv_fluc_TR = get_correlations(stats_data_truth)
uu_fluc_DA, vv_fluc_DA, ww_fluc_DA, uv_fluc_DA = get_correlations(stats_data)
if flag_plotIG:
    uu_fluc_IG, vv_fluc_IG, ww_fluc_IG, uv_fluc_IG = get_correlations(stats_data_IG)
# === uu subplot ===
plot_with_styles(axes[0],
                 yuIp, uu_fluc_DA,
                 yuIp, uu_fluc_TR,
                 yuIp, uu_fluc_TR,
                 yp_dns, uu_dns,
                 colB1, markers[1], r'$uu^+$')

# === vv subplot ===
plot_with_styles(axes[1],
                 yvIp, vv_fluc_DA,
                 yvIp, vv_fluc_TR,
                 yvIp, vv_fluc_TR,
                 yp_dns, vv_dns,
                 colB1, markers[2], r'$vv^+$')

# === ww subplot ===
plot_with_styles(axes[2],
                 yuIp, ww_fluc_DA,
                 yuIp, ww_fluc_TR,
                 yuIp, ww_fluc_TR,
                 yp_dns, ww_dns,
                 colB1, markers[5], r'$ww^+$')

# === uv subplot ===
plot_with_styles(axes[3],
                 yvIp, -uv_fluc_DA,
                 yvIp, -uv_fluc_TR,
                 yvIp, -uv_fluc_TR,
                 yp_dns, -uv_dns,
                 colB1, markers[7], r'$-uv^+$')

# === Limits & Layout ===
# for ax in axes:
#     ax.set_xlim(0.20, 1E3)
  
ax=axes[0]
ax.set_ylim(-0.5, 10)
ax.yaxis.set_major_locator(ticker.MultipleLocator(2.0))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.4))

ax=axes[1]
ax.set_ylim(-0.1, 2)
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.10))

ax=axes[2]
ax.set_ylim(-0.1, 3.2)
ax.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.2))

ax=axes[3]
ax.set_ylim(-0.1, 1.4)
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))

plt.tight_layout(rect=[0,0,1,1])

plt.savefig(save_loc + 'correlations.png',
            bbox_inches='tight', dpi=dpi, transparent=False,
            facecolor='w', edgecolor='k', orientation='landscape')
plt.show()

# %%
<<<<<<< HEAD
# # ^ |> Pressure
# plt.close('all')
# fig = plt.figure(figsize=(8,6))

# markevery = (0.0, 0.02)    
#    # --> Get data
# pp_DA = stats_data['pp']
# pp_IG = stats_data_IG['pp']
# pp_TR = stats_data_truth['pp']

# plt.plot(yuIp, pp_TR,
#           '-', color=colB1,
#             markerfacecolor=colB1,markeredgecolor=colB1,
#             linewidth=0.00, marker=markers[9],
#               markersize=ms, markevery=markevery,
#                 markeredgewidth=0.0, alpha=alp_mark)

# plt.plot(yuIp,pp_DA, linestyle='-',color = colB1, linewidth=lw_line,alpha=alp_line)
# plt.plot(yuIp,pp_IG, linestyle=dash_tight,color = colB1, linewidth=lw_line,alpha=alp_line)

# #-- Plot DNS data
# plt.plot(yp_dns, pp_dns, linestyle=dashdot_loose, color = 'k', linewidth=lw_dns,alpha=alp_dns)

# # var = r'${p}^+_{\sf{rms}}$'
# var = r'${pp}^+$'
# plt.ylabel(var,fontsize=36,rotation=90,labelpad=16)
# plt.xlabel(r'$y^+$',rotation=0,labelpad=5,fontsize=ls)
# plt.gca().yaxis.set_label_coords(-0.10, 0.5)

# plt.xlim(0.40,5E2)
# # plt.ylim(0,14)

# plt.xscale('log')

# plt.tick_params(axis='both',which='major',direction='in',colors='k',labelsize=ts,width=1.5,length=6,labelcolor='k',pad=6,\
#                 bottom=True,top=True,left=True,right=True)
# plt.tick_params(axis='both',which='minor',direction='in',colors='k',labelsize=ts,width=1.0,length=3,labelcolor='k',pad=3,\
#                 bottom=True,top=True,left=True,right=True)

# plt.gca().yaxis.set_major_locator(ticker.MultipleLocator(3.0))
# plt.gca().yaxis.set_minor_locator(ticker.MultipleLocator(0.6))
# plt.gca().xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=10))
# # plt.gca().xaxis.set_minor_locator(ticker.LogLocator(base=10.0, numticks=10))
# plt.gca().xaxis.set_minor_locator(LogLocator(base=10.0, subs=range(2, 10), numticks=10))
# plt.grid(True,alpha=0.2,which='major')
# plt.grid(True,alpha=0.08,which='minor')

# plt.tight_layout()
# plt.savefig( save_loc + 'p_var.png', bbox_inches='tight',dpi=dpi,transparent=False,\
#             facecolor='w',edgecolor='k',orientation='landscape')
# plt.show()
=======
# ^ |> Pressure variance
plt.close()
fig = plt.figure(figsize=(8,6))

markevery = (0.0, 0.02)
   # --> Get data
yuIp         = stats_data_truth['yuIp']
pp_fluc_DA   = stats_data['pp']        - stats_data['p_mean']**2
if flag_plotIG:
    pp_fluc_IG = stats_data_IG['pp']   - stats_data_IG['p_mean']**2
pp_fluc_TR   = stats_data_truth['pp']  - stats_data_truth['p_mean']**2

plt.plot(yuIp, pp_fluc_DA, linestyle='-', color=colB1, linewidth=lw_line, alpha=alp_line, zorder=1)
if flag_plotIG:
    plt.plot(yuIp, pp_fluc_IG, linestyle=dash_tight, color=colB1, linewidth=lw_line, alpha=alp_line)
plt.plot(yuIp, pp_fluc_TR,
         '-', color=colB1,
         markerfacecolor=colB1, markeredgecolor=colB1,
         linewidth=0.3*lw_line, marker=markers[9],
         markersize=ms, markevery=markevery,
         markeredgewidth=0.0, alpha=alp_TR, zorder=2)

#-- Plot DNS data
plt.plot(yp_dns, pp_dns, linestyle=dashdot_loose, color='k', linewidth=lw_dns, alpha=alp_dns)

plt.ylabel(r'$pp^+$', fontsize=ls, rotation=90, labelpad=26)
plt.xlabel(r'$y^+$', rotation=0, labelpad=10, fontsize=ls)

plt.xlim(0.40, 6E2)
plt.ylim(0, 9)

plt.xscale('log')

plt.tick_params(axis='both', which='major', direction='in', colors='k', labelsize=ts,
                width=1.5, length=6, labelcolor='k', pad=6,
                bottom=True, top=True, left=True, right=True)
plt.tick_params(axis='both', which='minor', direction='in', colors='k', labelsize=ts,
                width=1.0, length=3, labelcolor='k', pad=3,
                bottom=True, top=True, left=True, right=True)

plt.gca().yaxis.set_major_locator(ticker.MultipleLocator(3.0))
plt.gca().yaxis.set_minor_locator(ticker.MultipleLocator(0.6))
plt.gca().xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=10))
plt.gca().xaxis.set_minor_locator(LogLocator(base=10.0, subs=range(2, 10), numticks=10))
plt.grid(True, alpha=0.2, which='major')
plt.grid(True, alpha=0.08, which='minor')

plt.tight_layout()
plt.savefig(save_loc + 'p_var.png', bbox_inches='tight', dpi=dpi, transparent=False,
            facecolor='w', edgecolor='k', orientation='landscape')
plt.show()
>>>>>>> c818371 (Commit from Local)

# %%
