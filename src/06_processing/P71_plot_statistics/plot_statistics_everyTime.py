# %% [markdown]
'''
> __Project__    : Data assimilation with experimental measurements of turbulent boundary layer \
> __Author__     : Sugan Velazhagan                                                              \
> __Date__       : 06/17/2025                                                                     \
> __Description__: This code is for plotting the velocity statistics at every time
  
'''
# %% 
# ~  LOAD THE FRAMEWORK
#~ ~~~~~~~~~~~~~~~~~~~~~
from framework import * 
ts=tt.time()

# ~ LOGISTICS
export_dir = save_loc + 'everyTime_stats/'
create_directory(export_dir)

#* Read rescaling factor
rSc = np.load('./rSc.npy')
yy_meas = np.loadtxt(base_data + 'yp_exp.txt')
yy_meas = yy_meas*rSc

#* Load Re_tau=2000 DNS statistics.
DNS_2000 = np.loadtxt('../../DNS_statistics_Retau_2000.dat', skiprows = 75)
#----------------------------------------------------------------------
#-   y/delta |  y^+ | u'u' |  v'v' | w'w' | u'v' | u'w' |  v'w' |  k  |         
#----------------------------------------------------------------------
yp_DNS = DNS_2000[:,1]
uu_DNS = DNS_2000[:,2]
vv_DNS = DNS_2000[:,3]
ww_DNS = DNS_2000[:,4]
uv_DNS = DNS_2000[:,5]
KE_DNS = DNS_2000[:,8]
DNS_vor_2000 = np.loadtxt('../../DNS_vorticity_statistics_Retau_2000.dat', skiprows = 75)
pp_DNS = DNS_vor_2000[:,8]
DNS_TKE_2000 = np.loadtxt('../../DNS_TKE_Retau_2000.dat', skiprows = 75)
TKE_P_DNS = DNS_TKE_2000[:,2]
DNS_mean_2000 = np.loadtxt('../../DNS_mean_Retau_2000.dat', skiprows = 72)
U_DNS = DNS_mean_2000[:,2]

# %% [markdown]
#* Functions to get the statistics from DA and Exp data
from expData_operations import get_filename,get_measurement_velocities_and_pressure
def get_measurement_statistics(frame):
    file_addr = root_expData + get_filename(frame)
    u_exp, v_exp, w_exp, P_exp = get_measurement_velocities_and_pressure(file_addr,nex,ney,nez)
    P_exp      = P_exp - np.mean(P_exp)
    u_mean_exp = np.mean(u_exp,axis=(0,2))
    v_mean_exp = np.mean(v_exp,axis=(0,2))
    w_mean_exp = np.mean(w_exp,axis=(0,2))
    pp_exp     = np.mean(P_exp**2,axis=(0,2)) 
    uu_fluc    = -u_mean_exp**2 + np.mean(u_exp**2,axis=(0,2))
    vv_fluc    = -v_mean_exp**2 + np.mean(v_exp**2,axis=(0,2))
    ww_fluc    = -w_mean_exp**2 + np.mean(w_exp**2,axis=(0,2))
    uv_fluc    = np.mean(u_exp*v_exp,axis=(0,2)) - u_mean_exp*v_mean_exp

    uu_fluc = uu_fluc/rSc**2
    vv_fluc = vv_fluc/rSc**2
    ww_fluc = ww_fluc/rSc**2
    uv_fluc = uv_fluc/rSc**2

    return (uu_fluc, vv_fluc, ww_fluc, uv_fluc)
    # return u_mean_exp

def load_frame_stats(npz_file):
    data = np.load(npz_file, allow_pickle=True)
    stats_dict = {key: data[key] for key in data.files}
    data.close()
    
    return stats_dict

def get_individual_arrays(stats_data, ind):
    u_mean_AF = (stats_data['u_mean'])[ind,:]  # Shape: (n_frames, n_y_points)
    v_mean_AF = (stats_data['v_mean'])[ind,:-1]  
    w_mean_AF = (stats_data['w_mean'])[ind,:]  
    uu_AF = (stats_data['uu'])[ind,:]  
    vv_AF = (stats_data['vv'])[ind,:-1]  
    ww_AF = (stats_data['ww'])[ind,:]  
    uv_AF = (stats_data['uv'])[ind,:]  
    pp_AF = (stats_data['pp'])[ind,:]  
    uu_fluc = uu_AF - u_mean_AF**2
    vv_fluc = vv_AF - v_mean_AF**2
    ww_fluc = ww_AF - w_mean_AF**2
    uv_fluc = uv_AF[:-1]-0.5*(u_mean_AF[1:]+u_mean_AF[:-1])*v_mean_AF[:]

    uu_fluc = uu_fluc/rSc**2
    vv_fluc = vv_fluc/rSc**2
    ww_fluc = ww_fluc/rSc**2
    uv_fluc = uv_fluc/rSc**2

    return uu_fluc, vv_fluc, ww_fluc, uv_fluc
    # return u_mean_AF

file_addr = root_out_mount + sim_name + '/' + save_dir + '/stats_everyFrame.npz'
stats_data = load_frame_stats(file_addr)
# sim_name_IG = 'FC_IG_L001'
file_addr = root_out_mount + sim_name_IG + '/' + save_dir + '/stats_everyFrame.npz'
stats_data_IG = load_frame_stats(file_addr)

yuIp = stats_data['yuIp']
yvIp = stats_data['yvIp']
yuIp = yuIp * rSc
yvIp = yvIp * rSc
# %% [markdown]
# ~ PLOT PARAMETERS
# ~ ~~~~~~~~~~~~~~~~~~~~~
    # --> Colors
colB1 = '#03045e'
colB2 = '#0077b6'
colB3 = '#0096c7'
colO = '#d95f02'  # Orange

dpi=150

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

ms = 8  # Marker size
alp_mark = 1.0
mew_exp  = 1.2 # Experimental marker edge width
freq_meas = (0.0, 0.02) # Experimental marker frequency

# Create a custom dash-dot linestyle with loose spacing
# Dash-dot: (dash, space, dot, space)
dashdot_loose = (0, (4, 8, 2, 8))  # Increase numbers for looser spacing
# Example: (8, 8, 2, 8) for even looser dash-dot
dash_tight = (0, (5, 2))  # Tight dash pattern

# %%
# ~~ PLOTTING
# ~~~~~~~~~~~~~~~~~~~~~
# ^ |> CORRELATIONS
def plot_correlations(frame_exp, ind):
    plt.close()
    plot_settings()

    # === Figure & Subplots ===
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), sharex=False)
    axes = axes.flatten()   # [0]=uu, [1]=vv, [2]=ww, [3]=uv

    uu_fluc_exp, vv_fluc_exp, ww_fluc_exp, uv_fluc_exp = get_measurement_statistics(frame_exp)
    uu_fluc, vv_fluc, ww_fluc, uv_fluc = get_individual_arrays(stats_data, ind)
    uu_fluc_IG, vv_fluc_IG, ww_fluc_IG, uv_fluc_IG = get_individual_arrays(stats_data_IG, ind)

    # === Helper function for consistent plotting ===
    # def plot_with_styles(ax, y_exp, R_exp, y_DA, R_DA, y_DA2, R_DA2, y_IG, R_IG, y_dns, R_dns,
    def plot_with_styles(ax, y_exp, R_exp, y_DA, R_DA, y_IG, R_IG, y_dns, R_dns,
                        color, marker, label):
        
        mew_exp = 0.0; ms=6
        freq_meas = (0.0, 0.030)
        # Experimental (rescaled)
        ax.plot(y_exp,R_exp,
                '-', color=color,
                markerfacecolor=color, markeredgecolor=color,
                linewidth=0.0, marker=marker,
                markersize=ms, markevery=freq_meas,
                markeredgewidth=mew_exp, alpha=alp_mark
        )

        # Assimilated
        ax.plot(y_DA,R_DA,
                '-', color=color, linewidth=lw_line, alpha=alp_line
                ,label=f'Model {label}'
                )

        # # Assimilated
        # ax.plot(y_DA2,R_DA2,
        #         '-', color=colB2, linewidth=0.5*lw_line, alpha=alp_line
        #         # ,label=f'Model {label}'
        #         )

        # Initial Guess
        # ax.plot(y_IG,R_IG,
        #         linestyle=dash_tight, color=color, linewidth=lw_line, alpha=alp_line
        #         # ,label=f'IG {label}'
        #         )

        # DNS reference
        ax.plot(y_dns, R_dns,
                linestyle=dashdot_loose, color='k', linewidth=lw_dns, alpha=alp_dns
                #, label=f'DNS {label}'
                )

        ax.set_xscale('log')
        ax.set_xlabel(r'$y^+$', fontsize=ls)
        ax.set_ylabel(label, fontsize=ls)

        ax.set_xlim(0.20, 1E3)

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

    # === uu subplot ===
    plot_with_styles(axes[0],
                    yy_meas, uu_fluc_exp,
                    yuIp, uu_fluc,
                    yuIp, uu_fluc_IG,
                    yp_DNS, uu_DNS,
                    colB1, markers[1], r'$uu^+$')

    # === vv subplot ===
    plot_with_styles(axes[1],
                    yy_meas, vv_fluc_exp,
                    yvIp, vv_fluc,
                    yvIp, vv_fluc_IG,
                    yp_DNS, vv_DNS,
                    colB1, markers[2], r'$vv^+$')

    # === ww subplot ===
    plot_with_styles(axes[2],
                    yy_meas, ww_fluc_exp,
                    yuIp, ww_fluc,
                    yuIp, ww_fluc_IG,
                    yp_DNS, ww_DNS,
                    colB1, markers[5], r'$ww^+$')

    # === uv subplot ===
    plot_with_styles(axes[3],
                    yy_meas, -uv_fluc_exp,
                    yvIp, -uv_fluc,
                    yvIp, -uv_fluc_IG,
                    yp_DNS, -uv_DNS,
                    colB1, markers[7], r'$-uv^+$')

    ax=axes[0]
    # ax.set_ylim(-0.5, 30)
    ax.set_ylim(-0.5, 10)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(2.0))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.4))

    ax=axes[1]
    # ax.set_ylim(-0.1, 7)
    ax.set_ylim(-0.1, 2)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.10))

    ax=axes[2]
    # ax.set_ylim(-0.1, 12)
    ax.set_ylim(-0.1, 3.2)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.2))

    ax=axes[3]
    # ax.set_ylim(-0.1, 7)
    ax.set_ylim(-0.1, 1.4)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))

    time_p = ind * dt_m_p
    plt.suptitle(r'$t^+$ = ' + f'{time_p:.2f}', fontsize=ts+2)

    plt.tight_layout(rect=[0,0,1,1])

    plt.savefig(export_dir + f"correlations_subplots_test_{frame_exp:05d}.png",
    # plt.savefig(export_dir + f"correlations_subplots_{frame_exp:05d}.png",
    # plt.savefig(export_dir + f"correlations_subplots_freeYlim_{frame_exp:05d}.png",
                    bbox_inches='tight', dpi=dpi, transparent=False,
                    facecolor='w', edgecolor='k', orientation='landscape')
    # plt.show()

frames_exp = frame0 + (sim_no-1) * (nOBS-nOLP) +  np.arange(nOLP,nOBS+1)
frames = np.arange(t_start,t_end+1,t_skip)
# !Testing
plot_correlations(frames_exp[0], 0)

# %%
#~ Statistics at every time
if frames.shape[0] == frames_exp.shape[0]:
    Parallel(n_jobs=num_workers)(
        delayed(plot_correlations)(frame_exp, ind)
        for ind, frame_exp in tqdm(list(enumerate(frames_exp)), total=len(frames_exp), desc='Writing output')
    )
# %%
print('done')

# %%
