#%%[markdown]
"""
## Framework Initialization - Compute Statistics Module (32_DA_CHL)

### This module handles common initialization tasks:
    - Import common modules and configuration
    - Load settings and simulation parameters from YAML
    - Setup base and root directories
    - Configure logging
    - Load flow parameters
    - Setup grid data
    - Import utility functions
"""
#%%
#- ═════════════════════════════════════════════════════════════════════════════
#- IMPORTS
#- ═════════════════════════════════════════════════════════════════════════════

import sys
from pathlib import Path

#* Add common directory to path for imports
common_dir = str(Path(__file__).parent.parent / 'common')
if common_dir not in sys.path:
    sys.path.insert(0, common_dir)

#* Import common modules (numpy, matplotlib, yaml, logging, etc.)
from modules import *

#* Import configuration I/O (base_home, base_work, username, cwd)
from config_io import *

#- ═════════════════════════════════════════════════════════════════════════════
#- CONFIGURATION LOADING
#- ═════════════════════════════════════════════════════════════════════════════

#* Load settings from YAML files
with open('./settings.yaml') as file:
    settings = yaml.safe_load(file)
    globals().update(settings)

#* Load simulation parameters
for yaml_file in [rsrc_loc + "simulation_parameters.yaml"]:
    with open(yaml_file) as f:
        globals().update(yaml.safe_load(f))


#- ═════════════════════════════════════════════════════════════════════════════
#- LOGGING CONFIGURATION
#- ═════════════════════════════════════════════════════════════════════════════

#* Configure logging
log = logging.getLogger('Y71:')

if flag_debug:
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    log.debug(f" Settings loaded:")
    log.debug(f" --------------------")
    log.debug(f" Debug mode is ON")
else:
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

#* Suppress verbose plotting library output
logging.getLogger("PIL").setLevel(logging.WARNING)
logging.getLogger("matplotlib").setLevel(logging.WARNING)


#- ═════════════════════════════════════════════════════════════════════════════
#- DIRECTORY SETUP
#- ═════════════════════════════════════════════════════════════════════════════

#* Base directories (from config_io: base_home, base_work)

#* Root directories
base_data       = base_home + project_name + '/data/'
base_data_mount = base_work + project_name + '/data/'
root_out        = base_home + project_name + '/results/'
root_out_mount  = base_work + project_name + '/results/'
root_grid       = base_data + 'DNS_grid/'
root_truth      = base_data_mount + truth_dir + '/'
root_module     = base_home + project_name + module_path

#* Add module to Python path
sys.path.append(root_module)


#- ═════════════════════════════════════════════════════════════════════════════
#- UTILITY IMPORTS
#- ═════════════════════════════════════════════════════════════════════════════

from utilities import *


#- ═════════════════════════════════════════════════════════════════════════════
#- LOGGING RECONFIGURATION
#- ═════════════════════════════════════════════════════════════════════════════

#* Reconfigure logging (if needed after utilities import)
log = logging.getLogger('Y71:')

if flag_debug:
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
else:
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

logging.getLogger('matplotlib').setLevel(logging.INFO)


#- ═════════════════════════════════════════════════════════════════════════════
#- FLOW PARAMETERS
#- ═════════════════════════════════════════════════════════════════════════════

#* Load flow parameters
params   = np.load(rsrc_loc + 'flow_params.npz')
Re_tau   = params['Re_tau']
Re_cen   = params['Re_cen']
u_tau    = params['u_tau']
delta_nu = params['delta_nu']
t_delta  = Re_cen / Re_tau**2
dt_m_p   = dt * ObsNt / t_delta
dt_m     = dt * ObsNt
GPxi     = 1.892255625E-06
gradx_Pi = GPxi / dt
norm_pres = u_tau**2

log.info('=' * 70)
log.info(f" System details:")
log.info(f" --------------------")
log.info(f"Re_cen     = {Re_cen:.2f}")
log.info(f"Re_tau     = {Re_tau:.2f}")
log.info(f"u_tau      = {u_tau:.4e}")
log.info(f"delta_nu   = {delta_nu:.4e}")
log.info(f"dt_m_p     = {dt_m_p:.4e} ")
log.info(f"dt_m       = {dt_m:.4e} ")
log.info(f"dt         = {dt:.4e} ")
log.info(f"GPxi       = {GPxi:.4e} ")

#- ═════════════════════════════════════════════════════════════════════════════
#- SIMULATION-SPECIFIC SETUP
#- ═════════════════════════════════════════════════════════════════════════════
if len(sys.argv) > 1:
    sim_suffix = sys.argv[1]

#* Simulation number and naming
sim_label         = DA_series + "{:03d}".format(sim_no) + sim_suffix  # Label for the simulation
sim_name          = sim_prefix + sim_label
frame0            = frame_ref + (sim_no - 1) * (nOBS - nOLP)  # Start frame
frames_DNS        = np.arange(t_start, t_end+1, t_skip)
frames            = frame0 + frames_DNS // ObsNt
time_p            = frame0 * dt_m_p + frames_DNS * dt_m_p / ObsNt
N_frames          = len(frames_DNS)

log.info(f"Window no: {sim_no}")
log.info(f"Simulation name: {sim_name}")
log.info(f"DNS frames start at: {frames_DNS[0]} to {frames_DNS[-1]}, in steps of {t_skip}")

#* Set relevant I/O directories
root_sim      = base_data_mount + sim_name + '/'
field_folder  = root_sim
field_folder_truth = root_truth
export_loc    = root_out_mount + sim_name + '/' + export_dir +  '/'
save_loc      = root_out + sim_name + '/' + save_dir +  '/'

if flag_output:
    create_directory(export_loc)

log.info(f'Root simulation directory: {root_sim}')
log.info(f'Export location: {export_loc}')
log.info('=' * 70)
#- ═════════════════════════════════════════════════════════════════════════════
#- GRID CONFIGURATION
#- ═════════════════════════════════════════════════════════════════════════════

#* Import grid operations
from grid_operations import read_grid_data, get_grid_center_wHalos, get_staggered_grid

#* Set relevant directories
grid_loc         = root_grid + grid_name + '/'
grid_loc_FD      = root_grid + grid_name_FD + '/'

#* Load the full channel grid parameters file
with open(grid_loc_FD + 'grid_parameters.yaml') as file:
    globals().update(yaml.safe_load(file))

#* Load the subdomain grid parameters file
with open(grid_loc + 'grid_parameters.yaml') as file:
    globals().update(yaml.safe_load(file))

#* Read grid.data of full channel
_, _, xx_FD, yy_FD, zz_FD = read_grid_data(grid_loc_FD + 'grid.data', ngx_FD, ngy_FD, ngz_FD, Lx_FD, Ly_FD, Lz_FD)
xx_FD = xx_FD - xx_FD[istart_meas] #^ Measurements start from istart_meas
zz_FD = zz_FD - zz_FD[kstart_meas] #^ Measurements start from kstart_meas

#* Read grid.data of subdomain
_, _, xx, yy, zz = read_grid_data(grid_loc + 'grid.data', ngx, ngy, ngz, Lx, Ly, Lz)
zz = zz + zStart  #^ z grid always starts from 0

#* Convert to wall units
xp = xx / delta_nu
yp = yy / delta_nu
zp = zz / delta_nu

#* Cell centers (with Halos)
xx_cc = get_grid_center_wHalos(xx)
yy_cc = get_grid_center_wHalos(yy)
zz_cc = get_grid_center_wHalos(zz)

#* Staggered grid (with Halos)
xu, yu, zu = get_staggered_grid(xx, yy, zz, 1)
xv, yv, zv = get_staggered_grid(xx, yy, zz, 2)
xw, yw, zw = get_staggered_grid(xx, yy, zz, 3)

#* Load measurement locations (cell-centers)
xp_ob = np.load(grid_loc + "xp_ob.npy")
yp_ob = np.load(grid_loc + "yp_ob.npy")
zp_ob = np.load(grid_loc + "zp_ob.npy")

#- ═════════════════════════════════════════════════════════════════════════════
#- MEAN PROFILE
#- ═════════════════════════════════════════════════════════════════════════════
#* Load the mean U profile (in wall units, on yy grid)
U_yy     = np.load(rsrc_loc + 'U_profile.npy')
V_yy     = np.zeros_like(U_yy)
W_yy     = np.zeros_like(U_yy)

#- ═════════════════════════════════════════════════════════════════════════════
#- INITIALIZATION COMPLETE
#- ═════════════════════════════════════════════════════════════════════════════

log.info('Framework loaded successfully. Ready for processing...')
log.info('=' * 70)

# %%
