#%%[markdown]
"""
## Framework Initialization - Compute Error with Truth Module 

### This module handles common initialization tasks:
    - Import common modules and configuration
    - Load settings and simulation parameters from YAML
    - Setup base and root directories
    - Configure logging
    - Load flow parameters
    - Setup grid data (subdomain + full domain)
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
log = logging.getLogger('W20:')

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
root_module     = base_home + project_name + module_path

base_data       = base_home + project_name + '/data/'
base_data_mount = base_work + project_name + '/data/'

base_out        = base_home + project_name + '/results/'
base_out_mount  = base_work + project_name + '/results/'

root_grid       = base_data + 'DNS_grid/'

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
log = logging.getLogger('W20:')

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


# %%
#- ═════════════════════════════════════════════════════════════════════════════
#- FLOW PARAMETERS
#- ═════════════════════════════════════════════════════════════════════════════

#* Flow parameters from simulation_parameters.yaml (Re_cen, Re_tau, dt already in globals)
u_tau    = 1.0                     #^ nondimensional: half-height = 1, u_tau = 1
delta_nu = 1.0 / Re_tau            #^ nu = 1/Re_tau, so Re_tau = u_tau / nu
t_delta  = Re_cen / Re_tau**2
dt_p     = dt  / t_delta
gradx_Pi = u_tau**2 / 1.0        #^ dP/dx = -u_tau^2 / h, with half-height h = 1
GPxi     = gradx_Pi * dt         #^ Code takes pressure gradient in multiples of 'dt'
norm_pres = u_tau**2

log.info('=' * 70)
log.info(f" System details:")
log.info(f" --------------------")
log.info(f"Re_cen     = {Re_cen:.2f}")
log.info(f"Re_tau     = {Re_tau:.2f}")
log.info(f"u_tau      = {u_tau:.4e}")
log.info(f"delta_nu   = {delta_nu:.4e}")
log.info(f"dt_p       = {dt_p:.4e} ")
log.info(f"GPxi       = {GPxi:.4e} ")

# %%
#- ═════════════════════════════════════════════════════════════════════════════
#- SIMULATION-SPECIFIC SETUP
#- ═════════════════════════════════════════════════════════════════════════════
if len(sys.argv) > 1:
    inp_1 = sys.argv[1]

#*  Root of the simulation directory
root_sim          = base_data_mount + sim_series + sim_name + '/'
root_truth        = base_data_mount + sim_series + sim_name_truth + '/'
frames_DNS        = np.arange(t_start, t_end+1, t_skip)
time_p            = (frame0  + frames_DNS) * dt_p
N_frames          = len(frames_DNS)

#* Set relevant I/O directories
root_out_mount = base_out_mount + sim_series + sim_name + '/' + out_dir + '/'
root_out       = base_out + sim_series + sim_name + '/' + out_dir + '/'

if flag_output:
    create_directory(root_out_mount)
    create_directory(root_out)

log.info(f"Simulation Directory: {root_sim}")
log.info(f"DNS frames start at: {frames_DNS[0]} to {frames_DNS[-1]}, in steps of {t_skip}")
log.info(f'Output location: {root_out}')
log.info('=' * 70)

# %%
#- ═════════════════════════════════════════════════════════════════════════════
#- GRID CONFIGURATION
#- ═════════════════════════════════════════════════════════════════════════════
#* Import grid operations
from grid_operations import read_grid_data, get_grid_center_wHalos, get_staggered_grid

#* Set relevant directories
grid_loc         = root_grid + grid_name + '/'

#* Load the grid parameters file
with open(grid_loc + 'grid_parameters.yaml') as file:
    globals().update(yaml.safe_load(file))

#* Read grid.data
_, _, xx, yy, zz = read_grid_data(grid_loc + 'grid.data', ngx, ngy, ngz, Lx, Ly, Lz)
zz = zz + z0  #^ z grid always starts from 0

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

#- ═════════════════════════════════════════════════════════════════════════════
#- INITIALIZATION COMPLETE
#- ═════════════════════════════════════════════════════════════════════════════

log.info('Framework loaded successfully. Ready for processing...')
log.info('=' * 70)

# %%
