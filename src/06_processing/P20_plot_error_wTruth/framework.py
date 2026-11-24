#%%[markdown]
"""
## Framework Initialization - Plot Error with Truth Module

### This module handles common initialization tasks:
    - Import common modules and configuration
    - Load settings and simulation parameters from YAML
    - Setup base and root directories
    - Configure logging
    - Load flow parameters
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
log = logging.getLogger('P20:')

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
log = logging.getLogger('P20:')

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

logging.getLogger('matplotlib').setLevel(logging.ERROR)


# %%
#- ═════════════════════════════════════════════════════════════════════════════
#- FLOW PARAMETERS
#- ═════════════════════════════════════════════════════════════════════════════

#* Flow parameters from simulation_parameters.yaml (Re_cen, Re_tau, dt already in globals)
u_tau    = (Re_tau/Re_cen)          #^ nondimensional: half-height = 1, u_tau = 1
delta_nu = 1.0 / Re_tau            #^ nu = 1/Re_tau, so Re_tau = u_tau / nu
t_delta  = Re_cen / Re_tau**2
dt_p     = dt  / t_delta
norm_pres = u_tau**2

log.info('=' * 70)
log.info(f" System details:")
log.info(f" --------------------")
log.info(f"Re_cen     = {Re_cen:.2f}")
log.info(f"Re_tau     = {Re_tau:.2f}")
log.info(f"u_tau      = {u_tau:.4e}")
log.info(f"delta_nu   = {delta_nu:.4e}")
log.info(f"dt_p       = {dt_p:.4e} ")

# %%
#- ═════════════════════════════════════════════════════════════════════════════
#- SIMULATION-SPECIFIC SETUP
#- ═════════════════════════════════════════════════════════════════════════════

#* Set relevant I/O directories
data_loc = base_out_mount + sim_series + sim_name + '/' + out_dir  + '/'  #^ Error data (from W20)
save_loc = base_out + sim_series + sim_name + '/' + plot_dir + '/'  #^ Plots

if flag_output:
    create_directory(save_loc)

log.info(f'Data location: {data_loc}')
log.info(f'Save location: {save_loc}')
log.info('=' * 70)


#- ═════════════════════════════════════════════════════════════════════════════
#- INITIALIZATION COMPLETE
#- ═════════════════════════════════════════════════════════════════════════════

log.info('Framework loaded successfully. Ready for plotting...')
log.info('=' * 70)

# %%
