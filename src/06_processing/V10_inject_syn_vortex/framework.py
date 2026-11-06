#%%[markdown]
"""
## Framework Initialization - Inject Synthetic Vortex Module (100_BTK_CHL)

### This module handles common initialization tasks:
    - Import common modules and configuration
    - Load settings from YAML
    - Setup base and root directories
    - Configure logging
    - Load flow parameters (outer units: h = 1, U_b = 1)
    - Setup grid data (full channel, staggered with halos)
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
with open(str(Path(__file__).parent / 'settings.yaml')) as file:
    settings = yaml.safe_load(file)
    globals().update(settings)


#- ═════════════════════════════════════════════════════════════════════════════
#- LOGGING CONFIGURATION
#- ═════════════════════════════════════════════════════════════════════════════

#* Configure logging
log = logging.getLogger('V10:')

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
root_out        = base_home + project_name + '/results/'
root_sim        = base_data + sim_name + '/'
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
log = logging.getLogger('V10:')

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

#* Load the grid/flow parameters of the source package (ngx, ngy, ngz, Lx, Ly, Lz, Re_tau, Re_bulk)
with open(root_sim + 'grid_parameters.yaml') as file:
    globals().update(yaml.safe_load(file))

#* Cross-check the Reynolds number with the solver input file
Re_input = None
with open(root_sim + 'input_forward') as file:
    for line in file:
        tokens = line.split()
        if len(tokens) > 1 and tokens[1] == 'Re':
            Re_input = float(tokens[0])
            break
if Re_input is not None and not np.isclose(Re_input, Re_bulk):
    log.warning(f"Re in input_forward ({Re_input}) differs from Re_bulk ({Re_bulk}). Using input_forward.")
    Re_bulk = Re_input

#* Outer units: h = 1, U_b = 1
nu       = 1.0 / Re_bulk                     #^ Kinematic viscosity
u_tau    = Re_tau / Re_bulk                  #^ Friction velocity (= Re_tau * nu / h)
delta_nu = nu / u_tau                        #^ Viscous length scale
t_delta  = delta_nu / u_tau                  #^ Viscous time scale

log.info('=' * 70)
log.info(f" System details:")
log.info(f" --------------------")
log.info(f"Re_bulk    = {Re_bulk:.2f}")
log.info(f"Re_tau     = {Re_tau:.2f}")
log.info(f"nu         = {nu:.4e}")
log.info(f"u_tau      = {u_tau:.4e}")
log.info(f"delta_nu   = {delta_nu:.4e}")
log.info(f"t_delta    = {t_delta:.4e}")


#- ═════════════════════════════════════════════════════════════════════════════
#- SIMULATION-SPECIFIC SETUP
#- ═════════════════════════════════════════════════════════════════════════════
if len(sys.argv) > 1 and sys.argv[1] in ("lamb_oseen", "burgers", "burgers_packet", "hill", "taylor_green", "abc"):
    vortex_choice = sys.argv[1]

#* Case naming
case_name     = 'type_' + vortex_choice + case_suffix

#* Set relevant I/O directories
field_folder  = root_sim
save_loc      = base_data + save_dir + '/' + case_name + '/'
fig_loc       = save_loc + fig_dir + '/'

log.info(f"Vortex choice: {vortex_choice}")
log.info(f'Source simulation directory: {field_folder}')
log.info(f'Save location: {save_loc}')
log.info('=' * 70)


#- ═════════════════════════════════════════════════════════════════════════════
#- GRID CONFIGURATION
#- ═════════════════════════════════════════════════════════════════════════════

#* Import grid operations
from grid_operations import read_grid_data, get_grid_wHalos, get_grid_center_wHalos, get_staggered_grid

#* Read grid.data of full channel
_, _, xx, yy, zz = read_grid_data(root_sim + 'grid.data', ngx, ngy, ngz, Lx, Ly, Lz)

#* Wall locations
y_wall_bot = yy[0]
y_wall_top = yy[-1]

#* Convert to wall units
xp = xx / delta_nu
yp = (yy - y_wall_bot) / delta_nu
zp = zz / delta_nu

#* Cell faces and centers (with Halos)
xx_fc = get_grid_wHalos(xx)
yy_fc = get_grid_wHalos(yy)
zz_fc = get_grid_wHalos(zz)
xx_cc = get_grid_center_wHalos(xx)
yy_cc = get_grid_center_wHalos(yy)
zz_cc = get_grid_center_wHalos(zz)

#* Staggered grid (with Halos)
xu, yu, zu = get_staggered_grid(xx, yy, zz, 1)
xv, yv, zv = get_staggered_grid(xx, yy, zz, 2)
xw, yw, zw = get_staggered_grid(xx, yy, zz, 3)

log.info(f"Grid: ngx = {ngx}, ngy = {ngy}, ngz = {ngz}")
log.info(f"Domain: Lx = {Lx:.4f}, Ly = {Ly:.4f}, Lz = {Lz:.4f}")
log.info(f"dy+_wall = {yp[1]:.3f}, dx+ = {xp[1]-xp[0]:.3f}, dz+ = {zp[1]-zp[0]:.3f}")


#- ═════════════════════════════════════════════════════════════════════════════
#- INITIALIZATION COMPLETE
#- ═════════════════════════════════════════════════════════════════════════════

log.info('Framework loaded successfully. Ready for processing...')
log.info('=' * 70)

# %%
