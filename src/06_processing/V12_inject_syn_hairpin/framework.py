#%%[markdown]
"""
## Framework Initialization - Inject Synthetic Hairpin Module (Burgers core along a hairpin centreline)

### This module handles common initialization tasks:
    - Import common modules and configuration
    - Load settings and simulation parameters from YAML
    - Setup base and root directories
    - Configure logging
    - Load flow parameters
    - Setup grid data (full channel, staggered with halos)
    - Derive the Burgers core and hairpin centreline parameters
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
log = logging.getLogger('V12:')

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
log = logging.getLogger('V12:')

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
#^ Outer units: half-height h = 1, bulk velocity U_b = 1
nu       = 1.0 / Re_cen            #^ Kinematic viscosity
u_tau    = (Re_tau/Re_cen)         #^ Friction velocity in outer units
delta_nu = 1.0 / Re_tau            #^ Viscous length scale (= nu / u_tau)
t_delta  = Re_cen / Re_tau**2      #^ Viscous time scale  (= delta_nu / u_tau)
dt_p     = dt  / t_delta

log.info('=' * 70)
log.info(f" System details:")
log.info(f" --------------------")
log.info(f"Re_cen     = {Re_cen:.2f}")
log.info(f"Re_tau     = {Re_tau:.2f}")
log.info(f"nu         = {nu:.4e}")
log.info(f"u_tau      = {u_tau:.4e}")
log.info(f"delta_nu   = {delta_nu:.4e}")
log.info(f"dt_p       = {dt_p:.4e} ")

# %%
#- ═════════════════════════════════════════════════════════════════════════════
#- SIMULATION-SPECIFIC SETUP
#- ═════════════════════════════════════════════════════════════════════════════
# if len(sys.argv) > 1:
#     sim_name = sys.argv[1]   #^ Name of the simulation

#* Root of the simulation directories
root_sim     = base_data + sim_series + sim_name + '/'
root_package = base_data + package_name + '/'               #^ Template package (inputs, run scripts)

#* Set relevant I/O directories
root_out_mount = base_out_mount + sim_series + sim_name + '/' + out_dir + '/'
root_out       = base_out + sim_series + sim_name + '/' + out_dir + '/'

if flag_output:
    create_directory(root_sim)
    create_directory(root_out_mount)
    create_directory(root_out)

log.info(f"Simulation Directory: {root_sim} (frame {frame_IC})")
log.info(f"Template Package    : {root_package}")
log.info(f'Output location     : {root_out}')
log.info(f'Output (mount)      : {root_out_mount}')
log.info('=' * 70)

# %%
#- ═════════════════════════════════════════════════════════════════════════════
#- GRID CONFIGURATION
#- ═════════════════════════════════════════════════════════════════════════════
#* Import grid operations
from grid_operations import read_grid_data, get_grid_wHalos, get_grid_center_wHalos, get_staggered_grid

#* Set relevant directories
grid_loc         = root_grid + grid_name + '/'

#* Load the grid parameters file
with open(grid_loc + 'grid_parameters.yaml') as file:
    globals().update(yaml.safe_load(file))

#* Read grid.data
_, _, xx, yy, zz = read_grid_data(grid_loc + 'grid.data', ngx, ngy, ngz, Lx, Ly, Lz)
z0 = globals().get('z0', 0.0)  #^ Not every grid_parameters.yaml has z0
zz = zz + z0  #^ z grid always starts from 0

#* Wall locations
y_wall_bot = yy[0]
y_wall_top = yy[-1]

#* Convert to wall units
xp = xx / delta_nu
yp = (yy - y_wall_bot) / delta_nu
zp = zz / delta_nu

#* Cell faces (with Halos) -> edges of the vector potential
xx_fc = get_grid_wHalos(xx)
yy_fc = get_grid_wHalos(yy)
zz_fc = get_grid_wHalos(zz)

#* Cell centers (with Halos)
xx_cc = get_grid_center_wHalos(xx)
yy_cc = get_grid_center_wHalos(yy)
zz_cc = get_grid_center_wHalos(zz)

#* Staggered grid (with Halos)
xu, yu, zu = get_staggered_grid(xx, yy, zz, 1)
xv, yv, zv = get_staggered_grid(xx, yy, zz, 2)
xw, yw, zw = get_staggered_grid(xx, yy, zz, 3)

log.info(f"Grid: {grid_name}, ngx = {ngx}, ngy = {ngy}, ngz = {ngz}")
log.info(f"dx+ = {xp[1]-xp[0]:.2f}, dy+_wall = {yp[1]:.3f}, dz+ = {zp[1]-zp[0]:.2f}")

# %%
#- ═════════════════════════════════════════════════════════════════════════════
#- HAIRPIN VORTEX PARAMETERS
#- ═════════════════════════════════════════════════════════════════════════════
#@ Centreline: X(q) = X_head - (H sin(alpha), s_w H cos(alpha), 0) [1 - exp(-q^2/2)] + (0, 0, B q)
#@ Core:       omega = Gamma g(q) / (pi delta^2) exp(-r^2/delta^2) T(q),  g(q) = exp(-(q/q_m)^6)
#@             max(u_theta) = U_ref at r = 1.1209 delta (slender limit)

#* Peak of (1 - exp(-eta^2)) / eta  ->  max(u_theta) = C_LO * Gamma / (2 pi delta)
_eta  = np.linspace(0.5, 2.0, 200001)
C_LO  = np.max((1.0 - np.exp(-_eta**2)) / _eta)            #^ ~ 0.6382

#* Strength
U_ref    = sense * strength * u_tau                          #^ Signed peak swirl velocity
Gamma    = 2.0 * np.pi * spread * U_ref / C_LO               #^ Circulation of the core
omega_0  = Gamma / (np.pi * spread**2)                       #^ Peak vorticity on the centreline
a_strain = 4.0 * nu / spread**2                              #^ Equivalent Burgers strain rate (delta = sqrt(4 nu / a))

#* Head location (x, z default: middle of the domain)
x_c = 0.5 * (xx[0] + xx[-1]) if x_core is None else x_core
z_c = 0.5 * (zz[0] + zz[-1]) if z_core is None else z_core
y_c = y_wall_bot + y_core if wall == 'bottom' else y_wall_top - y_core
y_wall_near = y_wall_bot if wall == 'bottom' else y_wall_top
s_wall      = 1.0 if wall == 'bottom' else -1.0              #^ Legs point towards the near wall

#* Hairpin geometry
alpha   = np.deg2rad(incline)
kappa_0 = height / span**2                                   #^ Head curvature
y_base  = y_core - height * np.cos(alpha)                    #^ Height of the leg ends from the wall (q -> infinity)

log.info(f" Hairpin vortex (Burgers core):")
log.info(f" --------------------")
log.info(f"wall       = {wall}, sense = {sense}")
log.info(f"U_ref      = {U_ref:.4e}  ({sense * strength:.2f} u_tau)")
log.info(f"Gamma      = {Gamma:.4e}")
log.info(f"omega_0    = {omega_0:.4e}  (omega_0+ = {omega_0 * nu / u_tau**2:.3f})")
log.info(f"spread     = {spread:.4f}  ({spread / delta_nu:.1f} wall units)")
log.info(f"y_core     = {y_core:.4f}  ({y_core / delta_nu:.1f} wall units, head)")
log.info(f"head       = ({x_c:.4f}, {y_c:.4f}, {z_c:.4f})")
log.info(f"H, B       = {height:.4f}, {span:.4f}  ({height / delta_nu:.1f}, {span / delta_nu:.1f} wall units)")
log.info(f"alpha, q_m = {incline:.1f} deg, {q_max:.2f}")
log.info(f"kappa_0    = {kappa_0:.4f}  (R_0 = {1.0 / kappa_0:.4f}, kappa_0 * delta = {kappa_0 * spread:.3f})")
log.info(f"y_base     = {y_base:.4f}  ({y_base / delta_nu:.1f} wall units, leg ends)")
if kappa_0 * spread > 0.5:
    log.warning(f"kappa_0 * delta = {kappa_0 * spread:.2f} > 0.5: core not slender, increase span or reduce height")
if y_base < spread:
    log.warning(f"Leg ends at {y_base:.3f} < delta from the wall: legs overlap the wall damping, raise y_core or reduce height")

#- ═════════════════════════════════════════════════════════════════════════════
#- INITIALIZATION COMPLETE
#- ═════════════════════════════════════════════════════════════════════════════

log.info('Framework loaded successfully. Ready for processing...')
log.info('=' * 70)

# %%
