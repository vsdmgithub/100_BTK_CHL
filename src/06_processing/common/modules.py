# %% [markdown]
#& ═════════════════════════════════════════════════════════════════════════════
#& MODULE IMPORTS - Organized for HPC-CFD Scientific Computing
#& ═════════════════════════════════════════════════════════════════════════════


#- ═════════════════════════════════════════════════════════════════════════════
#- SYSTEM & CORE MODULES
#- ═════════════════════════════════════════════════════════════════════════════
import os                # Operating system interface - file paths, environment variables
import sys               # System-specific parameters and functions - path manipulation
import time as tt        # Time access and conversions - performance timing
import importlib         # Dynamic module importing - reload modules during development
# import gc                # Garbage collector interface - memory management for large datasets
# import warnings          # Warning control - suppress non-critical warnings
from pathlib import Path # Object-oriented filesystem paths - modern path handling

#* Environment Configuration (NumExpr threading)
os.environ["NUMEXPR_MAX_THREADS"] = "32"  # Max threads for NumExpr operations
os.environ["NUMEXPR_NUM_THREADS"] = "32"  # Active threads for NumExpr
#^ Note: Adjust based on your HPC node specifications for optimal performance

#! WARNING: Path cleanup - remove unwanted paths before importing custom modules
unwanted = '/home/vulcan/OneDrive/NUMERICAL_COMPUTATION/channel_subdomain_solver_2025/6_DA_subdomain_turbulent_channel/src/12_processing/subdomain_channel'
if unwanted in sys.path:
    sys.path.remove(unwanted)

#- ═════════════════════════════════════════════════════════════════════════════
#- MATHEMATICAL & NUMERICAL COMPUTATION
#- ═════════════════════════════════════════════════════════════════════════════

#* Core Numerical Libraries
import numpy as np                                    # Array operations, linear algebra, FFT
import scipy as sc                                    # Scientific computing algorithms
# from scipy import linalg                              # Advanced linear algebra operations
# from scipy import optimize                            # Optimization and root finding
# from scipy import integrate                           # Numerical integration and ODE solvers
# from scipy import signal                              # Signal processing, convolution, filtering
# from scipy import fft                                 # Fast Fourier Transform operations
# from scipy import sparse                              # Sparse matrix operations for large systems

#* Interpolation Tools (Critical for CFD post-processing)
import scipy.interpolate as si                        # Interpolation module namespace
from scipy.interpolate import interpn                 # N-dimensional interpolation
from scipy.interpolate import UnivariateSpline        # 1D spline interpolation with smoothing
# from scipy.interpolate import RegularGridInterpolator # Fast structured grid interpolation
# from scipy.interpolate import LinearNDInterpolator   # Unstructured N-D linear interpolation
# from scipy.interpolate import RBFInterpolator        # Radial basis function interpolation
# from scipy.interpolate import BSpline as bsp        # B-spline basis functions
# from scipy.interpolate import CubicSpline           # Cubic spline with boundary conditions
from scipy.optimize import curve_fit                  # Non-linear curve fitting

#* Statistical Analysis
# from scipy import stats                               # Statistical functions and distributions

#* Symbolic Mathematics (Optional - for analytical work)
# import sympy as sym                                 # Symbolic mathematics, equation solving


#- ═════════════════════════════════════════════════════════════════════════════
#- DATA STRUCTURES & MANIPULATION
#- ═════════════════════════════════════════════════════════════════════════════
# import pandas as pd                                   # DataFrames for tabular data, time series
# import xarray as xr                                 # N-D labeled arrays (NetCDF, HDF5) - excellent for CFD grids
# import dask                                         # Parallel computing with larger-than-memory arrays
# import dask.array as da                             # Distributed NumPy arrays for HPC
# import h5py                                         # HDF5 binary data format - standard for large CFD datasets


#- ═════════════════════════════════════════════════════════════════════════════
#- VISUALIZATION & PLOTTING
#- ═════════════════════════════════════════════════════════════════════════════

#* Matplotlib Core
import matplotlib as mpl                              # Matplotlib base configuration
import matplotlib.pyplot as plt                       # MATLAB-like plotting interface
from matplotlib import ticker                         # Axis tick customization
from matplotlib.ticker import (                       # Fine-grained tick control
    MultipleLocator,      # Set tick spacing at multiples of base
    FormatStrFormatter,   # Custom tick label formatting
    AutoMinorLocator,     # Automatic minor tick placement
    ScalarFormatter,      # Scalar number formatting
    LogLocator,           # Logarithmic axis tick placement
    FuncFormatter         # Custom function-based formatting
)

#* Advanced Layout Tools
from mpl_toolkits.axes_grid1 import make_axes_locatable # Colorbar positioning, axis division
from matplotlib.gridspec import GridSpec              # Advanced subplot grid layouts
# from mpl_toolkits.mplot3d import Axes3D             # 3D plotting capabilities

#* Colormaps and Styling
from matplotlib import cm                             # Colormap utilities
# import seaborn as sns                               # Statistical visualization, better default styles
# sns.set_theme()                                     # Apply seaborn styling globally

#* Advanced Visualization (CFD-specific)
# import pyvista as pv                                # 3D visualization for CFD meshes/surfaces
# import mayavi.mlab as mlab                          # Advanced 3D scientific visualization
# import plotly.graph_objects as go                   # Interactive web-based plots
# import holoviews as hv                              # High-level declarative visualization


#- ═════════════════════════════════════════════════════════════════════════════
#- FILE I/O & DATA FORMATS
#- ═════════════════════════════════════════════════════════════════════════════
import yaml                                           # YAML file reading/writing - config files
# import json                                           # JSON data interchange format
# import pickle                                       # Python object serialization
# import toml                                         # TOML configuration files
# from netCDF4 import Dataset                         # NetCDF file format - standard in atmospheric CFD
# import vtk                                          # VTK file format for ParaView/VisIt

#* CFD-Specific File Formats
import tecplot as tp                                # Tecplot data format and scripting
from tecplot.constant import *                      # Tecplot constants
# import meshio                                       # Universal mesh I/O (supports 20+ formats)
# import pyFoam                                       # OpenFOAM case manipulation and parsing


#- ═════════════════════════════════════════════════════════════════════════════
#- LOGGING & DEBUGGING
#- ═════════════════════════════════════════════════════════════════════════════
import logging                                        # Flexible logging framework
import logging.config                                 # Logging configuration from files/dicts
# from loguru import logger                           # Modern logging with better defaults

#^ Logging Configuration Example
# logging.basicConfig(
#     level=logging.INFO,
#     format='%(asctime)s | %(levelname)-8s | %(message)s',
#     datefmt='%Y-%m-%d %H:%M:%S'
# )


#- ═════════════════════════════════════════════════════════════════════════════
#- PARALLEL COMPUTING & PERFORMANCE
#- ═════════════════════════════════════════════════════════════════════════════
from joblib import Parallel, delayed                  # Simple parallel for-loops, disk caching
from tqdm import tqdm                                 # Progress bars for loops and iterators
# from tqdm.auto import tqdm                          # Auto-selects notebook vs terminal version

#* Advanced Parallelization
# from multiprocessing import Pool, cpu_count         # Process-based parallelism
# from concurrent.futures import ThreadPoolExecutor   # Thread-based parallelism for I/O
# import mpi4py                                       # MPI bindings for distributed HPC computing
# from mpi4py import MPI                              # Message Passing Interface

#* Performance Profiling
# import cProfile                                     # Performance profiling
# import line_profiler                                # Line-by-line profiling
# import memory_profiler                              # Memory usage profiling

#* JIT Compilation (Speed up Python code)
# import numba                                        # JIT compiler for NumPy code
# from numba import jit, njit, prange                 # Decorators for automatic compilation


#- ═════════════════════════════════════════════════════════════════════════════
#- UTILITY LIBRARIES
#- ═════════════════════════════════════════════════════════════════════════════
# from datetime import datetime, timedelta              # Date and time manipulation
# from collections import defaultdict, Counter          # Enhanced dictionary types
# from itertools import product, combinations           # Iterator tools for combinatorics
# import re                                             # Regular expressions for text parsing
import shutil                                       # High-level file operations (copy, move, delete)
# from glob import glob                               # Unix-style pathname pattern expansion
# import argparse                                     # Command-line argument parsing


#- ═════════════════════════════════════════════════════════════════════════════
#- CFD-SPECIFIC LIBRARIES (Commonly Used)
#- ═════════════════════════════════════════════════════════════════════════════

#* Turbulence Analysis
# import pyLOM                                        # Low-order modeling, POD, DMD analysis
# import modred                                       # Model reduction methods (POD, DMD, BPOD)
# from pyDMD import DMD, HODMD                        # Dynamic Mode Decomposition variants

#* Grid Generation & Manipulation
# import pygmsh                                       # Python interface for Gmsh mesh generator
# import meshio                                       # Universal mesh file converter
# from scipy.spatial import Delaunay, KDTree          # Computational geometry tools

#* Post-Processing Utilities
# import pyevtk                                       # Export VTK files for ParaView
# from paraview.simple import *                       # ParaView Python scripting
# import visit_utils                                  # VisIt Python interface

#* Machine Learning for CFD (Data-Driven Approaches)
# import sklearn                                      # Scikit-learn for ML algorithms
# from sklearn.decomposition import PCA               # Principal Component Analysis
# import tensorflow as tf                             # Deep learning for RANS closure models
# import torch                                        # PyTorch for neural network training


#- ═════════════════════════════════════════════════════════════════════════════
#- CONFIGURATION & BEST PRACTICES
#- ═════════════════════════════════════════════════════════════════════════════

#^ NumPy Print Options (for better array display)
np.set_printoptions(
    precision=6,      # Decimal precision for floating point output
    suppress=True,    # Suppress scientific notation for small numbers
    linewidth=120,    # Characters per line before wrapping
    threshold=50      # Total array elements before summarizing
)

#^ Matplotlib Configuration (for publication-quality figures)
# mpl.rcParams['figure.dpi'] = 100                      # Display resolution
# mpl.rcParams['savefig.dpi'] = 300                     # Save resolution for publications
# mpl.rcParams['font.size'] = 11                        # Base font size
# mpl.rcParams['font.family'] = 'serif'                 # Font family (serif, sans-serif)
# mpl.rcParams['text.usetex'] = True                  # Use LaTeX for text rendering
# mpl.rcParams['axes.labelsize'] = 12                   # Axis label font size
# mpl.rcParams['axes.titlesize'] = 13                   # Title font size
# mpl.rcParams['legend.fontsize'] = 10                  # Legend font size
# mpl.rcParams['xtick.labelsize'] = 10                  # X-tick label size
# mpl.rcParams['ytick.labelsize'] = 10                  # Y-tick label size
# mpl.rcParams['lines.linewidth'] = 1.5                 # Default line width
# mpl.rcParams['lines.markersize'] = 6                  # Default marker size

#! Warning Control (suppress common numerical warnings if needed)
# warnings.filterwarnings('ignore', category=RuntimeWarning)
# warnings.filterwarnings('ignore', category=FutureWarning)
