"""
modules.py
===========
Centralized Python module imports for CFD simulations, post-processing, data analysis,
and visualization.

✅ Categories:
  - Math / Numerical
  - Plotting / Visualization
  - System / OS / Logging
  - CFD-specific Tools
  - Parallel Processing / HPC
  - Data I/O
  - Optional Machine Learning & Post-processing Tools
"""

# ==================================================
# 🟢 MATH / NUMERICAL
# ==================================================
import numpy as np                     # Core numerical array operations

# Interpolation tools
import scipy as sc
import scipy.interpolate as si
from scipy.interpolate import interpn
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import SmoothBivariateSpline
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import griddata
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import BSpline as bsp, make_interp_spline

# FFT and signal processing
from scipy.fft import fftn, ifftn
from scipy.signal import butter, filtfilt, find_peaks

# Linear algebra & sparse solvers
from numpy.linalg import inv, solve, eig
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve

# ==================================================
# 🟢 PLOTTING / VISUALIZATION
# ==================================================
import matplotlib.pyplot as plt        # Standard plotting
import matplotlib.cm as cm             # Colormap utilities
import matplotlib.ticker as ticker     # Tick control
from matplotlib import ticker
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator
from matplotlib.ticker import ScalarFormatter

from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns                      # Statistical visualization

# Interactive and 3D Visualization
import pyvista as pv                   # 3D visualization for CFD post-processing
import vedo                            # Lightweight mesh/volume rendering (VTK, OBJ files)

# Animations
import matplotlib.animation as animation
import plotly.express as px
import plotly.graph_objects as go

# ==================================================
# 🟢 SYSTEM / OS / LOGGING
# ==================================================
import os
import sys
import importlib
import time as tt                     # Timing utilities
from pathlib import Path              # Modern path handling

# Logging setup
import logging as lg
import logging.config

# Parallel execution
from tqdm import tqdm                 # Progress bars
import multiprocessing                # Multi-core processing
from joblib import Parallel, delayed  # Easy parallel loops
from mpi4py import MPI                # MPI for distributed computing
import numba                          # JIT compilation
import cupy                           # GPU-accelerated array computations (NVIDIA CUDA)

# ==================================================
# 🟢 DATA I/O
# ==================================================
import h5py                           # HDF5 file support
import netCDF4                        # netCDF support
import pandas as pd                   # Tabular data / CSV
import csv                            # CSV read/write
import json                           # JSON read/write
import yaml                           # YAML file parsing
import pickle                         # Python object serialization

# ==================================================
# 🟢 CFD-SPECIFIC TOOLS
# ==================================================
import tecplot as tp                  # Tecplot Python API (Post-processing)
from tecplot.constant import *        # Tecplot constants (for zones, variables, etc.)

# Mesh & Geometry Processing
import meshio                         # Mesh file I/O (VTK, GMSH, etc.)
import pygmsh                         # Mesh generation
import fenics                         # Finite Element PDE solving
import dolfin                         # Core of FEniCS
import firedrake                      # Alternative FEM solver

# ==================================================
# 🟢 MACHINE LEARNING & DATA ANALYSIS
# ==================================================
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import tensorflow as tf
import torch
import keras

# ==================================================
# 🟢 IMAGE / VIDEO PROCESSING
# ==================================================
from PIL import Image
import imageio
import cv2                            # OpenCV for videos and image processing

# ==================================================
# 🟢 WARNINGS / EXCEPTIONS
# ==================================================
import warnings
warnings.filterwarnings("ignore")     # Disable warnings (use with caution)

