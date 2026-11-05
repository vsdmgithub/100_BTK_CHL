import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib import ticker
from seaborn import color_palette
from cycler import cycler

# ============================================================
# 1️⃣ BASIC — Minimal style (like your original version)
# ============================================================
import matplotlib.pyplot as plt

def plot_settings():
    """
    Configure matplotlib settings for consistent, publication-quality plots.
    Based on W01_compute_utau/main.py
    """
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "font.weight": "bold",
        "axes.labelweight": "bold",
        "axes.titleweight": "bold",
        "axes.linewidth": 1.5,
        "text.latex.preamble": r"\usepackage{sfmath}"
    })

def plot_settings_JFM():
    """
    Configure matplotlib settings for Journal of Fluid Mechanics (JFM)
    publication-quality figures.

    JFM figure specifications:
    - Single column width : 84 mm  (3.31 in)
    - Double column width : 174 mm (6.85 in)
    - Body font           : Times Roman (serif), ~10 pt
    - Caption font        : ~8 pt
    - Ticks               : inward, mirrored on all four sides
    - Grid                : off (JFM convention)
    - Save dpi            : >= 300 (600 recommended for line art)

    LaTeX packages: newtxtext + newtxmath give Times-compatible
    text and maths that match the JFM typeset font.
    """
    plt.rcParams.update({
        # -- LaTeX rendering (Times-compatible serif) --------------------
        "text.usetex"          : True,
        "font.family"          : "serif",
        "font.serif"           : ["Times", "Times New Roman",
                                  "Computer Modern Roman"],
        "text.latex.preamble"  : r"\usepackage{newtxtext,newtxmath}\usepackage{mathrsfs}",

        # -- Figure geometry (single-column default) ---------------------
        "figure.figsize"       : (3.31, 2.56),   # 84 mm × 65 mm
        "figure.dpi"           : 150,
        "savefig.dpi"          : 600,
        "savefig.bbox"         : "tight",
        "savefig.pad_inches"   : 0.02,

        # -- Font sizes (JFM body ~10 pt) --------------------------------
        "font.size"            : 9,
        "axes.labelsize"       : 10,
        "axes.titlesize"       : 10,
        "legend.fontsize"      : 8,
        "xtick.labelsize"      : 9,
        "ytick.labelsize"      : 9,

        # -- Axes --------------------------------------------------------
        "axes.linewidth"       : 0.75,
        "axes.labelweight"     : "bold",

        # -- Ticks (inward, mirrored — JFM/fluid-mechanics standard) ----
        "xtick.direction"      : "in",
        "ytick.direction"      : "in",
        "xtick.top"            : True,
        "ytick.right"          : True,
        "xtick.major.size"     : 3.5,
        "ytick.major.size"     : 3.5,
        "xtick.minor.size"     : 2.0,
        "ytick.minor.size"     : 2.0,
        "xtick.major.width"    : 0.75,
        "ytick.major.width"    : 0.75,
        "xtick.minor.width"    : 0.5,
        "ytick.minor.width"    : 0.5,

        # -- Lines and markers -------------------------------------------
        "lines.linewidth"      : 1.0,
        "lines.markersize"     : 4.0,

        # -- Legend ------------------------------------------------------
        "legend.frameon"       : False,
        "legend.handlelength"  : 1.5,
        "legend.borderpad"     : 0.3,
        "legend.labelspacing"  : 0.3,

        # -- Grid (off by default for JFM) -------------------------------
        "axes.grid"            : False,
    })


def plot_settings_woLatex():
    """
    Configure matplotlib settings for consistent, publication-quality plots.
    Cluster-compatible version without LaTeX dependency.
    """
    plt.rcParams.update({
        "text.usetex": False,  # Disabled LaTeX for cluster compatibility
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        # "font.weight": "bold",
        # "axes.labelweight": "bold",
        # "axes.titleweight": "bold",
        "axes.linewidth": 1.5,
        "mathtext.fontset": "dejavusans"  # Use mathtext instead of LaTeX
    })

def apply_plot_style():
    """
    Minimal serif-font style setup for cluster plotting (no LaTeX dependency).
    """
    plt.rcParams.update({
        "text.usetex": False,                 # Disable LaTeX
        "font.family": "serif",               # Use built-in serif
        "font.serif": ["Times New Roman", "DejaVu Serif", "Times"],
        "mathtext.fontset": "dejavuserif",    # Use math symbols that match serif
        "axes.labelsize": 14,
        "axes.titlesize": 15,
        "font.size": 13,
        "legend.fontsize": 12,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "figure.dpi": 150,
    })

def apply_plot_style_wLatex1():
     """
     Minimal LaTeX-compatible setup for quick plotting.
     """
     plt.rcParams.update({
         "text.usetex": True,
         "font.family": "serif",
         "font.serif": ["Computer Modern Roman"],
         "text.latex.preamble": r"\usepackage{newtxmath}",
     })

def apply_plot_style_wLatex2():
   """Journal-style with LaTeX fonts."""
   plt.rcParams.update({
       "text.usetex": True,
       "font.family": "serif",
       "font.serif": ["Times"],
      "text.latex.preamble": r"\usepackage{siunitx}\usepackage{newtxtext,newtxmath}",
        "axes.labelsize": 12,
       "font.size": 12,
       "legend.fontsize": 10,
       "xtick.labelsize": 10,
       "ytick.labelsize": 10,
       "axes.linewidth": 1,
   })

# ============================================================
# 2️⃣ LINE PLOT STYLE — For profile, convergence, or error plots
# ============================================================

def apply_lineplot_style():
    """
    Line plot style — best for profiles, time series, spectra, etc.
    """
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
        "text.latex.preamble": r"\usepackage{newtxtext,newtxmath}",
        
        # Figure setup
        "figure.dpi": 200,
        "figure.figsize": (6, 4),
        "savefig.dpi": 300,
        "figure.autolayout": True,

        # Axes
        "axes.labelsize": 14,
        "axes.titlesize": 14,
        "axes.linewidth": 1.0,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "grid.linestyle": "--",
        "grid.linewidth": 0.5,

        # Ticks
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.size": 6,
        "ytick.major.size": 6,
        "xtick.minor.size": 3,
        "ytick.minor.size": 3,

        # Legend
        "legend.fontsize": 12,
        "legend.frameon": False,
        "legend.handlelength": 2.0,

        # Lines
        "lines.linewidth": 1.8,
        "lines.markersize": 5,

        # Save
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
    })

    # Consistent color cycle
    plt.rcParams["axes.prop_cycle"] = cycler(color=[
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
        "#9467bd", "#8c564b", "#e377c2", "#7f7f7f"
    ])


# ============================================================
# 3️⃣ CONTOUR / IMAGE PLOTS — For fields, λ₂, or u' visualizations
# ============================================================

def apply_contour_style():
    """
    Contour/image style — better color scaling, no grids, LaTeX labels.
    """
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
        "text.latex.preamble": r"\usepackage{newtxtext,newtxmath}",

        "figure.dpi": 200,
        "figure.figsize": (6, 3),
        "savefig.dpi": 300,

        # Axes & ticks
        "axes.labelsize": 14,
        "axes.titlesize": 14,
        "axes.linewidth": 1.0,
        "axes.grid": False,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,

        # Image/contour
        "image.cmap": "viridis",
        "image.interpolation": "nearest",

        # Colorbar
        "axes.titlepad": 6,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
    })


# ============================================================
# 4️⃣ CUSTOM STYLE — Advanced, configurable
# ============================================================

def apply_custom_style(
    usetex=True,
    fontsize=14,
    fontfamily="serif",
    colormap="viridis",
    figsize=(6, 4),
    grid=True,
    linewidth=1.8,
    markersize=5,
    dpi=200,
    tex_preamble=r"\usepackage{newtxtext,newtxmath}"
):
    """
    Advanced configurable style for publication-quality figures.
    Example:
        apply_custom_style(fontfamily='sans-serif', colormap='plasma')
    """
    plt.rcParams.update({
        "text.usetex": usetex,
        "font.family": fontfamily,
        "font.serif": ["Computer Modern Roman"],
        "text.latex.preamble": tex_preamble,

        "figure.dpi": dpi,
        "figure.figsize": figsize,
        "figure.autolayout": True,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,

        "axes.labelsize": fontsize,
        "axes.titlesize": fontsize,
        "axes.linewidth": 1.0,
        "axes.grid": grid,

        "xtick.labelsize": fontsize - 2,
        "ytick.labelsize": fontsize - 2,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,

        "legend.fontsize": fontsize - 2,
        "legend.frameon": False,
        "legend.handlelength": 2.0,

        "lines.linewidth": linewidth,
        "lines.markersize": markersize,
        "image.cmap": colormap,
    })

    plt.rcParams["axes.prop_cycle"] = cycler(color=[
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
        "#9467bd", "#8c564b", "#e377c2", "#7f7f7f"
    ])

def visualize_cartesian_grid(x, y, xlabel='x', ylabel='y',skip=1,gName=''):
    # Create a grid of points
    """
    Visualize a Cartesian grid with specified x and y coordinates, and save the grid plot as an image.

    Parameters:
    x (array-like): Array of x-coordinates for the grid lines.
    y (array-like): Array of y-coordinates for the grid lines.
    xlabel (str, optional): Label for the x-axis. Default is 'x'.
    ylabel (str, optional): Label for the y-axis. Default is 'y'.
    skip (int, optional): Interval for skipping points when plotting grid lines. Default is 1 (no skipping).
    gName (str, optional): Suffix for the saved image file name. Default is an empty string.

    The function creates a grid by plotting vertical and horizontal lines at specified intervals, 
    adds labels to the axes, and saves the plot as a PNG image with a resolution of 300 dpi.
    """
    plt.close()
    fig,ax=plt.subplots()
    # xx, yy = np.meshgrid(x, y)

    lw=0.25
    alpha=1.0
    # Plot the grid lines
    for i in range(0,len(x),skip):
        ax.plot([x[i]]*len(y), y, 'k-', lw=lw,alpha=alpha) # Vertical lines
    for j in range(0,len(y),skip):
        ax.plot(x,[y[j]]*len(x), 'k-', lw=lw,alpha=alpha) # Horizontal lines

    # Add labels to the axes
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    if (skip != 1):
        plt.title('Cartesian Grid, skipping every {} pts'.format(skip))
    else:
        plt.title('Cartesian Grid')

    # Show the plot
    plt.savefig('./Grid'+gName+'.png',dpi=300,bbox_inches='tight')
    plt.show()

def plot_contour(data, X, Y, cmap='viridis', 
                 xmin=None, xmax=None, ymin=None, ymax=None, 
                 vmin=None, vmax=None, xlabel='X-axis', ylabel='Y-axis', plot_name='contour_plot',save_dir = './'):
  """
  Function to plot a contour plot from 2D data with customizable options.

  Parameters:
  - data : 2D array
      The data to be plotted.
  - X, Y : 2D arrays
      Grids corresponding to the data.
  - cmap : str or Colormap, optional
      Colormap for the contour plot (default: 'viridis').
  - xmin, xmax : float, optional
      Minimum and maximum limits for the x-axis (default: None, determined from data).
  - ymin, ymax : float, optional
      Minimum and maximum limits for the y-axis (default: None, determined from data).
  - vmin, vmax : float, optional
      Minimum and maximum values for the colormap (default: None, determined from data).
  - xlabel, ylabel : str, optional
      Labels for the x-axis and y-axis (default: 'X-axis', 'Y-axis').
  - plot_name : str, optional
      Name of the plot (default: 'contour_plot').
  - save_dir : str, optional
      Name of the directory to save (default: './').
  """
  fig, ax = plt.subplots()
  # Contour plot with options

  lev = 400
  cont =0.5
  if vmin is None and vmax is None:
      vmax=cont * np.amax(data)
      vmin=cont * np.amin(data)

  contour = ax.contourf(X, Y, data, levels=lev,
                        cmap=cmap,
                        # extend='both',
                        origin='lower', 
                        vmin=vmin, vmax=vmax,
                        )

  # Adding a colorbar with options
  # cbar = fig.colorbar(contour, ax=ax, orientation='vertical',
  #                     fraction=0.046, pad=0.04, 
  #                     shrink=1.0, aspect=20)
  cbar = fig.colorbar(contour, ticks=np.linspace(vmin,vmax,5))

  # Colorbar labels and title
  # cbar.set_label('Intensity', rotation=270, labelpad=15)
  # cbar.ax.set_title('Title', fontsize=10)
  # Customize colorbar ticks for readability
  # cbar.ax.yaxis.set_major_locator(plt.MaxNLocator(nbins=4))  # Adjust the number of ticks as needed
  # cbar.ax.yaxis.set_major_locator(plt.LinearLocator(numticks=4))  # Adjust the number of ticks as needed

  # Axis labels and title 
  ax.set_xlabel(xlabel)
  ax.set_ylabel(ylabel)
  # ax.set_title(' ')

  # Axis limits
  if xmin is not None and xmax is not None:
      ax.set_xlim(xmin, xmax)
  if ymin is not None and ymax is not None:
      ax.set_ylim(ymin, ymax)

  # Customize ticks
  ax.xaxis.set_major_locator(plt.MaxNLocator(4))
  ax.yaxis.set_major_locator(plt.MaxNLocator(4))

  # Grid
  # ax.grid(True, linestyle='--', linewidth=0.5)

  # Adding annotations
  # ax.annotate('Peak', xy=(0, 0), xytext=(1, 1), 
              # arrowprops=dict(facecolor='black', shrink=0.05))

  # Display the plot
  plt.show()

  # Save the plot as a high-resolution PNG file
  plt.savefig(save_dir+'{}.png'.format(plot_name), dpi=600, bbox_inches='tight')
