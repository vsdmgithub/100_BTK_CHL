import os
import numpy as np
import logging
import sys
import seaborn as sns
import matplotlib.pyplot as plt

# Define logging levels for different modes
LOG_LEVELS = {
    "normal": logging.INFO,
    "debug": logging.DEBUG,
    "production": logging.WARNING
}

def setup_logging(mode="normal"):
    """Set up logging to print INFO to stdout and ERROR to stderr."""
    log_level = LOG_LEVELS.get(mode, logging.INFO)

    # Create separate handlers
    stdout_handler = logging.StreamHandler(sys.stdout)  # INFO & DEBUG to stdout
    stderr_handler = logging.StreamHandler(sys.stderr)  # WARNING & ERROR to stderr

    # Set levels for handlers
    stdout_handler.setLevel(logging.DEBUG)   # Captures DEBUG & INFO
    stderr_handler.setLevel(logging.WARNING) # Captures WARNING, ERROR, CRITICAL

    # Define log format
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    stdout_handler.setFormatter(formatter)
    stderr_handler.setFormatter(formatter)

    # Get the root logger and configure it
    logger = logging.getLogger()
    logger.setLevel(log_level)
    logger.handlers = []  # Clear existing handlers
    logger.addHandler(stdout_handler)
    logger.addHandler(stderr_handler)

    return logger

if __name__ == "__main__":
    # Example: Choose mode
    logger = setup_logging(mode="debug")

    # Example log messages
    logger.debug("This is a debug message.")    # Goes to stdout
    logger.info("This is an info message.")     # Goes to stdout
    logger.warning("This is a warning message.")# Goes to stderr
    logger.error("This is an error message.")   # Goes to stderr
    logger.critical("This is a critical message.") # Goes to stderr

def create_directory(directory_path):
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)
        print(f"Directory '{directory_path}' created successfully.")
    else:
        print(f"Directory '{directory_path}' already exists.")

def print_progress(i, total_iterations, interval=10):
    """
    Print the progress of a loop at specified intervals.

    Parameters:
    i (int): Current iteration index.
    total_iterations (int): Total number of iterations.
    interval (int): Percentage interval for status update. Default is 10.
    """
    if total_iterations == 0:
        return  # Avoid division by zero
    progress_percentage = ((i + 1) / total_iterations) * 100
    
    if progress_percentage % interval == 0 or (i + 1) == total_iterations:
        print(f'Progress: {progress_percentage:.0f}%')

def print_line(status):
# ✅ Task completed
# ❌ Error in processing data
# 📝 Note: Remember to update this function
# ⚠️ Warning: This might break the code
# ☠️ Danger: Do not modify this line
# 🧪Testing in progress...",
# ⏳ In Progress: Loading data
# 🔴 Stop execution if input is invalid
# 🟢 Ready to proceed
# 🟡 Caution: Check the values before running
    """Prints a formatted line based on a single-word input status."""
    symbols = {
        "tick": "✅ Task completed",
        "wrong": "❌ Error detected",
        "stickynote": "🗒️ Important note",
        "note": "📝 Reminder: Update this section",
        "warning": "⚠️ Warning: Check your input",
        "danger": "☠️ Critical Issue: Proceed with caution!",
        "progress": "⏳ In Progress...",
        "testing": "🧪 Testing in progress...",
        "red": "🔴 Stop: Review required",
        "green": "🟢 All good to proceed",
        "yellow": "🟡 Caution: Verify before continuing",
    }

    # Print the corresponding line or a default message
    print(symbols.get(status.lower(), "ℹ️ Unknown status: No matching symbol found"))

def plot_palette(palette_name, n_colors=8):
    palette = sns.color_palette(palette_name, n_colors)
    sns.palplot(palette)
    plt.title(palette_name)
    plt.show()

# Try different palettes
if __name__ == "__main__":
    for p in ["colorblind", "muted", "deep", "dark", "cubehelix", "viridis"]:
        plot_palette(p)

def get_bcname(flag):
    if flag == 1:
        return 'bot'
    elif flag == 2:
        return 'top'
    elif flag == 3:
        return 'lft'
    elif flag == 4:
        return 'rht'
    elif flag == 5:
        return 'bck'
    elif flag == 6:
        return 'frt'
    else:
        print('Invalid flag')
        exit()