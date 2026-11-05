"""
Configuration module for setting up home and work directory paths.

This module automatically detects the current machine/cluster environment
and sets appropriate base paths for home and work directories.

Supported environments:
    - Personal PC (username: sugan, vulcan)
    - Cluster 1: Rockfish (username: smuruga4, path contains 'tzaki1')
    - Cluster 2: HPCMP - Warhawk & Jean (username: smuruga4, path contains '/p/')
"""

import os
from pathlib import Path
from typing import Tuple


def get_environment_paths() -> Tuple[str, str]:
    """
    Determine base_home and base_work paths based on the current environment.

    Returns:
        Tuple[str, str]: A tuple containing (base_home, base_work) paths.

    Raises:
        RuntimeError: If the username or cluster environment cannot be determined.
        ValueError: If required environment variables are not set.
    """
    username = os.getenv("USER") or os.getenv("USERNAME")

    if not username:
        raise ValueError(
            "Unable to determine username. Neither USER nor USERNAME "
            "environment variables are set."
        )

    cwd = str(Path.cwd())

    # Configuration mapping for different environments
    if username == "sugan":
        base_home = "/home/sugan/OneDrive/NUMERICAL_COMPUTATION/"
        base_work = "/home/sugan/Rockfish_scratch4/"

    elif username == "vulcan":
        base_home = "/home/vulcan/OneDrive/NUMERICAL_COMPUTATION/"
        base_work = "/home/vulcan/Rockfish_scratch4/"

    elif username == "smuruga4":
        # Cluster 1: Rockfish scratch4
        if "tzaki1" in cwd:
            base_home = "/scratch4/tzaki1/smuruga4/"
            base_work = "/scratch4/tzaki1/smuruga4/"
        # Cluster 2: HPCMP Warhawk
        elif "/p/" in cwd:
            base_home = "/p/home/smuruga4/"
            base_work = "/p/work1/smuruga4/"
        else:
            raise RuntimeError(
                f"Unknown cluster environment for user '{username}'. "
                f"Current working directory: {cwd}\n"
                f"Expected 'tzaki1' or '/p/' in path."
            )
    else:
        raise RuntimeError(
            f"Unknown user: '{username}'. "
            f"Supported users are: sugan, vulcan, smuruga4"
        )

    return base_home, base_work


def validate_paths(base_home: str, base_work: str, strict: bool = False) -> None:
    """
    Validate that the configured paths exist.

    Args:
        base_home: Path to the home directory.
        base_work: Path to the work directory.
        strict: If True, raise an error if paths don't exist.
                If False, only print warnings.

    Raises:
        FileNotFoundError: If strict=True and any path doesn't exist.
    """
    for path_name, path_value in [("base_home", base_home), ("base_work", base_work)]:
        path_obj = Path(path_value)
        if not path_obj.exists():
            message = f"Warning: {path_name} does not exist: {path_value}"
            if strict:
                raise FileNotFoundError(message)
            else:
                print(f"[WARNING] {message}")
        elif not path_obj.is_dir():
            message = f"{path_name} exists but is not a directory: {path_value}"
            if strict:
                raise NotADirectoryError(message)
            else:
                print(f"[WARNING] {message}")


# Initialize paths
try:
    username = os.getenv("USER") or os.getenv("USERNAME")
    cwd = str(Path.cwd())
    base_home, base_work = get_environment_paths()

    # Optional: Validate paths (set strict=True to raise errors instead of warnings)
    validate_paths(base_home, base_work, strict=False)

except (RuntimeError, ValueError) as e:
    print(f"[ERROR] Failed to initialize configuration: {e}")
    raise


# Export variables for backward compatibility
__all__ = ["username", "cwd", "base_home", "base_work", "get_environment_paths", "validate_paths"]
