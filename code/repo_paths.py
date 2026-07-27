"""
repo_paths.py

The logs and CSVs under results/ are committed alongside the code, so any path a
script prints has to read the same whatever machine produced it. Everything these
scripts write lives inside the repository, so printing the path relative to the
repository root is both shorter and reproducible.
"""

import os

ROOT = os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))


def rel(path):
    """Path relative to the repository root, e.g. 'results/emergence_sweep.csv'.

    Falls back to the normalised absolute path for anything outside the
    repository, which nothing here currently writes.
    """
    full = os.path.abspath(path)
    if os.path.commonpath([full, ROOT]) != ROOT:
        return os.path.normpath(full)
    return os.path.relpath(full, ROOT).replace(os.sep, "/")
