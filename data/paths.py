import os

"""
Centralised dataset paths.

This file contains NO logic and NO data loading.
It only defines where data lives on disk.
"""

# Root directory of the project data
BASE_DIR = os.path.expanduser("~/Documents/iBSC/Project")

# Subdirectories following dataset convention
IMAGES_DIR = os.path.join(BASE_DIR, "imagesTr")
LABELS_DIR = os.path.join(BASE_DIR, "labelsTr")
ZONES_DIR = os.path.join(BASE_DIR, "zonesTr")
