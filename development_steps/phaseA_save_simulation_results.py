"""
Phase A — Save simulation outputs

This script runs the simulation pipeline and saves
lesion-level results to a structured CSV file.

Output
------
outputs/simulation_results.csv
"""

import os
import sys
import numpy as np
import pandas as pd


# ---------------------------------------------------------
# Ensure project root is on path
# ---------------------------------------------------------

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(PROJECT_ROOT)


# ---------------------------------------------------------
# Import pipeline
# ---------------------------------------------------------

from pipelines.run_simulation_pipeline import run_simulation_pipeline
from pipelines.run_simulation_pipeline import load_full_dataset


# ---------------------------------------------------------
# Output path
# ---------------------------------------------------------

OUTPUT_DIR = os.path.join(PROJECT_ROOT, "outputs")
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "simulation_results.csv")

os.makedirs(OUTPUT_DIR, exist_ok=True)


# ---------------------------------------------------------
# Run simulation
# ---------------------------------------------------------

print("Running simulation pipeline...")

results = run_simulation_pipeline()

print("\nSimulation complete. Converting results to table...")


# ---------------------------------------------------------
# Load dataset again (for lesion metadata)
# ---------------------------------------------------------

dataset = load_full_dataset()

dataset.verify_all_geometry()
dataset.extract_all_lesions()
dataset.compute_all_lesion_volumes()
dataset.compute_all_lesion_centroids()


# ---------------------------------------------------------
# Build table
# ---------------------------------------------------------

rows = []

for patient_result in results:

    patient_id = patient_result["patient_id"]
    lesion_stats = patient_result["lesion_results"]

    patient = next(p for p in dataset if p.id == patient_id)

    for lesion in patient.lesions:

        stats = lesion_stats[lesion.id]

        expected_positive_cores = np.mean(
            stats["distribution_positive_core_counts"]
        )

        rows.append({
            "patient_id": patient_id,
            "lesion_id": lesion.id,
            "lesion_volume_mm3": lesion.volume,
            "hit_probability": stats["hit_probability"],
            "mean_percent_positive": stats["mean_percentage_positive"],
            "expected_positive_cores": expected_positive_cores
        })


df = pd.DataFrame(rows)

df.to_csv(OUTPUT_FILE, index=False)

print(f"\n✓ Saved simulation results to:\n{OUTPUT_FILE}")
print(f"\nTotal rows: {len(df)}")