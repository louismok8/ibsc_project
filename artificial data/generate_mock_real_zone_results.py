# analysis/generate_mock_real_zone_results.py

import pandas as pd
import numpy as np
import os

# -------------------------------
# Parameters
# -------------------------------

n_patients = 435
zones = ["left_anterior", "left_posterior", "right_anterior", "right_posterior"]
output_file = "data/real_zone_results.csv"

# -------------------------------
# Generate patient IDs
# -------------------------------

patient_ids = [str(i+1).zfill(3) for i in range(n_patients)]

# -------------------------------
# Generate mock biopsy results
# -------------------------------

np.random.seed(42)  # reproducible

rows = []

for pid in patient_ids:
    for zone in zones:
        # Randomly assign positive/negative outcome (~30% positive rate)
        biopsy_positive = np.random.binomial(1, 0.3)
        rows.append({
            "patient": pid,
            "zone": zone,
            "biopsy_positive": biopsy_positive
        })

# -------------------------------
# Save CSV
# -------------------------------

df = pd.DataFrame(rows)

# Ensure directory exists
os.makedirs("data", exist_ok=True)

df.to_csv(output_file, index=False)

print(f"✓ Mock real_zone_results.csv created with {n_patients*len(zones)} rows")