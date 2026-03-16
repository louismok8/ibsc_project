# analysis/generate_mock_real_patient_results.py

import pandas as pd
import numpy as np

# Number of patients
n_patients = 435

# Generate patient IDs: P001, P002, ..., P435
patient_ids = [f"{str(i+1).zfill(3)}" for i in range(n_patients)]

# Generate mock positive/negative outcomes (random 0/1)
# You can adjust probability if desired
np.random.seed(42)  # reproducible
positives = np.random.binomial(1, 0.3, size=n_patients)  # ~30% positive rate

# Create DataFrame
df = pd.DataFrame({
    "patient_id": patient_ids,
    "positive": positives
})

# Save CSV
df.to_csv("data/real_patient_results.csv", index=False)
print("✓ Mock real_patient_results.csv created with 435 patients")