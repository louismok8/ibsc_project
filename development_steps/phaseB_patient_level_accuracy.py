# analysis/step21_patient_level_accuracy.py

import pandas as pd
from sklearn.metrics import confusion_matrix, accuracy_score, precision_score, recall_score, f1_score


# ---------------------------------------------------------
# File paths
# ---------------------------------------------------------
SIMULATION_CSV = "outputs/simulation_results.csv"
REAL_BIOPSY_CSV = "data/real_patient_results.csv"
PATIENT_PREDICTIONS_CSV = "outputs/patient_predictions.csv"

# ---------------------------------------------------------
# Load data
# ---------------------------------------------------------
sim_df = pd.read_csv(SIMULATION_CSV, dtype={"patient_id": str})
real_df = pd.read_csv(REAL_BIOPSY_CSV, dtype={"patient_id": str})
sim_df["patient_id"] = sim_df["patient_id"].str.zfill(3)
real_df["patient_id"] = real_df["patient_id"].str.zfill(3)

# ---------------------------------------------------------
# Aggregate lesion-level simulation to patient-level
# ---------------------------------------------------------
# Example aggregation:
# - patient_hit_probability = max(hit_probability across lesions)
# - patient_expected_positive_cores = sum(expected_positive_cores across lesions)

patient_agg = sim_df.groupby("patient_id").agg(
    patient_hit_probability=("hit_probability", "max"),
    patient_expected_positive_cores=("expected_positive_cores", "sum")
).reset_index()

# ---------------------------------------------------------
# Define simulation cut-off
# ---------------------------------------------------------
# Option 1: hit_probability > 0.5
patient_agg["predicted_positive"] = (patient_agg["patient_hit_probability"] > 0.5).astype(int)

# Option 2: expected_positive_cores >= 1
# Uncomment if you prefer this metric
# patient_agg["predicted_positive"] = (patient_agg["patient_expected_positive_cores"] >= 1).astype(int)

# ---------------------------------------------------------
# Merge with real biopsy outcomes
# ---------------------------------------------------------
merged = patient_agg.merge(real_df, on="patient_id", how="left")

# Ensure y_true and y_pred are integers
merged["positive"] = merged["positive"].fillna(0).astype(int)  # assume missing = 0
merged["predicted_positive"] = merged["predicted_positive"].astype(int)

y_true = merged["positive"].values
y_pred = merged["predicted_positive"].values

# ---------------------------------------------------------
# Compute confusion matrix and metrics
# ---------------------------------------------------------
y_true = merged["positive"].values
y_pred = merged["predicted_positive"].values

cm = confusion_matrix(y_true, y_pred)
accuracy = accuracy_score(y_true, y_pred)
precision = precision_score(y_true, y_pred)
recall = recall_score(y_true, y_pred)  # Sensitivity
f1 = f1_score(y_true, y_pred)

print("Confusion Matrix:")
print(cm)
print(f"\nAccuracy: {accuracy:.3f}")
print(f"Precision: {precision:.3f}")
print(f"Sensitivity (Recall): {recall:.3f}")
print(f"F1 Score: {f1:.3f}")

# ---------------------------------------------------------
# Save patient-level predictions
# ---------------------------------------------------------
merged.to_csv(PATIENT_PREDICTIONS_CSV, index=False)
print(f"\n✓ Patient-level predictions saved to {PATIENT_PREDICTIONS_CSV}")