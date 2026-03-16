"""
Master pipeline for MRI-targeted prostate biopsy simulation.

This script replaces step1–step18 development scripts and runs
the full simulation workflow in a single reproducible pipeline.

Pipeline stages
---------------
1. Load dataset
2. Verify geometry
3. Extract lesions
4. Compute lesion geometry (volume + centroid)
5. Build biopsy template
6. Restrict template to prostate
7. Select target holes
8. Define ideal biopsy needles
9. Discretise needles
10. Run Monte Carlo biopsy simulation
11. Aggregate and print results
"""

import os
import sys
import numpy as np
import nibabel as nib
import pandas as pd

from sklearn.metrics import (
    confusion_matrix,
    accuracy_score,
    precision_score,
    recall_score,
    f1_score
)

# ---------------------------------------------------------
# Ensure project root is on Python path
# ---------------------------------------------------------

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(PROJECT_ROOT)

# ---------------------------------------------------------
# Imports from project modules
# ---------------------------------------------------------

from data.paths import IMAGES_DIR, LABELS_DIR, ZONES_DIR
from models.dataset import Dataset
from models.patients import Patient
from models.template import BiopsyTemplate
from models.simulation import BiopsySimulation


# ---------------------------------------------------------
# Utility: Load NIfTI file
# ---------------------------------------------------------

def load_nifti(path):
    """
    Load NIfTI image and return data, affine, spacing.
    """
    img = nib.load(path)
    data = img.get_fdata(dtype=np.float32)
    affine = img.affine
    spacing = img.header.get_zooms()[:3]

    return data, affine, spacing


# ---------------------------------------------------------
# Step 1: Load full dataset
# ---------------------------------------------------------

def load_full_dataset():
    """
    Load all patients from the dataset directory structure.
    """

    patients = []

    for fname in sorted(os.listdir(LABELS_DIR)):

        if not fname.endswith(".nii.gz"):
            continue

        patient_id = fname.replace(".nii.gz", "")

        print(f"Loading patient {patient_id}")

        t2, t2_aff, spacing = load_nifti(
            os.path.join(IMAGES_DIR, f"{patient_id}_0000.nii.gz")
        )

        lesion, lesion_aff, _ = load_nifti(
            os.path.join(LABELS_DIR, f"{patient_id}.nii.gz")
        )

        prostate, prostate_aff, _ = load_nifti(
            os.path.join(ZONES_DIR, f"{patient_id}.nii.gz")
        )

        patient = Patient(
            patient_id=patient_id,
            t2=t2,
            lesion_mask=lesion,
            prostate_mask=prostate,
            affines={
                "t2": t2_aff,
                "lesion": lesion_aff,
                "prostate": prostate_aff,
            },
            spacing=spacing,
        )

        patients.append(patient)

    dataset = Dataset(patients)

    print(f"\n✓ Loaded dataset with {len(dataset)} patients\n")

    return dataset


# ---------------------------------------------------------
# Main simulation pipeline
# ---------------------------------------------------------

def run_simulation_pipeline(
    n_simulations=1000,
    n_cores=5,
    sigma_max_mm=5.0,
    step_mm=1.0,
):

    # -----------------------------------------------------
    # Load dataset
    # -----------------------------------------------------

    dataset = load_full_dataset()

    # -----------------------------------------------------
    # Geometry validation
    # -----------------------------------------------------

    print("Verifying geometry...")
    dataset.verify_all_geometry()
    print("✓ Geometry valid\n")

    # -----------------------------------------------------
    # Lesion extraction
    # -----------------------------------------------------

    print("Extracting lesions...")
    dataset.extract_all_lesions()

    total_lesions = sum(len(p.lesions) for p in dataset)
    print(f"✓ Extracted {total_lesions} lesions\n")

    # -----------------------------------------------------
    # Lesion geometry
    # -----------------------------------------------------

    print("Computing lesion geometry...")
    dataset.compute_all_lesion_volumes()
    dataset.compute_all_lesion_centroids()
    print("✓ Lesion geometry computed\n")

    # -----------------------------------------------------
    # Run simulation for each patient
    # -----------------------------------------------------

    all_results = []

    for p_idx, patient in enumerate(dataset, start=1):

        print(f"Patient {p_idx}/{len(dataset)}: {patient.id}")
        print("=" * 60)
        print(f"Running simulation for patient {patient.id}")
        print("=" * 60)

        # Build biopsy template centred on prostate
        prostate_center = patient.compute_prostate_centroid()

        template = BiopsyTemplate(
            grid_size=19,
            spacing=5.0,
            origin=prostate_center,
            direction=(0.0, 0.0, 1.0),
        )

        # Create simulation object
        sim = BiopsySimulation(patient, template)

        # Restrict template holes
        valid_holes = sim.restrict_to_prostate()
        print(f"Valid template holes: {len(valid_holes)}")

        # Target lesions
        sim.select_target_holes()

        # Construct ideal needles
        sim.define_ideal_needles(core_length_mm=20.0)

        # Discretise needles
        sim.discretise_needles(step_mm=step_mm)

        # -------------------------------------------------
        # Monte Carlo simulation
        # -------------------------------------------------

        print("Running Monte Carlo simulation...")

        results = sim.run_monte_carlo(
            n_simulations=n_simulations,
            n_cores=n_cores,
            sigma_max_mm=sigma_max_mm,
            step_mm=step_mm,
        )

        print("✓ Monte Carlo complete\n")

        # Store results
        all_results.append({
            "patient_id": patient.id,
            "lesion_results": results
        })

        # -------------------------------------------------
        # Print summary
        # -------------------------------------------------

        for lesion_id, stats in results.items():

            mean_cores = np.mean(stats["distribution_positive_core_counts"])

            print(f"Lesion {lesion_id}")
            print(f"  Hit probability: {stats['hit_probability']:.3f}")
            print(f"  Mean % positive: {100*stats['mean_percentage_positive']:.1f}%")
            print(f"  Mean positive cores: {mean_cores:.2f}")
            print("-" * 40)

        print()

    print("\n✓ Simulation pipeline complete")

    return all_results



# ---------------------------------------------------------
# Evaluation pipeline (Phase A–C)
# ---------------------------------------------------------

def run_full_pipeline():

    print("\n==============================")
    print("Running FULL simulation + evaluation pipeline")
    print("==============================\n")

    results = run_simulation_pipeline()

    dataset = load_full_dataset()

    dataset.verify_all_geometry()
    dataset.extract_all_lesions()
    dataset.compute_all_lesion_volumes()
    dataset.compute_all_lesion_centroids()

    # -----------------------------------------------------
    # Phase A — Save lesion simulation results
    
    # OUTPUT: outputs/simulation_results.py : per-lesion simulation results
    # -----------------------------------------------------

    print("\nPhase A — Saving simulation results")

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

    sim_df = pd.DataFrame(rows)

    os.makedirs(os.path.join(PROJECT_ROOT, "outputs"), exist_ok=True)

    sim_path = os.path.join(PROJECT_ROOT, "outputs/simulation_results.csv")
    sim_df.to_csv(sim_path, index=False)

    print(f"✓ Saved simulation results → {sim_path}")

    # -----------------------------------------------------
    # Phase B — Patient-level accuracy

    # OUTPUT: outputs/patient_predictions.csv : Aggregates lesion-level results to patient-level predictions
    # max(hit_probability) threshold >= 0.5 to determine if patient is positive
    # Compare with real data
    # Compute confusion matrix
    # -----------------------------------------------------

    print("\nPhase B — Patient-level detection accuracy")

    real_path = os.path.join(PROJECT_ROOT, "data/real_patient_results.csv")

    real_df = pd.read_csv(real_path, dtype={"patient_id": str})
    sim_df["patient_id"] = sim_df["patient_id"].astype(str).str.zfill(3)
    real_df["patient_id"] = real_df["patient_id"].astype(str).str.zfill(3)

    patient_agg = sim_df.groupby("patient_id").agg(
        patient_hit_probability=("hit_probability", "max"),
        patient_expected_positive_cores=("expected_positive_cores", "sum")
    ).reset_index()

    patient_agg["predicted_positive"] = (
        patient_agg["patient_hit_probability"] > 0.5
    ).astype(int)

    merged = patient_agg.merge(real_df, on="patient_id", how="left")

    merged["positive"] = merged["positive"].fillna(0).astype(int)

    y_true = merged["positive"].values
    y_pred = merged["predicted_positive"].values

    cm = confusion_matrix(y_true, y_pred)

    print("\nConfusion Matrix")
    print(cm)

    print("\nMetrics")
    print(f"Accuracy:  {accuracy_score(y_true,y_pred):.3f}")
    print(f"Precision: {precision_score(y_true,y_pred):.3f}")
    print(f"Recall:    {recall_score(y_true,y_pred):.3f}")
    print(f"F1 Score:  {f1_score(y_true,y_pred):.3f}")

    pred_path = os.path.join(PROJECT_ROOT, "outputs/patient_predictions.csv")
    merged.to_csv(pred_path, index=False)

    print(f"\n✓ Saved patient predictions → {pred_path}")

    # -----------------------------------------------------
    # Phase C — Lesion zones

    # OUTPUT: outputs/lesion_zones.csv : assigns a lesion to a zone within prostate based on centroid
    # -----------------------------------------------------

    def assign_lesion_zones(dataset):

        rows = []

        for patient in dataset:

            prostate_centroid = patient.compute_prostate_centroid()
            px, py, pz = prostate_centroid

            for lesion in patient.lesions:

                lx, ly, lz = lesion.centroid

                side = "left" if lx < px else "right"
                region = "posterior" if ly < py else "anterior"

                zone = f"{side}_{region}"

                rows.append({
                    "patient_id": patient.id,
                    "lesion_id": lesion.id,
                    "centroid_x": lx,
                    "centroid_y": ly,
                    "centroid_z": lz,
                    "zone": zone
                })

        df = pd.DataFrame(rows)

        output_file = os.path.join(PROJECT_ROOT, "outputs/lesion_zones.csv")

        df.to_csv(output_file, index=False)

        print(f"✓ Saved lesion zones → {output_file}")

    print("\nPhase C — Assigning lesion zones")

    assign_lesion_zones(dataset)

    print("\n✓ Evaluation pipeline complete")


    # -----------------------------------------------------
    # Phase D — Zone-level simulation predictions

    # OUTPUT: outputs/zone_level_predictions.csv : Merges lesion-level zone assignments (from Phase C) with simulation hit probabilities (Phase A)
    # -----------------------------------------------------

    print("\nPhase D — Zone-level simulation predictions")

    # Merge lesion zones with lesion simulation results
    zone_df = pd.read_csv(os.path.join(PROJECT_ROOT, "outputs/lesion_zones.csv"))
    sim_df = pd.read_csv(os.path.join(PROJECT_ROOT, "outputs/simulation_results.csv"))

    # Ensure correct types
    zone_df["patient_id"] = zone_df["patient_id"].astype(str).str.zfill(3)
    sim_df["patient_id"] = sim_df["patient_id"].astype(str).str.zfill(3)

    # Merge to get hit probabilities per lesion with their zones
    merged_zone = zone_df.merge(
        sim_df[["patient_id", "lesion_id", "hit_probability"]],
        on=["patient_id", "lesion_id"],
        how="left"
    )

    # For each lesion, define a "lesion_hit" if hit_probability > threshold (e.g., 0.5)
    merged_zone["lesion_hit"] = (merged_zone["hit_probability"] > 0.5).astype(int)

    # Aggregate to zone level: if any lesion in zone hit → zone predicted positive
    zone_preds = merged_zone.groupby(["patient_id", "zone"])["lesion_hit"].max().reset_index()
    zone_preds.rename(columns={"lesion_hit": "predicted_positive"}, inplace=True)

    # Save
    # Only include rows for zones that contain at least one lesion (zone-level prediction)
    zone_path = os.path.join(PROJECT_ROOT, "outputs/zone_level_predictions.csv")
    zone_preds.to_csv(zone_path, index=False)

    print(f"✓ Zone-level predictions saved → {zone_path}")


    # -----------------------------------------------------
    # Phase E — Zone-level accuracy metrics

    # OUTPUT: outputs/zone_accuracy_metrics.csv : Compares predicted zone-level results (Phase D) with real biopsy zone-level results
    # Confirms how well the simulation predicts positives/negatives zone-wise across all patients
    # -----------------------------------------------------

    print("\nPhase E — Zone-level accuracy metrics")

    # Load the real biopsy zone-level results
    real_zone_path = os.path.join(PROJECT_ROOT, "data/real_zone_results.csv")
    real_zone_df = pd.read_csv(real_zone_path, dtype={"patient": str, "zone": str})

    # Ensure consistent formatting
    zone_df = zone_preds.copy()
    zone_df["patient_id"] = zone_df["patient_id"].astype(str).str.zfill(3)
    real_zone_df["patient"] = real_zone_df["patient"].astype(str).str.zfill(3)

    # Merge simulation predictions with real data
    merged_zone = pd.merge(
        real_zone_df,
        zone_df,
        left_on=["patient", "zone"],
        right_on=["patient_id", "zone"],
        how="left"
    )

    # Fill missing predictions with 0 (zones with no lesions assumed negative)
    merged_zone["predicted_positive"] = merged_zone["predicted_positive"].fillna(0).astype(int)
    merged_zone["biopsy_positive"] = merged_zone["biopsy_positive"].astype(int)

    # Compute per-zone metrics
    zone_metrics = []

    for zone_name in merged_zone["zone"].unique():
        zone_data = merged_zone[merged_zone["zone"] == zone_name]

        y_true = zone_data["biopsy_positive"].values
        y_pred = zone_data["predicted_positive"].values

        TP = np.sum((y_true == 1) & (y_pred == 1))
        TN = np.sum((y_true == 0) & (y_pred == 0))
        FP = np.sum((y_true == 0) & (y_pred == 1))
        FN = np.sum((y_true == 1) & (y_pred == 0))

        sensitivity = TP / (TP + FN) if (TP + FN) > 0 else np.nan
        specificity = TN / (TN + FP) if (TN + FP) > 0 else np.nan
        accuracy = (TP + TN) / (TP + TN + FP + FN)

        zone_metrics.append({
            "zone": zone_name,
            "TP": TP,
            "TN": TN,
            "FP": FP,
            "FN": FN,
            "sensitivity": sensitivity,
            "specificity": specificity,
            "accuracy": accuracy
        })

    zone_metrics_df = pd.DataFrame(zone_metrics)

    # Save to outputs
    zone_metrics_file = os.path.join(PROJECT_ROOT, "outputs/zone_accuracy_metrics.csv")
    zone_metrics_df.to_csv(zone_metrics_file, index=False)

    print(f"✓ Zone-level metrics saved → {zone_metrics_file}")
    print("\nZone-level metrics summary:")
    print(zone_metrics_df)






# ---------------------------------------------------------
# Entry point
# ---------------------------------------------------------

if __name__ == "__main__":

    run_full_pipeline()