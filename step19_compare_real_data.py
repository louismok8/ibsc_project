import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from step1_load_data import load_full_dataset
from models.template import BiopsyTemplate
from models.simulation import BiopsySimulation
from scipy import stats

# -----------------------------
# 1. Load dataset and run Monte Carlo if not saved
# -----------------------------
dataset = load_full_dataset()
dataset.verify_all_geometry()
dataset.extract_all_lesions()
dataset.compute_all_lesion_volumes()
dataset.compute_all_lesion_centroids()

patient = next(iter(dataset))

template = BiopsyTemplate(
    grid_size=19,
    spacing=5.0,
    origin=patient.compute_prostate_centroid(),
    direction=(0.0, 0.0, 1.0),
)

sim = BiopsySimulation(patient, template)
sim.restrict_to_prostate()
sim.select_target_holes()
sim.define_ideal_needles(core_length_mm=20.0)
sim.discretise_needles(step_mm=1.0)

# Run Monte Carlo (or load precomputed results)
results = sim.run_monte_carlo(
    n_simulations=1000,
    n_cores=5,
    sigma_max_mm=3.0,
    step_mm=1.0,
)

# -----------------------------
# 2. Load real biopsy data
# -----------------------------
# Expecting a CSV or Excel with columns: lesion_id, positive_cores, percent_positive
real_data_path = "data/real_biopsy_stats.csv"
real_df = pd.read_csv(real_data_path)

# -----------------------------
# 3. Aggregate simulated metrics per lesion
# -----------------------------
sim_data = []
for lesion_id, stats in results.items():
    sim_data.append({
        "lesion_id": lesion_id,
        "mean_positive_cores": np.mean(stats["distribution_positive_core_counts"]),
        "mean_percent_positive": stats["mean_percentage_positive"] * 100,
    })

sim_df = pd.DataFrame(sim_data)

# -----------------------------
# 4. Merge real and simulated data
# -----------------------------
merged_df = pd.merge(real_df, sim_df, on="lesion_id", how="inner", suffixes=("_real", "_sim"))

# Kolmogorov-Smirnov test for distribution similarity
ks_stat, ks_p = stats.ks_2samp(
    merged_df["percent_positive_real"],
    merged_df["mean_percent_positive"]
)

print(f"\nKS test p-value (% positive): {ks_p:.4f}")

# Mann-Whitney for core counts
u_stat, u_p = stats.mannwhitneyu(
    merged_df["positive_cores_real"],
    merged_df["mean_positive_cores"],
    alternative="two-sided"
)

print(f"Mann-Whitney p-value (# cores): {u_p:.4f}")

# Correlation
corr, corr_p = stats.pearsonr(
    merged_df["percent_positive_real"],
    merged_df["mean_percent_positive"]
)

print(f"Pearson correlation: r={corr:.3f}, p={corr_p:.4f}")

print("✓ Comparison table (first 5 lesions):")
print(merged_df.head())

# -----------------------------
# 5. Visual comparison
# -----------------------------
def plot_comparison(real_col, sim_col, ylabel, title):
    plt.figure(figsize=(6,4))
    plt.hist(merged_df[real_col], bins=10, alpha=0.5, label="Observed")
    plt.hist(merged_df[sim_col], bins=10, alpha=0.5, label="Simulated")
    plt.xlabel(ylabel)
    plt.ylabel("Number of lesions")
    plt.title(title)
    plt.legend()
    plt.show()

plot_comparison("percent_positive_real", "mean_percent_positive", "% Positive Core", "Observed vs Simulated % Positive Core")
plot_comparison("positive_cores_real", "mean_positive_cores", "Number of Positive Cores", "Observed vs Simulated # Positive Cores")

# -----------------------------
# 6. Summary statistics
# -----------------------------
for col_real, col_sim in [("percent_positive_real", "mean_percent_positive"),
                          ("positive_cores_real", "mean_positive_cores")]:
    diff = merged_df[col_sim] - merged_df[col_real]
    print(f"{col_real} vs {col_sim}:")
    print(f"  Mean difference: {diff.mean():.2f}")
    print(f"  Std deviation:  {diff.std():.2f}")
    print(f"  Min / Max diff: {diff.min():.2f} / {diff.max():.2f}")
    print("-"*40)

print("✓ Step 19 complete: simulated vs observed comparison done.")

# Get the file for data/real_biopsy_stats.csv