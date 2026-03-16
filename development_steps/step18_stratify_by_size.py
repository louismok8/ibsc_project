import matplotlib.pyplot as plt
import numpy as np
from step1_load_data import load_full_dataset
from models.template import BiopsyTemplate
from models.simulation import BiopsySimulation
from scipy import stats

def compute_size_bins(volumes):
    q1, q2 = np.percentile(volumes, [33.3, 66.7])
    return q1, q2

def main():

    dataset = load_full_dataset()
    dataset.verify_all_geometry()
    dataset.extract_all_lesions()
    dataset.compute_all_lesion_volumes()
    dataset.compute_all_lesion_centroids()

    patient = next(iter(dataset))

    # Extract all lesion volumes
    volumes = np.array([lesion.volume for lesion in patient.lesions])
    q1, q2 = compute_size_bins(volumes)

    # Assign lesions to size bins
    size_groups = {"Small": [], "Medium": [], "Large": []}
    for lesion in patient.lesions:
        if lesion.volume < q1:
            size_groups["Small"].append(lesion)
        elif lesion.volume < q2:
            size_groups["Medium"].append(lesion)
        else:
            size_groups["Large"].append(lesion)

    # Run simulation (Monte Carlo)
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

    results = sim.run_monte_carlo(
        n_simulations=1000,
        n_cores=5,
        sigma_max_mm=3.0,
        step_mm=1.0,
    )

    # Aggregate statistics per size bin
    hit_probs_dict = {}
    expected_cores_dict = {}
    percent_positive_dict = {}

    print("✓ Step 18 complete: Stratified lesion statistics\n")

    for size_label, lesions in size_groups.items():
        if not lesions:
            print(f"{size_label}: No lesions in this bin")
            hit_probs_dict[size_label] = []
            expected_cores_dict[size_label] = []
            percent_positive_dict[size_label] = []
            continue

        hit_probs = [results[lesion.id]["hit_probability"] for lesion in lesions]
        expected_positives = [
            np.mean(results[lesion.id]["distribution_positive_core_counts"])
            for lesion in lesions
        ]
        mean_percent_positive = [
            results[lesion.id]["mean_percentage_positive"]
            for lesion in lesions
        ]

        hit_probs_dict[size_label] = hit_probs
        expected_cores_dict[size_label] = expected_positives
        percent_positive_dict[size_label] = mean_percent_positive

        print(f"Size group: {size_label}")
        print(f"  # Lesions: {len(lesions)}")
        print(f"  Mean hit probability: {np.mean(hit_probs):.3f}")
        print(f"  Mean expected positive cores: {np.mean(expected_positives):.2f}")
        print(f"  Mean % positive core: {100*np.mean(mean_percent_positive):.1f}%")
        print("-"*50)

        # -----------------------------
        # Statistical comparison example: Small vs Large
        # -----------------------------
        if (len(hit_probs_dict["Small"]) > 1 and 
            len(hit_probs_dict["Large"]) > 1):

            stat, p = stats.mannwhitneyu(
                hit_probs_dict["Small"],
                hit_probs_dict["Large"],
                alternative="two-sided"
            )

            print("\nStatistical test (Small vs Large hit probability)")
            print(f"  Mann-Whitney U p-value: {p:.4f}")


    # ===== Visualization =====
    labels = ["Small", "Medium", "Large"]

    def plot_metric(metric_dict, ylabel, title):
        means = [np.mean(metric_dict[label]) if metric_dict[label] else 0 for label in labels]
        plt.figure()
        plt.bar(labels, means, color=["skyblue", "orange", "green"])
        plt.ylabel(ylabel)
        plt.title(title)
        plt.ylim(0, 1.05 if "probability" in ylabel.lower() else None)
        plt.show()

    plot_metric(hit_probs_dict, "Mean Hit Probability", "Hit Probability by Lesion Size")
    plot_metric(expected_cores_dict, "Mean Expected Positive Cores", "Expected Positive Cores by Lesion Size")
    plot_metric(percent_positive_dict, "Mean % Positive Core / 100", "% Positive Core by Lesion Size")

if __name__ == "__main__":
    main()