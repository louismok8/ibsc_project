from step1_load_data import load_full_dataset
from models.template import BiopsyTemplate
from models.simulation import BiopsySimulation
import numpy as np


def main():

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

    print("Running Monte Carlo simulation...\n")

    results = sim.run_monte_carlo(
        n_simulations=1000,
        n_cores=5,
        sigma_max_mm=3.0,
        step_mm=1.0,
    )

    print("✓ Phase 8 Monte Carlo complete\n")

    for lesion_id, stats in results.items():

        print(f"Lesion {lesion_id}")
        print(f"  Hit probability: {stats['hit_probability']:.3f}")
        print(f"  Mean % positive: {100*stats['mean_percentage_positive']:.1f}%")
        print(f"  Mean # positive cores: "
              f"{np.mean(stats['distribution_positive_core_counts']):.2f}")
        print("-" * 40)


if __name__ == "__main__":
    main()