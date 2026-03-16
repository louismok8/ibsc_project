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

    realised = sim.apply_error_to_needles(
        sigma_max_mm=3.0,
        random_seed=42
    )

    print("✓ Step 13 complete\n")

    for lesion in patient.lesions:
        ideal = lesion.simulation_results["needle_points"]
        real = lesion.simulation_results["realised_needle_points"]

        print(f"Patient: {patient.id}")
        print(f"Lesion: {lesion.id}")

        print("  First ideal point:    ", ideal[0])
        print("  First realised point: ", real[0])
        print()
        print("  Tip ideal:            ", ideal[-1])
        print("  Tip realised:         ", real[-1])
        print("-" * 50)

    tip_displacement = np.linalg.norm(real[-1] - ideal[-1])
    print(f"  Tip displacement magnitude: {tip_displacement:.2f} mm")


if __name__ == "__main__":
    main()