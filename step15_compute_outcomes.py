from step1_load_data import load_full_dataset
from models.template import BiopsyTemplate
from models.simulation import BiopsySimulation


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
    sim.apply_error_to_needles(
        sigma_max_mm=3.0,
        random_seed=42
    )

    sim.intersect_needles_with_lesions()
    outcomes = sim.compute_needle_outcomes(step_mm=1.0)

    print("✓ Step 15 complete\n")

    for lesion in patient.lesions:
        o = lesion.simulation_results["outcomes"]
        print(f"Patient: {patient.id}")
        print(f"Lesion: {lesion.id}")
        print(f"  Hit flag: {o['hit_flag']}")
        print(f"  Positive length: {o['positive_length_mm']:.2f} mm")
        print(f"  Percentage positive: {100*o['percentage_positive']:.1f}%")
        print("-" * 40)


if __name__ == "__main__":
    main()