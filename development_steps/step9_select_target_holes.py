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

    print(f"Patient ID: {patient.id}")
    print(f"Number of lesions: {len(patient.lesions)}")

    prostate_center = patient.compute_prostate_centroid()

    template = BiopsyTemplate(
        grid_size=19,
        spacing=5.0,
        origin=prostate_center,
        direction=(0.0, 0.0, 1.0),
    )

    sim = BiopsySimulation(patient, template)
    sim.restrict_to_prostate()

    selected = sim.select_target_holes()

    print("✓ Step 9 complete")
    print(f"Lesions: {len(selected)}")

    for lesion in patient.lesions[:5]:
        d = lesion.simulation_results["hole_distance_mm"]
        print(
            f"  Lesion {lesion.id}: "
            f"distance to selected hole = {d:.2f} mm"
        )


if __name__ == "__main__":
    main()