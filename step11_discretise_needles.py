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

    sampled = sim.discretise_needles(step_mm=1.0)

    print("✓ Step 11 complete\n")

    for lesion in patient.lesions[:3]:
        pts = lesion.simulation_results["needle_points"]
        print(f"Lesion {lesion.id}")
        print(f"  Number of sampled points: {len(pts)}")
        print(f"  First point: {pts[0]}")
        print(f"  Last point:  {pts[-1]}")
        print()

    # Testing Error
    for lesion in patient.lesions[:1]:
        errors = sim.generate_error_field(lesion, sigma_max_mm=3.0)

        print("Error shape:", errors.shape)
        print("First point error:", errors[0])
        print("Tip error:", errors[-1])


if __name__ == "__main__":
    main()