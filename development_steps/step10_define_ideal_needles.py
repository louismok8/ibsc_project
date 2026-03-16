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

    print("✓ Step 10 complete\n")

    for lesion in patient.lesions:
        needle = lesion.simulation_results["ideal_needle"]

        print(f"Lesion {lesion.id}")
        print(f"  Start:  {needle['start']}")
        print(f"  End:    {needle['end']}")
        print(f"  Length: {needle['length_mm']} mm\n")


if __name__ == "__main__":
    main()