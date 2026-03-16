from step1_load_data import load_full_dataset
from models.template import BiopsyTemplate
from models.simulation import BiopsySimulation


def main():
    dataset = load_full_dataset()
    dataset.verify_all_geometry()

    patient = next(iter(dataset))

    prostate_center = patient.compute_prostate_centroid()

    template = BiopsyTemplate(
        grid_size=19,
        spacing=5.0,
        origin=prostate_center,   # 🔑 FIX
        direction=(0.0, 0.0, 1.0),
    )

    sim = BiopsySimulation(patient, template)
    valid_holes = sim.restrict_to_prostate()

    print(f"Total template holes: {len(template)}")
    print(f"Valid prostate-intersecting holes: {len(valid_holes)}")


if __name__ == "__main__":
    main()
