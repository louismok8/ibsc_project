from step1_load_data import load_full_dataset

def main():
    dataset = load_full_dataset()

    dataset.verify_all_geometry()
    dataset.extract_all_lesions()
    dataset.compute_all_lesion_volumes()

    volumes = [
        lesion.volume
        for patient in dataset
        for lesion in patient.lesions
    ]

    print(f"✓ Step 4 complete: computed volumes for {len(volumes)} lesions")
    print(f"  Example volumes (mm³): {volumes[:5]}")

if __name__ == "__main__":
    main()
