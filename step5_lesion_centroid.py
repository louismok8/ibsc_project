from step1_load_data import load_full_dataset

def main():
    dataset = load_full_dataset()

    dataset.verify_all_geometry()
    dataset.extract_all_lesions()
    dataset.compute_all_lesion_volumes()
    dataset.compute_all_lesion_centroids()

    examples = [(lesion.patient_id, lesion.id, lesion.centroid)
        for patient in dataset
        for lesion in patient.lesions
    ][:5]

    print("✓ Step 5 complete: computed lesion centroids")
    for pid, lid, c in examples:
        print(f"  Patient {pid}, Lesion {lid}: centroid = {c} (mm)")

if __name__ == "__main__":
    main()
