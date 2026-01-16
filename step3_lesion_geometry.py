from step1_load_data import load_full_dataset


def main():
    dataset = load_full_dataset()

    dataset.verify_all_geometry()

    dataset.extract_all_lesions()

    total_lesions = sum(len(p.lesions) for p in dataset)
    print(f"✓ Step 3 complete: extracted {total_lesions} lesions.")


if __name__ == "__main__":
    main()



"""
“Each lesion is now an independent object with its own binary mask, 
ready for geometric characterisation and simulation.”
"""