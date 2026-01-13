from step1_load_data import load_full_dataset


def main():
    dataset = load_full_dataset()

    # Always validate geometry before extracting lesions
    dataset.verify_all_geometry()

    # Step 3: extract lesions
    dataset.extract_all_lesions()

    # Optional sanity output
    total_lesions = sum(len(p.lesions) for p in dataset)
    print(f"✓ Step 3 complete: extracted {total_lesions} lesions.")


if __name__ == "__main__":
    main()



"""
“Each lesion is now an independent object with its own binary mask, 
ready for geometric characterisation and simulation.”
"""