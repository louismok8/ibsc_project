from step1_load_data import load_full_dataset
import matplotlib.pyplot as plt
import numpy as np

def main():
    dataset = load_full_dataset()

    dataset.verify_all_geometry()
    dataset.extract_all_lesions()
    dataset.compute_all_lesion_volumes()
    dataset.compute_all_lesion_centroids()

    volumes = np.array([
        lesion.volume
        for patient in dataset
        for lesion in patient.lesions
    ])

    print(f"Total lesions: {len(volumes)}")
    print(f"Min volume: {volumes.min():.1f} mm³")
    print(f"Median volume: {np.median(volumes):.1f} mm³")
    print(f"Max volume: {volumes.max():.1f} mm³\n")

    q1, q2 = np.percentile(volumes, [33.3, 66.7])

    print(f"Suggested size bins (quantile-based):")
    print(f"  Small   : < {q1:.1f} mm³")
    print(f"  Medium  : {q1:.1f} – {q2:.1f} mm³")
    print(f"  Large   : > {q2:.1f} mm³")

    # Plot global histogram
    plt.figure()
    plt.hist(volumes, bins=20)
    plt.xlabel("Lesion volume (mm³)")
    plt.ylabel("Number of lesions")
    plt.title("Global lesion volume distribution")
    plt.show()


if __name__ == "__main__":
    main()
