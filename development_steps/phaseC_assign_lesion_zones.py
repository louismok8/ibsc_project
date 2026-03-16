import pandas as pd

OUTPUT_FILE = "outputs/lesion_zones.csv"


def assign_lesion_zones(dataset):

    rows = []

    for patient in dataset:

        prostate_centroid = patient.compute_prostate_centroid()

        px, py, pz = prostate_centroid

        for lesion in patient.lesions:

            lx, ly, lz = lesion.centroid

            side = "left" if lx < px else "right"
            region = "posterior" if ly < py else "anterior"

            zone = f"{side}_{region}"

            rows.append({
                "patient_id": patient.id,
                "lesion_id": lesion.id,
                "centroid_x": lx,
                "centroid_y": ly,
                "centroid_z": lz,
                "zone": zone
            })

    df = pd.DataFrame(rows)

    df.to_csv(OUTPUT_FILE, index=False)

    print(f"✓ Saved lesion zones → {OUTPUT_FILE}")


    