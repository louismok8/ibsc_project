import os
import sys
from pathlib import Path

import nibabel as nib
import numpy as np

PHASE_ROOT = Path(__file__).resolve().parent
if str(PHASE_ROOT) not in sys.path:
    sys.path.insert(0, str(PHASE_ROOT))

from data.paths import IMAGES_DIR, LABELS_DIR, ZONES_DIR
from models.patients import Patient
from models.dataset import Dataset


def load_nifti(path):
    img = nib.load(path)
    data = img.get_fdata(dtype=np.float32)
    return data, img.affine, img.header.get_zooms()[:3]


def load_full_dataset():
    patients = []

    for fname in sorted(os.listdir(LABELS_DIR)):
        if not fname.endswith(".nii.gz"):
            continue

        patient_id = fname.replace(".nii.gz", "")

        t2, t2_aff, spacing = load_nifti(
            os.path.join(IMAGES_DIR, f"{patient_id}_0000.nii.gz")
        )
        lesion, lesion_aff, _ = load_nifti(
            os.path.join(LABELS_DIR, f"{patient_id}.nii.gz")
        )
        prostate, prostate_aff, _ = load_nifti(
            os.path.join(ZONES_DIR, f"{patient_id}.nii.gz")
        )

        patients.append(
            Patient(
                patient_id=patient_id,
                t2=t2,
                lesion_mask=lesion,
                prostate_mask=prostate,
                affines={
                    "t2": t2_aff,
                    "lesion": lesion_aff,
                    "prostate": prostate_aff,
                },
                spacing=spacing,
            )
        )

    return Dataset(patients)


if __name__ == "__main__":
    dataset = load_full_dataset()
    print(f"Loaded {len(dataset)} patients.")