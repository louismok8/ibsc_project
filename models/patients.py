import numpy as np
from models.lesion import Lesion



class Patient:
    """
    Represents a single patient and all associated MRI-space data.

    This object owns:
    - T2-weighted MRI
    - lesion label mask
    - prostate binary mask
    - spatial metadata (affines, spacing)

    All geometry checks and lesion-level operations belong here.
    """

    def __init__(self, patient_id, t2, lesion_mask, prostate_mask, affines, spacing):
        self.id = patient_id

        self.t2 = t2
        self.lesion_mask = lesion_mask
        self.prostate_mask = prostate_mask

        self.affines = affines
        self.spacing = spacing

        # Filled later in Step 3
        self.lesions = []

    # ============================
    # STEP 2: Geometry verification
    # ============================

    def verify_geometry(self):
        """
        Run all spatial consistency checks for this patient.
        """
        self._check_shapes()
        self._check_affines()
        self._check_masks()

    def _check_shapes(self):
        if self.t2.shape != self.lesion_mask.shape:
            raise ValueError(f"[{self.id}] T2 / lesion shape mismatch")
        if self.t2.shape != self.prostate_mask.shape:
            raise ValueError(f"[{self.id}] T2 / prostate shape mismatch")

    def _check_affines(self, tol=1e-5):
        ref = self.affines["t2"]
        for name, affine in self.affines.items():
            if not np.allclose(ref, affine, atol=tol):
                raise ValueError(f"[{self.id}] Affine mismatch: {name}")

    def _check_masks(self):
        # Lesion mask must be integer-labeled
        lesion_vals = np.unique(self.lesion_mask)
        if not np.all(np.equal(lesion_vals, lesion_vals.astype(int))):
            raise ValueError(f"[{self.id}] Non-integer lesion labels")

        # Prostate mask must be binary
        prostate_vals = np.unique(self.prostate_mask)
        if not set(prostate_vals).issubset({0, 1}):
            raise ValueError(f"[{self.id}] Prostate mask not binary")

    # ============================
    # STEP 3: Lesion extraction
    # ============================

    def extract_lesions(self):
        """
        Extract individual lesion masks from the multi-label lesion mask.

        For lesion indices ℓ = 1 ... L:
            lesion_ℓ = (lesion_mask == ℓ)

        Populates:
            self.lesions : list[Lesion]

        Returns
        -------
        list[Lesion]
            List of extracted lesion objects.
        """

        lesion_labels = np.unique(self.lesion_mask)

        # Exclude background label (0)
        lesion_labels = lesion_labels[lesion_labels != 0]

        lesions = []

        for label in lesion_labels:
            binary_mask = (self.lesion_mask == label)

            lesion = Lesion(
                lesion_id=int(label),
                binary_mask=binary_mask,
                patient_id=self.id
            )

            lesions.append(lesion)

        self.lesions = lesions
        return lesions