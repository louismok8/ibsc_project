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
        self.lesion_mask = lesion_mask # one 3D array of lesions with labels of 0(background), 1(tumour 1), 2+(tumour 2+)
        self.prostate_mask = prostate_mask

        self.affines = affines
        self.spacing = spacing

        # Filled later in Step 3
        self.lesions = []


    # ============================
    # STEP 2: Geometry verification
    #   Ensuring that all arrays are in the same voxel grid, same physical space, all masks obey same rules
    # ============================
    def verify_geometry(self): # Public orchestration method to be called by upstream Dataset method or patient.verify_geometry()
        self._check_shapes()
        self._check_affines()
        self._check_masks()

    def _check_shapes(self):
        # t2 shape must equal lesion mask shape (including background so should be same)
        if self.t2.shape != self.lesion_mask.shape:
            raise ValueError(f"[{self.id}] T2 / lesion shape mismatch")
        # t2 shape must equal prostate mask shape (including background so should be same)
        if self.t2.shape != self.prostate_mask.shape:
            raise ValueError(f"[{self.id}] T2 / prostate shape mismatch")

    def _check_affines(self, tol=1e-5): # checks for actual physical spatial coordination of voxels across all 3 masks: t2, lesion, prostate
        ref = self.affines["t2"] # select t2 as the reference coordinate space
        for name, affine in self.affines.items():
            if not np.allclose(ref, affine, atol=tol): # checks reference affine to actual affine of "t2": t2_affine, "lesion": lesion_affine, "prostate": prostate_affine with specified tolerance
                raise ValueError(f"[{self.id}] Affine mismatch: {name}")

    def _check_masks(self): 
        # Lesion mask must be integer-labeled
        lesion_vals = np.unique(self.lesion_mask) # Extract all unique lesion mask labels (0,1,2+)
        if not np.all(np.equal(lesion_vals, lesion_vals.astype(int))): # Ensure all lesion mask labels are integers, as floats can cause errors
            raise ValueError(f"[{self.id}] Non-integer lesion labels")

        # Prostate mask must be binary
        prostate_vals = np.unique(self.prostate_mask) # Extract all unique prostate mask labels (0,1)
        if not set(prostate_vals).issubset({0, 1}): # Ensure all lesion mask labels are integers, as floats can cause errors
            raise ValueError(f"[{self.id}] Non-integer lesion labels")


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
        
        # Create 1D array of each lesion's labels
        lesion_labels = np.unique(self.lesion_mask)

        # Exclude background label (0) via a boolean mask on the array lesion_labels
        lesion_labels = lesion_labels[lesion_labels != 0]

        # Stores Lesion class objects/instances eg: lesions = [Lesions(id=1),Lesions(id=2)]
        lesions = []

        for label in lesion_labels:
            # compares every voxel to label (1,2) in lesion labels, returns Boolean array of same shape as lesion
            binary_mask = (self.lesion_mask == label)

            lesion = Lesion(
                lesion_id=int(label),
                binary_mask=binary_mask,
                patient_id=self.id
            )

            lesions.append(lesion)
        
        # Constructor assignment
        self.lesions = lesions
        return lesions
    

    def compute_lesion_volumes(self):
        """
        Compute volumes for all lesions in this patient.
        """
        for lesion in self.lesions:
            lesion.compute_volume(self.spacing)


    def compute_lesion_centroids(self):
        """
        Compute world-coordinate centroids for all lesions.
        """
        affine = self.affines["t2"]

        for lesion in self.lesions:
            lesion.compute_centroid(affine)
