import numpy as np


class Lesion:
    """
    Represents a single MRI-defined lesion.

    This object owns:
    - a binary lesion mask
    - a lesion index (label)
    - a reference to the parent patient

    Geometry (volume, centroid) is computed in later steps.
    """

    def __init__(self, lesion_id, binary_mask, patient_id):
        """
        Parameters
        ----------
        lesion_id : int
            Lesion label index (ℓ >= 1)

        binary_mask : np.ndarray
            Boolean or {0,1} array indicating lesion voxels.

        patient_id : str
            Identifier of the parent patient (for traceability).
        """
        self.id = lesion_id
        self.mask = binary_mask.astype(np.uint8)
        self.patient_id = patient_id

        # Placeholders for later steps
        self.volume = None
        self.centroid = None
        self.simulation_results = {}

    
    def compute_volume(self, spacing):
        """
        Compute lesion volume in mm^3.
        """
        voxel_volume = float(spacing[0] * spacing[1] * spacing[2])
        voxel_count = int(self.mask.sum())

        self.volume = float(voxel_count * voxel_volume)
        return self.volume
    
    import numpy as np


    def compute_centroid(self, affine):
        """
        Compute lesion centroid in world coordinates (mm).

        Parameters
        ----------
        affine : np.ndarray, shape (4, 4)
            Voxel-to-world affine matrix.
        """
        coords = np.argwhere(self.mask == 1)

        if coords.size == 0:
            raise ValueError(
                f"Lesion {self.id} in patient {self.patient_id} has empty mask"
            )

        centroid_voxel = coords.mean(axis=0)

        # Convert to homogeneous coordinates
        centroid_h = np.append(centroid_voxel, 1.0)

        centroid_world = affine @ centroid_h

        # Store only (x, y, z)
        self.centroid = centroid_world[:3].astype(float)
        return self.centroid



