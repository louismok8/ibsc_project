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
