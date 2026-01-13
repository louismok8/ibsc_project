"""
Dataset-level container and orchestration.
"""

class Dataset:
    """
    Represents the entire patient cohort.

    Responsibilities:
    - Store Patient objects
    - Apply operations across patients
    - Provide cohort-level utilities
    """

    def __init__(self, patients):
        """
        Parameters
        ----------
        patients : list[Patient]
            List of fully constructed Patient objects.
        """
        self.patients = patients

    def verify_all_geometry(self):
        """
        Run spatial consistency checks for all patients.

        Raises
        ------
        ValueError
            If any patient fails geometry validation.
        """
        for patient in self.patients:
            patient.verify_geometry()

    def __len__(self):
        """Return number of patients in the dataset."""
        return len(self.patients)

    def __iter__(self):
        """Allow iteration: for patient in dataset"""
        return iter(self.patients)
    
    def extract_all_lesions(self):
        """
        Extract lesions for all patients in the dataset.
        """
        for patient in self.patients:
            patient.extract_lesions()

