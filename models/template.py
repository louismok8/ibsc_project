import numpy as np


class BiopsyTemplate:
    """
    Represents an idealised transperineal biopsy template grid.

    The template is defined in MRI world coordinates (mm)
    and is independent of any patient anatomy.
    """

    def __init__(
        self,
        grid_size=19,
        spacing=5.0,
        origin=(0.0, 0.0, 0.0),
        direction=(0.0, 0.0, 1.0),
    ):
        """
        Parameters
        ----------
        grid_size : int
            Number of holes per axis (default 19 → 19×19 grid)

        spacing : float
            Distance between adjacent holes (mm)

        origin : tuple[float, float, float]
            World-space center of the grid (mm)

        direction : tuple[float, float, float]
            Needle insertion direction (world space)
        """
        self.grid_size = grid_size
        self.spacing = float(spacing)
        self.origin = np.asarray(origin, dtype=float) # Convert origin to Numpy Array
        d = np.asarray(direction, dtype=float) # Convert direction to Numpy Array
        self.direction = d / np.linalg.norm(d) # Normalise direction vector to be a unit vector

        self.holes = self._build_grid()


    # Helper function
    def _build_grid(self):
        """
        Build grid hole origins in a plane.
        """
        half = self.grid_size // 2  # 19 // 2 = 9 (meaning grid will run from -9 to 9)
        coords = [] # List storage for every hole in the template

        for i in range(-half, half + 1): # horizontal movement
            for j in range(-half, half + 1): # vertical movement
                offset = np.array([ # calculates hole position relative to center
                    i * self.spacing,
                    j * self.spacing,
                    0.0
                ])

                point = self.origin + offset # point world coordinate of each hole

                coords.append({
                    "origin": point,
                    "direction": self.direction
                })

        return coords

    def __len__(self):
        return len(self.holes)