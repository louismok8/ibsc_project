import numpy as np
from scipy import stats
from scipy.stats import norm


class BiopsySimulation:
    """
    Handles patient-specific biopsy feasibility checks.
    """

    def __init__(self, patient, template):
        """
        Parameters
        ----------------
        patient : Patient
            Patient object with prostate mask and affine
        template : BiopsyTemplate
            Idealised biopsy template
        """
        self.patient = patient
        self.template = template

        self.valid_holes = []

    def restrict_to_prostate(
        self,
        step_mm=1.0,
        max_depth_mm=150.0,
    ):
        """
        Keep only template holes whose needle trajectory
        intersects the prostate mask.

        Parameters
        ----------
        step_mm : float
            Sampling step along needle (mm)

        max_depth_mm : float
            Maximum needle insertion depth (mm)
        """
        affine = self.patient.affines["t2"]
        inv_affine = np.linalg.inv(affine)

        prostate = self.patient.prostate_mask
        shape = prostate.shape

        valid = []

        for hole in self.template.holes:
            if self._hole_hits_prostate(
                hole,
                prostate,
                shape,
                inv_affine,
                step_mm,
                max_depth_mm,
            ):
                valid.append(hole)

        self.valid_holes = valid
        return valid

    def _hole_hits_prostate(
        self,
        hole,
        prostate_mask,
        shape,
        inv_affine,
        step_mm,
        max_depth_mm,
    ):
        """
        Check whether a single hole intersects the prostate.
        """
        origin = hole["origin"]
        direction = hole["direction"]

        n_steps = int(max_depth_mm / step_mm)

        for i in range(n_steps):
            point_world = origin + i * step_mm * direction

            point_h = np.append(point_world, 1.0)
            voxel = inv_affine @ point_h
            voxel = voxel[:3]

            idx = np.round(voxel).astype(int)

            if np.any(idx < 0) or np.any(idx >= shape):
                continue

            if prostate_mask[tuple(idx)] > 0:
                return True

        return False
    

    def select_target_holes(self):
        """
        For each lesion, select the closest valid template hole
        based on Euclidean distance in world coordinates.

        Stores result inside each lesion:
            lesion.simulation_results["selected_hole"]

        Returns
        -------
        dict
            Mapping (lesion_id -> selected hole dict)
        """

        if not self.valid_holes:
            raise RuntimeError("Must run restrict_to_prostate() first.")

        selected = {}

        for lesion in self.patient.lesions:

            centroid = lesion.centroid
            min_dist = np.inf
            best_hole = None

            for hole in self.valid_holes:
                hole_origin = hole["origin"]

                dist = np.linalg.norm(hole_origin - centroid)

                if dist < min_dist:
                    min_dist = dist
                    best_hole = hole

            lesion.simulation_results["selected_hole"] = best_hole
            lesion.simulation_results["hole_distance_mm"] = float(min_dist)

            selected[lesion.id] = best_hole

        return selected
    

    def define_ideal_needles(self, core_length_mm=20.0):
        """
        Construct ideal (error-free) biopsy core segments
        centred at each lesion centroid.

        Parameters
        ----------
        core_length_mm : float
            Length of biopsy core (mm)

        Stores in each lesion:
            lesion.simulation_results["ideal_needle"]
                {
                    "start": np.ndarray (3,),
                    "end": np.ndarray (3,),
                    "length_mm": float
                }

        Returns
        -------
        dict
            Mapping lesion_id -> needle dict
        """

        results = {}

        for lesion in self.patient.lesions:

            if "selected_hole" not in lesion.simulation_results:
                raise RuntimeError(
                    "Must run select_target_holes() before defining needles."
                )

            centroid = lesion.centroid
            direction = lesion.simulation_results["selected_hole"]["direction"]

            half_length = core_length_mm / 2.0

            start = centroid - half_length * direction
            end = centroid + half_length * direction

            needle = {
                "start": start.astype(float),
                "end": end.astype(float),
                "length_mm": float(core_length_mm),
            }

            lesion.simulation_results["ideal_needle"] = needle
            results[lesion.id] = needle

        return results
    

    def discretise_needles(self, step_mm=1.0):
        """
        Discretise each ideal needle into regularly spaced 3D points.

        Parameters
        ----------
        step_mm : float
            Sampling distance along the needle (mm)

        Stores in each lesion:
            lesion.simulation_results["needle_points"]
                np.ndarray shape (N, 3)

        Returns
        -------
        dict
            Mapping lesion_id -> array of sampled points
        """

        results = {}

        for lesion in self.patient.lesions:

            if "ideal_needle" not in lesion.simulation_results:
                raise RuntimeError(
                    "Must run define_ideal_needles() before discretisation."
                )

            needle = lesion.simulation_results["ideal_needle"]

            start = needle["start"]
            end = needle["end"]

            # Direction vector
            vec = end - start
            length = np.linalg.norm(vec)

            if length == 0:
                raise ValueError("Needle has zero length.")

            direction = vec / length

            # Number of sampling steps
            n_steps = int(np.floor(length / step_mm)) + 1

            points = []

            for i in range(n_steps):
                point = start + i * step_mm * direction
                points.append(point)

            points = np.array(points, dtype=float)

            lesion.simulation_results["needle_points"] = points
            results[lesion.id] = points

        return results
    

    def generate_error_field(
        self,
        lesion,
        sigma_max_mm=3.0,
    ):
        """
        Generate Gaussian spatial perturbations along a discretised needle.

        Error increases linearly from template end to tip.

        Parameters
        ----------
        lesion : Lesion
            Lesion object with discretised needle

        sigma_max_mm : float
            Maximum standard deviation (mm) at needle tip

        Returns
        -------
        np.ndarray
            Array of shape (N, 3) representing XYZ perturbations
            to apply to each sampled needle point
        """

        if "needle_points" not in lesion.simulation_results:
            raise RuntimeError(
                "Must run discretise_needles() before generating error."
            )

        points = lesion.simulation_results["needle_points"]

        start = points[0]
        end = points[-1]

        vec = end - start
        length = np.linalg.norm(vec)

        if length == 0:
            raise ValueError("Needle length is zero.")

        direction = vec / length

        # Compute distance of each point from template end
        distances = np.linalg.norm(points - start, axis=1)

        # Linear variance growth
        sigmas = sigma_max_mm * (distances / length)

        # Sample Gaussian noise
        errors = norm.rvs(
            loc=0.0,
            scale=sigmas[:, None],
            size=points.shape
        )

        return errors
    

    def apply_error_to_needles(
        self,
        sigma_max_mm=3.0,
        random_seed=None,
    ):
        """
        Apply stochastic placement error to all discretised needles.

        Parameters
        ----------
        sigma_max_mm : float
            Maximum Gaussian standard deviation at needle tip (mm)

        random_seed : int or None
            Optional seed for reproducibility

        Stores in each lesion:
            lesion.simulation_results["realised_needle_points"]

        Returns
        -------
        dict
            Mapping lesion_id -> realised needle points array
        """

        if random_seed is not None:
            np.random.seed(random_seed)

        results = {}

        for lesion in self.patient.lesions:

            if "needle_points" not in lesion.simulation_results:
                raise RuntimeError(
                    "Must run discretise_needles() before applying error."
                )

            ideal_points = lesion.simulation_results["needle_points"]

            # Generate spatial perturbation field
            errors = self.generate_error_field(
                lesion,
                sigma_max_mm=sigma_max_mm
            )

            realised_points = ideal_points + errors

            lesion.simulation_results["realised_needle_points"] = realised_points

            results[lesion.id] = realised_points

        return results


    def intersect_needles_with_lesions(self):
        """
        Determine which realised needle points intersect each lesion mask.

        For every lesion:
            - Convert realised needle points from world → voxel coordinates
            - Sample lesion mask
            - Store boolean array of intersection flags

        Stores in each lesion:
            lesion.simulation_results["lesion_intersections"]
                np.ndarray of shape (N,) of dtype bool

        Returns
        -------
        dict
            Mapping lesion_id -> boolean intersection array
        """

        affine = self.patient.affines["t2"]
        inv_affine = np.linalg.inv(affine)

        results = {}

        for lesion in self.patient.lesions:

            if "realised_needle_points" not in lesion.simulation_results:
                raise RuntimeError(
                    "Must run apply_error_to_needles() before hit testing."
                )

            points_world = lesion.simulation_results["realised_needle_points"]
            lesion_mask = lesion.mask
            shape = lesion_mask.shape

            intersections = []

            for point_world in points_world:

                # Convert world → voxel
                point_h = np.append(point_world, 1.0)
                voxel = inv_affine @ point_h
                voxel = voxel[:3]

                idx = np.round(voxel).astype(int)

                # Check bounds
                if np.any(idx < 0) or np.any(idx >= shape):
                    intersections.append(False)
                    continue

                # Check lesion mask
                if lesion_mask[tuple(idx)] > 0:
                    intersections.append(True)
                else:
                    intersections.append(False)

            intersections = np.array(intersections, dtype=bool)

            lesion.simulation_results["lesion_intersections"] = intersections
            results[lesion.id] = intersections

        return results
    

    def compute_needle_outcomes(self, step_mm=1.0):
        """
        Compute biopsy outcome metrics for each lesion.

        Metrics:
            - hit_flag (0 or 1)
            - positive_length_mm
            - percentage_positive (0–1)

        Parameters
        ----------
        step_mm : float
            Sampling spacing used during discretisation.

        Stores in each lesion:
            lesion.simulation_results["outcomes"]

        Returns
        -------
        dict
            Mapping lesion_id -> outcome dict
        """

        results = {}

        for lesion in self.patient.lesions:

            if "lesion_intersections" not in lesion.simulation_results:
                raise RuntimeError(
                    "Must run intersect_needles_with_lesions() first."
                )

            intersections = lesion.simulation_results["lesion_intersections"]

            # Hit flag
            hit_flag = int(np.any(intersections))

            # Positive length
            n_positive = int(np.sum(intersections))
            positive_length = float(n_positive * step_mm)

            # Total core length
            needle = lesion.simulation_results["ideal_needle"]
            total_length = needle["length_mm"]

            percentage_positive = (
                positive_length / total_length
                if total_length > 0 else 0.0
            )

            outcomes = {
                "hit_flag": hit_flag,
                "positive_length_mm": positive_length,
                "percentage_positive": float(percentage_positive),
            }

            lesion.simulation_results["outcomes"] = outcomes
            results[lesion.id] = outcomes

        return results
    

    def run_monte_carlo(
        self,
        n_simulations=1000,
        n_cores=5,
        sigma_max_mm=3.0,
        step_mm=1.0,
    ):
        """
        Run Monte Carlo biopsy simulation per lesion.

        Parameters
        ----------
        n_simulations : int
            Number of stochastic repeats

        n_cores : int
            Number of targeted cores per lesion

        sigma_max_mm : float
            Maximum placement error at needle tip

        step_mm : float
            Sampling spacing

        Returns
        -------
        dict
            Per-lesion Monte Carlo statistics
        """

        results = {}

        for lesion in self.patient.lesions:

            hit_flags = []
            positive_percentages = []
            positive_core_counts = []

            for sim_i in range(n_simulations):

                lesion_hits = 0
                lesion_percentages = []

                for core_i in range(n_cores):

                    # Apply stochastic error
                    self.apply_error_to_needles(
                        sigma_max_mm=sigma_max_mm
                    )

                    self.intersect_needles_with_lesions()
                    self.compute_needle_outcomes(step_mm=step_mm)

                    outcome = lesion.simulation_results["outcomes"]

                    lesion_hits += outcome["hit_flag"]
                    lesion_percentages.append(outcome["percentage_positive"])

                # Record per simulation
                hit_flags.append(int(lesion_hits > 0))
                positive_percentages.append(
                    float(np.mean(lesion_percentages))
                )
                positive_core_counts.append(int(lesion_hits))

                if sim_i % 100 == 0 and sim_i > 0:
                    current_estimate = np.mean(hit_flags)
                    print(
                        f"[Lesion {lesion.id}] "
                        f"Iteration {sim_i}: "
                        f"Hit probability ≈ {current_estimate:.3f}"
                    )

            # Store aggregated statistics
            hit_mean = np.mean(hit_flags)

            # 95% CI for hit probability (binomial)
            ci_low, ci_high = stats.binom.interval(
                0.95,
                n=n_simulations,
                p=hit_mean
            )

            ci_low /= n_simulations
            ci_high /= n_simulations

            mean_percent = np.mean(positive_percentages)

            # 95% CI for mean percentage positive (t-based)
            ci_percent = stats.t.interval(
                0.95,
                df=len(positive_percentages)-1,
                loc=mean_percent,
                scale=stats.sem(positive_percentages)
            )

            results[lesion.id] = {
                "hit_probability": float(hit_mean),
                "hit_probability_ci": (float(ci_low), float(ci_high)),
                "mean_percentage_positive": float(mean_percent),
                "mean_percentage_positive_ci": tuple(ci_percent),
                "distribution_percentage_positive": np.array(positive_percentages),
                "distribution_positive_core_counts": np.array(positive_core_counts),
            }

        return results