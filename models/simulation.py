import numpy as np
from scipy import stats # statistics module (CI, binomial stats, t-distributions)
from scipy.stats import norm # normal distribution object within stats module (Gaussian needle placement errors)


class BiopsySimulation:
    """
    Handles patient-specific biopsy simulations
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

        self.valid_holes = [] # empty list to store valid template holes that can reach the prostate depending on patient


    def restrict_to_prostate(   # filters template holes to find those which intersect the prostate
        self,
        step_mm=1.0,    # distance between sampled points along the needle (sample every 1mm)
        max_depth_mm=150.0,     # max needle insertion depth
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
        affine = self.patient.affines["t2"]     # load MRI affine transformation matrix
        inv_affine = np.linalg.inv(affine)      # computes inverse matrix (world -> voxel)

        prostate = self.patient.prostate_mask   # load prostate mask (3d binary mask)
        shape = prostate.shape                  # load prostate mask shape (valid voxel indices)

        valid = []  # Temporary list of holes that intersect the prostate

        for hole in self.template.holes:    # iterate through every hole in template
            if self._hole_hits_prostate(    # calls helper function (defined below) to check validity
                hole,
                prostate,
                shape,
                inv_affine,
                step_mm,
                max_depth_mm,
            ):
                valid.append(hole)          # if valid, append to temporary valid list

        self.valid_holes = valid            # Save temporary valid list elements into main list of valid holes
        return valid

    def _hole_hits_prostate(                # HELPER FUNCTION: Check whether a single hole intersects the prostate
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
        origin = hole["origin"]     # Extract hole origin
        direction = hole["direction"]   # Extract hole direction

        n_steps = int(max_depth_mm / step_mm)   # compute sampling steps (150/1=150 samples)

        for i in range(n_steps):    # iterates along needle path
            point_world = origin + i * step_mm * direction  # traces needle line 

            point_h = np.append(point_world, 1.0)   # convert to homogenous coordinates, add 1 to end of coord for affine transformation
            voxel = inv_affine @ point_h    # matrix multiplication to convert world -> MRI voxel coords
            voxel = voxel[:3] # keep only x,y,z

            idx = np.round(voxel).astype(int) # round voxel coord to integer values

            if np.any(idx < 0) or np.any(idx >= shape):     # ensure voxel index is inside the mask
                continue

            if prostate_mask[tuple(idx)] > 0:   # if continue (after verififying voxel index is inside the mask, test if voxel coord correspond to being within the prostate, if so True, if not False
                return True

        return False
    

    def select_target_holes(self):  # HELPER FUNCTION: Choose the best hole for each lesion
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

        if not self.valid_holes:    # Ensure filtering has been done first
            raise RuntimeError("Must run restrict_to_prostate() first.")    # Stops execution with error message

        selected = {}   # Prepare output dictionary which maps lesions ID: selected best hole

        for lesion in self.patient.lesions:     # iterate through each lesion for a specific patient

            centroid = lesion.centroid      # obtain lesion centroid
            min_dist = np.inf   # start with infinite distance
            best_hole = None    # temporarily preset best hole to be no holes

            for hole in self.valid_holes:   # iterate through each valid hole
                hole_origin = hole["origin"]    # extract hole origin

                dist = np.linalg.norm(hole_origin - centroid)   # calculate distance between centroid of lesion and hole origin

                if dist < min_dist: # this if block ensures hole with origin closest to lesion centroid becomes best hole
                    min_dist = dist
                    best_hole = hole

            lesion.simulation_results["selected_hole"] = best_hole  # save result in Lesion Class Constructor
            lesion.simulation_results["hole_distance_mm"] = float(min_dist) # save distance result in Lesion Class Constructor

            selected[lesion.id] = best_hole # append to output dictionary

        return selected # output dictionary is output
    

    def define_ideal_needles(self, core_length_mm=20.0):    # HELPER FUNCITON: Defines ideal needles
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

        results = {}    # Prepare output dictionary which maps lesions ID: ideal needle 

        for lesion in self.patient.lesions: # iterate through each lesion for a specific patient

            if "selected_hole" not in lesion.simulation_results:    # ensure a best hole has already been selected
                raise RuntimeError(
                    "Must run select_target_holes() before defining needles."
                )

            centroid = lesion.centroid  # obtain lesion centroid
            direction = lesion.simulation_results["selected_hole"]["direction"] # get direction of best hole

            half_length = core_length_mm / 2.0

            # generate a 20mm core by going 10mm front and back from lesion centroid
            start = centroid - half_length * direction
            end = centroid + half_length * direction

            needle = {  # generate the needle object
                "start": start.astype(float),
                "end": end.astype(float),
                "length_mm": float(core_length_mm),
            }

            lesion.simulation_results["ideal_needle"] = needle  # add a new entry ideal needle: needle to Lesion Class Constructor
            results[lesion.id] = needle    # append to output dictionary 

        return results
    

    def discretise_needles(self, step_mm=1.0):  # Turns EACH needle into sample points (1mm along needle)
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

        results = {}    # Prepare output dictionary which maps lesions ID: points

        for lesion in self.patient.lesions: # iterate through each lesion for a specific patient

            if "ideal_needle" not in lesion.simulation_results: # ensure an ideal needle has already been selected
                raise RuntimeError(
                    "Must run define_ideal_needles() before discretisation."
                )

            needle = lesion.simulation_results["ideal_needle"]

            start = needle["start"]
            end = needle["end"]

            # Direction vector
            vec = end - start   # compute direction vector
            length = np.linalg.norm(vec)    # compute length of idrection vector

            if length == 0:
                raise ValueError("Needle has zero length.")

            direction = vec / length    # normalise direction 

            # Number of sampling steps
            n_steps = int(np.floor(length / step_mm)) + 1   # define number of sampling steps + 1

            points = []     # prepare sample points list

            for i in range(n_steps):
                point = start + i * step_mm * direction
                points.append(point)

            points = np.array(points, dtype=float)  # turns points list to an array

            lesion.simulation_results["needle_points"] = points     # add a new entry needle points: points to Lesion Class Constructor
            results[lesion.id] = points # append to output dictionary

        return results
    

    def generate_error_field(   # HELPER FUNCTION: Generates ERROR FIELD
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

        if "needle_points" not in lesion.simulation_results:    # ensure needle points has already been selected 
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
        sigmas = sigma_max_mm * (distances / length)    # error grows linearly along the needle

        # Sample Gaussian noise
        errors = norm.rvs(      # 
            loc=0.0,
            scale=sigmas[:, None],
            size=points.shape
        )

        return errors
    

    def apply_error_to_needles(     # HELPER FUNCTION: Applies error field to needles
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

            realised_points = ideal_points + errors     # Apply error to each needle points

            lesion.simulation_results["realised_needle_points"] = realised_points   # add a new entry realised needle points: points to Lesion Class Constructor

            results[lesion.id] = realised_points    # append to output dictionary

        return results


    def intersect_needles_with_lesions(self):       # HELPER FUNCITON: Test if needle point intersect with lesion mask
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

            if "realised_needle_points" not in lesion.simulation_results:   # ensure realised needle points has already been made
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

            lesion.simulation_results["lesion_intersections"] = intersections   # add a new entry lesion intersections: True or False to Lesion Class Constructor
            results[lesion.id] = intersections

        return results
    

    def compute_needle_outcomes(self, step_mm=1.0):     # HELPER FUNCTION: Compute needle outcome metrics such as hit_flag, positive_length_mm and percentage_positive
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

            if "lesion_intersections" not in lesion.simulation_results: # ensure lesion_intersections have alreayd been made
                raise RuntimeError(
                    "Must run intersect_needles_with_lesions() first."
                )

            intersections = lesion.simulation_results["lesion_intersections"]

            # Hit flag
            hit_flag = int(np.any(intersections))   # 1 if there is an intersection with lesion

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

            lesion.simulation_results["outcomes"] = outcomes    # add a new entry outocmes: needle outcome metrics for each lesion to Lesion Class Constructor
            results[lesion.id] = outcomes

        return results
    

    def run_monte_carlo(        # HELPER FUNCTION: Monte Carlo simulation
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

            for sim_i in range(n_simulations):  # iterate through each simulation for n simulations

                lesion_hits = 0
                lesion_percentages = []

                for core_i in range(n_cores):   # Foe each simuation, iterate through its 5 core samples

                    # Apply stochastic error
                    self.apply_error_to_needles(
                        sigma_max_mm=sigma_max_mm
                    )

                    self.intersect_needles_with_lesions()
                    self.compute_needle_outcomes(step_mm=step_mm)

                    outcome = lesion.simulation_results["outcomes"]

                    lesion_hits += outcome["hit_flag"]
                    lesion_percentages.append(outcome["percentage_positive"])

                # Record per SIMULATION (5 core summary)
                hit_flags.append(int(lesion_hits > 0))
                positive_percentages.append(
                    float(np.mean(lesion_percentages))
                )
                positive_core_counts.append(int(lesion_hits))

                if sim_i % 100 == 0 and sim_i > 0:  # Give summary stats every 100 simulations
                    current_estimate = np.mean(hit_flags)
                    print(
                        f"[Lesion {lesion.id}] "
                        f"Iteration {sim_i}: "
                        f"Hit probability ≈ {current_estimate:.3f}"
                    )

            # Store aggregated statistics (after simulating all simulations)

            # Final hit mean (hitflag mean)
            hit_mean = np.mean(hit_flags)

            # Generate 95% CI for hit probability / hitflag mean (binomial)
            ci_low, ci_high = stats.binom.interval(
                0.95,
                n=n_simulations,
                p=hit_mean
            )

            ci_low /= n_simulations
            ci_high /= n_simulations

            # Final postive percetage mean
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