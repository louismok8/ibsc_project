# MRI-Targeted Prostate Biopsy Simulation Project

## Project Overview

This project aims to **computationally estimate the true hit rate of MRI-targeted transperineal prostate biopsy**, accounting for realistic needle placement uncertainty. The goal is to understand how often a biopsy needle that is *intended* to target an MRI-visible lesion actually intersects that lesion, given known sources of targeting imprecision.

The work is motivated by the clinical reality that MRI-targeted biopsies are often performed using **cognitive targeting**: clinicians visually map MRI findings onto 2D transrectal or transperineal ultrasound images and manually guide the biopsy needle. This process is inherently operator-dependent and subject to multiple sources of error.

Rather than modelling biopsy outcomes at the patient level, this project explicitly performs a **lesion-level analysis**, treating each MRI-defined lesion as an independent target.

---

## Clinical Motivation

In current clinical practice:

- MRI provides high-resolution 3D images with:
  - Prostate segmentation
  - One or more MRI-visible suspicious lesions
- Biopsies are typically performed under ultrasound guidance
- Clinicians cognitively map MRI findings onto ultrasound images
- Targeting is affected by:
  - Operator judgment
  - Tissue deformation
  - Needle bending (especially at the tip)
  - Template constraints
  - Depth control errors

As a result:
- A biopsy core may be positive without having sampled the *intended* lesion
- Observed positive core rates may overestimate true targeting accuracy

This project attempts to **separate intended targeting from actual sampling success**.

---

## Core Research Question

> Given a biopsy needle that is *perfectly aimed* at the centre of an MRI-defined lesion, what is the probability that, after realistic targeting errors, the needle actually intersects that lesion?

This probability is referred to as the **true lesion hit rate**.

---

## Key Design Principles

### 1. Lesion-Level Analysis (Not Patient-Level)

- **Unit of analysis**: individual MRI-defined lesion
- A patient with multiple lesions contributes multiple independent simulations
- Lesions are *not* merged or collapsed into a single patient outcome

This enables:
- Lesion-size–dependent analysis
- Accurate estimation of targeting performance per lesion

---

### 2. MRI as Ground Truth for Target Definition

- MRI lesion masks define the *intended targets*
- Pathology is **not** used to define targeting success
- A hit is defined purely geometrically

---

### 3. Hit Definition

A biopsy core is considered a **hit** if:

- **Any part** of the needle core intersects the MRI-defined lesion mask

Notes:
- Partial intersections count as hits
- Full and partial overlaps are treated equally for hit-rate estimation
- Percentage of positive core length is also recorded as a secondary metric

---

### 4. Error Modelling Philosophy

The project **does not explicitly model biomechanics**.

Instead:
- A stochastic needle placement error model is used
- This implicitly captures the combined effect of:
  - Cognitive targeting error
  - Tissue deformation
  - Needle bending
  - Depth uncertainty

The error is applied to the *final realised needle trajectory*.

---

## Data Description

For each patient, the following files are provided (in `.nii.gz` format):

- **T2-weighted MRI**
- **Lesion mask**
  - Integer labels: `0 … L`
  - `0` = background
  - Each `ℓ = 1 … L` corresponds to one lesion
- **Prostate mask**
  - Binary: `1 = prostate`, `0 = background`

Patients may have **multiple lesions**.

---

## Simulation Framework

### Biopsy Geometry

- Transperineal template biopsy
- 19 × 19 grid
- 5 mm spacing
- Fixed insertion direction
- Needle core length ≈ 15–20 mm

Only grid holes that intersect the prostate are considered valid.

---

### Simulation Workflow (Per Lesion)

1. Extract binary mask for the lesion
2. Compute lesion volume and centroid
3. Select the closest valid template grid hole to the lesion centroid
4. Define an ideal (error-free) needle trajectory:
   - Core centred on lesion centroid
5. Discretise the needle core into sampled points
6. Apply stochastic placement error to generate realised trajectories
7. Test needle–lesion intersection
8. Record outcomes

---

### Monte Carlo Simulation

For each lesion:

- A fixed number of targeted needles is simulated (e.g. 5)
- Simulation is repeated **1,000–100,000 times**
- Statistics are collected until convergence

Recorded metrics include:
- Hit probability
- Number of positive needles per simulation
- Percentage positive core length per needle

---

## Lesion Size Stratification

Lesions are grouped by volume into predefined size bins:

- **Small lesions**
- **Medium lesions**
- **Large lesions**

Important:
- The simulation is identical for all lesions
- Size groups are used **only for analysis and aggregation**

Statistics are computed **separately for each size group**.

---

## Outputs

Primary outputs include:

- Per-lesion true hit rate
- Distribution of number of positive cores per lesion
- Distribution of positive core percentages
- Size-stratified hit-rate statistics

These outputs can be compared against:
- Observed clinical biopsy positivity rates
- Without assuming that observed positivity equals correct targeting

---

## Interpretation Goals

The project aims to:

- Quantify how lesion size affects sampling success
- Explain why clinical positive core rates may overestimate targeting accuracy
- Provide a principled estimate of true MRI-targeted biopsy performance

---

## Scope Limitations

This project explicitly does **not**:

- Model ultrasound imaging
- Perform MRI–ultrasound registration
- Model tissue biomechanics explicitly
- Infer cancer grade or pathology extent

These are outside scope and intentionally abstracted into the error model.

---

## Summary

This project provides a **lesion-level, MRI-grounded, simulation-based estimate** of prostate biopsy targeting accuracy under realistic uncertainty. It bridges the gap between ideal MRI targeting and real-world biopsy outcomes, offering insight into the reliability—and limitations—of current clinical practice.
