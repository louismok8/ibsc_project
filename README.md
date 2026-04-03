# MRI-Targeted Prostate Biopsy Simulation Project

## Project Overview

This project aims to **computationally estimate the true hit rate of MRI-targeted transperineal prostate biopsy**, accounting for realistic needle placement uncertainty. The goal is to understand how often a biopsy needle that is *intended* to target an MRI-visible lesion actually intersects that lesion, given known sources of targeting imprecision.

The work is motivated by the clinical reality that MRI-targeted biopsies are often performed using **cognitive targeting**: clinicians visually map MRI findings onto 2D transrectal or transperineal ultrasound images and manually guide the biopsy needle. This process is inherently operator-dependent and subject to multiple sources of error.

Rather than modelling biopsy outcomes at the patient level, this project explicitly performs a **lesion-level analysis**, treating each MRI-defined lesion as an independent target.

---

## Project structure

```bash
iBSC/Project/
│
├── data/                  
│   ├── paths.py           # Paths
├── imagesTr/              # MRI scans
├── labelsTr/              # Lesion masks
├── zonesTr/               # Prostate masks
│
├── models/                # Core classes
│   ├── dataset.py
│   ├── patients.py
│   ├── lesion.py
│   ├── template.py
│   └── simulation.py
│
├── pipelines/
│   └── run_simulation_pipeline.py
│
├── artificial data/              # Mock data generation
├── outputs/               # Results
```

---

## Core Research Question

> Given a biopsy needle that is *perfectly aimed* at the centre of an MRI-defined lesion, what is the probability that, after realistic targeting errors, the needle actually intersects that lesion?

This probability is referred to as the **true lesion hit rate**.

---

## Outputs
Lesion-level → simulation_results.csv

Patient-level → patient_predictions.csv

Zone assignment → lesion_zones.csv

Zone predictions → zone_level_predictions.csv

Zone metrics → zone_accuracy_metrics.csv

---

## How to Run
### Generate mock data
python analysis/generate_mock_real_patient_results.py
python analysis/generate_mock_real_zone_results.py

### Run full pipeline
python pipelines/run_simulation_pipeline.py

---

## Key Parameters
n_simulations → number of Monte Carlo runs (default: 1000)

n_cores → biopsy cores per lesion (default: 5)

sigma_max_mm → max needle error (3mm for MRI-targeted, 5mm for standard TRUS)

step_mm → sampling resolution (1mm along the needle)

---

## Summary

This project provides a **lesion-level, MRI-grounded, simulation-based estimate** of prostate biopsy targeting accuracy under realistic uncertainty. It bridges the gap between ideal MRI targeting and real-world biopsy outcomes, offering insight into the reliability and limitations of current clinical practice.
