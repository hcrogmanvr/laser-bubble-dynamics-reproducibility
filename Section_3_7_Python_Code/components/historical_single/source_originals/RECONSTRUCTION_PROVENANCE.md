# Laser-Driven Multi-Bubble Collapse and Plasma-Mediated Optical Flashes

Reproducibility code for the study **"Laser-Driven Multi-Bubble Collapse and Plasma-Mediated Optical Flashes: A Coupled Dynamical Model"**.

## Scope
This repository implements the reduced-order model described in the manuscript: laser gating, phenomenological optical darkening, plasma-induced absorption quenching, internal-energy evolution, an EOS-capped pressure closure, reduced Rayleigh-Plesset bubble dynamics, collapse-triggered bursts, camera response, mean-field multi-bubble coupling, and an impulse-based membrane-damage proxy.

> Important: the original source scripts were not attached to this chat. This package is a clean reproducibility implementation reconstructed from the equations and model description in the manuscript. Phenomenological closures are clearly labeled and should be checked against the authors' archived original scripts before publication.

## Repository layout
- `src/model.py` — core equations and single/multi-bubble solvers
- `src/io_utils.py` — YAML parameter loading
- `configs/default.yaml` — baseline model parameters
- `scripts/run_single.py` — single-bubble simulation and figures
- `scripts/run_multi.py` — multi-bubble ensemble summary
- `scripts/calibrate_growth.py` — R^2-space growth calibration
- `outputs/` — generated CSV and PNG outputs

## Installation
```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Run
```bash
python scripts/run_single.py
python scripts/run_multi.py
python scripts/calibrate_growth.py
```

## Model provenance
Established physics used in the manuscript includes Rayleigh-Plesset bubble dynamics, laser-induced breakdown/cavitation literature, and plasma/sonoluminescence context. Several closures used here (transmission saturation, ionization quenching, collapse threshold/reset, camera mapping, membrane damage proxy) are phenomenological components of the study and are not claimed to be first-principles laws.

## Citation
See `CITATION.cff`.

## License
MIT for code. Manuscript text, figures, and experimental source data remain subject to their original publication rights.
