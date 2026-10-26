# Compact upload edition

This edition fits below the requested 25 MB upload limit. It preserves all simulation code, original experimental inputs, figure-specific parameters, manuscript images, native formulas, calculation summaries and numerical arrays required by the default plotting and verification commands.

Removed files are generated plots/coordinates and duplicate component previews, plus eight regenerable coupling/optical/phase-sensitivity trajectory archives. Their control-summary records and complete generating code remain included. No experimental measurements or simulation parameters were removed.

    python -m pip install -r requirements-lock.txt
    python reproduce.py --mode verify

This regenerates every manuscript plot, plotted coordinates, tables and verification outputs.

    python reproduce.py --mode full

This additionally regenerates the omitted trajectory archives by recalculating the full numerical suite. The root README contains the full scientific scope and provenance qualifications. The compact ZIP does not contain Git history; it is ready to extract and upload as repository contents.

RELEASE_SHA256SUMS.json records the actual compact files. Earlier input-package manifests in docs/input_package_manifests/ remain historical provenance records.
