FIGURE 8 - ORIGINAL EXPERIMENTAL GROWTH ANALYSIS (NO BAR GRAPHS)

Run:
  pip install numpy scipy matplotlib openpyxl
  python figure8.py --input "The four Graphs Fig-2 Data (2).xlsx" --output .

Figures:
  Figure_8A_Red_Growth.png/.pdf              representative red growth-phase D² fit and diameter
  Figure_8B_Four_Traces_and_Slope_CI.png/.pdf four growth-phase D² fits and dot/whisker 95% intervals (NO BAR GRAPH)
  Figure_8C_Full_Expansion_Contraction.png/.pdf original complete measured diameter time-series

Tables:
  Figure_8_Results_Table.csv and .xlsx
  In workbook: regression summary, notes, and original source measurements.

Method:
  Time in source is minutes; convert to seconds by multiplying by 60.
  Diameter is interpreted as micrometers because original paper Figure 2 labels µm,
  despite the source spreadsheet labeling the diameter column as mm.
  For each trace, identify first observed maximum diameter and fit from the first
  measurement through that maximum (included). Ordinary least-squares linear
  regression is performed on diameter squared D² versus time. 95% slope CI uses
  Student t with n-2 degrees of freedom and regression standard error.
  Figure display time is 0-60 s, but fit includes maxima at 60.54s (red) and
  61.26s (green); no fit data are silently excluded.

CAUTIONS:
  The Y-error=2 column's semantics are undocumented; these bars are not shown.
  Dot/whisker intervals refer to regression slope uncertainty, not raw Y-error.
  The 30-80 um²/s published range is descriptive, not a probabilistic CI.
  Curves are empirical experimental fits, NOT independent simulation predictions.
  The green and red maxima beyond 60 s are included in regressions, but not
  shown in the 0-60 s figure panels.

Reference: Erenso et al. (2024), European Biophysics Journal 53, 123-131.
DOI https://doi.org/10.1007/s00249-024-01701-3
