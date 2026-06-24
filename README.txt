START HERE — REPRODUCIBLE HEATMAP PACKAGE
==========================================

1. Double-click run_heatmap_windows.bat on Windows.
   It installs Python packages and creates the PNG and PDF heatmaps.

2. To change the plotted data, edit centered_matrix.csv in Excel.
   Only change values in sample columns N1 through HF10.

3. To change color intensity, open make_heatmap.py and edit:
       COLOR_LIMIT = 3.5
   Smaller numbers show stronger red/blue colors. Larger numbers show softer colors.

4. The script checks for all 35 expected people and exactly one TMEM94 row.

What the values mean:
Each value is already a person's normalized expression minus the mean normalized
expression of the nine Normal samples for the same gene.

Blue = lower than Normal-group mean; white = near Normal-group mean; red = higher.

The heatmap is descriptive. Do not treat individual colored cells as statistical
significance results.
