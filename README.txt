DNA Methylation Analysis Shiny App
==================================

Overview
This Shiny app allows users to upload DNA methylation data and perform a complete downstream analysis, including visualization, differential methylation, pathway enrichment, power analysis, and classification via Random Forest.
File Requirements
- Input must be a .csv file.
- The file should include:
•	A 'Gene' column (probe/gene identifiers)
•	A 'Region' column (e.g., Island, Shore, Shelf)
•	Multiple sample columns (e.g., Tumor1, Tumor2, Normal1, etc.)
Usage Instructions
1. Upload your methylation data using the file input.
2. Select the sample groups (detected automatically).
3. Use the buttons in the sidebar to:
    - Generate boxplots and heatmaps
    - Run differential methylation (DE) analysis
    - Perform pathway enrichment (via enrichR)
    - Conduct power analysis
    - Visualize sample clustering (PCA / UMAP)
    - Train and evaluate a Random Forest classifier
Random Forest Tab
- Trains a classifier using probes with adjusted p-value ≤ 0.05 from DE results.
- Includes three subtabs:
•	Prediction: actual vs predicted group labels
•	Metrics: accuracy, sensitivity, specificity, AUC
•	Feature Importance: top features ranked by importance
- Each subtab includes a download button.

Tips
- Choose exactly two groups for DE, power, and Random Forest.
- Use consistent sample name prefixes (e.g., Tumor1, Tumor2).
- All plots are interactive (via Plotly) where applicable.
- Download buttons are provided for results in each section.

Contact
jcaperella@gmail.com

Note
If you encounter issues, double-check your input formatting, or check for NA values in key columns.
