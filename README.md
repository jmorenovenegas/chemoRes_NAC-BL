# chemoRes_NAC-BL 
WCGNA pre NAC samples. We look for correlation with pCR in basal-like samples.

Analysis steps:
- Data loading, preprocessing, normalization and cleaning.
- Network construction and module detection.
- Relate modules to clinical chemoresistance(-pCR).
- WGCNA relevant modules enrichment.
- Module expansion using DIAMOnD.py (user must set path to python source with the required packages).
- Subnetwork building and enrichment of expanded modules.
- Link communities analysis.

Complete analysis in chemoRes_NAC-BL.R script.
Functions used in RCode directory.
Interactomes used in DIAMOnD in Interactomes directory. 
RCode/interactomes_processing.R is the script used to process the interactomes. More information is available in the main script.

R version 3.5.1



