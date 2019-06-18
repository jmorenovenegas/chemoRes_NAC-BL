# chemoRes_NAC-BL
WCGNA pre NAC samples. We look for correlation with pCR in basal-like samples.  
Analysis steps:  
(chemoRes_NAC-BL.R)  
- Data loading, preprocessing, normalization and cleaning.
- Network construction and module detection.
- Relate modules to clinical chemoresistance(-pCR).
- WGCNA relevant modules enrichment.
- Module expansion using DIAMOnD.py (user must set path to python source with the required packages).
- Subnetwork building and enrichment of expanded modules.
- Link communities analysis subnetworks preprocessing.  

(link_communities.R)  
- Link communities analysis  

Intructions to perform the analysis:
1. Install python(2.7 or higher) and required libraries
2. Set path to python source in chemoRes_NAC-BL.R script(pybin variable)
3. Run chemoRes_NAC-BL.R (Recommended in RStudio)
4. Run link_communities.R (recommended in RStudio to adjust parameters in the script). You must pass .gml files obtained in the previous step as input.

Used functions available in RCode directory.

Interactomes used in DIAMOnD available in Interactomes directory.

RCode/interactomes_processing.R is the script used to process the interactomes. 

More information is available in the main script.

R version 3.5.1 in RStudio
- NanoStringQCPro
- dplyr
- WGCNA
- biomaRt
- clusterProfiler
- linkcomm
- ggplot2
- igraph
- org.Hs.eg.db
- biomaRt
- pathview
- ReactomePA

Python version 2.7.15
- networkx
- numpy
- scipy



