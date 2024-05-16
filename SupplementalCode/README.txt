Files used to reproduce empirical analyses in Kortessis et al. "Increasing environmental fluctuations can dampen variability of endogenously cycling populations"



1. Folder list
  A.  BeanBeetle
    Short description: Contains .R files and data to analyze Desharnais (1979) bean beetle experimental data.
    
  B. LynxHare
    Short description: Contains .R files and data to analyze Lynx dataset.

2. File list
   A. Filename: /BeanBeetle/BB_MLE.R    
      Short description: Run this file to reproduce analysis of bean beetle dataset. Estiamtes model parameters and estimates variance components.

   B. Filename: /BeanBeetle/ModelFuncs.R
      Short description: Support functions for BB_MLE.R

   C. Filename: /BeanBeetle/BeanBeatle_LPA.csv
      Short description: Dataset from Desharnais (1979).

   D. Filename: /LynxHare/LV_MLE.R    
      Short description: Run this file to reproduce analysis of Lynx dataset. Estimates model parameters and estimates variance components.

   E. Filename: /LynxHare/ModelFuncs.R
      Short description: Support functions for LV_MLE.R

   F. Filename: /LynxHare/LynxHare.txt
      Short description: Dataset from http://people.whitman.edu/~hundledr/courses/M250F03/.
        

2. Relationship between files:        

To run the bean beetle analysis run /BeanBeetle/BB_MLE.R
To run the Lynx analysis run /LynxHare/LV_MLE.R
