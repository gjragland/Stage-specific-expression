# Stage-specific-expression

## Analysis of RNAseq data and RNAi phenotype for Freda et al., xx

Transcript mapping count data and code analyzing stage-specific patterns of gene expression during cold exposure in Drosophila melanogaster

AND

Analysis of RNAi line phenotypes

### Data (mapping counts output from RSEM:

Count output - 'geneCounts.tar.gz'

### Code:

#### RNAseq analysis

* Main block including glm (EdgeR) modeling - 'EdgeR_modelsAndTrajectories.R'

  * Calls
  
    * 'EdgeR_LineTrajectoryContrasts.r' (generate line trajectories for line phenotype plotting)
    
    * 'tissueSpecificityInFredaEtAlData.R' (analyze relationship between DE in this study and tissue-specific expression from flyatlas2 data)
    
      * Calls - 'calculateTauForTissueSpecificity.R' (calculate tau values from flyatlas2 data)
 
 #### RNAi analysis 
 
 csv data file available in supplemental materials to publication
 
 * 'RNAiStatsv6.R' (apply glm models to RNAi phenotypes)
 * 'RNAiPlotsv6.R' (plot RNAi phenotypes)
