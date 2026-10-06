# Batch Correction Assessment Discussion - Text Only

## Original Problem Statement

The DASCombat Shiny operator is used for technical batch correction and currently relies solely on visual PCA inspection for assessment. The goal is to implement literature-based approaches to automatically determine:

1. Whether technical batch effects are present (so batch correction could be useful)
2. Whether batch correction was successful

## Current Operator Analysis

### Current Implementation
The operator uses a two-step PCA visualization approach:
- Calculates PCA on original data (before correction)
- Calculates PCA on corrected data (after correction)  
- Displays side-by-side PC1 vs PC2 plots colored by batch
- Relies on user visual inspection to determine correction effectiveness

### Current Limitations
- Assessment is purely subjective and visual
- No quantitative batch effect metrics
- biological groups (covariates of interest) are not taken into account
- Limited to PC1/PC2 visualization only
- No automated recommendations or thresholds
- No detection of overcorrection or failed corrections
- Difficult to document methodology decisions

## Context: Dataset Characteristics

**Feature Space**: 140-180 features (moderate dimensional, not high-dimensional) 
**Sample Size**: most of the times 12-36 samples total. 
**Technical batches**: runs (12 sample / run) and chips (4 samples / barcode) 
**Group (covariate of interest) Structure**: Typically 3 samples per condition, 2-12 conditions 
**Outlier Presence**: Large outliers are common due to small sample sizes


## Implications for Method Selection

**PCA Sufficiency**: With more features than samples, PCA will capture most variance in the first 2-3 components. Complex non-linear methods like t-SNE or UMAP become unnecessary and potentially counterproductive.

**K-means Inappropriateness**: With only 3 samples per group, k-means clustering becomes statistically meaningless. Any clustering approach requiring predetermined group numbers should be avoided.

**Outlier Impact**: Large outliers have disproportionate influence on all distance-based and clustering methods. Robust, median-based approaches become essential.

**Statistical Power Limitations**: Traditional mixed-effects models for PVCA become unreliable with such small sample sizes. Simpler ANOVA-based approaches are more appropriate.


## Revised Methodology for Small Datasets

### 1. Primary method: Simple Variance Decomposition
ANOVA on PC1 scores. This approach:
- Uses PC1 which captures most variance in small datasets
- Applies simple two-way ANOVA (PC1 ~ batch + biology) when biological groups are known
- Provides clear variance partitioning without problematic variance component estimation
- Remains statistically valid with small sample sizes

### 2. Robust Distance-Based Separation  
Median-based distance ratios:
- Calculates median within-batch vs between-batch distances in PC space
- Uses median statistics for outlier resistance
- Provides interpretable separation ratios without clustering assumptions
- Compares batch separation strength to biological separation when possible
- Excellent resistance to outliers through median-based calculations
- Provides quality control for primary method

### 3. Simple Geometric Overlap
Basic geometric overlap analysis:
- Examines range overlap between batches in PC1-PC2 space
- Provides intuitive geometric interpretation

### 4. Conservative Outlier Flagging
- Flag when single samples dramatically alter overall results
- Monitor for excessive variance reduction (>90%)
- Compare biological group coherence before/after when groups are known
- Provide warnings rather than automated decisions

### Keep Visual PCA Display
- Maintains user familiarity and confidence
- Allows expert override of automated recommendations  
- Provides visual validation of quantitative results
- Serves as backup when quantitative methods are uncertain

## Implementation Philosophy

**Conservative Approach**: Use multiple metrics but only provide strong recommendations when methods agree. Flag borderline cases for user decision rather than making potentially incorrect automated choices.

**Transparent Uncertainty**: Always display raw metric values alongside interpreted recommendations. Make method limitations explicit to users.

**Context Awareness**: Adjust interpretation based on dataset characteristics (sample size, balance, data type) rather than using universal thresholds.

**User Override**: Always allow expert users to override automated recommendations, especially important in research contexts where domain knowledge matters.

### Risk Mitigation

**Conservative Thresholds**: Set decision boundaries that err on the side of caution, flagging uncertain cases for user review rather than making potentially incorrect automated decisions.

**Multiple Method Validation**: Use complementary approaches to cross-validate findings. When methods disagree, recommend user review rather than automated decision.

**Explicit Limitations**: Clearly communicate method assumptions and limitations to users. Provide guidance on when automated recommendations might be unreliable.

**Preserving User Control**: Maintain visual PCA display and allow easy override of automated recommendations. The goal is decision support, not replacement of user judgment.

**Dataset-Specific Adaptation**: Adjust method parameters and thresholds based on detected dataset characteristics (sample size, feature count, balance) rather than using one-size-fits-all approaches.

This comprehensive approach provides objective, quantitative assessment while acknowledging the inherent challenges and limitations of batch effect detection in small, complex biological datasets.





### General Statistical Concerns
**Threshold Arbitrariness**: Many recommended thresholds (20% batch variance, 0.5 silhouette cutoff, etc.) lack strong empirical foundation and may not apply across different data types or experimental designs.

**Multiple Testing**: When applying multiple assessment methods, the chance of false positive recommendations increases without appropriate correction.

**Context Insensitivity**: Methods don't account for experimental context - sometimes apparent "batch effects" represent legitimate technical differences that shouldn't be removed.




