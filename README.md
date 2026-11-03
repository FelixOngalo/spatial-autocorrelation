# Spatial Autocorrelation in R: Mathematical Foundations, Implementation, and Statistical Verification

## Project Overview

This project provides a detailed investigation of **spatial autocorrelation**, focusing on the mathematical principles, statistical methods, and computational techniques used to identify spatial patterns in geographic and spatially indexed data. It examines whether observations located near one another tend to exhibit similar values, contrasting values, or spatially random distributions.

The project implements **Global Moran’s I and Local Moran’s I from scratch in R**, translating their mathematical definitions into computational procedures rather than relying exclusively on existing statistical packages. It explores the construction and interpretation of spatial weights matrices, the role of neighbourhood relationships, and the distinction between global spatial dependence and local spatial clustering.

To establish statistical significance, the analysis incorporates **permutation-based hypothesis testing**. This approach generates reference distributions by repeatedly rearranging observed values across spatial locations, allowing the observed Moran statistics to be evaluated against what might be expected under spatial randomness.

The project also addresses the multiple-testing problem associated with local spatial statistics. Since Local Moran’s I produces a separate statistical test for each observation, conducting many tests simultaneously increases the likelihood of false-positive findings. The analysis therefore applies **False Discovery Rate (FDR) correction** to control the expected proportion of false discoveries among the identified significant locations.

Finally, the statistical procedures are evaluated using simulated spatial datasets with known underlying hotspot and coldspot structures. By comparing detected patterns with the predefined ground truth, the project assesses the accuracy of the methods and demonstrates how mathematical statistics can be used to distinguish genuine spatial clusters from patterns that may arise by chance.

## 1. Research Objectives

The primary objective is to understand, implement, and evaluate statistical methods for measuring spatial autocorrelation and identifying statistically significant local spatial patterns.

Specific objectives include:

- Establish the mathematical foundations of spatial autocorrelation and spatial dependence.
- Construct spatial weights matrices to represent relationships between neighbouring observations.
- Derive and implement Global Moran’s I using its mathematical definition.
- Implement Local Moran’s I to identify local clusters and spatial outliers.
- Develop permutation-based procedures for assessing statistical significance.
- Apply FDR correction to account for multiple simultaneous hypothesis tests.
- Identify spatial hotspots, coldspots, and spatial outliers.
- Evaluate statistical performance against simulated data with known spatial structures.
- Examine the limitations of spatial statistical inference and the sensitivity of results to neighbourhood definitions.

## 2. Mathematical Foundations

### 2.1 Spatial Autocorrelation

Spatial autocorrelation describes the degree to which observations at different geographic locations are statistically related according to their spatial proximity or neighbourhood relationships.

Three broad patterns can occur:

- **Positive spatial autocorrelation:** Nearby observations tend to have similar values.
- **Negative spatial autocorrelation:** Nearby observations tend to have dissimilar values.
- **Spatial randomness:** No systematic spatial relationship is detected.

These patterns are important because conventional statistical methods often assume that observations are independent. Spatial dependence can violate this assumption and affect statistical inference, standard errors, and model interpretation.

### 2.2 Spatial Weights Matrices

Spatial relationships are represented mathematically through a spatial weights matrix:

\[
W=(w_{ij})
\]

where \(w_{ij}\) describes the spatial relationship between observations \(i\) and \(j\).

The matrix determines which observations are considered neighbours and how strongly their values contribute to one another's spatial statistics.

Depending on the implementation, neighbourhood relationships may be defined using:

- Contiguity-based neighbourhoods.
- Distance-based neighbourhoods.
- Fixed-distance thresholds.
- K-nearest-neighbour relationships.

The project examines how the specification of spatial weights influences the resulting measures of spatial autocorrelation.

### 2.3 Global Moran’s I

Global Moran’s I measures the overall degree of spatial autocorrelation across a dataset.

Its mathematical expression is:

\[
I=\frac{n}{S_0}
\frac{\sum_i\sum_j w_{ij}(x_i-\bar{x})(x_j-\bar{x})}
{\sum_i(x_i-\bar{x})^2}
\]

where:

- \(n\) is the number of observations.
- \(x_i\) and \(x_j\) are observed values.
- \(\bar{x}\) is the sample mean.
- \(w_{ij}\) represents the spatial weight between observations.
- \(S_0=\sum_i\sum_j w_{ij}\) is the sum of all spatial weights.

The statistic provides an overall measure of spatial dependence. Positive values generally indicate spatial clustering of similar values, negative values suggest spatial dispersion of dissimilar values, and values near the expected value under the null hypothesis indicate limited evidence of global spatial autocorrelation.

### 2.4 Local Moran’s I

While Global Moran’s I describes the overall spatial structure, Local Moran’s I identifies where spatial clustering or spatial outliers occur.

A commonly used formulation is:

\[
I_i=\frac{x_i-\bar{x}}{m_2}
\sum_j w_{ij}(x_j-\bar{x})
\]

where:

\[
m_2=\frac{1}{n}\sum_i(x_i-\bar{x})^2
\]

The local statistic identifies four important spatial configurations:

- **High–High (HH):** A high-value observation surrounded by high-value neighbours.
- **Low–Low (LL):** A low-value observation surrounded by low-value neighbours.
- **High–Low (HL):** A high-value observation surrounded by low-value neighbours.
- **Low–High (LH):** A low-value observation surrounded by high-value neighbours.

High–High and Low–Low configurations represent spatial clusters, while High–Low and Low–High configurations indicate potential spatial outliers.

## 3. Computational Implementation in R

The statistical procedures are implemented from their mathematical definitions in R.

The implementation focuses on:

1. Preparing spatially indexed observations.
2. Constructing the spatial weights matrix.
3. Calculating Global Moran’s I.
4. Calculating Local Moran’s I for individual observations.
5. Generating permutation distributions.
6. Calculating permutation-based significance levels.
7. Applying FDR correction.
8. Classifying and visualising significant spatial patterns.

This approach makes the computational steps transparent and provides an opportunity to examine the relationship between mathematical definitions and their numerical implementation.

## 4. Permutation-Based Hypothesis Testing

The observed Moran statistics alone do not establish whether spatial autocorrelation is statistically significant.

The project therefore uses permutation testing to construct an empirical reference distribution under the null hypothesis of spatial randomness.

The procedure involves:

1. Calculating the observed Global or Local Moran’s I statistic.
2. Randomly permuting the observed values across spatial locations.
3. Recalculating the statistic for each permutation.
4. Repeating the procedure to generate an empirical null distribution.
5. Comparing the observed statistic with the permutation distribution.
6. Estimating a permutation-based p-value.

For a two-sided permutation test, a commonly used Monte Carlo p-value is:

\[
p=\frac{1+\sum_{b=1}^{B}
\mathbf{1}(|I_b|\geq |I_{\mathrm{obs}}|)}
{B+1}
\]

where:

- \(B\) is the number of permutations.
- \(I_b\) is the statistic obtained from permutation \(b\).
- \(I_{\mathrm{obs}}\) is the observed statistic.
- \(\mathbf{1}(\cdot)\) is an indicator function.

The addition of one to both the numerator and denominator avoids reporting a zero Monte Carlo p-value.

## 5. Multiple Testing and False Discovery Rate Correction

Local Moran’s I produces a separate statistic and significance test for every spatial observation. Consequently, testing many locations simultaneously increases the probability of identifying false-positive hotspots or coldspots.

To address this issue, the project implements **False Discovery Rate correction**, using the Benjamini–Hochberg procedure.

Given ordered p-values:

\[
p_{(1)}\leq p_{(2)}\leq \cdots \leq p_{(m)}
\]

the procedure compares each ordered p-value with:

\[
\frac{i}{m}q
\]

where:

- \(m\) is the number of simultaneous tests.
- \(i\) is the rank of the ordered p-value.
- \(q\) is the selected FDR significance level.

This approach controls the expected proportion of false discoveries among rejected hypotheses under its applicable assumptions.

The corrected significance values are used to distinguish statistically significant local spatial patterns from observations that do not provide sufficient evidence of spatial clustering.

## 6. Hotspot and Coldspot Detection

After calculating Local Moran’s I and applying FDR correction, the project classifies observations according to their local spatial configurations.

The classification process combines:

- The sign of the centred observation.
- The sign of its spatially lagged value.
- The statistical significance of the local Moran statistic.
- The adjusted p-value after FDR correction.

This produces interpretable spatial categories:

| Classification | Interpretation |
|---|---|
| High–High | Statistically significant hotspot |
| Low–Low | Statistically significant coldspot |
| High–Low | Potential high-value spatial outlier |
| Low–High | Potential low-value spatial outlier |
| Not significant | No statistically significant local association detected |

The distinction between a spatial cluster and a spatial outlier is important. A high-valued observation is not automatically a hotspot; it must also exhibit an appropriate relationship with its neighbours and satisfy the statistical significance criterion.

## 7. Simulation and Validation Against Known Ground Truth

An important component of the project is the use of simulated spatial data with predefined spatial structures.

Unlike real-world datasets, where the true underlying hotspot locations are generally unknown, simulated datasets allow the actual spatial configurations to be specified before statistical analysis.

The simulated data provide known examples of:

- High-value spatial clusters.
- Low-value spatial clusters.
- Spatially random observations.
- Potential spatial outliers.

The implemented procedures are then evaluated by comparing their detected classifications against the known underlying patterns.

This provides a controlled environment for investigating whether the statistical methods correctly identify spatial clusters and how often they incorrectly classify locations.

### Evaluation considerations

The comparison against ground truth can be used to examine:

- **True positives:** Correctly identified hotspots or coldspots.
- **False positives:** Locations incorrectly identified as significant clusters.
- **True negatives:** Locations correctly identified as non-clusters.
- **False negatives:** Genuine clusters that the procedure fails to detect.

Additional performance measures may include:

\[
\text{Precision}=
\frac{TP}{TP+FP}
\]

\[
\text{Recall}=
\frac{TP}{TP+FN}
\]

\[
\text{F1 Score}=
2\left(\frac{\text{Precision}\times\text{Recall}}
{\text{Precision}+\text{Recall}}\right)
\]

These measures provide a more complete assessment of detection performance than statistical significance alone.

## 8. Visualisation and Interpretation

Visualisation supports the interpretation of spatial dependence and local clustering.

The project uses graphical representations to examine:

- The distribution of simulated spatial observations.
- Global spatial autocorrelation results.
- Local Moran’s I statistics.
- Permutation reference distributions.
- Significant High–High and Low–Low clusters.
- Spatial outliers.
- Detected classifications compared with the known ground truth.

These visualisations help connect numerical statistical results to their spatial interpretation.

## 9. Limitations and Methodological Considerations

Several considerations are important when interpreting the results:

- **Spatial weights specification:** Different neighbourhood definitions can produce different spatial autocorrelation results.
- **Permutation uncertainty:** Monte Carlo p-values depend on the number of permutations. A larger number generally provides finer p-value resolution but requires more computation.
- **Multiple testing:** FDR correction reduces the expected proportion of false discoveries but does not eliminate all false positives.
- **Spatial dependence:** Local tests may be dependent on one another because neighbouring locations share observations.
- **Simulation assumptions:** Results obtained from simulated spatial structures may not generalise directly to complex real-world geographic processes.
- **Global versus local inference:** A significant global statistic does not imply that every individual location contains a significant local cluster.
- **Implementation verification:** From-scratch implementations should be checked against mathematical properties and, where appropriate, established statistical software.

## 10. Learning Outcomes

Through this project, I explore the intersection of mathematical statistics, computational implementation, and spatial data analysis.

The project demonstrates the ability to:

- Translate mathematical statistical definitions into computational algorithms.
- Implement spatial statistics from first principles.
- Construct and interpret spatial weights matrices.
- Conduct Monte Carlo permutation tests.
- Apply multiple-testing corrections.
- Distinguish global spatial dependence from local spatial association.
- Evaluate statistical procedures against known simulated outcomes.
- Interpret the limitations of statistical inference in spatially dependent data.

## Conclusion

This project demonstrates how mathematical statistics can be translated into practical computational methods for understanding spatial relationships. By implementing Global and Local Moran’s I from scratch, developing permutation-based hypothesis tests, and applying FDR correction, it provides a structured examination of spatial autocorrelation and hotspot detection.

The use of simulated data with known ground truth adds an important verification component, making it possible to assess not only whether spatial patterns are statistically significant but also whether the implemented methods correctly recover the underlying spatial structures.

Overall, the project combines mathematical reasoning, statistical inference, algorithm development, and reproducible R programming to investigate one of the fundamental problems in spatial data analysis.
