# MAGICAL Methodology and Circuit Interpretation

`pymagical` is based on the MAGICAL (Multiome Accessibility Gene Integration Calling and Looping) framework described in *Chen et al., Nature Computational Science 2023*. This methodology leverages a hierarchical Bayesian approach to map regulatory circuits consisting of Transcription Factors (TFs), cis-regulatory chromatin sites (Peaks), and target Genes.

## Hierarchical Bayesian Model

MAGICAL models the coordinated variation in chromatin accessibility and gene expression across conditions and samples. It introduces hidden variables to explicitly model signal and noise in both data types.

### 1. Chromatin Accessibility Model
The chromatin activity $A$ for a peak $p$ in a cell $k$ is modeled as:
$$A_{pk} = \sum_{m} B_{pm} T_{mk} + N_{A,pk}$$
*   **$B_{pm}$:** TF–peak binding confidence (Weight) between TF $m$ and peak $p$.
*   **$T_{mk}$:** Hidden TF activity of TF $m$ in cell $k$. This represents the regulatory capacity (protein level) and is distinct from TF expression.
*   **$N_A$:** Data noise in the ATAC-seq modality.

### 2. Gene Expression Model
The gene expression $R$ for a gene $g$ in a cell $k$ is modeled as:
$$R_{gk} = \sum_{p} L_{gp} \left( \sum_{m} B_{pm} T_{mk} \right) + N_{R,gk}$$
*   **$L_{gp}$:** Peak–gene looping confidence (Weight) between peak $p$ and gene $g$.
*   **$N_R$:** Data noise in the RNA-seq modality.

## Continuous Weights vs. Posterior Probabilities

During Gibbs sampling, MAGICAL iteratively estimates two types of variables for the linkages ($B$ and $L$):

### 1. Continuous Weights (Regression Coefficients)
The weights $B_{pm}$ and $L_{gp}$ are sampled from zero-mean Gaussian distributions. They represent the linear effect size and can be positive or negative.
*   **Positive Weight:** Indicates a positive correlation (e.g., TF opens the peak, or open peak increases gene expression).
*   **Negative Weight:** Indicates a negative correlation (e.g., TF closes the peak, or open peak decreases gene expression).

### 2. Binary States and Sampling Frequency
To account for the uncertainty of whether a regulatory link actually exists, the model uses binary indicator variables ($B_{state}$ and $L_{state}$). In each iteration, a state is drawn (0 or 1) based on the posterior probability. 

The **final posterior probability** reported in the output is the **sampling frequency** of the binary state being 1 across all iterations:
$$P(state = 1) = \frac{\sum_{n=1}^N state^{(n)}}{N}$$

## Biological Directionality (Activators vs. Repressors)

While the original paper focuses on the existence of circuits (the posterior probabilities), `pymagical` extends this by interpreting the biological directionality using the continuous weights.

The overall regulatory effect of a TF on a target gene through a specific peak is determined by the product of the average weights:
$$\text{Overall Effect Sign} = \text{sign}(Average(B) \times Average(L))$$

### Circuit Classification
*   **Activator (+):** The TF has a positive overall effect on gene expression.
    *   `[+,+]`: TF opens peak (+), and open peak increases gene expression (+).
    *   `[-,-]`: TF closes peak (-), and closed peak increases gene expression (meaning open peak decreases it, -).
*   **Repressor (-):** The TF has a negative overall effect on gene expression.
    *   `[+,-]`: TF opens peak (+), but open peak decreases gene expression (-).
    *   `[-,+]`: TF closes peak (-), and open peak would have increased gene expression (+).

## Output Notation

`pymagical` uses the following notation in the results edge list:
`TF_Name (Confidence_Probability, Overall_Effect [L_dir, B_dir])`

*   **Confidence_Probability:** The sampling frequency of the circuit (0.0 to 1.0).
*   **Overall_Effect:** `+` for Activation, `-` for Repression.
*   **L_dir:** Direction of Peak $\to$ Gene linkage weight.
*   **B_dir:** Direction of TF $\to$ Peak linkage weight.

*Example:* `KLF5 (0.90, - [-,+])`
This indicates a 90% confidence circuit where KLF5 opens the peak (`B_dir = +`), but the peak's accessibility is negatively correlated with the target gene's expression (`L_dir = -`), resulting in overall repression.
