# CONSTELLATION — Method Notes

Working notes for the Methods section of the manuscript.

---

## Statistical Testing Framework

### Test Statistic

For a given ligand-receptor pair (L, R) in a spatial context (within a cell type or between a sender-receiver pair), the test statistic is:

$$T = \sum_i L_i \cdot w_i = \mathbf{L}^\top \mathbf{w}$$

where:
- $L_i = \log(1 + \text{ligand expression in cell } i)$
- $w_i = \sum_{j \in \mathcal{N}(i)} K(d_{ij}) \cdot R_j$ is the spatially-weighted receptor signal around cell $i$
- $K(d) = \exp(-d / \tau)$ is an exponential decay kernel with characteristic distance $\tau$
- $\mathcal{N}(i)$ are the spatial neighbors of cell $i$ (from KNN graph)

For **within-type** analysis, both $i$ and $j$ belong to the same cell type. For **between-type** analysis, $i$ indexes sender cells (expressing ligand) and $j$ indexes receiver cells (expressing receptor).

### Null Hypothesis

$H_0$: The spatial co-localization of ligand and receptor expression is no greater than expected by chance.

Under the null, we consider the distribution of $T$ when the ligand vector $\mathbf{L}$ is randomly permuted among the $n$ cells in the group, while the spatial structure (encoded in $\mathbf{w}$) is held fixed. This preserves the marginal distributions of both L and R, as well as the spatial graph topology, and only disrupts the spatial co-expression pattern.

### Analytical Permutation Null

Rather than performing Monte Carlo permutations, we compute the exact first two moments of $T$ under the permutation null in closed form.

For $T = \mathbf{L}_\pi^\top \mathbf{w}$ where $\pi$ is a random permutation of $\{1, \ldots, n\}$:

$$E[T] = \frac{\sum_i L_i \cdot \sum_i w_i}{n}$$

$$\text{Var}[T] = \frac{SS_L \cdot SS_w}{n - 1}$$

where:
- $SS_L = \sum_i L_i^2 - \frac{(\sum_i L_i)^2}{n}$ (sum of squared deviations for ligand)
- $SS_w = \sum_i w_i^2 - \frac{(\sum_i w_i)^2}{n}$ (sum of squared deviations for weighted receptor)

**Derivation sketch.** These follow from the moments of a sum of products under sampling without replacement. For any two fixed vectors $\mathbf{a}$ and $\mathbf{b}$, if $\pi$ is a uniform random permutation:

$$E\left[\sum_i a_i b_{\pi(i)}\right] = \frac{\sum a_i \cdot \sum b_i}{n}$$

$$\text{Var}\left[\sum_i a_i b_{\pi(i)}\right] = \frac{1}{n-1}\left(\sum a_i^2 - \frac{(\sum a_i)^2}{n}\right)\left(\sum b_i^2 - \frac{(\sum b_i)^2}{n}\right)$$

This is a classical result (see e.g., Hoeffding 1952; Lehmann 1975, "Nonparametrics: Statistical Methods Based on Ranks").

### Z-score and P-value

We compute a standardized z-score:

$$z = \frac{T_{\text{obs}} - E[T]}{\sqrt{\text{Var}[T]}}$$

and obtain a two-sided p-value from the standard normal distribution:

$$p = 2\left(1 - \Phi(|z|)\right)$$

**Normal approximation.** The permutation distribution of $T$ converges to a Gaussian under mild conditions (combinatorial CLT). This approximation is accurate when $n$ is not too small and neither $\mathbf{L}$ nor $\mathbf{w}$ is dominated by a few extreme values. We require $n \geq 50$ cells for within-type tests and $n_{\text{sender}} \geq 10$ for between-type tests.

**What is exact vs approximate:**
- The null mean $E[T]$ and variance $\text{Var}[T]$ are **exact** — they are the true moments of the finite permutation distribution, not estimates.
- The p-value is **approximate**, relying on the normal approximation to the permutation distribution. For $n \geq 50$, this approximation is very tight (validated empirically: Pearson $r > 0.999$ between analytical z-scores and 1000-permutation z-scores across 220 tests on real data with 67K–98K cells per group).

**Note:** Higher-order corrections are possible — the skewness and kurtosis of the permutation distribution also have closed forms (see Pesarin & Salmaso 2010) — but are unnecessary for the cell counts typical of spatial transcriptomics datasets.

### Validation

Comparison of analytical vs. Monte Carlo permutation (1000 permutations) on a 708,647-cell human lymph node Xenium dataset:

| Test type    | N tests | Pearson r (z-scores) | Mean |Δz| | Max |Δz| |
|-------------|---------|---------------------|-----------|----------|
| Within-type  | 100     | 0.9997              | 0.041     | 0.32     |
| Between-type | 120     | 0.9996              | 0.037     | 0.24     |

The small residual differences reflect Monte Carlo sampling noise in the permutation approach, not error in the analytical formula.

The analytical approach provides a ~400x speedup over 1000 permutations per test, enabling genome-wide testing across all cell-type pairs without computational bottleneck.

### Multiple Testing Correction

P-values across all tested (cell-type context, LR pair) combinations are corrected jointly using Benjamini-Hochberg FDR. The total testing burden includes both within-type and between-type tests, filtered by expression fraction thresholds (default: both ligand and receptor expressed in ≥5% of cells in the relevant population).

---

## Kernel Parameter $\tau$

The decay parameter $\tau$ controls the effective interaction range: weight $= \exp(-d/\tau)$, so interactions beyond $\sim 3\tau$ contribute negligibly.

**Default:** $\tau$ is set to the median nearest-neighbor distance (NN1) across all cells, which approximates the cell diameter. For typical Xenium data at subcellular resolution, this is ~5 µm, meaning interactions are weighted most strongly within ~15 µm (~3 cell diameters).

**Rationale:** Cell-cell signaling (ligand secretion → receptor binding) occurs predominantly between adjacent or nearby cells. Using the cell diameter as the characteristic scale captures direct contact and short-range paracrine signaling while down-weighting distant cells.

---

## Computational Optimizations

Two caching strategies reduce redundant computation:

1. **Receptor caching:** The weighted receptor vector $\mathbf{w}$ depends only on the receptor gene and the spatial context (cell type, graph). Since a receptor may pair with multiple ligands, we compute $\mathbf{w}$ once per unique receptor and reuse it across all paired ligands. The summary statistics $(sum(\mathbf{w}), SS_w)$ are also cached.

2. **Ligand caching:** Similarly, the ligand vector $\mathbf{L}$ and its summary statistics $(sum(\mathbf{L}), SS_L)$ are computed once per unique ligand and reused across all paired receptors.

This reduces the per-test cost to a single dot product plus O(1) arithmetic, with gene extraction amortized across all pairs sharing that gene.
