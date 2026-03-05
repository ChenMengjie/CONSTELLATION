---
title: Methods
layout: default
nav_order: 3
---

# Statistical Methods
{: .no_toc }

<details open markdown="block">
  <summary>Table of contents</summary>
  {: .text-delta }
- TOC
{:toc}
</details>

---

## Test statistic

For a given ligand-receptor pair (L, R) in a spatial context (within a cell type or between a sender-receiver pair), the test statistic is:

$$T = \sum_i L_i \cdot w_i = \mathbf{L}^\top \mathbf{w}$$

where:
- $$L_i = \log(1 + \text{ligand expression in cell } i)$$
- $$w_i = \sum_{j \in \mathcal{N}(i)} K(d_{ij}) \cdot R_j$$ is the spatially-weighted receptor signal around cell $$i$$
- $$K(d) = \exp(-d / \tau)$$ is an exponential decay kernel with characteristic distance $$\tau$$
- $$\mathcal{N}(i)$$ are the spatial neighbors of cell $$i$$ (from KNN graph)

For **within-type** analysis, both $$i$$ and $$j$$ belong to the same cell type. For **between-type** analysis, $$i$$ indexes sender cells (expressing ligand) and $$j$$ indexes receiver cells (expressing receptor).

## Null hypothesis

$$H_0$$: The spatial co-localization of ligand and receptor expression is no greater than expected by chance.

Under the null, we consider the distribution of $$T$$ when the ligand vector $$\mathbf{L}$$ is randomly permuted among the $$n$$ cells in the group, while the spatial structure (encoded in $$\mathbf{w}$$) is held fixed. This preserves the marginal distributions of both L and R, as well as the spatial graph topology, and only disrupts the spatial co-expression pattern.

## Analytical permutation null

Rather than performing Monte Carlo permutations, we compute the exact first two moments of $$T$$ under the permutation null in closed form.

For $$T = \mathbf{L}_\pi^\top \mathbf{w}$$ where $$\pi$$ is a random permutation of $$\{1, \ldots, n\}$$:

$$E[T] = \frac{\sum_i L_i \cdot \sum_i w_i}{n}$$

$$\text{Var}[T] = \frac{SS_L \cdot SS_w}{n - 1}$$

where:
- $$SS_L = \sum_i L_i^2 - \frac{(\sum_i L_i)^2}{n}$$ (sum of squared deviations for ligand)
- $$SS_w = \sum_i w_i^2 - \frac{(\sum_i w_i)^2}{n}$$ (sum of squared deviations for weighted receptor)

These follow from the moments of a sum of products under sampling without replacement (Hoeffding, 1952; Lehmann, 1975).

### What is exact vs approximate

- The null mean $$E[T]$$ and variance $$\text{Var}[T]$$ are **exact** --- they are the true moments of the finite permutation distribution, not estimates.
- The p-value is **approximate**, relying on the normal approximation to the permutation distribution. For $$n \geq 50$$, this approximation is very tight (validated empirically: Pearson $$r > 0.999$$ between analytical z-scores and 1000-permutation z-scores across 220 tests on real data with 67K--98K cells per group).

## Z-score and p-value

We compute a standardized z-score:

$$z = \frac{T_{\text{obs}} - E[T]}{\sqrt{\text{Var}[T]}}$$

and obtain a two-sided p-value from the standard normal distribution:

$$p = 2\left(1 - \Phi(|z|)\right)$$

The analytical approach provides a ~400x speedup over 1000-permutation Monte Carlo, enabling genome-wide testing across all cell-type pairs without computational bottleneck.

## Fold enrichment

$$\text{FE} = \frac{T_{\text{obs}}}{E[T]}$$

A result is called significant if $$p_{\text{adj}} < 0.05$$ and $$z > 0$$ (i.e., $$\text{FE} > 1$$).

## Multiple testing correction

P-values across all tested (cell-type context, LR pair) combinations are corrected jointly using Benjamini-Hochberg FDR. The total testing burden includes both within-type and between-type tests, filtered by expression fraction thresholds (default: both ligand and receptor expressed in $$\geq 5\%$$ of cells in the relevant population).

When the total number of testable combinations exceeds the number actually performed (e.g., targeted testing of a single LR pair), dummy p-values of 1.0 are appended before BH correction to maintain the full testing burden.

## Kernel parameter tau

The decay parameter $$\tau$$ controls the effective interaction range: weight $$= \exp(-d/\tau)$$, so interactions beyond $$\sim 3\tau$$ contribute negligibly.

**Default:** $$\tau$$ is set to the median nearest-neighbor distance (NN1) across all cells, which approximates the cell diameter. For typical Xenium data at subcellular resolution, this is ~5 um, meaning interactions are weighted most strongly within ~15 um (~3 cell diameters).

**Rationale:** Cell-cell signaling (ligand secretion to receptor binding) occurs predominantly between adjacent or nearby cells. Using the cell diameter as the characteristic scale captures direct contact and short-range paracrine signaling while down-weighting distant cells.

**Spatial range characterization:** Sweeping $$\tau$$ across a range of values classifies LR interactions as contact-dependent (signal decays as $$\tau$$ increases beyond the contact range) or secreted (signal remains stable across $$\tau$$ values, reflecting diffusion-mediated signaling).

## Minimum cell and edge requirements

| Analysis type | Min cells | Min edges | Min expression fraction |
|--------------|-----------|-----------|------------------------|
| Within-type (autocrine) | 50 cells | 50 within-type edges | 5% of cells in type |
| Between-type (paracrine) | 10 sender + 10 receiver | 100 cross-type edges | 5% of sender (ligand) or receiver (receptor) |

## Validation

Comparison of analytical vs. Monte Carlo permutation (1000 permutations) on a 708,647-cell human lymph node Xenium dataset:

| Test type | N tests | Pearson r (z-scores) | Mean \|delta z\| | Max \|delta z\| |
|-----------|---------|---------------------|-----------|----------|
| Within-type | 100 | 0.9997 | 0.041 | 0.32 |
| Between-type | 120 | 0.9996 | 0.037 | 0.24 |

The small residual differences reflect Monte Carlo sampling noise in the permutation approach, not error in the analytical formula.
