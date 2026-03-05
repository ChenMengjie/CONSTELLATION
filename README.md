# CONSTELLATION

CONSTELLATION is a spatial-conditioned permutation testing framework for ligand-receptor (LR) interaction inference in spatial transcriptomics data. It tests whether LR co-expression is spatially enriched beyond what is expected under random cell labeling, using an analytical permutation null that provides exact p-values without Monte Carlo sampling.

## Key features

- **Analytical permutation null** — exact mean and variance of the test statistic under permutation, ~400x faster than Monte Carlo
- **Distance-weighted kernel** — exponential decay `K(d) = exp(-d/tau)` captures spatial range of interactions
- **Cell-type and compartment analysis** — test LR pairs within/between cell types, or across spatial compartments
- **Spatial range characterization** — sweep over tau to classify contact-dependent vs secreted signals
- **Built-in LR databases** — CellPhoneDB, CellTalkDB, ConnectomeDB, LIANA consensus (~4,700 pairs)
- **Cell Ontology integration** — hierarchical cell type grouping via EBI OLS API

## Installation

```bash
git clone https://github.com/ChenMengjie/CONSTELLATION.git
cd constellation
pip install -e .
```

**Requirements:** Python >= 3.9, numpy, scipy, pandas, scanpy, scikit-learn, matplotlib, seaborn, statsmodels, requests, tqdm.

## Quick start

### Cell-type level analysis

```python
import constellation as cst
import scanpy as sc
import numpy as np

# Load data
adata = sc.read_h5ad('data.h5ad')

# IMPORTANT: CONSTELLATION applies log1p internally.
# If your data is already log-normalized, reverse it:
adata.X = np.expm1(adata.X.toarray())

# Load LR database
lr_df, lr_pairs = cst.load_lr_resource('consensus', adata=adata)

# Build spatial graph
indices, distances = cst.build_spatial_graph(adata, k=10)

# Run full analysis (within-type + between-type, FDR correction)
results = cst.run_celltype_analysis(
    adata, cell_type_col='cell_type',
    lr_pairs=lr_pairs,
    indices=indices, distances=distances,
    tau=5.0,
)

# View results
cst.print_top_results(results, n=20)
```

### Compartment-level analysis

```python
import constellation as cst

# Prepare expression dictionary (raw counts, not log-transformed)
expr_dict = {gene: expression_vector for gene, expression_vector in ...}

# Run compartment analysis
results = cst.run_compartment_analysis(
    expr_dict=expr_dict,
    compartments=compartments,       # array of compartment labels
    coords=coords,                   # (n_cells, 2) coordinates
    lr_pairs=lr_pairs,
    k=10,
    tau=7.0,
)
```

### Lineage-level analysis

CONSTELLATION tests at the finest cell type resolution by default, which enables **specificity assessment** — determining whether an LR interaction is restricted to a particular subtype or broadly detected across related subtypes. The `run_lineage_analysis()` function organizes fine cell types into lineage groups via a user-supplied grouping dictionary (e.g., informed by Cell Ontology), enabling systematic comparison of results across subtypes within each lineage.

```python
import constellation as cst
from constellation.ontology import get_lymph_node_refined_groups

# Group fine cell types into lineages for specificity assessment
grouping = get_lymph_node_refined_groups()
# e.g. {"CD4_T": ["Naive_CD4_T", "Treg", "Tfh_GC"], "B_cell": ["Naive_B", "GC_DZ", "GC_LZ"], ...}

results = cst.run_lineage_analysis(
    adata, cell_type_col='cell_type',
    grouping=grouping,
    lr_pairs=lr_pairs,
    indices=indices, distances=distances,
    tau=5.0,
)
```

### Targeted testing

```python
# Test a specific LR pair across all cell type contexts
results = cst.test_lr_pair_between_types(
    ligand='CXCL13', receptor='CXCR5',
    cell_types_array=cell_types,
    X=adata.X, gene_names=gene_names,
    indices=indices, distances=distances,
    tau=5.0,
)

# Test all receptors for a given ligand
results = cst.test_ligand_all_receptors(
    ligand='CD40LG', lr_pairs=lr_pairs,
    cell_types_array=cell_types,
    X=adata.X, gene_names=gene_names,
    indices=indices, distances=distances,
)
```

## Method

### Test statistic

For a given LR pair in a spatial context, the test statistic is:

```
T = L^T * w
```

where:
- `L_i = log(1 + ligand expression in cell i)`
- `w_i = sum_j K(d_ij) * R_j` is the spatially-weighted receptor signal around cell `i`
- `K(d) = exp(-d / tau)` is the exponential decay kernel
- `N(i)` are the spatial k-nearest neighbors of cell `i`

### Analytical permutation null

Under the null hypothesis (random permutation of L while holding the spatial structure w fixed):

```
E[T] = sum(L) * sum(w) / n
Var[T] = SS_L * SS_w / (n - 1)
```

where `SS_L` and `SS_w` are sums of squared deviations. The z-score `z = (T - E[T]) / sqrt(Var[T])` follows a standard normal distribution for `n >= 50` cells (validated: Pearson r > 0.999 vs 1000-permutation Monte Carlo on 708K-cell data).

### Multiple testing correction

P-values are corrected using Benjamini-Hochberg FDR across all tested (context, LR pair) combinations. Expression fraction thresholds (default >= 5%) filter untestable pairs before correction.

### Kernel parameter tau

`tau` controls the effective interaction range. Default: median nearest-neighbor distance (~5 um for Xenium). Interactions beyond ~3*tau contribute negligibly. Sweeping tau classifies contact-dependent (signal decays with tau) vs secreted (signal stable across tau) interactions.

## Modules

| Module | Description |
|--------|-------------|
| `constellation.testing` | Core statistical testing: within-type, between-type, compartment, lineage |
| `constellation.io` | Data loading, spatial graph, LR databases, FDR correction |
| `constellation.visualization` | Heatmaps, dot plots, networks, spatial maps, compartment plots |
| `constellation.ontology` | Cell Ontology integration, hierarchical cell type grouping |

See [API.md](API.md) for the complete function reference.

## Input data requirements

| Input | Format | Notes |
|-------|--------|-------|
| Expression matrix | AnnData `.X` or `expr_dict` | **Raw or normalized counts, NOT log-transformed** |
| Cell coordinates | `adata.obsm['spatial']` or numpy array | x, y in microns |
| Cell type labels | `adata.obs[cell_type_col]` | Categorical column |
| LR pairs | List of `(ligand, receptor)` tuples | From `load_lr_resource()` |

CONSTELLATION applies `log1p()` internally. If your data is already log-normalized (e.g., from `sc.pp.log1p()`), you must reverse it with `np.expm1()` before passing to CONSTELLATION.

## Tutorials

Step-by-step notebooks demonstrating CONSTELLATION on three Xenium spatial transcriptomics datasets.

| Dataset | Notebook | Description |
|---------|----------|-------------|
| **Lymph Node** (Xenium Prime 5K, 708K cells, 22 types) | [`tutorial_lymphnode.ipynb`](tutorials/lymph_node/tutorial_lymphnode.ipynb) | Spatial range characterization, tau sweep, multi-method benchmarking (CD40LG–CD40, CXCL13–CXCR5) |
| **Colon Cancer** (Xenium V1, 308K cells, 3 compartments) | [`tutorial_crc_constellation.ipynb`](tutorials/crc/tutorial_crc_constellation.ipynb) | Per-compartment analysis, tumor microenvironment LR characterization |
| **Ovarian Cancer** (Xenium Prime 5K, 265K cells, 17 types) | [`tutorial_oc.ipynb`](tutorials/oc/tutorial_oc.ipynb) | Per-compartment analysis, 5-method benchmarking, CXCL12–CXCR4 specificity |

Data preprocessing and cell type annotation notebooks: [lymph node](tutorials/lymph_node/tutorial_annotate_lymphnode.ipynb) | [colon cancer](tutorials/crc/tutorial_annotate_compartments.ipynb) | [ovarian cancer](tutorials/oc/tutorial_annotate_oc.ipynb)

## Citation

If you use CONSTELLATION in your research, please cite:

> [Citation to be added upon publication]

## License

MIT
