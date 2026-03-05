# CONSTELLATION Changelog

## [0.4.0] - 2025-02-10

### New Features
- **`run_lineage_analysis()`**: Organizes fine cell types into lineage groups for specificity assessment — testing whether an LR interaction is restricted to a particular subtype or broadly detected across related subtypes. Remaps via a `grouping` dict and delegates to `run_celltype_analysis()`. Compatible with ontology functions (`get_lymph_node_refined_groups()`, etc.).
- **Distance-weighted compartment test**: `run_compartment_analysis()` now uses the same kernel-weighted test statistic (`T = dot(L, w)`, `K = exp(-d/tau)`) as the cell-type level analysis. Replaces the previous binary neighbor test. Added `tau` parameter (default 5.0).
- **Analytical null for compartments**: Compartment tests use the same closed-form permutation null (`E[T] = sum(L)*sum(w)/n`, `Var[T] = SS_L*SS_w/(n-1)`) as cell-type tests. No permutations needed.
- **Within-compartment edges**: Compartment analysis restricts KNN edges to pairs where both endpoints belong to the same compartment, preventing cross-compartment signal leakage.

### Code Quality
- **Extracted `_analytical_test()` helper**: Consolidates the analytical permutation null computation (3 inline copies → 1 function) used by within-type, between-type, and compartment tests.
- **Extracted `_apply_burden_correction()` helper**: Consolidates FDR correction with burden padding (5 inline copies → 1 function). Handles `fdr_method='none'` without crashing.
- **Bug fix: `fdr_method='none'`**: Previously crashed with `ValueError` in `multipletests`. Now correctly skips correction and passes through raw p-values.
- **Removed unused `ct_sets` parameter**: Accepted by `test_between_type_lr` and targeted testing functions but never used. Removed from all signatures.
- **Parameterized filtering thresholds**: `min_cells` and `min_edges` in `test_within_type_lr` (default 50/50) and `test_between_type_lr` (default 10/100) are now explicit parameters instead of hardcoded.
- **`apply_fdr_correction()` returns a copy**: No longer modifies the input DataFrame in place.
- **`load_cell_type_mapping()` parameterized**: Accepts `source_col` and `target_col` instead of hardcoded column names.
- **Standardized spatial indexing**: `run_compartment_analysis` now uses `build_spatial_graph_from_coords` instead of direct `cKDTree`.

### Documentation
- **README.md**: Comprehensive rewrite with method overview, usage examples for cell-type, compartment, and lineage analysis, input requirements, and parameter guidance.
- **API.md**: Complete function reference for all 70+ exported functions across 4 modules.

## [0.3.0] - 2025-01-28

### Performance Improvements
- **Vectorized edge counting**: Replaced nested Python loops with NumPy boolean mask operations in `compute_testable_burden_within` and `compute_testable_burden_between`. Results in 16-19x speedup for burden computation.
- **Pre-computed cell type masks**: Between-type burden computation now pre-computes boolean masks for all cell types, avoiding redundant array allocations.

### New Features
- **Progress bars**: Added tqdm progress bars to all long-running functions:
  - `compute_testable_burden_within`
  - `compute_testable_burden_between`
  - `test_lr_pair_within_types`
  - `test_lr_pair_between_types`
  - `test_ligand_all_receptors`
  - `test_receptor_all_ligands`
  - `test_custom_lr_set`

### Dependencies
- Added `tqdm` as a required dependency

## [0.2.0] - 2025-01-27

### New Features
- **Targeted testing functions**: Added functions for testing specific LR pairs with full burden FDR correction:
  - `test_lr_pair_within_types`
  - `test_lr_pair_between_types`
  - `test_ligand_all_receptors`
  - `test_receptor_all_ligands`
  - `test_custom_lr_set`

- **Testable burden computation**: Added functions to pre-compute multiple testing burden:
  - `compute_testable_burden_within`
  - `compute_testable_burden_between`
  - `save_testable_burden`
  - `load_testable_burden`

- **Cell size statistics**: Added functions for spatial cell size analysis:
  - `compute_cell_sizes`
  - `compute_cell_sizes_with_area`
  - `summarize_cell_sizes`
  - `print_cell_size_table`

### Performance Improvements
- **Ligand grouping optimization**: LR pairs sharing the same ligand now reuse the permuted ligand matrix, reducing redundant computation.
- **CSC matrix support**: Sparse expression matrices are converted to CSC format for efficient column slicing.

## [0.1.0] - 2025-01-26

### Initial Release
- Core spatial LR testing functions:
  - `test_within_type_lr`: Test LR pairs within a single cell type (autocrine)
  - `test_between_type_lr`: Test LR pairs between cell types (paracrine)
- I/O functions for loading data and saving results
- Cell Ontology integration for hierarchical cell type analysis
- Visualization functions for results
