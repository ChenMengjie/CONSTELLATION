# CONSTELLATION API Reference

## Table of contents

- [High-level analysis](#high-level-analysis)
- [Compartment-level analysis](#compartment-level-analysis)
- [Batch testing](#batch-testing)
- [Targeted testing](#targeted-testing)
- [Data loading and I/O](#data-loading-and-io)
- [Spatial graph](#spatial-graph)
- [LR databases](#lr-databases)
- [FDR and burden](#fdr-and-burden)
- [Cell size statistics](#cell-size-statistics)
- [Visualization — cell-type interactions](#visualization--cell-type-interactions)
- [Visualization — LR-level plots](#visualization--lr-level-plots)
- [Visualization — compartment plots](#visualization--compartment-plots)
- [Cell Ontology](#cell-ontology)

---

## High-level analysis

### `run_celltype_analysis`

```python
cst.run_celltype_analysis(
    adata, cell_type_col, lr_pairs, indices, distances,
    expr_frac=None, tau=5.0, min_expr_frac=0.05,
    min_cells=50, fdr_method='fdr_bh', fdr_threshold=0.05,
    test_within=True, test_between=True, X_csc=None,
    lr_annotations=None, verbose=True,
)
```

Run full CONSTELLATION analysis: filter testable pairs, run within-type and between-type tests, apply FDR correction, and return a combined results DataFrame.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `adata` | AnnData | — | Annotated data matrix (raw counts, NOT log-transformed) |
| `cell_type_col` | str | — | Column in `adata.obs` with cell type labels |
| `lr_pairs` | list[tuple] | — | List of `(ligand, receptor)` tuples |
| `indices` | ndarray | — | KNN indices from `build_spatial_graph()` |
| `distances` | ndarray | — | KNN distances from `build_spatial_graph()` |
| `expr_frac` | dict | None | Pre-computed expression fractions (from `compute_expression_fractions`) |
| `tau` | float | 5.0 | Distance decay parameter (um) for kernel `K(d) = exp(-d/tau)` |
| `min_expr_frac` | float | 0.05 | Minimum fraction of cells expressing a gene to include |
| `min_cells` | int | 50 | Minimum cells in a type for within-type testing |
| `fdr_method` | str | 'fdr_bh' | Multiple testing correction method |
| `fdr_threshold` | float | 0.05 | FDR significance threshold |
| `test_within` | bool | True | Run within-type (autocrine) tests |
| `test_between` | bool | True | Run between-type (paracrine) tests |
| `X_csc` | sparse | None | Pre-computed CSC format of `adata.X` |
| `lr_annotations` | dict | None | LR pair annotations (category, etc.) |
| `verbose` | bool | True | Print progress |

**Returns:** `pd.DataFrame` with columns: `cell_type` / `sender` + `receiver`, `ligand`, `receptor`, `lr_pair`, `z_score`, `p_value`, `p_adj`, `fold_enrichment`, `significant`, `test_type` ('within' or 'between').

---

### `run_lineage_analysis`

```python
cst.run_lineage_analysis(
    adata, cell_type_col, grouping, lr_pairs, indices, distances,
    expr_frac=None, tau=5.0, min_expr_frac=0.05,
    min_cells=50, fdr_method='fdr_bh', fdr_threshold=0.05,
    test_within=True, test_between=True, X_csc=None,
    lr_annotations=None, verbose=True,
)
```

Run CONSTELLATION at the lineage level for specificity assessment. Remaps fine cell types to lineage groups via a grouping dictionary, then delegates to `run_celltype_analysis()`. This enables testing whether an LR interaction is specific to a particular subtype or broadly detected across related subtypes within a lineage.

**Parameters:** Same as `run_celltype_analysis`, plus:

| Parameter | Type | Description |
|-----------|------|-------------|
| `grouping` | dict | Maps lineage name to list of fine cell types. E.g., `{"CD4_T": ["Naive_CD4_T", "Treg", "Tfh"], ...}` |

**Returns:** Same DataFrame as `run_celltype_analysis`, with `sender`/`receiver` columns containing lineage names.

---

### `validate_inputs`

```python
cst.validate_inputs(
    adata, cell_type_col, lr_pairs, indices, distances,
    tau=5.0, verbose=True,
)
```

Validate inputs for CONSTELLATION analysis. Checks data dimensions, cell type labels, gene overlap with LR pairs, spatial graph consistency, and tau relative to NN distances.

**Returns:** dict with validation results and summary statistics.

---

### `scan_celltype_pairs`

```python
cst.scan_celltype_pairs(
    adata, cell_type_col, lr_pairs, indices, distances,
    min_expr_frac=0.05, min_cells=50, verbose=True,
)
```

Scan all cell-type pairs and report the number of testable LR pairs, edge counts, and testing burden. Useful for previewing the analysis scope before running.

**Returns:** `pd.DataFrame` with testable pair counts per cell-type context.

---

## Compartment-level analysis

### `run_compartment_analysis`

```python
cst.run_compartment_analysis(
    expr_dict, compartments, coords, lr_pairs,
    k=10, tau=5.0, indices=None, distances=None,
    min_cells=10, min_expr_frac=0.01,
    fdr_method='fdr_bh', fdr_threshold=0.05, verbose=True,
)
```

Run CONSTELLATION at the compartment level. Uses the same distance-weighted kernel test (`T = dot(L, w)`) as the cell-type level, with edges restricted to within-compartment pairs.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expr_dict` | dict | — | Maps gene name to expression vector (raw counts). E.g., `{"CXCL13": np.array([...])}` |
| `compartments` | array | — | Compartment labels per cell |
| `coords` | ndarray | — | (n_cells, 2) spatial coordinates |
| `lr_pairs` | list[tuple] | — | LR pairs to test |
| `k` | int | 10 | Number of nearest neighbors |
| `tau` | float | 5.0 | Distance decay parameter (um) |
| `indices` | ndarray | None | Pre-computed KNN indices (skips graph building if provided) |
| `distances` | ndarray | None | Pre-computed KNN distances |
| `min_cells` | int | 10 | Minimum cells in a compartment to test |
| `min_expr_frac` | float | 0.01 | Minimum expression fraction |
| `fdr_method` | str | 'fdr_bh' | Multiple testing correction |
| `fdr_threshold` | float | 0.05 | FDR threshold |
| `verbose` | bool | True | Print progress |

**Returns:** `pd.DataFrame` with columns: `lr_pair`, `ligand`, `receptor`, `compartment`, `n_cells`, `n_edges`, `S_obs`, `null_mean`, `null_std`, `fold_enrichment`, `z_score`, `p_value`, `p_adj`, `significant`.

---

### `scan_compartments`

```python
cst.scan_compartments(
    expr_dict, compartments, lr_pairs,
    min_expr_frac=0.01, min_cells=10, verbose=True,
)
```

Preview compartment analysis: report testable LR pairs and expression statistics per compartment.

**Returns:** `pd.DataFrame` with expression fractions and testability flags.

---

### `compute_distance_profile`

```python
cst.compute_distance_profile(
    expr_dict, ligand, receptor, distance_values, indices,
    bin_width=15, min_cells=100,
)
```

Compute LR colocalization metrics as a function of distance from a boundary.

**Returns:** dict with keys: `bin_centers`, `coexpr_frac`, `fold_enrichment`, `n_ligand_pos`.

---

## Batch testing

### `test_within_type_lr`

```python
cst.test_within_type_lr(
    cell_type, cell_types, X, gene_names, indices, distances,
    lr_pairs_subset, lr_annotations=None, tau=5.0, X_csc=None,
)
```

Test all LR pairs within a single cell type (autocrine signaling). Both ligand- and receptor-expressing cells belong to the same type.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `cell_type` | str | Cell type to test |
| `cell_types` | array | Cell type labels for all cells |
| `X` | ndarray/sparse | Expression matrix (cells x genes, raw counts) |
| `gene_names` | list[str] | Gene names matching columns of X |
| `indices` | ndarray | KNN indices |
| `distances` | ndarray | KNN distances |
| `lr_pairs_subset` | list[tuple] | LR pairs to test |
| `tau` | float | Distance decay parameter |
| `X_csc` | sparse | Pre-computed CSC format |

**Returns:** list of dicts, each with: `cell_type`, `ligand`, `receptor`, `z_score`, `p_value`, `fold_enrichment`, `n_cells`, `n_edges`.

---

### `test_between_type_lr`

```python
cst.test_between_type_lr(
    sender_type, receiver_type, cell_types, X, gene_names,
    indices, distances, lr_pairs_subset,
    ct_indices=None, ct_sets=None,
    lr_annotations=None, tau=5.0, X_csc=None,
)
```

Test all LR pairs between a sender and receiver cell type (paracrine signaling). Ligand expression in sender cells, receptor in receiver cells.

**Returns:** list of dicts, each with: `sender`, `receiver`, `ligand`, `receptor`, `z_score`, `p_value`, `fold_enrichment`, `n_sender`, `n_receiver`, `n_edges`.

---

### `compute_expression_fractions`

```python
cst.compute_expression_fractions(adata, cell_type_col, min_expr_frac=0.01)
```

Compute the fraction of cells expressing each gene, per cell type. Used to filter testable LR pairs.

**Returns:** dict of dicts: `{cell_type: {gene: fraction}}`.

---

### `filter_testable_lr_pairs`

```python
cst.filter_testable_lr_pairs(lr_pairs, expr_frac, cell_type, min_expr_frac=0.01)
```

Filter LR pairs to those where both genes are expressed above threshold in the given cell type.

**Returns:** list of `(ligand, receptor)` tuples.

---

### `filter_testable_lr_pairs_between`

```python
cst.filter_testable_lr_pairs_between(lr_pairs, expr_frac, sender, receiver, min_expr_frac=0.01)
```

Filter LR pairs for between-type analysis: ligand expressed in sender, receptor in receiver.

**Returns:** list of `(ligand, receptor)` tuples.

---

## Targeted testing

These functions test specific LR pairs with proper full-burden FDR correction.

### `test_lr_pair_within_types`

```python
cst.test_lr_pair_within_types(
    ligand, receptor, cell_types_array, X, gene_names,
    indices, distances, expr_frac=None, burden_within=None,
    burden_between=None, tau=5.0, min_expr_frac=0.05,
    min_cells=50, fdr_method='fdr_bh', X_csc=None,
)
```

Test a specific LR pair across all qualifying cell types (within-type). FDR correction accounts for the full testing burden.

**Returns:** `pd.DataFrame` with results per cell type.

---

### `test_lr_pair_between_types`

```python
cst.test_lr_pair_between_types(
    ligand, receptor, cell_types_array, X, gene_names,
    indices, distances, expr_frac=None, burden_within=None,
    burden_between=None, tau=5.0, min_expr_frac=0.05,
    min_cells=10, fdr_method='fdr_bh', X_csc=None,
)
```

Test a specific LR pair across all qualifying sender-receiver pairs. FDR correction accounts for the full testing burden.

**Returns:** `pd.DataFrame` with results per sender-receiver pair.

---

### `test_ligand_all_receptors`

```python
cst.test_ligand_all_receptors(
    ligand, lr_pairs, cell_types_array, X, gene_names,
    indices, distances, expr_frac=None, burden_within=None,
    burden_between=None, tau=5.0, min_expr_frac=0.05,
    min_cells_within=50, min_cells_between=10,
    fdr_method='fdr_bh', X_csc=None,
)
```

Test all receptors paired with a given ligand, across all cell type contexts.

**Returns:** `pd.DataFrame`.

---

### `test_receptor_all_ligands`

```python
cst.test_receptor_all_ligands(
    receptor, lr_pairs, cell_types_array, X, gene_names,
    indices, distances, expr_frac=None, burden_within=None,
    burden_between=None, tau=5.0, min_expr_frac=0.05,
    min_cells_within=50, min_cells_between=10,
    fdr_method='fdr_bh', X_csc=None,
)
```

Test all ligands paired with a given receptor, across all cell type contexts.

**Returns:** `pd.DataFrame`.

---

### `test_custom_lr_set`

```python
cst.test_custom_lr_set(
    lr_subset, cell_types_array, X, gene_names,
    indices, distances, expr_frac=None, burden_within=None,
    burden_between=None, tau=5.0, min_expr_frac=0.05,
    min_cells_within=50, min_cells_between=10,
    fdr_method='fdr_bh', X_csc=None,
)
```

Test a custom set of LR pairs with full-burden FDR correction.

**Returns:** `pd.DataFrame`.

---

## Data loading and I/O

### `load_spatial_data`

```python
cst.load_spatial_data(h5ad_path, cell_type_col='majority_voting', refined_mapping_path=None)
```

Load spatial transcriptomics data from an h5ad file.

**Returns:** tuple of `(adata, cell_types, gene_names)`.

---

### `load_lr_pairs`

```python
cst.load_lr_pairs(lr_path)
```

Load LR pairs with annotations from a CSV file.

**Returns:** tuple of `(lr_df, lr_pairs, lr_annotations)`.

---

### `save_results`

```python
cst.save_results(results_df, output_path, sig_output_path=None, fdr_threshold=0.05)
```

Save results DataFrame to CSV. Optionally save significant results to a separate file.

---

### `print_results_summary`

```python
cst.print_results_summary(results_df, fdr_threshold=0.05)
```

Print a summary of results: total tests, significant tests, by test type.

---

### `print_top_results`

```python
cst.print_top_results(results_df, n=20, ascending=False)
```

Print the top `n` results ranked by z-score.

---

## Spatial graph

### `build_spatial_graph`

```python
cst.build_spatial_graph(adata, k=10)
```

Build a k-nearest neighbor spatial graph from AnnData spatial coordinates (`adata.obsm['spatial']`).

**Returns:** tuple of `(indices, distances)` — both `(n_cells, k)` arrays.

---

### `build_spatial_graph_from_coords`

```python
cst.build_spatial_graph_from_coords(coords, k=10, verbose=True)
```

Build a k-nearest neighbor spatial graph from raw coordinates.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `coords` | ndarray | (n_cells, 2) array of x, y coordinates |
| `k` | int | Number of nearest neighbors |

**Returns:** tuple of `(indices, distances)`.

---

## LR databases

### `load_lr_resource`

```python
cst.load_lr_resource(resource='consensus', adata=None, verbose=True)
```

Load ligand-receptor pairs from LIANA's curated databases.

**Available resources:**

| Resource | Pairs | Description |
|----------|-------|-------------|
| `consensus` | ~4,700 | LIANA consensus (union of all databases) |
| `cellphonedb` | ~1,200 | CellPhoneDB v4 |
| `celltalkdb` | ~3,300 | CellTalkDB (literature-curated) |
| `connectomedb` | ~2,200 | Connectome DB 2020 |

If `adata` is provided, automatically filters to genes in the expression matrix.

**Returns:** tuple of `(lr_df, lr_pairs)`. `lr_pairs` is a list of `(ligand, receptor)` tuples.

---

### `show_lr_resources`

```python
cst.show_lr_resources()
```

Print available LR databases with pair counts.

---

### `filter_lr_pairs_by_genes`

```python
cst.filter_lr_pairs_by_genes(lr_pairs, gene_set, verbose=True)
```

Filter LR pairs to those where both ligand and receptor are in `gene_set`.

**Returns:** list of `(ligand, receptor)` tuples.

---

## FDR and burden

### `apply_fdr_correction`

```python
cst.apply_fdr_correction(results_df, method='fdr_bh')
```

Apply Benjamini-Hochberg (or other) FDR correction to p-values.

**Returns:** `pd.DataFrame` with `p_adj` column added.

---

### `compute_testable_burden_within`

```python
cst.compute_testable_burden_within(
    adata, cell_type_col, lr_pairs, indices, distances,
    min_cells=50, min_edges=50, min_expr_frac=0.01,
)
```

Pre-compute the total number of testable within-type tests. Used for conservative FDR correction that accounts for the full testing space.

**Returns:** dict with burden counts and details.

---

### `compute_testable_burden_between`

```python
cst.compute_testable_burden_between(
    adata, cell_type_col, lr_pairs, indices, distances,
    min_cells=10, min_edges=100, min_expr_frac=0.01,
)
```

Pre-compute the total number of testable between-type tests.

**Returns:** dict with burden counts and details.

---

### `save_testable_burden` / `load_testable_burden`

```python
cst.save_testable_burden(burden_within, burden_between, output_path)
cst.load_testable_burden(path)
```

Save/load pre-computed burden information to/from disk.

---

## Cell size statistics

### `compute_cell_sizes`

```python
cst.compute_cell_sizes(adata, cell_type_col, indices, distances)
```

Compute cell size statistics (nearest-neighbor distances) per cell type.

**Returns:** `pd.DataFrame` with median, mean, std of NN distances per type.

---

### `compute_cell_sizes_with_area`

```python
cst.compute_cell_sizes_with_area(adata, cell_type_col, area_col=None)
```

Compute cell size statistics including cell area if available in `adata.obs`.

**Returns:** `pd.DataFrame`.

---

### `summarize_cell_sizes`

```python
cst.summarize_cell_sizes(adata, cell_type_col, indices, distances, area_col=None)
```

Comprehensive cell size summary combining NN distances and area.

**Returns:** `pd.DataFrame`.

---

### `print_cell_size_table`

```python
cst.print_cell_size_table(df, top_n=None)
```

Print cell size statistics in a readable table format.

---

## Visualization — cell-type interactions

### `plot_celltype_pair_heatmap`

```python
cst.plot_celltype_pair_heatmap(
    results_df, value_col='z_score', agg_func='mean',
    sig_only=True, fdr_threshold=0.05, cmap='RdBu_r',
    center=0, figsize=(12, 10), title=None, save_path=None,
)
```

Heatmap of aggregated LR interaction scores between cell type pairs.

---

### `plot_combined_heatmap`

```python
cst.plot_combined_heatmap(
    results_df, value_col='z_score', agg_func='count',
    sig_only=True, fdr_threshold=0.05, cmap='YlOrBr',
    cell_type_order=None, cluster=False,
    figsize=None, title=None, save_path=None,
)
```

Combined heatmap with within-type on diagonal and between-type off-diagonal. Shows sender (rows) x receiver (columns).

---

### `plot_interaction_dotplot`

```python
cst.plot_interaction_dotplot(
    results_df, sig_only=True, fdr_threshold=0.05,
    color_col='z_score', color_agg='mean',
    min_count=1, max_dot_size=400, cmap='YlOrBr',
    cell_type_order=None, cluster=False,
    figsize=None, title=None, save_path=None,
)
```

Dot plot of cell-type interactions. Dot size = count of significant LR pairs, color = aggregated z-score.

---

### `plot_interaction_network`

```python
cst.plot_interaction_network(
    results_df, sig_only=True, fdr_threshold=0.05,
    edge_metric='count', node_size_metric='total_interactions',
    node_colors=None, show_self_loops=True,
    figsize=(10, 10), title=None, save_path=None,
)
```

Circular network diagram of cell-type interactions. Nodes = cell types, edges = significant LR pairs.

---

### `plot_cell_lineage_tree`

```python
cst.plot_cell_lineage_tree(
    lineage_dict, counts_dict,
    figsize=(10, 8), save_path=None,
)
```

Plot a hierarchical lineage tree of cell types with abundance counts.

---

### `plot_celltype_barplot`

```python
cst.plot_celltype_barplot(
    counts_dict, order=None, colors=None,
    figsize=(10, 5), title=None, save_path=None,
)
```

Bar chart of cell type abundances.

---

## Visualization — LR-level plots

### `plot_lr_dotplot`

```python
cst.plot_lr_dotplot(
    results_df, lr_pairs=None, celltype_pairs=None,
    top_n_lr=20, top_n_pairs=20,
    size_col='p_adj', color_col='z_score',
    sig_only=True, fdr_threshold=0.05,
    cmap='RdBu_r', center=0,
    figsize=None, title=None, save_path=None,
)
```

CellPhoneDB-style dot plot: LR pairs (rows) x cell-type pairs (columns). Size = significance, color = z-score.

---

### `plot_spatial_interactions`

```python
cst.plot_spatial_interactions(
    adata, cell_type_col, results_df=None,
    ligand=None, receptor=None, sender=None, receiver=None,
    indices=None, distances=None, tau=5.0,
    mode='cell_type', sender_color='#E8B87D', receiver_color='#7BA3A8',
    figsize=(12, 10), title=None, save_path=None,
)
```

Spatial scatter plot showing sender and receiver cells with interaction edges overlaid on tissue coordinates.

---

### `plot_lr_interaction_summary`

```python
cst.plot_lr_interaction_summary(
    results_df, group_col='cell_type', sig_only=True,
    fdr_threshold=0.05, top_n=20,
    figsize=(10, 6), save_path=None,
)
```

Bar chart summarizing the number of significant interactions per cell type.

---

### `plot_lr_category_heatmap`

```python
cst.plot_lr_category_heatmap(
    results_df, category_col='ligand_category_major',
    sig_only=True, fdr_threshold=0.05,
    figsize=(14, 8), save_path=None,
)
```

Heatmap of LR interactions grouped by functional category.

---

## Visualization — compartment plots

### `plot_compartment_heatmap`

```python
cst.plot_compartment_heatmap(
    results_df, value_col='significant', compartment_col='compartment',
    lr_col='lr_pair', compartment_order=None, cmap='RdYlGn',
    figsize=None, title='CONSTELLATION Compartment Detection',
    show_values=True, save_path=None,
)
```

Heatmap of LR detection across spatial compartments.

---

### `plot_compartment_spatial`

```python
cst.plot_compartment_spatial(
    coords, compartments, compartment_colors=None,
    point_size=0.5, figsize=(12, 10),
    title='Spatial Compartments', xlim=None, ylim=None, save_path=None,
)
```

Spatial scatter plot colored by compartment labels.

**Returns:** tuple of `(fig, ax)`.

---

### `plot_distance_profile`

```python
cst.plot_distance_profile(
    profile_data, ligand=None, receptor=None,
    boundary_pos=0, interface_pos=None,
    color='#5a7aaa', figsize=(8, 4),
    xlabel='Distance from Boundary (um)',
    ylabel='Coexpression Fold Enrichment',
    show_reference=True, save_path=None,
)
```

Plot LR colocalization as a function of distance from a boundary.

---

### `plot_boundary_profile`

```python
cst.plot_boundary_profile(
    profile_data, ligand=None, receptor=None,
    boundary_pos=0, interface_width=50,
    xlim=None, zone_labels=True,
    figsize=(8, 8), color='#5a7aaa', save_path=None,
)
```

Three-row boundary profile: coexpression fraction, cell counts, and fold enrichment vs distance.

**Returns:** tuple of `(fig, axes)`.

---

## Cell Ontology

### `search_cell_ontology`

```python
cst.search_cell_ontology(query, exact=False, max_results=10)
```

Search the Cell Ontology (CL) for terms matching a query string via the EBI OLS API.

**Returns:** list of dicts with `id`, `label`, `description`.

---

### `get_term_info`

```python
cst.get_term_info(cl_id)
```

Get detailed information about a Cell Ontology term.

**Returns:** dict with term details.

---

### `get_ancestors` / `get_children`

```python
cst.get_ancestors(cl_id, include_self=False)
cst.get_children(cl_id)
```

Traverse the ontology hierarchy upward (ancestors) or downward (children).

**Returns:** list of dicts with `id`, `label`.

---

### `map_annotation_to_ontology`

```python
cst.map_annotation_to_ontology(annotation, custom_mappings=None)
```

Map a cell type annotation string to a Cell Ontology ID.

**Returns:** dict with `id`, `label`, `score`.

---

### `build_annotation_mapping`

```python
cst.build_annotation_mapping(annotations, custom_mappings=None, verbose=True)
```

Map a list of cell type annotations to Cell Ontology IDs.

**Returns:** dict mapping annotation strings to ontology terms.

---

### `get_lymph_node_refined_groups`

```python
cst.get_lymph_node_refined_groups()
```

Get predefined refined groupings for lymph node cell types (17 lineage groups from 29 fine types). Compatible with `run_lineage_analysis()`.

**Returns:** dict mapping group name to list of cell types.

---

### `get_lymph_node_coarse_groups`

```python
cst.get_lymph_node_coarse_groups()
```

Get coarser groupings (major lineages: B, T, Myeloid, Stromal) for lymph node cells.

**Returns:** dict mapping lineage name to list of cell types.

---

### `build_hierarchy_for_plotting`

```python
cst.build_hierarchy_for_plotting(refined_groups, coarse_groups, colors=None)
```

Build a nested hierarchy dictionary suitable for `plot_cell_lineage_tree()`.

**Returns:** dict with hierarchical structure.

---

### `print_hierarchy`

```python
cst.print_hierarchy(hierarchy, indent=0)
```

Pretty-print a hierarchy dictionary.
