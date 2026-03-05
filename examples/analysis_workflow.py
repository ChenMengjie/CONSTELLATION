"""
CONSTELLATION — Cell-Type-Only Testing Workflow
================================================

Pseudocode showing the cell-type-only analysis pipeline for a new
spatial transcriptomics dataset. This mode tests LR interactions
conditioned on cell type identity alone (no compartment information).

A separate compartment-aware testing mode will incorporate spatial
compartment/zone information for context-specific interaction inference.

Usage:
    import constellation as cst
"""

import pandas as pd
import constellation as cst

# =============================================================================
# Step 1: Load data
# =============================================================================
# Load spatial transcriptomics data (h5ad with spatial coordinates)
adata, cell_types, gene_names = cst.load_spatial_data(
    "new_dataset.h5ad",
    cell_type_col="cell_type"       # column in adata.obs with cell type labels
)

# Load ligand-receptor database
lr_df, lr_pairs, lr_ann = cst.load_lr_pairs("lr_database.csv")
# lr_pairs: list of (ligand, receptor) tuples
# lr_ann:   dict {(lig, rec): {category info}}

X = adata.X                          # expression matrix (cells × genes)

# =============================================================================
# Step 2: Build spatial graph (k-NN)
# =============================================================================
indices, distances = cst.build_spatial_graph(adata, k=10)
# indices:   (n_cells, k) — neighbor cell indices
# distances: (n_cells, k) — neighbor distances in microns

# =============================================================================
# Step 3: Compute expression fractions per cell type
# =============================================================================
# Used to filter LR pairs to those actually expressed
expr_frac = cst.compute_expression_fractions(adata, cell_type_col="cell_type")
# expr_frac: {cell_type: {gene: fraction_expressing}}

# =============================================================================
# Step 4: Compute testable burden (for conservative FDR correction)
# =============================================================================
# Counts total number of tests that *would* be performed in a full analysis,
# so that targeted queries can use the same FDR correction.

burden_within = cst.compute_testable_burden_within(
    adata, "cell_type", lr_pairs, indices, distances,
    expr_frac=expr_frac
)

burden_between = cst.compute_testable_burden_between(
    adata, "cell_type", lr_pairs, indices, distances,
    expr_frac=expr_frac
)

# Save/load burden for later targeted queries
cst.save_testable_burden(burden_within, burden_between, "burden.json")
# burden_within, burden_between = cst.load_testable_burden("burden.json")

# =============================================================================
# Step 5a: Batch within-type testing (autocrine signaling)
# =============================================================================
# For each cell type, test all expressed LR pairs for spatial co-expression.
# Permutation reuse: permutations generated once per cell type, reused across
# all LR pairs → major speedup.

unique_types = sorted(set(cell_types))
within_results = []

for cell_type in unique_types:
    # Filter to LR pairs where both ligand & receptor are expressed
    testable = cst.filter_testable_lr_pairs(lr_pairs, expr_frac, cell_type)

    res = cst.test_within_type_lr(
        cell_type=cell_type,
        cell_types=cell_types,
        X=X,
        gene_names=gene_names,
        indices=indices,
        distances=distances,
        lr_pairs_subset=testable,
        lr_annotations=lr_ann,
        n_perm=200,
        tau=5.0                     # kernel decay distance (microns)
    )
    within_results.extend(res)

within_df = pd.DataFrame(within_results)

# =============================================================================
# Step 5b: Batch between-type testing (paracrine signaling)
# =============================================================================
# For each sender→receiver pair, test LR pairs where the ligand is expressed
# in the sender and receptor in the receiver.

between_results = []

for sender in unique_types:
    for receiver in unique_types:
        # Filter: ligand expressed in sender, receptor in receiver
        testable = cst.filter_testable_lr_pairs_between(
            lr_pairs, expr_frac, sender, receiver
        )

        res = cst.test_between_type_lr(
            sender_type=sender,
            receiver_type=receiver,
            cell_types=cell_types,
            X=X,
            gene_names=gene_names,
            indices=indices,
            distances=distances,
            lr_pairs_subset=testable,
            lr_annotations=lr_ann,
            n_perm=200,
            tau=5.0
        )
        between_results.extend(res)

between_df = pd.DataFrame(between_results)

# =============================================================================
# Step 6: FDR correction & save
# =============================================================================
within_df = cst.apply_fdr_correction(within_df, method='fdr_bh')
between_df = cst.apply_fdr_correction(between_df, method='fdr_bh')

cst.save_results(within_df, "within_type_results.csv")
cst.save_results(between_df, "between_type_results.csv")

# =============================================================================
# Step 7: Targeted queries (with full testable burden FDR)
# =============================================================================
# These functions apply FDR correction using the FULL testable burden
# (from Step 4), so results are directly comparable to the batch analysis.

# 7a. Test a specific LR pair across all cell type pairs
res_pair = cst.test_lr_pair_between_types(
    "CXCL13", "CXCR5",
    cell_types, X, gene_names, indices, distances,
    burden_info=burden_between
)

# 7b. Test all receptors for a given ligand
res_ligand = cst.test_ligand_all_receptors(
    "CXCL13", lr_pairs,
    cell_types, X, gene_names, indices, distances,
    analysis_type="between",
    burden_info=burden_between
)

# 7c. Test all ligands for a given receptor
res_receptor = cst.test_receptor_all_ligands(
    "CXCR5", lr_pairs,
    cell_types, X, gene_names, indices, distances,
    analysis_type="between",
    burden_info=burden_between
)

# 7d. Test a custom set of LR pairs
res_custom = cst.test_custom_lr_set(
    [("CCL19", "CCR7"), ("CXCL12", "CXCR4"), ("CCL21", "CCR7")],
    cell_types, X, gene_names, indices, distances,
    analysis_type="between",
    burden_info=burden_between
)

# =============================================================================
# Step 7e: Lineage-based testing (bottom-up, up to 2 generations)
# =============================================================================
# Generalized lineage testing using Cell Ontology (CL).
#
# Strategy:
#   1. Map each annotated cell type to a Cell Ontology ID (CL:XXXXXXX)
#   2. Query CL for parent/grandparent terms → build hierarchy automatically
#   3. Merge sibling cell types at each generation and re-test
#
#   Generation 0 (leaf):        Annotated cell types as-is (done in Steps 5a/5b)
#   Generation 1 (parent):      Merge siblings sharing a CL parent
#   Generation 2 (grandparent): Merge one more level up in CL
#
# This is fully general — works for any tissue as long as cell type
# annotations can be mapped to Cell Ontology terms.

# --- Step 7e.1: Map annotations to Cell Ontology IDs ---
unique_types = sorted(set(cell_types))

# Automatic mapping via OLS API (with optional user overrides)
cl_mapping = cst.build_annotation_mapping(
    unique_types,
    custom_mappings={                        # override ambiguous names
        "FDC": None,                         # stromal, not in CL immune branch
        "FRC": None,
        "Endothelial": None,
    },
    verbose=True
)
# cl_mapping: {"Naive CD4 T": "CL:0000895", "Treg": "CL:0000815", ...}

# --- Step 7e.2: Build hierarchy from CL parent/grandparent terms ---
# For each mapped cell type, query ancestors up to 2 generations
# to discover shared parents and grandparents.

leaf_to_parent = {}     # {leaf_label: parent_label}
leaf_to_grandparent = {} # {leaf_label: grandparent_label}

for cell_type_label, cl_id in cl_mapping.items():
    if cl_id is None:
        continue  # skip unmapped types (e.g. stromal)

    ancestors = cst.get_ancestors(cl_id)
    # ancestors: [{"id": "CL:...", "label": "parent"}, {"id": "CL:...", "label": "grandparent"}, ...]

    if len(ancestors) >= 1:
        leaf_to_parent[cell_type_label] = ancestors[0]["label"]
    if len(ancestors) >= 2:
        leaf_to_grandparent[cell_type_label] = ancestors[1]["label"]

# Build generation groups: {group_label: [leaf_labels]}
from collections import defaultdict

gen1_groups = defaultdict(list)  # parent level
for leaf, parent in leaf_to_parent.items():
    gen1_groups[parent].append(leaf)
gen1_groups = dict(gen1_groups)

gen2_groups = defaultdict(list)  # grandparent level
for leaf, grandparent in leaf_to_grandparent.items():
    gen2_groups[grandparent].append(leaf)
gen2_groups = dict(gen2_groups)

# Example output for lymph node:
#   gen1_groups = {
#       "CD4-positive, alpha-beta T cell": ["Naive CD4 T", "Memory CD4 T", "Treg", "Tfh"],
#       "CD8-positive, alpha-beta T cell": ["Naive CD8 T", "Effector CD8 T"],
#       "B cell":     ["Naive B", "GC B", "Memory B", "Plasma"],
#       "dendritic cell": ["cDC", "pDC"],
#       ...
#   }
#   gen2_groups = {
#       "T cell":     ["Naive CD4 T", ..., "Naive CD8 T", ...],
#       "lymphocyte": ["Naive B", ..., "NK"],
#       "myeloid cell": ["cDC", "pDC", "Macrophage"],
#       ...
#   }

# --- Step 7e.3: Relabel cells at each generation & re-test ---
import numpy as np

for gen_label, gen_groups in [("gen1", gen1_groups), ("gen2", gen2_groups)]:

    # Build leaf → merged group mapping
    leaf_to_group = {}
    for group_name, leaves in gen_groups.items():
        for leaf in leaves:
            leaf_to_group[leaf] = group_name

    # Relabel cell types (cells not in hierarchy keep original label)
    cell_types_merged = np.array([
        leaf_to_group.get(ct, ct) for ct in cell_types
    ])

    unique_merged = sorted(set(cell_types_merged))

    # Recompute expression fractions at merged level
    # (temporarily assign merged labels to adata.obs)
    adata.obs[f'cell_type_{gen_label}'] = cell_types_merged
    expr_frac_merged = cst.compute_expression_fractions(
        adata, cell_type_col=f'cell_type_{gen_label}'
    )

    # Between-type testing at this generation level
    merged_results = []
    for sender in unique_merged:
        for receiver in unique_merged:
            testable = cst.filter_testable_lr_pairs_between(
                lr_pairs, expr_frac_merged, sender, receiver
            )
            res = cst.test_between_type_lr(
                sender_type=sender,
                receiver_type=receiver,
                cell_types=cell_types_merged,
                X=X,
                gene_names=gene_names,
                indices=indices,
                distances=distances,
                lr_pairs_subset=testable,
                lr_annotations=lr_ann,
                n_perm=200,
                tau=5.0
            )
            merged_results.extend(res)

    merged_df = pd.DataFrame(merged_results)
    merged_df = cst.apply_fdr_correction(merged_df)
    merged_df['generation'] = gen_label
    cst.save_results(merged_df, f"between_type_{gen_label}_results.csv")

# =============================================================================
# Step 8: Summarize & visualize
# =============================================================================
cst.print_results_summary(between_df)
cst.print_top_results(between_df, n=20)

# Heatmap of z-scores across sender→receiver pairs
cst.plot_celltype_pair_heatmap(
    between_df, value_col='z_score', sig_only=True,
    save_path="celltype_heatmap.pdf"
)

# Category-stratified heatmaps
cst.plot_lr_category_heatmap(
    between_df, category_col='ligand_category_major',
    save_path="category_heatmaps.pdf"
)

# Summary bar chart of significant interactions per cell type
cst.plot_lr_interaction_summary(
    between_df, group_col='sender',
    save_path="interaction_summary.pdf"
)
