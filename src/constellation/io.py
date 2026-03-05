"""
Spatial LR I/O Functions

Functions for loading data and saving results.
"""

import pandas as pd
import scanpy as sc
from sklearn.neighbors import NearestNeighbors
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm


def load_spatial_data(h5ad_path, cell_type_col='majority_voting',
                      refined_mapping_path=None):
    """
    Load spatial transcriptomics data.

    Parameters
    ----------
    h5ad_path : str
        Path to h5ad file
    cell_type_col : str
        Column name for cell type annotation
    refined_mapping_path : str, optional
        Path to refined cell type mapping CSV

    Returns
    -------
    adata : AnnData
        Annotated data object
    cell_types : array
        Cell type labels
    gene_names : list
        Gene names
    """
    print(f"Loading data from {h5ad_path}...")
    adata = sc.read_h5ad(h5ad_path)

    # Apply refined mapping if provided
    if refined_mapping_path is not None:
        mapping = pd.read_csv(refined_mapping_path)
        ct_to_refined = dict(zip(mapping['majority_voting'], mapping['refined_group']))
        adata.obs['refined_group'] = adata.obs['majority_voting'].map(ct_to_refined)
        cell_type_col = 'refined_group'

    cell_types = adata.obs[cell_type_col].values
    gene_names = list(adata.var_names)

    print(f"Loaded: {adata.n_obs:,} cells, {adata.n_vars:,} genes")
    print(f"Cell types: {len(set(cell_types))}")

    return adata, cell_types, gene_names


def load_lr_pairs(lr_path):
    """
    Load LR pairs with annotations.

    Parameters
    ----------
    lr_path : str
        Path to LR pairs CSV

    Returns
    -------
    lr_df : DataFrame
        LR pairs dataframe
    lr_pairs : list of tuples
        List of (ligand, receptor) pairs
    lr_annotations : dict
        Annotations {(lig, rec): {...}}
    """
    lr_df = pd.read_csv(lr_path)
    lr_pairs = list(zip(lr_df['ligand'], lr_df['receptor']))

    # Build annotations dict for fast lookup
    lr_annotations = {}
    for _, row in lr_df.iterrows():
        lr_annotations[(row['ligand'], row['receptor'])] = {
            'ligand_category': row.get('ligand_category'),
            'ligand_category_major': row.get('ligand_category_major', 'Other'),
            'receptor_category': row.get('receptor_category')
        }

    print(f"Loaded {len(lr_pairs)} LR pairs")
    return lr_df, lr_pairs, lr_annotations


def load_lr_resource(resource='consensus', adata=None, verbose=True):
    """
    Load ligand-receptor pairs from LIANA's curated databases.

    Requires the ``liana`` package (``pip install liana``).

    Parameters
    ----------
    resource : str
        Name of the LR resource. Available options:
        'consensus' (default), 'cellphonedb', 'cellchatdb', 'celltalkdb',
        'connectomedb2020', 'icellnet', 'baccin2019', 'cellcall',
        'cellinker', 'embrace', 'guide2pharma', 'hpmr', 'italk',
        'kirouac2010', 'lrdb', 'mouseconsensus', 'ramilowski2015'.
    adata : AnnData, optional
        If provided, filter to LR pairs where both genes are in adata.var_names.
    verbose : bool
        Print summary information (default: True).

    Returns
    -------
    lr_df : DataFrame
        LR pairs dataframe with 'ligand' and 'receptor' columns.
    lr_pairs : list of tuples
        List of (ligand, receptor) pairs.

    Examples
    --------
    >>> import constellation as cst
    >>> # List available resources
    >>> cst.show_lr_resources()
    >>> # Load default (consensus) resource
    >>> lr_df, lr_pairs = cst.load_lr_resource()
    >>> # Load CellPhoneDB, filtered to genes in your data
    >>> lr_df, lr_pairs = cst.load_lr_resource('cellphonedb', adata=adata)
    """
    try:
        import liana as li
    except ImportError:
        raise ImportError(
            "liana is required to load LR resources. "
            "Install with: pip install liana")

    available = li.rs.show_resources()
    if resource not in available:
        raise ValueError(
            f"Unknown resource '{resource}'. "
            f"Available: {available}")

    lr_df = li.rs.select_resource(resource).copy()
    lr_df = lr_df[['ligand', 'receptor']].drop_duplicates().reset_index(drop=True)
    n_total = len(lr_df)

    if adata is not None:
        gene_set = set(adata.var_names)
        mask = lr_df['ligand'].isin(gene_set) & lr_df['receptor'].isin(gene_set)
        lr_df = lr_df[mask].reset_index(drop=True)

    lr_pairs = list(zip(lr_df['ligand'], lr_df['receptor']))

    if verbose:
        print(f"Resource: {resource} ({n_total} pairs total)")
        if adata is not None:
            print(f"Filtered to genes in adata: {len(lr_pairs)} testable pairs")
        else:
            print(f"Loaded {len(lr_pairs)} LR pairs")

    return lr_df, lr_pairs


def show_lr_resources():
    """
    Print available LR resources from LIANA.

    Requires the ``liana`` package (``pip install liana``).

    Returns
    -------
    list of str
        Available resource names.
    """
    try:
        import liana as li
    except ImportError:
        raise ImportError(
            "liana is required. Install with: pip install liana")

    resources = li.rs.show_resources()
    print("Available LR resources:")
    for r in sorted(resources):
        df = li.rs.select_resource(r)
        print(f"  {r:<25s} {len(df):>5,} pairs")
    return resources


def build_spatial_graph(adata, k=10):
    """
    Build k-nearest neighbor spatial graph.

    Parameters
    ----------
    adata : AnnData
        Annotated data object with spatial coordinates in obsm['spatial']
    k : int
        Number of neighbors

    Returns
    -------
    indices : array
        KNN indices (n_cells x k)
    distances : array
        KNN distances (n_cells x k)
    """
    print(f"Building spatial graph with k={k}...")
    coords = adata.obsm['spatial']

    nn = NearestNeighbors(n_neighbors=k+1, algorithm='ball_tree')
    nn.fit(coords)
    distances, indices = nn.kneighbors(coords)

    # Remove self
    distances = distances[:, 1:]
    indices = indices[:, 1:]

    print(f"Graph: {len(coords):,} cells, k={k} neighbors")
    return indices, distances


def apply_fdr_correction(results_df, method='fdr_bh'):
    """
    Apply multiple testing correction.

    Returns a copy of the input DataFrame with a 'p_adj' column added.
    The original DataFrame is not modified.

    Parameters
    ----------
    results_df : DataFrame
        Results with 'p_value' column
    method : str
        Correction method (default: Benjamini-Hochberg)

    Returns
    -------
    DataFrame
        Copy of results_df with 'p_adj' column added
    """
    df = results_df.copy()
    _, pvals_adj, _, _ = multipletests(df['p_value'], method=method)
    df['p_adj'] = pvals_adj
    return df


def save_results(results_df, output_path, sig_output_path=None, fdr_threshold=0.05):
    """
    Save results to files.

    Parameters
    ----------
    results_df : DataFrame
        Results dataframe
    output_path : str
        Path for full results (parquet)
    sig_output_path : str, optional
        Path for significant results (CSV)
    fdr_threshold : float
        FDR threshold for significant results
    """
    # Save full results
    results_df.to_parquet(output_path, index=False)
    print(f"Results saved to: {output_path}")

    # Save significant results
    if sig_output_path is not None:
        sig_df = results_df[results_df['p_adj'] < fdr_threshold].copy()
        sig_df = sig_df.sort_values('z_score', ascending=False)
        sig_df.to_csv(sig_output_path, index=False)
        print(f"Significant results: {sig_output_path}")

    return results_df


def load_cell_type_mapping(mapping_path, source_col='majority_voting',
                           target_col='refined_group'):
    """
    Load cell type mapping from CSV.

    Parameters
    ----------
    mapping_path : str
        Path to mapping CSV
    source_col : str
        Column name for original cell type labels (default: 'majority_voting')
    target_col : str
        Column name for mapped cell type labels (default: 'refined_group')

    Returns
    -------
    dict
        {original_type: refined_type}
    """
    mapping = pd.read_csv(mapping_path)
    return dict(zip(mapping[source_col], mapping[target_col]))


def print_results_summary(results_df, fdr_threshold=0.05):
    """
    Print summary of results.

    Parameters
    ----------
    results_df : DataFrame
        Results with 'p_value' and 'p_adj' columns
    fdr_threshold : float
        FDR threshold
    """
    n_total = len(results_df)
    n_sig_raw = (results_df['p_value'] < 0.05).sum()
    n_sig_adj = (results_df['p_adj'] < fdr_threshold).sum()

    print(f"\nTotal tests: {n_total:,}")
    print(f"Significant at p<0.05: {n_sig_raw:,} ({100*n_sig_raw/n_total:.1f}%)")
    print(f"Significant at FDR<{fdr_threshold}: {n_sig_adj:,} ({100*n_sig_adj/n_total:.1f}%)")


def print_top_results(results_df, n=20, ascending=False):
    """
    Print top results by z-score.

    Parameters
    ----------
    results_df : DataFrame
        Results dataframe
    n : int
        Number of results to show
    ascending : bool
        If True, show lowest z-scores (negative interactions)
    """
    if ascending:
        top = results_df[results_df['z_score'] < 0].nsmallest(n, 'z_score')
    else:
        top = results_df[results_df['z_score'] > 0].nlargest(n, 'z_score')

    for _, r in top.iterrows():
        sig = "***" if r['p_adj'] < 0.001 else "**" if r['p_adj'] < 0.01 else "*" if r['p_adj'] < 0.05 else ""

        if 'sender' in r:  # Between-type
            print(f"  {r['ligand']:10} → {r['receptor']:10} | "
                  f"{r['sender'][:15]:15} → {r['receiver'][:15]:15} | z={r['z_score']:6.1f} {sig}")
        else:  # Within-type
            print(f"  {r['ligand']:10} → {r['receptor']:10} | "
                  f"{r['cell_type']:20} | z={r['z_score']:6.1f} {sig}")


# =============================================================================
# TESTABLE BURDEN COMPUTATION
# =============================================================================

def compute_testable_burden_within(adata, cell_type_col, lr_pairs, indices, distances,
                                    min_cells=50, min_edges=50, min_expr_frac=0.01):
    """
    Compute the total number of testable within-type tests.

    This determines the multiple testing burden for FDR correction.
    Also computes recommended tau from median nearest neighbor distance.

    Parameters
    ----------
    adata : AnnData
        Annotated data object
    cell_type_col : str
        Column name for cell type
    lr_pairs : list of tuples
        All LR pairs [(ligand, receptor), ...]
    indices : array
        KNN indices from spatial graph
    distances : array
        KNN distances from spatial graph
    min_cells : int
        Minimum cells per type
    min_edges : int
        Minimum within-type edges
    min_expr_frac : float
        Minimum expression fraction

    Returns
    -------
    dict
        {
            'total_tests': int,
            'testable_cell_types': list,
            'testable_lr_per_type': {cell_type: [(lig, rec), ...]},
            'n_tests_per_type': {cell_type: int},
            'expr_frac': {cell_type: {gene: frac}},
            'tau': float (recommended tau from median NN1 distance)
        }
    """
    import numpy as np
    from scipy.sparse import issparse

    cell_types = adata.obs[cell_type_col].values
    gene_names = list(adata.var_names)
    X = adata.X
    unique_types = list(set(cell_types))

    # Compute expression fractions
    expr_frac = {}
    for ct in unique_types:
        idx = np.where(cell_types == ct)[0]
        if issparse(X):
            X_ct = X[idx, :].toarray()
        else:
            X_ct = X[idx, :]
        frac = (X_ct > 0).sum(axis=0) / len(idx)
        expr_frac[ct] = dict(zip(gene_names, np.array(frac).flatten()))

    # Find testable cell types and LR pairs
    testable_cell_types = []
    testable_lr_per_type = {}
    n_tests_per_type = {}

    # Pre-compute total cells for vectorized edge counting
    n_total = len(cell_types)

    for ct in tqdm(unique_types, desc="Computing within-type burden"):
        ct_idx = np.where(cell_types == ct)[0]
        if len(ct_idx) < min_cells:
            continue

        # Vectorized edge counting (O(n_ct * k) instead of nested loops)
        is_ct = np.zeros(n_total, dtype=bool)
        is_ct[ct_idx] = True
        nbr_indices = indices[ct_idx, :]  # (n_ct, k)
        mask = is_ct[nbr_indices]         # boolean lookup O(1) per element
        n_edges = mask.sum()

        if n_edges < min_edges:
            continue

        # Find testable LR pairs for this type
        testable_lr = []
        ct_frac = expr_frac[ct]
        for lig, rec in lr_pairs:
            if (ct_frac.get(lig, 0) >= min_expr_frac and
                ct_frac.get(rec, 0) >= min_expr_frac):
                testable_lr.append((lig, rec))

        if testable_lr:
            testable_cell_types.append(ct)
            testable_lr_per_type[ct] = testable_lr
            n_tests_per_type[ct] = len(testable_lr)

    total_tests = sum(n_tests_per_type.values())

    # Compute recommended tau from median NN1 distance (cell diameter proxy)
    nn1_distances = distances[:, 0]
    tau = float(np.median(nn1_distances))

    print(f"Within-type testable burden:")
    print(f"  Cell types: {len(testable_cell_types)}")
    print(f"  Total tests: {total_tests:,}")
    print(f"  Recommended tau (median NN1): {tau:.2f} μm")

    return {
        'total_tests': total_tests,
        'testable_cell_types': testable_cell_types,
        'testable_lr_per_type': testable_lr_per_type,
        'n_tests_per_type': n_tests_per_type,
        'expr_frac': expr_frac,
        'tau': tau,
    }


def compute_testable_burden_between(adata, cell_type_col, lr_pairs, indices, distances,
                                     min_cells=10, min_edges=100, min_expr_frac=0.01):
    """
    Compute the total number of testable between-type tests.

    Also computes recommended tau from median nearest neighbor distance.

    Parameters
    ----------
    adata : AnnData
        Annotated data object
    cell_type_col : str
        Column name for cell type
    lr_pairs : list of tuples
        All LR pairs
    indices : array
        KNN indices
    distances : array
        KNN distances
    min_cells : int
        Minimum cells per type
    min_edges : int
        Minimum edges between types
    min_expr_frac : float
        Minimum expression fraction

    Returns
    -------
    dict
        {
            'total_tests': int,
            'testable_pairs': [(sender, receiver, n_edges, n_lr), ...],
            'testable_lr_per_pair': {(sender, receiver): [(lig, rec), ...]},
            'expr_frac': {cell_type: {gene: frac}},
            'tau': float (recommended tau from median NN1 distance)
        }
    """
    import numpy as np
    from scipy.sparse import issparse

    cell_types = adata.obs[cell_type_col].values
    gene_names = list(adata.var_names)
    X = adata.X
    unique_types = list(set(cell_types))

    # Compute expression fractions
    expr_frac = {}
    ct_indices = {}
    ct_sets = {}

    for ct in unique_types:
        idx = np.where(cell_types == ct)[0]
        ct_indices[ct] = idx
        ct_sets[ct] = set(idx)

        if issparse(X):
            X_ct = X[idx, :].toarray()
        else:
            X_ct = X[idx, :]
        frac = (X_ct > 0).sum(axis=0) / len(idx)
        expr_frac[ct] = dict(zip(gene_names, np.array(frac).flatten()))

    # Find testable cell type pairs
    testable_pairs = []
    testable_lr_per_pair = {}

    # Pre-compute for vectorized edge counting
    n_total = len(cell_types)

    # Pre-compute boolean masks for all cell types (for fast lookup)
    ct_masks = {}
    for ct in unique_types:
        mask = np.zeros(n_total, dtype=bool)
        mask[ct_indices[ct]] = True
        ct_masks[ct] = mask

    # Filter to valid senders (enough cells)
    valid_senders = [s for s in unique_types if len(ct_indices[s]) >= min_cells]

    for sender in tqdm(valid_senders, desc="Computing between-type burden"):
        sender_idx = ct_indices[sender]

        # Get neighbor indices for all sender cells once
        sender_nbr_indices = indices[sender_idx, :]  # (n_sender, k)

        for receiver in unique_types:
            if sender == receiver:
                continue

            receiver_idx = ct_indices[receiver]
            if len(receiver_idx) < min_cells:
                continue

            # Vectorized edge counting: sender → receiver
            is_receiver = ct_masks[receiver]
            mask = is_receiver[sender_nbr_indices]  # boolean lookup O(1) per element
            n_edges = mask.sum()

            if n_edges < min_edges:
                continue

            # Find testable LR pairs
            testable_lr = []
            for lig, rec in lr_pairs:
                if (expr_frac[sender].get(lig, 0) >= min_expr_frac and
                    expr_frac[receiver].get(rec, 0) >= min_expr_frac):
                    testable_lr.append((lig, rec))

            if testable_lr:
                testable_pairs.append((sender, receiver, n_edges, len(testable_lr)))
                testable_lr_per_pair[(sender, receiver)] = testable_lr

    total_tests = sum(info[3] for info in testable_pairs)

    # Compute recommended tau from median NN1 distance
    nn1_distances = distances[:, 0]
    tau = float(np.median(nn1_distances))

    print(f"Between-type testable burden:")
    print(f"  Cell type pairs: {len(testable_pairs)}")
    print(f"  Total tests: {total_tests:,}")
    print(f"  Recommended tau (median NN1): {tau:.2f} μm")

    return {
        'total_tests': total_tests,
        'testable_pairs': testable_pairs,
        'testable_lr_per_pair': testable_lr_per_pair,
        'expr_frac': expr_frac,
        'ct_indices': ct_indices,
        'ct_sets': ct_sets,
        'tau': tau,
    }


def save_testable_burden(burden_within, burden_between, output_path):
    """
    Save testable burden information to a file.

    Parameters
    ----------
    burden_within : dict
        Output from compute_testable_burden_within
    burden_between : dict
        Output from compute_testable_burden_between
    output_path : str
        Path to save (pickle format)
    """
    import pickle

    data = {
        'within': burden_within,
        'between': burden_between,
    }

    with open(output_path, 'wb') as f:
        pickle.dump(data, f)

    print(f"Saved testable burden to: {output_path}")


def load_testable_burden(path):
    """
    Load testable burden information.

    Parameters
    ----------
    path : str
        Path to burden file

    Returns
    -------
    dict
        {'within': {...}, 'between': {...}}
    """
    import pickle

    with open(path, 'rb') as f:
        data = pickle.load(f)

    print(f"Within-type total tests: {data['within']['total_tests']:,}")
    print(f"Between-type total tests: {data['between']['total_tests']:,}")

    return data


# =============================================================================
# CELL SIZE / SPATIAL STATISTICS
# =============================================================================

def compute_cell_sizes(adata, cell_type_col, indices, distances):
    """
    Compute cell size statistics for each cell type based on nearest neighbor distances.

    Uses the distance to nearest neighbors as a proxy for cell size/spacing.
    The median nearest neighbor distance approximates cell diameter.

    Parameters
    ----------
    adata : AnnData
        Annotated data object
    cell_type_col : str
        Column name for cell type
    indices : array
        KNN indices (n_cells x k)
    distances : array
        KNN distances (n_cells x k)

    Returns
    -------
    DataFrame
        Cell size statistics per cell type:
        - n_cells: number of cells
        - nn1_median: median distance to 1st nearest neighbor (cell diameter proxy)
        - nn1_mean: mean distance to 1st NN
        - nn1_std: std of distance to 1st NN
        - nn1_q25, nn1_q75: quartiles
        - nn_avg_median: median of average distance to all k neighbors
    """
    import numpy as np
    import pandas as pd

    cell_types = adata.obs[cell_type_col].values
    unique_types = sorted(set(cell_types))

    # Distance to 1st nearest neighbor (best proxy for cell diameter)
    nn1_dist = distances[:, 0]

    # Average distance to all k neighbors
    nn_avg_dist = distances.mean(axis=1)

    results = []

    for ct in unique_types:
        ct_mask = cell_types == ct
        n_cells = ct_mask.sum()

        if n_cells < 10:
            continue

        nn1_ct = nn1_dist[ct_mask]
        nn_avg_ct = nn_avg_dist[ct_mask]

        results.append({
            'cell_type': ct,
            'n_cells': n_cells,
            'nn1_median': np.median(nn1_ct),
            'nn1_mean': np.mean(nn1_ct),
            'nn1_std': np.std(nn1_ct),
            'nn1_q25': np.percentile(nn1_ct, 25),
            'nn1_q75': np.percentile(nn1_ct, 75),
            'nn1_min': np.min(nn1_ct),
            'nn1_max': np.max(nn1_ct),
            'nn_avg_median': np.median(nn_avg_ct),
            'nn_avg_mean': np.mean(nn_avg_ct),
        })

    df = pd.DataFrame(results)
    df = df.sort_values('nn1_median')

    return df


def compute_cell_sizes_with_area(adata, cell_type_col, area_col=None):
    """
    Compute cell size statistics including cell area if available.

    Parameters
    ----------
    adata : AnnData
        Annotated data object
    cell_type_col : str
        Column name for cell type
    area_col : str, optional
        Column name for cell area (e.g., 'cell_area', 'nucleus_area')
        If None, searches for common area column names.

    Returns
    -------
    DataFrame
        Cell size statistics including area if available
    """
    import numpy as np
    import pandas as pd

    cell_types = adata.obs[cell_type_col].values
    unique_types = sorted(set(cell_types))

    # Try to find area column
    if area_col is None:
        area_candidates = ['cell_area', 'nucleus_area', 'area', 'Cell_Area', 'Nucleus_Area']
        for col in area_candidates:
            if col in adata.obs.columns:
                area_col = col
                print(f"Found area column: {area_col}")
                break

    has_area = area_col is not None and area_col in adata.obs.columns

    results = []

    for ct in unique_types:
        ct_mask = cell_types == ct
        n_cells = ct_mask.sum()

        if n_cells < 10:
            continue

        row = {
            'cell_type': ct,
            'n_cells': n_cells,
        }

        if has_area:
            areas = adata.obs.loc[ct_mask, area_col].values
            row.update({
                'area_median': np.median(areas),
                'area_mean': np.mean(areas),
                'area_std': np.std(areas),
                'area_q25': np.percentile(areas, 25),
                'area_q75': np.percentile(areas, 75),
                # Estimate diameter from area (assuming circular)
                'diameter_from_area': 2 * np.sqrt(np.median(areas) / np.pi),
            })

        results.append(row)

    df = pd.DataFrame(results)

    if has_area:
        df = df.sort_values('area_median')
    else:
        print("No area column found. Use compute_cell_sizes() with NN distances instead.")

    return df


def summarize_cell_sizes(adata, cell_type_col, indices, distances, area_col=None):
    """
    Comprehensive cell size summary combining NN distances and area (if available).

    Parameters
    ----------
    adata : AnnData
        Annotated data object
    cell_type_col : str
        Column name for cell type
    indices : array
        KNN indices
    distances : array
        KNN distances
    area_col : str, optional
        Column name for cell area

    Returns
    -------
    DataFrame
        Combined cell size statistics
    dict
        Summary statistics across all cells
    """
    import numpy as np
    import pandas as pd

    # Get NN-based sizes
    df_nn = compute_cell_sizes(adata, cell_type_col, indices, distances)

    # Try to get area-based sizes
    df_area = compute_cell_sizes_with_area(adata, cell_type_col, area_col)

    # Merge if area is available
    if 'area_median' in df_area.columns:
        df = df_nn.merge(df_area[['cell_type', 'area_median', 'area_mean', 'diameter_from_area']],
                         on='cell_type', how='left')
    else:
        df = df_nn

    # Global summary
    nn1_all = distances[:, 0]
    summary = {
        'total_cells': len(adata),
        'n_cell_types': len(df),
        'global_nn1_median': np.median(nn1_all),
        'global_nn1_mean': np.mean(nn1_all),
        'global_nn1_std': np.std(nn1_all),
        'recommended_tau': np.median(nn1_all),  # tau ≈ cell diameter
    }

    print(f"\nGlobal Cell Size Summary:")
    print(f"  Total cells: {summary['total_cells']:,}")
    print(f"  Median NN1 distance (cell diameter proxy): {summary['global_nn1_median']:.2f} μm")
    print(f"  Recommended tau: {summary['recommended_tau']:.2f} μm")

    return df, summary


def print_cell_size_table(df, top_n=None):
    """
    Print cell size table in a readable format.

    Parameters
    ----------
    df : DataFrame
        Output from compute_cell_sizes or summarize_cell_sizes
    top_n : int, optional
        Only show top N cell types by count
    """
    if top_n is not None:
        df = df.nlargest(top_n, 'n_cells')

    print(f"\n{'Cell Type':<20} {'N cells':>10} {'NN1 median':>12} {'NN1 mean':>10} {'NN1 std':>10}")
    print("-" * 65)

    for _, row in df.iterrows():
        print(f"{row['cell_type']:<20} {row['n_cells']:>10,} {row['nn1_median']:>12.2f} "
              f"{row['nn1_mean']:>10.2f} {row['nn1_std']:>10.2f}")

    if 'diameter_from_area' in df.columns:
        print(f"\n{'Cell Type':<20} {'Area median':>12} {'Diameter':>10}")
        print("-" * 45)
        for _, row in df.iterrows():
            if pd.notna(row.get('area_median')):
                print(f"{row['cell_type']:<20} {row['area_median']:>12.1f} {row['diameter_from_area']:>10.2f}")


def build_spatial_graph_from_coords(coords, k=10, verbose=True):
    """
    Build k-nearest neighbor spatial graph from raw coordinates.

    This is useful when working with non-AnnData formats (e.g., parquet + h5,
    Xenium raw outputs, etc.).

    Parameters
    ----------
    coords : array-like
        Spatial coordinates, shape (n_cells, 2) or (n_cells, 3)
    k : int
        Number of neighbors (default: 10)
    verbose : bool
        Print summary information (default: True)

    Returns
    -------
    indices : array
        KNN indices (n_cells x k)
    distances : array
        KNN distances (n_cells x k)

    Examples
    --------
    >>> import constellation as cst
    >>> coords = meta[['x_centroid', 'y_centroid']].values
    >>> indices, distances = cst.build_spatial_graph_from_coords(coords, k=10)
    """
    import numpy as np

    coords = np.asarray(coords)

    nn = NearestNeighbors(n_neighbors=k + 1, algorithm='ball_tree')
    nn.fit(coords)
    distances, indices = nn.kneighbors(coords)

    # Remove self (first neighbor)
    distances = distances[:, 1:]
    indices = indices[:, 1:]

    if verbose:
        print(f"Spatial graph: {len(coords):,} cells, k={k} neighbors")
        print(f"  Median NN1 distance: {np.median(distances[:, 0]):.2f}")

    return indices, distances


def filter_lr_pairs_by_genes(lr_pairs, gene_set, verbose=True):
    """
    Filter LR pairs to those where both genes are in a gene set.

    Parameters
    ----------
    lr_pairs : list of tuples
        List of (ligand, receptor) pairs
    gene_set : set or list
        Set of available gene names
    verbose : bool
        Print summary (default: True)

    Returns
    -------
    list of tuples
        Filtered LR pairs

    Examples
    --------
    >>> # Filter LR pairs to genes in Xenium panel
    >>> gene_set = set(gene_names)
    >>> lr_pairs_filtered = cst.filter_lr_pairs_by_genes(lr_pairs, gene_set)
    """
    gene_set = set(gene_set)

    filtered = [(lig, rec) for lig, rec in lr_pairs
                if lig in gene_set and rec in gene_set]

    if verbose:
        print(f"LR pairs: {len(lr_pairs)} total -> {len(filtered)} in gene set")

    return filtered
