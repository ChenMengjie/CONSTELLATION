"""
Spatial LR Testing Functions

Core functions for testing ligand-receptor interactions in spatial transcriptomics data.
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.sparse import issparse
from tqdm import tqdm


# =============================================================================
# INTERNAL HELPERS
# =============================================================================

def _apply_burden_correction(df, total_burden, n_performed, fdr_method):
    """
    Apply FDR correction with optional burden padding.

    If total_burden > n_performed, pads with p=1.0 dummy values so that
    FDR correction accounts for the full testing burden. Handles
    fdr_method='none' for deferred correction in composite functions.

    Parameters
    ----------
    df : DataFrame
        Results with 'p_value' column. Modified in place.
    total_burden : int or None
        Total number of tests in full analysis.
    n_performed : int
        Number of tests actually performed.
    fdr_method : str
        Correction method ('fdr_bh', 'none', etc.).

    Returns
    -------
    DataFrame
        Input df with 'p_adj', 'total_burden', 'n_tests_performed',
        'fdr_method' columns added.
    """
    from statsmodels.stats.multitest import multipletests

    if fdr_method == 'none':
        df['p_adj'] = df['p_value']
        df['total_burden'] = n_performed
    elif total_burden is not None and total_burden > n_performed:
        pvals_padded = list(df['p_value']) + [1.0] * (total_burden - n_performed)
        _, pvals_adj_padded, _, _ = multipletests(pvals_padded, method=fdr_method)
        df['p_adj'] = pvals_adj_padded[:n_performed]
        df['total_burden'] = total_burden
    else:
        if n_performed > 1:
            _, pvals_adj, _, _ = multipletests(df['p_value'], method=fdr_method)
            df['p_adj'] = pvals_adj
        else:
            df['p_adj'] = df['p_value']
        df['total_burden'] = n_performed

    df['n_tests_performed'] = n_performed
    df['fdr_method'] = fdr_method
    return df


def _analytical_test(L, w, n):
    """
    Compute analytical permutation test for T = dot(L, w).

    Under random permutation of L, E[T] = sum(L)*sum(w)/n and
    Var[T] = SS_L * SS_w / (n-1), where SS_x = sum(x^2) - (sum(x))^2/n.

    Parameters
    ----------
    L : array
        Ligand expression vector (log1p-transformed).
    w : array
        Distance-weighted receptor signal vector.
    n : int
        Number of cells.

    Returns
    -------
    dict
        Keys: S_obs, null_mean, null_std, z_score, p_value, fold_enrichment.
    """
    S_obs = float(np.dot(L, w))
    sum_L = L.sum()
    sum_w = w.sum()
    SS_L = (L * L).sum() - sum_L * sum_L / n
    SS_w = (w * w).sum() - sum_w * sum_w / n

    null_mean = sum_L * sum_w / n
    null_var = SS_L * SS_w / (n - 1) if n > 1 else 0.0
    null_std = np.sqrt(null_var) if null_var > 0 else 0.0

    if null_std < 1e-10:
        z = 0.0
        p = 1.0
    else:
        z = (S_obs - null_mean) / null_std
        p = 2 * (1 - stats.norm.cdf(abs(z)))

    if null_mean > 1e-10:
        fold_enrichment = S_obs / null_mean
    else:
        fold_enrichment = np.nan

    return {
        'S_obs': float(S_obs),
        'null_mean': float(null_mean),
        'null_std': float(null_std),
        'z_score': z,
        'p_value': p,
        'fold_enrichment': float(fold_enrichment),
    }


# =============================================================================
# CORE TESTING FUNCTIONS
# =============================================================================

def test_within_type_lr(cell_type, cell_types, X, gene_names, indices, distances,
                        lr_pairs_subset, lr_annotations=None, tau=5.0,
                        X_csc=None, min_cells=50, min_edges=50):
    """
    Test LR pairs within a cell type (autocrine signaling).

    Parameters
    ----------
    cell_type : str
        Cell type to test
    cell_types : array
        Cell type labels for all cells
    X : array or sparse matrix
        Expression matrix (cells x genes)
    gene_names : list
        Gene names
    indices : array
        KNN indices from spatial graph
    distances : array
        KNN distances from spatial graph
    lr_pairs_subset : list of tuples
        List of (ligand, receptor) pairs to test
    lr_annotations : dict, optional
        LR pair annotations {(lig, rec): {'ligand_category': ..., ...}}
    tau : float
        Kernel decay distance (microns)
    min_cells : int
        Minimum cells required for this type (default: 50)
    min_edges : int
        Minimum within-type edges required (default: 50)

    Returns
    -------
    list of dict
        Test results for each LR pair (z_score, p_value from analytical null)
    """
    # Get cell indices for this type
    ct_idx = np.where(cell_types == cell_type)[0]

    if len(ct_idx) < min_cells:
        return []

    # --- Vectorized edge finding (boolean lookup) ---
    n_total = len(cell_types)
    is_ct = np.zeros(n_total, dtype=bool)
    is_ct[ct_idx] = True

    nbr_indices = indices[ct_idx, :]      # (n_cells_ct, K)
    nbr_dists = distances[ct_idx, :]      # (n_cells_ct, K)
    mask = is_ct[nbr_indices]             # O(1) per element
    row_pos, col_pos = np.where(mask)

    n_edges = len(row_pos)
    if n_edges < min_edges:
        return []

    edge_dists = nbr_dists[row_pos, col_pos]

    # Kernel weights
    K = np.exp(-edge_dists / tau)

    # Local indices: row_pos is already local (index into ct_idx)
    local_edges_i = row_pos
    # For j: map global neighbor IDs to local position in ct_idx
    edges_j_global = nbr_indices[row_pos, col_pos]
    local_edges_j = np.searchsorted(ct_idx, edges_j_global)

    # --- Gene name lookup dict ---
    gene_name_to_idx = {g: i for i, g in enumerate(gene_names)}

    # --- Validate LR pairs ---
    valid_pairs = []
    for lig, rec in lr_pairs_subset:
        li = gene_name_to_idx.get(lig)
        ri = gene_name_to_idx.get(rec)
        if li is not None and ri is not None:
            valid_pairs.append((lig, rec, li, ri))

    if not valid_pairs:
        return []

    # --- Use CSC for fast column access ---
    if X_csc is None and issparse(X):
        X_csc = X.tocsc()

    n_ct = len(ct_idx)
    if lr_annotations is None:
        lr_annotations = {}

    # --- Cache receptor w-vectors (receptor may pair with multiple ligands) ---
    receptor_cache = {}  # rec_gi -> w
    unique_receptors = set((rec, rec_gi) for _, rec, _, rec_gi in valid_pairs)

    for rec, rec_gi in unique_receptors:
        if X_csc is not None:
            R_ct = np.log1p(np.asarray(X_csc[ct_idx, rec_gi].todense()).flatten())
        else:
            R_ct = np.log1p(X[ct_idx, rec_gi])
        KR = K * R_ct[local_edges_j]
        w = np.zeros(n_ct)
        np.add.at(w, local_edges_i, KR)
        receptor_cache[rec_gi] = w

    # --- Cache ligand vectors ---
    ligand_cache = {}  # lig_gi -> L_ct
    unique_ligands = set((lig, lig_gi) for lig, _, lig_gi, _ in valid_pairs)

    for lig, lig_gi in unique_ligands:
        if X_csc is not None:
            L_ct = np.log1p(np.asarray(X_csc[ct_idx, lig_gi].todense()).flatten())
        else:
            L_ct = np.log1p(X[ct_idx, lig_gi])
        ligand_cache[lig_gi] = L_ct

    # --- Test each LR pair using analytical permutation null ---
    results = []

    for lig, rec, lig_gi, rec_gi in valid_pairs:
        L_ct = ligand_cache[lig_gi]
        w = receptor_cache[rec_gi]

        test_result = _analytical_test(L_ct, w, n_ct)
        ann = lr_annotations.get((lig, rec), {})

        results.append({
            'ligand': lig,
            'receptor': rec,
            'cell_type': cell_type,
            'n_cells': len(ct_idx),
            'n_edges': n_edges,
            **test_result,
            'ligand_category': ann.get('ligand_category'),
            'ligand_category_major': ann.get('ligand_category_major', 'Other'),
            'receptor_category': ann.get('receptor_category')
        })

    return results


def test_between_type_lr(sender_type, receiver_type, cell_types, X, gene_names,
                         indices, distances, lr_pairs_subset,
                         ct_indices=None,
                         lr_annotations=None, tau=5.0,
                         X_csc=None, min_cells=10, min_edges=100):
    """
    Test LR pairs between sender and receiver cell types.

    Parameters
    ----------
    sender_type : str
        Sender cell type (expresses ligand)
    receiver_type : str
        Receiver cell type (expresses receptor)
    cell_types : array
        Cell type labels for all cells
    X : array or sparse matrix
        Expression matrix (cells x genes)
    gene_names : list
        Gene names
    indices : array
        KNN indices from spatial graph
    distances : array
        KNN distances from spatial graph
    lr_pairs_subset : list of tuples
        List of (ligand, receptor) pairs to test
    ct_indices : dict, optional
        Pre-computed cell type indices {cell_type: array of indices}
    lr_annotations : dict, optional
        LR pair annotations
    tau : float
        Kernel decay distance (microns)
    min_cells : int
        Minimum cells per type (default: 10)
    min_edges : int
        Minimum sender-receiver edges required (default: 100)

    Returns
    -------
    list of dict
        Test results for each LR pair (z_score, p_value from analytical null)
    """
    # Get cell indices
    if ct_indices is not None:
        sender_idx = ct_indices[sender_type]
        receiver_idx = ct_indices[receiver_type]
    else:
        sender_idx = np.where(cell_types == sender_type)[0]
        receiver_idx = np.where(cell_types == receiver_type)[0]

    if len(sender_idx) < min_cells or len(receiver_idx) < min_cells:
        return []

    # --- Vectorized edge finding (boolean lookup) ---
    n_total = len(cell_types)
    is_receiver = np.zeros(n_total, dtype=bool)
    is_receiver[receiver_idx] = True

    nbr_indices = indices[sender_idx, :]     # (n_sender, K)
    nbr_dists = distances[sender_idx, :]     # (n_sender, K)
    mask = is_receiver[nbr_indices]          # O(1) per element
    row_pos, col_pos = np.where(mask)

    n_edges = len(row_pos)
    if n_edges < min_edges:
        return []

    edge_dists = nbr_dists[row_pos, col_pos]
    K = np.exp(-edge_dists / tau)

    # Local indices
    local_edges_i = row_pos   # already local index into sender_idx
    edges_j_global = nbr_indices[row_pos, col_pos]
    # receiver_idx comes from np.where so is sorted
    local_edges_j = np.searchsorted(receiver_idx, edges_j_global)

    # --- Gene name lookup dict ---
    gene_name_to_idx = {g: i for i, g in enumerate(gene_names)}

    # --- Validate LR pairs ---
    valid_pairs = []
    for lig, rec in lr_pairs_subset:
        li = gene_name_to_idx.get(lig)
        ri = gene_name_to_idx.get(rec)
        if li is not None and ri is not None:
            valid_pairs.append((lig, rec, li, ri))

    if not valid_pairs:
        return []

    # --- Use CSC for fast column access ---
    if X_csc is None and issparse(X):
        X_csc = X.tocsc()

    n_sender = len(sender_idx)
    if lr_annotations is None:
        lr_annotations = {}

    # --- Cache receptor w-vectors (receptor may pair with multiple ligands) ---
    receptor_cache = {}  # rec_gi -> (w, sum_w, SS_w)
    unique_receptors = set((rec, rec_gi) for _, rec, _, rec_gi in valid_pairs)

    for rec, rec_gi in unique_receptors:
        if X_csc is not None:
            R_receiver = np.log1p(np.asarray(X_csc[receiver_idx, rec_gi].todense()).flatten())
        else:
            R_receiver = np.log1p(X[receiver_idx, rec_gi])
        KR = K * R_receiver[local_edges_j]
        w = np.zeros(n_sender)
        np.add.at(w, local_edges_i, KR)
        receptor_cache[rec_gi] = w

    # --- Cache ligand vectors ---
    ligand_cache = {}  # lig_gi -> L_sender
    unique_ligands = set((lig, lig_gi) for lig, _, lig_gi, _ in valid_pairs)

    for lig, lig_gi in unique_ligands:
        if X_csc is not None:
            L_sender = np.log1p(np.asarray(X_csc[sender_idx, lig_gi].todense()).flatten())
        else:
            L_sender = np.log1p(X[sender_idx, lig_gi])
        ligand_cache[lig_gi] = L_sender

    # --- Test each LR pair using analytical permutation null ---
    results = []

    for lig, rec, lig_gi, rec_gi in valid_pairs:
        L_sender = ligand_cache[lig_gi]
        w = receptor_cache[rec_gi]

        test_result = _analytical_test(L_sender, w, n_sender)
        ann = lr_annotations.get((lig, rec), {})

        results.append({
            'ligand': lig,
            'receptor': rec,
            'sender': sender_type,
            'receiver': receiver_type,
            'n_edges': n_edges,
            **test_result,
            'ligand_category': ann.get('ligand_category'),
            'ligand_category_major': ann.get('ligand_category_major', 'Other'),
            'receptor_category': ann.get('receptor_category')
        })

    return results


def compute_expression_fractions(adata, cell_type_col, min_expr_frac=0.01):
    """
    Compute expression fractions for each gene in each cell type.

    Parameters
    ----------
    adata : AnnData
        Annotated data object
    cell_type_col : str
        Column name for cell type labels
    min_expr_frac : float
        Minimum expression fraction threshold

    Returns
    -------
    dict
        {cell_type: {gene: fraction}}
    """
    cell_types = adata.obs[cell_type_col].values
    unique_types = list(set(cell_types))
    X = adata.X
    gene_names = list(adata.var_names)

    expr_frac = {}
    for ct in unique_types:
        idx = np.where(cell_types == ct)[0]
        if issparse(X):
            X_ct = X[idx, :].toarray()
        else:
            X_ct = X[idx, :]
        frac = (X_ct > 0).sum(axis=0) / len(idx)
        expr_frac[ct] = dict(zip(gene_names, np.array(frac).flatten()))

    return expr_frac


def filter_testable_lr_pairs(lr_pairs, expr_frac, cell_type, min_expr_frac=0.01):
    """
    Filter LR pairs to those expressed in a cell type.

    Parameters
    ----------
    lr_pairs : list of tuples
        List of (ligand, receptor) pairs
    expr_frac : dict
        Expression fractions {cell_type: {gene: fraction}}
    cell_type : str
        Cell type to check
    min_expr_frac : float
        Minimum expression fraction

    Returns
    -------
    list of tuples
        Filtered LR pairs
    """
    testable = []
    ct_frac = expr_frac[cell_type]
    for lig, rec in lr_pairs:
        if ct_frac.get(lig, 0) >= min_expr_frac and ct_frac.get(rec, 0) >= min_expr_frac:
            testable.append((lig, rec))
    return testable


def filter_testable_lr_pairs_between(lr_pairs, expr_frac, sender, receiver, min_expr_frac=0.01):
    """
    Filter LR pairs for between-type analysis.

    Parameters
    ----------
    lr_pairs : list of tuples
        List of (ligand, receptor) pairs
    expr_frac : dict
        Expression fractions
    sender : str
        Sender cell type
    receiver : str
        Receiver cell type
    min_expr_frac : float
        Minimum expression fraction

    Returns
    -------
    list of tuples
        Filtered LR pairs where ligand expressed in sender, receptor in receiver
    """
    testable = []
    for lig, rec in lr_pairs:
        if (expr_frac[sender].get(lig, 0) >= min_expr_frac and
            expr_frac[receiver].get(rec, 0) >= min_expr_frac):
            testable.append((lig, rec))
    return testable


# =============================================================================
# TARGETED TESTING FUNCTIONS (with full burden correction)
# =============================================================================

def test_lr_pair_within_types(ligand, receptor, cell_types_array, X, gene_names,
                               indices, distances,
                               total_burden=None,
                               testable_cell_types=None,
                               burden_info=None,
                               tau=None, fdr_method='fdr_bh'):
    """
    Test a specific LR pair across cell types (within-type analysis).

    FDR correction uses the FULL testable burden for conservative correction.

    Parameters
    ----------
    ligand : str
        Ligand gene name
    receptor : str
        Receptor gene name
    cell_types_array : array
        Cell type labels for all cells
    X : array or sparse matrix
        Expression matrix (cells x genes)
    gene_names : list
        Gene names
    indices : array
        KNN indices from spatial graph
    distances : array
        KNN distances from spatial graph
    total_burden : int, optional
        Total number of tests in full analysis (for FDR correction).
        If None and burden_info provided, uses burden_info['total_tests'].
    testable_cell_types : list, optional
        Cell types to test. If None and burden_info provided, uses burden_info['testable_cell_types'].
    burden_info : dict, optional
        Output from compute_testable_burden_within. Provides total_burden, testable_cell_types, and tau.
    tau : float, optional
        Kernel decay distance (microns). Priority order:
        1. User-supplied value (if provided)
        2. burden_info['tau'] (if burden_info provided)
        3. Computed from median NN1 distance (data-driven default)
    fdr_method : str
        Multiple testing correction method

    Returns
    -------
    DataFrame
        Results with p_adj corrected using full burden
    """
    if ligand not in gene_names or receptor not in gene_names:
        raise ValueError(f"Gene not found: {ligand if ligand not in gene_names else receptor}")

    # Extract from burden_info if provided (only if user didn't supply values)
    if burden_info is not None:
        if total_burden is None:
            total_burden = burden_info.get('total_tests')
        if testable_cell_types is None:
            testable_cell_types = burden_info.get('testable_cell_types')
        if tau is None:  # User-supplied tau takes precedence
            tau = burden_info.get('tau')

    # Default tau from data if still None
    if tau is None:
        tau = float(np.median(distances[:, 0]))

    # Determine cell types to test
    unique_types = list(set(cell_types_array))
    if testable_cell_types is None:
        testable_cell_types = [ct for ct in unique_types
                               if (cell_types_array == ct).sum() >= 50]

    results = []

    for cell_type in tqdm(testable_cell_types, desc=f"Testing {ligand}-{receptor}"):
        ct_results = test_within_type_lr(
            cell_type=cell_type,
            cell_types=cell_types_array,
            X=X,
            gene_names=gene_names,
            indices=indices,
            distances=distances,
            lr_pairs_subset=[(ligand, receptor)],
            tau=tau
        )
        results.extend(ct_results)

    if not results:
        return pd.DataFrame()

    df = pd.DataFrame(results)
    n_performed = len(df)

    # Apply FDR correction
    df = _apply_burden_correction(df, total_burden, n_performed, fdr_method)

    return df.sort_values('z_score', ascending=False)


def test_lr_pair_between_types(ligand, receptor, cell_types_array, X, gene_names,
                                indices, distances,
                                total_burden=None,
                                testable_pairs=None,
                                ct_indices=None,
                                burden_info=None,
                                tau=None, fdr_method='fdr_bh'):
    """
    Test a specific LR pair across cell type pairs (between-type analysis).

    FDR correction uses the FULL testable burden for conservative correction.

    Parameters
    ----------
    ligand : str
        Ligand gene name
    receptor : str
        Receptor gene name
    cell_types_array : array
        Cell type labels for all cells
    X : array or sparse matrix
        Expression matrix
    gene_names : list
        Gene names
    indices : array
        KNN indices
    distances : array
        KNN distances
    total_burden : int, optional
        Total number of tests in full analysis (for FDR correction).
        If None and burden_info provided, uses burden_info['total_tests'].
    testable_pairs : list of tuples, optional
        List of (sender, receiver) pairs to test.
        If None and burden_info provided, extracts from burden_info.
    ct_indices : dict, optional
        Pre-computed {cell_type: indices}
    burden_info : dict, optional
        Output from compute_testable_burden_between. Provides total_burden, testable_pairs, ct_indices, and tau.
    tau : float, optional
        Kernel decay distance (microns). Priority order:
        1. User-supplied value (if provided)
        2. burden_info['tau'] (if burden_info provided)
        3. Computed from median NN1 distance (data-driven default)
    fdr_method : str
        Multiple testing correction method

    Returns
    -------
    DataFrame
        Results with p_adj corrected using full burden
    """
    if ligand not in gene_names or receptor not in gene_names:
        raise ValueError(f"Gene not found: {ligand if ligand not in gene_names else receptor}")

    # Extract from burden_info if provided (only if user didn't supply values)
    if burden_info is not None:
        if total_burden is None:
            total_burden = burden_info.get('total_tests')
        if testable_pairs is None:
            testable_pairs = [(s, r) for s, r, _, _ in burden_info.get('testable_pairs', [])]
        if ct_indices is None:
            ct_indices = burden_info.get('ct_indices')
        if tau is None:  # User-supplied tau takes precedence
            tau = burden_info.get('tau')

    # Default tau from data if still None
    if tau is None:
        tau = float(np.median(distances[:, 0]))

    # Build cell type indices if not provided
    unique_types = list(set(cell_types_array))
    if ct_indices is None:
        ct_indices = {ct: np.where(cell_types_array == ct)[0] for ct in unique_types}

    # Determine pairs to test
    if testable_pairs is None:
        type_counts = {ct: len(ct_indices[ct]) for ct in unique_types}
        valid_types = [ct for ct in unique_types if type_counts[ct] >= 10]
        testable_pairs = [(s, r) for s in valid_types for r in valid_types if s != r]

    results = []

    for sender, receiver in tqdm(testable_pairs, desc=f"Testing {ligand}-{receptor}"):
        ct_results = test_between_type_lr(
            sender_type=sender,
            receiver_type=receiver,
            cell_types=cell_types_array,
            X=X,
            gene_names=gene_names,
            indices=indices,
            distances=distances,
            lr_pairs_subset=[(ligand, receptor)],
            ct_indices=ct_indices,
            tau=tau
        )
        results.extend(ct_results)

    if not results:
        return pd.DataFrame()

    df = pd.DataFrame(results)
    n_performed = len(df)

    # Apply FDR correction
    df = _apply_burden_correction(df, total_burden, n_performed, fdr_method)

    return df.sort_values('z_score', ascending=False)


def test_ligand_all_receptors(ligand, lr_pairs, cell_types_array, X, gene_names,
                               indices, distances, analysis_type='within',
                               total_burden=None,
                               burden_info=None,
                               tau=None, fdr_method='fdr_bh'):
    """
    Test all receptor pairs for a given ligand.

    FDR correction uses the FULL testable burden for conservative correction.

    Parameters
    ----------
    ligand : str
        Ligand gene name
    lr_pairs : list of tuples
        All available LR pairs [(ligand, receptor), ...]
    cell_types_array : array
        Cell type labels
    X : array or sparse matrix
        Expression matrix
    gene_names : list
        Gene names
    indices : array
        KNN indices
    distances : array
        KNN distances
    analysis_type : str
        'within' or 'between'
    total_burden : int, optional
        Total tests in full analysis. If None, extracted from burden_info.
    burden_info : dict, optional
        Output from compute_testable_burden_within/between.
        Provides testable cell types/pairs and total burden.
    tau : float
        Kernel decay distance
    fdr_method : str
        Multiple testing correction method

    Returns
    -------
    DataFrame
        Results with p_adj corrected using full burden
    """
    # Find all receptors for this ligand
    receptors = [rec for lig, rec in lr_pairs if lig == ligand and rec in gene_names]

    if not receptors:
        print(f"No receptors found for ligand {ligand}")
        return pd.DataFrame()

    print(f"Testing {ligand} with {len(receptors)} receptors: {receptors}")

    # Extract info from burden_info if provided
    if burden_info is not None:
        if total_burden is None:
            total_burden = burden_info['total_tests']
        if tau is None:
            tau = burden_info.get('tau')
        if analysis_type == 'within':
            testable_cell_types = burden_info.get('testable_cell_types')
        else:
            testable_pairs = [(s, r) for s, r, _, _ in burden_info.get('testable_pairs', [])]
            ct_indices = burden_info.get('ct_indices')
    else:
        testable_cell_types = None
        testable_pairs = None
        ct_indices = None

    # Default tau from data if still None
    if tau is None:
        tau = float(np.median(distances[:, 0]))

    all_results = []

    for receptor in tqdm(receptors, desc=f"Testing {ligand} → receptors"):
        if analysis_type == 'within':
            df = test_lr_pair_within_types(
                ligand, receptor, cell_types_array, X, gene_names,
                indices, distances,
                total_burden=None,  # Don't correct yet
                testable_cell_types=testable_cell_types,
                tau=tau, fdr_method='none'
            )
        else:
            df = test_lr_pair_between_types(
                ligand, receptor, cell_types_array, X, gene_names,
                indices, distances,
                total_burden=None,
                testable_pairs=testable_pairs,
                ct_indices=ct_indices,
                tau=tau, fdr_method='none'
            )

        if len(df) > 0:
            all_results.append(df)

    if not all_results:
        return pd.DataFrame()

    combined = pd.concat(all_results, ignore_index=True)
    n_performed = len(combined)

    # Apply FDR correction
    combined = _apply_burden_correction(combined, total_burden, n_performed, fdr_method)

    return combined.sort_values('z_score', ascending=False)


def test_receptor_all_ligands(receptor, lr_pairs, cell_types_array, X, gene_names,
                               indices, distances, analysis_type='within',
                               total_burden=None,
                               burden_info=None,
                               tau=None, fdr_method='fdr_bh'):
    """
    Test all ligand pairs for a given receptor.

    FDR correction uses the FULL testable burden for conservative correction.

    Parameters
    ----------
    receptor : str
        Receptor gene name
    lr_pairs : list of tuples
        All available LR pairs
    cell_types_array : array
        Cell type labels
    X : array or sparse matrix
        Expression matrix
    gene_names : list
        Gene names
    indices : array
        KNN indices
    distances : array
        KNN distances
    analysis_type : str
        'within' or 'between'
    total_burden : int, optional
        Total tests in full analysis
    burden_info : dict, optional
        Output from compute_testable_burden_within/between
    tau : float
        Kernel decay distance
    fdr_method : str
        Multiple testing correction method

    Returns
    -------
    DataFrame
        Results with p_adj corrected using full burden
    """
    # Find all ligands for this receptor
    ligands = [lig for lig, rec in lr_pairs if rec == receptor and lig in gene_names]

    if not ligands:
        print(f"No ligands found for receptor {receptor}")
        return pd.DataFrame()

    print(f"Testing {receptor} with {len(ligands)} ligands: {ligands}")

    # Extract info from burden_info
    if burden_info is not None:
        if total_burden is None:
            total_burden = burden_info['total_tests']
        if tau is None:
            tau = burden_info.get('tau')
        if analysis_type == 'within':
            testable_cell_types = burden_info.get('testable_cell_types')
        else:
            testable_pairs = [(s, r) for s, r, _, _ in burden_info.get('testable_pairs', [])]
            ct_indices = burden_info.get('ct_indices')
    else:
        testable_cell_types = None
        testable_pairs = None
        ct_indices = None

    # Default tau from data if still None
    if tau is None:
        tau = float(np.median(distances[:, 0]))

    all_results = []

    for ligand in tqdm(ligands, desc=f"Testing ligands → {receptor}"):
        if analysis_type == 'within':
            df = test_lr_pair_within_types(
                ligand, receptor, cell_types_array, X, gene_names,
                indices, distances,
                total_burden=None,
                testable_cell_types=testable_cell_types,
                tau=tau, fdr_method='none'
            )
        else:
            df = test_lr_pair_between_types(
                ligand, receptor, cell_types_array, X, gene_names,
                indices, distances,
                total_burden=None,
                testable_pairs=testable_pairs,
                ct_indices=ct_indices,
                tau=tau, fdr_method='none'
            )

        if len(df) > 0:
            all_results.append(df)

    if not all_results:
        return pd.DataFrame()

    combined = pd.concat(all_results, ignore_index=True)
    n_performed = len(combined)

    if total_burden is not None and total_burden > n_performed:
        pvals_padded = list(combined['p_value']) + [1.0] * (total_burden - n_performed)
        _, pvals_adj_padded, _, _ = multipletests(pvals_padded, method=fdr_method)
        combined['p_adj'] = pvals_adj_padded[:n_performed]
        combined['total_burden'] = total_burden
    else:
        if n_performed > 1:
            _, pvals_adj, _, _ = multipletests(combined['p_value'], method=fdr_method)
            combined['p_adj'] = pvals_adj
        else:
            combined['p_adj'] = combined['p_value']
        combined['total_burden'] = n_performed

    combined['n_tests_performed'] = n_performed
    combined['fdr_method'] = fdr_method

    return combined.sort_values('z_score', ascending=False)


def test_custom_lr_set(lr_subset, cell_types_array, X, gene_names,
                        indices, distances, analysis_type='within',
                        total_burden=None,
                        burden_info=None,
                        tau=None, fdr_method='fdr_bh'):
    """
    Test a custom set of LR pairs with full burden correction.

    FDR correction uses the FULL testable burden for conservative correction.

    Parameters
    ----------
    lr_subset : list of tuples
        LR pairs to test [(ligand, receptor), ...]
    cell_types_array : array
        Cell type labels
    X : array or sparse matrix
        Expression matrix
    gene_names : list
        Gene names
    indices : array
        KNN indices
    distances : array
        KNN distances
    analysis_type : str
        'within' or 'between'
    total_burden : int, optional
        Total tests in full analysis
    burden_info : dict, optional
        Output from compute_testable_burden_within/between
    tau : float
        Kernel decay distance
    fdr_method : str
        Multiple testing correction method

    Returns
    -------
    DataFrame
        Results with p_adj corrected using full burden
    """
    # Filter to valid pairs
    valid_pairs = [(l, r) for l, r in lr_subset
                   if l in gene_names and r in gene_names]

    if not valid_pairs:
        print("No valid LR pairs found")
        return pd.DataFrame()

    print(f"Testing {len(valid_pairs)} LR pairs")

    # Extract info from burden_info
    if burden_info is not None:
        if total_burden is None:
            total_burden = burden_info['total_tests']
        if tau is None:
            tau = burden_info.get('tau')
        if analysis_type == 'within':
            testable_cell_types = burden_info.get('testable_cell_types')
        else:
            testable_pairs = [(s, r) for s, r, _, _ in burden_info.get('testable_pairs', [])]
            ct_indices = burden_info.get('ct_indices')
    else:
        testable_cell_types = None
        testable_pairs = None
        ct_indices = None

    # Default tau from data if still None
    if tau is None:
        tau = float(np.median(distances[:, 0]))

    all_results = []

    for ligand, receptor in tqdm(valid_pairs, desc="Testing custom LR set"):
        if analysis_type == 'within':
            df = test_lr_pair_within_types(
                ligand, receptor, cell_types_array, X, gene_names,
                indices, distances,
                total_burden=None,
                testable_cell_types=testable_cell_types,
                tau=tau, fdr_method='none'
            )
        else:
            df = test_lr_pair_between_types(
                ligand, receptor, cell_types_array, X, gene_names,
                indices, distances,
                total_burden=None,
                testable_pairs=testable_pairs,
                ct_indices=ct_indices,
                tau=tau, fdr_method='none'
            )

        if len(df) > 0:
            all_results.append(df)

    if not all_results:
        return pd.DataFrame()

    combined = pd.concat(all_results, ignore_index=True)
    n_performed = len(combined)

    print(f"Tests performed: {n_performed}, Total burden: {total_burden or n_performed}")

    if total_burden is not None and total_burden > n_performed:
        pvals_padded = list(combined['p_value']) + [1.0] * (total_burden - n_performed)
        _, pvals_adj_padded, _, _ = multipletests(pvals_padded, method=fdr_method)
        combined['p_adj'] = pvals_adj_padded[:n_performed]
        combined['total_burden'] = total_burden
    else:
        if n_performed > 1:
            _, pvals_adj, _, _ = multipletests(combined['p_value'], method=fdr_method)
            combined['p_adj'] = pvals_adj
        else:
            combined['p_adj'] = combined['p_value']
        combined['total_burden'] = n_performed

    combined['n_tests_performed'] = n_performed
    combined['fdr_method'] = fdr_method

    return combined.sort_values('z_score', ascending=False)


# =============================================================================
# INPUT VALIDATION
# =============================================================================

def validate_inputs(adata, cell_type_col, lr_pairs, indices, distances,
                    tau=5.0, verbose=True):
    """
    Validate inputs for CONSTELLATION analysis.

    Checks data format, expression matrix, spatial coordinates, LR pair
    overlap, spatial graph structure, and tau parameter. Raises ValueError
    for critical issues and prints warnings for potential problems.

    Parameters
    ----------
    adata : AnnData
        Annotated data with expression matrix, cell type annotations, and
        spatial coordinates.
    cell_type_col : str
        Column in adata.obs containing cell type labels.
    lr_pairs : list of tuples
        List of (ligand, receptor) pairs.
    indices : array, shape (n_cells, k)
        KNN indices from spatial graph.
    distances : array, shape (n_cells, k)
        KNN distances from spatial graph.
    tau : float
        Kernel decay distance in microns (default: 5.0).
    verbose : bool
        Print validation summary (default: True).

    Returns
    -------
    dict
        Validation summary with keys: n_cells, n_genes, n_cell_types,
        n_lr_testable, median_nn1, coord_range, is_valid.

    Raises
    ------
    ValueError
        If critical issues are found (missing fields, wrong shapes, etc.).
    """
    import warnings

    n_cells = adata.n_obs
    info = {'is_valid': True}
    warn_msgs = []

    # ── 1. AnnData structure ──
    # Expression matrix
    if adata.X is None:
        raise ValueError("adata.X is None. Provide an expression matrix.")
    if not (isinstance(adata.X, np.ndarray) or issparse(adata.X)):
        raise ValueError(
            f"adata.X has unexpected type {type(adata.X).__name__}. "
            "Expected numpy ndarray or scipy sparse matrix.")

    # Cell type column
    if cell_type_col not in adata.obs.columns:
        raise ValueError(
            f"Column '{cell_type_col}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}")
    ct_values = adata.obs[cell_type_col]
    if ct_values.isna().any():
        n_na = ct_values.isna().sum()
        raise ValueError(
            f"{n_na} NaN values in adata.obs['{cell_type_col}']. "
            "Remove or fill NaN cell type labels before running.")
    unique_types = ct_values.unique()
    if len(unique_types) < 2:
        raise ValueError(
            f"Only {len(unique_types)} unique cell type(s) found. "
            "Need at least 2 for meaningful interaction testing.")
    info['n_cell_types'] = len(unique_types)

    # Gene names
    if len(adata.var_names) == 0:
        raise ValueError("adata.var_names is empty. No genes found.")
    info['n_genes'] = len(adata.var_names)
    info['n_cells'] = n_cells

    # Spatial coordinates
    if 'spatial' not in adata.obsm:
        raise ValueError(
            "'spatial' not found in adata.obsm. "
            "Add spatial coordinates: adata.obsm['spatial'] = coords")
    coords = np.asarray(adata.obsm['spatial'])
    if coords.shape[0] != n_cells:
        raise ValueError(
            f"Spatial coordinates have {coords.shape[0]} rows but "
            f"adata has {n_cells} cells.")
    if coords.shape[1] not in (2, 3):
        raise ValueError(
            f"Spatial coordinates have {coords.shape[1]} columns. "
            "Expected 2 (x, y) or 3 (x, y, z).")

    # ── 2. Expression matrix checks ──
    if issparse(adata.X):
        x_data = adata.X.data
    else:
        x_data = adata.X.ravel()

    if len(x_data) == 0 or np.all(x_data == 0):
        raise ValueError("Expression matrix is all zeros.")
    if np.any(x_data < 0):
        warn_msgs.append(
            "WARNING: Negative values detected in expression matrix. "
            "CONSTELLATION expects raw counts (>= 0).")

    # Check for log-transformed data
    x_max = np.max(x_data) if len(x_data) > 0 else 0
    if x_max < 20:
        # Sample to check for integers
        sample = x_data[:min(100000, len(x_data))]
        frac_integer = np.mean(np.equal(np.mod(sample, 1), 0))
        if frac_integer < 0.5:
            warn_msgs.append(
                f"WARNING: Expression values look log-transformed "
                f"(max={x_max:.1f}, {frac_integer:.0%} integers). "
                "CONSTELLATION applies log1p internally — input should be "
                "raw counts.")

    # ── 3. Spatial coordinates ──
    if np.any(np.isnan(coords)) or np.any(np.isinf(coords)):
        raise ValueError(
            "NaN or Inf values in spatial coordinates. Clean coordinates first.")

    # Compute NN1 distance from provided distances
    nn1_dists = distances[:, 0]
    median_nn1 = float(np.median(nn1_dists[nn1_dists > 0]))
    info['median_nn1'] = median_nn1

    coord_min = coords.min(axis=0)
    coord_max = coords.max(axis=0)
    info['coord_range'] = {
        'min': coord_min.tolist(),
        'max': coord_max.tolist(),
        'extent': (coord_max - coord_min).tolist(),
    }

    if median_nn1 < 1.0:
        warn_msgs.append(
            f"WARNING: Median nearest-neighbor distance is {median_nn1:.4f}. "
            "This is very small — coordinates may be normalized or in pixels "
            "instead of microns. CONSTELLATION expects coordinates in µm.")

    # ── 4. LR pair overlap ──
    gene_set = set(adata.var_names)
    all_ligands = set(lig for lig, _ in lr_pairs)
    all_receptors = set(rec for _, rec in lr_pairs)
    found_lig = all_ligands & gene_set
    found_rec = all_receptors & gene_set

    testable_pairs = [(l, r) for l, r in lr_pairs
                      if l in gene_set and r in gene_set]
    info['n_lr_testable'] = len(testable_pairs)
    info['n_lr_total'] = len(lr_pairs)
    info['n_ligands_found'] = len(found_lig)
    info['n_receptors_found'] = len(found_rec)

    if len(testable_pairs) == 0:
        raise ValueError(
            "No LR pairs have both ligand and receptor in adata.var_names. "
            f"Ligands found: {len(found_lig)}/{len(all_ligands)}, "
            f"Receptors found: {len(found_rec)}/{len(all_receptors)}. "
            "Check that gene names match (e.g., HUGO symbols).")

    lig_frac = len(found_lig) / max(len(all_ligands), 1)
    rec_frac = len(found_rec) / max(len(all_receptors), 1)
    if lig_frac < 0.1 or rec_frac < 0.1:
        warn_msgs.append(
            f"WARNING: Low LR gene overlap — "
            f"ligands: {len(found_lig)}/{len(all_ligands)} ({lig_frac:.0%}), "
            f"receptors: {len(found_rec)}/{len(all_receptors)} ({rec_frac:.0%}). "
            "Check that gene name format matches (e.g., HUGO symbols vs Ensembl IDs).")

    # ── 5. Spatial graph ──
    indices = np.asarray(indices)
    distances = np.asarray(distances)

    if indices.shape[0] != n_cells:
        raise ValueError(
            f"indices has {indices.shape[0]} rows but adata has {n_cells} cells.")
    if distances.shape[0] != n_cells:
        raise ValueError(
            f"distances has {distances.shape[0]} rows but adata has {n_cells} cells.")
    if indices.shape != distances.shape:
        raise ValueError(
            f"indices shape {indices.shape} != distances shape {distances.shape}.")
    if np.any(indices < 0) or np.any(indices >= n_cells):
        raise ValueError(
            f"indices contain values outside [0, {n_cells}). "
            "Spatial graph indices must reference valid cell positions.")
    if np.any(np.isnan(distances)) or np.any(np.isinf(distances)):
        raise ValueError("NaN or Inf values in distances array.")
    if np.any(distances < 0):
        raise ValueError("Negative values in distances array.")

    info['k_neighbors'] = indices.shape[1]

    # ── 6. Tau sanity ──
    if tau <= 0:
        raise ValueError(f"tau must be positive, got {tau}.")
    if median_nn1 > 0:
        if tau < median_nn1 * 0.1:
            warn_msgs.append(
                f"WARNING: tau={tau} is much smaller than median NN1 distance "
                f"({median_nn1:.1f}). Effective range 3*tau={3*tau:.1f} may be "
                "too small to capture neighbor interactions.")
        if tau > median_nn1 * 100:
            warn_msgs.append(
                f"WARNING: tau={tau} is much larger than median NN1 distance "
                f"({median_nn1:.1f}). Kernel weights will be nearly uniform — "
                "spatial structure may be washed out.")

    # ── Print summary ──
    if warn_msgs:
        info['is_valid'] = False  # warnings present

    if verbose:
        print(f'\n{"=" * 50}')
        print(f'CONSTELLATION INPUT VALIDATION')
        print(f'{"=" * 50}')
        print(f'Cells:          {n_cells:,}')
        print(f'Genes:          {info["n_genes"]:,}')
        print(f'Cell types:     {info["n_cell_types"]} '
              f'({", ".join(sorted(str(t) for t in unique_types[:5]))}{"..." if len(unique_types) > 5 else ""})')
        print(f'Neighbors (k):  {info["k_neighbors"]}')
        print(f'Median NN1:     {median_nn1:.2f}')
        print(f'Tau:            {tau} (effective range ≈ {3*tau:.0f})')
        print(f'LR pairs:       {len(testable_pairs)}/{len(lr_pairs)} testable '
              f'({len(found_lig)} ligands, {len(found_rec)} receptors in data)')
        extent = info['coord_range']['extent']
        print(f'Coord range:    {extent[0]:.0f} x {extent[1]:.0f}')

        if warn_msgs:
            print(f'\n⚠ {len(warn_msgs)} warning(s):')
            for msg in warn_msgs:
                print(f'  {msg}')
        else:
            print(f'\nAll checks passed.')
        print(f'{"=" * 50}\n')

    return info


# =============================================================================
# CELL-TYPE PAIR SCANNING
# =============================================================================

def scan_celltype_pairs(adata, cell_type_col, lr_pairs, indices, distances,
                        min_expr_frac=0.05, min_cells=50, verbose=True):
    """
    Scan all cell-type pairs and report LR gene overlap and testing burden.

    Use this before running the full analysis to understand the scope of
    testable interactions and estimate computational cost.

    Parameters
    ----------
    adata : AnnData
        Annotated data with expression matrix and cell type annotations.
    cell_type_col : str
        Column in adata.obs containing cell type labels.
    lr_pairs : list of tuples
        List of (ligand, receptor) pairs.
    indices : array, shape (n_cells, k)
        KNN indices from spatial graph.
    distances : array, shape (n_cells, k)
        KNN distances from spatial graph.
    min_expr_frac : float
        Minimum fraction of cells expressing a gene to count as testable
        (default: 0.05).
    min_cells : int
        Minimum number of cells for a cell type to be included (default: 50).
    verbose : bool
        Print summary tables (default: True).

    Returns
    -------
    summary : dict with keys:
        'cell_types' : DataFrame
            Per-cell-type info: n_cells, nn1_median, n_ligands, n_receptors,
            n_lr_within.
        'within' : DataFrame
            Within-type testable pairs per cell type.
        'between' : DataFrame
            Between-type testable pairs per (sender, receiver).
        'burden_within' : int
            Total within-type tests.
        'burden_between' : int
            Total between-type tests.
        'burden_total' : int
            Grand total testing burden.
        'expr_frac' : dict
            Expression fractions {cell_type: {gene: frac}}.
        'recommended_tau' : float
            Recommended tau based on median NN1 distance (cell diameter).
        'median_nn1' : float
            Global median nearest-neighbor distance.
        'cell_sizes' : dict
            Per-cell-type NN1 stats {cell_type: {nn1_median, nn1_q25, nn1_q75}}.

    Examples
    --------
    >>> import constellation as cst
    >>> summary = cst.scan_celltype_pairs(
    ...     adata, 'cell_type', lr_pairs, indices, distances)
    >>> # Use recommended tau
    >>> results = cst.run_celltype_analysis(
    ...     adata, 'cell_type', lr_pairs, indices, distances,
    ...     tau=summary['recommended_tau'],
    ...     expr_frac=summary['expr_frac'])
    >>> # Or override with custom tau
    >>> results = cst.run_celltype_analysis(
    ...     adata, 'cell_type', lr_pairs, indices, distances,
    ...     tau=10.0)
    """
    import pandas as pd

    cell_types = adata.obs[cell_type_col].values.astype(str)
    unique_types = sorted(set(cell_types))
    gene_set = set(adata.var_names)
    n_cells = len(cell_types)

    # LR genes in data
    all_ligands = set(lig for lig, _ in lr_pairs)
    all_receptors = set(rec for _, rec in lr_pairs)
    lr_in_data = [(l, r) for l, r in lr_pairs if l in gene_set and r in gene_set]

    # ── Cell size statistics ──
    nn1_dist = distances[:, 0]
    median_nn1_global = float(np.median(nn1_dist[nn1_dist > 0]))

    # Per-cell-type NN1 stats
    ct_nn1 = {}
    for ct in unique_types:
        ct_mask = cell_types == ct
        if ct_mask.sum() >= 10:
            nn1_ct = nn1_dist[ct_mask]
            ct_nn1[ct] = {
                'nn1_median': float(np.median(nn1_ct)),
                'nn1_q25': float(np.percentile(nn1_ct, 25)),
                'nn1_q75': float(np.percentile(nn1_ct, 75)),
            }

    # Recommended tau = median NN1 (approximates cell diameter)
    recommended_tau = round(median_nn1_global, 1)

    # Compute expression fractions
    expr_frac = compute_expression_fractions(adata, cell_type_col, min_expr_frac=0.0)

    # ── Per-cell-type stats ──
    ct_rows = []
    for ct in unique_types:
        n = int((cell_types == ct).sum())
        if n < min_cells:
            continue
        ct_frac = expr_frac[ct]
        # Count ligands/receptors expressed above threshold in this type
        lig_expressed = sum(1 for g in all_ligands & gene_set
                           if ct_frac.get(g, 0) >= min_expr_frac)
        rec_expressed = sum(1 for g in all_receptors & gene_set
                           if ct_frac.get(g, 0) >= min_expr_frac)
        # Within-type testable pairs
        n_within = sum(1 for l, r in lr_in_data
                       if ct_frac.get(l, 0) >= min_expr_frac
                       and ct_frac.get(r, 0) >= min_expr_frac)
        row = {
            'cell_type': ct,
            'n_cells': n,
            'n_ligands': lig_expressed,
            'n_receptors': rec_expressed,
            'n_lr_within': n_within,
        }
        if ct in ct_nn1:
            row['nn1_median'] = ct_nn1[ct]['nn1_median']
        ct_rows.append(row)

    ct_df = pd.DataFrame(ct_rows).sort_values('n_cells', ascending=False)
    valid_types = ct_df['cell_type'].tolist()

    # ── Within-type burden ──
    within_rows = []
    for _, row in ct_df.iterrows():
        if row['n_lr_within'] > 0:
            within_rows.append({
                'cell_type': row['cell_type'],
                'n_cells': row['n_cells'],
                'n_testable_lr': row['n_lr_within'],
            })
    within_df = pd.DataFrame(within_rows)
    burden_within = int(within_df['n_testable_lr'].sum()) if len(within_df) > 0 else 0

    # ── Between-type burden ──
    between_rows = []
    for sender in valid_types:
        s_frac = expr_frac[sender]
        for receiver in valid_types:
            if sender == receiver:
                continue
            r_frac = expr_frac[receiver]
            n_testable = sum(1 for l, r in lr_in_data
                             if s_frac.get(l, 0) >= min_expr_frac
                             and r_frac.get(r, 0) >= min_expr_frac)
            if n_testable > 0:
                between_rows.append({
                    'sender': sender,
                    'receiver': receiver,
                    'n_testable_lr': n_testable,
                })

    between_df = pd.DataFrame(between_rows)
    burden_between = int(between_df['n_testable_lr'].sum()) if len(between_df) > 0 else 0
    burden_total = burden_within + burden_between

    # ── Print summary ──
    if verbose:
        print(f'\n{"=" * 65}')
        print(f'CELL-TYPE PAIR SCAN')
        print(f'{"=" * 65}')
        print(f'LR database:    {len(lr_pairs)} pairs '
              f'({len(lr_in_data)} with both genes in data)')
        print(f'Cell types:     {len(valid_types)} '
              f'(>= {min_cells} cells, min_expr_frac={min_expr_frac})')
        print(f'')

        # Cell size section
        print(f'Cell spacing (NN1 distance as cell diameter proxy):')
        print(f'  Global median:  {median_nn1_global:.2f} um')
        nn1_vals = [ct_nn1[ct]['nn1_median'] for ct in valid_types if ct in ct_nn1]
        if nn1_vals:
            print(f'  Range by type:  {min(nn1_vals):.2f} - {max(nn1_vals):.2f} um')
        print(f'')

        # Tau recommendation
        print(f'Recommended tau:  {recommended_tau} um  '
              f'(= median NN1, approx. cell diameter)')
        print(f'  tau controls the spatial kernel: '
              f'weight = exp(-dist / tau)')
        print(f'  Effective interaction range ~ 3*tau = '
              f'{3*recommended_tau:.0f} um')
        print(f'  To test shorter-range interactions, decrease tau.')
        print(f'  To test longer-range interactions, increase tau.')
        print(f'')

        # Cell type table
        print(f'Per-cell-type summary:')
        print(f'{"Cell Type":<30s} {"Cells":>7s} {"NN1":>6s} '
              f'{"Lig":>5s} {"Rec":>5s} {"LR":>5s}')
        print(f'{"-"*30} {"-"*7} {"-"*6} {"-"*5} {"-"*5} {"-"*5}')
        for _, row in ct_df.iterrows():
            nn1_str = (f'{row["nn1_median"]:.1f}'
                       if 'nn1_median' in row and not pd.isna(row['nn1_median'])
                       else '  -')
            print(f'{row["cell_type"]:<30s} {row["n_cells"]:>7,} '
                  f'{nn1_str:>6s} '
                  f'{row["n_ligands"]:>5} {row["n_receptors"]:>5} '
                  f'{row["n_lr_within"]:>5}')

        print(f'')
        print(f'Testing burden:')
        print(f'  Within-type:  {burden_within:>8,} tests '
              f'({len(within_df)} cell types)')
        print(f'  Between-type: {burden_between:>8,} tests '
              f'({len(between_df)} directed pairs)')
        print(f'  Total:        {burden_total:>8,} tests')

        # Top between-type pairs
        if len(between_df) > 0:
            top = between_df.nlargest(10, 'n_testable_lr')
            print(f'\nTop between-type pairs:')
            for _, row in top.iterrows():
                print(f'  {row["sender"]:<25s} -> {row["receiver"]:<25s} '
                      f'{row["n_testable_lr"]:>5} LR pairs')

        print(f'{"=" * 65}\n')

    return {
        'cell_types': ct_df,
        'within': within_df,
        'between': between_df,
        'burden_within': burden_within,
        'burden_between': burden_between,
        'burden_total': burden_total,
        'expr_frac': expr_frac,
        'recommended_tau': recommended_tau,
        'median_nn1': median_nn1_global,
        'cell_sizes': ct_nn1,
    }


# =============================================================================
# HIGH-LEVEL ANALYSIS FUNCTION
# =============================================================================

def run_celltype_analysis(adata, cell_type_col, lr_pairs, indices, distances,
                 expr_frac=None, tau=5.0, min_expr_frac=0.05,
                 min_cells=50, fdr_method='fdr_bh', fdr_threshold=0.05,
                 test_within=True, test_between=True, X_csc=None,
                 lr_annotations=None, verbose=True):
    """
    Run full CONSTELLATION analysis: filter → test → FDR correct → summarize.

    This is the main entry point. For each testable (cell-type pair, LR pair)
    combination, it:
      1. Filters testable pairs by expression fraction thresholds
      2. Computes analytical z-scores (within-type and/or between-type)
      3. Applies genome-wide FDR correction (Benjamini-Hochberg)
      4. Prints a summary of significant results

    Parameters
    ----------
    adata : AnnData
        Annotated data with expression matrix and cell type annotations.
    cell_type_col : str
        Column in adata.obs containing cell type labels.
    lr_pairs : list of tuples
        List of (ligand, receptor) pairs to test.
    indices : array, shape (n_cells, k)
        KNN indices from spatial graph.
    distances : array, shape (n_cells, k)
        KNN distances from spatial graph.
    expr_frac : dict, optional
        Pre-computed expression fractions {cell_type: {gene: fraction}}.
        If None, computed automatically.
    tau : float
        Kernel decay distance in microns (default: 5.0).
        Weight = exp(-distance / tau). Effective range ≈ 3*tau.
    min_expr_frac : float
        Minimum fraction of cells expressing a gene to be testable (default: 0.05).
    min_cells : int
        Minimum number of cells for a cell type to be tested (default: 50).
    fdr_method : str
        Multiple testing correction method (default: 'fdr_bh').
        See statsmodels.stats.multitest.multipletests for options.
    fdr_threshold : float
        FDR significance threshold for summary (default: 0.05).
    test_within : bool
        Whether to test within-type (autocrine) interactions (default: True).
    test_between : bool
        Whether to test between-type (paracrine) interactions (default: True).
    X_csc : sparse matrix, optional
        Pre-computed CSC format of expression matrix for faster column access.
    lr_annotations : dict, optional
        LR pair annotations {(lig, rec): {'ligand_category': ..., ...}}.
    verbose : bool
        Print progress information (default: True).

    Returns
    -------
    results_df : DataFrame
        All test results with columns: ligand, receptor, sender, receiver,
        interaction_type, z_score, p_value, p_adj, fold_enrichment, n_edges,
        S_obs, null_mean, null_std, and any annotation columns.
        Sorted by p_adj ascending.

    Examples
    --------
    >>> import constellation as cst
    >>> from scipy.spatial import cKDTree
    >>>
    >>> # Load data and LR pairs
    >>> adata = sc.read_h5ad('spatial_data.h5ad')
    >>> lr_db = pd.read_csv('lr_pairs.csv')
    >>> lr_pairs = list(zip(lr_db['ligand'], lr_db['receptor']))
    >>>
    >>> # Build spatial graph
    >>> tree = cKDTree(adata.obsm['spatial'])
    >>> distances, indices = tree.query(adata.obsm['spatial'], k=11)
    >>> indices, distances = indices[:, 1:], distances[:, 1:]
    >>>
    >>> # Run full analysis (filter → test → FDR → summarize)
    >>> results = cst.run_celltype_analysis(
    ...     adata, cell_type_col='cell_type', lr_pairs=lr_pairs,
    ...     indices=indices, distances=distances,
    ...     tau=5.0
    ... )
    >>>
    >>> # Get significant results
    >>> sig = results[(results['p_adj'] < 0.05) & (results['z_score'] > 0)]
    """
    from statsmodels.stats.multitest import multipletests
    import pandas as pd

    # ── Step 0: Validate inputs ──
    if verbose:
        validate_inputs(adata, cell_type_col, lr_pairs, indices, distances,
                        tau=tau, verbose=True)

    cell_types = adata.obs[cell_type_col].values.astype(str)
    gene_names = adata.var_names.tolist()
    X = adata.X
    unique_types = sorted(set(cell_types))

    # Pre-compute CSC for sparse matrices
    if X_csc is None and issparse(X):
        if verbose:
            print('Converting to CSC format for fast column access...')
        X_csc = X.tocsc()

    # ── Step 1: Compute expression fractions ──
    if expr_frac is None:
        if verbose:
            print('Step 1/4: Computing expression fractions...')
        expr_frac = compute_expression_fractions(adata, cell_type_col, min_expr_frac=0.0)
    elif verbose:
        print('Step 1/4: Using pre-computed expression fractions.')

    # Filter cell types by minimum cell count
    ct_counts = {ct: (cell_types == ct).sum() for ct in unique_types}
    valid_types = [ct for ct in unique_types if ct_counts[ct] >= min_cells]
    if verbose:
        print(f'  Cell types: {len(valid_types)}/{len(unique_types)} '
              f'(>= {min_cells} cells)')

    # ── Step 2: Filter testable pairs ──
    if verbose:
        print(f'\nStep 2/4: Filtering testable combinations '
              f'(min_expr_frac={min_expr_frac})...')

    n_within_tests = 0
    n_between_tests = 0
    within_plan = {}   # {cell_type: [lr_pairs]}
    between_plan = {}  # {(sender, receiver): [lr_pairs]}

    if test_within:
        for ct in valid_types:
            testable = filter_testable_lr_pairs(
                lr_pairs, expr_frac, ct, min_expr_frac=min_expr_frac)
            if testable:
                within_plan[ct] = testable
                n_within_tests += len(testable)

    if test_between:
        for sender in valid_types:
            for receiver in valid_types:
                if sender == receiver:
                    continue
                testable = filter_testable_lr_pairs_between(
                    lr_pairs, expr_frac, sender, receiver,
                    min_expr_frac=min_expr_frac)
                if testable:
                    between_plan[(sender, receiver)] = testable
                    n_between_tests += len(testable)

    total_tests = n_within_tests + n_between_tests
    if verbose:
        print(f'  Within-type: {n_within_tests:,} tests '
              f'({len(within_plan)} cell types)')
        print(f'  Between-type: {n_between_tests:,} tests '
              f'({len(between_plan)} cell-type pairs)')
        print(f'  Total testable: {total_tests:,}')

    if total_tests == 0:
        if verbose:
            print('No testable combinations found.')
        return pd.DataFrame()

    # ── Step 3: Run analytical tests ──
    if verbose:
        print(f'\nStep 3/4: Computing analytical z-scores (tau={tau})...')

    all_results = []

    # Within-type
    if within_plan:
        if verbose:
            print(f'\n  Within-type ({n_within_tests:,} tests)...')
        iterator = tqdm(within_plan.items(), desc='  Within-type',
                        total=len(within_plan)) if verbose else within_plan.items()
        for ct, testable in iterator:
            res_list = test_within_type_lr(
                cell_type=ct, cell_types=cell_types, X=X,
                gene_names=gene_names, indices=indices, distances=distances,
                lr_pairs_subset=testable, lr_annotations=lr_annotations,
                tau=tau, X_csc=X_csc,
            )
            if res_list:
                for r in res_list:
                    r['sender'] = ct
                    r['receiver'] = ct
                    r['interaction_type'] = 'within'
                all_results.extend(res_list)

    # Between-type
    if between_plan:
        if verbose:
            print(f'\n  Between-type ({n_between_tests:,} tests)...')
        iterator = tqdm(between_plan.items(), desc='  Between-type',
                        total=len(between_plan)) if verbose else between_plan.items()
        for (sender, receiver), testable in iterator:
            res_list = test_between_type_lr(
                sender_type=sender, receiver_type=receiver,
                cell_types=cell_types, X=X, gene_names=gene_names,
                indices=indices, distances=distances,
                lr_pairs_subset=testable, lr_annotations=lr_annotations,
                tau=tau, X_csc=X_csc,
            )
            if res_list:
                for r in res_list:
                    r['interaction_type'] = 'between'
                all_results.extend(res_list)

    if not all_results:
        if verbose:
            print('No results generated.')
        return pd.DataFrame()

    # ── Step 4: FDR correction and summary ──
    results_df = pd.DataFrame(all_results)

    if verbose:
        print(f'\nStep 4/4: FDR correction ({fdr_method})...')
        print(f'  Tests performed: {len(results_df):,}')

    _, pvals_adj, _, _ = multipletests(results_df['p_value'], method=fdr_method)
    results_df['p_adj'] = pvals_adj
    results_df = results_df.sort_values('p_adj')

    # ── Summary ──
    sig_mask = (results_df['p_adj'] < fdr_threshold) & (results_df['z_score'] > 0)
    n_sig = sig_mask.sum()

    if verbose:
        print(f'\n{"=" * 50}')
        print(f'RESULTS SUMMARY')
        print(f'{"=" * 50}')
        print(f'Total tests:  {len(results_df):,}')
        print(f'Significant:  {n_sig:,} '
              f'(p_adj < {fdr_threshold}, z > 0)')

        if n_sig > 0:
            sig_df = results_df[sig_mask]
            n_within_sig = (sig_df['interaction_type'] == 'within').sum()
            n_between_sig = (sig_df['interaction_type'] == 'between').sum()
            print(f'  Within-type:  {n_within_sig:,}')
            print(f'  Between-type: {n_between_sig:,}')

            # Top cell-type pairs
            pair_counts = (sig_df.groupby(['sender', 'receiver'])
                           .size().sort_values(ascending=False))
            print(f'\nTop cell-type pairs ({min(10, len(pair_counts))} shown):')
            for (s, r), n in pair_counts.head(10).items():
                mean_fe = sig_df[(sig_df['sender'] == s) &
                                 (sig_df['receiver'] == r)]['fold_enrichment'].mean()
                print(f'  {s} -> {r}: {n} pairs (mean FE={mean_fe:.2f})')

    return results_df


# =============================================================================
# LINEAGE-LEVEL ANALYSIS
# =============================================================================

def run_lineage_analysis(adata, cell_type_col, grouping, lr_pairs, indices, distances,
                         expr_frac=None, tau=5.0, min_expr_frac=0.05,
                         min_cells=50, fdr_method='fdr_bh', fdr_threshold=0.05,
                         test_within=True, test_between=True, X_csc=None,
                         lr_annotations=None, verbose=True):
    """
    Run CONSTELLATION analysis at the lineage level for specificity assessment.

    Remaps fine-grained cell types to lineage groups, then delegates
    to ``run_celltype_analysis()``. This enables testing whether an LR
    interaction is specific to a particular subtype or broadly detected
    across related subtypes within a lineage.

    Parameters
    ----------
    adata : AnnData
        Annotated data with expression matrix and cell type annotations.
    cell_type_col : str
        Column in ``adata.obs`` containing fine-grained cell type labels.
    grouping : dict
        Mapping of lineage name to list of fine cell types.
        Example: ``{"CD4_T": ["Naive_CD4_T", "Treg", "Tfh_GC"], ...}``
    lr_pairs : list of tuples
        List of (ligand, receptor) pairs to test.
    indices : ndarray
        Pre-computed KNN indices from ``build_spatial_graph()``.
    distances : ndarray
        Pre-computed KNN distances from ``build_spatial_graph()``.
    expr_frac : dict, optional
        Pre-computed expression fractions. If None, computed internally.
        Note: must correspond to lineage-level groups, not fine types.
    tau : float, default=5.0
        Distance decay parameter for the exponential kernel.
    min_expr_frac : float, default=0.05
        Minimum expression fraction to consider a gene expressed in a group.
    min_cells : int, default=50
        Minimum cells in a lineage group to include it.
    fdr_method : str, default='fdr_bh'
        FDR correction method.
    fdr_threshold : float, default=0.05
        FDR significance threshold.
    test_within : bool, default=True
        Test within-lineage interactions.
    test_between : bool, default=True
        Test between-lineage interactions.
    X_csc : sparse matrix, optional
        Pre-computed CSC matrix.
    lr_annotations : dict, optional
        Annotations for LR pairs (e.g., category labels).
    verbose : bool, default=True
        Print progress and summary.

    Returns
    -------
    pd.DataFrame
        Same format as ``run_celltype_analysis()``, with sender/receiver
        columns containing lineage names.
    """
    import warnings

    # Invert grouping: {fine_type: lineage_name}
    fine_to_lineage = {}
    for lineage, fine_types in grouping.items():
        for ft in fine_types:
            if ft in fine_to_lineage:
                warnings.warn(
                    f"'{ft}' appears in multiple lineage groups: "
                    f"'{fine_to_lineage[ft]}' and '{lineage}'. "
                    f"Using '{lineage}'.")
            fine_to_lineage[ft] = lineage

    # Map fine types to lineage
    lineage_col = '_lineage'
    fine_types_in_data = adata.obs[cell_type_col].values.astype(str)
    lineage_labels = pd.Series(fine_types_in_data).map(fine_to_lineage)

    n_unmapped = lineage_labels.isna().sum()
    if n_unmapped > 0:
        unmapped_types = set(fine_types_in_data[lineage_labels.isna().values])
        if verbose:
            print(f'Warning: {n_unmapped:,} cells ({len(unmapped_types)} types) '
                  f'not in any lineage group: {sorted(unmapped_types)}')
        lineage_labels = lineage_labels.fillna('Unmapped')

    adata.obs[lineage_col] = lineage_labels.values

    if verbose:
        lineage_counts = adata.obs[lineage_col].value_counts()
        print(f'\nLineage groups ({len(grouping)}):')
        for lin in sorted(grouping.keys()):
            count = lineage_counts.get(lin, 0)
            n_subtypes = len(grouping[lin])
            print(f'  {lin}: {count:,} cells ({n_subtypes} subtypes)')
        if n_unmapped > 0:
            print(f'  Unmapped: {n_unmapped:,} cells (excluded)')

    try:
        # expr_frac must be recomputed for lineage grouping
        results = run_celltype_analysis(
            adata, cell_type_col=lineage_col, lr_pairs=lr_pairs,
            indices=indices, distances=distances, expr_frac=expr_frac,
            tau=tau, min_expr_frac=min_expr_frac, min_cells=min_cells,
            fdr_method=fdr_method, fdr_threshold=fdr_threshold,
            test_within=test_within, test_between=test_between,
            X_csc=X_csc, lr_annotations=lr_annotations, verbose=verbose,
        )
    finally:
        # Always clean up the temporary column
        if lineage_col in adata.obs.columns:
            del adata.obs[lineage_col]

    return results


# =============================================================================
# COMPARTMENT-LEVEL ANALYSIS
# =============================================================================

def run_compartment_analysis(
    expr_dict,
    compartments,
    coords,
    lr_pairs,
    k=10,
    tau=5.0,
    indices=None,
    distances=None,
    min_cells=10,
    fdr_method='fdr_bh',
    fdr_threshold=0.05,
    verbose=True,
):
    """
    Run CONSTELLATION analysis at the compartment level.

    Tests ligand-receptor colocalization within spatial compartments
    (e.g., Tumor, Interface, Tissue) using the same distance-weighted
    kernel test as the cell-type level analysis.

    For each compartment, edges are restricted to pairs where both cells
    belong to the compartment. The test statistic is T = dot(L, w), where
    L = log1p(ligand expression) and w_i = sum_j K(d_ij) * R_j with
    K = exp(-d/tau) and R = log1p(receptor expression).

    Parameters
    ----------
    expr_dict : dict
        Dictionary mapping gene names to expression arrays.
        Each array should have length n_cells.
    compartments : array-like
        Compartment labels for each cell (length n_cells).
    coords : array-like
        Spatial coordinates, shape (n_cells, 2).
    lr_pairs : list of tuples
        List of (ligand, receptor) pairs to test.
    k : int, default=10
        Number of neighbors for spatial graph.
    tau : float, default=5.0
        Distance decay parameter for the exponential kernel (microns).
    indices : array, optional
        Pre-computed KNN indices. If None, will be computed from coords.
    distances : array, optional
        Pre-computed KNN distances. If None, will be computed from coords.
    min_cells : int, default=10
        Minimum number of cells in a compartment to run test.
    fdr_method : str, default='fdr_bh'
        Method for FDR correction (passed to statsmodels.multipletests).
    fdr_threshold : float, default=0.05
        FDR threshold for significance.
    verbose : bool, default=True
        Print progress and summary.

    Returns
    -------
    pd.DataFrame
        Results with columns: lr_pair, ligand, receptor, compartment,
        n_cells, n_edges, S_obs, null_mean, null_std, fold_enrichment,
        z_score, p_value, p_adj, significant.

    Examples
    --------
    >>> # Prepare expression dictionary
    >>> expr_dict = {gene: get_expr(adata, gene) for gene in lr_genes}
    >>> compartments = adata.obs['compartment'].values
    >>> coords = adata.obsm['spatial']
    >>> lr_pairs = [('CXCL9', 'CXCR3'), ('CD40LG', 'CD40')]
    >>>
    >>> results = cst.run_compartment_analysis(
    ...     expr_dict, compartments, coords, lr_pairs,
    ...     k=10, tau=5.0, verbose=True
    ... )
    """
    from statsmodels.stats.multitest import multipletests

    compartments = np.asarray(compartments)
    coords = np.asarray(coords)
    n_cells = len(compartments)
    unique_compartments = sorted(set(compartments))

    if verbose:
        print('=' * 60)
        print('CONSTELLATION Compartment Analysis')
        print('=' * 60)
        print(f'Cells: {n_cells:,}')
        print(f'Compartments: {unique_compartments}')
        print(f'LR pairs: {len(lr_pairs)}')
        print(f'tau: {tau}')

    # Build spatial graph if not provided
    if indices is None or distances is None:
        if verbose:
            print(f'\nBuilding spatial graph (k={k})...')
        from .io import build_spatial_graph_from_coords
        indices, distances = build_spatial_graph_from_coords(
            coords, k=k, verbose=verbose)

    # --- Pre-compute within-compartment edge structure ---
    # For each compartment: find edges where both endpoints are in the compartment,
    # compute kernel weights, and cache local index mappings.
    comp_edge_data = {}

    for comp in unique_compartments:
        comp_mask = compartments == comp
        comp_idx = np.where(comp_mask)[0]
        n_comp = len(comp_idx)

        if n_comp < min_cells:
            continue

        # Boolean lookup for compartment membership
        is_comp = np.zeros(n_cells, dtype=bool)
        is_comp[comp_idx] = True

        # Find within-compartment edges
        nbr_indices = indices[comp_idx, :]     # (n_comp, K)
        nbr_dists = distances[comp_idx, :]     # (n_comp, K)
        mask = is_comp[nbr_indices]
        row_pos, col_pos = np.where(mask)

        n_edges = len(row_pos)
        if n_edges < min_cells:
            continue

        edge_dists = nbr_dists[row_pos, col_pos]
        K = np.exp(-edge_dists / tau)

        # Local indices
        local_edges_i = row_pos
        edges_j_global = nbr_indices[row_pos, col_pos]
        local_edges_j = np.searchsorted(comp_idx, edges_j_global)

        comp_edge_data[comp] = {
            'comp_idx': comp_idx,
            'n_comp': n_comp,
            'n_edges': n_edges,
            'K': K,
            'local_edges_i': local_edges_i,
            'local_edges_j': local_edges_j,
        }

    if verbose:
        print(f'\nCompartments with enough cells/edges: {len(comp_edge_data)}/{len(unique_compartments)}')
        for comp, data in comp_edge_data.items():
            print(f'  {comp}: {data["n_comp"]:,} cells, {data["n_edges"]:,} edges')

    # --- Filter valid LR pairs ---
    valid_lr = [(lig, rec) for lig, rec in lr_pairs
                if lig in expr_dict and rec in expr_dict]

    if verbose:
        print(f'\nTesting {len(valid_lr)} LR pairs x {len(comp_edge_data)} compartments...')

    # --- Run tests ---
    results = []
    iterator = tqdm(valid_lr, desc='LR pairs') if verbose else valid_lr

    for ligand, receptor in iterator:
        lr_name = f'{ligand}-{receptor}'
        l_expr_all = np.asarray(expr_dict[ligand], dtype=np.float64)
        r_expr_all = np.asarray(expr_dict[receptor], dtype=np.float64)

        for comp, edata in comp_edge_data.items():
            comp_idx = edata['comp_idx']
            n_comp = edata['n_comp']
            n_edges = edata['n_edges']
            K = edata['K']
            local_edges_i = edata['local_edges_i']
            local_edges_j = edata['local_edges_j']

            # Extract compartment expression and apply log1p
            L_comp = np.log1p(l_expr_all[comp_idx])
            R_comp = np.log1p(r_expr_all[comp_idx])

            # Compute w: distance-weighted receptor signal
            KR = K * R_comp[local_edges_j]
            w = np.zeros(n_comp)
            np.add.at(w, local_edges_i, KR)

            # Analytical permutation test
            test_result = _analytical_test(L_comp, w, n_comp)

            results.append({
                'lr_pair': lr_name,
                'ligand': ligand,
                'receptor': receptor,
                'compartment': comp,
                'n_cells': n_comp,
                'n_edges': n_edges,
                **test_result,
            })

    if not results:
        if verbose:
            print('No results generated.')
        return pd.DataFrame()

    results_df = pd.DataFrame(results)

    # FDR correction
    _, p_adj, _, _ = multipletests(results_df['p_value'], method=fdr_method)
    results_df['p_adj'] = p_adj
    results_df['significant'] = (results_df['p_adj'] < fdr_threshold) & (results_df['z_score'] > 0)
    results_df = results_df.sort_values('p_adj')

    # Summary
    n_sig = results_df['significant'].sum()

    if verbose:
        print(f'\n{"=" * 60}')
        print('RESULTS SUMMARY')
        print(f'{"=" * 60}')
        print(f'Total tests: {len(results_df):,}')
        print(f'Significant: {n_sig:,} (p_adj < {fdr_threshold}, z > 0)')

        if n_sig > 0:
            print(f'\nBy compartment:')
            for comp in unique_compartments:
                n_comp_sig = results_df[
                    (results_df['compartment'] == comp) &
                    results_df['significant']
                ].shape[0]
                print(f'  {comp}: {n_comp_sig}')

    return results_df


def scan_compartments(
    expr_dict,
    compartments,
    lr_pairs,
    min_expr_frac=0.01,
    min_cells=10,
    verbose=True,
):
    """
    Scan compartments to report testable LR pairs and expression statistics.

    This function provides a summary of which LR pairs are testable in each
    compartment based on expression thresholds, before running the full analysis.

    Parameters
    ----------
    expr_dict : dict
        Dictionary mapping gene names to expression arrays.
        Each array should have length n_cells.
    compartments : array-like
        Compartment labels for each cell (length n_cells).
    lr_pairs : list of tuples
        List of (ligand, receptor) pairs to check.
    min_expr_frac : float, default=0.01
        Minimum fraction of cells expressing a gene to consider it testable.
    min_cells : int, default=10
        Minimum number of ligand-positive cells in a compartment to be testable.
    verbose : bool, default=True
        Print summary report.

    Returns
    -------
    dict
        Dictionary with keys:
        - 'summary': DataFrame with per-compartment testable counts
        - 'testable_lr': dict mapping compartment -> list of testable (lig, rec) pairs
        - 'expr_frac': dict mapping compartment -> gene -> expression fraction
        - 'total_testable': total number of testable (LR pair, compartment) combinations
        - 'compartment_counts': dict mapping compartment -> n_cells

    Examples
    --------
    >>> report = cst.scan_compartments(expr_dict, compartments, lr_pairs)
    >>> print(f"Total testable: {report['total_testable']}")
    """
    import pandas as pd

    compartments = np.asarray(compartments)
    unique_compartments = sorted(set(compartments))
    n_cells = len(compartments)

    # Count cells per compartment
    compartment_counts = {}
    for comp in unique_compartments:
        compartment_counts[comp] = (compartments == comp).sum()

    # Compute expression fractions per compartment
    expr_frac = {}
    for comp in unique_compartments:
        comp_mask = compartments == comp
        n_comp = comp_mask.sum()
        expr_frac[comp] = {}
        for gene, expr in expr_dict.items():
            frac = (expr[comp_mask] > 0).sum() / n_comp
            expr_frac[comp][gene] = frac

    # Find testable LR pairs per compartment
    testable_lr = {}
    summary_rows = []

    for comp in unique_compartments:
        comp_mask = compartments == comp
        n_comp = compartment_counts[comp]
        testable = []

        for lig, rec in lr_pairs:
            if lig not in expr_dict or rec not in expr_dict:
                continue

            l_frac = expr_frac[comp].get(lig, 0)
            r_frac = expr_frac[comp].get(rec, 0)

            # Check if both genes are expressed above threshold
            if l_frac < min_expr_frac or r_frac < min_expr_frac:
                continue

            # Check minimum ligand-positive cells
            l_pos = (expr_dict[lig][comp_mask] > 0).sum()
            if l_pos < min_cells:
                continue

            testable.append((lig, rec))

        testable_lr[comp] = testable

        # Compute gene coverage
        n_lig_expressed = sum(1 for g in set(l for l, r in lr_pairs if l in expr_dict)
                              if expr_frac[comp].get(g, 0) >= min_expr_frac)
        n_rec_expressed = sum(1 for g in set(r for l, r in lr_pairs if r in expr_dict)
                              if expr_frac[comp].get(g, 0) >= min_expr_frac)

        summary_rows.append({
            'compartment': comp,
            'n_cells': n_comp,
            'pct_cells': 100 * n_comp / n_cells,
            'n_ligands_expressed': n_lig_expressed,
            'n_receptors_expressed': n_rec_expressed,
            'n_testable_lr': len(testable),
        })

    summary_df = pd.DataFrame(summary_rows)
    total_testable = sum(len(v) for v in testable_lr.values())

    if verbose:
        print('=' * 70)
        print('CONSTELLATION Compartment Scan')
        print('=' * 70)
        print(f'Total cells: {n_cells:,}')
        print(f'Compartments: {len(unique_compartments)}')
        print(f'Input LR pairs: {len(lr_pairs)}')
        print(f'Genes in expr_dict: {len(expr_dict)}')
        print(f'\nExpression threshold: {min_expr_frac:.1%}')
        print(f'Min ligand+ cells: {min_cells}')

        print(f'\n{"Compartment":<20} {"Cells":>10} {"Pct":>8} {"Lig":>6} {"Rec":>6} {"Testable LR":>12}')
        print('-' * 70)
        for _, row in summary_df.iterrows():
            print(f'{row["compartment"]:<20} {row["n_cells"]:>10,} {row["pct_cells"]:>7.1f}% '
                  f'{row["n_ligands_expressed"]:>6} {row["n_receptors_expressed"]:>6} '
                  f'{row["n_testable_lr"]:>12}')
        print('-' * 70)
        print(f'{"TOTAL":<20} {n_cells:>10,} {100:>7.1f}% {"":>6} {"":>6} {total_testable:>12}')
        print('=' * 70)

    return {
        'summary': summary_df,
        'testable_lr': testable_lr,
        'expr_frac': expr_frac,
        'total_testable': total_testable,
        'compartment_counts': compartment_counts,
    }


def compute_distance_profile(
    expr_dict,
    ligand,
    receptor,
    distance_values,
    indices,
    bin_edges=None,
    bin_width=15,
    min_cells=100,
):
    """
    Compute LR colocalization as a function of distance from a boundary.

    Parameters
    ----------
    expr_dict : dict
        Dictionary mapping gene names to expression arrays.
    ligand : str
        Ligand gene name.
    receptor : str
        Receptor gene name.
    distance_values : array
        Distance of each cell from the boundary (e.g., tumor edge).
        Negative values = inside boundary, positive = outside.
    indices : array
        KNN indices from spatial graph.
    bin_edges : array, optional
        Custom bin edges. If None, uses bin_width to create bins.
    bin_width : float, default=15
        Width of distance bins (microns).
    min_cells : int, default=100
        Minimum cells per bin to compute statistics.

    Returns
    -------
    dict
        Dictionary with keys:
        - bin_centers: array of bin centers
        - coexpr_fold: fold enrichment of L+R+ coexpression per bin
        - l_fraction: fraction of L+ cells per bin
        - r_fraction: fraction of R+ cells per bin
        - n_cells: number of cells per bin
        - coexpr_frac: fraction of L+ cells that are also R+ per bin
        - n_l_pos: count of L+ cells per bin
        - spatial_fold: fraction of L+ cells with an R+ neighbor /
          expected under independence (1-(1-r_global)^k). Requires
          ``indices`` to be provided.
    """
    l_expr = expr_dict[ligand]
    r_expr = expr_dict[receptor]
    distance_values = np.asarray(distance_values)

    if bin_edges is None:
        d_min, d_max = distance_values.min(), distance_values.max()
        bin_edges = np.arange(d_min, d_max + bin_width, bin_width)

    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Precompute boolean positivity arrays
    l_pos_all = l_expr > 0
    r_pos_all = r_expr > 0

    # Global receptor-positive fraction for spatial expected value
    r_pos_global = r_pos_all.mean()

    coexpr_fold = []
    l_fraction = []
    r_fraction = []
    n_cells_list = []
    coexpr_frac_list = []
    n_l_pos_list = []
    spatial_fold_list = []

    for lo, hi in zip(bin_edges[:-1], bin_edges[1:]):
        mask = (distance_values >= lo) & (distance_values < hi)
        n = mask.sum()

        if n < min_cells:
            coexpr_fold.append(np.nan)
            l_fraction.append(np.nan)
            r_fraction.append(np.nan)
            n_cells_list.append(n)
            coexpr_frac_list.append(np.nan)
            n_l_pos_list.append(0)
            spatial_fold_list.append(np.nan)
            continue

        l_pos_bin = l_pos_all[mask]
        r_pos_bin = r_pos_all[mask]

        l_frac = l_pos_bin.mean()
        r_frac = r_pos_bin.mean()
        lr_frac = (l_pos_bin & r_pos_bin).mean()
        expected = l_frac * r_frac
        fold = lr_frac / expected if expected > 0 else np.nan

        coexpr_fold.append(fold)
        l_fraction.append(l_frac)
        r_fraction.append(r_frac)
        n_cells_list.append(n)

        # Coexpression fraction: L+R+ / L+
        n_l = int(l_pos_bin.sum())
        n_l_pos_list.append(n_l)
        if n_l > 0:
            coexpr_frac_list.append((l_pos_bin & r_pos_bin).sum() / n_l)
        else:
            coexpr_frac_list.append(np.nan)

        # Spatial fold: fraction of L+ cells with at least one R+ neighbor
        # divided by expected under independence
        if indices is not None and n_l > 0:
            bin_indices = np.where(mask)[0]
            l_pos_idx = bin_indices[l_pos_bin]
            k = indices.shape[1]
            # For each L+ cell, check if any of its k neighbors are R+
            neighbor_r_pos = r_pos_all[indices[l_pos_idx]]  # (n_l, k)
            obs_frac = neighbor_r_pos.any(axis=1).mean()
            # Expected: probability that at least one of k neighbors is R+
            expected_spatial = 1.0 - (1.0 - r_pos_global) ** k
            spatial_fold_list.append(obs_frac / expected_spatial if expected_spatial > 0 else np.nan)
        else:
            spatial_fold_list.append(np.nan)

    return {
        'bin_centers': bin_centers,
        'coexpr_fold': np.array(coexpr_fold),
        'l_fraction': np.array(l_fraction),
        'r_fraction': np.array(r_fraction),
        'n_cells': np.array(n_cells_list),
        'coexpr_frac': np.array(coexpr_frac_list),
        'n_l_pos': np.array(n_l_pos_list),
        'spatial_fold': np.array(spatial_fold_list),
    }
