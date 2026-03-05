"""
Spatial LR Visualization Functions

Functions for visualizing LR analysis results.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns


def plot_celltype_pair_heatmap(results_df, value_col='z_score', agg_func='mean',
                                sig_only=True, fdr_threshold=0.05,
                                cmap='RdBu_r', center=0, figsize=(12, 10),
                                title=None, save_path=None):
    """
    Plot heatmap of cell type pairs.

    Parameters
    ----------
    results_df : DataFrame
        Between-type results with 'sender', 'receiver', and value columns
    value_col : str
        Column to aggregate ('z_score', 'p_adj', etc.)
    agg_func : str
        Aggregation function ('mean', 'count', 'max')
    sig_only : bool
        Only include significant results
    fdr_threshold : float
        FDR threshold for significance
    cmap : str
        Colormap
    center : float or None
        Center value for diverging colormap
    figsize : tuple
        Figure size
    title : str, optional
        Plot title
    save_path : str, optional
        Path to save figure
    """
    df = results_df.copy()
    if sig_only and 'p_adj' in df.columns:
        df = df[df['p_adj'] < fdr_threshold]

    if len(df) == 0:
        print("No significant results to plot")
        return

    # Create pivot table
    if agg_func == 'count':
        pivot = df.groupby(['sender', 'receiver']).size().unstack(fill_value=0)
    else:
        pivot = df.groupby(['sender', 'receiver'])[value_col].agg(agg_func).unstack(fill_value=0)

    # Plot
    fig, ax = plt.subplots(figsize=figsize)

    sns.heatmap(pivot, ax=ax, cmap=cmap, center=center,
                annot=False, cbar_kws={'label': value_col})

    if title:
        ax.set_title(title, fontsize=12, fontweight='bold')
    ax.set_xlabel('Receiver Cell Type')
    ax.set_ylabel('Sender Cell Type')
    ax.tick_params(axis='x', rotation=90, labelsize=9)
    ax.tick_params(axis='y', rotation=0, labelsize=9)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved: {save_path}")

    return fig, ax


def plot_lr_category_heatmap(results_df, category_col='ligand_category_major',
                              sig_only=True, fdr_threshold=0.05,
                              figsize=(14, 12), save_path=None):
    """
    Plot heatmaps for different LR categories.

    Parameters
    ----------
    results_df : DataFrame
        Between-type results
    category_col : str
        Category column to split by
    sig_only : bool
        Only include significant results
    fdr_threshold : float
        FDR threshold
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure
    """
    df = results_df.copy()
    if sig_only and 'p_adj' in df.columns:
        df = df[df['p_adj'] < fdr_threshold]

    categories = df[category_col].value_counts().head(4).index.tolist()

    fig, axes = plt.subplots(2, 2, figsize=figsize)

    for idx, cat in enumerate(categories):
        ax = axes[idx // 2, idx % 2]
        cat_df = df[df[category_col] == cat]

        if len(cat_df) > 0:
            pivot = cat_df.groupby(['sender', 'receiver']).size().unstack(fill_value=0)

            sns.heatmap(pivot, ax=ax, cmap='YlOrRd', annot=False,
                        cbar_kws={'label': '# Significant'})
            ax.set_title(f'{cat} (n={len(cat_df)})', fontsize=11, fontweight='bold')
        else:
            ax.set_title(f'{cat} (n=0)')

        ax.set_xlabel('Receiver')
        ax.set_ylabel('Sender')
        ax.tick_params(axis='x', rotation=90, labelsize=7)
        ax.tick_params(axis='y', rotation=0, labelsize=7)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved: {save_path}")

    return fig, axes


def plot_cell_lineage_tree(lineage_dict, counts_dict, figsize=(16, 10), save_path=None):
    """
    Plot cell lineage hierarchy.

    Parameters
    ----------
    lineage_dict : dict
        Nested dict of lineage hierarchy
    counts_dict : dict
        Cell counts per group
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Colors by major lineage
    lineage_colors = {
        'T cell': '#E41A1C',
        'B cell': '#377EB8',
        'Myeloid': '#4DAF4A',
        'Other': '#984EA3',
    }

    y_pos = 0.9
    x_start = 0.05

    for major_name, major_info in lineage_dict.items():
        color = major_info.get('color', '#999999')
        children = major_info.get('children', {})
        n_children = len(children)

        major_width = 0.2

        # Draw major category box
        rect = patches.FancyBboxPatch(
            (x_start, y_pos - 0.08), major_width, 0.12,
            boxstyle="round,pad=0.01",
            facecolor=color, edgecolor='black', linewidth=2
        )
        ax.add_patch(rect)
        ax.text(x_start + major_width/2, y_pos - 0.02, major_name,
                ha='center', va='center', fontsize=11, fontweight='bold')

        # Draw children
        if n_children > 0:
            child_y = y_pos - 0.25
            child_width = major_width / n_children - 0.01

            for i, (child_name, child_color) in enumerate(children.items()):
                child_x = x_start + i * (major_width / n_children)

                # Draw line
                ax.plot([x_start + major_width/2, child_x + child_width/2],
                        [y_pos - 0.08, child_y + 0.06], 'k-', linewidth=1)

                # Draw child box
                rect = patches.FancyBboxPatch(
                    (child_x, child_y - 0.04), child_width, 0.1,
                    boxstyle="round,pad=0.01",
                    facecolor=child_color, edgecolor='black', linewidth=1
                )
                ax.add_patch(rect)
                ax.text(child_x + child_width/2, child_y + 0.01, child_name,
                        ha='center', va='center', fontsize=8)

        x_start += major_width + 0.05

    ax.set_xlim(0, 1)
    ax.set_ylim(0.5, 1)
    ax.axis('off')
    ax.set_title('Cell Lineage Hierarchy', fontsize=14, fontweight='bold')

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight', facecolor='white')
        print(f"Saved: {save_path}")

    return fig, ax


def plot_celltype_barplot(counts_dict, order=None, colors=None,
                          figsize=(12, 8), title=None, save_path=None):
    """
    Plot bar chart of cell type abundances.

    Parameters
    ----------
    counts_dict : dict
        {cell_type: count}
    order : list, optional
        Order of cell types
    colors : list, optional
        Colors for each cell type
    figsize : tuple
        Figure size
    title : str, optional
        Plot title
    save_path : str, optional
        Path to save figure
    """
    if order is None:
        order = list(counts_dict.keys())

    counts = [counts_dict[g] for g in order]

    fig, ax = plt.subplots(figsize=figsize)

    y_pos = np.arange(len(order))

    if colors is None:
        colors = ['#377EB8'] * len(order)

    bars = ax.barh(y_pos, counts, color=colors, edgecolor='black', linewidth=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(order)
    ax.invert_yaxis()
    ax.set_xlabel('Number of Cells', fontsize=12)

    if title:
        ax.set_title(title, fontsize=14, fontweight='bold')

    # Add count labels
    for bar, count in zip(bars, counts):
        ax.text(bar.get_width() + max(counts)*0.01, bar.get_y() + bar.get_height()/2,
                f'{count:,}', va='center', fontsize=9)

    ax.set_xlim(0, max(counts) * 1.15)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved: {save_path}")

    return fig, ax


def plot_lr_interaction_summary(results_df, group_col='cell_type',
                                 sig_only=True, fdr_threshold=0.05,
                                 top_n=20, figsize=(10, 8), save_path=None):
    """
    Plot summary of significant interactions per cell type.

    Parameters
    ----------
    results_df : DataFrame
        Results dataframe
    group_col : str
        Column to group by
    sig_only : bool
        Only count significant
    fdr_threshold : float
        FDR threshold
    top_n : int
        Number of top groups to show
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure
    """
    df = results_df.copy()
    if sig_only and 'p_adj' in df.columns:
        df = df[df['p_adj'] < fdr_threshold]

    counts = df.groupby(group_col).size().sort_values(ascending=False).head(top_n)

    fig, ax = plt.subplots(figsize=figsize)

    bars = ax.barh(range(len(counts)), counts.values, color='steelblue', edgecolor='black')
    ax.set_yticks(range(len(counts)))
    ax.set_yticklabels(counts.index)
    ax.invert_yaxis()
    ax.set_xlabel('Number of Significant Interactions')
    ax.set_title(f'Significant LR Interactions by {group_col}', fontweight='bold')

    for bar, count in zip(bars, counts.values):
        ax.text(bar.get_width() + max(counts.values)*0.01,
                bar.get_y() + bar.get_height()/2,
                f'{count}', va='center', fontsize=9)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved: {save_path}")

    return fig, ax


def create_dotplot_markers(adata, markers_dict, groupby, categories_order=None,
                           standard_scale='var', cmap='Reds',
                           swap_axes=False, save_path=None, title=None):
    """
    Create scanpy dotplot for canonical markers.

    Parameters
    ----------
    adata : AnnData
        Annotated data object
    markers_dict : dict
        {category: [genes]} for markers
    groupby : str
        Column to group by
    categories_order : list, optional
        Order of categories
    standard_scale : str
        Scaling method ('var' or 'group')
    cmap : str
        Colormap
    swap_axes : bool
        If True, swap x and y axes
    save_path : str, optional
        Path suffix for saving
    title : str, optional
        Plot title
    """
    import scanpy as sc

    # Filter to genes in panel
    filtered_markers = {}
    for cat, genes in markers_dict.items():
        genes_in_panel = [g for g in genes if g in adata.var_names]
        if genes_in_panel:
            filtered_markers[cat] = genes_in_panel

    sc.pl.dotplot(
        adata,
        var_names=filtered_markers,
        groupby=groupby,
        categories_order=categories_order,
        standard_scale=standard_scale,
        cmap=cmap,
        swap_axes=swap_axes,
        save=save_path,
        show=False,
        title=title
    )

    if save_path:
        print(f"Saved dotplot")


# =============================================================================
# CELL-TYPE PAIR-WISE INTERACTION VISUALIZATIONS
# =============================================================================

def _prepare_interaction_data(results_df, sig_only=True, fdr_threshold=0.05,
                               min_count=1):
    """
    Internal helper to filter and aggregate results for visualization.

    Parameters
    ----------
    results_df : DataFrame
        Combined results from run_celltype_analysis().
    sig_only : bool
        Filter to significant results (p_adj < fdr_threshold and z_score > 0).
    fdr_threshold : float
        Significance threshold.
    min_count : int
        Minimum significant LR pairs to include a cell-type pair.

    Returns
    -------
    grouped : DataFrame
        Aggregated per (sender, receiver) with columns:
        sender, receiver, count, mean_z, max_z, mean_fe, is_within.
    all_types : list
        Sorted list of all cell types.
    """
    df = results_df.copy()
    if sig_only and 'p_adj' in df.columns:
        df = df[(df['p_adj'] < fdr_threshold) & (df['z_score'] > 0)]

    if len(df) == 0:
        return pd.DataFrame(), []

    grouped = df.groupby(['sender', 'receiver']).agg(
        count=('z_score', 'size'),
        mean_z=('z_score', 'mean'),
        max_z=('z_score', 'max'),
        mean_fe=('fold_enrichment', 'mean'),
    ).reset_index()

    grouped['is_within'] = grouped['sender'] == grouped['receiver']
    grouped = grouped[grouped['count'] >= min_count]

    all_types = sorted(set(grouped['sender']) | set(grouped['receiver']))

    return grouped, all_types


def plot_combined_heatmap(results_df, value_col='z_score', agg_func='count',
                           sig_only=True, fdr_threshold=0.05,
                           cmap='YlOrBr', center=None, annot=True,
                           cell_type_order=None, cluster=False,
                           diag_linewidth=2.0, figsize=None,
                           title=None, save_path=None):
    """
    Plot combined heatmap with within-type on diagonal and between-type off-diagonal.

    Parameters
    ----------
    results_df : DataFrame
        Combined results from run_celltype_analysis() containing both
        'within' and 'between' interaction_type rows.
    value_col : str
        Column to aggregate ('z_score', 'fold_enrichment'). Ignored when
        agg_func='count'.
    agg_func : str
        Aggregation function: 'mean', 'max', 'count'.
    sig_only : bool
        Filter to significant results (p_adj < fdr_threshold, z_score > 0).
    fdr_threshold : float
        Significance threshold.
    cmap : str
        Colormap for the heatmap.
    center : float or None
        Center value for diverging colormaps. Use 0 with 'RdBu_r'.
    annot : bool
        Annotate cells with values.
    cell_type_order : list or None
        Order of cell types. If None, alphabetical or clustered.
    cluster : bool
        If True and cell_type_order is None, order by hierarchical clustering.
    diag_linewidth : float
        Border thickness for diagonal (within-type) cells.
    figsize : tuple or None
        Figure size. Auto-computed if None.
    title : str or None
        Plot title.
    save_path : str or None
        Path to save figure.

    Returns
    -------
    fig, ax
    """
    df = results_df.copy()
    if sig_only and 'p_adj' in df.columns:
        df = df[(df['p_adj'] < fdr_threshold) & (df['z_score'] > 0)]

    if len(df) == 0:
        print("No significant results to plot")
        return None, None

    # Separate within and between
    within_df = df[df['interaction_type'] == 'within'] if 'interaction_type' in df.columns else df[df['sender'] == df['receiver']]
    between_df = df[df['interaction_type'] == 'between'] if 'interaction_type' in df.columns else df[df['sender'] != df['receiver']]

    # Get all cell types
    all_types = sorted(set(df['sender']) | set(df['receiver']))
    if cell_type_order is not None:
        all_types = [ct for ct in cell_type_order if ct in set(df['sender']) | set(df['receiver'])]

    n = len(all_types)
    type_to_idx = {ct: i for i, ct in enumerate(all_types)}

    # Build matrix
    matrix = np.full((n, n), np.nan)

    # Aggregate function
    def _agg(sub_df):
        if agg_func == 'count':
            return sub_df.groupby(['sender', 'receiver']).size()
        else:
            return sub_df.groupby(['sender', 'receiver'])[value_col].agg(agg_func)

    # Fill off-diagonal
    if len(between_df) > 0:
        agg_between = _agg(between_df)
        for (s, r), val in agg_between.items():
            if s in type_to_idx and r in type_to_idx:
                matrix[type_to_idx[s], type_to_idx[r]] = val

    # Fill diagonal
    if len(within_df) > 0:
        agg_within = _agg(within_df)
        for (s, r), val in agg_within.items():
            if s in type_to_idx:
                matrix[type_to_idx[s], type_to_idx[s]] = val

    # Hierarchical clustering
    if cluster and cell_type_order is None and n > 2:
        from scipy.cluster.hierarchy import linkage, leaves_list
        from scipy.spatial.distance import squareform
        matrix_fill = np.nan_to_num(matrix, nan=0.0)
        combined = matrix_fill + matrix_fill.T
        Z = linkage(combined, method='average')
        order = leaves_list(Z)
        all_types = [all_types[i] for i in order]
        matrix = matrix[np.ix_(order, order)]

    # Auto figsize
    if figsize is None:
        side = max(8, n * 0.55 + 3)
        figsize = (side, side * 0.9)

    # Plot
    fig, ax = plt.subplots(figsize=figsize)

    # Mask NaN for gray background
    masked = np.ma.masked_invalid(matrix)
    cmap_obj = plt.cm.get_cmap(cmap).copy()
    cmap_obj.set_bad('#E8E8E8')

    # Determine vmin/vmax
    valid = matrix[~np.isnan(matrix)]
    if len(valid) == 0:
        print("No data to plot")
        return None, None

    if center is not None:
        vmax = max(abs(valid.max() - center), abs(valid.min() - center))
        vmin_plot = center - vmax
        vmax_plot = center + vmax
    else:
        vmin_plot = valid.min()
        vmax_plot = valid.max()

    im = ax.imshow(masked, cmap=cmap_obj, aspect='equal',
                   vmin=vmin_plot, vmax=vmax_plot)

    # Colorbar
    label_map = {'count': '# Significant LR pairs',
                 'mean': f'Mean {value_col}',
                 'max': f'Max {value_col}'}
    cbar_label = label_map.get(agg_func, f'{agg_func}({value_col})')
    cbar = fig.colorbar(im, ax=ax, shrink=0.8, label=cbar_label)

    # Annotations
    if annot:
        fmt = 'd' if agg_func == 'count' else '.1f'
        for i in range(n):
            for j in range(n):
                val = matrix[i, j]
                if np.isnan(val):
                    continue
                # Determine text color (white on dark, black on light)
                norm_val = (val - vmin_plot) / (vmax_plot - vmin_plot) if vmax_plot != vmin_plot else 0.5
                text_color = 'white' if norm_val > 0.7 else 'black'
                text = f'{int(val)}' if agg_func == 'count' else f'{val:.1f}'
                fontweight = 'bold' if i == j else 'normal'
                fontsize = 8 if n > 20 else (9 if n > 10 else 10)
                ax.text(j, i, text, ha='center', va='center',
                        color=text_color, fontsize=fontsize,
                        fontweight=fontweight)

    # Diagonal borders
    for i in range(n):
        rect = patches.Rectangle(
            (i - 0.5, i - 0.5), 1, 1,
            fill=False, edgecolor='black',
            linewidth=diag_linewidth, clip_on=True, zorder=3
        )
        ax.add_patch(rect)

    # Labels
    ax.set_xticks(range(n))
    ax.set_xticklabels(all_types, rotation=45, ha='right', fontsize=10)
    ax.set_yticks(range(n))
    ax.set_yticklabels(all_types, fontsize=10)
    ax.set_xlabel('Receiver Cell Type', fontsize=11)
    ax.set_ylabel('Sender Cell Type', fontsize=11)

    if title is None:
        title = f'Cell-Type Interactions ({cbar_label})'
    ax.set_title(title, fontsize=12, fontweight='bold', pad=10)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved: {save_path}")

    return fig, ax


def plot_interaction_dotplot(results_df, sig_only=True, fdr_threshold=0.05,
                              color_col='z_score', color_agg='mean',
                              min_count=1, max_dot_size=400, min_dot_size=20,
                              cmap='YlOrBr', edgecolor='black', linewidth=0.8,
                              cell_type_order=None, cluster=False,
                              show_within=True, annotate_threshold=None,
                              figsize=None, title=None, save_path=None):
    """
    Plot dot plot of cell-type interactions.

    Dot size = number of significant LR pairs.
    Dot color = mean z-score (or other metric).

    Parameters
    ----------
    results_df : DataFrame
        Combined results from run_celltype_analysis().
    sig_only : bool
        Filter to significant results.
    fdr_threshold : float
        Significance threshold.
    color_col : str
        Column for dot color: 'z_score', 'fold_enrichment'.
    color_agg : str
        Aggregation for color: 'mean', 'max'.
    min_count : int
        Minimum significant pairs to show a dot.
    max_dot_size, min_dot_size : float
        Dot size range in points^2.
    cmap : str
        Colormap for dot color.
    edgecolor : str
        Dot edge color.
    linewidth : float
        Dot edge width.
    cell_type_order : list or None
        Order of cell types.
    cluster : bool
        Cluster cell types if cell_type_order is None.
    show_within : bool
        Show within-type results on diagonal (diamond markers).
    annotate_threshold : int or None
        Annotate dots with count if count >= threshold.
    figsize : tuple or None
        Figure size.
    title : str or None
        Plot title.
    save_path : str or None
        Path to save figure.

    Returns
    -------
    fig, ax
    """
    grouped, all_types = _prepare_interaction_data(
        results_df, sig_only=sig_only, fdr_threshold=fdr_threshold,
        min_count=min_count)

    if len(grouped) == 0:
        print("No significant results to plot")
        return None, None

    # Cell type ordering
    if cell_type_order is not None:
        all_types = [ct for ct in cell_type_order if ct in all_types]
    elif cluster and len(all_types) > 2:
        from scipy.cluster.hierarchy import linkage, leaves_list
        n = len(all_types)
        matrix = np.zeros((n, n))
        type_to_idx = {ct: i for i, ct in enumerate(all_types)}
        for _, row in grouped.iterrows():
            si = type_to_idx.get(row['sender'])
            ri = type_to_idx.get(row['receiver'])
            if si is not None and ri is not None:
                matrix[si, ri] = row['count']
        combined = matrix + matrix.T
        Z = linkage(combined, method='average')
        order = leaves_list(Z)
        all_types = [all_types[i] for i in order]

    n_types = len(all_types)
    type_to_idx = {ct: i for i, ct in enumerate(all_types)}

    # Color aggregation column
    color_map = {'mean': 'mean_z', 'max': 'max_z'}
    if color_col == 'fold_enrichment':
        color_map = {'mean': 'mean_fe', 'max': 'mean_fe'}
    color_key = color_map.get(color_agg, 'mean_z')

    # Auto figsize
    if figsize is None:
        side = max(8, n_types * 0.55 + 3)
        figsize = (side + 2, side)

    fig, ax = plt.subplots(figsize=figsize)

    # Background grid
    for i in range(n_types):
        ax.axhline(i, color='#E8E8E8', linewidth=0.5, zorder=0)
        ax.axvline(i, color='#E8E8E8', linewidth=0.5, zorder=0)

    ax.set_facecolor('#FAFAF8')

    # Separate within and between
    between = grouped[~grouped['is_within']]
    within = grouped[grouped['is_within']]

    # Size scaling
    all_counts = grouped['count'].values
    count_min, count_max = all_counts.min(), all_counts.max()

    def _scale_size(counts):
        if count_max > count_min:
            return min_dot_size + (counts - count_min) / (count_max - count_min) * (max_dot_size - min_dot_size)
        return np.full(len(counts), (max_dot_size + min_dot_size) / 2)

    # Color range
    all_colors = grouped[color_key].values
    vmin_c, vmax_c = all_colors.min(), all_colors.max()

    # Plot between-type (circles)
    if len(between) > 0:
        x_bt = [type_to_idx[r] for r in between['receiver']]
        y_bt = [type_to_idx[s] for s in between['sender']]
        sizes_bt = _scale_size(between['count'].values)
        colors_bt = between[color_key].values

        sc_bt = ax.scatter(x_bt, y_bt, s=sizes_bt, c=colors_bt,
                           cmap=cmap, vmin=vmin_c, vmax=vmax_c,
                           edgecolors=edgecolor, linewidths=linewidth,
                           marker='o', zorder=3)

    # Plot within-type (diamonds)
    if show_within and len(within) > 0:
        x_wt = [type_to_idx[r] for r in within['receiver']]
        y_wt = [type_to_idx[s] for s in within['sender']]
        sizes_wt = _scale_size(within['count'].values)
        colors_wt = within[color_key].values

        sc_wt = ax.scatter(x_wt, y_wt, s=sizes_wt, c=colors_wt,
                           cmap=cmap, vmin=vmin_c, vmax=vmax_c,
                           edgecolors=edgecolor, linewidths=linewidth,
                           marker='D', zorder=4)

    # Annotations
    if annotate_threshold is not None:
        for _, row in grouped.iterrows():
            if row['count'] >= annotate_threshold:
                xi = type_to_idx.get(row['receiver'])
                yi = type_to_idx.get(row['sender'])
                if xi is not None and yi is not None:
                    ax.text(xi, yi, f"{int(row['count'])}", ha='center',
                            va='center', fontsize=7, fontweight='bold',
                            color='white', zorder=5)

    # Colorbar
    color_label = f'{"Mean" if color_agg == "mean" else "Max"} {color_col}'
    ref_scatter = sc_bt if len(between) > 0 else sc_wt
    cbar = fig.colorbar(ref_scatter, ax=ax, shrink=0.7, label=color_label)

    # Size legend
    if count_max > count_min:
        legend_counts = np.unique(np.round(
            np.linspace(count_min, count_max, 4)).astype(int))
    else:
        legend_counts = [int(count_min)]
    legend_sizes = _scale_size(np.array(legend_counts, dtype=float))

    legend_elements = []
    for cnt, sz in zip(legend_counts, legend_sizes):
        legend_elements.append(
            plt.scatter([], [], s=sz, c='gray', edgecolors='black',
                        linewidths=0.5, marker='o', label=f'{int(cnt)}')
        )

    # Marker shape legend
    if show_within and len(within) > 0:
        legend_elements.append(
            plt.scatter([], [], s=80, c='gray', edgecolors='black',
                        linewidths=0.5, marker='o', label='Between-type'))
        legend_elements.append(
            plt.scatter([], [], s=80, c='gray', edgecolors='black',
                        linewidths=0.5, marker='D', label='Within-type'))

    leg = ax.legend(handles=legend_elements, title='# Sig. pairs',
                    loc='upper left', bbox_to_anchor=(1.15, 1.0),
                    framealpha=0.9, edgecolor='gray', fontsize=9,
                    title_fontsize=9)

    # Axis labels
    ax.set_xticks(range(n_types))
    ax.set_xticklabels(all_types, rotation=45, ha='right', fontsize=10)
    ax.set_yticks(range(n_types))
    ax.set_yticklabels(all_types, fontsize=10)
    ax.set_xlabel('Receiver Cell Type', fontsize=11)
    ax.set_ylabel('Sender Cell Type', fontsize=11)
    ax.set_xlim(-0.5, n_types - 0.5)
    ax.set_ylim(n_types - 0.5, -0.5)

    if title is None:
        title = 'Cell-Type Interaction Dot Plot'
    ax.set_title(title, fontsize=12, fontweight='bold', pad=10)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved: {save_path}")

    return fig, ax


def plot_interaction_network(results_df, sig_only=True, fdr_threshold=0.05,
                              edge_metric='count', min_count=1,
                              node_size_metric='total_interactions',
                              node_colors=None, default_node_color='#E8B87D',
                              edge_cmap='YlOrBr', edge_width_range=(0.5, 5.0),
                              edge_alpha=0.6, show_self_loops=True,
                              self_loop_scale=0.12,
                              cell_type_order=None,
                              node_label_fontsize=10,
                              figsize=(10, 10), title=None, save_path=None):
    """
    Plot circular network diagram of cell-type interactions.

    Pure matplotlib implementation (no networkx). Nodes placed on a circle,
    directed edges drawn as curved arrows.

    Parameters
    ----------
    results_df : DataFrame
        Combined results from run_celltype_analysis().
    sig_only : bool
        Filter to significant results.
    fdr_threshold : float
        Significance threshold.
    edge_metric : str
        What edge width encodes: 'count' or 'mean_z'.
    min_count : int
        Minimum significant pairs to draw an edge.
    node_size_metric : str
        Node sizing: 'total_interactions' or 'uniform'.
    node_colors : dict or None
        {cell_type: color} for custom node colors.
    default_node_color : str
        Default node fill color.
    edge_cmap : str
        Colormap for edge colors.
    edge_width_range : tuple
        (min_width, max_width) for edge line widths.
    edge_alpha : float
        Edge transparency.
    show_self_loops : bool
        Show within-type interactions as self-loops.
    self_loop_scale : float
        Radius of self-loop arcs relative to plot radius.
    cell_type_order : list or None
        Order of nodes around the circle.
    node_label_fontsize : int
        Font size for node labels.
    figsize : tuple
        Figure size.
    title : str or None
        Plot title.
    save_path : str or None
        Path to save figure.

    Returns
    -------
    fig, ax
    """
    from matplotlib.patches import FancyArrowPatch, Arc

    grouped, all_types = _prepare_interaction_data(
        results_df, sig_only=sig_only, fdr_threshold=fdr_threshold,
        min_count=min_count)

    if len(grouped) == 0:
        print("No significant results to plot")
        return None, None

    # Cell type ordering
    if cell_type_order is not None:
        all_types = [ct for ct in cell_type_order if ct in all_types]

    n = len(all_types)

    # Circular layout — start from top, go clockwise
    radius = 1.0
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    angles = np.pi / 2 - angles  # start from top

    node_positions = {}
    node_angles = {}
    for i, ct in enumerate(all_types):
        node_positions[ct] = (radius * np.cos(angles[i]),
                              radius * np.sin(angles[i]))
        node_angles[ct] = angles[i]

    # Node sizes
    if node_size_metric == 'total_interactions':
        sender_counts = grouped.groupby('sender')['count'].sum()
        receiver_counts = grouped.groupby('receiver')['count'].sum()
        total = sender_counts.add(receiver_counts, fill_value=0)
        min_node, max_node = 0.06, 0.14
        t_min, t_max = total.min(), total.max()
        node_radii = {}
        for ct in all_types:
            t = total.get(ct, 0)
            if t_max > t_min:
                node_radii[ct] = min_node + (t - t_min) / (t_max - t_min) * (max_node - min_node)
            else:
                node_radii[ct] = (min_node + max_node) / 2
    else:
        node_radii = {ct: 0.10 for ct in all_types}

    # Edge metrics
    between_edges = grouped[~grouped['is_within']].copy()
    within_edges = grouped[grouped['is_within']].copy()

    metric_col = 'count' if edge_metric == 'count' else 'mean_z'
    min_w, max_w = edge_width_range

    if len(between_edges) > 0:
        e_vals = between_edges[metric_col].values
        if len(within_edges) > 0:
            all_e_vals = np.concatenate([e_vals, within_edges[metric_col].values])
        else:
            all_e_vals = e_vals
    elif len(within_edges) > 0:
        all_e_vals = within_edges[metric_col].values
    else:
        all_e_vals = np.array([1.0])

    e_min, e_max = all_e_vals.min(), all_e_vals.max()

    def _edge_width(val):
        if e_max > e_min:
            return min_w + (val - e_min) / (e_max - e_min) * (max_w - min_w)
        return (min_w + max_w) / 2

    # Edge colormap
    cmap_obj = plt.cm.get_cmap(edge_cmap)

    def _edge_color(val):
        if e_max > e_min:
            return cmap_obj((val - e_min) / (e_max - e_min))
        return cmap_obj(0.5)

    # Plot
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_aspect('equal')
    ax.axis('off')

    # Track bidirectional edges for curvature direction
    seen_pairs = set()

    # Draw between-type edges
    for _, row in between_edges.iterrows():
        s_pos = np.array(node_positions[row['sender']])
        r_pos = np.array(node_positions[row['receiver']])
        s_rad = node_radii[row['sender']]
        r_rad = node_radii[row['receiver']]

        # Shorten arrow to node edges
        dx = r_pos - s_pos
        dist = np.linalg.norm(dx)
        if dist < 1e-6:
            continue
        direction = dx / dist
        start = s_pos + s_rad * direction
        end = r_pos - r_rad * direction

        # Curvature: opposite for bidirectional
        pair_key = frozenset([row['sender'], row['receiver']])
        if pair_key in seen_pairs:
            rad = -0.2
        else:
            rad = 0.2
            seen_pairs.add(pair_key)

        width = _edge_width(row[metric_col])
        color = _edge_color(row[metric_col])

        arrow = FancyArrowPatch(
            posA=tuple(start), posB=tuple(end),
            connectionstyle=f"arc3,rad={rad}",
            arrowstyle="->,head_length=6,head_width=4",
            color=color, linewidth=width,
            alpha=edge_alpha, zorder=2
        )
        ax.add_patch(arrow)

    # Draw self-loops (within-type)
    if show_self_loops:
        for _, row in within_edges.iterrows():
            ct = row['sender']
            x, y = node_positions[ct]
            angle = node_angles[ct]
            nr = node_radii[ct]

            width = _edge_width(row[metric_col])
            color = _edge_color(row[metric_col])

            # Arc center outside the node
            loop_cx = x + (nr + self_loop_scale * 0.6) * np.cos(angle)
            loop_cy = y + (nr + self_loop_scale * 0.6) * np.sin(angle)

            arc = Arc(
                (loop_cx, loop_cy),
                self_loop_scale, self_loop_scale,
                angle=np.degrees(angle),
                theta1=30, theta2=330,
                color=color, linewidth=width,
                alpha=edge_alpha, zorder=2
            )
            ax.add_patch(arc)

    # Draw nodes
    for ct in all_types:
        x, y = node_positions[ct]
        nr = node_radii[ct]
        color = default_node_color
        if node_colors and ct in node_colors:
            color = node_colors[ct]

        circle = plt.Circle(
            (x, y), nr,
            facecolor=color, edgecolor='black',
            linewidth=1.5, zorder=5
        )
        ax.add_patch(circle)

    # Draw labels outside nodes
    label_pad = 0.06
    for i, ct in enumerate(all_types):
        x, y = node_positions[ct]
        nr = node_radii[ct]
        angle = angles[i]

        lx = (radius + nr + label_pad) * np.cos(angle)
        ly = (radius + nr + label_pad) * np.sin(angle)

        # Determine alignment
        cos_a = np.cos(angle)
        if cos_a > 0.1:
            ha = 'left'
        elif cos_a < -0.1:
            ha = 'right'
        else:
            ha = 'center'

        # Rotate text to follow circle tangent
        rot = np.degrees(angle)
        if rot > 90:
            rot -= 180
        elif rot < -90:
            rot += 180

        ax.text(lx, ly, ct, fontsize=node_label_fontsize,
                fontweight='bold', ha=ha, va='center',
                rotation=rot, rotation_mode='anchor', zorder=6)

    # Set limits
    pad = 0.5
    ax.set_xlim(-radius - pad, radius + pad)
    ax.set_ylim(-radius - pad, radius + pad)

    # Edge colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap_obj,
                                norm=plt.Normalize(vmin=e_min, vmax=e_max))
    sm.set_array([])
    edge_label = '# Significant LR pairs' if edge_metric == 'count' else 'Mean z-score'
    cbar = fig.colorbar(sm, ax=ax, shrink=0.6, label=edge_label, pad=0.02)

    if title is None:
        title = 'Cell-Type Interaction Network'
    ax.set_title(title, fontsize=12, fontweight='bold', pad=20)

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight', facecolor='white')
        print(f"Saved: {save_path}")

    return fig, ax


def plot_lr_dotplot(results_df, lr_pairs=None, celltype_pairs=None,
                    top_n_lr=20, top_n_pairs=20,
                    size_col='p_adj', color_col='z_score',
                    sig_only=True, fdr_threshold=0.05,
                    cluster_rows=False, cluster_cols=False,
                    cmap='RdBu_r', center=0,
                    max_dot_size=300, min_dot_size=10,
                    figsize=None, title=None, save_path=None):
    """
    CellPhoneDB-style dot plot of LR pairs × cell-type pairs.

    Rows = LR pairs, columns = cell-type pairs (sender → receiver).
    Dot size = -log10(p_adj), dot color = z_score or fold_enrichment.

    Parameters
    ----------
    results_df : DataFrame
        Results from run_celltype_analysis().
    lr_pairs : list of tuples, optional
        Specific (ligand, receptor) pairs to show. If None, auto-selects top_n_lr.
    celltype_pairs : list of tuples, optional
        Specific (sender, receiver) pairs to show. If None, auto-selects top_n_pairs.
    top_n_lr : int
        Number of top LR pairs to show (by minimum p_adj).
    top_n_pairs : int
        Number of top cell-type pairs to show (by count of significant hits).
    size_col : str
        Column for dot size. 'p_adj' uses -log10(p_adj).
    color_col : str
        Column for dot color: 'z_score', 'fold_enrichment'.
    sig_only : bool
        Filter to significant results for auto-selection.
    fdr_threshold : float
        Significance threshold.
    cluster_rows, cluster_cols : bool
        Hierarchically cluster rows/columns.
    cmap : str
        Colormap for dot color.
    center : float or None
        Center value for diverging colormaps.
    max_dot_size, min_dot_size : float
        Dot size range in points^2.
    figsize : tuple or None
        Figure size. Auto-computed if None.
    title : str or None
        Plot title.
    save_path : str or None
        Path to save figure.

    Returns
    -------
    fig, ax
    """
    df = results_df.copy()

    # Filter to significant for auto-selection
    if sig_only and 'p_adj' in df.columns:
        sig_df = df[(df['p_adj'] < fdr_threshold) & (df['z_score'] > 0)]
    else:
        sig_df = df

    if len(sig_df) == 0:
        print("No significant results to plot")
        return None, None

    # Create pair labels
    sig_df = sig_df.copy()
    sig_df['lr_label'] = sig_df['ligand'] + ' → ' + sig_df['receptor']
    sig_df['ct_label'] = sig_df.apply(
        lambda r: r['sender'] if r['sender'] == r['receiver']
        else f"{r['sender']} → {r['receiver']}", axis=1)

    # Select LR pairs
    if lr_pairs is not None:
        lr_labels = [f"{l} → {r}" for l, r in lr_pairs]
        sig_df = sig_df[sig_df['lr_label'].isin(lr_labels)]
    else:
        # Top LR pairs by minimum p_adj
        best_p = sig_df.groupby('lr_label')['p_adj'].min().sort_values()
        lr_labels = best_p.head(top_n_lr).index.tolist()
        sig_df = sig_df[sig_df['lr_label'].isin(lr_labels)]

    # Select cell-type pairs
    if celltype_pairs is not None:
        ct_labels = []
        for s, r in celltype_pairs:
            ct_labels.append(s if s == r else f"{s} → {r}")
        sig_df = sig_df[sig_df['ct_label'].isin(ct_labels)]
    else:
        # Top cell-type pairs by number of significant hits
        ct_counts = sig_df.groupby('ct_label').size().sort_values(ascending=False)
        ct_labels = ct_counts.head(top_n_pairs).index.tolist()
        sig_df = sig_df[sig_df['ct_label'].isin(ct_labels)]

    if len(sig_df) == 0:
        print("No results matching selection criteria")
        return None, None

    # Get unique labels
    lr_labels = sorted(sig_df['lr_label'].unique())
    ct_labels = sorted(sig_df['ct_label'].unique())

    # Hierarchical clustering
    if (cluster_rows or cluster_cols) and len(lr_labels) > 2 and len(ct_labels) > 2:
        from scipy.cluster.hierarchy import linkage, leaves_list

        # Build pivot matrix for clustering
        pivot = sig_df.pivot_table(
            index='lr_label', columns='ct_label',
            values=color_col, aggfunc='first', fill_value=0)
        pivot = pivot.reindex(index=lr_labels, columns=ct_labels, fill_value=0)

        if cluster_rows and len(lr_labels) > 2:
            Z = linkage(pivot.values, method='average')
            order = leaves_list(Z)
            lr_labels = [lr_labels[i] for i in order]

        if cluster_cols and len(ct_labels) > 2:
            Z = linkage(pivot.values.T, method='average')
            order = leaves_list(Z)
            ct_labels = [ct_labels[i] for i in order]

    n_lr = len(lr_labels)
    n_ct = len(ct_labels)
    lr_to_idx = {l: i for i, l in enumerate(lr_labels)}
    ct_to_idx = {c: i for i, c in enumerate(ct_labels)}

    # Compute size values: -log10(p_adj)
    sig_df = sig_df.copy()
    sig_df['neg_log_p'] = -np.log10(sig_df['p_adj'].clip(lower=1e-300))

    # Size scaling
    size_vals = sig_df['neg_log_p'].values
    s_min, s_max = size_vals.min(), size_vals.max()

    def _scale_size(vals):
        if s_max > s_min:
            return min_dot_size + (vals - s_min) / (s_max - s_min) * (max_dot_size - min_dot_size)
        return np.full(len(vals), (max_dot_size + min_dot_size) / 2)

    # Color range
    color_vals = sig_df[color_col].values
    if center is not None:
        vmax = max(abs(color_vals.max() - center), abs(color_vals.min() - center))
        vmin_c, vmax_c = center - vmax, center + vmax
    else:
        vmin_c, vmax_c = color_vals.min(), color_vals.max()

    # Auto figsize
    if figsize is None:
        w = max(8, n_ct * 0.7 + 4)
        h = max(6, n_lr * 0.35 + 2)
        figsize = (w, h)

    fig, ax = plt.subplots(figsize=figsize)
    ax.set_facecolor('#FAFAF8')

    # Grid
    for i in range(n_lr):
        ax.axhline(i, color='#E8E8E8', linewidth=0.5, zorder=0)
    for j in range(n_ct):
        ax.axvline(j, color='#E8E8E8', linewidth=0.5, zorder=0)

    # Plot dots
    x_pos = [ct_to_idx[c] for c in sig_df['ct_label'] if c in ct_to_idx]
    y_pos = [lr_to_idx[l] for l in sig_df['lr_label'] if l in lr_to_idx]

    # Filter to rows with valid indices
    valid_mask = sig_df['ct_label'].isin(ct_to_idx) & sig_df['lr_label'].isin(lr_to_idx)
    plot_df = sig_df[valid_mask]

    x_pos = [ct_to_idx[c] for c in plot_df['ct_label']]
    y_pos = [lr_to_idx[l] for l in plot_df['lr_label']]
    sizes = _scale_size(plot_df['neg_log_p'].values)
    colors = plot_df[color_col].values

    sc = ax.scatter(x_pos, y_pos, s=sizes, c=colors,
                    cmap=cmap, vmin=vmin_c, vmax=vmax_c,
                    edgecolors='black', linewidths=0.5,
                    zorder=3)

    # Colorbar
    cbar = fig.colorbar(sc, ax=ax, shrink=0.7, label=color_col, pad=0.02)

    # Size legend
    if s_max > s_min:
        ref_vals = np.linspace(s_min, s_max, 4)
    else:
        ref_vals = [s_min]
    ref_sizes = _scale_size(ref_vals)

    legend_elements = []
    for val, sz in zip(ref_vals, ref_sizes):
        legend_elements.append(
            plt.scatter([], [], s=sz, c='gray', edgecolors='black',
                        linewidths=0.5, label=f'{val:.1f}'))

    leg = ax.legend(handles=legend_elements, title='-log10(p_adj)',
                    loc='upper left', bbox_to_anchor=(1.15, 1.0),
                    framealpha=0.9, edgecolor='gray', fontsize=9,
                    title_fontsize=9)

    # Axis labels
    ax.set_xticks(range(n_ct))
    ax.set_xticklabels(ct_labels, rotation=45, ha='right', fontsize=9)
    ax.set_yticks(range(n_lr))
    ax.set_yticklabels(lr_labels, fontsize=9)
    ax.set_xlabel('Cell-Type Pair', fontsize=11)
    ax.set_ylabel('Ligand → Receptor', fontsize=11)
    ax.set_xlim(-0.5, n_ct - 0.5)
    ax.set_ylim(n_lr - 0.5, -0.5)

    if title is None:
        title = 'LR Interaction Dot Plot'
    ax.set_title(title, fontsize=12, fontweight='bold', pad=10)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved: {save_path}")

    return fig, ax


def plot_spatial_interactions(adata, cell_type_col, results_df=None,
                              ligand=None, receptor=None,
                              sender=None, receiver=None,
                              indices=None, distances=None, tau=5.0,
                              mode='cell_type',
                              sender_color='#E8B87D', receiver_color='#7BA3A8',
                              edge_color='#999999', edge_alpha=0.3,
                              bg_color='#E8E8E8', point_size=0.5,
                              xlim=None, ylim=None,
                              figsize=(12, 10), title=None, save_path=None):
    """
    Overlay LR interactions on tissue spatial coordinates.

    Parameters
    ----------
    adata : AnnData
        Annotated data with spatial coordinates in adata.obsm['spatial'].
    cell_type_col : str
        Column in adata.obs with cell type labels.
    results_df : DataFrame, optional
        Results from run_celltype_analysis(). Required for mode='overview'.
    ligand, receptor : str, optional
        Specific LR pair to visualize. Required for mode='cell_type' and 'score'.
    sender, receiver : str, optional
        Cell types. Required for mode='cell_type' and 'score'.
    indices, distances : array, optional
        KNN graph. Required for mode='cell_type' (edge drawing).
    tau : float
        Kernel decay distance for mode='score'.
    mode : str
        Visualization mode:
        - 'cell_type': highlight sender/receiver cells with edges
        - 'score': color sender cells by interaction score contribution
        - 'overview': bin-level heatmap of total significant interactions
    sender_color, receiver_color : str
        Colors for sender/receiver cells.
    edge_color : str
        Color for interaction edges.
    edge_alpha : float
        Edge transparency.
    bg_color : str
        Background cell color.
    point_size : float
        Size of cell dots.
    xlim, ylim : tuple, optional
        (min, max) for axis limits (zoom).
    figsize : tuple
        Figure size.
    title : str, optional
        Plot title.
    save_path : str, optional
        Path to save figure.

    Returns
    -------
    fig, ax
    """
    from scipy.sparse import issparse
    from matplotlib.collections import LineCollection

    coords = adata.obsm['spatial']
    cell_types = adata.obs[cell_type_col].values.astype(str)

    fig, ax = plt.subplots(figsize=figsize)

    if mode == 'cell_type':
        # --- Mode 1: highlight sender/receiver with edges ---
        if sender is None or receiver is None:
            raise ValueError("sender and receiver required for mode='cell_type'")

        sender_mask = cell_types == sender
        receiver_mask = cell_types == receiver
        other_mask = ~sender_mask & ~receiver_mask

        # Background cells
        ax.scatter(coords[other_mask, 0], coords[other_mask, 1],
                   s=point_size * 0.3, c=bg_color, alpha=0.3,
                   rasterized=True, zorder=1)

        # If within-type, just highlight that type
        if sender == receiver:
            ax.scatter(coords[sender_mask, 0], coords[sender_mask, 1],
                       s=point_size, c=sender_color, alpha=0.7,
                       rasterized=True, zorder=2, label=f'{sender}')
        else:
            ax.scatter(coords[sender_mask, 0], coords[sender_mask, 1],
                       s=point_size, c=sender_color, alpha=0.7,
                       rasterized=True, zorder=2, label=f'{sender} (sender)')
            ax.scatter(coords[receiver_mask, 0], coords[receiver_mask, 1],
                       s=point_size, c=receiver_color, alpha=0.7,
                       rasterized=True, zorder=2, label=f'{receiver} (receiver)')

        # Draw edges between sender-receiver spatial neighbors
        if indices is not None:
            sender_idx = np.where(sender_mask)[0]
            receiver_set = set(np.where(receiver_mask)[0])

            # Collect edges
            edges = []
            for si in sender_idx:
                for ni in indices[si]:
                    if ni in receiver_set:
                        edges.append([coords[si], coords[ni]])

            if edges:
                # Subsample if too many edges
                max_edges = 50000
                if len(edges) > max_edges:
                    rng = np.random.RandomState(42)
                    idx = rng.choice(len(edges), max_edges, replace=False)
                    edges = [edges[i] for i in idx]

                lc = LineCollection(edges, colors=edge_color,
                                    alpha=edge_alpha, linewidths=0.3,
                                    zorder=1)
                ax.add_collection(lc)
                edge_info = f'{len(edges):,} edges shown'
            else:
                edge_info = 'no edges'
        else:
            edge_info = 'no graph provided'

        ax.legend(loc='upper right', fontsize=10, framealpha=0.9,
                  markerscale=5, edgecolor='gray')

        if title is None:
            lr_str = f' ({ligand}→{receptor})' if ligand and receptor else ''
            title = f'{sender} → {receiver}{lr_str}\n{edge_info}'

    elif mode == 'score':
        # --- Mode 2: color sender cells by interaction score ---
        if sender is None or receiver is None or ligand is None or receptor is None:
            raise ValueError("sender, receiver, ligand, receptor required for mode='score'")
        if indices is None or distances is None:
            raise ValueError("indices and distances required for mode='score'")

        sender_mask = cell_types == sender
        receiver_mask = cell_types == receiver
        sender_idx = np.where(sender_mask)[0]
        receiver_idx = np.where(receiver_mask)[0]

        # Gene indices
        gene_names = list(adata.var_names)
        lig_gi = gene_names.index(ligand)
        rec_gi = gene_names.index(receptor)

        X = adata.X
        if issparse(X):
            X_csc = X.tocsc()
            L_sender = np.log1p(np.asarray(X_csc[sender_idx, lig_gi].todense()).flatten())
            R_receiver = np.log1p(np.asarray(X_csc[receiver_idx, rec_gi].todense()).flatten())
        else:
            L_sender = np.log1p(X[sender_idx, lig_gi])
            R_receiver = np.log1p(X[receiver_idx, rec_gi])

        # Compute w vector for each sender cell
        n_total = len(cell_types)
        is_receiver = np.zeros(n_total, dtype=bool)
        is_receiver[receiver_idx] = True

        nbr_indices = indices[sender_idx, :]
        nbr_dists = distances[sender_idx, :]
        mask = is_receiver[nbr_indices]
        row_pos, col_pos = np.where(mask)

        edge_dists = nbr_dists[row_pos, col_pos]
        K = np.exp(-edge_dists / tau)

        local_edges_j = np.searchsorted(receiver_idx, nbr_indices[row_pos, col_pos])
        KR = K * R_receiver[local_edges_j]

        w = np.zeros(len(sender_idx))
        np.add.at(w, row_pos, KR)

        # Interaction score per sender cell = L * w
        scores = L_sender * w

        # Background
        other_mask = ~sender_mask & ~receiver_mask
        ax.scatter(coords[other_mask, 0], coords[other_mask, 1],
                   s=point_size * 0.3, c=bg_color, alpha=0.3,
                   rasterized=True, zorder=1)

        # Receiver cells
        ax.scatter(coords[receiver_mask, 0], coords[receiver_mask, 1],
                   s=point_size * 0.5, c=receiver_color, alpha=0.3,
                   rasterized=True, zorder=1, label=f'{receiver} (receiver)')

        # Sender cells colored by score
        sc = ax.scatter(coords[sender_idx, 0], coords[sender_idx, 1],
                        s=point_size, c=scores, cmap='YlOrRd',
                        vmin=0, rasterized=True, zorder=3)

        cbar = fig.colorbar(sc, ax=ax, shrink=0.6,
                            label=f'{ligand}→{receptor} score', pad=0.02)

        ax.legend(loc='upper right', fontsize=10, framealpha=0.9,
                  markerscale=5, edgecolor='gray')

        if title is None:
            title = f'{sender}: {ligand}→{receptor} interaction score'

    elif mode == 'overview':
        # --- Mode 3: spatial heatmap of interaction density ---
        if results_df is None:
            raise ValueError("results_df required for mode='overview'")

        sig = results_df[(results_df['p_adj'] < 0.05) & (results_df['z_score'] > 0)]

        # Count significant interactions per cell type
        ct_sig_count = {}
        for _, row in sig.iterrows():
            s, r = row['sender'], row['receiver']
            ct_sig_count[s] = ct_sig_count.get(s, 0) + 1
            if s != r:
                ct_sig_count[r] = ct_sig_count.get(r, 0) + 1

        # Assign score to each cell
        cell_scores = np.zeros(len(cell_types))
        for ct, count in ct_sig_count.items():
            mask = cell_types == ct
            cell_scores[mask] = count

        # Plot
        nonzero = cell_scores > 0
        ax.scatter(coords[~nonzero, 0], coords[~nonzero, 1],
                   s=point_size * 0.3, c=bg_color, alpha=0.3,
                   rasterized=True, zorder=1)
        sc = ax.scatter(coords[nonzero, 0], coords[nonzero, 1],
                        s=point_size, c=cell_scores[nonzero],
                        cmap='YlOrRd', rasterized=True, zorder=2)

        cbar = fig.colorbar(sc, ax=ax, shrink=0.6,
                            label='# Significant interactions', pad=0.02)

        if title is None:
            title = f'Spatial overview: {len(sig):,} significant interactions'

    else:
        raise ValueError(f"Unknown mode: {mode}. Use 'cell_type', 'score', or 'overview'.")

    # Common formatting
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)

    ax.invert_yaxis()
    ax.set_aspect('equal')
    ax.set_xlabel('X (μm)', fontsize=11)
    ax.set_ylabel('Y (μm)', fontsize=11)
    ax.set_title(title, fontsize=12, fontweight='bold', pad=10)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight', facecolor='white')
        print(f"Saved: {save_path}")

    return fig, ax


# =============================================================================
# COMPARTMENT-LEVEL VISUALIZATION
# =============================================================================

def plot_compartment_heatmap(
    results_df,
    value_col='significant',
    compartment_col='compartment',
    lr_col='lr_pair',
    compartment_order=None,
    cmap='RdYlGn',
    figsize=None,
    title='CONSTELLATION Compartment Detection',
    show_values=True,
    save_path=None,
):
    """
    Plot a heatmap of LR detection across compartments.

    Parameters
    ----------
    results_df : pd.DataFrame
        Results from run_compartment_analysis().
    value_col : str, default='significant'
        Column to use for heatmap values. Options:
        - 'significant': boolean detection (True/False)
        - 'fold_enrichment': fold enrichment values
        - 'z_score': z-score values
    compartment_col : str, default='compartment'
        Column containing compartment labels.
    lr_col : str, default='lr_pair'
        Column containing LR pair names.
    compartment_order : list, optional
        Order of compartments (columns). If None, uses sorted order.
    cmap : str, default='RdYlGn'
        Colormap for heatmap.
    figsize : tuple, optional
        Figure size. If None, auto-computed based on data.
    title : str
        Plot title.
    show_values : bool, default=True
        Show text annotations on cells.
    save_path : str, optional
        Path to save figure.

    Returns
    -------
    fig, ax : matplotlib figure and axes
    """
    import pandas as pd

    # Pivot data
    pivot = results_df.pivot_table(
        index=lr_col,
        columns=compartment_col,
        values=value_col,
        aggfunc='first'
    )

    if compartment_order is not None:
        pivot = pivot[[c for c in compartment_order if c in pivot.columns]]

    # Figure size
    if figsize is None:
        n_rows, n_cols = pivot.shape
        figsize = (max(4, n_cols * 1.5), max(6, n_rows * 0.4))

    fig, ax = plt.subplots(figsize=figsize)

    # Determine value range
    if value_col == 'significant':
        vmin, vmax = 0, 1
        data = pivot.values.astype(float)
    else:
        data = pivot.values
        vmin, vmax = None, None

    im = ax.imshow(data, cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)

    # Labels
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, fontsize=11)
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index, fontsize=10)
    ax.set_xlabel('Compartment', fontsize=12)
    ax.set_ylabel('LR Pair', fontsize=12)
    ax.set_title(title, fontsize=13, fontweight='bold')

    # Annotations
    if show_values:
        for i in range(len(pivot.index)):
            for j in range(len(pivot.columns)):
                val = pivot.values[i, j]
                if value_col == 'significant':
                    text = '+' if val else '-'
                    color = 'white' if val else 'gray'
                else:
                    text = f'{val:.1f}' if not np.isnan(val) else ''
                    color = 'white' if not np.isnan(val) and val > (vmax or data.max()) / 2 else 'black'
                ax.text(j, i, text, ha='center', va='center',
                        color=color, fontsize=11, fontweight='bold')

    # Colorbar (only for non-boolean)
    if value_col != 'significant':
        cbar = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
        cbar.set_label(value_col.replace('_', ' ').title(), fontsize=11)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight', facecolor='white')
        print(f"Saved: {save_path}")

    return fig, ax


def plot_distance_profile(
    profile_data,
    ligand=None,
    receptor=None,
    boundary_pos=0,
    interface_pos=None,
    color='#5a7aaa',
    figsize=(10, 5),
    title=None,
    xlabel='Distance from Boundary (μm)',
    ylabel='Coexpression Fold Enrichment',
    show_reference=True,
    save_path=None,
):
    """
    Plot LR colocalization as a function of distance from a boundary.

    Parameters
    ----------
    profile_data : dict
        Output from compute_distance_profile() with keys:
        - bin_centers: array of bin centers
        - coexpr_fold: fold enrichment values
        - l_fraction: ligand fraction (optional)
        - r_fraction: receptor fraction (optional)
    ligand : str, optional
        Ligand name for title.
    receptor : str, optional
        Receptor name for title.
    boundary_pos : float, default=0
        Position of the main boundary (e.g., tumor edge).
    interface_pos : float, optional
        Position of secondary boundary (e.g., interface edge).
    color : str, default='#5a7aaa'
        Line color.
    figsize : tuple
        Figure size.
    title : str, optional
        Plot title. If None, auto-generated from ligand/receptor.
    xlabel : str
        X-axis label.
    ylabel : str
        Y-axis label.
    show_reference : bool, default=True
        Show horizontal line at y=1 (no enrichment).
    save_path : str, optional
        Path to save figure.

    Returns
    -------
    fig, ax : matplotlib figure and axes
    """
    bin_centers = profile_data['bin_centers']
    coexpr_fold = profile_data['coexpr_fold']

    fig, ax = plt.subplots(figsize=figsize)

    # Main line
    ax.plot(bin_centers, coexpr_fold, 'o-', color=color,
            linewidth=2, markersize=4, label='Coexpression fold')

    # Reference line
    if show_reference:
        ax.axhline(y=1, color='gray', linestyle='--', alpha=0.5, label='Expected')

    # Boundary markers
    ax.axvline(x=boundary_pos, color='red', linestyle='--', alpha=0.4,
               label='Boundary')
    if interface_pos is not None:
        ax.axvline(x=interface_pos, color='orange', linestyle='--', alpha=0.4,
                   label='Interface')

    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)

    if title is None and ligand and receptor:
        title = f'{ligand}-{receptor} Colocalization vs Distance'
    if title:
        ax.set_title(title, fontsize=13, fontweight='bold')

    ax.legend(loc='upper right', fontsize=10)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight', facecolor='white')
        print(f"Saved: {save_path}")

    return fig, ax


def plot_compartment_spatial(
    coords,
    compartments,
    compartment_colors=None,
    point_size=0.5,
    figsize=(10, 8),
    title='Spatial Compartments',
    xlim=None,
    ylim=None,
    save_path=None,
):
    """
    Plot spatial map colored by compartment.

    Parameters
    ----------
    coords : array-like
        Spatial coordinates, shape (n_cells, 2).
    compartments : array-like
        Compartment labels for each cell.
    compartment_colors : dict, optional
        Mapping of compartment names to colors.
        If None, uses default colors.
    point_size : float
        Size of points.
    figsize : tuple
        Figure size.
    title : str
        Plot title.
    xlim : tuple, optional
        X-axis limits (min, max).
    ylim : tuple, optional
        Y-axis limits (min, max).
    save_path : str, optional
        Path to save figure.

    Returns
    -------
    fig, ax : matplotlib figure and axes
    """
    coords = np.asarray(coords)
    compartments = np.asarray(compartments)
    unique_comps = sorted(set(compartments))

    # Default colors
    if compartment_colors is None:
        default_colors = ['#d4a0a0', '#d4c090', '#a0b8d4', '#a0d4a0',
                          '#d4a0d4', '#a0d4d4', '#d4d4a0', '#b0b0b0']
        compartment_colors = {c: default_colors[i % len(default_colors)]
                              for i, c in enumerate(unique_comps)}

    fig, ax = plt.subplots(figsize=figsize)

    for comp in unique_comps:
        mask = compartments == comp
        color = compartment_colors.get(comp, '#808080')
        ax.scatter(coords[mask, 0], coords[mask, 1],
                   c=color, s=point_size, alpha=0.5, label=comp, rasterized=True)

    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)

    ax.set_aspect('equal')
    ax.set_xlabel('X (μm)', fontsize=11)
    ax.set_ylabel('Y (μm)', fontsize=11)
    ax.set_title(title, fontsize=13, fontweight='bold')
    ax.legend(markerscale=20, loc='upper right', fontsize=10)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight', facecolor='white')
        print(f"Saved: {save_path}")

    return fig, ax


def plot_boundary_profile(profile_data, ligand=None, receptor=None,
                          boundary_pos=0, interface_width=50,
                          xlim=None, zone_labels=True,
                          figsize=(8, 8), color='#5a7aaa', save_path=None):
    """
    Plot a 3-row boundary profile from compute_distance_profile() output.

    Row 1: Coexpression fraction (L+R+ / L+) vs distance.
    Row 2: Number of ligand-positive cells per bin (bar chart).
    Row 3: Spatial fold enrichment vs distance with y=1 reference line.

    Parameters
    ----------
    profile_data : dict
        Output of compute_distance_profile(). Must contain keys:
        bin_centers, coexpr_frac, n_l_pos, spatial_fold.
    ligand : str, optional
        Ligand gene name (for labels).
    receptor : str, optional
        Receptor gene name (for labels).
    boundary_pos : float, default=0
        Distance value of the tumor boundary (for vertical line).
    interface_width : float, default=50
        Width of the interface zone (for shading).
    xlim : tuple, optional
        (xmin, xmax) limits for the x-axis. If None, uses data range.
    zone_labels : bool, default=True
        Show Tumor / Interface / Tissue labels at the top of the plot.
    figsize : tuple, default=(8, 8)
        Figure size.
    color : str, default='#5a7aaa'
        Line/bar color.
    save_path : str, optional
        Path to save figure.

    Returns
    -------
    fig, axes : matplotlib figure and axes array
    """
    centers = profile_data['bin_centers']
    coexpr_frac = profile_data['coexpr_frac']
    n_l_pos = profile_data['n_l_pos']
    spatial_fold = profile_data['spatial_fold']

    pair_label = ''
    if ligand and receptor:
        pair_label = f'{ligand}\u2013{receptor}'

    fig, axes = plt.subplots(3, 1, figsize=figsize, sharex=True,
                             gridspec_kw={'height_ratios': [2, 1, 2]})

    # Helper: add boundary line, interface shading, and zone labels
    def _decorate(ax, is_top=False):
        ax.axvline(boundary_pos, color='red', ls='--', lw=1, alpha=0.7)
        ax.axvspan(boundary_pos, boundary_pos + interface_width,
                   color='orange', alpha=0.08)
        if zone_labels and is_top:
            # Place zone labels at the top of the axes
            x_lo, x_hi = ax.get_xlim() if xlim is None else xlim
            trans = ax.get_xaxis_transform()
            # Tumor label (left of boundary)
            tumor_x = (x_lo + boundary_pos) / 2
            if tumor_x >= x_lo:
                ax.text(tumor_x, 1.05, 'Tumor', transform=trans,
                        ha='center', va='bottom', fontsize=10, color='#c0392b',
                        fontweight='bold')
            # Interface label (in the shaded zone)
            iface_x = boundary_pos + interface_width / 2
            ax.text(iface_x, 1.05, 'Interface', transform=trans,
                    ha='center', va='bottom', fontsize=10, color='#e67e22',
                    fontweight='bold')
            # Tissue label (right of interface)
            tissue_x = (boundary_pos + interface_width + x_hi) / 2
            if tissue_x <= x_hi:
                ax.text(tissue_x, 1.05, 'Tissue', transform=trans,
                        ha='center', va='bottom', fontsize=10, color='#2980b9',
                        fontweight='bold')

    # Row 1: Coexpression fraction
    ax = axes[0]
    valid = ~np.isnan(coexpr_frac)
    ax.plot(centers[valid], coexpr_frac[valid], '-o', color=color,
            markersize=4, lw=1.5)
    lig_label = ligand if ligand else 'L'
    rec_label = receptor if receptor else 'R'
    ax.set_ylabel(f'{lig_label}+{rec_label}+ / {lig_label}+', fontsize=11)
    ax.set_title(f'Boundary Profile: {pair_label}' if pair_label
                 else 'Boundary Profile', fontsize=13, fontweight='bold',
                 pad=20 if zone_labels else 6)
    ax.grid(axis='y', alpha=0.3)
    if xlim:
        ax.set_xlim(xlim)
    _decorate(ax, is_top=True)

    # Row 2: n Ligand+ cells (bar chart)
    ax = axes[1]
    width = np.median(np.diff(centers)) * 0.8
    ax.bar(centers, n_l_pos, width=width, color=color, alpha=0.6,
           edgecolor='white', linewidth=0.5)
    _decorate(ax)
    ax.set_ylabel(f'n {lig_label}+ cells', fontsize=11)
    ax.grid(axis='y', alpha=0.3)

    # Row 3: Spatial fold enrichment
    ax = axes[2]
    valid_sf = ~np.isnan(spatial_fold)
    ax.plot(centers[valid_sf], spatial_fold[valid_sf], '-o', color=color,
            markersize=4, lw=1.5)
    ax.axhline(1.0, color='gray', ls=':', lw=1, alpha=0.7)
    _decorate(ax)
    ax.set_ylabel('Spatial fold enrichment', fontsize=11)
    ax.set_xlabel('Distance from tumor boundary (\u00b5m)', fontsize=11)
    ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight', facecolor='white')
        print(f"Saved: {save_path}")

    return fig, axes
