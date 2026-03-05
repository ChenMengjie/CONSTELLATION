"""
CONSTELLATION: CONtext-Specific spaTial cELLulAr inTeractions In OrgaNs

A spatial-conditioned permutation testing framework for ligand-receptor
interaction inference in spatial transcriptomics data.
"""

from .testing import (
    # High-level entry point
    run_celltype_analysis,
    run_lineage_analysis,
    validate_inputs,
    scan_celltype_pairs,
    # Compartment-level analysis
    run_compartment_analysis,
    scan_compartments,
    compute_distance_profile,
    # Batch testing
    test_within_type_lr,
    test_between_type_lr,
    compute_expression_fractions,
    filter_testable_lr_pairs,
    filter_testable_lr_pairs_between,
    # Targeted testing functions
    test_lr_pair_within_types,
    test_lr_pair_between_types,
    test_ligand_all_receptors,
    test_receptor_all_ligands,
    test_custom_lr_set,
)

from .io import (
    load_spatial_data,
    load_lr_pairs,
    load_lr_resource,
    show_lr_resources,
    build_spatial_graph,
    build_spatial_graph_from_coords,
    filter_lr_pairs_by_genes,
    apply_fdr_correction,
    save_results,
    load_cell_type_mapping,
    print_results_summary,
    print_top_results,
    # Testable burden computation
    compute_testable_burden_within,
    compute_testable_burden_between,
    save_testable_burden,
    load_testable_burden,
    # Cell size statistics
    compute_cell_sizes,
    compute_cell_sizes_with_area,
    summarize_cell_sizes,
    print_cell_size_table,
)

from .visualization import (
    plot_celltype_pair_heatmap,
    plot_lr_category_heatmap,
    plot_cell_lineage_tree,
    plot_celltype_barplot,
    plot_lr_interaction_summary,
    create_dotplot_markers,
    # Cell-type pair-wise interaction plots
    plot_combined_heatmap,
    plot_interaction_dotplot,
    plot_interaction_network,
    # LR-level and spatial plots
    plot_lr_dotplot,
    plot_spatial_interactions,
    # Compartment-level plots
    plot_compartment_heatmap,
    plot_distance_profile,
    plot_compartment_spatial,
    plot_boundary_profile,
)

from .ontology import (
    # OLS API functions
    search_cell_ontology,
    get_term_info,
    get_ancestors,
    get_children,
    # Mapping functions
    map_annotation_to_ontology,
    build_annotation_mapping,
    # Hierarchy functions
    get_immune_cell_hierarchy,
    define_refined_groups,
    get_group_at_level,
    # Lymph node specific
    get_lymph_node_refined_groups,
    get_lymph_node_coarse_groups,
    build_hierarchy_for_plotting,
    print_hierarchy,
    # Constants
    IMMUNE_CELL_ONTOLOGY,
)

__version__ = '0.4.0'
__package_name__ = 'CONSTELLATION'
__all__ = [
    # High-level
    'run_celltype_analysis',
    'run_lineage_analysis',
    'validate_inputs',
    'scan_celltype_pairs',
    # Compartment-level
    'run_compartment_analysis',
    'scan_compartments',
    'compute_distance_profile',
    # Testing - batch
    'test_within_type_lr',
    'test_between_type_lr',
    'compute_expression_fractions',
    'filter_testable_lr_pairs',
    'filter_testable_lr_pairs_between',
    # Testing - targeted
    'test_lr_pair_within_types',
    'test_lr_pair_between_types',
    'test_ligand_all_receptors',
    'test_receptor_all_ligands',
    'test_custom_lr_set',
    # I/O
    'load_spatial_data',
    'load_lr_pairs',
    'load_lr_resource',
    'show_lr_resources',
    'build_spatial_graph',
    'build_spatial_graph_from_coords',
    'filter_lr_pairs_by_genes',
    'apply_fdr_correction',
    'save_results',
    'load_cell_type_mapping',
    'print_results_summary',
    'print_top_results',
    # Testable burden
    'compute_testable_burden_within',
    'compute_testable_burden_between',
    'save_testable_burden',
    'load_testable_burden',
    # Cell size statistics
    'compute_cell_sizes',
    'compute_cell_sizes_with_area',
    'summarize_cell_sizes',
    'print_cell_size_table',
    # Visualization
    'plot_celltype_pair_heatmap',
    'plot_lr_category_heatmap',
    'plot_cell_lineage_tree',
    'plot_celltype_barplot',
    'plot_lr_interaction_summary',
    'create_dotplot_markers',
    'plot_combined_heatmap',
    'plot_interaction_dotplot',
    'plot_interaction_network',
    'plot_lr_dotplot',
    'plot_spatial_interactions',
    # Compartment visualization
    'plot_compartment_heatmap',
    'plot_distance_profile',
    'plot_compartment_spatial',
    'plot_boundary_profile',
    # Cell Ontology
    'search_cell_ontology',
    'get_term_info',
    'get_ancestors',
    'get_children',
    'map_annotation_to_ontology',
    'build_annotation_mapping',
    'get_immune_cell_hierarchy',
    'define_refined_groups',
    'get_group_at_level',
    'get_lymph_node_refined_groups',
    'get_lymph_node_coarse_groups',
    'build_hierarchy_for_plotting',
    'print_hierarchy',
    'IMMUNE_CELL_ONTOLOGY',
]
