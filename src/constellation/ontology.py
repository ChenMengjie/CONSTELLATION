"""
Cell Ontology Functions

Functions for querying Cell Ontology (CL) and building cell lineage hierarchies.
Uses the EBI Ontology Lookup Service (OLS) API for ontology queries.
"""

import requests
import time
from typing import Dict, List, Optional, Tuple, Set
from collections import defaultdict


# =============================================================================
# ONTOLOGY LOOKUP SERVICE (OLS) API
# =============================================================================

OLS_BASE_URL = "https://www.ebi.ac.uk/ols4/api"
CELL_ONTOLOGY_ID = "cl"


def search_cell_ontology(query: str, exact: bool = False,
                         max_results: int = 10) -> List[Dict]:
    """
    Search Cell Ontology for terms matching a query.

    Parameters
    ----------
    query : str
        Search term (e.g., "T cell", "B lymphocyte")
    exact : bool
        If True, only return exact matches
    max_results : int
        Maximum number of results to return

    Returns
    -------
    list of dict
        List of matching terms with 'id', 'label', 'description', 'synonyms'
    """
    url = f"{OLS_BASE_URL}/search"
    params = {
        "q": query,
        "ontology": CELL_ONTOLOGY_ID,
        "rows": max_results,
        "exact": str(exact).lower(),
    }

    try:
        response = requests.get(url, params=params, timeout=10)
        response.raise_for_status()
        data = response.json()

        results = []
        for doc in data.get("response", {}).get("docs", []):
            results.append({
                "id": doc.get("obo_id", doc.get("short_form")),
                "label": doc.get("label"),
                "description": doc.get("description", [""])[0] if doc.get("description") else "",
                "synonyms": doc.get("synonym", []),
                "iri": doc.get("iri"),
            })
        return results

    except requests.RequestException as e:
        print(f"Error querying OLS: {e}")
        return []


def get_term_info(cl_id: str) -> Optional[Dict]:
    """
    Get detailed information about a Cell Ontology term.

    Parameters
    ----------
    cl_id : str
        Cell Ontology ID (e.g., "CL:0000084" for T cell)

    Returns
    -------
    dict or None
        Term information including label, definition, synonyms
    """
    # Convert CL:0000084 to CL_0000084 for URL
    term_id = cl_id.replace(":", "_")
    url = f"{OLS_BASE_URL}/ontologies/{CELL_ONTOLOGY_ID}/terms"
    params = {"short_form": term_id}

    try:
        response = requests.get(url, params=params, timeout=10)
        response.raise_for_status()
        data = response.json()

        terms = data.get("_embedded", {}).get("terms", [])
        if terms:
            term = terms[0]
            return {
                "id": cl_id,
                "label": term.get("label"),
                "description": term.get("description", [""])[0] if term.get("description") else "",
                "synonyms": term.get("synonyms", []),
                "iri": term.get("iri"),
            }
        return None

    except requests.RequestException as e:
        print(f"Error getting term info: {e}")
        return None


def get_ancestors(cl_id: str, include_self: bool = False) -> List[Dict]:
    """
    Get ancestor terms (parents, grandparents, etc.) in the ontology.

    Parameters
    ----------
    cl_id : str
        Cell Ontology ID
    include_self : bool
        Whether to include the term itself

    Returns
    -------
    list of dict
        Ancestor terms ordered from immediate parent to root
    """
    term_id = cl_id.replace(":", "_")
    url = f"{OLS_BASE_URL}/ontologies/{CELL_ONTOLOGY_ID}/terms/http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2F{term_id}/ancestors"

    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        data = response.json()

        ancestors = []
        for term in data.get("_embedded", {}).get("terms", []):
            obo_id = term.get("obo_id")
            if obo_id and obo_id.startswith("CL:"):
                ancestors.append({
                    "id": obo_id,
                    "label": term.get("label"),
                })

        if include_self:
            info = get_term_info(cl_id)
            if info:
                ancestors.insert(0, {"id": cl_id, "label": info["label"]})

        return ancestors

    except requests.RequestException as e:
        print(f"Error getting ancestors: {e}")
        return []


def get_children(cl_id: str) -> List[Dict]:
    """
    Get direct children of a Cell Ontology term.

    Parameters
    ----------
    cl_id : str
        Cell Ontology ID

    Returns
    -------
    list of dict
        Direct child terms
    """
    term_id = cl_id.replace(":", "_")
    url = f"{OLS_BASE_URL}/ontologies/{CELL_ONTOLOGY_ID}/terms/http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2F{term_id}/children"

    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        data = response.json()

        children = []
        for term in data.get("_embedded", {}).get("terms", []):
            obo_id = term.get("obo_id")
            if obo_id and obo_id.startswith("CL:"):
                children.append({
                    "id": obo_id,
                    "label": term.get("label"),
                })
        return children

    except requests.RequestException as e:
        print(f"Error getting children: {e}")
        return []


# =============================================================================
# COMMON CELL TYPE MAPPINGS (for lymph node / immune cells)
# =============================================================================

# Pre-defined mappings for common immune cell types
IMMUNE_CELL_ONTOLOGY = {
    # T cells
    "T cell": "CL:0000084",
    "CD4+ T cell": "CL:0000624",
    "CD8+ T cell": "CL:0000625",
    "naive T cell": "CL:0000898",
    "naive CD4+ T cell": "CL:0000895",
    "naive CD8+ T cell": "CL:0000900",
    "memory T cell": "CL:0000813",
    "effector T cell": "CL:0000911",
    "regulatory T cell": "CL:0000815",
    "Treg": "CL:0000815",
    "T follicular helper cell": "CL:0002038",
    "Tfh": "CL:0002038",
    "T follicular regulatory cell": "CL:0002039",
    "Tfr": "CL:0002039",

    # B cells
    "B cell": "CL:0000236",
    "naive B cell": "CL:0000788",
    "memory B cell": "CL:0000787",
    "germinal center B cell": "CL:0000844",
    "GC B cell": "CL:0000844",
    "plasma cell": "CL:0000786",
    "plasmablast": "CL:0000980",

    # Myeloid cells
    "monocyte": "CL:0000576",
    "macrophage": "CL:0000235",
    "dendritic cell": "CL:0000451",
    "conventional dendritic cell": "CL:0000990",
    "cDC": "CL:0000990",
    "plasmacytoid dendritic cell": "CL:0000784",
    "pDC": "CL:0000784",

    # NK and other
    "NK cell": "CL:0000623",
    "natural killer cell": "CL:0000623",
    "mast cell": "CL:0000097",
    "innate lymphoid cell": "CL:0001065",
    "ILC": "CL:0001065",
}


# =============================================================================
# LINEAGE HIERARCHY DEFINITIONS
# =============================================================================

def get_immune_cell_hierarchy() -> Dict:
    """
    Get a predefined immune cell lineage hierarchy.

    Returns
    -------
    dict
        Nested dictionary representing the hierarchy
    """
    return {
        "Lymphoid": {
            "cl_id": "CL:0000542",
            "children": {
                "T cell": {
                    "cl_id": "CL:0000084",
                    "children": {
                        "CD4+ T cell": {
                            "cl_id": "CL:0000624",
                            "children": {
                                "naive CD4+ T cell": {"cl_id": "CL:0000895"},
                                "memory CD4+ T cell": {"cl_id": "CL:0000897"},
                                "Treg": {"cl_id": "CL:0000815"},
                                "Tfh": {"cl_id": "CL:0002038"},
                                "Tfr": {"cl_id": "CL:0002039"},
                            }
                        },
                        "CD8+ T cell": {
                            "cl_id": "CL:0000625",
                            "children": {
                                "naive CD8+ T cell": {"cl_id": "CL:0000900"},
                                "memory CD8+ T cell": {"cl_id": "CL:0000909"},
                                "effector CD8+ T cell": {"cl_id": "CL:0001050"},
                            }
                        },
                    }
                },
                "B cell": {
                    "cl_id": "CL:0000236",
                    "children": {
                        "naive B cell": {"cl_id": "CL:0000788"},
                        "GC B cell": {
                            "cl_id": "CL:0000844",
                            "children": {
                                "dark zone GC B cell": {"cl_id": "CL:0000846"},
                                "light zone GC B cell": {"cl_id": "CL:0000847"},
                            }
                        },
                        "memory B cell": {"cl_id": "CL:0000787"},
                        "plasma cell": {"cl_id": "CL:0000786"},
                    }
                },
                "NK cell": {"cl_id": "CL:0000623"},
                "ILC": {"cl_id": "CL:0001065"},
            }
        },
        "Myeloid": {
            "cl_id": "CL:0000763",
            "children": {
                "monocyte": {"cl_id": "CL:0000576"},
                "macrophage": {"cl_id": "CL:0000235"},
                "dendritic cell": {
                    "cl_id": "CL:0000451",
                    "children": {
                        "cDC": {"cl_id": "CL:0000990"},
                        "pDC": {"cl_id": "CL:0000784"},
                    }
                },
                "mast cell": {"cl_id": "CL:0000097"},
            }
        },
    }


# =============================================================================
# ANNOTATION MAPPING FUNCTIONS
# =============================================================================

def map_annotation_to_ontology(annotation: str,
                                custom_mappings: Dict[str, str] = None) -> Optional[str]:
    """
    Map a cell type annotation string to a Cell Ontology ID.

    Parameters
    ----------
    annotation : str
        Cell type annotation (e.g., "Naive CD4 T", "GC-Tfh-SAP")
    custom_mappings : dict, optional
        Custom mappings {annotation: cl_id}

    Returns
    -------
    str or None
        Cell Ontology ID if found
    """
    # Normalize annotation
    ann_lower = annotation.lower().strip()

    # Check custom mappings first
    if custom_mappings:
        for key, cl_id in custom_mappings.items():
            if key.lower() == ann_lower or key.lower() in ann_lower:
                return cl_id

    # Check predefined immune cell mappings
    for key, cl_id in IMMUNE_CELL_ONTOLOGY.items():
        if key.lower() == ann_lower:
            return cl_id

    # Try partial matching
    for key, cl_id in IMMUNE_CELL_ONTOLOGY.items():
        if key.lower() in ann_lower or ann_lower in key.lower():
            return cl_id

    # Try OLS search as fallback
    results = search_cell_ontology(annotation, max_results=1)
    if results:
        return results[0]["id"]

    return None


def build_annotation_mapping(annotations: List[str],
                             custom_mappings: Dict[str, str] = None,
                             verbose: bool = True) -> Dict[str, Optional[str]]:
    """
    Map a list of annotation labels to Cell Ontology IDs.

    Parameters
    ----------
    annotations : list of str
        List of cell type annotations
    custom_mappings : dict, optional
        Custom annotation -> CL ID mappings
    verbose : bool
        Print progress and unmapped terms

    Returns
    -------
    dict
        {annotation: cl_id or None}
    """
    mapping = {}
    unmapped = []

    for ann in annotations:
        cl_id = map_annotation_to_ontology(ann, custom_mappings)
        mapping[ann] = cl_id
        if cl_id is None:
            unmapped.append(ann)
        time.sleep(0.1)  # Rate limiting for API

    if verbose:
        print(f"Mapped {len(annotations) - len(unmapped)}/{len(annotations)} annotations")
        if unmapped:
            print(f"Unmapped: {unmapped}")

    return mapping


# =============================================================================
# REFINED GROUPING FUNCTIONS
# =============================================================================

def define_refined_groups(annotations: List[str],
                          grouping_rules: Dict[str, List[str]]) -> Dict[str, str]:
    """
    Create a mapping from annotations to refined groups.

    Parameters
    ----------
    annotations : list of str
        List of original annotation labels
    grouping_rules : dict
        {refined_group: [list of original annotations]}

    Returns
    -------
    dict
        {original_annotation: refined_group}
    """
    mapping = {}

    for group_name, members in grouping_rules.items():
        for member in members:
            if member in annotations:
                mapping[member] = group_name

    # Check for unmapped annotations
    unmapped = [a for a in annotations if a not in mapping]
    if unmapped:
        print(f"Warning: {len(unmapped)} unmapped annotations: {unmapped[:5]}...")

    return mapping


def get_group_at_level(hierarchy: Dict, target_level: int,
                       current_level: int = 0) -> Dict[str, List[str]]:
    """
    Extract groups at a specific level of the hierarchy.

    Parameters
    ----------
    hierarchy : dict
        Nested hierarchy dictionary
    target_level : int
        Level to extract (0 = root, 1 = first children, etc.)
    current_level : int
        Current recursion level

    Returns
    -------
    dict
        {group_name: [leaf_descendants]}
    """
    groups = {}

    def get_leaves(node: Dict) -> List[str]:
        """Recursively get all leaf names under a node."""
        leaves = []
        if "children" in node:
            for child_name, child_node in node["children"].items():
                leaves.extend(get_leaves(child_node))
        else:
            return [node.get("label", "unknown")]
        return leaves

    def traverse(node: Dict, name: str, level: int):
        if level == target_level:
            # This is our target level - collect all leaves under it
            if "children" in node:
                leaves = []
                for child_name, child_node in node["children"].items():
                    leaves.extend(get_leaves(child_node) if "children" in child_node else [child_name])
                groups[name] = leaves
            else:
                groups[name] = [name]
        elif "children" in node:
            for child_name, child_node in node["children"].items():
                traverse(child_node, child_name, level + 1)

    for name, node in hierarchy.items():
        traverse(node, name, current_level)

    return groups


# =============================================================================
# LYMPH NODE SPECIFIC GROUPINGS
# =============================================================================

def get_lymph_node_refined_groups() -> Dict[str, List[str]]:
    """
    Get predefined refined groupings for lymph node cell types.
    Based on biological function and spatial organization.

    Returns
    -------
    dict
        {refined_group: [original_annotations]}
    """
    return {
        # T cell subsets
        "Naive_CD4_T": ["Naive"],
        "Memory_CD4_T": ["CM Pre-non-Tfh", "CM PreTfh", "T-Trans-Mem", "T-Eff-Mem"],
        "Treg": ["Eff-Tregs", "Eff-Tregs-IL32"],
        "Tfr": ["Tfr"],
        "Tfh_GC": ["GC-Tfh-SAP", "GC-Tfh-OX40", "Tfh-LZ-GC"],
        "Tfh_border": ["Tfh T:B border", "Tfh-Mem"],
        "CD8_T": ["Naive CD8 T", "RM CD8 T", "SCM CD8 T", "cycling T"],

        # B cell subsets
        "Naive_B": ["NBC", "NBC CD229+", "NBC early activation", "NBC IFN-activated"],
        "Pre_GC_B": ["Early GC-commited NBC", "GC-commited NBC", "preGC"],
        "GC_DZ": ["DZ early Sphase", "DZ late Sphase", "DZ cell cycle exit", "DZ non proliferative"],
        "GC_LZ": ["LZ_DZ transition", "DZ_LZ transition", "LZ_DZ reentry commitment",
                  "PC committed Light Zone GCBC"],
        "Memory_B_ncs": ["ncsMBC"],
        "Memory_B_cs": ["csMBC", "MBC FCRL5+", "MBC derived PC precursor",
                        "Reactivated proliferative MBCs"],
        "Plasma": ["IgG+ PC precursor", "preMature IgG+ PC", "MBC derived IgG+ PC", "PB"],

        # Myeloid
        "Monocyte_Mac": ["Monocytes", "C1Q Slan-like", "ITGAX Slan-like"],
        "cDC": ["DC1 mature", "DC2", "DC5", "aDC1"],
        "pDC": ["PDC"],

        # Other
        "NK": ["CD16-CD56+ NK", "CD16+CD56- NK"],
        "Mast": ["Mast"],
    }


def get_lymph_node_coarse_groups() -> Dict[str, List[str]]:
    """
    Get coarser groupings (major lineages) for lymph node cells.

    Returns
    -------
    dict
        {coarse_group: [refined_groups]}
    """
    return {
        "CD4_T": ["Naive_CD4_T", "Memory_CD4_T", "Treg", "Tfr", "Tfh_GC", "Tfh_border"],
        "CD8_T": ["CD8_T"],
        "B_cell": ["Naive_B", "Pre_GC_B", "GC_DZ", "GC_LZ", "Memory_B_ncs", "Memory_B_cs", "Plasma"],
        "Myeloid": ["Monocyte_Mac", "cDC", "pDC"],
        "Other": ["NK", "Mast"],
    }


# =============================================================================
# HIERARCHY VISUALIZATION DATA
# =============================================================================

def build_hierarchy_for_plotting(refined_groups: Dict[str, List[str]],
                                  coarse_groups: Dict[str, List[str]],
                                  colors: Dict[str, str] = None) -> Dict:
    """
    Build a hierarchy dictionary suitable for plotting.

    Parameters
    ----------
    refined_groups : dict
        {refined_group: [original_annotations]}
    coarse_groups : dict
        {coarse_group: [refined_groups]}
    colors : dict, optional
        {group_name: color}

    Returns
    -------
    dict
        Nested dict for plot_cell_lineage_tree()
    """
    default_colors = {
        "CD4_T": "#E41A1C",
        "CD8_T": "#FF7F00",
        "B_cell": "#377EB8",
        "Myeloid": "#4DAF4A",
        "Other": "#984EA3",
    }

    if colors is None:
        colors = default_colors

    hierarchy = {}

    for coarse_name, refined_list in coarse_groups.items():
        color = colors.get(coarse_name, "#999999")
        children = {}

        for refined_name in refined_list:
            # Generate slightly different shade for children
            children[refined_name] = color

        hierarchy[coarse_name] = {
            "color": color,
            "children": children,
        }

    return hierarchy


def print_hierarchy(hierarchy: Dict, indent: int = 0):
    """
    Pretty print a hierarchy dictionary.

    Parameters
    ----------
    hierarchy : dict
        Nested hierarchy
    indent : int
        Current indentation level
    """
    for name, node in hierarchy.items():
        prefix = "  " * indent
        if isinstance(node, dict):
            cl_id = node.get("cl_id", "")
            print(f"{prefix}{name} ({cl_id})")
            if "children" in node:
                print_hierarchy(node["children"], indent + 1)
        else:
            print(f"{prefix}{name}")
