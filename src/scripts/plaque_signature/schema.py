"""
schema.py

Defines canonical column groups for plaque signature modeling.
"""

# Columns always present across all mice
SPATIAL_COLS = ["x_centroid", "y_centroid"]
MORPH_COLS = ["cell_area", "nucleus_area"]
CLUSTER_COL = "cluster_leiden"

# Plaque-related columns (present only in Tg 17.9-month)
PLAQUE_COLS = [
    "distance_to_plaque",
    "nearest_plaque_center_dist",
    "inside_any_plaque",
    "nearest_plaque_id",
    "nearest_plaque_area",
]

TARGET_COL = "distance_to_plaque"

# Function to detect gene columns given a dataframe
def get_gene_columns(df):
    exclude = set(
        ["cell_id"]
        + SPATIAL_COLS
        + MORPH_COLS
        + [CLUSTER_COL]
        + PLAQUE_COLS
        + [
            "transcript_counts",
            "control_probe_counts",
            "control_codeword_counts",
            "unassigned_codeword_counts",
            "total_counts",
        ]
    )
    return [c for c in df.columns if c not in exclude]
