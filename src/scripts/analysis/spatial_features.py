"""
spatial_features.py

Functions for computing spatial features related to plaques, including geometry features
like perimeter, major/minor axis lengths, and orientation.
"""

from __future__ import annotations

import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import KDTree
from shapely.geometry import MultiPolygon, Polygon

from src.utils.logging_utils import logger


def compute_plaque_geometry_features(plaques_df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute geometry features for each plaque polygon.

    For each plaque, computes:
    - plaque_area (already exists, but included for completeness)
    - plaque_perimeter
    - plaque_major_axis_length (from minimum rotated rectangle)
    - plaque_minor_axis_length (from minimum rotated rectangle)
    - plaque_orientation (angle of major axis in degrees, 0-180)

    Args:
        plaques_df: DataFrame with columns:
            - 'plaque_id': unique identifier for each plaque
            - 'geometry': Shapely Polygon or MultiPolygon objects

    Returns:
        DataFrame with plaque_id and geometry feature columns:
            - plaque_id
            - plaque_area
            - plaque_perimeter
            - plaque_major_axis_length
            - plaque_minor_axis_length
            - plaque_orientation (degrees, 0-180)
    """
    if "geometry" not in plaques_df.columns:
        raise ValueError("plaques_df must contain a 'geometry' column")
    if "plaque_id" not in plaques_df.columns:
        raise ValueError("plaques_df must contain a 'plaque_id' column")

    logger.info(f"Computing geometry features for {len(plaques_df)} plaques")

    results = []

    for _idx, row in plaques_df.iterrows():
        geom = row["geometry"]
        plaque_id = row["plaque_id"]

        # Handle MultiPolygon by using the largest polygon
        if isinstance(geom, MultiPolygon):
            # Use the polygon with the largest area
            geom = max(geom.geoms, key=lambda p: p.area)

        if not isinstance(geom, Polygon) or not geom.is_valid:
            logger.warning(f"⚠️ Skipping invalid geometry for plaque_id {plaque_id}")
            results.append(
                {
                    "plaque_id": plaque_id,
                    "plaque_area": np.nan,
                    "plaque_perimeter": np.nan,
                    "plaque_major_axis_length": np.nan,
                    "plaque_minor_axis_length": np.nan,
                    "plaque_orientation": np.nan,
                }
            )
            continue

        # Compute area (should already exist, but compute for consistency)
        area = geom.area

        # Compute perimeter
        perimeter = geom.length

        # Compute oriented bounding box (minimum rotated rectangle)
        # This gives us the smallest rectangle that can contain the polygon
        mrr = geom.minimum_rotated_rectangle

        # Get the coordinates of the minimum rotated rectangle
        coords = np.array(mrr.exterior.coords[:-1])  # Exclude last duplicate point

        # Compute edge lengths and store with indices
        edge_data = []
        for i in range(len(coords)):
            p1 = coords[i]
            p2 = coords[(i + 1) % len(coords)]
            edge_length = np.hypot(p2[0] - p1[0], p2[1] - p1[1])
            edge_data.append((edge_length, i))

        # Sort by edge length (descending)
        edge_data.sort(reverse=True, key=lambda x: x[0])

        # Major axis is the longest edge, minor axis is the second longest
        # (For rectangles, these will be the two distinct edge lengths)
        major_axis = edge_data[0][0] if len(edge_data) > 0 else 0.0
        minor_axis = edge_data[1][0] if len(edge_data) > 1 else 0.0

        # Compute orientation: angle of the major axis
        # Use the first longest edge (or any edge of that length)
        major_edge_idx = edge_data[0][1]
        p1 = coords[major_edge_idx]
        p2 = coords[(major_edge_idx + 1) % len(coords)]
        dx = p2[0] - p1[0]
        dy = p2[1] - p1[1]

        # Angle in radians, then convert to degrees
        # atan2 gives angle from positive x-axis, range [-π, π]
        angle_rad = math.atan2(dy, dx)
        angle_deg = math.degrees(angle_rad)

        # Normalize to 0-180 degrees (since orientation is symmetric)
        # We want the angle of the major axis, regardless of direction
        if angle_deg < 0:
            angle_deg += 180
        angle_deg = angle_deg % 180

        results.append(
            {
                "plaque_id": plaque_id,
                "plaque_area": area,
                "plaque_perimeter": perimeter,
                "plaque_major_axis_length": major_axis,
                "plaque_minor_axis_length": minor_axis,
                "plaque_orientation": angle_deg,
            }
        )

    geometry_features_df = pd.DataFrame(results)

    logger.info(
        f"✅ Computed geometry features for {len(geometry_features_df)} plaques. "
        f"Mean area: {geometry_features_df['plaque_area'].mean():.2f}, "
        f"Mean perimeter: {geometry_features_df['plaque_perimeter'].mean():.2f}"
    )

    return geometry_features_df


def add_plaque_geometry_features(
    cells_df: pd.DataFrame,
    plaques_df: pd.DataFrame,
    nearest_plaque_id_col: str = "nearest_plaque_id",
) -> pd.DataFrame:
    """
    Add plaque geometry features to cells based on their nearest plaque.

    For each cell, adds columns:
    - nearest_plaque_area (may already exist)
    - nearest_plaque_perimeter
    - nearest_plaque_major_axis
    - nearest_plaque_minor_axis
    - nearest_plaque_orientation

    Args:
        cells_df: DataFrame with cell data, must contain nearest_plaque_id_col
        plaques_df: DataFrame with plaque data, must contain:
            - 'plaque_id': unique identifier
            - 'geometry': Shapely Polygon or MultiPolygon objects
        nearest_plaque_id_col: Name of the column in cells_df that contains
            the nearest plaque ID (default: 'nearest_plaque_id')

    Returns:
        DataFrame with added geometry feature columns
    """
    if nearest_plaque_id_col not in cells_df.columns:
        raise ValueError(
            f"cells_df must contain column '{nearest_plaque_id_col}'. "
            f"Available columns: {cells_df.columns.tolist()}"
        )

    logger.info(
        f"Adding plaque geometry features to {len(cells_df)} cells "
        f"using {len(plaques_df)} plaques"
    )

    # Compute geometry features for all plaques
    geometry_features_df = compute_plaque_geometry_features(plaques_df)

    # Merge geometry features into cells based on nearest_plaque_id
    # Rename columns to add 'nearest_' prefix
    geometry_features_renamed = geometry_features_df.rename(
        columns={
            "plaque_area": "nearest_plaque_area",
            "plaque_perimeter": "nearest_plaque_perimeter",
            "plaque_major_axis_length": "nearest_plaque_major_axis",
            "plaque_minor_axis_length": "nearest_plaque_minor_axis",
            "plaque_orientation": "nearest_plaque_orientation",
        }
    )

    # Merge with cells_df
    # If some columns already exist (like nearest_plaque_area), we need to handle them
    # by dropping them first or using suffixes
    existing_cols = [
        col
        for col in geometry_features_renamed.columns
        if col in cells_df.columns and col != "plaque_id"
    ]

    # Temporarily drop existing columns to avoid merge conflicts
    cells_df_temp = cells_df.drop(columns=existing_cols) if existing_cols else cells_df

    # Merge the geometry features
    cells_with_features = cells_df_temp.merge(
        geometry_features_renamed,
        left_on=nearest_plaque_id_col,
        right_on="plaque_id",
        how="left",
    )

    # Drop the plaque_id column from the merge (it's redundant with nearest_plaque_id)
    if "plaque_id" in cells_with_features.columns:
        cells_with_features = cells_with_features.drop(columns=["plaque_id"])

    # Restore other columns from original dataframe that weren't geometry features
    for col in cells_df.columns:
        if col not in cells_with_features.columns:
            cells_with_features[col] = cells_df[col]

    # Check how many cells got matched
    n_matched = cells_with_features["nearest_plaque_perimeter"].notna().sum()
    logger.info(
        f"✅ Added geometry features. "
        f"{n_matched}/{len(cells_with_features)} cells matched to plaques "
        f"({100 * n_matched / len(cells_with_features):.1f}%)"
    )

    return cells_with_features


def compute_multi_plaque_proximity_features(
    cells_df: pd.DataFrame,
    plaques_df: pd.DataFrame,
    dist_threshold: float = 30.0,
    x_col: str | None = None,
    y_col: str | None = None,
) -> pd.DataFrame:
    """
    Compute multi-plaque proximity features for each cell.

    For each cell, computes:
    - n_plaques_within_R: count of plaques within radius R (default 30 µm)
    - mean_dist_to_plaques_within_R: mean distance to plaques within radius R

    Uses KDTree on plaque centroids for fast spatial queries.

    Args:
        cells_df: DataFrame with cell coordinates. Must contain either:
            - 'coord_X' and 'coord_Y', or
            - 'x_centroid' and 'y_centroid'
        plaques_df: DataFrame with plaque data, must contain:
            - 'geometry': Shapely Polygon or MultiPolygon objects, or
            - 'centroid_x' and 'centroid_y': plaque centroid coordinates
        dist_threshold: Radius in microns for counting nearby plaques (default: 30.0)
        x_col: Name of x-coordinate column in cells_df (auto-detected if None)
        y_col: Name of y-coordinate column in cells_df (auto-detected if None)

    Returns:
        DataFrame with added columns:
            - n_plaques_within_R (where R is dist_threshold)
            - mean_dist_to_plaques_within_R
    """
    # Auto-detect coordinate columns
    if x_col is None or y_col is None:
        if "coord_X" in cells_df.columns and "coord_Y" in cells_df.columns:
            x_col = "coord_X"
            y_col = "coord_Y"
        elif "x_centroid" in cells_df.columns and "y_centroid" in cells_df.columns:
            x_col = "x_centroid"
            y_col = "y_centroid"
        else:
            raise ValueError(
                "Could not auto-detect coordinate columns. "
                "Expected 'coord_X'/'coord_Y' or 'x_centroid'/'y_centroid'. "
                f"Available columns: {cells_df.columns.tolist()}"
            )

    if x_col not in cells_df.columns or y_col not in cells_df.columns:
        raise ValueError(
            f"cells_df must contain columns '{x_col}' and '{y_col}'. "
            f"Available columns: {cells_df.columns.tolist()}"
        )

    logger.info(
        f"Computing multi-plaque proximity features for {len(cells_df)} cells "
        f"using {len(plaques_df)} plaques with radius={dist_threshold} µm"
    )

    # Extract plaque centroids
    if "centroid_x" in plaques_df.columns and "centroid_y" in plaques_df.columns:
        plaque_centroids = plaques_df[["centroid_x", "centroid_y"]].values
    elif "geometry" in plaques_df.columns:
        # Compute centroids from geometry
        plaque_centroids = np.array(
            [
                [geom.centroid.x, geom.centroid.y]
                for geom in plaques_df["geometry"]
                if geom is not None and hasattr(geom, "centroid")
            ]
        )
        if len(plaque_centroids) != len(plaques_df):
            logger.warning(
                f"⚠️ Could not extract centroids for all plaques. "
                f"Expected {len(plaques_df)}, got {len(plaque_centroids)}"
            )
    else:
        raise ValueError(
            "plaques_df must contain either 'centroid_x'/'centroid_y' columns "
            "or 'geometry' column with Shapely geometries"
        )

    if len(plaque_centroids) == 0:
        logger.warning("⚠️ No plaque centroids found. Setting all proximity features to 0.")
        cells_out = cells_df.copy()
        cells_out[f"n_plaques_within_{dist_threshold:.0f}um"] = 0
        cells_out[f"mean_dist_to_plaques_within_{dist_threshold:.0f}um"] = np.nan
        return cells_out

    # Build KDTree on plaque centroids
    tree = KDTree(plaque_centroids)

    # Extract cell coordinates
    cell_coords = cells_df[[x_col, y_col]].values

    # Query all cells at once for efficiency
    n_plaques_list = []
    mean_dist_list = []

    for cell_xy in cell_coords:
        # Find all plaques within radius
        indices = tree.query_ball_point(cell_xy, dist_threshold)

        if len(indices) == 0:
            n_plaques_list.append(0)
            mean_dist_list.append(np.nan)
        else:
            # Count plaques
            n_plaques_list.append(len(indices))

            # Compute mean distance to plaques within radius
            nearby_centroids = plaque_centroids[indices]
            distances = np.linalg.norm(nearby_centroids - cell_xy, axis=1)
            mean_dist_list.append(float(np.mean(distances)))

    # Create output dataframe
    cells_out = cells_df.copy()
    col_name_count = f"n_plaques_within_{dist_threshold:.0f}um"
    col_name_mean = f"mean_dist_to_plaques_within_{dist_threshold:.0f}um"

    cells_out[col_name_count] = n_plaques_list
    cells_out[col_name_mean] = mean_dist_list

    # Log statistics
    n_zero = (cells_out[col_name_count] == 0).sum()
    n_nonzero = (cells_out[col_name_count] > 0).sum()
    mean_n_plaques = cells_out[col_name_count].mean()
    median_n_plaques = cells_out[col_name_count].median()

    logger.info(
        f"✅ Computed multi-plaque proximity features. "
        f"Cells with 0 plaques: {n_zero} ({100*n_zero/len(cells_out):.1f}%), "
        f"Cells with >0 plaques: {n_nonzero} ({100*n_nonzero/len(cells_out):.1f}%). "
        f"Mean plaques per cell: {mean_n_plaques:.2f}, "
        f"Median: {median_n_plaques:.1f}"
    )

    return cells_out


def plot_plaque_proximity_distribution(
    cells_df: pd.DataFrame,
    dist_threshold: float = 30.0,
    figsize: tuple[int, int] = (10, 6),
    save_path: str | None = None,
) -> None:
    """
    Create diagnostic plots showing the distribution of n_plaques_within_R.

    This helps determine if the chosen radius threshold is informative.
    If almost all cells see 0 plaques (or all see many), the radius may need adjustment.

    Args:
        cells_df: DataFrame with multi-plaque proximity features (from
            compute_multi_plaque_proximity_features)
        dist_threshold: The radius threshold used (for labeling)
        figsize: Figure size (width, height) in inches
        save_path: Optional path to save the figure. If None, displays the plot.
    """
    col_name = f"n_plaques_within_{dist_threshold:.0f}um"

    if col_name not in cells_df.columns:
        raise ValueError(
            f"Column '{col_name}' not found in cells_df. "
            "Run compute_multi_plaque_proximity_features first."
        )

    fig, axes = plt.subplots(1, 2, figsize=figsize)

    # Plot 1: Histogram
    ax1 = axes[0]
    counts = cells_df[col_name].value_counts().sort_index()
    ax1.bar(counts.index, counts.values, alpha=0.7, edgecolor="black")
    ax1.set_xlabel(f"Number of plaques within {dist_threshold:.0f} µm")
    ax1.set_ylabel("Number of cells")
    ax1.set_title(f"Distribution of n_plaques_within_{dist_threshold:.0f}um")
    ax1.grid(True, alpha=0.3)

    # Add statistics text
    mean_val = cells_df[col_name].mean()
    median_val = cells_df[col_name].median()
    n_zero = (cells_df[col_name] == 0).sum()
    pct_zero = 100 * n_zero / len(cells_df)

    stats_text = (
        f"Mean: {mean_val:.2f}\n"
        f"Median: {median_val:.1f}\n"
        f"Cells with 0 plaques: {n_zero} ({pct_zero:.1f}%)"
    )
    ax1.text(
        0.98,
        0.98,
        stats_text,
        transform=ax1.transAxes,
        verticalalignment="top",
        horizontalalignment="right",
        bbox={"boxstyle": "round", "facecolor": "wheat", "alpha": 0.5},
        fontsize=9,
    )

    # Plot 2: Cumulative distribution
    ax2 = axes[1]
    sorted_counts = np.sort(cells_df[col_name].values)
    cumulative = np.arange(1, len(sorted_counts) + 1) / len(sorted_counts)
    ax2.plot(sorted_counts, cumulative, linewidth=2)
    ax2.set_xlabel(f"Number of plaques within {dist_threshold:.0f} µm")
    ax2.set_ylabel("Cumulative proportion of cells")
    ax2.set_title(f"Cumulative distribution of n_plaques_within_{dist_threshold:.0f}um")
    ax2.grid(True, alpha=0.3)
    ax2.axvline(
        median_val,
        color="red",
        linestyle="--",
        alpha=0.7,
        label=f"Median: {median_val:.1f}",
    )
    ax2.legend()

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        logger.info(f"✅ Saved diagnostic plot to {save_path}")
    else:
        plt.show()


def suggest_radius_threshold(
    cells_df: pd.DataFrame,
    distance_col: str = "distance_to_plaque",
    percentile: float = 50.0,
) -> float:
    """
    Suggest a radius threshold based on the median distance to closest plaque.

    This can be used as an alternative to a fixed threshold (e.g., 30 µm).

    Args:
        cells_df: DataFrame with distance to nearest plaque
        distance_col: Name of the column containing distance to nearest plaque
        percentile: Percentile to use for threshold suggestion (default: 50.0 for median)

    Returns:
        Suggested radius threshold in microns
    """
    if distance_col not in cells_df.columns:
        raise ValueError(
            f"Column '{distance_col}' not found in cells_df. "
            f"Available columns: {cells_df.columns.tolist()}"
        )

    # Filter out NaN values
    distances = cells_df[distance_col].dropna()

    if len(distances) == 0:
        logger.warning("⚠️ No valid distances found. Returning default threshold of 30.0")
        return 30.0

    suggested_threshold = float(np.percentile(distances, percentile))

    logger.info(
        f"💡 Suggested radius threshold (based on {percentile}th percentile "
        f"of {distance_col}): {suggested_threshold:.2f} µm"
    )

    return suggested_threshold


def compute_pig_neighborhood_means(
    cells_df: pd.DataFrame,
    pig_genes: list[str],
    k_neighbors: int = 100,
) -> pd.DataFrame:
    """
    Compute neighborhood mean expression of PIG genes for each cell.

    Builds a KDTree on (coord_X, coord_Y).
    For each cell:
      - Finds its k nearest neighbors (excluding itself).
      - For each gene in pig_genes:
          computes the mean expression of that gene among neighbors.

    Adds columns:
      - neigh_mean_{gene} for each gene in pig_genes

    Args:
        cells_df: DataFrame with cell data. Must contain:
            - Either 'coord_X'/'coord_Y' or 'x_centroid'/'y_centroid': cell coordinates
            - Columns for each gene in pig_genes: gene expression values
        pig_genes: List of PIG gene names (e.g., ['Gfap', 'Apoe', 'Cst3', ...])
        k_neighbors: Number of nearest neighbors to consider (default: 100)

    Returns:
        DataFrame with added columns:
            - neigh_mean_{gene} for each gene in pig_genes
    """
    # Auto-detect coordinate columns
    x_col = None
    y_col = None

    if "coord_X" in cells_df.columns and "coord_Y" in cells_df.columns:
        x_col = "coord_X"
        y_col = "coord_Y"
    elif "x_centroid" in cells_df.columns and "y_centroid" in cells_df.columns:
        x_col = "x_centroid"
        y_col = "y_centroid"
    else:
        raise ValueError(
            "Could not auto-detect coordinate columns. "
            "Expected 'coord_X'/'coord_Y' or 'x_centroid'/'y_centroid'. "
            f"Available columns: {cells_df.columns.tolist()[:30]}..."
        )

    # Filter to genes that are actually present in the dataframe
    present_genes = [g for g in pig_genes if g in cells_df.columns]
    missing_genes = [g for g in pig_genes if g not in cells_df.columns]

    if not present_genes:
        raise ValueError(
            f"None of the specified PIG genes are present in cells_df. "
            f"Requested: {pig_genes}, Available columns: {cells_df.columns.tolist()[:20]}..."
        )

    if missing_genes:
        logger.warning(
            f"⚠️ Some PIG genes not found in cells_df: {missing_genes}. "
            f"Computing neighborhood means only for present genes: {present_genes}"
        )

    logger.info(
        f"Computing neighborhood PIG expression means for {len(cells_df)} cells "
        f"using k={k_neighbors} neighbors and {len(present_genes)} PIG genes"
    )

    # Extract coordinates using detected column names
    coords = cells_df[[x_col, y_col]].values

    # Build KDTree
    tree = KDTree(coords)

    # Initialize output dataframe (copy to avoid modifying input)
    cells_out = cells_df.copy()

    # Initialize arrays to store neighborhood means for each gene
    neigh_means = {gene: np.full(len(cells_df), np.nan) for gene in present_genes}

    # Query neighbors for each cell
    # We query k_neighbors + 1 to account for the cell itself
    k_query = min(k_neighbors + 1, len(cells_df))

    for i in range(len(cells_df)):
        # Query k_neighbors + 1 nearest neighbors
        # tree.query returns (distance, index) where index can be scalar or array
        if k_query == 1:
            # Only one result (the cell itself), skip
            continue

        distances, indices = tree.query(coords[i], k=k_query)

        # Remove the first index (the cell itself, distance should be ~0)
        # indices is an array when k > 1
        neighbor_indices = indices[1:] if len(indices) > 1 else []

        if len(neighbor_indices) == 0:
            # No neighbors found, leave as NaN
            continue

        # For each PIG gene, compute mean expression among neighbors
        for gene in present_genes:
            # neighbor_indices are integer positions in the dataframe
            neighbor_expressions = cells_df.iloc[neighbor_indices][gene].values
            # Compute mean, handling NaN values
            mean_expr = np.nanmean(neighbor_expressions)
            neigh_means[gene][i] = mean_expr

    # Add columns to output dataframe
    for gene in present_genes:
        col_name = f"neigh_mean_{gene}"
        cells_out[col_name] = neigh_means[gene]

    # Log statistics
    n_valid = sum(cells_out[f"neigh_mean_{gene}"].notna().sum() for gene in present_genes) / len(
        present_genes
    )
    logger.info(
        f"✅ Computed neighborhood PIG expression means. "
        f"Mean valid values per gene: {n_valid:.0f}/{len(cells_df)} "
        f"({100 * n_valid / len(cells_df):.1f}%)"
    )

    return cells_out
