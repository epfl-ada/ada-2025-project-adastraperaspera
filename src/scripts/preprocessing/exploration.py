from __future__ import annotations

from typing import Any, Literal

import pandas as pd
from pandas import DataFrame, Series


def summarize_missing_by_column(
    df: DataFrame,
    *,
    sort_by: Literal["pct_missing", "n_missing"] = "pct_missing",
    ascending: bool = False,
    include_non_missing: bool = True,
    serialize: bool = False,
) -> DataFrame | list[dict[str, Any]]:
    """
    Summarize missing (NA) values for each column.

    This uses `DataFrame.isna()` to detect missing entries (covering NaN, None, and NaT),
    counts them per column, computes proportions, and returns a table sorted by the
    most missing columns first by default.

    Parameters
    ----------
    df : pandas.DataFrame
        The input DataFrame.
    sort_by : {"pct_missing", "n_missing"}, default "pct_missing"
        Which statistic to sort by.
    ascending : bool, default False
        Sort order (descending by default to surface the most NA-heavy columns).
    include_non_missing : bool, default True
        If True, include the non-missing count as a convenience column.
    serialize : bool, default False
        If True, return a JSON-serializable list of dicts instead of a DataFrame.

    Returns
    -------
    pandas.DataFrame or list[dict]
        Columns:
          - column: str
          - dtype: str
          - n_missing: int
          - pct_missing: float in [0, 1]
          - n_non_missing: int (optional)
        If `serialize=True`, returns `list[dict]` with the same fields.
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError("df must be a pandas DataFrame")

    n_rows = len(df)
    na_mask = df.isna()

    n_missing = na_mask.sum(axis=0).astype("int64")
    pct_missing = (n_missing / n_rows).astype("float64")
    result = pd.DataFrame(
        {
            "column": df.columns,
            "dtype": df.dtypes.astype(str).values,
            "n_missing": n_missing.values,
            "pct_missing": pct_missing.values,
        }
    )

    if include_non_missing:
        result["n_non_missing"] = (n_rows - n_missing).values

    result = result.sort_values(
        by=sort_by, ascending=ascending, kind="mergesort"
    ).reset_index(drop=True)

    return result.to_dict(orient="records") if serialize else result


def summarize_missing_by_row(
    df: DataFrame,
    *,
    ddof: int = 0,
) -> tuple[Series, Series]:
    """
    Compute per-row missing counts and dataset-level aggregates.

    Parameters
    ----------
    df : pandas.DataFrame
        The input DataFrame.
    ddof : int, default 0
        Delta degrees of freedom for the standard deviation of per-row missing counts.
        Use 0 for population std (default) or 1 for sample std.

    Returns
    -------
    missing_per_row : pandas.Series
        Number of missing values in each row (index aligns with `df.index`).
    aggregates : pandas.Series
        Dataset-level statistics over the row-wise missing counts, including:
          - rows, columns
          - mean_n_missing, std_n_missing, min_n_missing, median_n_missing,
            q1_n_missing (25th pct), q3_n_missing (75th pct), max_n_missing
          - mean_pct_missing_per_row
          - pct_rows_no_missing, pct_rows_any_missing, pct_rows_all_missing
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError("df must be a pandas DataFrame")

    n_rows, n_cols = df.shape
    na_mask = df.isna()

    missing_per_row = na_mask.sum(axis=1).astype("int64")
    pct_missing_per_row = (missing_per_row / n_cols).astype("float64")

    aggregates = pd.Series(
        {
            "rows": int(n_rows),
            "columns": int(n_cols),
            "mean_n_missing": float(missing_per_row.mean()),
            "std_n_missing": float(missing_per_row.std(ddof=ddof)),
            "min_n_missing": int(missing_per_row.min()),
            "q1_n_missing": float(missing_per_row.quantile(0.25)),
            "median_n_missing": float(missing_per_row.median()),
            "q3_n_missing": float(missing_per_row.quantile(0.75)),
            "max_n_missing": int(missing_per_row.max()),
            "mean_pct_missing_per_row": float(pct_missing_per_row.mean()),
            "pct_rows_no_missing": float((missing_per_row == 0).mean()),
            "pct_rows_any_missing": float((missing_per_row > 0).mean()),
            "pct_rows_all_missing": float((missing_per_row == n_cols).mean()),
        }
    )

    return missing_per_row, aggregates
