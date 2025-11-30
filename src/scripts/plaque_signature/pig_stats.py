import pandas as pd
from scipy.stats import linregress

def compute_glial_pig_df(pig_z_df, glial_clusters, age_map, disease_map):
    """
    Restrict PIG scores to glial clusters and add age/disease annotations.
    """
    pig_glial = pig_z_df[pig_z_df["cluster"].isin(glial_clusters)].copy()
    pig_glial["age"] = pig_glial["mouse"].map(age_map)
    pig_glial["disease"] = pig_glial["mouse"].map(disease_map)
    return pig_glial


def compute_age_progression(pig_glial_df, glial_clusters):
    rows = []
    for c in glial_clusters:
        df_c = pig_glial_df[pig_glial_df["cluster"] == c]
        tg = df_c[df_c["disease"]=="TG"]
        if tg["age"].nunique() < 2:
            continue
        slope, intercept, r, p, stderr = linregress(tg["age"], tg["mean_pig"])
        rows.append({"cluster": c, "slope": slope, "r": r, "p": p})
    return pd.DataFrame(rows)


def compute_disease_effect(pig_glial_df, glial_clusters):
    rows = []
    for c in glial_clusters:
        df = pig_glial_df[pig_glial_df["cluster"] == c]
        tg_mean = df[df["disease"]=="TG"]["mean_pig"].mean()
        wt_mean = df[df["disease"]=="WT"]["mean_pig"].mean()
        rows.append({
            "cluster": c,
            "TG_minus_WT": tg_mean - wt_mean,
            "TG_mean": tg_mean,
            "WT_mean": wt_mean
        })
    return pd.DataFrame(rows)


def build_summary_table(disease_effect_df, age_df):
    summary = disease_effect_df.merge(age_df, on="cluster", how="left")
    summary["abs_slope"] = summary["slope"].abs()
    return summary.set_index("cluster")
