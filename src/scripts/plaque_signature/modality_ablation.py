import numpy as np
import pandas as pd
from .preprocess import ModalityScaler
from .train_signature import run_ablation_experiments


MODALITY_SETS = [
    ["genes"],
    ["morph"],
    ["spatial"],
    ["cluster"],
    ["genes", "morph"],
    ["genes", "cluster"],
    ["genes", "spatial"],
    ["genes", "morph", "cluster"],
    ["genes", "morph", "spatial", "cluster"],  # full model
]


def run_modality_ablation(data, train_idx, test_idx, modality_sets=MODALITY_SETS):
    """
    Runs ablation experiments over different combinations of modalities.

    Parameters
    ----------
    data : dict
        Output from load_mouse_csv().
    train_idx, test_idx : array-like
        Train/test splits on data["df"].
    modality_sets : list of list of str
        Which combinations of modalities to test.

    Returns
    -------
    DataFrame
        Results of the ablation across all modality sets and models.
    """
    
    all_results = []

    for mods in modality_sets:
        print(f"\n=== Running modalities: {mods} ===")

        # Prepare modality blocks
        blocks_train = {}
        blocks_test = {}

        for m in ["genes", "morph", "spatial", "cluster"]:
            if m in mods:
                blocks_train[m] = data[m].iloc[train_idx]
                blocks_test[m] = data[m].iloc[test_idx]
            else:
                blocks_train[m] = None
                blocks_test[m] = None

        # Scaling
        scaler = ModalityScaler()
        scaler.fit(blocks_train)

        X_train = scaler.transform(blocks_train).astype(np.float32)
        X_test = scaler.transform(blocks_test).astype(np.float32)

        y_train = data["target"].iloc[train_idx].values.astype(np.float32)
        y_test = data["target"].iloc[test_idx].values.astype(np.float32)

        # Run ablations (7 models)
        df = run_ablation_experiments(X_train, y_train, X_test, y_test)
        df["Modalities"] = ", ".join(mods)

        all_results.append(df)

    return pd.concat(all_results, ignore_index=True)
