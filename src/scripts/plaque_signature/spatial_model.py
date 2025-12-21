import numpy as np
from sklearn.metrics import r2_score

from .load_data import load_mouse_csv
from .preprocess import ModalityScaler
from .train_signature import train_hgb


def train_spatial_only(data_mouse, train_idx, test_idx):
    """
    Train a spatial-only HGB regression model on a single reference mouse.

    Returns:
        model, scaler, y_test, y_pred_test
    """

    blocks_train = {
        "genes": None,
        "morph": None,
        "cluster": None,
        "spatial": data_mouse["spatial"].iloc[train_idx],
    }

    blocks_test = {
        "genes": None,
        "morph": None,
        "cluster": None,
        "spatial": data_mouse["spatial"].iloc[test_idx],
    }

    scaler = ModalityScaler()
    scaler.fit(blocks_train)

    X_train = scaler.transform(blocks_train).astype("float32")
    X_test = scaler.transform(blocks_test).astype("float32")

    y_train = data_mouse["target"].iloc[train_idx].values.astype("float32")
    y_test = data_mouse["target"].iloc[test_idx].values.astype("float32")

    model = train_hgb(X_train, y_train)
    y_pred_test = model.predict(X_test)

    print("HGB Spatial-only R²:", r2_score(y_test, y_pred_test))

    return model, scaler, y_test, y_pred_test


def predict_spatial_only(model, scaler, mouse_paths):
    """
    Apply a trained spatial-only model to all mice.

    Parameters:
        model : trained HGB model
        scaler : spatial-only ModalityScaler
        mouse_paths : dict, name → csv path

    Returns:
        dict {mouse_name: {"df": df_mouse, "pred": predictions}}
    """
    all_preds = {}

    for name, path in mouse_paths.items():
        print(f"Processing {name}...")

        data_mouse = load_mouse_csv(path)
        df_mouse = data_mouse["df"]

        blocks_mouse = {
            "genes": None,
            "morph": None,
            "cluster": None,
            "spatial": data_mouse["spatial"],
        }

        X_mouse = scaler.transform(blocks_mouse).astype("float32")
        preds = model.predict(X_mouse)

        all_preds[name] = {"df": df_mouse, "pred": preds}

    print("Done.")
    return all_preds
