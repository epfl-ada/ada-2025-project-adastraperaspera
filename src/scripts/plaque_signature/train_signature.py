"""
train_signature.py

Training and evaluating plaque-distance prediction models with modality ablations.
Models included:
- Ridge Regression
- Lasso Regression
- PLS Regression (n_components=10)
- XGBoost
- LightGBM
- CatBoost
- HistGradientBoostingRegressor (HGB)

Also includes:
- unified evaluator
- signature extraction
- full modality ablation runner
"""

import numpy as np
import pandas as pd
import time

from sklearn.linear_model import Ridge, Lasso, LinearRegression
from sklearn.metrics import r2_score
from sklearn.cross_decomposition import PLSRegression
from sklearn.ensemble import HistGradientBoostingRegressor

# LightGBM
import lightgbm as lgb

# XGBoost
try:
    from xgboost import XGBRegressor
    XGB_AVAILABLE = True
except ImportError:
    XGB_AVAILABLE = False

# CatBoost
try:
    from catboost import CatBoostRegressor
    CATBOOST_AVAILABLE = True
except ImportError:
    CATBOOST_AVAILABLE = False


# ---------------------------------------------------------
#                   MODEL TRAINERS
# ---------------------------------------------------------

def train_linear(X_train, y_train):
    """Ordinary least squares regression."""
    model = LinearRegression()
    model.fit(X_train, y_train)
    return model

def train_ridge(X_train, y_train):
    model = Ridge(alpha=1.0)
    model.fit(X_train, y_train)
    return model


def train_lasso(X_train, y_train, alpha=0.1):
    model = Lasso(alpha=alpha, max_iter=5000)
    model.fit(X_train, y_train)
    return model


def train_pls(X_train, y_train, n_components=10):
    n_features = X_train.shape[1]
    k = min(n_components, n_features)
    pls = PLSRegression(n_components=k)
    pls.fit(X_train, y_train)
    return pls


def train_xgboost(X_train, y_train):
    if not XGB_AVAILABLE:
        raise ImportError("XGBoost is not installed.")

    model = XGBRegressor(
        n_estimators=250,
        learning_rate=0.05,
        max_depth=6,
        subsample=0.8,
        colsample_bytree=0.8,
        reg_lambda=1.0,
        random_state=42    )
    model.fit(X_train, y_train)
    return model


def train_lightgbm(X_train, y_train):
    train_data = lgb.Dataset(X_train, label=y_train)
    params = {
        "objective": "regression",
        "metric": "rmse",
        "learning_rate": 0.05,
        "num_leaves": 31,
        "feature_fraction": 0.8,
        "bagging_fraction": 0.8,
        "bagging_freq": 1,
        "verbose": -1,
    }
    model = lgb.train(params, train_data, num_boost_round=200)
    return model


def train_catboost(X_train, y_train):
    if not CATBOOST_AVAILABLE:
        raise ImportError("CatBoost is not installed.")

    model = CatBoostRegressor(
        depth=6,
        iterations=400,
        learning_rate=0.05,
        loss_function="RMSE",
        random_seed=42,
        verbose=False,
    )
    model.fit(X_train, y_train)
    return model


def train_hgb(X_train, y_train):
    model = HistGradientBoostingRegressor(
        learning_rate=0.05,
        max_depth=6,
        max_iter=400,
        random_state=42,
    )
    model.fit(X_train, y_train)
    return model


# ---------------------------------------------------------
#                   EVALUATION
# ---------------------------------------------------------

def evaluate_model(model, X_train, y_train, X_test, y_test):
    """Return train & test R²."""
    y_pred_train = model.predict(X_train)
    y_pred_test = model.predict(X_test)
    return {
        "train_r2": r2_score(y_train, y_pred_train),
        "test_r2": r2_score(y_test, y_pred_test),
    }


# ---------------------------------------------------------
#                 SIGNATURE EXTRACTION
# ---------------------------------------------------------

def extract_plaque_signature(pls_model, X):
    """Return first latent PLS axis."""
    return pls_model.transform(X)[:, 0]


# ---------------------------------------------------------
#                 FULL ABLATION RUNNER
# ---------------------------------------------------------

def run_ablation_experiments(X_train, y_train, X_test, y_test):
    """
    Run all models on the provided feature matrix for a SINGLE modality set.
    """

    model_trainers = {
        "Ridge": train_ridge,
        "Lasso": train_lasso,
        "PLS": lambda X, y: train_pls(X, y, n_components=10),
        "LightGBM": train_lightgbm,
        "HGB": train_hgb,
    }

    if XGB_AVAILABLE:
        model_trainers["XGBoost"] = train_xgboost
    if CATBOOST_AVAILABLE:
        model_trainers["CatBoost"] = train_catboost

    results = []

    for name, trainer in model_trainers.items():
        try:
            start = time.time()
            model = trainer(X_train, y_train)
            elapsed = time.time() - start

            scores = evaluate_model(model, X_train, y_train, X_test, y_test)
            scores["Model"] = name
            scores["Time"] = elapsed

            results.append(scores)

        except Exception as e:
            results.append({
                "Model": name,
                "train_r2": None,
                "test_r2": None,
                "Time": None,
                "Error": str(e),
            })

    return pd.DataFrame(results)
