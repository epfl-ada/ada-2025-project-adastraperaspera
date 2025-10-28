import pandas as pd
import numpy as np
from sklearn.linear_model import LassoCV, ElasticNetCV
from sklearn.ensemble import RandomForestRegressor
from sklearn.neural_network import MLPRegressor
from xgboost import XGBRegressor

def run_lasso(X_train, X_test, y_train, y_test, gene_cols):
    model = LassoCV(cv=5, n_jobs=-1, random_state=42)
    model.fit(X_train, y_train)
    results = {
        "model": "LASSO",
        "train_r2": model.score(X_train, y_train),
        "test_r2": model.score(X_test, y_test),
        "alpha": model.alpha_
    }
    imp = pd.DataFrame({"gene": gene_cols, "importance": np.abs(model.coef_)})
    return model, results, imp.sort_values("importance", ascending=False)

def run_elasticnet(X_train, X_test, y_train, y_test, gene_cols):
    model = ElasticNetCV(cv=5, n_jobs=-1, random_state=42, l1_ratio=[0.1,0.5,0.9])
    model.fit(X_train, y_train)
    results = {
        "model": "ElasticNet",
        "train_r2": model.score(X_train, y_train),
        "test_r2": model.score(X_test, y_test),
        "alpha": model.alpha_,
        "l1_ratio": model.l1_ratio_
    }
    imp = pd.DataFrame({"gene": gene_cols, "importance": np.abs(model.coef_)})
    return model, results, imp.sort_values("importance", ascending=False)

def run_rf(X_train, X_test, y_train, y_test, gene_cols):
    model = RandomForestRegressor(
        n_estimators=300, max_depth=15, min_samples_leaf=100,
        n_jobs=-1, random_state=42)
    model.fit(X_train, y_train)
    results = {
        "model": "RandomForest",
        "train_r2": model.score(X_train, y_train),
        "test_r2": model.score(X_test, y_test)
    }
    imp = pd.DataFrame({"gene": gene_cols, "importance": model.feature_importances_})
    return model, results, imp.sort_values("importance", ascending=False)

def run_xgb(X_train, X_test, y_train, y_test, gene_cols):
    model = XGBRegressor(
        n_estimators=500, learning_rate=0.05, max_depth=6,
        subsample=0.8, colsample_bytree=0.8, reg_lambda=1.0,
        n_jobs=-1, random_state=42, verbosity=0)
    model.fit(X_train, y_train)
    results = {
        "model": "XGBoost",
        "train_r2": model.score(X_train, y_train),
        "test_r2": model.score(X_test, y_test)
    }
    imp = pd.DataFrame({"gene": gene_cols, "importance": model.feature_importances_})
    return model, results, imp.sort_values("importance", ascending=False)


