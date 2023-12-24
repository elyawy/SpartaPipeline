import numpy as np
from sklearn.linear_model import SGDRegressor, Lasso
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.model_selection import GridSearchCV
from lightgbm import LGBMRegressor


def performance_metric(pearson_correlations):
    sorted_correlations = sorted(pearson_correlations)
    print(sorted_correlations)
    print(len(pearson_correlations)//10)
    percentile_pearson = sorted_correlations[len(pearson_correlations)//20]
    return percentile_pearson

REG_TYPES = ["lasso", "lgbm"]#, "sgd"]


def get_regression_methods(reg_type):
    if reg_type=="lasso":
        reg = Lasso()
        parameters = {'alpha':np.logspace(-7,4,20)}
        clf_lassocv = GridSearchCV(estimator = reg,
                                    param_grid = parameters, cv=3,
                                    scoring = 'neg_mean_squared_error')
        reg_pipline = Pipeline([('scaler', StandardScaler()),('reg', clf_lassocv)])
        return reg_pipline
    elif reg_type=="lgbm":
        reg_pipline = Pipeline([('reg', LGBMRegressor())])
        return reg_pipline
    elif reg_type=="sgd":
        return Pipeline([('scaler', StandardScaler()),('reg', SGDRegressor())])
