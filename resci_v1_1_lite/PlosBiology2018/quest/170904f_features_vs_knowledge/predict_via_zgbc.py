import argparse
import os
import random
import time

import numpy as np
import pandas as pd

from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.pipeline import Pipeline

from joblib import Parallel, delayed
import multiprocessing


def _make_and_save_prediction(j, params):

    orig_index_values = params['orig_index_values']
    rows = params['rows']
    fraction_to_train = params['fraction_to_train']
    features = params['features'].astype(float)
    target = params['target'].astype(float)
    estimators = params['estimators']
    predictor_name = params['predictor_name']
    i = params['i']   # input dataset

    df_out = pd.DataFrame(
        index=orig_index_values,
        columns=['predicted'])
    nums = [x for x in range(rows)]
    random.shuffle(nums)
    half = int(np.floor(rows * fraction_to_train))
    first = nums[:half]
    last = nums[half:]

    clf = Pipeline([
        ('zscore', StandardScaler()),
        ('classification', GradientBoostingClassifier(
            n_estimators=estimators))
    ])
    clf.fit(features[first, :], target[first])
    tp = clf.predict(features[last, :])

    df_out.iloc[last, 0] = tp
    df_out = df_out.dropna()

    p_out = os.path.join(
        i,
        predictor_name,
        'prediction_{}_{}.csv.gz'.format(
            str(int(round(time.time() * 1000))),
            str(j))
    )
    _ensure_presence_of_directory(p_out)
    df_out.to_csv(p_out, compression='gzip')


def predict_via_zgbc(i, p, e, n):
    """
    Prediction via zgbc (Z-score, Gradient-Boosting-Classifier)

    Input:
        i       str: path to dataset; e.g.:
                    170827_human_ranked_biophysics_experiments_log_attention
        p       int: percentiles to use for training
                    recommended: 90
        n       int: number of iterations (sequential independent runs)
                    e.g.: 100
        e       int: estimators for gradient-boost-regression
                    processing time ~linearly with e
                    while 100 is default, higher values can increse
                    prediciton further, if many features
                    recommended: 300 unless many more features
    """

    p = int(np.round(p))   # percentiles to use to train
    estimators = int(np.round(e))
    n = int(n)

    predictor_name = 'zgbc_p{}_e{}'.format(p, estimators)
    index_label = 'gene_ncbi'
    target_label = 'prediction_target'

    df_features = pd.read_csv(
        os.path.join(
            i,
            'input',
            'features.csv.gz'),
        index_col=index_label).replace(
            'False', 0).replace('True', 1)

    df_target = pd.read_csv(
        os.path.join(
            i,
            'input',
            'target.csv.gz'),
        index_col=index_label).replace(
            'False', 0).replace('True', 1)

    if not all(df_features.index == df_target.index):
        raise EnvironmentError("Features and training indicies don't match")
    else:
        params = dict()
        params['features'] = df_features.values
        params['target'] = df_target.loc[:, target_label].values
        params['fraction_to_train'] = p / 100
        params['rows'] = params['features'].shape[0]
        params['orig_index_values'] = df_features.index
        params['estimators'] = estimators
        params['predictor_name'] = predictor_name
        params['i'] = i

        num_cores = multiprocessing.cpu_count()
        Parallel(n_jobs=num_cores)(delayed(_make_and_save_prediction)
                                   (j, params) for j in np.arange(n))


# def _add_date_in_milliseconds_and_iteration_(p, j):
#     """
#     Adds date in milliseconds, and iteration number; the latter
#     hopefully avoids some writing mistakes happening at
#     same node

#     Input:
#         p       str: e.g.: path to file
#         j       int interation
#     """
#     [fo, fn] = os.path.split(p)
#     [fb, ext] = os.path.splitext(fn)
#     dt = str(int(round(time.time() * 1000)))
#     p = os.path.join(fo, fb + '_' + dt + '_' + str(j) + ext)
#     return p


def _ensure_presence_of_directory(directory_path=None, ):
    '''
    Ensure that the directory of exists. Creates dictionary with cognate
    name in case that directory does not exist. Anticipates that files have
    extensions separated by '.' symbol (can not create directory with . in
    its name); If file does not have an extension separated by '.' a folder
    will with its filname will be created, a behavior that can be avoided
    by calling os.path.dirname prior this function.

    Input:
        directory_path      str; Name of a directory or the full path of a file
    '''
    if directory_path is None:
        raise ValueError('No input specfied for ensure_presence_of_directory')

    directory_path_n, ext = os.path.split(directory_path)

    if '.' in ext:
        directory_path = directory_path_n

    if not os.path.exists(directory_path):
        os.makedirs(directory_path)


def main():
    pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Runs zgbc for a given prediction set.')
    parser.add_argument('-i', type=str,
                        help='path/to/prediction/base')
    parser.add_argument('-n', type=int,
                        help='number/of/iterations')
    parser.add_argument('-p', type=int,
                        help='percentiles/to/train')
    parser.add_argument('-e', type=int,
                        help='estimators/for/model')

    args = parser.parse_args()

    predict_via_zgbc(
        os.path.abspath(args.i),    # input directory
        args.p,                     # percentiles to train
        args.e,                     # estimators
        args.n,                     # number of iterations
    )
