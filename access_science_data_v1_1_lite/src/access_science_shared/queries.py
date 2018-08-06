import glob
import os
import sqlite3
import tempfile
import time

import numpy as np
import pandas as pd


def sql_w_cache(sql, p_database, refresh=False):
    """
    Performs an sql query, where the result is cached locally
    on disc within the systems's temporary folder for up
    to one week.

    Thus this function may not be used for queries where
    the database is undergoing changes.

    Input:
        sql         SQL string
        p_database  str path do to database
        refresh     default: False
                    (set to True to reload from database)

    Output:
        df_data     dataframe, as queried
    """

    # Internal Settings
    cache_format = 'msgpack'
    maximal_age_in_sec = 60 * 60 * 24 * 7

    # Set general environment of function
    current_time = time.time()

    p_cache_folder = os.path.join(
        tempfile.gettempdir(),
        'access_science/shared/sql_cache')

    if not os.path.exists(p_cache_folder):
        os.makedirs(p_cache_folder)

    # Helper function
    def _filename_from_query_id(query_id):
        fn = os.path.join(
            p_cache_folder, 'cached_query_{}.{}'.format(
                query_id, cache_format))
        return fn

    # Get information from register
    p_register = os.path.join(
        p_cache_folder,
        'register.csv')

    if os.path.exists(p_register):

        # read register
        df_register = pd.read_csv(p_register)

        # remove outdated files
        is_outdated = current_time - df_register.loc[
            :, 'timestamp'] > maximal_age_in_sec
        for ts in df_register.loc[is_outdated, 'query_id'].values:
            fn = _filename_from_query_id(int(ts))
            os.remove(fn)
        df_register = df_register.loc[~is_outdated, :]

        # Check for presence
        is_present = (df_register.loc[:, 'sql_query'] == sql
                      ) & (df_register.loc[:, 'db_location'] == p_database)

        if refresh:

            # Initiate Settings for data retreival
            if df_register.shape[0] > 0:
                current_query_id = df_register.iloc[
                    -1, :]['query_id'] + 1
                current_query_id = int(current_query_id)
            else:  # if old entries had been removed
                current_query_id = 1
            query_database = True

            # Update Register
            df_register = df_register.loc[~is_present, :]

        else:
            pr = np.count_nonzero(is_present)
            if pr == 1:
                # Initiate Settings for data retreival
                current_query_id = df_register.loc[
                    is_present, 'query_id']
                current_query_id = int(current_query_id)
                query_database = False
            elif pr == 0:
                if df_register.shape[0] > 0:
                    current_query_id = df_register.iloc[
                        -1, :]['query_id'] + 1
                    current_query_id = int(current_query_id)
                else:  # if old entries had been removed
                    current_query_id = 1
                query_database = True
            elif pr > 1:
                raise EnvironmentError(
                    'There appear multiple copies of' +
                    'prior entry. This indicates logical ' +
                    'mistake in underlying code' +
                    'Check with author (Thomas Stoeger)')

    else:  # no register

        # Safety: remove existing files
        pattern = os.path.join(
            p_cache_folder, '*.{}'.format(cache_format))
        existing_files = glob.glob(pattern)
        for fi in existing_files:
            os.remove(fi)

        # Initiate register
        df_register = pd.DataFrame(
            columns=[
                'query_id',
                'sql_query',
                'db_location',
                'timestamp'])

        # Initiate Settings for data retreival
        current_query_id = 1
        query_database = True

    # Adress cache or database
    p_cache = _filename_from_query_id(
        current_query_id)

    if query_database:
        con = sqlite3.connect(p_database)
        df_data = pd.read_sql_query(sql, con)
        con.close()
        df_data.to_msgpack(p_cache)

        try:
            df_register = df_register.append({
                'query_id': current_query_id,
                'sql_query': sql,
                'db_location': p_database,
                'timestamp': current_time},
                ignore_index=True)

        except Exception:
            os.remove(p_cache)
            raise EnvironmentError(
                'Could not update cache register.\n' +
                'Removing cachged file from cache.')

    else:  # Get from Cache
        if os.path.exists(p_cache):
            df_data = pd.read_msgpack(p_cache)
        else:
            print('''
                Cache appears corrupted. Will delete everything in chache.
                Please rerun.''')

            print(p_cache)
            import shutil
            shutil.rmtree(p_cache_folder)

    df_register.to_csv(p_register, index=False)

    return df_data
