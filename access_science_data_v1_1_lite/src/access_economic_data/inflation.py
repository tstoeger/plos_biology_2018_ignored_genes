import os
import pandas as pd

from access_science_shared import inout


def consumer_price_adjustment_factor(ref_year):
    """
    Will load adjustment factor for consumer price, relative
    to ref_year; Source of data: urbaan consumer price index from
    usinflationcalculator.com

    Input:
        ref_year    int, year that shall be used as a reference

    Output:
        df_cpi

    """

    p_cpi = inout.get_path(
        'rbusa',
        os.path.join(
            'economics',
            'cpi',
            'usinflationcalculator.csv'))

    df_cpi = pd.read_csv(p_cpi, skiprows=1, usecols=['Year', 'Avg'])
    df_cpi = df_cpi.dropna()
    df_cpi = df_cpi.set_index('Year')
    df_cpi['adjustment_factor'] = df_cpi.loc[
        ref_year, 'Avg'] / df_cpi['Avg']
    df_cpi = df_cpi.reset_index().drop('Avg', axis=1)

    return df_cpi
