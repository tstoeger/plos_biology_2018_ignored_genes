import os
import pandas as pd

from access_science_shared import inout


def prj(cols_to_use=None):
    """
    Will load NIH project descriptions.
    """

    p_folder = inout.get_path(
        'rbusa',
        os.path.join(
            'nih_exporter',
            'prj'))

    era = range(1985, 2017)
    agg = []

    for year in era:
        p = os.path.join(
            p_folder,
            'RePORTER_PRJFUNDING_C_FY{}_REVIEWED.csv.gz'.format(
                year))

        dfy = pd.read_csv(
            p,
            low_memory=False,
            encoding='latin1',
            usecols=cols_to_use)
        agg.append(dfy)

    df = pd.concat(agg, axis=0)

    df['CORE_PROJECT_NUM'] = df['CORE_PROJECT_NUM'].astype(str)

    return df


def publnk():
    """
    Will retreive the linkage of pubmed_ids to nih
    grants
    """

    p_folder = inout.get_path(
        'rbusa',
        os.path.join(
            'nih_exporter',
            'publnk'))

    era = range(1980, 2017)
    agg = []

    for year in era:
        p = os.path.join(
            p_folder,
            'RePORTER_PUBLNK_C_{}.csv.gz'.format(
                year))

        dfy = pd.read_csv(p, low_memory=False)
        agg.append(dfy)

    df = pd.concat(agg, axis=0)
    df = df.rename(columns={'PMID': 'pubmed_id'})
    df = df.drop_duplicates()   # ~30 duplications

    df['PROJECT_NUMBER'] = df['PROJECT_NUMBER'].astype(str)

    return df
