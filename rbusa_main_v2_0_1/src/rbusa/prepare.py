import os
import re
import sqlite3

import numpy as np
import pandas as pd

from shutil import copytree

import rbusa.inout as io

from rbusa.settings import retreive_nih_exporter_links
from rbusa import downloader


def consumer_price_index():
    """
    Gets the consoumer price index from rbusa_manual
    """

    p_ext = 'economics/cpi/'

    p_in = io.get_rbusa_manual_data_path(p_ext)
    p_out = io.get_output_path(p_ext)

    copytree(p_in, p_out)


def diego_fregolente_on_nih():
    """
    Information on gender of NIH supported authors,
    as computattionally inferred by Diego Fregolente
    """

    p_ext = 'gender/diego_fregolente_on_nih/'

    p_in = io.get_rbusa_manual_data_path(p_ext)
    p_out = io.get_output_path(p_ext)

    copytree(p_in, p_out)


def nih_exporter():
    """
    Downloads all NIH funding information; Note that NIH yearly
    replaces monthly updates by annual data, thus requiring
    updates in the /cfg/nih_exporter_schemes.csv file
    """

    subdir_name = 'nih_exporter'
    out_p_main = io.get_output_path(subdir_name)
    io.ensure_presence_of_directory(out_p_main)

    internal_p_main = io.get_internal_path(subdir_name)
    io.ensure_presence_of_directory(internal_p_main)

    download_urls = retreive_nih_exporter_links()

    # Download Files
    for url_of_download in download_urls:

        # Direct files that need to be reconstitued
        # with patches into an intenral directory
        if '_PRJ_' in url_of_download:
            out_base = internal_p_main
        elif '_PRJFUNDING_' in url_of_download:
            out_base = internal_p_main
        elif '_DUNS_' in url_of_download:
            out_base = internal_p_main
        else:
            out_base = out_p_main

        if 'ALL.zip' in url_of_download:
            out_p = out_base
        else:
            # Organize by subfolders uniting one data category
            # (e.g: abastracts, or linkage to medline)
            category = re.findall(
                '.*RePORTER_([A-Z]*)', url_of_download)[0]
            category = category.lower()
            out_p = os.path.join(out_base, category)
            io.ensure_presence_of_directory(out_p)

        # Download and extract data and save as gzip
        _, filename = os.path.split(url_of_download)
        path_downloaded_file = os.path.join(out_p, filename)
        downloader.download(url_of_download, path_downloaded_file)
        downloader.unzip(path_downloaded_file, delete_zip=True)

        anticipated_filename_csv = filename.replace('.zip', '.csv')

        potential_subfolder = os.path.splitext(filename)[0]  # For recent

        p_legacy = os.path.join(
            out_p,
            anticipated_filename_csv)    # till 2016
        p_recent = os.path.join(
            out_p,
            potential_subfolder,
            anticipated_filename_csv)
        if os.path.exists(p_legacy):
            anticipated_csv = p_legacy
        elif os.path.exists(p_recent):
            anticipated_csv = p_recent

        p_final_location = os.path.join(
            out_p, anticipated_filename_csv + '.gz')

        df = pd.read_csv(anticipated_csv, low_memory=False, encoding='latin1')
        df.to_csv(p_final_location, index=False, compression='gzip')
        os.remove(anticipated_csv)

        p = os.path.join(out_p, potential_subfolder)

        if os.path.exists(p):
            os.rmdir(p)

    # reconstitute project files from year specific patches
    era = range(1985, 2000)
    for year in era:
        p_patch = os.path.join(
            internal_p_main,
            'prjfunding',
            'RePORTER_PRJFUNDING_C_FY{}.csv.gz'.format(year))

        p_main = os.path.join(
            internal_p_main,
            'prj',
            'RePORTER_PRJ_C_FY{}.csv.gz'.format(year))

        df_patch = pd.read_csv(
            p_patch, low_memory=False).set_index(
                'APPLICATION_ID', verify_integrity=True)
        df_main = pd.read_csv(
            p_main, low_memory=False).set_index(
                'APPLICATION_ID', verify_integrity=True)

        df = df_patch.combine_first(df_main)   # patch is prioritized
        df = df.reset_index()

        p_final = os.path.join(
            out_p_main,
            'prj',
            'RePORTER_PRJFUNDING_C_FY{}_REVIEWED.csv.gz'.format(year))
        io.ensure_presence_of_directory(p_final)
        df.to_csv(p_final, encoding='latin1', compression='gzip', index=False)

    era = range(2000, 2009)
    for year in era:
        p_patch = os.path.join(
            internal_p_main,
            'duns',
            'RePORTER_DUNS_C_FY{}.csv.gz'.format(year))

        p_main = os.path.join(
            internal_p_main,
            'prj',
            'RePORTER_PRJ_C_FY{}.csv.gz'.format(year))

        df_patch = pd.read_csv(
            p_patch, low_memory=False).set_index(
                'APPLICATION_ID', verify_integrity=True)
        df_main = pd.read_csv(
            p_main, low_memory=False).set_index(
                'APPLICATION_ID', verify_integrity=True)

        df = df_patch.combine_first(df_main)   # patch is prioritized
        df = df.reset_index()

        p_final = os.path.join(
            out_p_main,
            'prj',
            'RePORTER_PRJFUNDING_C_FY{}_REVIEWED.csv.gz'.format(year))
        df.to_csv(p_final, encoding='latin1', compression='gzip', index=False)

    era = range(2009, 2017)
    for year in era:
        p_main = os.path.join(
            internal_p_main,
            'prj',
            'RePORTER_PRJ_C_FY{}.csv.gz'.format(year))

        df_main = pd.read_csv(
            p_main, low_memory=False).set_index(
                'APPLICATION_ID', verify_integrity=True)

        df = df_main
        df = df.reset_index()

        p_final = os.path.join(
            out_p_main,
            'prj',
            'RePORTER_PRJFUNDING_C_FY{}_REVIEWED.csv.gz'.format(year))
        df.to_csv(p_final, encoding='latin1', compression='gzip', index=False)


def wos_dais():
    """
    Extracts WOS - DAIS information (as obtained from Diego
    Fregolente, Brian Uzi lab, according to NICO's access to
    Web of Sciences's disambiguation)

    Will split records into smaller batches that can be loaded
    without memory-overload on typical computes (and being able
    to be opened by pandas)

    Will additionally save a partial copy of these indiviual
    batches, that only contains those records that are indirectly
    mapped to an individual gene

    Note that the execution of this function might take 2-3h.
    (as loops and batches are used internally to avoid usage
    of too much memory)

    Output format:
        table with  WOS (ID), DAIS (author ID), and position in paper
    """

    # set output folder
    p_out = io.get_output_path('disambiguation/wos_dais/')
    io.ensure_presence_of_directory(p_out)

    # get WoS Identifiers, mapped to genes
    p_medline = io.get_medline_wos_path('db/medline_wos.sqlite')
    sql = '''
    SELECT DISTINCT
        ut2pmid.ut AS medline_mapped_wos
    FROM ut2pmid
    ORDER BY ut2pmid.ut ASC;
    '''
    con = sqlite3.connect(p_medline)
    df_data = pd.read_sql_query(sql, con)
    con.close()
    wos_mapped_to_medline = df_data['medline_mapped_wos'].values

    # load WOS - DAIS data
    p = io.get_rbusa_manual_data_path(
        'author_disambiguation/170607_paper_IDS_DAIS.txt.gz')
    df = pd.read_table(p)

    # prepare batches
    amount_of_records = df.shape[0]
    batch_size = 200000
    boundaries = np.arange(0, amount_of_records, batch_size)
    boundaries = np.append(boundaries, amount_of_records)

    # Process batch-wise
    amount_of_batches = len(boundaries) - 1
    for j in range(amount_of_batches):
        curr_batch = int(j + 1)
        print('Processing batch ', curr_batch, 'out of ', amount_of_batches)
        curr_start = boundaries[j]
        curr_end = boundaries[j + 1]

        df_batch = df.iloc[curr_start: curr_end, :].copy()

        df_batch['WOS'] = df_batch['ID|DAIS'].str.extract(
            '^WOS:(.*)\|.*', expand=False)
        df_batch['DAIS'] = df_batch['ID|DAIS'].str.extract(
            '^WOS:.*.*\|(.*)', expand=False)
        df_batch = df_batch.drop('ID|DAIS', axis=1)
        df_batch = df_batch.set_index('WOS', verify_integrity=True)

        df_batch = df_batch['DAIS']  # turn to series

        separator = '\,'
        has_multiple_authors = df_batch.str.contains(separator)

        df_single = df_batch.loc[~has_multiple_authors]
        df_single = df_single.reset_index()
        df_single.loc[:, 'position'] = 1

        df_multi = df_batch[has_multiple_authors]
        df_multi = df_multi.str.split(separator).apply(
            pd.Series, 1).stack().reset_index()

        df_multi = df_multi.rename(columns={'level_1': 'position', 0: 'DAIS'})
        df_multi['position'] = df_multi['position'] + 1
        df_out_all = pd.concat([df_single, df_multi], axis=0)

        f = df_out_all.loc[:, 'WOS'].isin(wos_mapped_to_medline)
        df_out_medline = df_out_all.loc[f, :]

        outfile_all = 'wos_dais_all_batch{}'.format(curr_batch)
        outfile_medline = 'wos_dais_gene_mapped_batch_{}'.format(curr_batch)

        def to_disk(df, filebase):
            p = os.path.join(
                p_out,
                filebase + '.csv.gz')
            df.to_csv(p, compression='gzip', index=False)

        to_disk(df_out_all, outfile_all)
        to_disk(df_out_medline, outfile_medline)

        del df_batch
        del df_out_all
        del df_out_medline
        del df_multi
        del df_single


def zeng_2016():
    """
    Downloads gender data from Zeng et al. 2016
    (Version 2 from Fighare), and saves data in a more
    memory efficient way

    https://figshare.com/articles/
    Publication_records_for_top_US_researchers_in_6_fields/3472616

    """

    out_folder = io.get_output_path('gender/zeng2016')
    io.ensure_presence_of_directory(out_folder)

    # Researchers' informatoin
    stored_location = downloader.download_data_set(
        'zeng_researchers',
        folder_contains_dots=False)[0]
    p = os.path.join(
        out_folder,
        'zeng2016_researchers.csv.gz')
    df = pd.read_csv(
        stored_location,
        encoding='ISO-8859-1')
    df.to_csv(p, compression='gzip', index=False)

    # Researchers' publications
    stored_location = downloader.download_data_set(
        'zeng_publications', folder_contains_dots=False)[0]
    p = os.path.join(
        out_folder,
        'zeng2016_publications.csv.gz')
    df = pd.read_csv(
        stored_location,
        encoding='ISO-8859-1')
    df.to_csv(p, compression='gzip', index=False)
