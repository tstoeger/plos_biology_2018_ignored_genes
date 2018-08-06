import os
import inspect
import pandas as pd


def get_data_path(requested='out'):
    mydir = os.path.dirname(
        os.path.abspath(inspect.getfile(inspect.currentframe())))
    data_dir = os.path.join(mydir, '..', '..', 'data')
    abs_path = os.path.abspath(data_dir)
    if requested == 'out':
        return os.path.join(abs_path, 'out')
    elif requested == 'internal':
        return os.path.join(abs_path, 'internal')


def retreive_general_download_settings():
    '''
    Loads settings file with individual download links,
    Excludes settings for genomes

    Output:
        df      DataFrame with setting file with 'name_of_dataset' as index
    '''

    file_path = _get_settings_file_path('download_links.txt')

    df = pd.read_table(
        file_path,
        sep='\t',
        header=0,
        index_col='name_of_dataset')

    return df


def retreive_nih_exporter_links():
    '''
    Infers the names of the NIH exporter files that should be
    present on NIH's server

    Output:
        download_urls   list of URLs of individual files

    '''

    file_path = _get_settings_file_path('nih_exporter_schemes.csv')
    download_settings = pd.read_csv(file_path)

    download_urls = []
    for r in download_settings.iterrows():

        # retreive settings
        url_scheme = r[1]['URL scheme (? Indicates Increment)']
        increment_start = r[1]['First Value of Increment']
        increment_end = r[1]['Last Value of Increment']

        # check how many characters (?) need to be replaced by increment
        if '?????' in url_scheme:
            raise ValueError("""
                nih_exporter_schemes.csv contains unanticipated format.
                Maximally four ? signs supported.""")
            characters_to_exchange = -1  # safety
        elif '????' in url_scheme:
            characters_to_exchange = 4
        elif '???' in url_scheme:
            characters_to_exchange = 3
        elif '??' in url_scheme:
            characters_to_exchange = 2
        elif '?' in url_scheme:
            characters_to_exchange = 1
        else:
            characters_to_exchange = 0

        # Make list of files that need to be downloaded
        if characters_to_exchange == 0:

            # If not substitution is indicated, include scheme as download
            download_url = url_scheme
            download_urls.append(download_url)

        elif characters_to_exchange > 1:

            # convert (values need to be defined)
            increment_start = int(increment_start)
            increment_end = int(increment_end)

            # add individual increments to download list
            increment_range = range(increment_start, increment_end + 1)
            for i in increment_range:
                pattern_to_replace = '?' * characters_to_exchange
                pattern_of_increment = '{0:0>' + str(
                    characters_to_exchange) + '}'
                pattern_that_replaces = pattern_of_increment.format(i)
                download_url = url_scheme.replace(
                    pattern_to_replace, pattern_that_replaces)
                download_urls.append(download_url)

        else:  # safety
            download_urls.append([])

    return download_urls


def _get_settings_file_path(filebase):
    """
    Gets absolute path, which sits ./../../cfg relative
    to this file

    Input:
        filebase  str; name for file (including extension)
        file_path str; full path to the settings file

    """

    mydir = os.path.dirname(
        os.path.abspath(
            inspect.getfile(inspect.currentframe())))
    file_path = os.path.join(
        mydir, '..', '..', 'cfg', filebase)

    return file_path
