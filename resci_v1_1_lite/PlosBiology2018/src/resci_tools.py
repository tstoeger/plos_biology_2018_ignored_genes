import datetime
import math
import os

import pandas as pd

"""
The scope of this .py are general tools for data handling that are independent
of any of the datasets contained in the analysis of science of science.

ensure_presence_of_directory    ensures that directory is present
export_image                    export image so that editable in illustrator
split_text_to_multiple_rows     pivots records contained in same pandas cell

"""

# Settings
p_material = '~/Dropbox/Work/manuscripts/genes/material'


def ensure_presence_of_directory(directory_path=None, ):
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


def export_image(p, insert_date_time=True):
    """
    Will export a current figure;
        - makes font edit-able in illustrator

    Input:
        insert_date_time    optional; default True: will insert
                                date, hour, minute before file
                                extension
    """

    import matplotlib as mpl
    mpl.rcParams['pdf.fonttype'] = 42  # edit-able in illustrator

    p = os.path.join(
        p_material,
        p)

    if insert_date_time:
        [fo, fn] = os.path.split(p)
        [fb, ext] = os.path.splitext(fn)
        dt = datetime.datetime.today().strftime('%y%m%d_%H%M')
        p = os.path.join(fo, fb + '_' + dt + ext)

    ensure_presence_of_directory(p)

    mpl.pyplot.savefig(p, bbox_inches='tight')


def export_raster_image(p, dpi, insert_date_time=True):
    """
    Will export a current figure with 600 dpi

    Input:
        p       str path to file
        dpi     int dots per inch
        insert_date_time    optional; default True: will insert
                                date, hour, minute before file
                                extension
    """

    import matplotlib as mpl
    mpl.rcParams['pdf.fonttype'] = 42  # edit-able in illustrator

    p = get_material_path(p, insert_date_time)

    mpl.pyplot.savefig(p, dpi=dpi, bbox_inches='tight')


def export_full_frame(p, df, insert_date_time=True, save_index=True):
    """
    Will export a dataframe to materials
        - makes font edit-able in illustrator

    Input:
        p                   subpath within material folder
        df                  dataframe to export
        insert_date_time    optional; default True: will insert
                                date, hour, minute before file
                                extension
        save_index          optional; default True: will also
                                export the index
    """

    p = os.path.join(
        p_material,
        p)

    if p.endswith('.csv.gz'):
        p = p[:-3]
        compress = True
        file_format = 'csv'
    elif p.endswith('.csv'):
        compress = False
        file_format = 'csv'
    elif p.endswith('.xlsx'):
        file_format = 'xlsx'
    else:
        raise EnvironmentError(
            'No support for preseent file type.')

    if insert_date_time:
        [fo, fn] = os.path.split(p)
        [fb, ext] = os.path.splitext(fn)
        dt = datetime.datetime.today().strftime('%y%m%d_%H%M')
        p = os.path.join(fo, fb + '_' + dt + ext)

    ensure_presence_of_directory(p)

    if file_format == 'csv':
        if compress:
            p = p + '.gz'
            df.to_csv(p, compression='gzip', index=save_index)
        else:
            df.to_csv(p, index=save_index)
    elif file_format == 'xlsx':
        df.to_excel(p, index=save_index)


def get_material_path(p, insert_date_time=False):
    """
    Takes extension path p, and makes it as a subfile
    within material folder
    """

    p = os.path.join(
        p_material,
        p)

    if insert_date_time:
        [fo, fn] = os.path.split(p)
        [fb, ext] = os.path.splitext(fn)
        dt = datetime.datetime.today().strftime('%y%m%d_%H%M')
        p = os.path.join(fo, fb + '_' + dt + ext)

    ensure_presence_of_directory(p)

    return p


def split_text_to_multiple_rows(df, column, separator):
    """
    Separates entries, that have been separted by separator
    within a given
    column of the dataframe df into separate rows

    Input:
        df          DataFrame
        column      String; Name of the column,
                        which contains records separated by separator
        separator   String: Regular Expression that
                        separates the records in column

    Output:
        df_stacked  DataFrame, where records in
                        column have been separated into individual rows

    """

    # Separate rows, that should be processed (to save overhead)
    f = df[column].str.contains(separator)
    dff_s = df.loc[f, :]    # _s  separate
    dff_ns = df.loc[~f, :]  # _ns not separate
    orig_column_order = df.columns  # track original order of columns

    # Separate records and place them in separate rows
    s = dff_s[column].str.split(separator).apply(pd.Series, 1).stack()
    s.index = s.index.droplevel(-1)
    s.name = column
    del dff_s[column]
    dff_s = dff_s.join(s)

    # Recreate DataFrame in original input format
    dff_s = dff_s.reindex(orig_column_order, axis=1)
    df = pd.concat([dff_s, dff_ns])
    df_stacked = df

    # Convert to columns
    # df.join(s.apply(lambda x: Series(x.split(':'))))

    return df_stacked


def rotate(origin, point, degrees):
    """
    Rotate a point counterclockwise by a given angle around a given origin.

    The angle should be given in degrees.
    """

    angle = math.radians(degrees)

    ox, oy = origin
    px, py = point

    qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
    qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
    return qx, qy
