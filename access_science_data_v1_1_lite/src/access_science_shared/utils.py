import pandas as pd


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

    return df_stacked
