import os
import socket

_thomas_defaut_support_path = '~/Dropbox/Work/dynamic/resci_support_data'

if any(socket.gethostname() == x for x in [
        'Thomass-MacBook-Pro.local',
        'Thomass-MBP',
        'Thomass-iMac.local',
]):

    _internal_main_folder = _thomas_defaut_support_path

elif os.path.exists(_thomas_defaut_support_path):
    _internal_main_folder = _thomas_defaut_support_path

else:
    raise ValueError(
        'Did not yet set reference paths for current machine, {}'.format(
            socket.gethostname()))


def get_number_of_files(dir_path, pattern):
    """
    Counts the number of files in a given directory.extension

    Input:
        dir_path    str, path to folder
        pattern     str, pattern that should be used for finding files

    Output:
        number_of_files     int, number of files in dir_path that match pattern

    """

    from fnmatch import filter

    number_of_files = len(filter(os.listdir(dir_path), pattern))
    return number_of_files


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


def get_internal_path(extension=None):
    '''
    Returns subfolder within internal part of rbusa.

    Input:
        extension   str, optional, subfolder
    Output:
        outpath     str, folder within internal part of rbusa
    '''

    if extension is not None:
        extension = str.replace(extension, '\\', os.path.sep)
        extension = str.replace(extension, '/', os.path.sep)

        outpath = os.path.join(_internal_main_folder, extension)
    else:
        outpath = _internal_main_folder

    return outpath


def ensure_absence_of_file(file_path):
    """
    Throws error, if path already exists.

    Input:
        file_path   str, path to file or folder
    """

    abs_path = os.path.abspath(file_path)
    if os.path.exists(abs_path):
        raise EnvironmentError('{} already exists'.format(abs_path))


def check_number_of_files_in_directory(dir_path, pattern):
    """
    Counts the number of files in a given directory.extension

    Input:
        dir_path    str, path to folder
        pattern     str, pattern that should be used for finding files

    Output:
        number_of_files     int, number of files in dir_path that match pattern

    """

    from fnmatch import filter

    number_of_files = len(filter(os.listdir(dir_path), pattern))
    return number_of_files
