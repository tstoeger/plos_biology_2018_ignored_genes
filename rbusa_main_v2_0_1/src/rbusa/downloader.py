import gzip
import os
import re
import shutil
import time
import urllib

import numpy as np

import rbusa.inout as io

from rbusa.settings import retreive_general_download_settings


def download_data_set(name_of_dataset, folder_contains_dots=False):
    """
    Gets the path of a download.

    Input:
        name_of_dataset     string, as defined in download settings file
                            (in /cfg/download_links.csv)
        folder_contains_dots default is False; set to True that folder will
                            be present even if it contains a dot
    Output:
        locations_of_storage    list, contains storage path of every file
    """

    df = retreive_general_download_settings()

    if (name_of_dataset in df.index) is False:
        raise ValueError(
            'Could not find {} in settings file.'.format(name_of_dataset))
    else:
        df = df.loc[[name_of_dataset], :]

    print('Initialize download of {}'.format(name_of_dataset))

    locations_of_storage = []

    for r in df.iterrows():

        info = r[1]
        location_to_download = info['location_on_server']

        u = urllib.parse.urlparse(location_to_download)
        after_domain = u[2]
        subpath = 'downloads/{}{}'.format(
            info['location_on_present_machine'], after_domain)
        location_to_store = io.get_internal_path(subpath)

        if not folder_contains_dots:
            io.ensure_presence_of_directory(
                os.path.dirname(            # embed for files lacking extension
                    location_to_store))
        else:
            io.ensure_presence_of_directory(location_to_store)

        if os.path.isfile(location_to_store):
            raise EnvironmentError(
                '{} has already been downloaded.'.format(name_of_dataset))
        else:
            download(location_to_download, location_to_store)
            locations_of_storage.append(location_to_store)

    return locations_of_storage


"""
    #################   SUPPORT FUNCTIONS #########################

    download     downloads
    isolate_valid_replicons        extracts valid replicons from gbff
    ungzip       decompresses gz file (gnu zip)
    unzip        decompresses zip file

"""


def download(location_to_download, location_to_store):
    """
    Tries to download a file, and throws error if it repeatedly fails

    Input:
        location_to_download    str, location on server
        location_to_store       str, location on client
    """

    # Attempt to download
    has_been_downloaded = False
    attempts = 0
    while attempts < 4:
        try:
            urllib.request.urlretrieve(
                location_to_download, location_to_store)
            print('Succesfully downloaded {} to {}'.format(
                location_to_download, location_to_store))
            has_been_downloaded = True
            break
        except:
            print('Failed to download {} . Will retry in two minutes'.format(
                location_to_download))
            time.sleep(120)
            attempts += 1

    if not has_been_downloaded:
        raise EnvironmentError('Failed to download {} .'.format(
            location_to_download))


# def isolate_valid_replicons(p_genomic_out, taxon_of_interest):
#     """
#     Will take a _genomic.gff.gz and isolate standard
#     replicons and safe their sequence in plain text
#     - Compatibility with biopython (which does not support .gz)
#     - ignore non standard-replicons (e.g.: fragments),
#         NOTE: ignoring non standard-replcions might lead
#         to a dismmissal of some valid genes (e.g. birds having
#         microchromosomes, which are still not properly assembled)

#     Input:
#         p_genomic_out       str, path to genomic .gbff.gz file

#     Output:
#         outfolder           str, path where output is stored

#     """

#     # Support function to write replicion, which should be synonymous
#     # with biological chromosome
#     def write_replicion(lines_to_write, current_replicon):
#         """
#         Subfunction to write cached lines to genbank format

#         Input:
#             lines_to_write      list; lines that should be stored
#             current_replicons  str; name of current replicon

#         """

#         out_path = os.path.join(outfolder, current_replicon + '.gbff')

#         if os.path.isfile(out_path):
#             print('Skipping file. Output already exists.', out_path)
#         else:
#             with open(out_path, 'x') as outfile:
#                 outfile.writelines(lines_to_write)

#     # Main code (processing of genome files) ########################

#     # Initialize
#     # Get settings
#     settings = retreive_genome_download_settings(taxon_of_interest)

#     # Ensure output folder without genbank
#     outfolder = p_genomic_out.replace('.gbff.gz', '')
#     if not os.path.exists(outfolder):
#             os.makedirs(outfolder)
#     n = io.check_number_of_files_in_directory(outfolder, '.gbff')
#     if n > 0:
#         raise EnvironmentError('genomic.gbff folder is not empty')

#     # For error check: obtain tracker of anticipated replicons
#     anticipated_replicons = [
#         re.findall(
#             '^([A-Z0-9_]*).*', r)[  # refseq without version
#             0] for r in settings['replicons'].values()]
#     anticipated_replicons = np.array(anticipated_replicons)
#     is_in_read_file = np.zeros(len(anticipated_replicons))

#     lines_to_write = []
#     current_chromosome = ''
#     with gzip.open(p_genomic_out, 'r') as infile:

#         is_in_allowed_chromosome = False

#         for line in infile:
#             line = line.decode('utf-8')

#             # Track, whether one is within a standard replicon
#             if line.startswith('LOCUS'):

#                 if len(lines_to_write) > 0:
#                     if is_in_allowed_chromosome:  # standard replicon
#                         write_replicion(lines_to_write, current_chromosome)

#                 current_chromosome = re.findall(
#                     '^LOCUS *([A-Z0-9_]*).*', line)[0]
#                 if any(current_chromosome):

#                     f = anticipated_replicons == current_chromosome
#                     if any(f):
#                         is_in_read_file[f] = is_in_read_file[f] + 1
#                         is_in_allowed_chromosome = True
#                     else:
#                         is_in_allowed_chromosome = False
#                 else:
#                     is_in_allowed_chromosome = False

#                 if is_in_allowed_chromosome:
#                     lines_to_write = []

#             if is_in_allowed_chromosome:
#                 lines_to_write.append(line)

#         # Write last contig
#         if is_in_allowed_chromosome:
#             write_replicion(lines_to_write, current_chromosome)

#         # Check for completeness
#         if any(is_in_read_file > 1):
#             raise ValueError(
#                 'One replicon appears present several times.'
#                 'Every replicon should only appear once.')
#         elif any(is_in_read_file < 1):
#             raise ValueError(
#                 'One replicon appears absents absent.'
#                 'Every replicon should appear once.')

#     return outfolder


def ungzip(source_file, return_path=False):
    """
    Will decompress a gnu zip file at its current location.
    Aniticpates, and will remove .gz file extension.

    Input:
        source_file     str, path to file
        return_path     bool, shall output path be returned?

    Output:
        p_out           str, path to decompressed file
    """

    anticipated_extension_pattern = '\.gz$'

    if len(re.findall(anticipated_extension_pattern, source_file)) == 0:
        raise ValueError('Could not find anticipated file extension')

    p_in = source_file
    p_out = re.sub('\.gz$', '', source_file)

    with gzip.open(p_in, 'rb') as f_in, open(p_out, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

    if return_path:
        return p_out


def unzip(source_file, delete_zip=False):
    """
    Will decompress a zip file into a subfolder at current location

    Input:
        source_file     str, path to file
    """

    import zipfile

    dest_dir, ext = os.path.split(source_file)

    with zipfile.ZipFile(source_file) as zf:
        zf.extractall(dest_dir)

    if delete_zip:
        os.remove(source_file)
