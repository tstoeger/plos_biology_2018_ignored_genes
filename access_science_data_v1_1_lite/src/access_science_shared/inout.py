import os

locations = [   # upper paths will be preferred
    '/Users/tstoeger/Dropbox/Work/resources',
    '/home/workspace/scibio/resources'
]

datasets = {
    'extra_geisen': 'geisen/extra_geisen',
    'geisen': 'geisen/geisen_main_v1_2',
    'gtx_atlas': 'gtx_atlas/gtx_atlas_170308',
    'medline_wos': 'medline_wos_main/medline_wos_main_v1_1_1',
    'rbusa': 'rbusa/rbusa_main_v1_1',
    'publications': 'publications',
    'drugbank': 'drugbank/drugbank_v5_0_7',
    'geisen_manual': 'geisen/geisen_manual_v1_1',
    'ebi_gwas': 'gwas/embl_gwas_1_0_x',
    'genome_rnai': 'genome_rnai/GenomeRNAi_v17_AllScreens',
}


def get_path(dataset=None, extension=None):
    '''
    Returns subfolder within internal part of geisen.

    Input:
        dataset     str, name of dataset, e.g.: geisen
        extension   str, optional, subfolder
    Output:
        outpath     str, folder within internal part of geisen
    '''

    is_found = False
    for l in locations:    # upper paths will be preferred
        if not is_found:
            p = os.path.join(l, datasets[dataset])
            if os.path.exists(p):
                ref_folder = p
                is_found = True
    if not is_found:
        raise EnvironmentError(
            'Could not find location of {} dataset'.format(
                dataset))

    if extension is not None:
        extension = str.replace(extension, '\\', os.path.sep)
        extension = str.replace(extension, '/', os.path.sep)

        outpath = os.path.join(ref_folder, extension)
    else:
        outpath = ref_folder

    return outpath
