# NOTE

This is a subset of the access_science_data framework used by Thomas Stoeger on 2018-July-26. It contains functions relevant to a Plos Biology manuscript (and only refers to the subset of curated bioogy contained within it).


# access_science_data
Repository to accessing science of science data

This repository contains functions for accessing scientific data, but not for processing it.

Its data sources can be heterogeneous, and not all target data sources are required for all functions of access_science_data

## Example Datasets:

### geisen_main   

A curated (and preprocessed) aggregation of all high-level biological information. The scope of that source are all information which relates directly to individual genes or gene products. Curator: Thomas Stoeger

### rbusa_main

A curated (and preprocessed) aggregation of societal and economic data that does not directly relate to genes, even if it can be of interest to meta-studies of science. Curator: Thomas Stoeger

### medline_wos_main

A linkage of gene-containing MedLine to WebOfScience. Curator: Martin Gerlach

### gtx_atlas

EMBL-EBI gene expression atlas of all public high-quality gene expression studies.


# Setup

a) Add the src folder to your path profile
(on MacOS:  the hidden text file .bash_profile in your home folder)

add the following line
export PYTHONPATH="/Users/tstoeger/Projects/access_science_data/src:$PYTHONPATH"

Now access_science_data can be used from other code.

b) Within /src/access_science_shared/ add the local paths of the data repostiories to the top of inout.py
