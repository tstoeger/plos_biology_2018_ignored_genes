# Large-scale investigation of the reasons why potentially important genes are ignored.

## Background

This repository contains copies of several repositories used for Stoeger, Gerlach, Morimoto, Amaral, Plos Biology, 2018.

## Content

An overview of the data resources and repositories is given in overview_of_datasets_lite.pdf .

Briefly, robusa (or rbusa for legacy reasons) collects information outside of the biological sciences, whereas geisen collects information from the biological sciences. Access_science_data is used to read information from data repositories, and resci contains high-level analyses for individual figure panels.

- rbusa_main v2.0.1  -> all
- geisen_main v1_2_1  -> all
- geisen_manual v1.0.1 -> all
- access_science_data v1.1 -> subset relevant to Stoeger et al., Plos Biology, 2018
- resci v1.1  -> subset relevant to Stoeger et al., Plos Biology, 2018

## Data Sources


Note that raw data of used data sources can not be shared for legal reasons. For data, which will not be automatically downloaded through the in-built downloading function of geisen, please see our mnanuscript for data sources, and ways to license them. 

## Setup

Location of data repositories (or output folders while generating them) can be set at the top of the following modules. 

- access_science_data_v1_1_lite/src/access_science_shared/inout.py
- geisen_main_v1_2_1/src/inout.py
- rbusa_main_v2_0_1/src/rbusa/inout.py
- resci_v_1_1_lite/PlosBiology2018/src/resci_tools.py
- resci_v_1_1_lite/PlosBiology2018/src/resci_inout.py

## Is this everything? I'd like to use or enhance this framework.

The primary role of this repository is to supplement the publication of Stoeger, Gerlach, Morimoto, Amaral, Plos Biology, 2018. As hinted at by the lite versions of individual repositories, we have created a framework to access essentially any gene-related information from Python. This framework exceeds the scope of the current publication, and it exceeds commonly used bioinformatic tools. If interested to contribute, or use parts of our full framework, please do not hesitate to contact thomas.stoeger@northwestern.edu
