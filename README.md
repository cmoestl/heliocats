# heliocats


This package contains codes used for the creation of catalogs of interplanetary coronal mass ejections and high speed solar wind streams, and their analysis. 

by [C. Möstl](https://www.iwf.oeaw.ac.at/en/user-site/christian-moestl/)

Current status (March 2020): **Work in progress!** 

This is a continuation of work done in the EU HELCATS project (2014-2017): 
[https://www.helcats-fp7.eu](https://www.helcats-fp7.eu), in particular concerning ICMECAT (working package 4).

If you want to use parts of this code for generating results for peer-reviewed scientific publications, please contact me per email (christian.moestl@oeaw.ac.at) or via https://twitter.com/chrisoutofspace .



## Installation 

Install python 3.7.6 with miniconda:

on Linux:

	  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
	  bash Miniconda3-latest-Linux-x86.sh

on MacOS:

	  curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
	  bash Miniconda3-latest-MacOSX-x86_64.sh

Create a conda environment:

	  conda env create -f environment.yml

	  conda activate helio

	  pip install -r requirements.txt
	  
go to a directory of your choice

	  git clone https://github.com/cmoestl/heliocats


## Usage

For running the jupyter notebook (files with .ipynb), first activate the helio environment:

    conda activate helio
    
and then use either (depending on your preference for lab or notebook):

    jupyter lab cme_statistics.ipynb
    jupyter notebook cme_statistics.ipynb


For the python scripts, activate the environment (*conda activate helio*) and run:

    python icmecat_maker.py

processes data into a normalized format and creates the ICMECAT catalog.

    python cme_statistics.py
    
makes CME statistics to get the results and plots for the paper Möstl et al. (2020, in preparation).

    python sc_positions_insitu.py

makes movies of spacecraft positions and in situ data.

    python data_update.py
    
makes real time downloads of various data sources.



Folder *data/* contains all data files needed and produced during the analysis. 
All outputs can be found in the folder *results/*

