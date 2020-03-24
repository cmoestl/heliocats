# heliocats


This package contains codes used for the creation of catalogs of interplanetary coronal mass ejections and high speed solar wind streams, and their analysis. 

by [C. Möstl](https://www.iwf.oeaw.ac.at/en/user-site/christian-moestl/)

Current status (March 2020): **Work in progress!** 

This is a continuation of work done in the EU HELCATS project (2014-2017): 
[https://www.helcats-fp7.eu](https://www.helcats-fp7.eu), in particular concerning 
the Interplanetary Coronal Mass ejection CATalog ICMECAT (working package 4).

If you want to use parts of this code for generating results for peer-reviewed scientific publications, please contact me per email (christian.moestl@oeaw.ac.at) or via https://twitter.com/chrisoutofspace .


---
## Usage

For running the jupyter notebooks (files with .ipynb) or the python scripts (.py), first activate the helio environment:

    conda activate helio

Folder *data/* contains all data files needed and produced during the analysis. 
All outputs can be found in the folder *results/* or subfolders. Jupyter notebooks can be converted to scripts by e.g. for icmecat.ipynb to icmecat.py:

    jupyter nbconvert --to script icmecat.ipynb
    

### ICME catalog 


**work in progress, don't use yet!**

Before running the icmecat scripts, you need to download 6 data files (in total 4.1 GB) from this 
figshare repository: [https://doi.org/10.6084/m9.figshare.11973693.v1](https://doi.org/10.6084/m9.figshare.11973693.v1)

and place them in the data/ folder.

    data/helcats_all_data_removed.p
    data/maven_2014_2018_removed_smoothed.p
    data/psp_2018_2019.p
    data/wind_2018_2020.p
    data/sta_2018_2019_beacon.p
    data/ulysses_1990_2009_helcats.p


This creates the ICMECAT catalog (also a jupyter notebook is available):

    python icmecat.py

The catalog is available in these formats: .p, .xlsx, .json, .csv, .html, .txt   

Load this catalog with 

    import pickle
    file='icmecat/HELCATS_ICMECAT_v20.p'
    ic=pickle.load( open(file, 'rb'))
    
"ic" is a pandas dataframe, the names of all parameters can be seen with 

    ic.keys()

and the data can be accessed by
    
    ic.icmecat_id
    ic.icme_start_time
    ic.mo_bmax
    ic.mo_sc_heliodistance
    ...

which works particularly well in ipython or a jupyter notebook.

All times (ic.icme_start_time, ic.mo_start_time, ic.mo_end_time) are python datetime objects. 
They can be used directly plotting with matplotlib, but can also easily 
converted into many formats by:

    from sunpy.time import parse_time
    parse_time(ic.icme_start_time).plot_date
    parse_time(ic.icme_start_time).iso
    ...



### ICME statistics

use either (depending on your preference for jupyter lab, notebook or a script):

    jupyter lab cme_statistics.ipynb
  
    jupyter notebook cme_statistics.ipynb

    python cme_statistics.py
    
These codes make CME statistics to get the results and plots for the paper Möstl et al. (2020, in preparation). 


### Spacecraft positions and data

    python sc_positions_insitu.py

makes movies of spacecraft positions and in situ data.


### Real-time data update

    python data_update.py
    
makes real time downloads and plots of various data sources.





---

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
	  

