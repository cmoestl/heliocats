# heliocats


This package contains codes used for the creation of catalogs of interplanetary coronal mass ejections and their analysis. [This is a link to a google colab notebook](https://colab.research.google.com/drive/1_zJMGJnX3XJx7FCHD04SC3Y0KCMdDGMz) for instructions how to read the catalogs produced with this package.


**Authors**: [C. Möstl](https://www.iwf.oeaw.ac.at/en/user-site/christian-moestl/), A. J. Weiss, R. L. Bailey, M. A. Reiss, C. L. Simon Wedlund, A. Isavnin, R. M. Winslow, D. Stansby.

Current status (October 2021): **work in progress** 

This is a continuation of work done in the EU HELCATS project (2014-2017): 
[https://www.helcats-fp7.eu](https://www.helcats-fp7.eu). This package is used for updates of the the Interplanetary Coronal Mass ejection CATalog (ICMECAT).

If you want to use parts of this code for generating results for **peer-reviewed scientific publications, 
please contact me per email** (christian.moestl@oeaw.ac.at) or via https://twitter.com/chrisoutofspace for co-authorships.




---
## Usage

For running the jupyter notebooks (files with .ipynb) or the python scripts (.py), first activate the helio environment (for installation instructions, see bottom of this readme):

    conda activate helio

A data folder location can be set in *config.py* (e.g. *data/*) that contains all data files needed during the analysis. Outputs can be found in the folder *results/* or subfolders therein, e.g. the files created for the ICMECAT are in folder *icmecat/*. Jupyter notebooks can be converted to scripts by e.g. for icmecat.ipynb to icmecat.py:

    jupyter nbconvert --to script icmecat.ipynb
    
---




### ICME catalog 

Before running the icmecat scripts, you need to download 10 data files for 8 spacecraft the we made 
(in total 7.5 GB) from this figshare repository: 
[https://doi.org/10.6084/m9.figshare.11973693](https://doi.org/10.6084/m9.figshare.11973693)
and place these files in the a folder e.g. named "data", 
(the name of this folder is set by the variable *data_path* in file config.py):

    data/wind_2018_2019_heeq.p
    data/wind_2007_2018_heeq_helcats.p
    data/psp_2018_2019_sceq.p
    data/stereoa_2007_2019_sceq.p
    data/stereoa_2019_2020_sceq_beacon.p
    data/stereob_2007_2014_sceq.p
    data/maven_2014_2018_removed_smoothed.p
    data/ulysses_1990_2009_rtn.p
    data/vex_2007_2014_sceq_removed.p
    data/messenger_2007_2015_sceq_removed.p
 

This creates the ICMECAT catalog (also a jupyter notebook *icmecat.ipynb* is available, which can be run by using jupyter notebook or jupyter lab):

    python icmecat.py

The catalog is available in these formats: .p (pandas dataframe or numpy record array), .xlsx, .json, .csv, .html, .txt, .h5 (hdf5)   

Load this catalog into python with 

    import pickle
    file='icmecat/HELCATS_ICMECAT_v20_pandas.p'
    [ic,header,parameters]=pickle.load( open(file, 'rb'))
    
    
"ic" is a pandas dataframe, the names of all parameters can be seen with 

    ic.keys()

and the data can be accessed by
    
    ic.icmecat_id
    ic.icme_start_time
    ic.mo_bmax
    ic.mo_sc_heliodistance
    ...

which works particularly well in ipython or a jupyter notebook. Further, 

    print(header)
    print(parameters)
    
gives the header description and the list of all parameters with units.     

Alternatively, you can load the ICMECAT with 

    import pickle
    file='icmecat/HELCATS_ICMECAT_v20_numpy.p'
    [ic_nprec,ic_np,header,parameters]=pickle.load( open(file, 'rb'))  

which returns a numpy record array (ic_nprec) or a numpy structured array (ic_np) consisting of strings and floats.

In the pandas dataframe, all times (ic.icme_start_time, ic.mo_start_time, ic.mo_end_time) are python datetime objects. 
In the numpy arrays, these times are given in matplotlib format (fractional days since year 0001 Jan 1, plus 1 day). 
Both formats can be used directly plotting with matplotlib, but can also easily 
converted into many other formats by (given you have sunpy installed):

    from sunpy.time import parse_time
    parse_time(ic.icme_start_time).plot_date
    parse_time(ic.icme_start_time).iso
    ...


### HI CME arrival catalog 

use either (depending on your preference for jupyter lab, notebook or a script):

    jupyter lab arrcat.ipynb
  
    jupyter notebook arrcat.ipynb

    python arrcat.py
    
These codes makes the HELCATS CME arrival catalog, see e.g. Möstl et al. (2017, Space Weather). The catalog is available in essentially the same formats as the ICMECAT, and can be used similarly to above.





### ICME statistics

use either (depending on your preference for jupyter lab, notebook or a script):

    jupyter lab cme_statistics.ipynb
  
    jupyter notebook cme_statistics.ipynb

    python cme_statistics.py
    
These codes make CME statistics to obtain results for planetary space weather studies.


### Spacecraft positions and data

    python sc_positions_insitu.py

makes movies of spacecraft positions and in situ data, all configuration variables are set in the script.


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

go to a directory of your choice

	  git clone https://github.com/cmoestl/heliocats
	  

Create a conda environment using the environment.yml and requirements.txt file in the heliocats root directory, and activate the environment in between:

	  conda env create -f environment.yml

	  conda activate helio

	  pip install -r requirements.txt
	  


