# heliocats

This package contains codes used for the creation of catalogs of interplanetary coronal mass ejections (ICMEs) and their analysis. [This is a link to a google colab notebook](https://colab.research.google.com/drive/1_zJMGJnX3XJx7FCHD04SC3Y0KCMdDGMz) for instructions how to read the catalogs produced with this package. The ICMECAT as the main result to be produced with this package is published at https://helioforecast.space/icmecat.


**Authors**: C. Möstl, E. E. Davies, E. Weiler, Austrian Space Weather Office, Geosphere Austria https:/helioforecast.space; contributions by A. J. Weiss, R. L. Bailey, M. A. Reiss, C. L. Simon Wedlund, A. Isavnin, R. M. Winslow.

Last update: April 2025.

This is a continuation of work done in the EU HELCATS project (2014-2017): [https://www.helcats-fp7.eu](https://www.helcats-fp7.eu). 

If you want to use parts of this code for generating results for **peer-reviewed scientific publications, 
please contact me per email** (chris.moestl@outlook.com) for co-authorships.


---
## Usage

For running the jupyter notebooks (files with .ipynb) or the python scripts (.py), first activate the helio4 environment (for installation instructions, see bottom of this readme):

    conda activate helio4

A data folder location can be set in *config.py* (e.g. *data/*) that contains all data files needed during the analysis. Outputs can be found in the folder *results/* or subfolders therein, e.g. the files created for the ICMECAT are in folder *icmecat/*. Jupyter notebooks are converted to scripts in their first cell.
    
For some scripts there should be a systemwide installation of ffmpeg available.
    
---



### ICME catalog (work in progress)

Before running the icmecat scripts, you need to download data files the we made (in total around 10 GB) from this figshare repository: 
[https://doi.org/10.6084/m9.figshare.11973693](https://doi.org/10.6084/m9.figshare.11973693)
and place these files in the a folder e.g. named "data", 
(the name of this folder is set by the variable *data_path* in file config.py).

This creates the ICMECAT catalog (also a jupyter notebook *icmecat.ipynb* is available, which can be run and edited by using jupyter notebook or jupyter lab):

    python icmecat.py

The catalog is available in these formats: .p (pandas dataframe or numpy record array), .xlsx, .json, .csv, .html, .txt, .h5 (hdf5)   

Load this catalog into python with 

    import pickle
    file='icmecat/HELCATS_ICMECAT_v23_pandas.p'
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
    file='icmecat/HELIO4CAST_ICMECAT_v23_numpy.p'
    [ic_nprec,ic_np,header,parameters]=pickle.load( open(file, 'rb'))  

which returns a numpy record array (ic_nprec) or a numpy structured array (ic_np) consisting of strings and floats.

In the pandas dataframe, all times (ic.icme_start_time, ic.mo_start_time, ic.mo_end_time) are python datetime objects. 
In the numpy arrays, these times are given in matplotlib format (fractional days since 1970-01-01). 
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




---

## Installation 

Install python with miniconda:

on Linux:

	  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
	  bash Miniconda3-latest-Linux-x86_64.sh

on MacOS:

	  curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
	  bash Miniconda3-latest-MacOSX-x86_64.sh

go to a directory of your choice

	  git clone https://github.com/cmoestl/heliocats
	  

Create a conda environment using the "envs/env_helio4.yml", and activate the environment:

	  conda env create -f env_helio4.yml

	  conda activate helio4


Some codes at one point may need to use ffmpeg, which can be downloaded for Mac and Linux from this site: https://ffmpeg.org/download.html

or install from the command line

on Mac 
	  brew install ffmpeg

on Linux
	  sudo apt install ffmpeg




MIT LICENSE
Copyright 2020-2025, Christian Moestl, Emma Davies, Eva Weiler
Permission is hereby granted, free of charge, to any person obtaining a copy of this 
software and associated documentation files (the "Software"), to deal in the Software
without restriction, including without limitation the rights to use, copy, modify, 
merge, publish, distribute, sublicense, and/or sell copies of the Software, and to 
permit persons to whom the Software is furnished to do so, subject to the following 
conditions:
The above copyright notice and this permission notice shall be included in all copies 
or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
