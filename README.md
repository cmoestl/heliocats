# heliocats


This package contains codes used for the creation of catalogs of interplanetary coronal mass ejections and high speed solar wind streams, and their analysis. 

by [C. Möstl](https://www.iwf.oeaw.ac.at/en/user-site/christian-moestl/)

Current status (January 2020): **Work in progress!** 

This is a continuation of work done in the EU HELCATS project (2014-2017): https://www.helcats-fp7.eu, in particular concerning ICMECAT (working package 4).

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


## Programs

Activate the environment (*conda activate helio*) and run with *python program_name.py*:

- icmecat_maker.py        
processes data into a normalized format and creates the ICMECAT catalog.

- sc_positions_insitu.py  
makes movies of spacecraft positions and in situ data.

- data_update.py          
makes real time downloads of various data sources.

- cme_stats.py            
makes CME statistics.