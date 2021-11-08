# -*- coding: utf-8 -*-

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="heliocats",
    version="1.2.1",
    author="Christian Moestl",
    author_email="christian.moestl@oeaw.ac.at",
    keywords=["geophysics", "heliophysics", "space weather"],
    description="Making heliospheric event catalogs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/cmoestl/heliocats",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)
