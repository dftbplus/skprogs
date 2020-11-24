#!/usr/bin/env python3
from distutils.core import setup

setup(
    name="sktools",
    version='20.2',
    description="Tools to create SK-parameters",
    author="DFTB+ developers",
    url="http://www.dftbplus.org",
    platforms="platform independent",
    package_dir={"": "src"},
    packages=["sktools", "sktools.hsd", "sktools.calculators", "sktools.skgen"],
    scripts=[
        "bin/skgen",
    ],
    classifiers=[
        "Programming Language :: Python",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering",
    ],
    long_description="""
Processing and converting data related to the DFTB+ package
-----------------------------------------------------------
A few scripts which should make the life of DFTB+ users easier, by providing
functions to process and convert various DFTB+ data formats.
""",
    requires=[ "numpy" ]
)
