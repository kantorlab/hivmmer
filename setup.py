#!/usr/bin/env python

from glob import glob
from setuptools import setup

setup(
    name="hivmmer",
    version="0.1.2",
    author="Mark Howison",
    author_email="mhowison@brown.edu",
    url="https://github.com/kantorlab/hivmmer",
    description="""
        An alignment and variant-calling pipeline for Illumina deep sequencing of
        HIV-1, based on the probabilistic aligner HMMER.""",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Topic :: Scientific/Engineering :: Bio-Informatics"],
    install_requires=["BioPython>=1.69", "numpy>=1.13.0", "pandas>=0.22.0"],
    scripts=glob("scripts/*")
)
