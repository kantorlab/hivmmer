from setuptools import setup, find_packages

with open("hivmmer/VERSION") as f:
    version = f.read().strip()

setup(
    name="hivmmer",
    version=version,
    author="Mark Howison",
    author_email="mhowison@ripl.org",
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
        "Programming Language :: Python :: 3.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics"],
    provides=["hivmmer"],
    install_requires=["BioPython>=1.69", "matplotlib>=3.1.1", "numpy>=1.13.0", "openpyxl", "pandas>=0.22.0", "xlrd"],
    packages=find_packages(),
    package_data={"hivmmer": ["VERSION", "*.csv", "*.hmm.*", "*.tsv"]},
    scripts=["scripts/hivmmer"],
    entry_points={
        "console_scripts": ["hivmmer-filter=hivmmer.filter:_run",
                            "hivmmer-translate=hivmmer.translate:_run"]
    }
)
