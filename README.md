![conda](https://img.shields.io/conda/v/kantorlab/hivmmer.svg) ![python](https://img.shields.io/badge/python-3.7-blue.svg)

# hivmmer

An alignment and variant-calling pipeline for Illumina deep sequencing of
HIV-1, based on the probabilistic aligner [HMMER](http://hmmer.org).

## Pipeline steps

1. Constructs an amino acid profile Hidden Markov Model (pHMM) from a multiple
   sequence alignment of all HIV-1 Group M amino acid sequences publicly
   available in the [Los Alamos HIV Sequence Database](http//www.hiv.lanl.gov) for
   the a given gene or region of the HIV genome.
2. Filters the NGS data, retaining reads with length >75 and mean PHRED quality
   score > 25, and consolidates duplicate sequences (with `hivmmer-filter`).
   The number of duplicates are tracked to enable correct inference of frequencies
   later in the pipeline.
3. Translates each de-duplicated sequence into all six possible frames (forward
   and reverse), retaining only the translated sequences that contain no stop
   codons (with `hivmmer-translate`).
4. Aligns the translated reads to the reference pHMM with hmmsearch, producing
   a multiple sequence alignment of translated reads.
5. Constructs a sample-specific amino acid pHMM from the multiple sequence
   alignment of translated reads.
6. Repeats the HMMER alignment against the sample-specific pHMM for increased
   sensitivity.
7. Maps the translated amino acid coordinates in the alignment to the original
   frame and coordinates in the nucleotide reads to construct a codon frequency
   table (adjusting the counts for duplicate reads; with `hivmmer-codons`).

## Usage

```
hivmmer --id ID --fq1 FASTQ1 --fq2 FASTQ2 --ref REFERENCE [--cpu N]
        [-h|--help] [-v|--version]
```

`ID` specifies a name for the analysis that will be used as the basename for
all output.

`FASTQ1` and `FASTQ2` are the forward and reverse Illumina reads.

Optionally, you can use `N` threads to speed-up the HMMER stages of the pipeline.

## Installation

**hivmmer requires Python 3.7**

### Quick install with Anaconda Python

On 64-bit Linux, it is also possible to install hivmmer using prebuilt
packages from the [kantorlab Anaconda channel](https://anaconda.org/kantorlab).

First, install the [Anaconda](https://www.continuum.io/anaconda-overview)
or Miniconda distribution of Python 3.

Once the `conda` command is in your PATH, hivmmer and all its dependencies can
be installed into its own isolated conda environment with the single command:

    conda create -c kantorlab -n hivmmer hivmmer

Once installed, activate the `hivmmer` conda environment with:

    source activate hivmmer

This will place hivmmer and all its dependencies in your PATH.

We have primarily tested hivmmer on CentOS 6.8, but in theory it should run on
any 64-bit Linux system with glibc >= 2.12.

All relevant conda recipes are available from the
[Kantor Lab's conda-recipes repository](https://github.com/kantorlab/conda-recipes).

### Quick install with Docker

On systems other than 64-bit Linux, you can run hivmmer via a
[Docker](https://www.docker.com) container.

First, visit the [Docker website](https://www.docker.com) to download and
install Docker for your host operating system.

Second, pull the pre-compiled hivmmer Docker image, which includes all dependencies,
from DockerHub:

    docker pull kantorlab/hivmmer

Each time you want to use Agalma, run the docker image with:

    docker run -it kantorlab/hivmmer

This will launch a new Docker container with hivmmer, and provide an
interactive prompt to access to the container.

### Manual installation

hivmmer can be installed with pip using the included setup.py, and has the
following dependencies on external tools (which must be in your PATH):

* HMMER 3.2.1
* PEAR 0.9.11

Note: PEAR source code is available under an academic license from [https://www.h-its.org/en/research/sco/software/#NextGenerationSequencingSequenceAnalysis](https://www.h-its.org/en/research/sco/software/#NextGenerationSequencingSequenceAnalysis).

### pHMM references

hivmmer comes with prepackaged amino acid profile Hidden Markov Models for the entire
HIV genome, based on curated multiple-sequence alignments downloaded from the
[Los Alamos HIV Sequence Database](http//www.hiv.lanl.gov).

## Funding

Development of hivmmer is made possible through funding from the National
Institutes of Health under awards R01AI108441, R01AI120792, R01AI147333, and
the Providence/Boston Center for AIDS Research (P30AI042853).

## Authors

Mark Howison (<mhowison@ripl.org>)

For bug reports and questions, please create an
[issue on Github](https://github.com/kantorlab/hivmmer/issues).

## License

Copyright 2018, Brown University, Providence, RI. All Rights Reserved.

Copyright 2019-2020, Innovative Policy Lab (d/b/a Research Improving People's Lives), Providence, RI. All Rights Reserved.

See LICENSE for full terms of use.

