# Validation for alignment score thresholding

## Setup

The required python packages can be installed with Anaconda Python using:

    conda create -n hivmmer-validation -c kantorlab hivmmer=0.1.3 sra-tools=2.9.1.1 scons=3.1.1 r-ggplot2 r-tidyverse
    conda activate hivmmer-validation

Next, replace hivmmer 0.1.3 with the current version on this branch, with:

    pip install -e ..

## Estimated Thresholds

The estimated threshold is the minimum of the lowest observed E-value across
the "invalid" data sets (mismatched reference sequences, random seequences, and
UniVec sequences).

| region | E-value threshold |
| ------ | ----------------- |
| env    | 0.00148 |
| gag    | 0.00180 |
| nef    | 0.00065 |
| pol    | 0.000464 |
| tat    | 0.000245 |
| vif    | 0.000543 |
| vpr    | 0.000935 |
| vpu    | 0.00110 |

