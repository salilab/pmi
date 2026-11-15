#!/bin/bash -e

# Set up an environment to run CI tests, e.g. with GitHub Actions or Travis

eval "$(conda shell.bash hook)"
conda activate test
conda install --yes -q pip scipy matplotlib imp-nightly libboost-devel gxx_linux-64 eigen cereal swig cmake
pip install pytest-cov coverage pytest-flake8 'pytest<9'
