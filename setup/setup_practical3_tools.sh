#!/bin/bash
set -e

echo "Creating and activating micromamba environment practical3-env..."
micromamba create -y -n practical3-env python=3.11
eval "$(micromamba shell hook -s bash)"
micromamba activate practical3-env

echo "Installing structure and clumpp..."
micromamba install -y -c conda-forge -c bioconda structure clumpp

echo "Attempting to install structureHarvester (via pip)..."
micromamba run -n practical3-env pip install structureHarvester || echo "structureHarvester not available on pip"

echo "Testing structure:"
structure --help || echo "structure not found"

echo "Testing clumpp:"
clumpp --help || echo "clumpp not found"

echo "Testing structureHarvester:"
structureHarvester --help || echo "structureHarvester not found"

echo "Done. Environment 'practical3-env' ready."