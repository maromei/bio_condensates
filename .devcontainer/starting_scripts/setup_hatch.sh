#!/bin/bash

# for some reason the first call always fails when building the environment
# --> redirect error and output to null. It will work fine once the ipykernel
# is installed.

hatch env create > /dev/null 2>&1
hatch run python -m ipykernel install --user --name bio_condensates --display-name "Python Bio-Condensates"
