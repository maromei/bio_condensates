#!/bin/bash

set -e # exit on error

for file in .devcontainer/starting_scripts/*.sh; do
  bash "$file"
done

