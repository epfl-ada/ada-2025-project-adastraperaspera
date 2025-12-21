#!/bin/bash
# Activate virtual environment and set environment variables
source ada-venv/bin/activate
export PYTHONPATH="$PWD/src"
export DYLD_LIBRARY_PATH="/opt/homebrew/opt/libomp/lib:$DYLD_LIBRARY_PATH"
echo "Virtual environment activated"
echo "PYTHONPATH set to: $PYTHONPATH"
echo "DYLD_LIBRARY_PATH set for XGBoost"
