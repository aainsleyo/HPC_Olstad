#!/bin/bash

set -e

echo "Project 3: running analysis..."

# go to the directory where the script lives
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "Working directory: $SCRIPT_DIR"

# run python script
python3 plotting.py > output.log 2>&1

echo "Done."
echo "Outputs:"
echo " - output.log (stats + messages)"
echo " - .pdf and .png plots (if generated)"
