#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Download files
wget
$DIR/src/kinfin.py "$@"
