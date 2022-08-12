#!/bin/bash

set -eo pipefail

[ -d venv ] || python3 -m venv venv

source venv/bin/activate

pip install --upgrade pip

pip install -e '.[dev]'
