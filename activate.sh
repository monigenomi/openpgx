#!/bin/bash

direxists() {
    [ -d "$1" ]
}

direxists venv || python3 -m venv venv

source venv/bin/activate

pip install -e '.[dev]'