#!/bin/bash
set -e
./configure --prefix=$PREFIX CXX="$CXX -std=c++98"
make
make install
