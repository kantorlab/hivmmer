#!/bin/bash
set -e
./configure --prefix=$PREFIX CXX="$CXX -std=c++98 -Wimplicit-fallthrough=0"
make
make install
