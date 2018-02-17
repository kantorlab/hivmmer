#!/bin/bash
set -e
sed -i~ 's/pear_CFLAG/#pear_CFLAG/' src/Makefile.am
autoreconf
./configure --prefix=$PREFIX
make
make install
