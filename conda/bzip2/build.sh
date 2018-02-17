#!/bin/bash
set -e
make install PREFIX=$PREFIX CC="$CC" CFLAGS="$CFLAGS"
make -f Makefile-libbz2_so CC="$CC" CFLAGS="$CFLAGS"
ln -s libbz2.so.$PKG_VERSION libbz2.so
cp -d libbz2.so* $PREFIX/lib/
