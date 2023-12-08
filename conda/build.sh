#!/usr/bin/env bash
mkdir -p $PREFIX/bin

PHG2_DIR=$PREFIX/share/phg2-$PKG_VERSION-$PKG_BUILDNUM/

mkdir -p $PHG2_DIR

mv $SRC_DIR/* $PHG2_DIR

# Soft symlink to "point" to phg script, as a hard symlink
# leads to being unable to find the jars in lib/
# Helpful: https://stackoverflow.com/a/29786294
ln -s $PHG2_DIR/bin/phg $PREFIX/bin/phg2
