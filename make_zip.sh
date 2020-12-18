#!/bin/sh

set -eu

DIR=$(mktemp -d)
cp submission.cpp $DIR/SuperMinesweeper.cpp
cd $DIR
zip submission.zip SuperMinesweeper.cpp
cd -
mv $DIR/submission.zip .
