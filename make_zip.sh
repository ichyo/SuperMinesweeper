#!/bin/sh

set -eu

g++ -std=gnu++11 -O3 -o /dev/null ./submission.cpp # make sure same compile option works

DIR=$(mktemp -d)
cp submission.cpp $DIR/SuperMinesweeper.cpp
cd $DIR
zip submission.zip SuperMinesweeper.cpp
cd -
mv $DIR/submission.zip .
