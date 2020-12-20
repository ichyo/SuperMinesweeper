#!/bin/zsh

set -eu

BASE=$(dirname "$0")

javac **/*.java
jar cfm $BASE/../tester/tester.jar Manifest.txt $BASE/**/*.class
