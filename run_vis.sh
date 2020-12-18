#!/bin/sh

set -eu

SEED=${1:-1}
EXEC=./submission

make

java -jar tester/tester.jar -exec $EXEC -seed $SEED -noSummary -debug -saveAll /tmp/vis/io
