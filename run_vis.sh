#!/bin/sh

set -eu

SEED=${1:-1}
EXEC=./submission

make

java -jar tester/tester.jar -exec $EXEC -seed $SEED -delay 3 -noSummary -debug -saveAll /tmp/vis/io -screen 2
