#!/bin/sh

set -eu

SEED=${1:-1}
EXEC=./submission

make

java -jar tester/tester.jar -exec $EXEC -seed $SEED -delay 0 -novis -noSummary -debug -saveAll /tmp/novis/io -screen 2
