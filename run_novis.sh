#!/bin/sh

set -eu

SEED=1
EXEC=./submission

while getopts "s:" "flag"; do
    case $flag in
        s) SEED="${OPTARG}";;
    esac
done
shift $((OPTIND-1))

make

java -jar tester/tester.jar -exec $EXEC -seed $SEED -delay 0 -novis -noSummary -debug -saveAll /tmp/novis/io -screen 2
