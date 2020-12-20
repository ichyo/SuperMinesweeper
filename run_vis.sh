#!/bin/sh

set -eu

SEED=1
EXEC=./submission

PAUSE=""

while getopts "ps:" "flag"; do
    case $flag in
        p) PAUSE="-pause";;
        s) SEED="${OPTARG}";;
    esac
done
shift $((OPTIND-1))

make

java -jar tester/tester.jar $PAUSE -exec $EXEC -seed $SEED -delay 3 -noSummary -saveAll /tmp/vis/io -screen 2
