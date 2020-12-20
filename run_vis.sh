#!/bin/sh

set -eu

SEED=1
EXEC=./submission

PAUSE=""

MORE_ARG=""

while getopts "ps:d:l:m:" "flag"; do
    case $flag in
        p) PAUSE="-pause";;
        s) SEED="${OPTARG}";;
        d) MORE_ARG="${MORE_ARG} -d ${OPTARG}";;
        l) MORE_ARG="${MORE_ARG} -n ${OPTARG}";;
        m) MORE_ARG="${MORE_ARG} -m ${OPTARG}";;
    esac
done
shift $((OPTIND-1))

make

java -jar tester/tester.jar $PAUSE -exec $EXEC -seed $SEED -delay 3 -noSummary -saveAll /tmp/vis/io -screen 2 $MORE_ARG
