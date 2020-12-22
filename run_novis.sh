#!/bin/sh

set -eu

TLE=5000
SEED=1
EXEC=./submission

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

java -jar tester/tester.jar -exec $EXEC -seed $SEED -delay 0 -novis -noSummary -saveAll /tmp/novis/io -screen 2 -timeLimit $TLE $MORE_ARG
