#!/bin/sh

set -eu

TL=7500 # msec
CASE=200

THREADS=6
BASE=1000000

SOURCE=./submission.cpp
EXEC=./submission

MORE_ARG=""

while getopts "n:t:d:" "flag"; do
    case $flag in
        n) CASE=${OPTARG};;
        t) THREADS=${OPTARG};;
        d) MORE_ARG="${MORE_ARG} -d ${OPTARG}";;
    esac
done
shift $((OPTIND-1))

OUTPUT=/tmp/MM/$(date '+%Y%m%d-%H%M%S')
[ -f "$OUTPUT" ] && { echo "$OUTPUT exists"; exit 1; }

make

mkdir -p "$OUTPUT"

echo "source_file: ${SOURCE}" | tee -a $OUTPUT/info.yml
echo "source_md5sum: $(md5sum "$SOURCE" | cut -d ' ' -f 1)" | tee -a $OUTPUT/info.yml
echo "executable_md5sum: $(md5sum "$EXEC" | cut -d ' ' -f 1)" | tee -a $OUTPUT/info.yml
echo "executable: ${EXEC}" | tee -a $OUTPUT/info.yml
echo "threads: ${THREADS}" | tee -a $OUTPUT/info.yml
echo "time_limit: ${TL}" | tee -a $OUTPUT/info.yml
echo "first_seed: ${BASE}" | tee -a $OUTPUT/info.yml
echo "case_num: ${CASE}" | tee -a $OUTPUT/info.yml
echo "start_time: $(date '+%Y-%m-%d %H:%M:%S')" | tee -a $OUTPUT/info.yml
echo "additional_args: ${MORE_ARG}"

echo ""
echo $OUTPUT

java -jar tester/tester.jar -exec $EXEC -seed $BASE+$CASE -novis -noSummary -timeLimit $TL -threads $THREADS -saveScores $OUTPUT/scores.txt -saveAll $OUTPUT/io $MORE_ARG > /dev/null

echo "end_time: $(date '+%Y-%m-%d %H:%M:%S')" | tee -a $OUTPUT/info.yml

python ./show_stats.py $OUTPUT
