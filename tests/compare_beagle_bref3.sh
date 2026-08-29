#!/bin/bash
# Validation harness for binary (.bref3) reference input: converts
# $DIR/ref.vcf.gz to bref3, then runs Java Beagle and rusty-beagle with
# ref=<the bref3 file> and diffs the decompressed output VCFs.
#
# A bref3 file's own block/sequence-coding boundaries are fixed at
# conversion time and are semantically relevant to marker clustering, so
# they need not match what a fresh read of the same data from VCF would
# choose -- imputing from ref.vcf.gz vs. imputing from its bref3 conversion
# can legitimately give different (equally valid) output in Java Beagle
# itself. The correct comparison is therefore same-file-vs-same-file (both
# programs reading the bref3 file this script builds), not the VCF suites'
# baseline. See docs/PORT_NOTES.md.
#
# Usage: compare_beagle_bref3.sh <beagle.jar> <bref3.jar> <work_dir> [extra beagle args...]
set -u
JAR=$1
BREF3_JAR=$2
DIR=$3
shift 3
EXTRA_ARGS=("$@")

RUSTY=${RUSTY:-target/release/rusty-beagle}

norm() {
  zcat "$1" | grep -v -e '^##filedate=' -e '^##source='
}

echo "--- converting ref.vcf.gz to bref3 ---"
zcat "$DIR/ref.vcf.gz" | java -jar "$BREF3_JAR" > "$DIR/ref.bref3" 2> "$DIR/bref3_stdout.txt"
BREF3_RC=$?
if [ $BREF3_RC -ne 0 ]; then
  echo "bref3 conversion FAILED (rc=$BREF3_RC)"; tail -20 "$DIR/bref3_stdout.txt"; exit 4
fi

echo "--- java beagle (ref=*.bref3) ---"
java -jar "$JAR" gt="$DIR/target.vcf.gz" ref="$DIR/ref.bref3" \
    out="$DIR/java_bref3_out" "${EXTRA_ARGS[@]}" > "$DIR/java_bref3_stdout.txt" 2>&1
JAVA_RC=$?
if [ $JAVA_RC -ne 0 ]; then
  echo "java beagle FAILED (rc=$JAVA_RC)"; tail -20 "$DIR/java_bref3_stdout.txt"; exit 2
fi

echo "--- rusty-beagle (ref=*.bref3) ---"
"$RUSTY" gt="$DIR/target.vcf.gz" ref="$DIR/ref.bref3" \
    out="$DIR/rust_bref3_out" "${EXTRA_ARGS[@]}" > "$DIR/rust_bref3_stdout.txt" 2>&1
RUST_RC=$?
if [ $RUST_RC -ne 0 ]; then
  echo "rusty-beagle FAILED (rc=$RUST_RC)"; tail -20 "$DIR/rust_bref3_stdout.txt"; exit 3
fi

norm "$DIR/java_bref3_out.vcf.gz" > "$DIR/java_bref3_out.vcf"
norm "$DIR/rust_bref3_out.vcf.gz" > "$DIR/rust_bref3_out.vcf"
if cmp -s "$DIR/java_bref3_out.vcf" "$DIR/rust_bref3_out.vcf"; then
  echo "IDENTICAL: $(wc -l < "$DIR/java_bref3_out.vcf") lines"
  exit 0
else
  echo "DIFFER:"
  diff "$DIR/java_bref3_out.vcf" "$DIR/rust_bref3_out.vcf" | head -30
  echo "..."
  diff "$DIR/java_bref3_out.vcf" "$DIR/rust_bref3_out.vcf" | wc -l
  exit 1
fi
