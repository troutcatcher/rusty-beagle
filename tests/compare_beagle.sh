#!/bin/bash
# Validation harness: run Java Beagle and rusty-beagle on the same inputs
# and diff the decompressed output VCFs (ignoring ##filedate/##source).
#
# Usage: compare_beagle.sh <beagle.jar> <work_dir> [extra beagle args...]
set -u
JAR=$1
DIR=$2
shift 2
EXTRA_ARGS=("$@")

RUSTY=${RUSTY:-target/release/rusty-beagle}

norm() {
  zcat "$1" | grep -v -e '^##filedate=' -e '^##source='
}

echo "--- java beagle ---"
java -jar "$JAR" gt="$DIR/target.vcf.gz" ref="$DIR/ref.vcf.gz" \
    out="$DIR/java_out" "${EXTRA_ARGS[@]}" > "$DIR/java_stdout.txt" 2>&1
JAVA_RC=$?
if [ $JAVA_RC -ne 0 ]; then
  echo "java beagle FAILED (rc=$JAVA_RC)"; tail -20 "$DIR/java_stdout.txt"; exit 2
fi

echo "--- rusty-beagle ---"
"$RUSTY" gt="$DIR/target.vcf.gz" ref="$DIR/ref.vcf.gz" \
    out="$DIR/rust_out" "${EXTRA_ARGS[@]}" > "$DIR/rust_stdout.txt" 2>&1
RUST_RC=$?
if [ $RUST_RC -ne 0 ]; then
  echo "rusty-beagle FAILED (rc=$RUST_RC)"; tail -20 "$DIR/rust_stdout.txt"; exit 3
fi

norm "$DIR/java_out.vcf.gz" > "$DIR/java_out.vcf"
norm "$DIR/rust_out.vcf.gz" > "$DIR/rust_out.vcf"
if cmp -s "$DIR/java_out.vcf" "$DIR/rust_out.vcf"; then
  echo "IDENTICAL: $(wc -l < "$DIR/java_out.vcf") lines"
  exit 0
else
  echo "DIFFER:"
  diff "$DIR/java_out.vcf" "$DIR/rust_out.vcf" | head -30
  echo "..."
  diff "$DIR/java_out.vcf" "$DIR/rust_out.vcf" | wc -l
  exit 1
fi
