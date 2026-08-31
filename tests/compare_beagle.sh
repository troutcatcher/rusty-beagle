#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
#
# rusty-beagle - a Rust port of Beagle 5.5 genotype phasing and imputation.
# Copyright (C) 2026 The rusty-beagle authors
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
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
