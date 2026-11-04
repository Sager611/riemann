#!/bin/sh
# Fetch the three data tables of the audited preprint (single files only) into $CIRC2D_SRC.
P='preprints/A-sharp-Fourier-certificate-for-planar-circle-packing-September-23-2026'
D="${CIRC2D_SRC:-$(cd "$(dirname "$0")/../../.." && pwd)/.scratch/circ2d/src}"
mkdir -p "$D/build/data"
for f in coefficients.tex W_plus.tex W_minus.tex; do
  curl -s -L -o "$D/build/data/$f" "https://raw.githubusercontent.com/openai/math/main/$P/build/data/$f"
done
