#!/bin/sh
# dbnb2: replay Raymond's Lambda chain (Zenodo 23752103) restarted from Re s > 0.87338196661320752.
# Needs: uv, curl, ~1.5 GB scratch, ~45 min of one core. Run from an empty scratch directory.
# EV = path of Riemann/Deep/DBNB2/evidence (this directory).
set -eu
EV=${EV:?set EV to Riemann/Deep/DBNB2/evidence}
PY="uv run --no-project --python 3.12 --with python-flint==0.8.0 python3"
curl -s -L -o Proof_and_Replay.zip "https://zenodo.org/api/records/23752103/files/Proof_and_Replay.zip/content"
echo "bf5c7b3dca9e17e84e146ecc5ca0c6f671b409613088212e3f0b8a0b15f6aaf6  Proof_and_Replay.zip" | shasum -a 256 -c
mkdir -p pr && (cd pr && unzip -q ../Proof_and_Replay.zip)
B=pr/Proof_and_Replay/Lambda_Peer_Review
cp "$EV"/make_flow_b2.py "$EV"/patch_verifier.py "$EV"/float_check.py .
mkdir -p out
# 0. self-check: b0 = 7/8 must reproduce the author's published trajectory exactly
$PY make_flow_b2.py --workspace $B/inputs/workspace --out out/FLOW_REPRO78.json.gz --b0 7/8
# 1. headline flow (A) and conservative flow (B: retained row-55 tube kept by upward resets)
$PY make_flow_b2.py --workspace $B/inputs/workspace --out out/FLOW_B2.json.gz
$PY make_flow_b2.py --workspace $B/inputs/workspace --out out/FLOW_B2T.json.gz --tube
# 2. author's arc verifier and compact-witness script, only old-run constants replaced
python3 patch_verifier.py A out/FLOW_B2.json.gz verifyA
python3 patch_verifier.py B out/FLOW_B2T.json.gz verifyB
mkdir -p ws/dbnb2
ln -sfn "$PWD/$B/inputs/workspace/lambda_round14" ws/lambda_round14
ln -sfn "$PWD/$B/inputs/workspace/lambda_round15" ws/lambda_round15
ln -sfn "$PWD/out/FLOW_B2.json.gz" ws/dbnb2/FLOW_B2.json.gz
ln -sfn "$PWD/out/FLOW_B2T.json.gz" ws/dbnb2/FLOW_B2T.json.gz
for V in verifyA verifyB; do
  $PY $V/verify_arcs_b2.py --workspace ws
  $PY $V/compact_arc_witnesses_b2.py
done
# 3. independent float model (from the TeX formulas)
python3 float_check.py verifyA/ACTUAL_67_ROW_ARC_INPUTS.json
# 4. full 15-stage replay with the new arcs stage swapped in
mkdir -p $B/inputs/workspace/dbnb2 && cp out/FLOW_B2.json.gz $B/inputs/workspace/dbnb2/
cp verifyA/verify_arcs_b2.py $B/verification/arcs/verify_arcs_fresh.py
cp verifyA/compact_arc_witnesses_b2.py $B/verification/arcs/compact_arc_witnesses.py
SP=$($PY -c "import flint,os;print(os.path.dirname(os.path.dirname(flint.__file__)))")
W="$PWD/replayA"
(cd $B && $PY replay/run.py --all --jobs 1 --runtime "$SP" --work "$W")
