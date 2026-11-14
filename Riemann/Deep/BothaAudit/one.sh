#!/bin/bash
# one.sh (bothaaudit): one float-scan point: gen.py (Python, Arb) -> deig (C) even, odd -> post.py. Each step under memcap 290 MB + /usr/bin/time -l.
# usage: one.sh <out.jsonl> <c> <N> <prec> [k_even k_odd its [L_botha]]
cd "$(dirname "$0")"; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
X=../../../.scratch/bothaaudit; MC="bash ../../../.scratch/memcap.sh 290 1700"
out=$1; c=$2; N=$3; prec=$4; ke=${5:-3}; ko=${6:-2}; its=${7:-4}; Lb=${8:-}; tag=$X/g_${c//\//o}_$N
UV="uv run --no-project --with python-flint --with mpmath python3"
/usr/bin/time -l $MC $UV gen.py $c $N $prec $tag 2> $tag.gen.err || { cat $tag.gen.err; exit 1; }
/usr/bin/time -l $MC $X/deig $tag.even $prec $ke $its > $tag.even.out 2> $tag.even.err || { cat $tag.even.err; exit 1; }
/usr/bin/time -l $MC $X/deig $tag.odd $prec $ko $its > $tag.odd.out 2> $tag.odd.err || { cat $tag.odd.err; exit 1; }
rss() { rg -o '[0-9]+  maximum resident' $1 | awk '{printf "%.0f", $1/1048576}'; }
pk="{\"maxrss_MB\": [$(rss $tag.gen.err), $(rss $tag.even.err), $(rss $tag.odd.err)], \"L_botha\": \"$Lb\"}"
$UV post.py $out $c $N $prec $tag.even.out $tag.odd.out "$pk"
rm -f $tag.even $tag.odd
