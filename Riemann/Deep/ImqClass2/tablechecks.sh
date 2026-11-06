#!/bin/bash
# tablechecks.sh: checks of the final tables (state.bin after the downloader stopped); a few seconds of python and gp.
cd "$(dirname "$0")"
K=$(od -A n -t u8 -j 16 -N 8 state.bin | tr -d ' ')
cp state.bin state_final.bin
./mjx --report state_final.bin final
{
echo "# tablechecks $(date '+%Y-%m-%d %H:%M'), tables k < $K"
echo "== rows N <= 5000 against ImqClass/lists_sha256.txt"
python3 rowcheck.py final_sha.txt final_perh.txt ../../Riemann/Deep/ImqClass/lists_sha256.txt
echo "== PARI: largest d of every row (rowcheck.gp)"
perl -e 'alarm 600; exec @ARGV' gp -q rowcheck.gp < /dev/null
echo "== PARI: largest d of 200 random h (perhcheck.gp)"
python3 perhcheck.py final_perh.txt && perl -e 'alarm 600; exec @ARGV' gp -q perhcheck.gp < /dev/null
echo "== Holmin-Jones-Kurlberg-McLeman-Petersen (GRH), odd h <= 9999"
python3 ../../Riemann/Deep/ImqClass2/hk_compare.py final_perh.txt class_group_orders.txt 9999
echo "== tier 4 of chunk x5 against the tables (overlap [2.8e10, K*2^28))"
python3 t4check.py runs/p_x5.t4 runs/x5.found dumps/tail.txt 10000 28000000000 30999999999 $K runs/append.txt
echo "== overlap: largest table d with h <= 10^4 and minimum h of the table files inside [2.8e10, K*2^28)"
awk 'NF==2 && $2<=10000' dumps/tail.txt | sort -n -k1,1 | tail -1
rg 'min h' mjx.log | awk -v K=$K '{split($1,a,"."); k=a[2]+0; sub(":","",a[2]); if (k>=105 && k<K) {for(i=1;i<=NF;i++) if($i=="min" && $(i+1)=="h") {if (m=="" || $(i+2)+0<m+0) {m=$(i+2); f=$1}}}} END {print "files k=105..", K-1, ": min h", m, "in", f}'
} > checks.txt 2>&1
cat checks.txt
