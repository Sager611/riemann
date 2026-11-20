#!/bin/bash
# download + extract arXiv e-prints for grhharv4; writes <id>/ and <id>.hits
cd /Users/adriansager/playground/riemann/.scratch/grhharv4
for i in $(cat ids.txt); do
  [ -d "$i" ] && continue
  curl -s -L -A 'grhharv4-research' -o "$i.tgz" "https://arxiv.org/e-print/$i"
  mkdir -p "$i"
  t=$(file -b "$i.tgz")
  case "$t" in
    gzip*) if tar -xzf "$i.tgz" -C "$i" 2>/dev/null; then :; else gunzip -c "$i.tgz" > "$i/main.tex"; fi ;;
    *PDF*) echo "$i PDFONLY" ;;
    *tar*) tar -xf "$i.tgz" -C "$i" ;;
    *) cp "$i.tgz" "$i/main.tex" ;;
  esac
  rm -f "$i.tgz"
  find "$i" -type f \( -name '*.pdf' -o -name '*.png' -o -name '*.jpg' -o -name '*.eps' \) -delete
  rg -n -i --no-heading 'GRH|ERH|Riemann|zero-free|zero free|quasi' -g '*.tex' "$i" > "$i.hits" 2>/dev/null
  n=$(wc -l < "$i.hits"); f=$(find "$i" -name '*.tex' | wc -l)
  echo "$i texfiles=$f hits=$n"
  sleep 3
done
