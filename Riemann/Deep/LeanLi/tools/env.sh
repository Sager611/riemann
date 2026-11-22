# source me: sets LEAN (toolchain binary) and LEAN_PATH without lake (`lake env true` alone peaks at 661 MB RSS).
R=$(git rev-parse --show-toplevel)
export LEAN="$HOME/.elan/toolchains/leanprover--lean4---$(sed 's/.*://' "$R/lean-toolchain")/bin/lean"
LEAN_PATH="$R/.lake/build/lib/lean"
for p in "$R"/.lake/packages/*; do LEAN_PATH="$LEAN_PATH:$p/.lake/build/lib/lean"; done
export LEAN_PATH
