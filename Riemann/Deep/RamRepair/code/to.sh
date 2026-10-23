#!/bin/bash
# usage: to.sh SECONDS cmd...   (hard timeout via perl alarm)
s=$1; shift
exec perl -e 'alarm shift; exec @ARGV' "$s" "$@"
