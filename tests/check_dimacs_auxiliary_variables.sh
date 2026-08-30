#!/usr/bin/env bash

binary=$1
instance=$2

output=$("$binary" \
    --vertices 2 \
    --no-SMS \
    --all-graphs \
    --hide-graphs \
    --dimacs "$instance" 2>&1)
status=$?

printf '%s\n' "$output"

if [[ $status -ne 20 ]]; then
    echo "expected SMS to exhaust the models (exit 20), got $status" >&2
    exit 1
fi

if ! grep -q "Number of graphs: 2" <<<"$output"; then
    echo "expected both assignments of the graph-edge variable" >&2
    exit 1
fi
