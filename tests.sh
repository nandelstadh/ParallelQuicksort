#!/bin/bash

OUT="results.csv"

# overwrite previous results
echo "n,d,t,result" > "$OUT"

N=(1000000 10000000 100000000 1000000000)
D=(u n e)
T=(1 2 4 8)

for d in "${D[@]}"; do
for n in "${N[@]}"; do
for t in "${T[@]}"; do

    echo "Running n=$n d=$d t=$t"

    result=$(./main "$n" "$d" "$t")

    echo "$n,$d,$t,$result" >> "$OUT"

done
done
done
