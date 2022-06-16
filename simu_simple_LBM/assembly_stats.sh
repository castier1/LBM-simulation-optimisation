#!/bin/bash

# Clean the repository and generate the assembly code
make clean
rm -f lbm.data lbm.assembly
make
objdump -d lbm > lbm.assembly

# Find all the occurency
operations=("add" "sub" "lea" "mov" "movb" "push" "pop" "addsd" "mulsd" "subsd" "movsd" "divsd" "addq" "movq" "callq" "callq" "retq" "jmpq" "movapd" "movupd" "movlpd" "movhpd" "addpd" "subpd" "mulpd" "divpd")

for op in "${operations[@]}"; do
output=$(grep -wo "$op" "lbm.assembly" | uniq -c)
if [ -n "$output" ]; then
        echo "$output" >> "lbm.data"
else
    echo "    0 $op" >> "lbm.data"
fi
done

gnuplot assembly_plot.gp
