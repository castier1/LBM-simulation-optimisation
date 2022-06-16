set terminal png size 800,500 enhanced font "Arial,12"
set output 'assembly.png'
set title "Histogramme du nombre d'occurrences par instructions"
set grid x
set grid y
set xtics ('add' 0, 'sub' 1, 'lea' 2, 'mov' 3, 'movb' 4, 'push' 5, 'pop' 6, 'addsd' 7, 'mulsd' 8, 'subsd' 9, 'movsd' 10, 'divsd' 11, 'addq' 12, 'movq' 13, 'callq' 14, 'callq' 15, 'retq' 16, 'jmpq' 17, 'movapd' 18, 'movupd' 19, 'movlpd' 20, 'movhpd' 21, 'addpd' 22, 'subpd' 23, 'mulpd' 24, 'divpd' 25) rotate by 45 right
set ylabel "Occurency"
set xlabel "Instructions"
#set linestyle 1 lt 2 lw 2
set style data histograms
plot "lbm0.data" with histogram title "v0", "lbm.data" with histogram title "v1"
