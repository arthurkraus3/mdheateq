set term png medium size 1280,1280 enhanced
set out 'OUT/MD_000050.png'
set grid xtics
set grid ytics
show grid
set title 'time: 5.000000e+00, T=8.150697e+00'
plot [0:3.000000e+01] [0:3.000000e+01] 'OUT/MD_000050.txt' u 1:2:(0.2) w circles fill solid t ''
set out
