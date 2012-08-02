#!/bin/bash

#min=$(for i in $(cat matrix.txt); do echo $i; done | sort -n | head -n100 | tail -n1)
#max=$(for i in $(cat matrix.txt); do echo $i; done | sort -n | tail -n2000 |head -n1)
min=$(for i in $(cat matrix.txt); do echo $i; done | sort -n | head -n50 | tail -n1)
max=$(for i in $(cat matrix.txt); do echo $i; done | sort -n | tail -n10 |head -n1)

echo "
set title \"Heat Map generated from a file containing Z values only\"
unset key
set tic scale 0

# Color runs from white to green
#set palette defined (0 0 1 0, 1 1 1 1)
set palette rgbformulae 7,5,15
#set palette model CMY rgbformulae 7,5,15
set cbrange [$min:$max]
set cblabel \"Score\"
unset cbtics

#set xrange [-0.5:50.0]
#set yrange [-0.5:100.0]

set yrange [-0.5:50.0]
set xrange [-0.5:100.0]


set view map
splot '-' matrix with image
$(cat matrix.txt)
e
e
" | gnuplot
