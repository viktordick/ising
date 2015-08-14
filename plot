#!/bin/bash
function plot {
(
echo "set key at graph 0.1,0.95 l" 
echo "set xr [0.419:0.451]"
if [[ $1 == M ]]; then
    echo "set yr [0:1]"
    echo "set label '\$M\$' at graph 0.01,0.92 l"
elif [[ $1 == chi ]]; then
    echo "set yr [0:300]"
    echo "set label '\$\\chi\$' at graph 0.01,0.95 l"
fi
echo "p 1/0 t ''"
echo 'set pointsize 0.5'
lc=1
lt=1
for i in $(ls result)
do
    test "$(cat result/$i | wc -l )" -gt 2 || continue
    extent=${i%.*}
    if [[ $1 == M ]]; then
        echo "rep 'result/$extent.txt' u 3:6:7 w e t '$extent' lc $lc lt $lt"
    elif [[ $1 == chi ]]; then
        echo "rep 'result/$extent.txt' u 3:(\$8/\$2**0):(\$9/\$2**0) w e t '$extent' lc $lc lt $lt"
    fi
    let lc++
    let lt++
done
) > tmp.gpi
gplpdf tmp.gpi > /dev/null
rm tmp.gpi
mv tmp.pdf $1.pdf
}

plot M
plot chi
