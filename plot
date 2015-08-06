#!/bin/bash
(
echo "set ytics nomirror"
echo "set y2tics"
echo "set key at graph 0.05,0.95 l" 
echo "set xr [0.419:0.451]"
# echo "set xr [0.295:0.455]"
# echo "set xr [0:1]"
echo "set yr [0:1]"
# echo "unset link y"
echo "p 1/0 t ''"
echo 'set pointsize 0.5'
echo "set label '\$\\chi\$' at graph 0.98,0.95 r"
lc=1
lt=1
for i in $(ls result)
do
    test "$(cat result/$i | wc -l )" -gt 2 || continue
    extent=${i%.*}
    #     echo "rep 'result/$extent.txt' u 3:6:7 w e t '$extent' lc $lc lt $lt"
    echo "rep 'result/$extent.txt' u 3:(\$8/\$2**0):(\$9/\$2**0) axes x1y2 w e t '$extent' lc $lc lt $lt"
    let lc++
    let lt++
    if [[ $lc == 9 ]]; then
        lc=1
        let lt++
    fi
done
# echo "rep 'result/016.txt' u 3:8:9 axes x1y1 w e t '016' lc 1 lt 1"

# echo "rep 'result/064.txt' u 3:8:9 axes x1y2 w e t '064'" lc 3
) >tmp.gpi 
gplpdf tmp.gpi > /dev/null
rm tmp.gpi
mv tmp.pdf plot.pdf
