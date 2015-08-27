#!/bin/bash
pow=-1.75
function plot {
(
echo "set grid"
# echo "set xr [0.42:0.46]"
echo "set xr [0:1]"
# echo "set xtics 0.0005"
if [[ $1 == M ]]; then
#     echo "set key at graph 0.92,0.5 r" 
    echo "set key bottom right"
    echo "set yr [0:1]"
    echo "set label '\$M\$' at graph 0.01,0.96 l"
elif [[ $1 == chi ]]; then
    echo "set key at graph 0.98,0.95 r" 
    echo "set yr [0:0.06]"
    echo "set label '\$L^{$pow}\\chi\$' at graph 0.005,0.94 l"
fi
echo "p 0 t ''"
echo 'set pointsize 0.5'
echo "f(x)=x"
# echo "f(x)=exp(-4*x)"
lc=1
lt=1
# for i in 1024.txt 
for i in $(ls -v result)
do
    test "$(ls result/$i | wc -l )" -gt 2 || continue
    if [[ $1 == M ]]; then
        echo "rep '< cat result/$i/*' u (f(\$3)):6:7 w e t '$i' lc $lc lt $lt"
        echo "rep '' u (f(\$3)):6 smooth csplines t '' lc $lc dashtype (5,10)"
    elif [[ $1 == chi ]]; then
        echo "rep '< cat result/$i/*' u (f(\$3)):(\$8*\$2**$pow):(\$9*\$2**$pow) w e t '$i' lc $lc lt $lt"
        echo "rep '' u (f(\$3)):(\$8*\$2**$pow) smooth csplines t '' lc $lc dashtype (5,10) "
#         echo "set yr [0:*]"
    fi
    let lc++
    let lt++
done
echo "###MAKEPDF font ',8'"
) > tmp.gpi
gplpdf tmp.gpi > /dev/null
rm tmp.gpi
mv tmp.pdf $1.pdf
}

plot M
plot chi
