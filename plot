#!/bin/bash
pow=-1.75
# pow=1
function plot {
(
echo "set samples 10000"
echo "betac = 0.5*atanh(1/sqrt(2))"
echo "set grid"
echo "set xr [0:1]"
echo "set xr [0.3:0.6]"
echo "set xr [0.42:0.45]"
# echo "set xr [0.436:0.443]"
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
echo 'set pointsize 0.3'
echo 'set bar 0.5'
echo "f(x)=x"
# echo "f(x)=exp(-4*x)"
lc=1
lt=1
for i in $(ls -v result)
# for i in 256
do
    test "$(ls result/$i | wc -l )" -gt 2 || continue
    if [[ $1 == M ]]; then
        echo "rep '< cat result/$i/*' u (f(\$3)):6:7 w e t '$i' lc $lc lt $lt"
#         echo "rep '' u (f(\$3)):6 smooth csplines t '' lc $lc dashtype (5,10)"
    elif [[ $1 == chi ]]; then
        echo "rep '< cat result/$i/*' u (f(\$3)):(\$8*\$2**$pow):(\$9*\$2**$pow) w e t '$i' lc $lc lt $lt"
#         echo "rep '' u (f(\$3)):(\$8*\$2**$pow) smooth csplines t '' lc $lc dashtype (5,10) "
#         echo "set yr [0:*]"
        echo "set arrow from betac,0 to betac,0.06 nohead dashtype (5,10) lc 0"
    fi
    let lc++
    let lt++
done
echo "###MAKEPDF font ',8'"
) > tmp-$1.gpi
gplpdf tmp-$1.gpi > /dev/null
rm tmp-$1.gpi
mv tmp-$1.pdf $1.pdf
}

plot M &
plot chi
