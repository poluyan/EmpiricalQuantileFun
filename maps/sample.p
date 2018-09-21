reset
#set termoption enhanced

set term pngcairo color dashed
set terminal png size 1000, 1000 font 'Open Sans,16'
set output 'sample.png'

##set term svg size 1000, 1000 font 'Roman,25' # empty background
#set term svg enhanced mouse size 320, 320 font 'Open Sans,14' # white backgraund
#set output '4.svg'

unset key

reset
xs = -0.2
ys = -0.7
f(x,y)= 3*((1-(x-xs))**2)*exp(-((x-xs)**2) - ((y-ys)+1)**2) - 10*((x-xs)/5 - (x-xs)**3 -  (y-ys)**5)*exp(-(x-xs)**2-(y-ys)**2) - 1/3*exp(-((x-xs)+1)**2 - (y-ys)**2) 
set xrange [-3:3]
set yrange [-3:3]
set isosample 1000, 1000
set table 'test.dat'
splot f(x,y)
unset table

set contour base
set cntrparam levels discrete 1.8
unset surface
set table 'cont.dat'
splot f(x,y)
unset table

reset
set xrange [-3:3]
set yrange [-3:3]
set border 3
unset key
set palette defined ( 1 1 1 1, 1 1 1 1 )
unset colorbox
set size square
set tics scale 0


set style line 104 lt 1 lc rgb 'red' pt 13 ps 1.5 lw 2.0

filter1(x,y) = f(x,y) > 1. ? x : 1/0
filter2(x,y) = f(x,y) > 1. ? y : 1/0

filter3(x,y) = f(x,y) > 1.8 ? x : 1/0
plot \
'cont1.dat' u (filter1($1,$2)):2 w filledcurves lc "gray",\
'cont2.dat' u 1:2 w filledcurves xy=0.6,1 lc "gray",\
'sample_explicit.dat' with points ls 104,\
'sampled_explicit.dat' with dots lc rgb 'black'


