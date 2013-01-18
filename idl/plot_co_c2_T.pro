file_name = 'co_c2_T.eps'
restore, 'template12.dat'

file1 = 'prop'
file2 = 'mass_fractions'
title = '$expand from T=6000(K), ' $ 
        + 'with n(C)=10^9(cm^{-3}) & n(O)=10^{10}(cm^{-3}), ' $
        + '\tau~100(days)$' 

buf = 1
position = [0.8, 0.33]
xlog = 0
ylog = 1

xrange = [6000,1000]
yrange = [1.e0,1.e7]

s1 = read_ascii( file1, template = template12 )
s2 = read_ascii( file2, template = template12 )

names = ['species', '$n(CO)$', '$10^5\timesn(C_2)$']

xtitle = '$Temperature (K)$'
ytitle = '$Number density (cm^{-3})$'

d_small = 1.e-30
na = 1.1e10

x1 = s1.c * 1.e9

y1 = na * s1.d / s1.d[0] * s2.d / 2. + d_small
y2 = 1.e5 * na * s1.d / s1.d[0] * s2.e / 2. + d_small

p1 = plot( x1, y1, name = names[1], 'r', $
           xlog = xlog, ylog = ylog, $
           ;title = title, 
           xtitle = xtitle, ytitle = ytitle, $
           xrange = xrange, yrange = yrange, buffer = buf )
p2 = plot( x1, y2, name = names[2], /overplot, 'b' )

l = legend( target = [p1,p2], $
            position = position ) 

p1.save, file_name

end

