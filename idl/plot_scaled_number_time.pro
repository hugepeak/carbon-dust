restore, 'template12.dat'

file1 = 'prop'
file2 = 'mass_fractions'
title = '$expand from T=6000(K) ' $
        + 'with n(C)=10^9(cm^{-3}) & n(O)=10^{10}(cm^{-3}), ' $
        + '\tau~215(days)$'

file_name = 'idl.pdf'
buf = 0
position = [0.9, 0.8]
xlog = 1
ylog = 1

xrange = []
yrange = [1.e3,1.e17]

s1 = read_ascii( file1, template = template12 )
s2 = read_ascii( file2, template = template12 )

names = ['species', '$O$', '$C$', '$10^3\timesCO$', '$10^6\timesC_2$', $
         '$10^{12}\timesC_3$', '$10^{14}\timesC_4$', $
         '$10^{16}\timesC_5$', '$10^{17}\timesC_6$', $
         '$10^{20}\timesC_7$', '$10^{22}\timesC_8^C$', $
         '$10^{16}\timesC_8^R$' ]

xtitle = '$t - t(6000K) (s)$'
ytitle = '$Number density (cm^{-3})$' 

d_small = 1.e-30
na = 1.1e10
tau = 2.e7

x1 = s1.b - tau

y1 = na * s1.d * s1.d[0]^(-1) * s2.b + d_small
y2 = na * s1.d * s1.d[0]^(-1) * s2.c + d_small
y3 = na * s1.d * s1.d[0]^(-1) * s2.d * 1.e3 + d_small
y4 = na * s1.d * s1.d[0]^(-1) * s2.e * 1.e6 + d_small
y5 = na * s1.d * s1.d[0]^(-1) * s2.f * 1.e13 + d_small
y6 = na * s1.d * s1.d[0]^(-1) * s2.g * 1.e14 + d_small
y7 = na * s1.d * s1.d[0]^(-1) * s2.h * 1.e17 + d_small
y8 = na * s1.d * s1.d[0]^(-1) * s2.i * 1.e18 + d_small
y9 = na * s1.d * s1.d[0]^(-1) * s2.j * 1.e21 + d_small
y10 = na * s1.d * s1.d[0]^(-1) * s2.k * 1.e22 + d_small
y11 = na * s1.d * s1.d[0]^(-1) * s2.l * 1.e16 + d_small

p1 = plot( x1, y1, name = names[1], 'c', $
           xlog = xlog, ylog = ylog, $
           title = title, xtitle = xtitle, ytitle = ytitle, $
           xrange = xrange, yrange = yrange, buffer = buf )
p2 = plot( x1, y2, name = names[2], /overplot, 'b' )
p3 = plot( x1, y3, name = names[3], /overplot, 'g' )
p4 = plot( x1, y4, name = names[4], /overplot, 'r' )
p5 = plot( x1, y5, name = names[5], /overplot, 'c--' )
p6 = plot( x1, y6, name = names[6], /overplot, 'm--' )
p7 = plot( x1, y7, name = names[7], /overplot, 'y--' )
p8 = plot( x1, y8, name = names[8], /overplot, 'b--' )
p9 = plot( x1, y9, name = names[9], /overplot, 'g--' )
p10 = plot( x1, y10, name = names[10], /overplot, 'm' )
p11 = plot( x1, y11, name = names[11], /overplot, '2k' )

l = legend( target = [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11], $
            position = position ) 

p1.save, file_name

end

