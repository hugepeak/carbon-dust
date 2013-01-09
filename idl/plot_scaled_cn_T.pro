file_name = 'scaled_cn_T.jpeg'
restore, 'template12.dat'

file1 = 'prop'
file2 = 'mass_fractions'
title = '$expand from T=6000(K) ' $
        + 'with n(C)=10^9(cm^{-3}) & n(O)=10^{10}(cm^{-3}), ' $
        + '\tau~100(days)$'

buf = 1
position = [0.34, 0.8]
xlog = 0
ylog = 1

xrange = [8000,1000]
yrange = [1.e-10, 1.e2]

s1 = read_ascii( file1, template = template12 )
s2 = read_ascii( file2, template = template12 )

names = ['species', '$O$', '$C$', '$10^3\timesCO$', '$10^{10}\timesC_2$', $
         '$10^{17}\timesC_3$', '$10^{18}\timesC_4$', $
         '$10^{20}\timesC_5$', '$10^{21}\timesC_6$', $
         '$10^{24}\timesC_7$', '$10^{26}\timesC_8^C$', $
         '$10^{21}\timesC_8^R$' ]

xtitle = '$Temperature (K)$'
ytitle = 'Abundances'

d_small = 1.e-30
na = 1.1e10

x1 = s1.c * 1.e9

y1 = s2.b + d_small
y2 = s2.c + d_small
y3 = s2.d * 1.e3 + d_small
y4 = s2.e * 1.e10 + d_small
y5 = s2.f * 1.e17 + d_small
y6 = s2.g * 1.e18 + d_small
y7 = s2.h * 1.e20 + d_small
y8 = s2.i * 1.e21 + d_small
y9 = s2.j * 1.e24 + d_small
y10 = s2.k * 1.e26 + d_small
y11 = s2.l * 1.e21 + d_small

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

