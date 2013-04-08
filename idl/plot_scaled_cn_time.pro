file_name = 'scaled_cn_time.eps'
restore, 'template12.dat'

file1 = 'prop'
file2 = 'mass_fractions'
title = '$expand from T=6000(K), ' $ 
        + 'with n(C)=10^9(cm^{-3}) & n(O)=10^{10}(cm^{-3})$'

buf = 1
position = [0.35, 0.75]
xlog = 1
ylog = 1

xrange = [1.e4,3.e8]
yrange = [1.e-8,1.e1]

s1 = read_ascii( file1, template = template12 )
s2 = read_ascii( file2, template = template12 )

names = ['species', 'O', 'C', '$10^3\timesCO$', $
         '$10^{10}\timesC_2$', '$10^{16}\timesC_3$', '$10^{15}\timesC_4$', $
         '$10^{17}\timesC_5$', '$10^{15}\timesC_6$', $
         '$10^{18}\timesC_7$', '$10^{23}\timesC_8^C$', '$10^{16}\timesC_8^R$' ]

subtitle = '$T=T_0\times(\tau/t)^{1.8}, \rho=\rho_0\times(\tau/t)^3$'
xtitle = '$t - t_{6000K} (s)$'
ytitle = 'Number density (1/cc)'
ytitle = 'Abundance'

d_small = 1.e-30
x1 = s1.b - s1.b[0] + d_small 
na = 1.1e10
fac = na * s1.d / s1.d[0]
fac = 1.

y1 = s2.b + d_small
y2 = s2.c + d_small
y3 = s2.d / 2. * 1.e3 + d_small
y4 = s2.e / 2. * 1.e10 + d_small
y5 = s2.f / 3. * 1.e16 + d_small
y6 = s2.g / 4. * 1.e15 + d_small
y7 = s2.h / 5. * 1.e17 + d_small
y8 = s2.i / 6. * 1.e15 + d_small
y9 = s2.j / 7. * 1.e18 + d_small
y10 = s2.k / 8. * 1.e23 + d_small
y11 = s2.l / 8. * 1.e16 + d_small

y1 *= fac
y2 *= fac
y3 *= fac
y4 *= fac
y5 *= fac
y6 *= fac
y7 *= fac
y8 *= fac
y9 *= fac
y10 *= fac
y11 *= fac

p1 = plot( x1, y1, name = names[1], 'c', $
           xlog = xlog, ylog = ylog, $
           ;title = title, $ 
           xtitle = xtitle, ytitle = ytitle, $
           xrange = xrange, yrange = yrange, buffer = buf )
p2 = plot( x1, y2, name = names[2], /overplot, 'b' )
p3 = plot( x1, y3, name = names[3], /overplot, 'g' )
p4 = plot( x1, y4, name = names[4], /overplot, 'r' )
p5 = plot( x1, y5, name = names[5], /overplot, 'c--' )
p6 = plot( x1, y6, name = names[6], /overplot, 'm--' )
p7 = plot( x1, y7, name = names[7], /overplot, 'k--' )
p8 = plot( x1, y8, name = names[8], /overplot, 'b--' )
p9 = plot( x1, y9, name = names[9], /overplot, 'g--' )
p10 = plot( x1, y10, name = names[10], /overplot, 'm' )
p11 = plot( x1, y11, name = names[11], /overplot, '2k' )

l = legend( target = [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11], $
            position = position ) 

p1.save, file_name

end

