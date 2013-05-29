file_name = 'cn_time.eps'
restore, 'template12.dat'

file1 = 'prop'
file2 = 'mass_fractions'
title = '$expand from T=6000(K), ' $ 
        + 'with n(C)=10^9(cm^{-3}) & n(O)=10^{10}(cm^{-3}), ' $
        + '\tau~100(days)$'

buf = 1
position = [0.3, 0.75]
xlog = 0
ylog = 1

xrange = [8.e6,1.4e7]
yrange = [1.e-15, 10]

s1 = read_ascii( file1, template = template12 )
s2 = read_ascii( file2, template = template12 )

names = ['species', '1', '10', '$10^2$', $
         '$10^3$', '$10^4$', '$10^5$', '$10^6$', $
         '$10^7$', '$10^8$', '$10^9$' ]

subtitle = '$T=T_0\times(\tau/t)^{1.8}, \rho=\rho_0\times(\tau/t)^3$'
xtitle = '$t (s)$'
ytitle = 'Mass fraction'

d_small = 1.e-30
x1 = s1.b 
na = 1.1e10

y1 = na * s1.d * s1.d[0]^(-1) * s2.b + d_small
y2 = na * s1.d * s1.d[0]^(-1) * s2.c + d_small
y3 = na * s1.d * s1.d[0]^(-1) * s2.d / 2. + d_small
y4 = na * s1.d * s1.d[0]^(-1) * s2.e / 2. + d_small
y5 = na * s1.d * s1.d[0]^(-1) * s2.f / 3. + d_small
y6 = na * s1.d * s1.d[0]^(-1) * s2.g / 4. + d_small
y7 = na * s1.d * s1.d[0]^(-1) * s2.h / 5. + d_small
y8 = na * s1.d * s1.d[0]^(-1) * s2.i / 6. + d_small
y9 = na * s1.d * s1.d[0]^(-1) * s2.j / 7. + d_small
y10 = na * s1.d * s1.d[0]^(-1) * s2.k / 8.+ d_small

y1 = s2.b + d_small
y2 = s2.c + d_small
y3 = s2.d + d_small
y4 = s2.e  + d_small
y5 = s2.f + d_small
y6 = s2.g + d_small
y7 = s2.h + d_small
y8 = s2.i + d_small
y9 = s2.j + d_small
y10 = s2.k+ d_small

p1 = plot( x1, y1, name = names[1], 'c', $
           xlog = xlog, ylog = ylog, $
           xtitle = xtitle, ytitle = ytitle, $
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
;p11 = plot( x1, y11, name = names[11], /overplot, '2k' )

l = legend( target = [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10], $
            position = position ) 

p1.save, file_name

end

