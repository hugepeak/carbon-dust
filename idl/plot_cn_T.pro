restore, 'template12.dat'

file1 = 'prop'
file2 = 'mass_fractions'
title = '$expand from T=6000(K), \rho=2.5\times 10^{-13}(g/cm^3), ' $
        + '\tau~215(days), ' $
        + 'with 10%C & 90%O$'

file_name = 'idl.pdf'
buf = 1
position = [0.3, 0.8]
xlog = 0
ylog = 1

xrange = [8000,1000]
yrange = [1.e-30, 1.e2]

s1 = read_ascii( file1, template = template12 )
s2 = read_ascii( file2, template = template12 )

names = ['species', 'O', 'C', '$CO$', $
         '$C_2$', '$C_3$', '$C_4$', '$C_5$', '$C_6$', $
         '$C_7$', '$C_8^C$', '$C_8^R$' ]

subtitle = '$T=T_0\times(\tau/t)^{1.8}, \rho=\rho_0\times(\tau/t)^3$'
xtitle = '$Temperature (K)$'
ytitle = '$Number fraction$'

d_small = 1.e-30

x1 = s1.c * 1.e9

y1 = s2.b + d_small
y2 = s2.c + d_small
y3 = s2.d / 2. + d_small
y4 = s2.e / 2. + d_small
y5 = s2.f / 3. + d_small
y6 = s2.g / 4. + d_small
y7 = s2.h / 5. + d_small
y8 = s2.i / 6. + d_small
y9 = s2.j / 7. + d_small
y10 = s2.k / 8.+ d_small
y11 = s2.l / 8.+ d_small

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

