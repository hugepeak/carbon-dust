file_name = 'c3.eps'
restore, 'template12.dat'

file1 = 'prop'
file2 = 'mass_fractions'
title = '$expand from T=6000(K), ' $ 
        + 'with n(C)=10^9(cm^{-3}) & n(O)=10^{10}(cm^{-3})$'
title = '$Network: C, O, CO, C_2 and C_3 without photodissociation$'

buf = 1
position = [0.85, 0.5]
xlog = 0
ylog = 1

xrange = [0,50]
yrange = [1.e-11,1.e1]

s1 = read_ascii( file1, template = template12 )
s2 = read_ascii( file2, template = template12 )

names = ['species', 'O', 'C', 'CO', $
         '$C_2$', '$C_3$', '$10^{15}\timesC_4$', $
         '$10^{17}\timesC_5$', '$10^{15}\timesC_6$', $
         '$10^{18}\timesC_7$', '$10^{23}\timesC_8^C$', '$10^{16}\timesC_8^R$' ]

subtitle = '$T=T_0\times(\tau/t)^{1.8}, \rho=\rho_0\times(\tau/t)^3$'
xtitle = '$t - t(6000K) (days)$'
ytitle = 'Abundance'

d_small = 1.e-30
x1 = s1.b - s1.b[0] + d_small 
x1 = x1/8.64e5
na = 1.1e10

y1 = s2.b + d_small
y2 = s2.c + d_small
y3 = s2.d / 2. + d_small
y4 = s2.e / 2. + d_small
y5 = s2.f / 3. + d_small

p1 = plot( x1, y1, name = names[1], 'k', $
           xlog = xlog, ylog = ylog, $
           ;title = title, $ 
           xtitle = xtitle, ytitle = ytitle, $
           xrange = xrange, yrange = yrange, buffer = buf )
p2 = plot( x1, y2, name = names[2], /overplot, '--' )
p3 = plot( x1, y3, name = names[3], /overplot, '-.' )
p4 = plot( x1, y4, name = names[4], /overplot, '-:' )
p5 = plot( x1, y5, name = names[5], /overplot, ':' )

l = legend( target = [p1,p2,p3,p4,p5], $
            position = position ) 

p1.save, file_name

end

