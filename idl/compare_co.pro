file_name = 'compare_co.jpeg'
restore, 'template12.dat'

file1 = 'prop'
file2 = 'mass_fractions'
title = '$expand from T=6000(K), ' $ 
        + 'with n(C)=10^9(cm^{-3}) & n(O)=10^{10}(cm^{-3}), ' $
        + '\tau~100(days)$' 

buf = 1
position = [0.5, 0.8]
xlog = 0
ylog = 1

xrange = [6000,3000]
yrange = [1.e0, 1.e10]

s1 = read_ascii( file1, template = template12 )
s2 = read_ascii( file2, template = template12 )

names = ['names', '$Network Calculation$', $
         '$n_{CO}^{\gamma\ eq} = k_{CO}n_Cn_O\tau_\gamma$']

xtitle = '$Temperature (K)$'
ytitle = '$n_{CO} (cm^{-3})$'

d_small = 1.e-30
na = 1.1e10

x1 = s1.c * 1.e9

y1 = s1.d / s1.d[0] * na * s2.d / 2. + d_small
;y2 = 800. * ( x1 / 5.e3 )^(-22.6)
;y3 = 6.359e-26 * (x1)^(4.5) * exp(1.287e5/x1)
x2 = [6000.,5500.,5000.,4500.,4000.,3500.,3000]
y2 = [1.32e1,6.27e1,4.24e2,4.61e3,9.67e4,5.26e6,1.21e9]

x3 = [3500,3500]
y3 = yrange

p1 = plot( x1, y1, name = names[1], 'r', $
           xlog = xlog, ylog = ylog, $
           ;title = title, 
           xtitle = xtitle, ytitle = ytitle, $
           xrange = xrange, yrange = yrange, buffer = buf )
p2 = plot( x2, y2, name = names[2], /overplot, 'ko', sym_filled = 1 )
p3 = plot( x3, y3, /overplot, 'k--' )

t1 = text( 0.73, 0.75, font_size = 20, '$\tau_\gamma$' )
t2 = text( 0.82, 0.75, font_size = 20, '$\tau_e$' )

l = legend( target = [p1,p2], $
            position = position ) 

p1.save, file_name

end

