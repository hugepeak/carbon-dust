file_name = 'T_time.jpeg'
restore, 'template12.dat'

buf = 1
position = [0.7, 0.5]
xlog = 1
ylog = 0

file1 = 'prop'
file2 = 'prop2'

s1 = read_ascii( file1, template = template12 )

names = ['name', '$T(t) = 3800\cdot(100d\cdott^{-1})$', $
         '$T(t) = 18500\cdot(100d\cdott^{-1})^{1.8}$' ]

title = '$expand from T=6000(K), \rho=2.857\times 10^{-13}(g/cm^3), ' $
        + '\tau~100(days)$ '
xtitle = '$t - t(6000K) (s)$'
ytitle = '$T (K)$'

xrange = [1.e0,1.e8]
yrange = [0,7000]

d_small = 1.e-30

x1 = s1.b - s1.b[0] + d_small
y1 = s1.c * 1.e9

p1 = plot( x1, y1, name = names[1], 'r', $
           xlog = xlog, ylog = ylog, $
           xtitle = xtitle, ytitle = ytitle, $
           xrange = xrange, yrange = yrange, buffer = buf )

p1.save, file_name

end

