file_name = 'c2_flows_T.jpeg'
restore, 'template12.dat'

buf = 1
position = [0.5, 0.85]
xlog = 0
ylog = 1

file1 = 'prop'
file2 = 'c2_flows'

s1 = read_ascii( file1, template = template12 )
s2 = read_ascii( file2, template = template12 )

names = ['name', '$C + C -> C_2 + \gamma$', '$CO + C -> C_2 + O$' ]

title = '$expand from T=6000(K) ' $
        + 'with n(C)=10^9(cm^{-3}) & n(O)=10^{10}(cm^{-3}), ' $
        + '\tau~100(days)$'
xtitle = '$t - t(6000K) (s)$'
xtitle = '$T(K)$'
ytitle = '$C_2 flows (reactions per second per atom)$'

xrange = [6000,1000]
yrange = [1.e-12, 1.e-8]

d_small = 1.e-30
tau = s2.a[0] 

x1 = s2.b * 1.e9  

y1 = s2.c + d_small
y2 = s2.d + d_small

p1 = plot( x1, y1, name = names[1], 'r', $
           xlog = xlog, ylog = ylog, $
           title = title, xtitle = xtitle, ytitle = ytitle, $
           xrange = xrange, yrange = yrange, buffer = buf )
p2 = plot( x1, y2, name = names[2], /overplot, 'b' )

l = legend( target = [p1,p2], $
            position = position ) 

p1.save, file_name

end

