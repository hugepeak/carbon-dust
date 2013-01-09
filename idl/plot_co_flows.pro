file_name = 'co_flows_T.jpeg'
restore, 'template12.dat'

buf = 1
position = [0.85, 0.8]
xlog = 0
ylog = 1

file1 = 'prop'
file2 = 'co_flows'

s1 = read_ascii( file1, template = template12 )
s2 = read_ascii( file2, template = template12 )

names = ['name', '$CO + e^- -> CO + e^-$', $
         '$CO + \gamma -> C + O$', $
         '$C + O -> CO + \gamma$']

title = '$expand from T=6000(K) ' $
        + 'with n(C)=10^9(cm^{-3}) & n(O)=10^{10}(cm^{-3}), ' $
        + '\tau~100(days)$'
xtitle = '$t - t(6000K) (s)$'
xtitle = '$T(K)$'
ytitle = '$Flows (reactions per second per atom)$'

xrange = [6000,1000]
yrange = [1.e-10, 1.e-7]

d_small = 1.e-30
tau = s2.a[0] 

x1 = s2.b * 1.e9

y1 = s2.c + d_small
y2 = s2.d + d_small
y3 = s2.e + d_small

p1 = plot( x1, y1, name = names[1], 'r', $
           xlog = xlog, ylog = ylog, $
           title = title, xtitle = xtitle, ytitle = ytitle, $
           xrange = xrange, yrange = yrange, buffer = buf )
p2 = plot( x1, y2, name = names[2], /overplot, 'b' )
p3 = plot( x1, y3, name = names[3], /overplot, 'g' )

l = legend( target = [p1,p2,p3], $
            position = position ) 

p1.save, file_name

end

