file_name = 'c2_flows_T.jpeg'
restore, 'template12.dat'

buf = 1
position = [0.9, 0.85]
xlog = 0
ylog = 1

file0 = 'prop'
file1 = 'c2_flows1'
file2 = 'c2_flows2'
file3 = 'c2_flows3'
file4 = 'c2_flows4'

s0 = read_ascii( file0, template = template12 )
s1 = read_ascii( file1, template = template12 )
s2 = read_ascii( file2, template = template12 )
s3 = read_ascii( file3, template = template12 )
s4 = read_ascii( file4, template = template12 )

names = ['name', '$C + C -> C_2 + \gamma$', '$CO + C -> C_2 + O$', $
         '$C_2 + \gamma -> C + C$', '$C_2 + O -> CO + C$', $
         'total production', 'total destruction' ]

title = '$expand from T=6000(K) ' $
        + 'with n(C)=10^9(cm^{-3}) & n(O)=10^{10}(cm^{-3})$'
xtitle = '$T(K)$'
ytitle = '$C_2 flows (reactions per second per atom)$'

xrange = [6000,1000]
yrange = [1.e-12, 1.e-9]

d_small = 1.e-30

x1 = s1.b * 1.e9  

y1 = s1.c + d_small
y2 = s2.c + d_small
y3 = s3.c + d_small
y4 = s4.c + d_small
y5 = y1 + y2
y6 = y3 + y4

p1 = plot( x1, y1, name = names[1], 'r', $
           xlog = xlog, ylog = ylog, $
           title = title, xtitle = xtitle, ytitle = ytitle, $
           xrange = xrange, yrange = yrange, buffer = buf )
p2 = plot( x1, y2, name = names[2], /overplot, 'm' )
p3 = plot( x1, y3, name = names[3], /overplot, 'g' )
p4 = plot( x1, y4, name = names[4], /overplot, 'b' )
p5 = plot( x1, y5, name = names[5], /overplot, 'k' )
p6 = plot( x1, y6, name = names[6], /overplot, '2.' )

l = legend( target = [p1,p2,p3,p4,p5,p6], $
            position = position ) 

p1.save, file_name

end

