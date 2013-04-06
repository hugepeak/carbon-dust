file_name = 'c3_flows_time.eps'
restore, 'template12.dat'

buf = 0
position = [0.9, 0.85]
xlog = 1
ylog = 1

xrange = [1.e6,1.e8]
yrange = [1.e-22,1.e-19]

file0 = 'prop'
file1 = 'c3_flows1'
file2 = 'c3_flows2'
file3 = 'c3_flows3'

s0 = read_ascii( file0, template = template12 )
s1 = read_ascii( file1, template = template12 )
s2 = read_ascii( file2, template = template12 )
s3 = read_ascii( file3, template = template12 )

names = ['name', '$C + C_3 -> C_4 + \gamma$', '$C_3 + \gamma -> C_2 + C$', $
         '$C_3 + O -> CO + C_2$' ] 

xtitle = '$t - t_{6000K} (s)$'
ytitle = '$C_3 flows (reactions per second per atom)$'

d_small = 1.e-30

x1 = s1.b * 1.e9  
x1 = s1.a - s1.a[0] + d_small

y1 = s1.c + d_small
y2 = s2.c + d_small
y3 = s3.c + d_small

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

