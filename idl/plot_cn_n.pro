file_name = 'cn_n.jpeg'
restore, 'template12.dat'

file1 = 'prop'
file2 = 'mass_fractions'

buf = 1
position = [0.8, 0.8]
xlog = 0
ylog = 0

xrange = []
yrange = []

s1 = read_ascii( file1, template = template12 )
s2 = read_ascii( file2, template = template12 )

names = ['names', 'T=2031(K)', 'T=1802(K)', '$CO$', $
         '$C_2$', '$C_3$', '$C_4$', '$C_5$', '$C_6$', $
         '$C_7$', '$C_8^C$', '$C_8^R$' ]

xtitle = '$Number of Carbon Atoms$'
ytitle = '$log Abundance$'

i_212 = 214
i_213 = 216

x1 = [3.,4.,5.,6.,7.,8.]
y1 = [s2.f(i_212),s2.g(i_212),s2.h(i_212),s2.i(i_212),s2.j(i_212),s2.k(i_212)]
y2 = [s2.f(i_213),s2.g(i_213),s2.h(i_213),s2.i(i_213),s2.j(i_213),s2.k(i_213)]

y1 /= x1
y2 /= x1

yy1 = alog10(y1)
yy2 = alog10(y2)

p1 = plot( x1, yy1, name = names[1], 'k', $
           xlog = xlog, ylog = ylog, $
           ;title = title, 
           xtitle = xtitle, ytitle = ytitle, $
           xminor = 0, $
           xrange = xrange, yrange = yrange, buffer = buf )
p2 = plot( x1, yy2, name = names[2], 'r', /overplot )

l = legend( target = [p1,p2], position = position ) 

p1.save, file_name

print, 'file ' + file_name + ' written to current directory' 

end

