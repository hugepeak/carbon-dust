file_name = 'cn_n.jpeg'
restore, 'template12.dat'

file = 'idl'
file1 = 'prop'
file2 = 'mass_fractions'

buf = 1
position = [0.8, 0.8]
xlog = 0
ylog = 0

xrange = [3,8]
yrange = [-35,-15]

s1 = read_ascii( file1, template = template12 )
s2 = read_ascii( file2, template = template12 )

names = ['names', 'T=1953(K)', 'T=1824(K)', '$CO$', $
         '$C_2$', '$C_3$', '$C_4$', '$C_5$', '$C_6$', $
         '$C_7$', '$C_8^C$', '$C_8^R$' ]

xtitle = '$Number of Carbon Atoms$'
ytitle = '$log Abundance$'

video_file = file + '.mp4'
video = idlffvideowrite(video_file)
framerate = 10
framedims = [640,512]
stream = video.addvideostream(framedims[0], framedims[1], framerate)

x1 = [3.,4.,5.,6.,7.,8.]

for i = 1, 214 do begin

  y1 = [s2.f(i),s2.g(i),s2.h(i),s2.i(i),s2.j(i),s2.k(i)]

  y1 /= x1
  yy1 = alog10(y1)

  p1 = plot( x1, yy1, name = names[1], 'k', $
             xlog = xlog, ylog = ylog, $
             ;title = title, 
             xtitle = xtitle, ytitle = ytitle, $
             xminor = 0, $
             xrange = xrange, yrange = yrange, buffer = buf )

  timestamp = video.put(stream, p1.copywindow())

endfor

video.cleanup

print, 'File "' + video_file + '" written to current directory.'

end

