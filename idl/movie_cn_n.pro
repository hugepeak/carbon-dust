pro movie_cn_n, file = file

s1 = read_abunds( file )

buf = 1
xlog = 1
ylog = 1

xrange = [1.e0, 1.e6]
yrange = [1.e-20, 1.e1]

xtitle = '$Number of Carbon Atoms in Molecule$'
ytitle = '$Abundance$'

video_file = file + '.mp4'
video = idlffvideowrite(video_file)
framerate = 10
framedims = [640,512]
stream = video.addvideostream(framedims[0], framedims[1], framerate)

d_small = 1.e-30

for i = 0, n_elements(s1) - 1 do begin

  x1 = s1[i].z
  y1 = s1[i].y + d_small

  title = 'time=' + string( s1[i].t, format='(e7.1)' ) $
          + ', dt=' + string( s1[i].dt, format='(e7.1)' ) $
          + ', $T_9$=' + string( s1[i].t9, format='(e7.1)' ) $
          + ', $\rho$=' + string( s1[i].rho, format='(e7.1)' ) + '(g/$cm^3$)' 

  p1 = plot( x1, y1, name = '1', 'k', $
             xlog = xlog, ylog = ylog, $
             title = title, $
             xtitle = xtitle, ytitle = ytitle, $
             xminor = 0, $
             xrange = xrange, yrange = yrange, buffer = buf )

  timestamp = video.put(stream, p1.copywindow())

endfor

video.cleanup

print, 'File "' + video_file + '" written to current directory.'

end

