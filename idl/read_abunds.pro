function read_abunds, input_file

openr, lu, /get, input_file

readf, lu, format = '(17x,i)', i_species

s_tmp =$
  create_struct($
    't', 0.d0,$
    'dt', 0.d0,$
    't9', 0.d0,$
    'rho', 0.d0,$
    'z', dblarr( i_species ),$
    'a', dblarr( i_species ),$
    'y', dblarr( i_species ),$
    'dy', dblarr( i_species ) $
  )

s = [s_tmp]

x = 0.d0

while not eof( lu ) do begin

  readf, lu, format = '(4x,D)', x
  s_tmp.t = x

  readf, lu, format = '(5x,D)', x
  s_tmp.dt = x

  readf, lu, format = '(5x,D)', x
  s_tmp.t9 = x

  readf, lu, format = '(13x,D)', x
  s_tmp.rho = x

  readf, lu, format = '(%"\n")'

  for i = 0, i_species - 1 do begin
    readf,$
      lu,$
      x1, x2, x3, x4
      s_tmp.z[i] = x1
      s_tmp.a[i] = x2
      s_tmp.y[i] = x3
      s_tmp.dy[i] = x4
  endfor

  readf, lu, format = '(%"\n")'
  readf, lu, format = '(%"\n")'
  readf, lu, format = '(%"\n")'

  s = [s,s_tmp]

endwhile

s = s[1:n_elements(s) - 1]

return, s

end
