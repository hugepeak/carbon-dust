function read_equil, input_file

openr, lu, /get, input_file

readf, lu, format = '(12x,i)', i_elements

i = 0

s_tmp =$
  create_struct($
    'yz_net', dblarr( i_elements + 1 ),$
    'yz_nse', dblarr( i_elements + 1 ),$
    'yz_qse', dblarr( i_elements + 1 ),$
    'yz_wse', dblarr( i_elements + 1 ) )

s = [s_tmp]

x = 0.d0

while not eof( lu ) do begin

  readf, lu, format = '(6x,D)', x
  s_tmp.time = x

  readf, lu, format = '(5x,D)', x
  s_tmp.t9 = x

  readf, lu, format = '(6x,D)', x
  s_tmp.rho = x

  readf, lu, format = '(5x,E)', x
  s_tmp.ye = x

  readf, lu, format = '(9x,E)', x
  s_tmp.ye_wse = x

  readf, lu, format = '(5x,E)', x
  s_tmp.yh = x

  readf, lu, format = '(9x,E)', x
  s_tmp.yh_nse = x

  readf, lu, format = '(9x,E)', x
  s_tmp.yh_wse = x

  for iz = 0, i_elements do begin
    readf,$
      lu,$
      i, x1, x2, x3, x4
      s_tmp.i[i] = i
      s_tmp.yz_net[i] = x1
      s_tmp.yz_nse[i] = x2
      s_tmp.yz_qse[i] = x3
      s_tmp.yz_wse[i] = x4
      ;print, i, x1, x2, x3, x4
  endfor

  readf, lu, format = '(%"\n")'

  s = [s,s_tmp]

endwhile

s = s[1:n_elements(s) - 1]

return, s

end
