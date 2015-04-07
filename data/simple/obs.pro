pro obs, w, obs1g, obs1gn, obs2g, obs2gn, unc, par1, par2

;Syntax
  if n_params() lt 8 then begin
    print, 'syntax: obs, w, obs1g, obs1gn, obs2g, obs2gn, unc, par1, par2'
    return
  endif

;NIRSpec G395H wavelength scale (approximately)
  wmin = 2.9d+00			;microns
  wmax = 5.0d+00			;microns
  disp = 0.000658d+00			;microns/pixel, Nyquist for R=3000
  nw = floor((wmax-wmin)/disp)
  w = wmin+disp*dindgen(nw)		;microns

;Gaussian 1 parameters
  amp1 = 1.0d+00			;arbitrary units
  cen1 = 4.0d+00			;microns, Brackett alpha
  sig1 = 4.0d-03			;microns, FWHM 

;Gaussian 2 parameters
  amp2 = 0.3d+00			;same units as amp1
  cen2 = cen1 + 0.008d+00		;microns, +0.008 micron offset
  sig2 = 0.006d+00			;microns, FWHM, 50% broader than sig1

;Extract wavelength around Gaussian
  iw = where(abs(w-cen1) le 10*sig1, nw)
  w = w[iw]

;Calculate 1 Gaussian
  obs1g = amp1*exp(-0.5d+00*((w-cen1)/sig1)^2)

;Calculate 1 Gaussian
  obs2g = obs1g + amp2*exp(-0.5d+00*((w-cen2)/sig2)^2)

;Add noise.
  unc = 1.0d-02					;same arbitrary units as amp
  obs1gn = obs1g + unc * randomn(1L, nw)	;fixed RNG seed
  obs2gn = obs2g + unc * randomn(2L, nw)	;fixed RNG seed

;Bundle parmeters for return argument.
  par1 = [amp1, cen1, sig1]
  par2 = [amp2, cen2, sig2]

;Write 1 Gaussian observation with no noise.
  parstr = ['amp','cen','sig'] + '=' + strtrim(string(par1, form='(f9.6)'),2)
  file = 'obs1g.dat'
  print, 'writing ' + file
  openw, unit, file, /get_lun
  printf, unit, '# Column 1: Wavelength [microns]'
  printf, unit, '# Column 2: Flux [arbitrary units]'
  printf, unit, '# No noise'
  printf, unit, '# Gaussian 1: ' + strjoin(parstr, ',  ')
  for iw=0L, nw-1 do begin
    printf, unit, form='(2f20.16)', w[iw], obs1g[iw]
  endfor
  free_lun, unit

;Write 1 Gaussian observation with noise.
  par1str = ['amp','cen','sig'] + '=' + strtrim(string(par1, form='(f9.6)'),2)
  file = 'obs1gn.dat'
  print, 'writing ' + file
  openw, unit, file, /get_lun
  printf, unit, '# Column 1: Wavelength [microns]'
  printf, unit, '# Column 2: Flux [arbitrary units]'
  printf, unit, '# Column 3: Flux uncertainty [same units as flux]'
  printf, unit, '# Gaussian 1: ' + strjoin(par1str, ',  ')
  for iw=0L, nw-1 do begin
    printf, unit, form='(3f20.16)', w[iw], obs1gn[iw], unc
  endfor
  free_lun, unit

;Write 2 Gaussian observation with no noise.
  par1str = ['amp','cen','sig'] + '=' + strtrim(string(par1, form='(f9.6)'),2)
  par2str = ['amp','cen','sig'] + '=' + strtrim(string(par2, form='(f9.6)'),2)
  file = 'obs2g.dat'
  print, 'writing ' + file
  openw, unit, file, /get_lun
  printf, unit, '# Column 1: Wavelength [microns]'
  printf, unit, '# Column 2: Flux [arbitrary units]'
  printf, unit, '# No noise'
  printf, unit, '# Gaussian 1: ' + strjoin(par1str, ',  ')
  printf, unit, '# Gaussian 2: ' + strjoin(par2str, ',  ')
  for iw=0L, nw-1 do begin
    printf, unit, form='(2f20.16)', w[iw], obs2g[iw]
  endfor
  free_lun, unit

;Write 2 Gaussian observation with noise.
  par1str = ['amp','cen','sig'] + '=' + strtrim(string(par1, form='(f9.6)'),2)
  par2str = ['amp','cen','sig'] + '=' + strtrim(string(par2, form='(f9.6)'),2)
  file = 'obs2gn.dat'
  print, 'writing ' + file
  openw, unit, file, /get_lun
  printf, unit, '# Column 1: Wavelength [microns]'
  printf, unit, '# Column 2: Flux [arbitrary units]'
  printf, unit, '# Column 3: Flux uncertainty [same units as flux]'
  printf, unit, '# Gaussian 1: ' + strjoin(par1str, ',  ')
  printf, unit, '# Gaussian 2: ' + strjoin(par2str, ',  ')
  for iw=0L, nw-1 do begin
    printf, unit, form='(2f20.16)', w[iw], obs2gn[iw]
  endfor
  free_lun, unit

end
