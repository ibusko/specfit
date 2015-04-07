pro fit

;Generate observations and return model parameters.
  obs, w, obs1g, obs1gn, obs2g, obs2gn, unc, par1, par2

;Default parameter settings.
  settings = {fixed:1, mpside:2}

;Open diagnostic output file.
  logfile = 'fit.log'
  openw, unit, logfile, /get_lun
  print, 'writing ' + logfile

;Loop through observations
  ids = ['obs1g','obs1gn','obs2g','obs2gn']
  nid = n_elements(ids)
  for iid=0, nid-1 do begin
    id = ids[iid]

;Get observation for current case
    if ~execute('obs='+id) then message, 'error getting observation for ' + id

;Parameter setup.
    if iid lt 2 then begin
      ng = 1
      parmod = par1
    endif else begin
      ng = 2
      parmod = [par1,par2]
    endelse
    noisy = iid mod 2 eq 1

;Loop through free/fixed parameter cases.
    for ifree=0, 2 do begin

;Offsets for initial parameter guess to force non-trivial solution.
      paroff = dblarr(3)
      paroff[0] = 0.1d+00 * parmod[0]
      if ifree ge 1 then paroff[1] = parmod[2]
      if ifree ge 2 then paroff[2] = 0.25d+00 * parmod[2]

;Initial parameter guess with free parameters offset.
      parini = parmod
      parini[0:ifree] += paroff[0:ifree]
      if ng eq 2 then parini[3:3+ifree] += paroff[0:ifree]

;Let some parameters be free.
      parinfo = replicate(settings, 3*ng)
      parinfo[0:ifree].fixed=0
      if ng eq 2 then parinfo[3:3+ifree].fixed=0
      fixstr = strjoin(strtrim(parinfo.fixed,2))	;for output

;Fit observation
      parfit = mpfitfun('func',w,obs,unc,parini $
                       , parinfo=parinfo,/quiet $
                       , perror=parunc,covar=covar $
                       , niter=niter,status=status)

;Print fit results
      if ifree eq 0 then begin
        if noisy then begin
          printf, unit
        endif else begin
          head = 'Case   It S  Fixed Row    Amp1      Cen1      Sig1' 
          if ng eq 2 then head +=    '      Amp2      Cen2      Sig2'
          if ng eq 1 then h='111' else h='111111'
          printf, unit
          printf, unit, head
          printf, unit, form='(a-6,a3,a2,a7,a4,6f10.6)' $
                      , 'model','  0',' x',h,'Mod',parmod
          printf, unit
        endelse
      endif 
      printf, unit, form='(a-6,i3,i2,a7,a4,6f10.6)' $
                  , id,niter,status,fixstr,'Ini',parini
      printf, unit, form='(18x,a4,6f10.6)','Fin',parfit
      if noisy then printf, unit, form='(18x,a4,6f10.6)','Unc',parunc

;End of loop through number of free parameters.
    endfor

;Print covariance array for case with all parameters floating.
    if noisy then begin
      pcor = covar / (parunc # parunc)
      for i=0, 3*ng-1 do begin
        printf, unit, form='(18x,a4,6f10.6)','Cor',pcor[*,i]
      endfor

;Plot fit through noisy data with all parameters free.
      psfile = id + '.ps'
      print, 'writing ' + psfile
      set_plot, 'ps'
      !p.font = 0
      device, file=psfile, bits=8, /color, /encap
      loadct, 12, /silent
      plot, w, obs, xsty=3, ysty=3, thick=5, psym=7 $
          , xtit='Wavelength (microns)', ytit='Intensity' $
          , xmarg=[7,1], ymarg=[3.5,0.5], xthi=5, ythi=5
      if ng eq 2 then begin
        oplot, w, func(w,parfit[0:2]),thi=5, line=5, col=10
        oplot, w, func(w,parfit[3:5]),thi=5, line=5, col=10
      endif
      oplot, w, func(w,parmod), thi=5, col=100
      oplot, w, func(w,parfit),thi=5, line=2, col=190
      xyouts, /norm, 0.12, 0.92, 'Model', col=100
      xyouts, /norm, 0.12, 0.88, 'Fit', col=190
      device, /close
      set_plot, 'x'
      !p.font = -1
    endif

;End of loop through observations.
  endfor
  free_lun, unit

end
