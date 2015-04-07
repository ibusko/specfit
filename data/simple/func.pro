function func, w, par

  if n_elements(par) mod 3 ne 0 then begin
    message, 'expected 3 parameters per Gaussian'
  endif
  ng = n_elements(par) / 3

  obs = 0.0d+00
  for ig=0, ng-1 do begin
    amp = par[3*ig]
    cen = par[3*ig+1]
    sig = par[3*ig+2]
    obs += amp*exp(-0.5d+00*((w-cen)/sig)^2)
  endfor

  return, obs

end
