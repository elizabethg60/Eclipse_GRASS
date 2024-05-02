pro specsurface, star, wave, spec, extcof, temperature=temperature, noplot=noplot, l0=l0, lstep=lstep, asymm=asymm, linemodel=linemodel, wavelength=wavelength, coeff=coeff, modwavint=modwavint, mumin=mumin, silent=silent, ipconf=ipconf

;procedure specsurface is part of the MODSTAR package
;
;+
; NAME:
;       specsurface
;
; PURPOSE:
;       This procedure computes a spectrum according to a stellar surface
;       defined in the structure "star". The structure follows the definitions
;       defined in the procedure createsurface.pro that belongs to the same
;       package. The procedure 
;
; INPUTS:
;       star:    structure defining stellar surface (see createsurface.pro)
;
; OUTPUTS:
;       wave:    wavelength array of output spectrum
;       spec:    output spectrum
;
; OPTIONAL INPUTS:
;       extcof: extinction coeffient - extinction will be applied, only
;                                      relevant if spatially resolved (solar)
;                                      observations are modelled
;       temperature:  temperature of star can be specified; usually not needed
;                     because the temperature is defined in star.temperature
;                                      
; KEYWORD PARAMETERS:
;       noplot   :  do not plot
;       silent   :  be quiet
;       l0       :  central wavelength in case a model lineprofile should be calculated
;       lstep    :  stepsize for model lineprofile
;       asymm    :  adds some asymmetry to model lineprofile
;       linemodel:  set keyword to 1 if a model lineprofile should be provided
;       wavelength: wavelength array for linemodel coefficients
;       coeff    :  linemodel coefficients (see cbspline2d for explanation)
;       ipconf   :  structure with interpolated spectral grid (see create_linesplines.pro)
;       modwavint:  set keyword if linemodel is used and coefficients are saved for index number x-axis instead of wavelengths
;       mumin    :  minimum for mu range for which the model spectra are defined
;
;  AR  17-Apr-2015 created
;  AR  05-Mar-2016 added line profile model
;  AR     Apr-2021 changed coeff logic, removed keyword mg
  
if n_params() lt 3 then begin
    message,/info,'specsurface, star, wave, spec, [temperature=, noplot=])'
    retall
endif

c = 299792.458d0

;set temperature if not defined
if keyword_set(temperature) then star.temperature=double(temperature) else temperature = 5780.d0

if keyword_set(linemodel) then begin
   tempwave = wavelength
   x = dindgen(n_elements(wavelength))
endif else begin
   lineprofile, l0, temperature, wave=tempwave, spec=tempspec, turb=5., step=lstep, noplot=noplot, depth=0.6 ;,/intnorm
   if keyword_set(asymm) then begin
      wingspec = (1. - shift(tempspec,n_elements(tempwave)/50.)) / 2.
      tempspec = tempspec - wingspec
   endif
   if not keyword_set(noplot) then plot,tempwave,tempspec
   tempspec = double(tempspec)
endelse

if not keyword_set(mumin) then mumin=0.

;calculate brightness: contrast, limb darkening, and projected surface
brightness = star.limbdark[*star.vis] * star.aproj[*star.vis]; * star.contrast[*star.vis] ;attenuation done in extspec

;initialisation
speccontrast = dblarr(n_elements(*star.vis), n_elements(tempwave))
temparrayleft = star.temperature[*star.vis] ;array of elements not yet calculated
ncind = 1

while ncind gt 0 do begin 
   ;next temperature
   tempdum = temparrayleft[0]
   ind = where(star.temperature[*star.vis] eq tempdum,nind)
   if not keyword_set(silent) then print,nind,' elements found with T(K) = ',tempdum

   ;find elements with other temperatures - bookkeeping
   dum = where(temparrayleft eq tempdum,complement=cind,ncomplement=ncind)
   if ncind gt 0 then temparrayleft = temparrayleft[cind] 
   
   if tempdum eq 0.d0 then begin
      speccontrast[ind,*] = 0.d0 
   endif else begin
      bb = planck(tempwave,tempdum)
      newwave = tempwave # (1.d0 - (star.vproj[*star.vis])[ind]/c)
      mugrid  = rebin(1#(star.mu[*star.vis])[ind] > mumin, n_elements(newwave[*,0]), nind) ; only use spectra for mu>mumin

      if keyword_set(linemodel) then begin
         if keyword_set(coeff) then begin
            ; presence of keyword coeff indicates spline model
            if keyword_set(modwavint) then dumwave = (newwave - tempwave[0]) / (tempwave[1] - tempwave[0]) else dumwave = newwave
            tempspec = cbspline2d(dumwave, mugrid, coeff, expol=1)
         endif else begin
            ; assume interpolated grid 
            x = tempwave
            x1 = newwave
            y1 = mugrid
            tempspec = interpolate(ipconf.zz, (x1-min(x))/(max(x)-min(x))*(n_elements(ipconf.zz[*,0])-1.),(y1-min(ipconf.mu))/(max(ipconf.mu)-min(ipconf.mu))*(n_elements(ipconf.zz[0,*])-1.))
         endelse
      endif else begin
         print, 'PANIC, only linemodel is implemented'
         retall
      endelse
      if keyword_set(extcof) then extspec,tempwave,tempspec,(star.airm[*star.vis])[ind],extcof ; extspec not yet tested with matrix
      speccontrast[ind,*] = transpose(tempspec * (bb#brightness[ind]))

   endelse
endwhile

;integration
spec = total(speccontrast, 1, /double)
wave = tempwave


END

