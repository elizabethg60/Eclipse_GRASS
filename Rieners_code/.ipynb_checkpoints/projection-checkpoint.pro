
pro projection, star, blsh_par=blsh_par, drot_par=drot_par, ld_par=ld_par, sun=sun

; procedure projection is part of the MODSTAR package
;
; calculates the variables in structure star that depend on the position of the observer
; the variables are: mu, aproj, vproj, vis
;
; INPUT PARAMETER
;            star:  structure for stellar grid information (created by createsurface.pro)
;
;
; OPTIONAL INPUT PARAMETERS
;            blsh_par:    convective blueshift
;            drot_par:    (differential) rotation
;            ld_par:      limb darkening
;
; OPTIONAL KEYWORDS:
;
;           sun:    if set, solar rotation law is used
;
; OUTPUT PARAMETER
;           star:   structure with updated stellar grid information 
;
;


Rsun = 6.95508d5              ;solar radius in km
clight = 299792.458

Rstar = star.radius * Rsun

;Some Constants
DDTOR = !DPI/180.D0 
cDincl = DOUBLE(COS(DDTOR*star.incl))
sDincl = DOUBLE(SIN(DDTOR*star.incl))

star.mu = cDincl*COS(DDTOR*(star.lat-90.D0))- $ ;recalc mu
  sDincl*SIN(DDTOR*(star.lat-90.D0))*COS(DDTOR*star.long)
star.aproj = star.mu * star.area  ;recalc aproj

;cartesian x and y-positions
phi    =  (star.long)*DDTOR
theta  =  (90.d0 - star.lat)  *DDTOR
rotang = -(90.d0 - star.incl) *DDTOR
;coordinates in system without inclination
xincl90 = SIN(theta)*SIN(phi)
yincl90 = COS(theta)
zincl90 = SIN(theta)*COS(phi)
;coordinates viewed under inclination ne 90.
star.xproj =  xincl90
star.yproj =  yincl90*COS(rotang) + zincl90*SIN(rotang)
star.zproj = -yincl90*SIN(rotang) + zincl90*COS(rotang)

;rotation
if keyword_set(sun) then begin
;solar differential rotation, Snodgrass, H., Ulrich, R. (1990), ApJ, 351, 309:
;for a summary of results, look into Beck, 2000
   if not keyword_set(drot_par) then begin
      mrad2deg = 86400.*180./(!Pi*1.e6)
;      drot_par = [2.972, -0.484, -0.361] * mrad2deg
      drot_par = [2.851, -0.343, -0.474] * mrad2deg
   endif
   A = drot_par[0]              ;     14.713
   B = drot_par[1]              ;     -2.396
   C = drot_par[2]              ;     -1.787
;convert sidereal to synodic rotation
   sidsyn = 0.9324
   A = A * sidsyn
   B = B * sidsyn
   C = C * sidsyn
;calculate P and v_eq
   star.P = 360.d0 / A
   star.veq = 2.*!Pi*RStar / star.P / 86400.d0
;differential rotation
   dlaw = (A + B*SIN(star.lat*DDTOR)^2. + C*SIN(star.lat*DDTOR)^4.)/A
endif else begin
   A = star.Omega
;calculate P and v_eq
   star.P = 360.d0 / A
   star.veq = 2.*!Pi*RStar / star.P / 86400.d0
;differential rotation
   dlaw = 1.d + star.alpha*SIN((star.lat)*DDTOR)^2.
endelse

star.vproj = star.veq*sDincl*dlaw*star.xproj
;relativistic effect
;star.vproj = star.vproj / sqrt(1.-(star.veq*sDincl*dlaw*COS(star.lat*DDTOR))^2./(clight*clight))

;convective blue shift
ind = where(star.temperature ge 0.99*median(star.temperature), nind, complement=cind, ncomplement=ncind)
if not keyword_set(blsh_par) then begin
   blsh_par = fltarr(4)
   blsh_par[0] = 0.;.35   ;amplitude
   blsh_par[1] = 0.; 0.5  ;mu symmetry center
   blsh_par[2] = 0.;10.   ;steepness
   blsh_par[3] = 0.;-0.05 ;y offset
endif
star.vproj[ind] = star.vproj[ind] + blsh_par[0] * (2.d0/(1.d0 + exp((star.mu[ind] - blsh_par[1])*blsh_par[2])) - 1.d0 ) + blsh_par[3]
;star.vproj[ind] = star.vproj[ind] - 0.15*star.mu[ind] + 0.105
;if nind gt 0 then star.vproj[ind] = star.vproj[ind] + .5


;Test:
;blueshift parametrized as 4th order
;polynomial - Cegla et al., 2014,
;parameters: priv. comm. M. Oshagh
;x = acos(star.mu[ind])/!Pi*180.
;vbs = ((4.564 + x*(0.2909 + x*(0.1076 + x*(-0.004107 + x*4.072e-5)))) / 1000. ) < 0.3
;star.vproj[ind] = star.vproj[ind] + vbs - 0.3


;meridional flow - prototype
;flowlat = -(abs(star.lat)<20.) * sign(star.lat) / 1000.
;star.vproj = star.vproj - flowlat * cos(star.long*DDTOR) * sin(star.lat*DDTOR)

;identify visible points
*star.vis = WHERE(star.mu GE 0.D0, numvis) ;remake the visible array
star.numvis = numvis              ;change numvis

;limb darkening
star.limbdark[*] = 0.d0
;star.limbdark[*star.vis] = (1.d0 - star.epsilon * (1.d0 - star.mu[*star.vis]))
if not keyword_set(ld_par) then begin
   ;Cox (2000) @ 600nm 
   ld_par = fltarr(3)
   ld_par[1] = .88 
   ld_par[2] =-.23 
   ld_par[0] = 1. - ld_par[1] - ld_par[2]
endif
star.limbdark[*star.vis] = poly(star.mu[*star.vis],ld_par)

end
