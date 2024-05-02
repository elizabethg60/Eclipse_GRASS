FUNCTION createsurface, nequator, incl, Omega, radius=radius

;+
; NAME:
;       createsurface
;
; PURPOSE:
;       This procedure creates a data structure containing all general information
;       about the star model like the coordinate grid of the stellar surface, the
;       mg and contrast intensity or the number of pixels on the surface. This
;       structure is then used by all of the following procedures.
;
; CALLING SEQUENCE:
;      star = createsurface(nequator, incl [, Omega, radius])
;      if Omega is not provided, parameters will be set to solar values (see projection.pro)
;
; INPUTS:
;       nequator:      Number of fields at the equator
;       incl:          Inclination angle
;       Omega:         Angular Velocity at Equator (deg/d)
;
; OPTIONAL INPUTS:
;       radius:        stellar radius (solar radii)
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       The function returns the coordinate grid of the surface in form of the following structure:
;       star={$
;            nequator:nequator       ,  $   ;# of points along equator defines resolution
;            incl:incl               ,  $   ;inclination angle under which the surface is seen
;            npoints:npoints         ,  $   ;# of points on the sphere (nearly the same size)
;                                           controlled by nequator
;            long:long               ,  $   ;Array for the longitude of each field
;            lat:lat                 ,  $   ;Array for the latitude of each field
;            area:area               ,  $   ;Array for the approximated area of each field
;            mu:mu                   ,  $   ;Array for the cos mu(norm with line of sight) of each field
;            aproj:aproj             ,  $   ;Array for the projected visible surface area of each field 
;                                           area * mu (negative if not visible)
;            numvis:numvis           ,  $   ;# of visible surface fields
;            vis:vis                 ,  $   ;Array for the indices of visible surface fields
;            dtheta:dtheta           ,  $   ;delta theta (180.D0/ntheta)
;            deltaphi:deltaphi       ,  $   ;
;            nphi:nphi               ,  $   ;
;            ntheta:ntheta           ,  $   ;# of areas on one side
;            mg:dblarr(npoints)      ,  $   ;
;            contrast:dblarr(npoints),  $   ;
;            dphi:dphi               ,  $   ;
;            veq:veq                 ,  $   ;rotation velocity at Eq.[km/s]
;            alpha:-0.2d0            ,  $   ;Differential rotating parameter for the model (default for sun)
;            epsilon:0.6d0           ,  $   ;default for sun
;            t0:5780.d0              ,  $   ;scalar value of temperature
;            dt:[500.d0, 5.d0]       ,  $   ;amplitude and scatter of T-variations 
;            temperature:5780.d0     ,  $   ;array of temperatures
;            evol:0.d0               ,  $   ;evolution parameter (see lightcurve.pro)
;            clong:0.d0              ,  $   ;face on longitude for map grid
;            rotation:0.d0           ,  $   ;to keep track of rotation of star
;            lamborig:679.0          ,  $   ;original lambda in nm (at which the continumm is calculated )
;            band:[400d0, 750d0]}           ;default band for Kepler in nanometers
;
; OPTIONAL OUTPUTS:
;
; PROCEDURE:
;
; EXAMPLE:
;       star = createsurface(10,80)
;       print, star.npoints
;
; MODIFICATION HISTORY:
;       2004-04-27:    Creation by AR
;       2006-04-26:    Adjusted by KS to create a structure
;       2020-12-10:    Modified by AR - added keyword radius
;-

; Check if the number of params is correct
If N_Params() LT 2 THEN BEGIN
    MESSAGE,/Info,'Syntax: star = createsurface(nequator, incl [, Omega, radius] )'
    RETAll
END

if not keyword_set(Omega) then begin
   Sun = 1
   Omega = 0
endif else begin
   Sun = 0
endelse

if not keyword_set(radius) then radius = 1.d

; Constants and some fields
DDTOR = !DPI/180.D0              ;double degree to radian
cDincl = DOUBLE(COS(DDTOR*incl)) ;cos of inclination
sDincl = DOUBLE(SIN(DDTOR*incl)) ;sin of inclination
ntheta = LONG(nequator/2)        ;# of areas on one side
thetacoord = DBLARR(ntheta)      ;create theta array
dtheta = 180.D0/ntheta           ;delta theta
nphi = INTARR(ntheta)           
dphi = DBLARR(ntheta)
phicoord = DBLARR(ntheta,nequator)

; Create a memory wasting grid
npoints = LONG(0)
FOR i=0, ntheta-1 DO BEGIN
    thetacoord[i] = -90.D0 + i*dtheta + 0.5*dtheta
    nphi[i] = FIX(nequator*SIN((thetacoord[i]+90.D0)*DDTOR))
    FOR j=0, nphi[i]-1 DO BEGIN
        npoints = LONG(npoints + 1)
        dphi[i] = 360.D0/nphi[i]
        phicoord[i,j] = j*dphi[i]
    ENDFOR
ENDFOR

; Now we do the real grid
lat  = DBLARR(npoints)
long = DBLARR(npoints)
area = DBLARR(npoints)
deltaphi = DBLARR(npoints)
k = LONG(0)

FOR i=0, ntheta-1 DO BEGIN
    FOR j=0, nphi[i]-1 DO BEGIN
        long[k]  = phicoord[i,j]
        lat[k]   = thetacoord[i]
        deltaphi[k]  = dphi[i]
        area[k]  = SIN((lat[k]+90.D0)*DDTOR)*dtheta*dphi[i]
        k = k+1
    ENDFOR
ENDFOR

;pointer because size may vary
vis = PTR_NEW(0,/allocate_heap)
numvis = LONG(0)

star={$
     nequator:nequator       ,  $   ;# points along equator defines resolution
     incl:incl               ,  $   ;inclination
     npoints:npoints         ,  $   ;# points on the sphere
     long:long               ,  $   ;longitudes
     lat:lat                 ,  $   ;latitudes
     area:area               ,  $   ;area of elements
     mu:DBLARR(npoints)      ,  $   ;cos(norm with line of sight)
     aproj:DBLARR(npoints)   ,  $   ;area * mu
     vproj:DBLARR(npoints)   ,  $   ;line-of-sight velocity
     xproj:DBLARR(npoints)   ,  $   ;projected x-coordinates
     yproj:DBLARR(npoints)   ,  $   ;projected y-coordinates
     zproj:DBLARR(npoints)   ,  $   ;projected z-coordinates
     ras:DBLARR(npoints)     ,  $   ;right ascension
     dec:DBLARR(npoints)     ,  $   ;declination
     vvec:DBLARR(npoints)    ,  $   ;barycentric velocity correction
     ecl:make_array(npoints,value=1.d)  ,$   ;intensity 0. for eclipsed stars
     airm:make_array(npoints,value=1.d) ,$   ;airmass over disc
     attn:make_array(npoints,value=1.d) ,$   ;attenuation from extinction over disc
     numvis:numvis           ,  $   ;number on visible side
     vis:vis                 ,  $   ;indices of visible elements
     dtheta:dtheta           ,  $   ;
     deltaphi:deltaphi       ,  $   ;
     nphi:nphi               ,  $   ;
     ntheta:ntheta           ,  $   ;
     mg:DBLARR(npoints)      ,  $   ;
     contrast:make_array(npoints,value=1.d) ,$   ;contrast can be used to create intensity patterns
     limbdark:DBLARR(npoints),  $   ;projected brightness due to limb darkening
     dphi:dphi               ,  $   ;
     veq:0.0                 ,  $   ;
     P:0.0                   ,  $   ;
     Omega:Omega             ,  $   ;
     Radius:Radius           ,  $   ;
     alpha:0.0               ,  $   ;
     epsilon:0.6d0           ,  $   ;default for sun
     t0:5780.d0              ,  $   ;default for sun
     dt:[500.d0, 5.d0]       ,  $   ;default 
     temperature:DBLARR(npoints),$  ;
     evol: 0.d0              ,  $   ;defaut is no evolution
     clong:   0.d0           ,  $   ;face on longitude for map grid
     rotation:0.d0           }      ;to keep track of rotation of star

;calculate projection
projection, star, sun=sun

RETURN, star
END
