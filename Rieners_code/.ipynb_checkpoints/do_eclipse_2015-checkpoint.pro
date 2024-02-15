
;goto, plot


;
;define star and setup
; loop over observation steps, not over rotation angle
;

;*** definition of parameters ****

;correction of computer time in seconds
;corrsec = 90.

;limb darkening
; ld_par_set = 0
;   ld_par_set  = [  0.30,   0.93, -0.23]       
;   ld_par_set = [ 0.131,  0.918, -0.049] ;where is this from ? close to 400nm in Cox (2000) ?
;   ld_par_set = [0.30505, 1.13123, -0.78604, 0.40560,  0.02297, -0.07880] ; Pierce & Slaughter, 1977 @ 5800nm
ld_par_set = [0.28392, 1.36896, -1.75998, 2.22154, -1.56074,  0.44630] ; Neckel & Labs, 1994 @ 5800 nm
;   ld_par_set = [0.20924, 1.30798, -1.20411, 1.12505, -0.67196,  0.14381] ; Neckel & Labs, 1994 @ 4930 nm
;   ld_par_set = [0.19571, 1.30551, -1.25845, 1.50626, -1.05472,  0.30570] ; Neckel & Labs, 1994 @ 4774 nm

;convective blueshift
;blsh_par_set = 0
;   blsh_par_set = [0.3714, 0.4711,10.327,-0.044]
;  blsh_par_set = [0.275, 0.375,5.,0.025] ;looks more or less like Fig.18 in Beeck et al., 2013
;   blsh_par_set = [0.27, 0.325, 12., 0.11] ;test
;   blsh_par_set = [0.275, 0.5, 30., 0.025]

;differential rotation
;   drot_par_set = 0
   drot_par_set = [14.713,-2.396,-1.787]


;*** Set 1 ***
;linemodel = 0
;corrsec = 90.
;blsh_par_set = [0.27, 0.325, 12., 0.11]

;*** Set 2 ***
linemodel = 1
corrsec = 0.
blsh_par_set = [0., 0., 0., 0.]


;*********************************


;hmi = 0
chunks = [0,45,124]
sigmaclip = 1.
lincorr = 1

ps = 1
newprofile  = linemodel
vradprofile = linemodel
barycalc    = 1
spot = 0
noplot = 1

addpath,'./eclipse_ephemeris'
c = 299792.458d
nn = 301 

; create spectra
;dv = 15.                        ; width of spectral window in km/s
dv = 30.
;wcen = 6175.0443
wcen = 6705.3
cd, '../FluxStar'
create_linesplines, wcen, dv=dv, /noplot, sunsave=sunsave
cd, '../Solar_Eclipse'

restore, '~/Data/playground/FluxStar/IAG_Limb_surface.sav'
ipconf = {x:x, mu:mu, zz:zz}
checkwave = x
xdata = x

;Sun template
;solar atlas
if not keyword_set(atlas_vis) then restore,'../Solar_Atlas/Solar_Atlas_VIS.sav'
wave_sun_vis  = 1.e8 / reform(atlas_vis.field1[0,*])
suntemplate_orig = reform(atlas_vis.field1[1,*])
ind = where(wave_sun_vis gt wcen*(1.-1.2*dv/c) and wave_sun_vis lt wcen*(1.+1.2*dv/c))
wave_sun_vis = wave_sun_vis[ind]
suntemplate_orig = suntemplate_orig[ind]
contf, suntemplate_orig, cont, sbin=6, nord=1, frac=0.01, /plot
suntemplate_norm = suntemplate_orig / cont

npe  = 100  ;number of points around equator
incl = 97.05  ;inclination
ecl  = 1                        ; ON/OFF simulate eclipse of 20-Mar-2015


star = createsurface(npe, incl)
star.contrast[*] = 1.d0
star.temperature[*] = 5780.d0

;put a spot on the surface
if keyword_set(spot) then spot,star


;initialize ephemiris
eph = createphemeris_eclipse_mar2015(corrsec)
maxecl = 0.

n = eph.nOPUS ; number of steps, defined by observation
velarrbary = dblarr(n)
velarrcheap = dblarr(n)
velarrgauss = dblarr(n)
velarrterra = dblarr(n)
velarrxcorr = dblarr(2,n)
stararr = replicate(star,n)
starvis = ptrarr(n,/allocate_heap)
cheapfluxarr = dblarr(n)

;extinction model
;extcof = replicate(0.8d1, n)
extcof = .4 - .3*(eph.obs_time-eph.obs_time[0])/(eph.obs_time[n-1]-eph.obs_time[0])
;extcof = extcof * 3.

if not keyword_set(barycalc) then begin
   print,'Using pre-tabulated barycentric corrections...'
   restore,'vvec_arr.sav'
endif else begin
   vvecarr = dblarr(n,star.npoints)
endelse
   
nwave=n_elements(checkwave)
if keyword_set(newprofile) then specarr = fltarr(n,nwave) else restore,'modelspectra.sav'

ld_par = ld_par_set
blsh_par = blsh_par_set
drot_par = drot_par_set
projection,star, blsh_par=blsh_par, drot_par=drot_par, ld_par=ld_par, /sun

if keyword_set(hmi) then begin
   hmi2surface,'hmi.V_720s.20150320_093600_TAI.1.Dopplergram.fits',hmivproj,star,hdr
   rsunrr = 963.499939d0 / 3600.d0 / 180.d0 * !DPi
   xv = 31642.8d0 & yv = 2421.3d0; - 30000.d0
   hmivproj = hmivproj + 1537.287744 + sin(rsunrr*star.xproj)*xv + sin(rsunrr*star.yproj)*yv
   hmi2surface,'hmi.Ic_720s.20150320_093600_TAI.1.continuum.fits',hmicont,star
   hmi2surface,'hmi.M_720s.20150320_093600_TAI.1.magnetogram.fits',hmimg,star
endif

;restore, '~/Kepler/MODSTAR/spec/lineprofile/G2_coeff.sav'
;restore, '~/Data/playground/Sun/IAG_Limb_6173.sav'
;restore, '~/Data/playground/Sun/ld_par_6713.sav' ; pfit

for i = 0, n-1 do begin

   ;initialize
   ld_par = ld_par_set
   blsh_par = blsh_par_set
   drot_par = drot_par_set

   ; model and standard parameters
;   projection,star, blsh_par=blsh_par, drot_par=drot_par, ld_par=ld_par, /sun
   ; use observed spectra that include blueshift, use standard limb darkening
   ;projection, star, /sun, blsh_par=replicate(0.d, 4), ld_par=ld_par, drot_par=drot_par
   ; use observed spectra that include blueshift, apply limb darkening fit
;   blsh_par = [0., 0., 0., 0.]
   projection, star, blsh_par=blsh_par, ld_par=pfit, drot_par=drot_par, /sun
;   projection, star, blsh_par=blsh_par, ld_par=ld_par, drot_par=drot_par, /sun

   if keyword_set(vignet) then vignetting, star, eph, i

   if keyword_set(hmi) then begin
      star.vproj[*] = !values.f_nan
      star.vproj = hmivproj/1000.
      star.limbdark[*] = !values.f_nan
      star.limbdark = (hmicont/max(hmicont,/nan))
      star.mg = hmimg
      ;indmag = where(abs(hmimg) lt 200,nmag,complement=complement)
      ;star.vproj[indmag] = star.vproj[indmag] - 0.0004*star.mg[indmag]
   endif

   if ecl eq 1 then begin
      eclipse, star, eph, i
      noneclvis = WHERE(star.mu GE 0.D0, numvis)
      fecl = (total((star.aproj)[noneclvis]) - total((star.aproj)[*star.vis]) ) / total((star.aproj)[noneclvis])
      if fecl gt maxecl then begin
         maxecl = fecl
         indmax = i
      endif
   endif

   extinction, star, eph, i, extcof[i]

   if keyword_set(barycalc) then begin
      bary_sun,star,eph,i 
      vvecarr[i,*] = star.vvec
   endif else begin
      star.vvec = reform(vvecarr[i,*])
   endelse
   star.vproj[*star.vis] = star.vproj[*star.vis] - star.vvec[*star.vis] ; from heliocentric to geocentric

   stararr[i,*] = star
   *starvis[i] = *star.vis

   if not keyword_set(noplot) then begin
      brightness = star.temperature*star.limbdark*star.attn
      
      wset,0
      plotsurface,star,star.vproj-median(star.vproj),lineco=!black,min=-veq,max=veq,colortable=37,/grid
      ;movie stuff
      image = tvrd(/true)
      if i eq 0 then images = fltarr(3, n_elements(image(0,*,0)), n_elements(image(0,0,*)), n)
      images(*,*,*,i) = image
      
      ;do this later to create movie:
      ;from IDL: image2pngs,images
      ;not from IDL:  ~/local/ffmpeg-git-20210528-i686-static/ffmpeg -i frame%4d.png -pix_fmt yuv420p output.mp4

      ;plotsurface,star,star.vproj,lineco=!black,min=-.02,max=.02,colortable=3,/grid
      wset,1
      plotsurface,star,brightness,lineco=!white,min=0.,max=max((brightness)[*star.vis]),colortable=3
      ;plotsurface,star,hmimg,lineco=!white,min=-10.,max=10.,colortable=1
      ;plotsurface,star,star.mu,lineco=!white,min=0.,max=1,colortable=37,/grid
      ;plotsurface,star,star.attn,lineco=!white,min=min((star.attn)[*star.vis]),max=max((star.attn)[*star.vis]),colortable=3;,/grid
   endif

   if keyword_set(newprofile) then begin
      ;profile calculation
      ;MuRaM model
      ;specsurface, star, wave, spec, extcof[i], /noplot, /linemodel, wavelength=wavelength[*,1], coeff=coeff[0,1], /modwavint
      ;IAG observations
      ;specsurface, star, wave, spec, extcof[i], /noplot, /linemodel, wavelength=checkwave, coeff=coeffiag
      specsurface, star, wave, spec, /noplot, /linemodel, wavelength=checkwave, ipconf=ipconf
      specarr[i,*] = spec
      plot, wave_sun_vis, suntemplate_norm, /nodata
      oplot, wave, spec/max(spec)
      oplot, wave_sun_vis, suntemplate_norm, col=!RED
   endif else begin
      spec = specarr[i,*]
   endelse

   
   if keyword_set(vradprofile) then begin
      spec=spec/max(spec)
      velwave = (wave-wcen)/wcen * c
      ;wait = get_kbrd()
      ;velocity calculation
      ;1. Gauss fit
;      a = [.5,.0,1.5]
;      fit = gaussfit(velwave,1.d0-spec,a,nterms=n_elements(a),estimates=a)
;      velarrgauss[i] = a[1]
      ;2. barycenter of the line
      velarrbary[i] = total(velwave*(1.d0-spec),/double) / total(1.d0-spec,/double)
      ;3. xcorr velocity
;      vrad, checkwave, checkspec, wave, spec, radvel, vchimax=vchimax, range=3.
;      velarrxcorr[0,i] = -radvel
;      velarrxcorr[1,i] = -vchimax
      ;4. terra
      terra, wave_sun_vis, suntemplate_norm, wave, spec/max(spec), radvel
      velarrterra[i] = radvel
   endif
   ;4. weighted velocity shifts
   cheapspec, star, vel, cheapflux
   velarrcheap[i] = vel
   cheapfluxarr[i] = cheapflux

   ;rotatesurface,star,rotang

endfor

save,stararr,starvis,filename='stararr.sav'
save,eph,filename='eph.sav'
if keyword_set(barycalc) then save,vvecarr,filename='vvec_arr.sav'

print
print,'Rotation period: ',star.P,' d',format='(17A,F5.2,A2)'
print,'vsini:           ',star.veq,' km/s',format='(17A,F5.2,A5)'
print

if keyword_set(newprofile) then save,wave,specarr,filename='modelspectra.sav'
if keyword_set(vradprofile) then save,velarrgauss,velarrbary,velarrxcorr,filename='modelrv.sav' else restore,'modelrv.sav'

if linemodel eq 0 then velarrplot = velarrcheap else velarrplot = velarrterra
;velarrplot = velarrgauss
;velarrplot = velarrbary
;velarrplot = velarrxcorr[0,*]
;velarrplot = velarrterra

;restore,'~/Data/playground/Solar_Eclipse/eclchisq_arr.sav'
;velcorr_2 = median(velcorrarr[indgen(8),0,*],dimension=1)

;restore,'~/Data/playground/Sun_RV/eclchisq_arr2015-03-20.sav'
restore,'~/Data/playground/Sun_RV/eclchisq_arr_bruker2015-03-20.sav'
rvord = [indgen(12),15,16]
;restore,'~/Data/playground/Sun_RV/6160_eclchis2015-03-20.sav'
;rvord = 4
velcorr_2 = median(velcorrarr[[rvord],*],dimension=1)

;restore,'~/Data/playground/Solar_Eclipse/eclchisq_arr_phx.sav'
;i = 19
;velcorr_2 = velcorrarr[i,0,*]
;velcorr_2 = velrvarr[i,0,*]
;velcorr_2 = median(velrvarr[indgen(7),0,*],dimension=1)

restore,'~/Data/playground/Solar_Eclipse/obs_flux.sav'
mflux = fltarr(n)
cflux = fltarr(n)
for i = 0,n-1 do mflux[i] = total(specarr[i,*])
gain2p5ind = 150+indgen(9) & flux[gain2p5ind] = flux[gain2p5ind]*2.5
for i = 0,n-1 do cflux[i] = total(cheapfluxarr[i,*])


indref = 125 + indgen(30)
setcolors,/systemvariables,/silent
fiftmin = 0.010416667d
tint = 2.d0*fiftmin
xr = 2457102.d0+[-10.*tint,1.*tint]+0.0001
xrflux = 2457102.d0+[-9.*tint,0.*tint]+0.0001
dataplot = velcorr_2

offdata = median(dataplot[indref])
offbary = median(eph.barycorr[indref])
offplot= median(velarrplot[indref])

offset =  offdata - offplot
print,'Indref observed:  ',offdata*1.e3,' m/s',format='(A21, F7.2,A4)'
print,'Indref JPL bary:  ',offbary*1.e3,' m/s',format='(A21, F7.2,A4)'
print,'Indref calc bary: ',offplot*1.e3,' m/s',format='(A21, F7.2,A4)'

residual = 1.e3*(reform(dataplot)-velarrplot - offset)
if keyword_set(lincorr) then begin
   wset,2
   chunkres,residual,eph.obs_time,chunks,linsol,sigmaclip
   wset,0
   linsol[chunks[1]:*] = 1.d0 ;only correct the part before 1st contact
   residual = residual - linsol
   dataplot = dataplot - linsol/1.e3
endif
if keyword_set(ressm) then residual = smooth(residual,ressm)

velarrplot = velarrplot + offset

plot:

;plot rotation, blueshift, and limb darkening
wset,2
;differential rotation
DDTOR = !DPI/180.D0
Ar = drot_par[0] & Br = drot_par[1] & Cr = drot_par[2]
sidsyn = 0.9324 & Ar = Ar * sidsyn & Br = Br * sidsyn & Cr = Cr * sidsyn
lat = findgen(1001)/1000. * 90.
dlaw = (Ar + Br*SIN(lat*DDTOR)^2. + Cr*SIN(lat*DDTOR)^4.)/Ar
plot,lat,dlaw,ps=0,position=[.1,.7,.95,.95],title='Differential Rotation',xtitle='Latitude',ytitle='rel. rot.'

;convective blueshift
mu = findgen(101)/100.
x = acos(mu)/!Pi*180.
vbs_ce = (4.564 + x*(0.2909 + x*(0.1076 + x*(-0.004107 + x*4.072e-5))))/1000. < 0.3 ;Cegla 2014
vbs_be = blsh_par[0] * (2.d0/(1.d0 + exp((mu - blsh_par[1])*blsh_par[2])) - 1.d0 ) + blsh_par[3]

plot,mu,vbs_be,ps=0,xr=[1,.05],position=[.1,.35,.95,.6],yr=[-.5,.4],/noerase,title='Convective Blueshift',xtitle='mu',ytitle='velocity [km/s]'
;oplot,mu,vbs_ce,ps=0,col=!RED

;limb darkening
ld = poly(mu,pfit)
plot,mu,ld,ps=0,xr=[1,.05],position=[.1,.07,.95,.25],/noerase,title='Limb Darkening',xtitle='mu',ytitle='rel. brightness'

wset, 0

restore,'vel_noeclipse.sav'
indtransit = where( abs(velarrplot-velarrplot_noeclipse) gt .035)
if ecl eq 0 then begin
   indtransit = indref
   velarrplot_noeclipse = velarrplot
   save,velarrplot_noeclipse,filename='vel_noeclipse.sav'
endif

ps_open,'ConvBlueshift',/encapsulated,/color,thick=3
plot,mu,vbs_be,ps=0,xr=[1,.05],yr=minmax(vbs_be)+[-.1,.1],/noerase,xtitle='!4l!3',ytitle='Convective blueshift [km/s]',position=[.12,.15,.95,.6]
ps_close
spawn,'epstopdf ConvBlueshift.eps'


if keyword_set(ps) then ps_open,'Eclipse',/color,/encapsulated,thick=3 else wset,2

dummy = LABEL_DATE(DATE_FORMAT=[' '])
;plot RV
plot,  eph.obs_time, velarrplot,ps=4,yrange = [-440.,990.],/ystyle,xticks=11,/xstyle,xr=xr,position=[.1,.4,.95,.95],ytitle='RV [m/s]',xtickformat='LABEL_DATE',/nodata
oplot, eph.obs_time, (velarrplot)*1.e3, ps=0,col=!RED
;oplot, eph.obs_time, (velarrxcorr[0,*]-median(velarrxcorr[0,indref]))*1.e3, ps=0,col=!RED,thick=2
oplot, eph.obs_time, dataplot*1.e3,ps=1
;oplot, eph.obs_time, eph.barycorr*1.e3,col=!RED,ps=0

;plot residuals
;plot,  eph.obs_time,residual,ps=1,yrange=[-39.,39.],/xstyle,xr=xr,position=[.1,.4,.95,.6],/noerase,/ystyle,ytitle='O-C [m/s]',xtickformat='LABEL_DATE',xticks=11
;if keyword_set(lincorr) then oplot, eph.obs_time,linsol,ps=0,col=!GREEN
;oplot,[0,1.e8],[0.,0.],col=!GREY,linestyle=3

;plot airmass
plot, eph.obs_time, eph.airm,ps=1,/xstyle,yr=[1.,4.],position=[.1,.3,.95,.4],xr=xr,/noerase,ytitle='airmass',xtickformat='LABEL_DATE',xticks=11,/ystyle,yticks=3
;plot water
;oplot, eph.obs_time, eqwh2o/median(eqwh2o)*2.5,ps=1,color=!BLUE


;flux = flux[1:*]
dummy = LABEL_DATE(DATE_FORMAT=['%H:%I'])
;plot flux
mf = mean( flux[135:140]) 
cf = mean(cflux[135:140])
plot,  eph.obs_time, flux/mf,ps=1,xtickformat='LABEL_DATE',xticks=11,/xstyle,xr=xr,xtitle='Fri Mar 20, 2015',yr=[0,1.5],/noerase,position=[.1,.1,.95,.3],/ystyle,ytitle='relative flux';, /nodata
oplot, eph.obs_time, cflux/cf,ps=0,col=!RED,thick=2
;oplot, eph.obs_time, mf*mflux/max(mflux),ps=0,col=!RED,thick=2
;oplot, eph.obs_time, mf*(flux/max(flux))/0.95 / (cflux/max(cflux)),ps=4,col=!GREY

if keyword_set(ps) then begin
   ps_close
   spawn,'epstopdf Eclipse.eps'
endif


;same plots but in two figures

rms = 2.5
if keyword_set(ps) then ps_open,'Eclipse_RV',/color,/encapsulated,thick=3,/portrait else wset,2

dummy = LABEL_DATE(DATE_FORMAT=[' '])
;plot RV
plot,  eph.obs_time, dataplot,ps=4,yrange = [-100,100]+minmax((dataplot)*1.e3),/ystyle,xticks=11,/xstyle,xr=xr,position=[.15,.5,.95,.95],ytitle='RV [m/s]',xtickformat='LABEL_DATE',/nodata,charsize=1
oplot, eph.obs_time, (velarrplot)*1.e3, ps=0,col=!RED
;oplot, eph.obs_time, (velarrxcorr[0,*]-median(velarrxcorr[0,indref]))*1.e3, ps=0,col=!RED,thick=2
oplot, eph.obs_time, (dataplot)*1.e3,ps=1
;oplot, eph.obs_time, eph.barycorr*1.e3,col=!RED,ps=0
oplot, eph.obs_time, velarrplot_noeclipse*1.e3,ps=0,col=!GREY,linestyle=2
oplot, eph.obs_time, eph.barycorr*1.e3,ps=0,col=!BLUE

dummy = LABEL_DATE(DATE_FORMAT=['%H:%I'])
;plot residuals
plot,  eph.obs_time,residual,ps=1,yrange=[-30,110],/xstyle,xr=xr,position=[.15,.3,.95,.5],/noerase,/ystyle,ytitle='O-C [m/s]',xtickformat='LABEL_DATE',xticks=11,charsize=1,xtitle='Fri Mar 20, 2015'
;if not keyword_set(lincorr) then oplot, eph.obs_time,linsol,ps=0,col=!GREEN
;oploterror,eph.obs_time,residual,replicate(rms,n_elements(residual)),/nohat,errcolor=!GREY,ps=1
;oplot,[0,1.e8],[0.,0.],col=!GREY,linestyle=3

;relres = residual/1000./abs(velarrplot-velarrplot_noeclipse)*100. ;relative residual
;errres = (residual+rms)/1000./abs(velarrplot-velarrplot_noeclipse)*100. ;relative residual
;dummy = LABEL_DATE(DATE_FORMAT=['%H:%I'])
;if ecl ne 0 then plot,  eph.obs_time[indtransit],relres[indtransit],ps=1,yrange=1.2*minmax(relres[indtransit]),/xstyle,xr=xr,position=[.15,.1,.95,.3],/noerase,/ystyle,ytitle='O-C [%]',xtickformat='LABEL_DATE',xticks=11,xtitle='Fri Mar 20, 2015',charsize=1
;oplot,[0,1.e8],[0.,0.],col=!GREY,linestyle=3


if keyword_set(ps) then begin
   ps_close
   spawn,'epstopdf Eclipse_RV.eps'
endif


if keyword_set(ps) then ps_open,'Eclipse_flux',/color,/encapsulated,thick=3 else wset,2

xrflux=xr
xticksflux = 11
;xticksflux = 9

dummy = LABEL_DATE(DATE_FORMAT=[' '])
;plot flux
mf = mean( flux[135:140]) 
cf = mean(cflux[135:140])
plotsym,0,thick=2
plot,  eph.obs_time, flux/mf,ps=8,xtickformat='LABEL_DATE',xticks=xticksflux,/xstyle,xr=xrflux,yr=[0,1.2],/noerase,position=[.1,.5,.95,.95],/ystyle,ytitle='relative flux'
plotsym,0,/fill
oplot, eph.obs_time, cflux/cf,ps=0,col=!RED,thick=2
;oplot, eph.obs_time, cf/mf*flux/cflux,ps=1,col=!BLUE,thick=2
;oplot, eph.obs_time, mf*(flux/max(flux))/0.95 / (cflux/max(cflux)),ps=4,col=!GREY
xyouts,xrflux[0]+tint/2.,.28,'Observation'
xyouts,xrflux[0]+tint/2.,.60,'Model',col=!RED
;xyouts,xr[0]+tint,.92,'Observation / Model',col=!BLUE

;flux = flux[1:*]
dummy = LABEL_DATE(DATE_FORMAT=['%H:%I'])
;plot airmass
plot, eph.obs_time, eph.airm,ps=1,/xstyle,yr=[1.,4.4],position=[.1,.3,.95,.5],xr=xrflux,/noerase,ytitle='airmass',xtickformat='LABEL_DATE',xticks=xticksflux,/ystyle,xtitle='Fri Mar 20, 2015',ytickname=['1',' ','2',' ','3',' ','4']
;plot water
;oplot, eph.obs_time, eqwh2o/median(eqwh2o)*2.5,ps=1,color=!BLUE

if keyword_set(ps) then begin
   ps_close
   spawn,'epstopdf Eclipse_flux.eps'
endif





;print information
;print,'Spectral synthesis:'
;print,'mean: ',mean((reform(dataplot)-reform(velarrxcorr[0,*]) + median(velarrxcorr[0,indref]))[goodind])*1.e3, ' m/s',format='(A5,F5.2,A4)'
;print,'rms: ',stdev((reform(dataplot)-reform(velarrxcorr[0,*]) + median(velarrxcorr[0,indref]))[goodind])*1.e3, ' m/s',format='(A5,F5.2,A4)'
print,'Weighted velocity:'
print,'mean: ',mean(residual), ' m/s',format='(A5,F5.2,A4)'
print,'rms: ',stdev(residual), ' m/s',format='(A5,F5.2,A4)'
print,'Before transit:'
print,'rms: ',stdev(residual[chunks[0]:chunks[1]]), ' m/s',format='(A5,F5.2,A4)'
print,'After transit:'
print,'rms: ',stdev(residual[chunks[2]:*]), ' m/s',format='(A5,F5.2,A4)'
print,'Before and after transit:'
print,'rms: ',stdev([residual[chunks[0]:chunks[1]],residual[chunks[2]:*]]), ' m/s',format='(A5,F5.2,A4)'
print,'Before and after transit 5min average:'
print,'rms: ',stdev(rebin([residual[chunks[0]:chunks[1]],residual[chunks[2]:*]],27)), ' m/s',format='(A5,F5.2,A4)'



if ecl eq 1 then print,'Maximum fraction eclipsed: ',maxecl*100.,'%',format='(A27, F6.1,A1)'

end


