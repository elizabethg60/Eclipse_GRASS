
pro extinction,star,eph,tid,extcof

;procedure to compute extinction on visible solar disk
;fills arrays star.attn and star.airm
  
DDTOR = !Pi/180.d0

rotangle = + (eph.sunrota)[tid] - (eph.pa)[tid]

horix =  cos(rotangle)*star.xproj - sin(rotangle)*star.yproj
horiy =  sin(rotangle)*star.xproj + cos(rotangle)*star.yproj

elevdisc = (eph.eleap)[tid] + horiy * (eph.Rs)[tid] 
z = 90.d0 - elevdisc
;star.airm = 1.d0/cos(DDTOR*z)
star.airm = 1.d0/(cos(DDTOR*z) + 0.025*exp(-11.*cos(DDTOR*z)))

dmag = extcof * star.airm
star.attn = 10.d0^(-dmag/2.5d0)


end
