pro cheapspec, star, vel, cheapflux

  ;similar to line weakening?
  ;star.contrast = (star.limbdark/max(star.limbdark,/nan))^0.1

  brightness = star.limbdark*star.attn*star.aproj
  cheapflux = total((brightness)[*star.vis],/double) 
  vel = total((star.contrast*star.vproj*brightness)[*star.vis],/double) / cheapflux

end
