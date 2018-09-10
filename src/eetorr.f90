Real Function eetorr(srt)
  Parameter (etam=0.5475, arho=0.77)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
  Save
  s2 = srt**2
  eetorr = 81.*(s2-4*arho**2)/(s2-4*etam**2)*rrtoee(srt)
  Return
End Function eetorr
