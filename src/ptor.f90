Real Function ptor(srt)
  Parameter (pimass=0.140, arho=0.77)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
  Save
  s2 = srt**2
  ptor = 9*(s2-4*arho**2)/(s2-4*pimass**2)*rtop(srt)
  Return
End Function ptor
