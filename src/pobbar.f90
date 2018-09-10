Real Function pobbar(srt)
  Parameter (pimass=0.140, arho=0.77, aomega=0.782)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
  Save
  sppb4p = xppbar(srt)*factr2(4)*ene**2/fsum
  pi2 = (s-(pimass+aomega)**2)*(s-(pimass-aomega)**2)/4/s
  pobbar = 4./9.*(sppb4p/2)/pi2*wtot
  Return
End Function pobbar
