Real Function robbar(srt)
  Parameter (pimass=0.140, arho=0.77, aomega=0.782)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
  Save
  sppb5p = xppbar(srt)*factr2(5)*ene**3/fsum
  pi2 = (s-(arho+aomega)**2)*(s-(arho-aomega)**2)/4/s
  robbar = 4./27.*sppb5p/pi2*wtot
  Return
End Function robbar
