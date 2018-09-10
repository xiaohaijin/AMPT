Real Function prbbar(srt)
  Parameter (pimass=0.140, arho=0.77, aomega=0.782)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
  Save
  sppb3p = xppbar(srt)*factr2(3)*ene/fsum
  pi2 = (s-(pimass+arho)**2)*(s-(pimass-arho)**2)/4/s
  prbbar = 4./27.*sppb3p/pi2*wtot
  Return
End Function prbbar
