Real Function ppbbar(srt)
  Parameter (pimass=0.140, arho=0.77, aomega=0.782)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
  Save
  sppb2p = xppbar(srt)*factr2(2)/fsum
  pi2 = (s-4*pimass**2)/4
  ppbbar = 4./9.*sppb2p/pi2*wtot
  Return
End Function ppbbar
