Real Function pptope(srt)
  Parameter (pimass=0.140, etam=0.5475)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
  Save
  s2 = srt**2
  pf2 = (s2-(pimass+etam)**2)*(s2-(pimass-etam)**2)/2/sqrt(s2)
  pi2 = (s2-4*pimass**2)*s2/2/sqrt(s2)
  pptope = 1./3.*pf2/pi2*petopp(srt)
  Return
End Function pptope
