Real Function xop2oe(srt)
  Parameter (pimass=0.140, etam=0.5475, aomega=0.782)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
  Save
  s2 = srt**2
  pf2 = (s2-(aomega+etam)**2)*(s2-(aomega-etam)**2)/2/sqrt(s2)
  pi2 = (s2-(aomega+pimass)**2)*(s2-(aomega-pimass)**2)/2/sqrt(s2)
  xop2oe = 1./3.*pf2/pi2*xoe2op(srt)
  Return
End Function xop2oe
