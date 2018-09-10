Real Function fkaon(p, pmax)
  Save
  fmax = 0.148
  If (pmax==0.) pmax = 0.000001
  fkaon = (1.-p/pmax)*(p/pmax)**2
  If (fkaon>fmax) fkaon = fmax
  fkaon = fkaon/fmax
  Return
End Function fkaon
