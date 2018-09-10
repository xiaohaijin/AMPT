Real Function akplam(pkaon)
  Save
  p = pkaon
  If (pkaon<0.2 .Or. pkaon>=10.0) sigma = 0.
  If (pkaon>=0.2 .And. pkaon<0.9) sigma = 50.*p**2 - 67.*p + 24.
  If (pkaon>=0.9 .And. pkaon<10.0) sigma = 3.0*pkaon**(-2.6)
  akplam = sigma
  Return
End Function akplam
