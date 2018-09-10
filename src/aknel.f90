Real Function aknel(pkaon)
  Save
  If (pkaon<0.5 .Or. pkaon>=4.0) sigma1 = 0.
  If (pkaon>=0.5 .And. pkaon<1.0) sigma1 = 20.*pkaon**2.74
  If (pkaon>=1.0 .And. pkaon<4.0) sigma1 = 20.*pkaon**(-1.8)
  aknel = sigma1
  Return
End Function aknel
