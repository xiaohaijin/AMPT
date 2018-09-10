Real Function aknsgm(pkaon)
  Save
  If (pkaon<0.5 .Or. pkaon>=6.0) sigma2 = 0.
  If (pkaon>=0.5 .And. pkaon<1.0) sigma2 = 1.2*pkaon**(-1.3)
  If (pkaon>=1.0 .And. pkaon<6.0) sigma2 = 1.2*pkaon**(-2.3)
  aknsgm = sigma2
  Return
End Function aknsgm
