Real Function akpsgm(pkaon)
  Save
  If (pkaon<0.2 .Or. pkaon>=1.5) sigma1 = 0.
  If (pkaon>=0.2 .And. pkaon<1.5) sigma1 = 0.6*pkaon**(-1.8)
  akpsgm = sigma1
  Return
End Function akpsgm
