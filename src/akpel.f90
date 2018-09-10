Real Function akpel(pkaon)
  Save
  If (pkaon<0.25 .Or. pkaon>=4.0) sigma2 = 0.
  If (pkaon>=0.25 .And. pkaon<4.0) sigma2 = 13.*pkaon**(-0.9)
  akpel = sigma2
  Return
End Function akpel
