Real Function aknpsg(pkaon)
  Save
  If (pkaon<=0.345) Then
    sigma1 = 0.624*pkaon**(-1.83)
  Else
    sigma1 = 0.7*pkaon**(-2.09)
  End If
  aknpsg = sigma1
  Return
End Function aknpsg
