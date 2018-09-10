Real Function pionpp(srt)
  Save
  pmass = 0.14
  pmass1 = 0.938
  pionpp = 0.00001
  If (srt<=1.22) Return
  plab = sqrt(((srt**2-pmass**2-pmass1**2)/(2.*pmass1))**2-pmass**2)
  pmin = 0.3
  pmax = 25.0
  If (plab>pmax) Then
    pionpp = 20./10.
    Return
  End If
  If (plab<pmin) Then
    pionpp = 0.
    Return
  End If
  a = 24.3
  b = -12.3
  c = 0.324
  an = -1.91
  d = -2.44
  pionpp = a + b*(plab**an) + c*(alog(plab))**2 + d*alog(plab)
  If (pionpp<=0) pionpp = 0
  pionpp = pionpp/10.
  Return
End Function pionpp
