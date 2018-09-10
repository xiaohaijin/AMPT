Real Function pipp1(srt)
  Save
  pmass = 0.14
  pmass1 = 0.938
  pipp1 = 0.0001
  If (srt<=1.22) Return
  plab = sqrt(((srt**2-pmass**2-pmass1**2)/(2.*pmass1))**2-pmass**2)
  pmin = 0.3
  pmax = 25.0
  If (plab>pmax) Then
    pipp1 = 20./10.
    Return
  End If
  If (plab<pmin) Then
    pipp1 = 0.
    Return
  End If
  a = 26.6
  b = -7.18
  c = 0.327
  an = -1.86
  d = -2.81
  pipp1 = a + b*(plab**an) + c*(alog(plab))**2 + d*alog(plab)
  If (pipp1<=0) pipp1 = 0
  pipp1 = pipp1/10.
  Return
End Function pipp1
