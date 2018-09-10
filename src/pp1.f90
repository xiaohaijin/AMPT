Real Function pp1(srt)
  Save
  pmass = 0.9383
  pp1 = 0.
  plab2 = ((srt**2-2.*pmass**2)/(2.*pmass))**2 - pmass**2
  If (plab2<=0) Return
  plab = sqrt(plab2)
  pmin = 0.968
  pmax = 2080
  If ((plab<pmin) .Or. (plab>pmax)) Then
     pp1 = 0.
     Return
  End If
  a = 30.9
  b = -28.9
  c = 0.192
  d = -0.835
  an = -2.46
  pp1 = a + b*(plab**an) + c*(alog(plab))**2
  If (pp1<=0) pp1 = 0.0
  Return
End Function pp1
