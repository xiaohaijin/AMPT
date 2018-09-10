Real Function pp2(srt)
  Save
  pmass = 0.9383
  scheck = ((srt**2-2.*pmass**2)/(2.*pmass))**2 - pmass**2
  If (scheck<0) Then
     Write (99, *) 'scheck33: ', scheck
     scheck = 0.
  End If
  plab = sqrt(scheck)
  pmin = 2.
  pmax = 2050
  If (plab>pmax) Then
     pp2 = 8.
     Return
  End If
  If (plab<pmin) Then
     pp2 = 25.
     Return
  End If
  a = 11.2
  b = 25.5
  c = 0.151
  d = -1.62
  an = -1.12
  pp2 = a + b*(plab**an) + c*(alog(plab))**2 + d*alog(plab)
  If (pp2<=0) pp2 = 0
  Return
End Function pp2
