Real Function pplpk(srt)
  Save
  pmass = 0.9383
  pplpk = 0.
  scheck = ((srt**2-2.*pmass**2)/(2.*pmass))**2 - pmass**2
  If (scheck<0) Then
     Write (99, *) 'scheck35: ', scheck
     scheck = 0.
  End If
  plab = sqrt(scheck)
  pmin = 2.82
  pmax = 25.0
  If (plab>pmax) Then
     pplpk = 0.036
     Return
  End If
  If (plab<pmin) Then
     pplpk = 0.
     Return
  End If
  a = 0.0654
  b = -3.16
  c = -0.0029
  an = -4.14
  pplpk = a + b*(plab**an) + c*(alog(plab))**2
  If (pplpk<=0) pplpk = 0
  Return
End Function pplpk
