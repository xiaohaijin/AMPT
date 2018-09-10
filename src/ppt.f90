Real Function ppt(srt)
  Save
  pmass = 0.9383
  scheck = ((srt**2-2.*pmass**2)/(2.*pmass))**2 - pmass**2
  If (scheck<0) Then
     Write (99, *) 'scheck34: ', scheck
     scheck = 0.
  End If
  plab = sqrt(scheck)
  pmin = 3.
  pmax = 2100
  If ((plab<pmin) .Or. (plab>pmax)) Then
     ppt = 55.
     Return
  End If
  a = 45.6
  b = 219.0
  c = 0.410
  d = -3.41
  an = -4.23
  ppt = a + b*(plab**an) + c*(alog(plab))**2 + d*alog(plab)
  If (ppt<=0) ppt = 0.0
  Return
End Function ppt
