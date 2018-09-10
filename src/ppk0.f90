Real Function ppk0(srt)
  Real xarray(7), earray(7)
  Save
  Data xarray/0.030, 0.025, 0.025, 0.026, 0.02, 0.014, 0.06/
  Data earray/3.67, 4.95, 5.52, 6.05, 6.92, 7.87, 10./
  pmass = 0.9383
  ppk0 = 0
  If (srt<=2.63) Return
  If (srt>4.54) Then
     ppk0 = 0.037
     Return
  End If
  plab = sqrt(((srt**2-2.*pmass**2)/(2.*pmass))**2-pmass**2)
  If (plab<earray(1)) Then
     ppk0 = xarray(1)
     Return
  End If
  Do ie = 1, 7
     If (earray(ie)==plab) Then
        ppk0 = xarray(ie)
        Goto 10
     Else If (earray(ie)>plab) Then
        ymin = alog(xarray(ie-1))
        ymax = alog(xarray(ie))
        xmin = alog(earray(ie-1))
        xmax = alog(earray(ie))
        ppk0 = exp(ymin+(alog(plab)-xmin)*(ymax-ymin)/(xmax-xmin))
        Goto 10
     End If
  End Do
10 Continue
  Return
End Function ppk0
