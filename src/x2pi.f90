Real Function x2pi(srt)
  Real xarray(15), earray(15)
  Save
  Data earray/2.23, 2.81, 3.67, 4.0, 4.95, 5.52, 5.97, 6.04, 6.6, 6.9, 7.87, 8.11, 10.01, 16.0, 19./
  Data xarray/1.22, 2.51, 2.67, 2.95, 2.96, 2.84, 2.8, 3.2, 2.7, 3.0, 2.54, 2.46, 2.4, 1.66, 1.5/
  pmass = 0.9383
  x2pi = 0.000001
  If (srt<=2.2) Return
  plab = sqrt(((srt**2-2.*pmass**2)/(2.*pmass))**2-pmass**2)
  If (plab<earray(1)) Then
     x2pi = xarray(1)
     Return
  End If
  Do ie = 1, 15
     If (earray(ie)==plab) Then
        x2pi = xarray(ie)
        Return
     Else If (earray(ie)>plab) Then
        ymin = alog(xarray(ie-1))
        ymax = alog(xarray(ie))
        xmin = alog(earray(ie-1))
        xmax = alog(earray(ie))
        x2pi = exp(ymin+(alog(plab)-xmin)*(ymax-ymin)/(xmax-xmin))
        Return
     End If
  End Do
  Return
End Function x2pi
