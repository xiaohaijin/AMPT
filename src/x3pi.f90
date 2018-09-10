Real Function x3pi(srt)
  Real xarray(12), earray(12)
  Save
  Data xarray/0.02, 0.4, 1.15, 1.60, 2.19, 2.85, 2.30, 3.10, 2.47, 2.60, 2.40, 1.70/
  Data earray/2.23, 2.81, 3.67, 4.00, 4.95, 5.52, 5.97, 6.04, 6.60, 6.90, 10.01, 19./
  pmass = 0.9383
  x3pi = 1.E-06
  If (srt<=2.3) Return
  plab = sqrt(((srt**2-2.*pmass**2)/(2.*pmass))**2-pmass**2)
  If (plab<earray(1)) Then
     x3pi = xarray(1)
     Return
  End If
  Do ie = 1, 12
     If (earray(ie)==plab) Then
        x3pi = xarray(ie)
        Return
     Else If (earray(ie)>plab) Then
        ymin = alog(xarray(ie-1))
        ymax = alog(xarray(ie))
        xmin = alog(earray(ie-1))
        xmax = alog(earray(ie))
        x3pi = exp(ymin+(alog(plab)-xmin)*(ymax-ymin)/(xmax-xmin))
        Return
     End If
  End Do
  Return
End Function x3pi
