Real Function x33pi(srt)
  Real xarray(12), earray(12)
  Save
  Data xarray/0.02, 0.22, 0.74, 1.10, 1.76, 1.84, 2.20, 2.40, 2.15, 2.60, 2.30, 1.70/
  Data earray/2.23, 2.81, 3.67, 4.00, 4.95, 5.52, 5.97, 6.04, 6.60, 6.90, 10.01, 19./
  pmass = 0.9383
  x33pi = 1.E-06
  If (srt<=2.3) Return
  plab = sqrt(((srt**2-2.*pmass**2)/(2.*pmass))**2-pmass**2)
  If (plab<earray(1)) Then
     x33pi = xarray(1)
     Return
  End If
  Do ie = 1, 12
     If (earray(ie)==plab) Then
        x33pi = xarray(ie)
        Return
     Else If (earray(ie)>plab) Then
        ymin = alog(xarray(ie-1))
        ymax = alog(xarray(ie))
        xmin = alog(earray(ie-1))
        xmax = alog(earray(ie))
        x33pi = exp(ymin+(alog(plab)-xmin)*(ymax-ymin)/(xmax-xmin))
        Return
     End If
  End Do
  Return
End Function x33pi
