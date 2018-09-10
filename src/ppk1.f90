Real Function ppk1(srt)
  Real xarray(7), earray(7)
  Save
  Data xarray/0.013, 0.025, 0.016, 0.012, 0.017, 0.029, 0.025/
  Data earray/3.67, 4.95, 5.52, 5.97, 6.05, 6.92, 7.87/
  pmass = 0.9383
  ppk1 = 0.
  If (srt<=2.63) Return
  If (srt>4.08) Then
     ppk1 = 0.025
     Return
  End If
  plab = sqrt(((srt**2-2.*pmass**2)/(2.*pmass))**2-pmass**2)
  If (plab<earray(1)) Then
     ppk1 = xarray(1)
     Return
  End If
  Do ie = 1, 7
     If (earray(ie)==plab) Then
        ppk1 = xarray(ie)
        Goto 10
     Else If (earray(ie)>plab) Then
        ymin = alog(xarray(ie-1))
        ymax = alog(xarray(ie))
        xmin = alog(earray(ie-1))
        xmax = alog(earray(ie))
        ppk1 = exp(ymin+(alog(plab)-xmin)*(ymax-ymin)/(xmax-xmin))
        Goto 10
     End If
  End Do
10 Continue
  Return
End Function ppk1
