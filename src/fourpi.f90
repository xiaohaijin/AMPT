Real Function fourpi(srt)
  Real xarray(10), earray(10)
  Save
  Data xarray/0.0001, 1.986597, 6.411932, 7.636956, 9.598362, 9.889740, 10.24317, 10.80138, 11.86988, 12.83925/
  Data earray/2.468, 2.718, 2.968, 0.322E+01, 0.347E+01, 0.372E+01, 0.397E+01, 0.422E+01, 0.447E+01, 0.472E+01/
  pmass = 0.14
  pmass1 = 0.938
  fourpi = 0.000001
  If (srt<=1.52) Return
  plab = srt
  If (plab<earray(1)) Then
    fourpi = 0.00001
    Return
  End If
  Do ie = 1, 10
    If (earray(ie)==plab) Then
      fourpi = xarray(ie)
      Return
    Else If (earray(ie)>plab) Then
      ymin = alog(xarray(ie-1))
      ymax = alog(xarray(ie))
      xmin = alog(earray(ie-1))
      xmax = alog(earray(ie))
      fourpi = exp(ymin+(alog(plab)-xmin)*(ymax-ymin)/(xmax-xmin))
      Return
    End If
  End Do
  Return
End Function fourpi
