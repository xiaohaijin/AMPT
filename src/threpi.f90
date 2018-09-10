Real Function threpi(srt)
  Real xarray(15), earray(15)
  Save
  Data xarray/8.0000000E-06, 6.1999999E-05, 1.881940, 5.025690, 11.80154, 13.92114, 15.07308, 11.79571, 11.53772, 10.01197, 9.792673, 9.465264, 8.970490, 7.944254, 6.886320/
  Data earray/0.122E+01, 0.147E+01, 0.172E+01, 0.197E+01, 0.222E+01, 0.247E+01, 0.272E+01, 0.297E+01, 0.322E+01, 0.347E+01, 0.372E+01, 0.397E+01, 0.422E+01, 0.447E+01, 0.472E+01/
  pmass = 0.14
  pmass1 = 0.938
  threpi = 0.000001
  If (srt<=1.36) Return
  plab = srt
  If (plab<earray(1)) Then
    threpi = 0.00001
    Return
  End If
  Do ie = 1, 15
    If (earray(ie)==plab) Then
      threpi = xarray(ie)
      Return
    Else If (earray(ie)>plab) Then
      ymin = alog(xarray(ie-1))
      ymax = alog(xarray(ie))
      xmin = alog(earray(ie-1))
      xmax = alog(earray(ie))
      threpi = exp(ymin+(alog(plab)-xmin)*(ymax-ymin)/(xmax-xmin))
      Return
    End If
  End Do
  Return
End Function threpi
