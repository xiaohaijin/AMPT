Real Function twopi(srt)
  Real xarray(19), earray(19)
  Save
  Data xarray/0.300E-05, 0.187E+01, 0.110E+02, 0.149E+02, 0.935E+01, 0.765E+01, 0.462E+01, 0.345E+01, 0.241E+01, 0.185E+01, 0.165E+01, 0.150E+01, 0.132E+01, 0.117E+01, 0.116E+01, 0.100E+01, 0.856E+00, 0.745E+00, 0.300E-05/
  Data earray/0.122E+01, 0.147E+01, 0.172E+01, 0.197E+01, 0.222E+01, 0.247E+01, 0.272E+01, 0.297E+01, 0.322E+01, 0.347E+01, 0.372E+01, 0.397E+01, 0.422E+01, 0.447E+01, 0.472E+01, 0.497E+01, 0.522E+01, 0.547E+01, 0.572E+01/
  pmass = 0.14
  pmass1 = 0.938
  twopi = 0.000001
  If (srt<=1.22) Return
  plab = srt
  If (plab<earray(1)) Then
    twopi = 0.00001
    Return
  End If
  Do ie = 1, 19
    If (earray(ie)==plab) Then
      twopi = xarray(ie)
      Return
    Else If (earray(ie)>plab) Then
      ymin = alog(xarray(ie-1))
      ymax = alog(xarray(ie))
      xmin = alog(earray(ie-1))
      xmax = alog(earray(ie))
      twopi = exp(ymin+(alog(plab)-xmin)*(ymax-ymin)/(xmax-xmin))
      Return
    End If
  End Do
  Return
End Function twopi
