Real Function xpp(srt)
  Real xarray(14), earray(14)
  Save
  Data earray/20., 30., 40., 60., 80., 100., 170., 250., 310., 350., 460., 560., 660., 800./
  Data xarray/150., 90., 80.6, 48.0, 36.6, 31.6, 25.9, 24.0, 23.1, 24.0, 28.3, 33.6, 41.5, 47/
  pmass = 0.9383
  ekin = 2000.*pmass*((srt/(2.*pmass))**2-1.)
  If (ekin<earray(1)) Then
    xpp = xarray(1)
    If (xpp>55) xpp = 55
    Return
  End If
  If (ekin>earray(14)) Then
    xpp = xarray(14)
    Return
  End If
  Do ie = 1, 14
    If (earray(ie)==ekin) Then
      xpp = xarray(ie)
      If (xpp>55) xpp = 55.
      Return
    End If
    If (earray(ie)>ekin) Then
      ymin = alog(xarray(ie-1))
      ymax = alog(xarray(ie))
      xmin = alog(earray(ie-1))
      xmax = alog(earray(ie))
      xpp = exp(ymin+(alog(ekin)-xmin)*(ymax-ymin)/(xmax-xmin))
      If (xpp>55) xpp = 55.
      Goto 50
    End If
  End Do
  50 Continue
  Return
End Function xpp
