Real Function xnp(srt)
  Real xarray(11), earray(11)
  Save
  Data earray/20., 30., 40., 60., 90., 135.0, 200., 300., 400., 600., 800./
  Data xarray/410., 270., 214.5, 130., 78., 53.5, 41.6, 35.9, 34.2, 34.3, 34.9/
  pmass = 0.9383
  ekin = 2000.*pmass*((srt/(2.*pmass))**2-1.)
  If (ekin<earray(1)) Then
    xnp = xarray(1)
    If (xnp>55) xnp = 55
    Return
  End If
  If (ekin>earray(11)) Then
    xnp = xarray(11)
    Return
  End If
  Do ie = 1, 11
    If (earray(ie)==ekin) Then
      xnp = xarray(ie)
      If (xnp>55) xnp = 55.
      Return
    End If
    If (earray(ie)>ekin) Then
      ymin = alog(xarray(ie-1))
      ymax = alog(xarray(ie))
      xmin = alog(earray(ie-1))
      xmax = alog(earray(ie))
      xnp = exp(ymin+(alog(ekin)-xmin)*(ymax-ymin)/(xmax-xmin))
      If (xnp>55) xnp = 55
      Goto 50
    End If
  End Do
  50 Continue
  Return
End Function xnp
