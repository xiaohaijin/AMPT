Real Function pipik(srt)
  Real xarray(5), earray(5)
  Save
  Data xarray/0.001, 0.7, 1.5, 1.7, 2.0/
  Data earray/1., 1.2, 1.6, 2.0, 2.4/
  pmass = 0.9383
  pipik = 0.
  If (srt<=1.) Return
  If (srt>2.4) Then
    pipik = 2.0/2.
    Return
  End If
  If (srt<earray(1)) Then
    pipik = xarray(1)/2.
    Return
  End If
  Do ie = 1, 5
    If (earray(ie)==srt) Then
      pipik = xarray(ie)
      Goto 10
    Else If (earray(ie)>srt) Then
      ymin = alog(xarray(ie-1))
      ymax = alog(xarray(ie))
      xmin = alog(earray(ie-1))
      xmax = alog(earray(ie))
      pipik = exp(ymin+(alog(srt)-xmin)*(ymax-ymin)/(xmax-xmin))
      Goto 10
    End If
  End Do
  10 pipik = pipik/2.
  Continue
  Return
End Function pipik
