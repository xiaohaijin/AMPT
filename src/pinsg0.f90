Real Function pinsg0(srt)
  Save
  srt0 = 0.938 + 2.*0.498
  If (srt<srt0) Then
    pinsg0 = 0.0
    Return
  End If
  ratio = srt0**2/srt**2
  pinsg0 = 1.121*(1.-ratio)**1.86*ratio**2
  Return
End Function pinsg0
